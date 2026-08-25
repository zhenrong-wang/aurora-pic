#!/usr/bin/env python3
"""Analyze the prospective AuroraPIC/native-eduPIC particle histories."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


HISTORY_COLUMNS = {
    "age_steps": "tail_mean_age_steps",
    "energetic_steps": "tail_mean_energetic_steps",
    "energetic_duty_fraction": "tail_mean_energetic_duty_fraction",
    "consecutive_energetic_steps":
        "tail_mean_consecutive_energetic_steps",
    "entries": "tail_mean_entries",
    "elastic_collisions": "tail_mean_elastic_collisions",
    "excitation_collisions": "tail_mean_excitation_collisions",
    "ionization_collisions": "tail_mean_ionization_collisions",
    "charge_exchange_collisions":
        "tail_mean_charge_exchange_collisions",
    "bgk_collisions": "tail_mean_bgk_collisions",
    "born_during_window_fraction": "tail_born_during_window_fraction",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def row_key(row: dict[str, str]) -> tuple[int, str]:
    return int(row["phase_bin"]), row["region"]


def native_rows(directory: Path) -> list[dict[str, str]]:
    moments = {row_key(row): row for row in read_rows(
        directory / "edupic_phase_eedf_moments.csv")}
    result = []
    for history in read_rows(directory / "edupic_phase_eedf_history.csv"):
        merged = dict(moments[row_key(history)])
        merged.update(history)
        result.append(merged)
    if len(result) != len(moments):
        raise ValueError("native moment/history row shapes differ")
    return result


def select(rows: list[dict[str, str]], regions: set[str], lower: float,
           upper: float) -> list[dict[str, str]]:
    selected = [row for row in rows if row["region"] in regions and
                lower <= float(row["phase_fraction"]) < upper]
    if not selected:
        raise ValueError("particle-history selection is empty")
    return selected


def aggregate(rows: list[dict[str, str]], regions: set[str], lower: float,
              upper: float) -> dict[str, float | int]:
    selected = select(rows, regions, lower, upper)
    represented_key = ("represented_observations" if
                       "represented_observations" in selected[0] else
                       "represented_observations_m-2")
    tail_key = ("tail_represented_observations" if
                "tail_represented_observations" in selected[0] else
                "tail_represented_observations_m-2")
    macro = [int(row["macro_observations"]) for row in selected]
    represented = [float(row[represented_key]) for row in selected]
    tail_represented = [float(row[tail_key]) for row in selected]
    macro_weights = [represented[index] / macro[index] if macro[index] else 0.0
                     for index in range(len(selected))]
    if "tail_macro_observations" in selected[0]:
        tail_macro = [int(row["tail_macro_observations"])
                      for row in selected]
        history_weights = [float(value) for value in tail_macro]
    else:
        tail_macro = [int(round(tail_represented[index] /
                                macro_weights[index]))
                      if macro_weights[index] else 0
                      for index in range(len(selected))]
        history_weights = tail_represented
    history_denominator = sum(history_weights)
    if history_denominator <= 0.0:
        raise ValueError("particle-history selection has no energetic tail")
    result: dict[str, float | int] = {
        "rows": len(selected),
        "minimum_macro_observations_per_region_phase_bin": min(macro),
        "tail_macro_observations": sum(tail_macro),
        "maximum_overflow_fraction": max(float(row["overflow_fraction"])
                                          for row in selected),
    }
    for name, column in HISTORY_COLUMNS.items():
        result[name] = (
            sum(float(row[column]) * weight
                for row, weight in zip(selected, history_weights)) /
            history_denominator)
    age = float(result["age_steps"])
    if age <= 0.0:
        raise ValueError("particle-history aggregate has zero mean age")
    result.update({
        "current_streak_fraction":
            float(result["consecutive_energetic_steps"]) / age,
        "tail_entries_per_1000_age_steps":
            1000.0 * float(result["entries"]) / age,
        "elastic_collisions_per_1000_age_steps":
            1000.0 * float(result["elastic_collisions"]) / age,
        "excitation_collisions_per_1000_age_steps":
            1000.0 * float(result["excitation_collisions"]) / age,
        "ionization_collisions_per_1000_age_steps":
            1000.0 * float(result["ionization_collisions"]) / age,
        "charge_exchange_collisions_per_1000_age_steps":
            1000.0 * float(result["charge_exchange_collisions"]) / age,
        "bgk_collisions_per_1000_age_steps":
            1000.0 * float(result["bgk_collisions"]) / age,
    })
    return result


def mean_members(members: list[dict[str, float | int]],
                 metrics: tuple[str, ...]) -> dict[str, float]:
    return {metric: sum(float(member[metric]) for member in members) /
            len(members) for metric in metrics}


def relative_ratio(candidate: float, reference: float) -> float | None:
    return candidate / reference if reference != 0.0 else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("aurorapic_moments", type=Path)
    parser.add_argument("native_directories", nargs=3, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rule = json.loads(args.rule.read_text(encoding="utf-8"))
    diagnostic = rule["diagnostic_contract"]
    regions = set(diagnostic["critical_regions"])
    lower, upper = map(float, diagnostic["critical_phase_fraction"])
    candidate_rows = read_rows(args.aurorapic_moments)
    members_rows = [native_rows(path) for path in args.native_directories]
    candidate = aggregate(candidate_rows, regions, lower, upper)
    members = [aggregate(rows, regions, lower, upper)
               for rows in members_rows]
    metrics = tuple(HISTORY_COLUMNS) + (
        "current_streak_fraction", "tail_entries_per_1000_age_steps",
        "elastic_collisions_per_1000_age_steps",
        "excitation_collisions_per_1000_age_steps",
        "ionization_collisions_per_1000_age_steps",
        "charge_exchange_collisions_per_1000_age_steps",
        "bgk_collisions_per_1000_age_steps",
    )
    native = mean_members(members, metrics)
    ratios = {metric: relative_ratio(float(candidate[metric]), native[metric])
              for metric in metrics}
    differences = {metric: float(candidate[metric]) - native[metric]
                   for metric in metrics}
    ranges = {metric: ((max(float(member[metric]) for member in members) -
                        min(float(member[metric]) for member in members)) /
                       max(abs(native[metric]), 1e-300))
              for metric in metrics}
    born_range = (max(float(member["born_during_window_fraction"])
                      for member in members) -
                  min(float(member["born_during_window_fraction"])
                      for member in members))

    persistence = (
        abs(differences["energetic_duty_fraction"]) >= 0.03 or
        ratios["current_streak_fraction"] is not None and
        abs(float(ratios["current_streak_fraction"]) - 1.0) >= 0.15)
    exposure_metrics = (
        "elastic_collisions_per_1000_age_steps",
        "excitation_collisions_per_1000_age_steps",
        "ionization_collisions_per_1000_age_steps",
    )
    exposure = any(ratios[name] is not None and
                   abs(float(ratios[name]) - 1.0) >= 0.15
                   for name in exposure_metrics)
    turnover = (
        ratios["tail_entries_per_1000_age_steps"] is not None and
        abs(float(ratios["tail_entries_per_1000_age_steps"]) - 1.0) >= 0.15 or
        abs(differences["born_during_window_fraction"]) >= 0.02)

    repeatability_primary = all(ranges[name] <= 0.15 for name in (
        "age_steps", "energetic_duty_fraction", "current_streak_fraction",
        "tail_entries_per_1000_age_steps",
        "elastic_collisions_per_1000_age_steps"))
    repeatability_inelastic = all(ranges[name] <= 0.30 for name in (
        "excitation_collisions_per_1000_age_steps",
        "ionization_collisions_per_1000_age_steps"))
    expected_checkpoints = {
        "13507": "7b88021958430ed2157cc0dba4b3d20fd23ec03593c1e3f29b4c66ba540b13cd",
        "24601": "8935e4313ddef5c171a55f4631c9aab2665a772e799a20edaa297ee9ec78c2fd",
        "35713": "92bd734fc3710139ae9c212bce23f2f7559bb7e5a17599dca21e0b9a4eacf267",
    }
    passivity = {directory.name.removeprefix("seed-"):
                  sha256(directory / "picdata.bin") == expected_checkpoints[
                      directory.name.removeprefix("seed-")]
                  for directory in args.native_directories}
    population = [candidate, *members]
    gates = {
        "critical_observation_population": all(
            int(value["minimum_macro_observations_per_region_phase_bin"]) >=
            int(diagnostic[
                "minimum_macro_observations_per_critical_region_phase_bin"])
            for value in population),
        "critical_tail_population": all(
            int(value["tail_macro_observations"]) >=
            int(diagnostic["minimum_tail_macro_observations_critical_aggregate"])
            for value in population),
        "histogram_overflow": all(
            float(value["maximum_overflow_fraction"]) <=
            float(diagnostic["maximum_overflow_fraction"])
            for value in population),
        "finite_metrics": all(math.isfinite(float(value[metric]))
                              for value in population for metric in metrics),
        "native_diagnostic_passivity": all(passivity.values()),
        "native_repeatability_primary": repeatability_primary,
        "native_repeatability_inelastic": repeatability_inelastic,
        "native_repeatability_born_fraction": born_range <= 0.01,
        "aurorapic_history_enabled": all(
            int(row.get("history_enabled", "0")) == 1
            for row in candidate_rows),
    }
    result = {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "native_edupic_to_aurorapic_energetic_particle_history_result",
        "rule_sha256": sha256(args.rule),
        "gates": gates,
        "all_measurement_gates_passed": all(gates.values()),
        "native_passivity": passivity,
        "critical_phase_0p125_to_0p5": {
            "aurorapic": candidate,
            "native_edupic_members": members,
            "native_edupic_ensemble_mean": native,
            "native_edupic_relative_range": ranges,
            "native_edupic_born_fraction_absolute_range": born_range,
            "aurorapic_to_native_ratio": ratios,
            "aurorapic_minus_native": differences,
        },
        "prospective_decision_outcome": {
            "energetic_persistence_supported": persistence,
            "collision_exposure_supported": exposure,
            "population_turnover_supported": turnover,
            "null_interpretation_selected": not (
                persistence or exposure or turnover),
            "interpretation_allowed": all(gates.values()),
        },
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "aurorapic_moments_sha256": sha256(args.aurorapic_moments),
            "native_history_sha256": [sha256(
                path / "edupic_phase_eedf_history.csv")
                for path in args.native_directories],
            "native_moments_sha256": [sha256(
                path / "edupic_phase_eedf_moments.csv")
                for path in args.native_directories],
        },
        "aggregation_note": "Per-row history means are weighted by energetic-tail observations. Derived event rates divide aggregate mean event counts by aggregate mean tracked age; all histories start at the bounded continuation origin.",
        "claim_boundary": rule["claim_boundary"],
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_measurement_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
