#!/usr/bin/env python3
"""Analyze the prospective AuroraPIC/native-eduPIC threshold ledger."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


COUNT_COLUMNS = (
    "electron_time_macro_observations",
    "energetic_time_macro_observations",
    "interstep_promotions", "interstep_demotions",
    "elastic_collision_promotions", "elastic_collision_demotions",
    "excitation_collision_promotions", "excitation_collision_demotions",
    "ionization_collision_promotions", "ionization_collision_demotions",
    "charge_exchange_collision_promotions",
    "charge_exchange_collision_demotions",
    "attachment_collision_promotions", "attachment_collision_demotions",
    "bgk_collision_promotions", "bgk_collision_demotions",
    "energetic_births", "subthreshold_births",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def aggregate(rows: list[dict[str, str]], regions: set[str], lower: float,
              upper: float) -> dict[str, float | int]:
    selected = [row for row in rows if row["region"] in regions and
                lower <= float(row["phase_fraction"]) < upper]
    if not selected:
        raise ValueError("threshold-crossing selection is empty")
    result: dict[str, float | int] = {
        "rows": len(selected),
        **{column: sum(int(row[column]) for row in selected)
           for column in COUNT_COLUMNS},
    }
    observations = int(result["electron_time_macro_observations"])
    if observations <= 0:
        raise ValueError("threshold-crossing selection has no electron time")
    result["energetic_fraction"] = (
        int(result["energetic_time_macro_observations"]) / observations)
    for column in COUNT_COLUMNS[2:]:
        result[f"{column}_per_million_electron_steps"] = (
            1.0e6 * int(result[column]) / observations)
    return result


RATE_METRICS = (
    "energetic_fraction",
    "interstep_promotions_per_million_electron_steps",
    "interstep_demotions_per_million_electron_steps",
    "elastic_collision_demotions_per_million_electron_steps",
    "excitation_collision_demotions_per_million_electron_steps",
    "ionization_collision_demotions_per_million_electron_steps",
    "subthreshold_births_per_million_electron_steps",
    "energetic_births_per_million_electron_steps",
)


def mean_members(members: list[dict[str, float | int]]) -> dict[str, float]:
    return {metric: sum(float(member[metric]) for member in members) /
            len(members) for metric in RATE_METRICS}


def ratio(candidate: float, reference: float) -> float | None:
    return candidate / reference if reference != 0.0 else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("aurorapic_crossings", type=Path)
    parser.add_argument("native_directories", nargs=3, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rule = json.loads(args.rule.read_text(encoding="utf-8"))
    diagnostic = rule["diagnostic_contract"]
    regions = set(diagnostic["critical_regions"])
    lower, upper = map(float, diagnostic["critical_phase_fraction"])
    candidate = aggregate(
        read_rows(args.aurorapic_crossings), regions, lower, upper)
    native_paths = [directory /
                    "edupic_phase_eedf_threshold_crossings.csv"
                    for directory in args.native_directories]
    members = [aggregate(read_rows(path), regions, lower, upper)
               for path in native_paths]
    native = mean_members(members)
    ratios = {metric: ratio(float(candidate[metric]), native[metric])
              for metric in RATE_METRICS}
    differences = {metric: float(candidate[metric]) - native[metric]
                   for metric in RATE_METRICS}
    ranges = {metric: ((max(float(member[metric]) for member in members) -
                        min(float(member[metric]) for member in members)) /
                       max(abs(native[metric]), 1e-300))
              for metric in RATE_METRICS}

    promotion_ratio = ratios[
        "interstep_promotions_per_million_electron_steps"]
    demotion_ratio = ratios[
        "interstep_demotions_per_million_electron_steps"]
    excitation_ratio = ratios[
        "excitation_collision_demotions_per_million_electron_steps"]
    ionization_ratio = ratios[
        "ionization_collision_demotions_per_million_electron_steps"]
    birth_ratio = ratios[
        "subthreshold_births_per_million_electron_steps"]
    energetic_birth_population = [int(candidate["energetic_births"]),
                                  *(int(member["energetic_births"])
                                    for member in members)]
    promotion_limited = promotion_ratio is not None and promotion_ratio <= 0.90
    demotion_enhanced = demotion_ratio is not None and demotion_ratio >= 1.10
    collision_loss = any(value is not None and abs(value - 1.0) >= 0.15
                         for value in (excitation_ratio, ionization_ratio))
    birth_supply = (
        birth_ratio is not None and abs(birth_ratio - 1.0) >= 0.15 or
        min(energetic_birth_population) >= 100 and
        ratios["energetic_births_per_million_electron_steps"] is not None and
        abs(float(ratios[
            "energetic_births_per_million_electron_steps"]) - 1.0) >= 0.15)

    primary_repeatability = all(ranges[metric] <= 0.15 for metric in (
        "energetic_fraction",
        "interstep_promotions_per_million_electron_steps",
        "interstep_demotions_per_million_electron_steps",
        "subthreshold_births_per_million_electron_steps"))
    collision_repeatability = all(ranges[metric] <= 0.30 for metric in (
        "elastic_collision_demotions_per_million_electron_steps",
        "excitation_collision_demotions_per_million_electron_steps",
        "ionization_collision_demotions_per_million_electron_steps"))
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
        "electron_time_population": all(
            int(value["electron_time_macro_observations"]) >=
            int(diagnostic[
                "minimum_electron_time_macro_observations_critical_aggregate"])
            for value in population),
        "promotion_population": all(
            int(value["interstep_promotions"]) >=
            int(diagnostic["minimum_interstep_promotions_critical_aggregate"])
            for value in population),
        "demotion_population": all(
            int(value["interstep_demotions"]) >=
            int(diagnostic["minimum_interstep_demotions_critical_aggregate"])
            for value in population),
        "finite_metrics": all(math.isfinite(float(value[metric]))
                              for value in population
                              for metric in RATE_METRICS),
        "native_diagnostic_passivity": all(passivity.values()),
        "native_repeatability_primary": primary_repeatability,
        "native_repeatability_collision_losses": collision_repeatability,
    }
    result = {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "native_edupic_to_aurorapic_unconditional_threshold_crossing_result",
        "rule_sha256": sha256(args.rule),
        "gates": gates,
        "all_measurement_gates_passed": all(gates.values()),
        "native_passivity": passivity,
        "critical_phase_0p125_to_0p5": {
            "aurorapic": candidate,
            "native_edupic_members": members,
            "native_edupic_ensemble_mean": native,
            "native_edupic_relative_range": ranges,
            "aurorapic_to_native_ratio": ratios,
            "aurorapic_minus_native": differences,
        },
        "prospective_decision_outcome": {
            "promotion_limited_supported": promotion_limited,
            "demotion_enhanced_supported": demotion_enhanced,
            "collision_loss_difference_supported": collision_loss,
            "birth_supply_difference_supported": birth_supply,
            "null_interpretation_selected": not (
                promotion_limited or demotion_enhanced or collision_loss or
                birth_supply),
            "interpretation_allowed": all(gates.values()),
        },
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "aurorapic_crossings_sha256": sha256(args.aurorapic_crossings),
            "native_crossings_sha256": [sha256(path)
                                         for path in native_paths],
        },
        "aggregation_note": "Raw macro events and electron-timesteps are summed over the predeclared phase/region rows before rates are formed. Interstep and accepted-collision events are separate views, not additive categories.",
        "claim_boundary": rule["claim_boundary"],
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_measurement_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
