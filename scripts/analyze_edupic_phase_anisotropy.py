#!/usr/bin/env python3
"""Analyze the prospective AuroraPIC/native-eduPIC anisotropy campaign."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


METRICS = (
    "tail_fraction",
    "temperature_x_fraction",
    "temperature_y_fraction",
    "temperature_z_fraction",
    "tail_longitudinal_energy_fraction",
    "tail_directional_population_imbalance",
    "tail_mean_velocity_x",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def aggregate(source: list[dict[str, str]], regions: set[str],
              lower: float, upper: float) -> dict[str, float | int]:
    selected = [row for row in source if row["region"] in regions and
                lower <= float(row["phase_fraction"]) < upper]
    if not selected:
        raise ValueError("anisotropy selection is empty")
    represented_key = ("represented_observations" if
                       "represented_observations" in selected[0] else
                       "represented_observations_m-2")
    tail_key = ("tail_represented_observations" if
                "tail_represented_observations" in selected[0] else
                "tail_represented_observations_m-2")
    represented = [float(row[represented_key]) for row in selected]
    tail = [float(row[tail_key]) for row in selected]
    macro = [int(row["macro_observations"]) for row in selected]
    weights = [represented[index] / macro[index] if macro[index] else 0.0
               for index in range(len(selected))]
    tail_macro = sum(tail[index] / weights[index]
                     for index in range(len(selected)) if weights[index] > 0.0)

    def weighted(column: str, values: list[float]) -> float:
        denominator = sum(values)
        return (sum(float(row[column]) * weight
                    for row, weight in zip(selected, values)) / denominator
                if denominator else 0.0)

    temperatures = [weighted(column, represented) for column in
                    ("temperature_x", "temperature_y", "temperature_z")]
    temperature_sum = sum(temperatures)
    return {
        "rows": len(selected),
        "minimum_macro_observations_per_region_phase_bin": min(macro),
        "represented_observations": sum(represented),
        "tail_represented_observations": sum(tail),
        "tail_macro_observations": int(round(tail_macro)),
        "tail_fraction": sum(tail) / sum(represented),
        "temperature_x": temperatures[0],
        "temperature_y": temperatures[1],
        "temperature_z": temperatures[2],
        "temperature_x_fraction": temperatures[0] / temperature_sum,
        "temperature_y_fraction": temperatures[1] / temperature_sum,
        "temperature_z_fraction": temperatures[2] / temperature_sum,
        # Version 1 emits a per-phase/region ratio but not its energy
        # denominator. Tail-count weighting is therefore the prospectively
        # retained aggregate estimator; individual octants are also reported.
        "tail_longitudinal_energy_fraction": weighted(
            "tail_longitudinal_energy_fraction", tail),
        "tail_directional_population_imbalance": weighted(
            "tail_directional_population_imbalance", tail),
        "tail_mean_velocity_x": weighted("tail_mean_velocity_x", tail),
    }


def mean_members(members: list[dict[str, float | int]]) -> dict[str, float]:
    return {metric: sum(float(member[metric]) for member in members) /
            len(members) for metric in METRICS}


def comparison(candidate: dict[str, float | int],
               members: list[dict[str, float | int]]) -> dict[str, object]:
    native = mean_members(members)
    ranges = {}
    for metric in METRICS:
        values = [float(member[metric]) for member in members]
        ranges[metric] = ((max(values) - min(values)) /
                          max(abs(native[metric]), 1e-300))
    return {
        "aurorapic": candidate,
        "native_edupic_members": members,
        "native_edupic_ensemble_mean": native,
        "native_edupic_relative_range": ranges,
        "aurorapic_to_native_ratio": {
            metric: float(candidate[metric]) / native[metric]
            if native[metric] != 0.0 else None for metric in METRICS},
        "aurorapic_minus_native": {
            metric: float(candidate[metric]) - native[metric]
            for metric in METRICS},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("aurorapic_moments", type=Path)
    parser.add_argument("native_directories", nargs=3, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rule = json.loads(args.rule.read_text(encoding="utf-8"))
    diagnostic = rule["diagnostic_contract"]
    critical = set(diagnostic["critical_regions"])
    lower, upper = map(float, diagnostic["critical_phase_fraction"])
    aurora_rows = rows(args.aurorapic_moments)
    native_paths = [directory / "edupic_phase_eedf_moments.csv"
                    for directory in args.native_directories]
    native_rows = [rows(path) for path in native_paths]
    critical_comparison = comparison(
        aggregate(aurora_rows, critical, lower, upper),
        [aggregate(member, critical, lower, upper)
         for member in native_rows])
    octants = []
    phase = lower
    while phase < upper - 1e-12:
        next_phase = min(upper, phase + 0.125)
        value = comparison(
            aggregate(aurora_rows, critical, phase, next_phase),
            [aggregate(member, critical, phase, next_phase)
             for member in native_rows])
        value["phase_fraction"] = [phase, next_phase]
        octants.append(value)
        phase = next_phase

    critical_diff = critical_comparison["aurorapic_minus_native"]
    longitudinal = abs(float(
        critical_diff["tail_longitudinal_energy_fraction"])) >= 0.05
    component = any(abs(float(critical_diff[key])) >= 0.05 for key in (
        "temperature_x_fraction", "temperature_y_fraction",
        "temperature_z_fraction"))
    directional = any(abs(float(octant["aurorapic_minus_native"]
                                    ["tail_directional_population_imbalance"]))
                      >= 0.05 for octant in octants)
    minimum_population = int(
        diagnostic["minimum_macro_observations_per_critical_region_phase_bin"])
    minimum_tail = int(
        diagnostic["minimum_tail_macro_observations_critical_aggregate"])
    all_aggregates = [critical_comparison["aurorapic"],
                      *critical_comparison["native_edupic_members"]]
    passivity_expected = {
        "13507": "7b88021958430ed2157cc0dba4b3d20fd23ec03593c1e3f29b4c66ba540b13cd",
        "24601": "8935e4313ddef5c171a55f4631c9aab2665a772e799a20edaa297ee9ec78c2fd",
        "35713": "92bd734fc3710139ae9c212bce23f2f7559bb7e5a17599dca21e0b9a4eacf267",
    }
    passivity = {}
    for directory in args.native_directories:
        seed = directory.name.removeprefix("seed-")
        passivity[seed] = sha256(directory / "picdata.bin") == \
            passivity_expected[seed]
    gates = {
        "critical_observation_population": all(
            int(value["minimum_macro_observations_per_region_phase_bin"])
            >= minimum_population for value in all_aggregates),
        "critical_tail_population": all(
            int(value["tail_macro_observations"]) >= minimum_tail
            for value in all_aggregates),
        "native_diagnostic_passivity": all(passivity.values()),
        "finite_metrics": all(math.isfinite(float(value[metric]))
                              for value in all_aggregates
                              for metric in METRICS),
    }
    result = {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "native_edupic_to_aurorapic_velocity_anisotropy_result",
        "rule_sha256": sha256(args.rule),
        "gates": gates,
        "all_measurement_gates_passed": all(gates.values()),
        "native_passivity": passivity,
        "critical_phase_0p125_to_0p5": critical_comparison,
        "critical_phase_octants": octants,
        "prospective_decision_outcome": {
            "longitudinal_mechanism_supported": longitudinal,
            "component_redistribution_supported": component,
            "phase_directionality_supported": directional,
            "null_interpretation_selected": not (
                longitudinal or component or directional),
        },
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "aurorapic_moments_sha256": sha256(args.aurorapic_moments),
            "native_moments_sha256": [sha256(path) for path in native_paths],
        },
        "aggregation_note": "Tail longitudinal energy fraction is tail-count-weighted across phase/region rows because diagnostic version 1 emits row ratios without their kinetic-energy denominators.",
        "claim_boundary": rule["claim_boundary"],
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_measurement_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
