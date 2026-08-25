#!/usr/bin/env python3
"""Analyze the locked two-microstate threshold-crossing replication."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_edupic_threshold_crossings import (
    RATE_METRICS, aggregate, mean_members, ratio, read_rows, sha256,
)


REPEATABILITY_METRICS = (
    "energetic_fraction",
    "interstep_promotions_per_million_electron_steps",
    "interstep_demotions_per_million_electron_steps",
    "excitation_collision_demotions_per_million_electron_steps",
    "ionization_collision_demotions_per_million_electron_steps",
    "subthreshold_births_per_million_electron_steps",
)

PROCESS_COLUMNS = {
    "excitation": "excitation_collision_demotions",
    "ionization": "ionization_collision_demotions",
}


def relative_range(members: list[dict[str, float | int]],
                   metric: str) -> float:
    values = [float(member[metric]) for member in members]
    mean = sum(values) / len(values)
    return (max(values) - min(values)) / max(abs(mean), 1.0e-300)


def analyze(rule_path: Path, candidate_paths: list[Path],
            candidate_report_paths: list[Path],
            native_directories: list[Path]) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    diagnostic = rule["diagnostic_contract"]
    regions = set(diagnostic["critical_regions"])
    lower, upper = map(float, diagnostic["critical_phase_fraction"])
    candidates = [aggregate(read_rows(path), regions, lower, upper)
                  for path in candidate_paths]
    candidate_reports = [json.loads(path.read_text(encoding="utf-8"))
                         for path in candidate_report_paths]
    native_paths = [directory / "edupic_phase_eedf_threshold_crossings.csv"
                    for directory in native_directories]
    natives = [aggregate(read_rows(path), regions, lower, upper)
               for path in native_paths]
    native_mean = mean_members(natives)
    candidate_ranges = {
        metric: relative_range(candidates, metric)
        for metric in REPEATABILITY_METRICS
    }
    member_ratios = [
        {metric: ratio(float(candidate[metric]), native_mean[metric])
         for metric in RATE_METRICS}
        for candidate in candidates
    ]

    locked = rule["locked_inputs"]
    native_hashes = [sha256(path) for path in native_paths]
    checkpoint_hashes = [sha256(directory / "picdata.bin")
                         for directory in native_directories]
    population = [*candidates, *natives]
    process_floor = int(
        diagnostic["minimum_process_demotions_for_relative_comparison"])
    process_qualified = {
        name: all(int(member[column]) >= process_floor
                  for member in population)
        for name, column in PROCESS_COLUMNS.items()
    }
    repeatability_limit = 0.08
    rule_hash = sha256(rule_path)
    locked_states = rule["locked_initial_states"]
    reports_linked = all(
        report.get("all_gates_passed") is True and
        report.get("rule_sha256") == rule_hash and
        report.get("inputs", {}).get("initial_state_id") == locked["id"] and
        report.get("inputs", {}).get("input_checkpoint_sha256") ==
            locked["checkpoint_sha256"] and
        report.get("inputs", {}).get("base_config_sha256") ==
            locked["base_config_sha256"] and
        report.get("inputs", {}).get("prior_report_sha256") ==
            locked["prior_report_sha256"] and
        report.get("inputs", {}).get("solver_sha256") ==
            locked["solver_sha256"] and
        report.get("output_hashes", {}).get(
            "phase_eedf_threshold_crossings.csv") == sha256(path)
        for report, locked, path in zip(
            candidate_reports, locked_states, candidate_paths, strict=True))
    gates = {
        "candidate_count": len(candidates) == 2,
        "native_count": len(natives) == 3,
        "native_crossing_hashes_locked":
            native_hashes == locked["native_crossing_sha256"],
        "native_checkpoint_hashes_locked":
            checkpoint_hashes == locked["native_checkpoint_sha256"],
        "candidate_runner_reports_linked_and_passing": reports_linked,
        "electron_time_population": all(
            int(member["electron_time_macro_observations"]) >= int(
                diagnostic[
                    "minimum_electron_time_macro_observations_critical_aggregate"])
            for member in population),
        "promotion_population": all(
            int(member["interstep_promotions"]) >= int(
                diagnostic["minimum_interstep_promotions_critical_aggregate"])
            for member in population),
        "demotion_population": all(
            int(member["interstep_demotions"]) >= int(
                diagnostic["minimum_interstep_demotions_critical_aggregate"])
            for member in population),
        "finite_metrics": all(
            math.isfinite(float(member[metric]))
            for member in population for metric in RATE_METRICS),
        "candidate_repeatability": all(
            value <= repeatability_limit
            for value in candidate_ranges.values()),
        "qualified_process_populations": all(process_qualified.values()),
    }
    promotion_metric = "interstep_promotions_per_million_electron_steps"
    demotion_metric = "interstep_demotions_per_million_electron_steps"
    promotion_limited = all(
        member[promotion_metric] is not None and
        float(member[promotion_metric]) <= 0.90
        for member in member_ratios)
    demotion_enhanced_rejected = all(
        member[demotion_metric] is not None and
        float(member[demotion_metric]) < 1.10
        for member in member_ratios)
    all_gates = all(gates.values())
    confirmed = all_gates and promotion_limited and demotion_enhanced_rejected
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "aurorapic_threshold_crossing_microstate_replication_result",
        "rule_sha256": rule_hash,
        "gates": gates,
        "all_measurement_and_repeatability_gates_passed": all_gates,
        "critical_phase_0p125_to_0p5": {
            "aurorapic_microstates": candidates,
            "aurorapic_relative_range": candidate_ranges,
            "native_edupic_members": natives,
            "native_edupic_ensemble_mean": native_mean,
            "aurorapic_microstate_to_native_mean_ratios": member_ratios,
            "process_population_qualified": process_qualified,
        },
        "prospective_decision_outcome": {
            "promotion_limited_in_both_microstates": promotion_limited,
            "demotion_enhanced_rejected_in_both_microstates":
                demotion_enhanced_rejected,
            "promotion_limited_mechanism_confirmed": confirmed,
            "interpretation_allowed": all_gates,
        },
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "aurorapic_crossings_sha256": [sha256(path)
                                             for path in candidate_paths],
            "aurorapic_runner_report_sha256": [sha256(path)
                                                 for path in candidate_report_paths],
            "native_crossings_sha256": native_hashes,
            "native_checkpoints_sha256": checkpoint_hashes,
        },
        "aggregation_note": "Counts are summed over the locked region/phase rows before rates are formed. Interstep transitions span the previous accepted collision and following push and are not a field-only decomposition.",
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("candidate_crossings", nargs=2, type=Path)
    parser.add_argument("candidate_reports", nargs=2, type=Path)
    parser.add_argument("native_directories", nargs=3, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = analyze(args.rule, args.candidate_crossings,
                     args.candidate_reports,
                     args.native_directories)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result[
        "all_measurement_and_repeatability_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
