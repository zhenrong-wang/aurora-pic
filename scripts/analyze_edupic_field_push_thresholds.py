#!/usr/bin/env python3
"""Analyze the locked AuroraPIC/native eduPIC field-push campaign."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from analyze_edupic_threshold_crossings import sha256


RATE_METRICS = (
    "field_push_promotions_per_million_pushes",
    "field_push_demotions_per_million_pushes",
)


def aggregate(path: Path, regions: set[str], lower: float,
              upper: float) -> dict[str, float | int]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = [row for row in csv.DictReader(stream)
                if row["region"] in regions and
                lower <= float(row["phase_fraction"]) < upper]
    if not rows:
        raise ValueError("field-push threshold selection is empty")
    observations = sum(int(row["field_push_macro_observations"])
                       for row in rows)
    promotions = sum(int(row["field_push_promotions"]) for row in rows)
    demotions = sum(int(row["field_push_demotions"]) for row in rows)
    if observations <= 0 or promotions > observations or demotions > observations:
        raise ValueError("field-push threshold counts are invalid")
    return {
        "rows": len(rows),
        "field_push_macro_observations": observations,
        "field_push_promotions": promotions,
        "field_push_demotions": demotions,
        "field_push_promotions_per_million_pushes":
            1.0e6 * promotions / observations,
        "field_push_demotions_per_million_pushes":
            1.0e6 * demotions / observations,
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
    candidates = [aggregate(path, regions, lower, upper)
                  for path in candidate_paths]
    native_paths = [directory /
                    "edupic_phase_eedf_field_push_thresholds.csv"
                    for directory in native_directories]
    natives = [aggregate(path, regions, lower, upper)
               for path in native_paths]
    native_mean = {
        metric: sum(float(member[metric]) for member in natives) / len(natives)
        for metric in RATE_METRICS
    }
    ratios = [
        {metric: float(member[metric]) / native_mean[metric]
         if native_mean[metric] else None
         for metric in RATE_METRICS}
        for member in candidates
    ]
    candidate_ranges = {
        metric: relative_range(candidates, metric) for metric in RATE_METRICS}
    native_ranges = {
        metric: relative_range(natives, metric) for metric in RATE_METRICS}

    rule_hash = sha256(rule_path)
    reports = [json.loads(path.read_text(encoding="utf-8"))
               for path in candidate_report_paths]
    locked_states = rule["locked_initial_states"]
    reports_linked = all(
        report.get("all_gates_passed") is True and
        report.get("rule_sha256") == rule_hash and
        report.get("inputs", {}).get("initial_state_id") == locked["id"] and
        report.get("inputs", {}).get("solver_sha256") ==
            locked["solver_sha256"] and
        report.get("inputs", {}).get("input_checkpoint_sha256") ==
            locked["checkpoint_sha256"] and
        report.get("output_hashes", {}).get(
            "phase_eedf_threshold_crossings.csv") == sha256(path)
        for report, locked, path in zip(
            reports, locked_states, candidate_paths, strict=True))
    expected_checkpoints = rule["locked_inputs"][
        "expected_passive_native_checkpoint_sha256"]
    native_passivity = {
        directory.name.removeprefix("seed-"):
            sha256(directory / "picdata.bin") == expected_checkpoints[
                directory.name.removeprefix("seed-")]
        for directory in native_directories
    }
    population = [*candidates, *natives]
    gates = {
        "candidate_runner_reports_linked_and_passing": reports_linked,
        "native_diagnostic_passivity": all(native_passivity.values()),
        "finite_metrics": all(
            math.isfinite(float(member[metric]))
            for member in population for metric in RATE_METRICS),
        "field_push_observation_population": all(
            int(member["field_push_macro_observations"]) >= int(
                diagnostic[
                    "minimum_field_push_observations_critical_aggregate"])
            for member in population),
        "field_push_promotion_population": all(
            int(member["field_push_promotions"]) >= int(
                diagnostic[
                    "minimum_field_push_promotions_critical_aggregate"])
            for member in population),
        "field_push_demotion_population": all(
            int(member["field_push_demotions"]) >= int(
                diagnostic["minimum_field_push_demotions_critical_aggregate"])
            for member in population),
        "candidate_repeatability": all(
            candidate_ranges[metric] <= 0.08 for metric in RATE_METRICS),
        "native_repeatability": all(
            native_ranges[metric] <= 0.15 for metric in RATE_METRICS),
    }
    all_gates = all(gates.values())
    promotion = "field_push_promotions_per_million_pushes"
    demotion = "field_push_demotions_per_million_pushes"
    promotion_deficit = all(
        member[promotion] is not None and float(member[promotion]) <= 0.90
        for member in ratios)
    promotion_parity = all(
        member[promotion] is not None and
        0.95 <= float(member[promotion]) <= 1.05 for member in ratios)
    demotion_enhanced = all(
        member[demotion] is not None and float(member[demotion]) >= 1.10
        for member in ratios)
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "native_edupic_to_aurorapic_field_push_threshold_result",
        "rule_sha256": rule_hash,
        "gates": gates,
        "all_measurement_repeatability_and_passivity_gates_passed": all_gates,
        "native_passivity": native_passivity,
        "critical_phase_0p125_to_0p5": {
            "aurorapic_microstates": candidates,
            "native_edupic_members": natives,
            "native_edupic_ensemble_mean": native_mean,
            "aurorapic_to_native_ratios": ratios,
            "aurorapic_relative_range": candidate_ranges,
            "native_relative_range": native_ranges,
        },
        "prospective_decision_outcome": {
            "field_push_promotion_deficit_supported":
                all_gates and promotion_deficit,
            "field_push_promotion_parity_supported":
                all_gates and promotion_parity,
            "field_push_demotion_enhancement_supported":
                all_gates and demotion_enhanced,
            "mixed_or_intermediate_result": all_gates and not (
                promotion_deficit or promotion_parity or demotion_enhanced),
            "interpretation_allowed": all_gates,
        },
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "aurorapic_crossings_sha256": [sha256(path)
                                             for path in candidate_paths],
            "aurorapic_runner_report_sha256": [sha256(path)
                                                 for path in candidate_report_paths],
            "native_field_push_sha256": [sha256(path)
                                           for path in native_paths],
            "native_checkpoints_sha256": [sha256(directory / "picdata.bin")
                                             for directory in native_directories],
        },
        "algorithm_contract_limit": rule["algorithm_contract_limit"],
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
                     args.candidate_reports, args.native_directories)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result[
        "all_measurement_repeatability_and_passivity_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
