#!/usr/bin/env python3
"""Analyze prospectively pooled matched-half-step threshold blocks."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_edupic_field_push_thresholds import (
    RATE_METRICS, aggregate, relative_range,
)
from analyze_edupic_threshold_crossings import sha256


def pool(blocks: list[dict[str, float | int]]) -> dict[str, float | int]:
    observations = sum(int(item["field_push_macro_observations"])
                       for item in blocks)
    promotions = sum(int(item["field_push_promotions"]) for item in blocks)
    demotions = sum(int(item["field_push_demotions"]) for item in blocks)
    if observations <= 0:
        raise ValueError("pooled field-push observation count is empty")
    return {
        "field_push_macro_observations": observations,
        "field_push_promotions": promotions,
        "field_push_demotions": demotions,
        "field_push_promotions_per_million_pushes":
            1.0e6 * promotions / observations,
        "field_push_demotions_per_million_pushes":
            1.0e6 * demotions / observations,
    }


def analyze(rule_path: Path, first_result_path: Path,
            first_crossings: list[Path], second_crossings: list[Path],
            first_reports: list[Path], second_reports: list[Path]
            ) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    diagnostic = rule["diagnostic_contract"]
    regions = set(diagnostic["critical_regions"])
    lower, upper = map(float, diagnostic["critical_phase_fraction"])
    first = [aggregate(path, regions, lower, upper)
             for path in first_crossings]
    second = [aggregate(path, regions, lower, upper)
              for path in second_crossings]
    pooled = [pool([a, b]) for a, b in zip(first, second, strict=True)]
    native = rule["locked_inputs"]["native_field_push_rates_per_million"]
    promotion = "field_push_promotions_per_million_pushes"
    demotion = "field_push_demotions_per_million_pushes"
    ratios = [{
        promotion: float(item[promotion]) / float(native["promotions"]),
        demotion: float(item[demotion]) / float(native["demotions"]),
    } for item in pooled]
    ranges = {metric: relative_range(pooled, metric)
              for metric in RATE_METRICS}
    rule_hash = sha256(rule_path)
    locked = rule["locked_continuation_states"]
    first_report_values = [json.loads(path.read_text(encoding="utf-8"))
                           for path in first_reports]
    second_report_values = [json.loads(path.read_text(encoding="utf-8"))
                            for path in second_reports]
    first_linked = all(
        sha256(crossing) == state["first_block_crossings_sha256"] and
        sha256(report_path) == state["prior_report_sha256"] and
        report.get("all_gates_passed") is True and
        report.get("inputs", {}).get("initial_state_id") == state["id"]
        for state, crossing, report_path, report in zip(
            locked, first_crossings, first_reports, first_report_values,
            strict=True))
    second_linked = all(
        report.get("all_gates_passed") is True and
        report.get("rule_sha256") == rule_hash and
        report.get("inputs", {}).get("initial_state_id") == state["id"] and
        report.get("inputs", {}).get("input_checkpoint_sha256") ==
            state["checkpoint_sha256"] and
        report.get("inputs", {}).get("prior_report_sha256") ==
            state["prior_report_sha256"] and
        report.get("output_hashes", {}).get(
            "phase_eedf_threshold_crossings.csv") == sha256(crossing) and
        report.get("algorithm_contract", {}).get(
            "collision_velocity_sampling") == "leapfrog_half_step" and
        report.get("algorithm_contract", {}).get(
            "random_stream_continued_from_checkpoint") is True
        for state, crossing, report in zip(
            locked, second_crossings, second_report_values, strict=True))
    first_result = json.loads(first_result_path.read_text(encoding="utf-8"))
    first_result_linked = (
        sha256(first_result_path) ==
            rule["locked_inputs"]["first_block_result_sha256"] and
        first_result.get("prospective_decision_outcome", {}).get(
            "interpretation_allowed") is False and
        first_result.get("gates", {}).get("candidate_repeatability") is False)
    population = [*first, *second]
    gates = {
        "first_block_result_linked": first_result_linked,
        "first_block_reports_and_crossings_linked": first_linked,
        "second_block_reports_and_crossings_linked": second_linked,
        "finite_metrics": all(
            math.isfinite(float(item[metric]))
            for item in pooled for metric in RATE_METRICS),
        "each_block_observation_population": all(
            int(item["field_push_macro_observations"]) >= int(
                diagnostic[
                    "minimum_field_push_observations_critical_aggregate"])
            for item in population),
        "each_block_promotion_population": all(
            int(item["field_push_promotions"]) >= int(
                diagnostic[
                    "minimum_field_push_promotions_critical_aggregate"])
            for item in population),
        "each_block_demotion_population": all(
            int(item["field_push_demotions"]) >= int(
                diagnostic[
                    "minimum_field_push_demotions_critical_aggregate"])
            for item in population),
        "pooled_candidate_repeatability": all(
            ranges[metric] <= 0.08 for metric in RATE_METRICS),
    }
    all_gates = all(gates.values())
    persistence = all(float(item[promotion]) <= 0.90 for item in ratios)
    parity = all(0.95 <= float(item[promotion]) <= 1.05 for item in ratios)
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "pooled_matched_half_step_threshold_replication_result",
        "rule_sha256": rule_hash,
        "gates": gates,
        "all_provenance_population_and_repeatability_gates_passed": all_gates,
        "critical_phase_0p125_to_0p5": {
            "first_blocks": first,
            "second_blocks": second,
            "pooled_microstates": pooled,
            "native_edupic_ensemble_mean": native,
            "pooled_aurorapic_to_native_ratios": ratios,
            "pooled_aurorapic_relative_range": ranges,
        },
        "prospective_decision_outcome": {
            "deficit_persists": all_gates and persistence,
            "staggering_explains_deficit": all_gates and parity,
            "partial_staggering_effect": all_gates and not (
                persistence or parity),
            "interpretation_allowed": all_gates,
        },
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "first_result_sha256": sha256(first_result_path),
            "first_crossings_sha256": [sha256(path)
                                           for path in first_crossings],
            "second_crossings_sha256": [sha256(path)
                                            for path in second_crossings],
            "first_report_sha256": [sha256(path) for path in first_reports],
            "second_report_sha256": [sha256(path) for path in second_reports],
        },
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("first_result", type=Path)
    parser.add_argument("first_crossings", nargs=2, type=Path)
    parser.add_argument("second_crossings", nargs=2, type=Path)
    parser.add_argument("first_reports", nargs=2, type=Path)
    parser.add_argument("second_reports", nargs=2, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = analyze(
        args.rule, args.first_result, args.first_crossings,
        args.second_crossings, args.first_reports, args.second_reports)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result[
        "all_provenance_population_and_repeatability_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
