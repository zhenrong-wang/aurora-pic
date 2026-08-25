#!/usr/bin/env python3
"""Analyze the prospective electrode half-cell Gauss control."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_edupic_field_push_thresholds import (
    RATE_METRICS, aggregate, relative_range,
)
from analyze_edupic_threshold_crossings import sha256


def decision(ratios: list[float], gates_pass: bool) -> dict[str, bool]:
    persistence = all(value <= 0.90 for value in ratios)
    parity = all(0.95 <= value <= 1.05 for value in ratios)
    return {
        "deficit_persists": gates_pass and persistence,
        "endpoint_correction_explains_deficit": gates_pass and parity,
        "partial_endpoint_effect": gates_pass and not (persistence or parity),
        "interpretation_allowed": gates_pass,
        "directional_persistence_signal_without_interpretation":
            persistence and not gates_pass,
    }


def analyze(rule_path: Path, crossing_paths: list[Path],
            report_paths: list[Path]) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    diagnostic = rule["diagnostic_contract"]
    regions = set(diagnostic["critical_regions"])
    lower, upper = map(float, diagnostic["critical_phase_fraction"])
    members = [aggregate(path, regions, lower, upper)
               for path in crossing_paths]
    native = rule["locked_inputs"]["native_field_push_rates_per_million"]
    baseline = rule["locked_inputs"][
        "baseline_pooled_candidate_rates_per_million"]
    promotion = "field_push_promotions_per_million_pushes"
    demotion = "field_push_demotions_per_million_pushes"
    ratios = [{
        promotion: float(member[promotion]) / float(native["promotions"]),
        demotion: float(member[demotion]) / float(native["demotions"]),
    } for member in members]
    baseline_ratios = [{
        promotion: float(member[promotion]) /
            float(baseline[state["id"]]["promotions"]),
        demotion: float(member[demotion]) /
            float(baseline[state["id"]]["demotions"]),
    } for member, state in zip(
        members, rule["locked_continuation_states"], strict=True)]
    ranges = {metric: relative_range(members, metric)
              for metric in RATE_METRICS}
    rule_hash = sha256(rule_path)
    reports = [json.loads(path.read_text(encoding="utf-8"))
               for path in report_paths]
    linked = all(
        report.get("all_gates_passed") is True and
        report.get("rule_sha256") == rule_hash and
        report.get("inputs", {}).get("initial_state_id") == state["id"] and
        report.get("inputs", {}).get("solver_sha256") ==
            state["solver_sha256"] and
        report.get("inputs", {}).get("input_checkpoint_sha256") ==
            state["checkpoint_sha256"] and
        report.get("inputs", {}).get("prior_report_sha256") ==
            state["prior_report_sha256"] and
        report.get("output_hashes", {}).get(
            "phase_eedf_threshold_crossings.csv") == sha256(path) and
        report.get("algorithm_contract", {}).get(
            "collision_velocity_sampling") == "leapfrog_half_step" and
        report.get("algorithm_contract", {}).get(
            "random_stream_continued_from_checkpoint") is True and
        report.get("algorithm_contract", {}).get(
            "electrode_half_cell_gauss_correction") is True
        for report, state, path in zip(
            reports, rule["locked_continuation_states"], crossing_paths,
            strict=True))
    gates = {
        "runner_reports_linked_and_passing": linked,
        "finite_metrics": all(
            math.isfinite(float(member[metric]))
            for member in members for metric in RATE_METRICS),
        "field_push_observation_population": all(
            int(member["field_push_macro_observations"]) >= int(
                diagnostic[
                    "minimum_field_push_observations_critical_aggregate"])
            for member in members),
        "field_push_promotion_population": all(
            int(member["field_push_promotions"]) >= int(
                diagnostic[
                    "minimum_field_push_promotions_critical_aggregate"])
            for member in members),
        "field_push_demotion_population": all(
            int(member["field_push_demotions"]) >= int(
                diagnostic[
                    "minimum_field_push_demotions_critical_aggregate"])
            for member in members),
        "candidate_repeatability": all(
            ranges[metric] <= 0.08 for metric in RATE_METRICS),
    }
    all_gates = all(gates.values())
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "electrode_half_cell_gauss_control_result",
        "rule_sha256": rule_hash,
        "gates": gates,
        "all_runner_population_provenance_and_repeatability_gates_passed":
            all_gates,
        "critical_phase_0p125_to_0p5": {
            "endpoint_corrected_aurorapic_microstates": members,
            "native_edupic_ensemble_mean": native,
            "endpoint_corrected_aurorapic_to_native_ratios": ratios,
            "endpoint_corrected_aurorapic_to_pooled_baseline_ratios":
                baseline_ratios,
            "endpoint_corrected_aurorapic_relative_range": ranges,
        },
        "prospective_decision_outcome": decision(
            [float(item[promotion]) for item in ratios], all_gates),
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "crossings_sha256": [sha256(path) for path in crossing_paths],
            "runner_report_sha256": [sha256(path) for path in report_paths],
        },
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("crossings", nargs=2, type=Path)
    parser.add_argument("reports", nargs=2, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = analyze(args.rule, args.crossings, args.reports)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result[
        "all_runner_population_provenance_and_repeatability_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
