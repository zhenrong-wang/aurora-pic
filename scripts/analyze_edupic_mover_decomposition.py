#!/usr/bin/env python3
"""Evaluate the locked AuroraPIC/eduPIC mover-decomposition comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


METRICS = (
    "mean_origin_energy_eV",
    "origin_longitudinal_energy_fraction",
    "mean_positive_linear_work_eV",
    "mean_quadratic_work_eV",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def means(members: list[dict[str, object]]) -> dict[str, float]:
    return {
        metric: sum(float(member[metric]) for member in members) / len(members)
        for metric in METRICS
    }


def relative_ranges(members: list[dict[str, object]]) -> dict[str, float]:
    result = {}
    for metric in METRICS:
        values = [float(member[metric]) for member in members]
        mean = sum(values) / len(values)
        result[metric] = (max(values) - min(values)) / max(abs(mean), 1e-300)
    return result


def classify(value: float, parity: tuple[float, float]) -> str:
    if value < parity[0]:
        return "aurorapic_lower"
    if value > parity[1]:
        return "aurorapic_higher"
    return "parity"


def analyze(rule_path: Path, native_path: Path,
            report_paths: list[Path]) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    native = json.loads(native_path.read_text(encoding="utf-8"))
    reports = [json.loads(path.read_text(encoding="utf-8"))
               for path in report_paths]
    rule_hash = sha256(rule_path)
    locked = rule["locked_inputs"]
    diagnostic = rule["diagnostic_contract"]
    states = locked["aurorapic_continuation_states"]
    aurora = [report["critical_scope"] for report in reports]
    native_members = native["members"]
    aurora_means = means(aurora)
    native_means = means(native_members)
    aurora_ranges = relative_ranges(aurora)
    native_ranges = relative_ranges(native_members)
    ratios = [
        {metric: float(member[metric]) / native_means[metric]
         for metric in METRICS}
        for member in aurora
    ]
    mean_ratios = {
        metric: aurora_means[metric] / native_means[metric]
        for metric in METRICS
    }
    parity = tuple(map(float,
                       rule["prospective_decision_rule"]["parity_interval"]))
    reports_linked = len(reports) == len(states) == 2 and all(
        report.get("all_gates_passed") is True
        and report.get("rule_sha256") == rule_hash
        and report.get("state_id") == state["id"]
        and report.get("inputs", {}).get("solver_sha256") ==
            locked["aurorapic_solver_sha256"]
        and report.get("inputs", {}).get("checkpoint_sha256") ==
            state["checkpoint_sha256"]
        and report.get("inputs", {}).get("prior_report_sha256") ==
            state["prior_report_sha256"]
        and report.get("inputs", {}).get("prior_deck_sha256") ==
            state["prior_deck_sha256"]
        for report, state in zip(reports, states, strict=True)
    )
    native_linked = (
        native.get("rule_sha256") == rule_hash
        and native.get(
            "all_native_integrity_population_repeatability_and_closure_gates_passed"
        ) is True
        and [member.get("seed") for member in native_members] ==
            locked["native_seeds"]
    )
    all_members = [*aurora, *native_members]
    tolerance = float(diagnostic["work_closure_relative_tolerance"])
    gates = {
        "aurorapic_member_count": len(aurora) == 2,
        "native_member_count": len(native_members) == 3,
        "aurorapic_reports_linked_and_passing": reports_linked,
        "native_result_linked_and_passing": native_linked,
        "finite_metrics": all(
            math.isfinite(float(member[metric]))
            for member in all_members for metric in METRICS),
        "population_floors": all(
            int(member["promotion_band_observations"]
                if "promotion_band_observations" in member
                else member["field_push_promotion_band_observations"]) >=
                int(diagnostic["minimum_band_observations_per_member"])
            and int(member["promotion_band_promotions"]
                    if "promotion_band_promotions" in member
                    else member["field_push_promotion_band_promotions"]) >=
                int(diagnostic["minimum_band_promotions_per_member"])
            for member in all_members),
        "origin_energy_partition": all(
            0.0 <= float(member[
                "origin_longitudinal_energy_fraction"]) <= 1.0
            for member in all_members),
        "aurorapic_closures": all(
            abs(float(member["linear_work_closure_residual_eV"])) <=
                tolerance * max(
                    1.0,
                    abs(float(member["positive_linear_work_sum_eV"])) +
                    abs(float(member["negative_linear_work_sum_eV"])))
            and abs(float(member[
                "total_work_decomposition_closure_residual_eV"])) <=
                tolerance * max(
                    1.0, abs(float(member["linear_work_sum_eV"])) +
                    abs(float(member["quadratic_work_sum_eV"])))
            for member in aurora),
        "aurorapic_repeatability": all(
            value <= 0.10 for value in aurora_ranges.values()),
        "native_repeatability": all(
            value <= 0.15 for value in native_ranges.values()),
    }
    all_gates = all(gates.values())
    field_strength = all(
        ratio["mean_quadratic_work_eV"] <= 0.90 for ratio in ratios)
    alignment = all(
        ratio["mean_positive_linear_work_eV"] <= 0.90 for ratio in ratios)
    longitudinal = all(
        ratio["origin_longitudinal_energy_fraction"] <= 0.90
        for ratio in ratios)
    origin = all(
        ratio["mean_origin_energy_eV"] <= 0.995 for ratio in ratios)
    supported = [
        name for name, value in (
            ("field_strength_deficit", field_strength),
            ("favorable_alignment_deficit", alignment),
            ("longitudinal_energy_partition_deficit", longitudinal),
            ("lower_origin_energy_distribution", origin),
        ) if all_gates and value
    ]
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "aurorapic_edupic_near_threshold_mover_decomposition_result",
        "rule_sha256": rule_hash,
        "gates": gates,
        "all_integrity_population_repeatability_and_closure_gates_passed":
            all_gates,
        "aurorapic_members": aurora,
        "aurorapic_ensemble_means": aurora_means,
        "aurorapic_relative_ranges": aurora_ranges,
        "native_edupic_ensemble_means": native_means,
        "native_edupic_relative_ranges": native_ranges,
        "aurorapic_member_to_native_mean_ratios": ratios,
        "aurorapic_ensemble_mean_to_native_mean_ratios": mean_ratios,
        "ensemble_mean_ratio_classification": {
            metric: classify(value, parity)
            for metric, value in mean_ratios.items()
        },
        "prospective_decision_outcome": {
            "interpretation_allowed": all_gates,
            "field_strength_deficit_supported": all_gates and field_strength,
            "favorable_alignment_deficit_supported": all_gates and alignment,
            "longitudinal_energy_partition_deficit_supported":
                all_gates and longitudinal,
            "lower_origin_energy_distribution_supported":
                all_gates and origin,
            "supported_associations": supported,
            "result": (
                "field_strength_and_favorable_alignment_deficits_supported_"
                "origin_energy_and_longitudinal_partition_deficits_not_supported"
                if all_gates and field_strength and alignment and
                   not longitudinal and not origin
                else "mixed_or_intermediate_near_threshold_mover_result"
            ),
        },
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "native_result_sha256": sha256(native_path),
            "aurorapic_runner_report_sha256": [
                sha256(path) for path in report_paths
            ],
        },
        "interpretation_note": (
            "The exact decomposition supports simultaneous sampled-field and "
            "favorable-alignment deficits in the locked scope. These are "
            "associations within one evolved discharge, not independent causal "
            "effects or identification of the first field divergence."
        ),
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("native_result", type=Path)
    parser.add_argument("aurorapic_reports", nargs=2, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = analyze(args.rule, args.native_result, args.aurorapic_reports)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result[
        "all_integrity_population_repeatability_and_closure_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
