#!/usr/bin/env python3
"""Evaluate the locked AuroraPIC/EDIPIC promotion-band work comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


METRICS = (
    "band_supply_fraction",
    "band_promotion_probability",
    "mean_signed_work_eV",
    "mean_positive_work_eV",
    "mean_negative_work_eV",
)
REPEATABILITY_METRICS = (
    "band_supply_fraction",
    "band_promotion_probability",
    "mean_positive_work_eV",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def means(members: list[dict[str, float]]) -> dict[str, float]:
    return {
        metric: sum(float(member[metric]) for member in members) / len(members)
        for metric in METRICS
    }


def relative_ranges(members: list[dict[str, float]]) -> dict[str, float]:
    result = {}
    for metric in METRICS:
        values = [float(member[metric]) for member in members]
        mean = sum(values) / len(values)
        result[metric] = (max(values) - min(values)) / max(abs(mean), 1.0e-300)
    return result


def classify_ratio(value: float, parity: tuple[float, float]) -> str:
    if value < parity[0]:
        return "aurorapic_lower"
    if value > parity[1]:
        return "aurorapic_higher"
    return "parity"


def analyze(rule_path: Path, native_path: Path,
            aurora_report_paths: list[Path]) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    native = json.loads(native_path.read_text(encoding="utf-8"))
    reports = [json.loads(path.read_text(encoding="utf-8"))
               for path in aurora_report_paths]
    rule_hash = sha256(rule_path)
    locked = rule["locked_inputs"]
    diagnostic = rule["diagnostic_contract"]
    decision = rule["prospective_decision_rule"]
    locked_states = locked["aurorapic_continuation_states"]

    aurora_members = [dict(report["critical_scope"]) for report in reports]
    native_members = [dict(member) for member in native["members"]]
    aurora_means = means(aurora_members)
    native_means = means(native_members)
    aurora_ranges = relative_ranges(aurora_members)
    native_ranges = relative_ranges(native_members)
    ratios = [
        {metric: float(member[metric]) / native_means[metric]
         for metric in METRICS}
        for member in aurora_members
    ]
    mean_ratios = {
        metric: aurora_means[metric] / native_means[metric]
        for metric in METRICS
    }
    parity = tuple(map(float, decision["parity_interval"]))

    reports_linked = len(reports) == len(locked_states) == 2 and all(
        report.get("all_gates_passed") is True
        and report.get("rule_sha256") == rule_hash
        and report.get("state_id") == state["id"]
        and report.get("inputs", {}).get("solver_sha256")
            == locked["aurorapic_solver_sha256"]
        and report.get("inputs", {}).get("checkpoint_sha256")
            == state["checkpoint_sha256"]
        and report.get("inputs", {}).get("prior_report_sha256")
            == state["prior_report_sha256"]
        and report.get("inputs", {}).get("prior_deck_sha256")
            == state["prior_deck_sha256"]
        for report, state in zip(reports, locked_states, strict=True)
    )
    native_linked = (
        native.get("rule_sha256") == rule_hash
        and native.get(
            "all_native_integrity_population_repeatability_and_work_closure_gates_passed"
        ) is True
        and [member.get("seed") for member in native_members]
            == locked["native_seeds"]
        and all(
            member.get("final_checkpoint_sha256")
                == locked["expected_passive_native_checkpoint_sha256"][str(member["seed"])]
            for member in native_members
        )
    )
    all_members = [*aurora_members, *native_members]
    gates = {
        "aurorapic_member_count": len(aurora_members) == 2,
        "native_member_count": len(native_members) == 3,
        "aurorapic_runner_reports_linked_and_passing": reports_linked,
        "native_result_linked_and_passing": native_linked,
        "finite_metrics": all(
            math.isfinite(float(member[metric]))
            for member in all_members for metric in METRICS
        ),
        "band_observation_population": all(
            int(member["promotion_band_observations"])
                >= int(diagnostic["minimum_band_observations_per_member"])
            for member in all_members
        ),
        "band_promotion_population": all(
            int(member["promotion_band_promotions"])
                >= int(diagnostic["minimum_band_promotions_per_member"])
            for member in all_members
        ),
        "work_closure": all(
            abs(float(member["work_closure_residual_eV"]))
                <= float(diagnostic["work_closure_relative_tolerance"])
                * max(abs(float(member["positive_macro_work_sum_eV"]))
                      if "positive_macro_work_sum_eV" in member else 1.0,
                      1.0)
            for member in all_members
        ),
        "aurorapic_repeatability": all(
            aurora_ranges[metric] <= 0.10
            for metric in REPEATABILITY_METRICS
        ),
        "native_repeatability": all(
            native_ranges[metric] <= 0.15
            for metric in REPEATABILITY_METRICS
        ),
    }
    all_gates = all(gates.values())

    supply_supported = all(
        member["band_supply_fraction"] <= 0.90 for member in ratios)
    positive_work_supported = all(
        member["mean_positive_work_eV"] <= 0.90 for member in ratios)
    conversion_supported = all(
        member["band_promotion_probability"] <= 0.90 for member in ratios)
    outcome = {
        "interpretation_allowed": all_gates,
        "band_supply_deficit_supported": all_gates and supply_supported,
        "positive_work_deficit_supported": all_gates and positive_work_supported,
        "promotion_conversion_deficit_supported":
            all_gates and conversion_supported,
        "result": (
            "positive_work_and_promotion_conversion_deficits_supported_"
            "band_supply_deficit_not_supported"
            if all_gates and positive_work_supported and conversion_supported
               and not supply_supported
            else "mixed_or_intermediate_near_threshold_work_result"
        ),
    }
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "aurorapic_edipic_promotion_band_work_result",
        "rule_sha256": rule_hash,
        "gates": gates,
        "all_integrity_population_repeatability_and_work_closure_gates_passed":
            all_gates,
        "critical_scope": {
            "region": diagnostic["critical_regions"][0],
            "phase_fraction": diagnostic["critical_phase_fraction"],
            "promotion_band_eV": diagnostic["promotion_band_eV"],
        },
        "aurorapic_members": aurora_members,
        "aurorapic_ensemble_means": aurora_means,
        "aurorapic_relative_ranges": aurora_ranges,
        "native_edipic_ensemble_means": native_means,
        "native_edipic_relative_ranges": native_ranges,
        "aurorapic_member_to_native_mean_ratios": ratios,
        "aurorapic_ensemble_mean_to_native_mean_ratios": mean_ratios,
        "ensemble_mean_ratio_classification": {
            metric: classify_ratio(value, parity)
            for metric, value in mean_ratios.items()
        },
        "prospective_decision_outcome": outcome,
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "native_result_sha256": sha256(native_path),
            "aurorapic_runner_report_sha256": [
                sha256(path) for path in aurora_report_paths
            ],
        },
        "interpretation_note": (
            "Promotion probability is a conditional outcome of the within-band "
            "energy distribution and field work; the two supported deficits are "
            "not asserted to be statistically or causally independent."
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
        "all_integrity_population_repeatability_and_work_closure_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
