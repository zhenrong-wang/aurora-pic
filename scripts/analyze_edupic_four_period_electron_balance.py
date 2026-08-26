#!/usr/bin/env python3
"""Resolve the five-member four-period electron-balance equivalence test."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_edupic_collision_enabled_common_state_ensemble import analyze, sha256


CONSERVATIVE_T4_95 = 2.131846786326649


def equivalence_interval(
        edupic_mean: float, edupic_sd: float,
        aurora_mean: float, aurora_sd: float,
        members_each: int = 5) -> dict[str, object]:
    difference = aurora_mean - edupic_mean
    standard_error = math.sqrt(
        edupic_sd * edupic_sd / members_each +
        aurora_sd * aurora_sd / members_each)
    half_width = CONSERVATIVE_T4_95 * standard_error
    return {
        "aurorapic_minus_edupic_mean_particles": difference,
        "standard_error_particles": standard_error,
        "critical_value": CONSERVATIVE_T4_95,
        "critical_value_contract": "t_0.95_with_4_df_conservative_lower_df",
        "confidence_level": 0.90,
        "lower_particles": difference - half_width,
        "upper_particles": difference + half_width,
    }


def resolve(rule_path: Path, lock_path: Path, report_path: Path,
            output_root: Path) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    result = analyze(rule_path, lock_path, report_path, output_root)
    parent_hash_matches = (
        report.get("parent_execution_sha256") ==
        rule["basis"]["parent_four_period_execution_sha256"])
    reused_verified = bool(report.get("reused_parent_members_verified", False))
    result["integrity"]["parent_execution_hash_matches"] = parent_hash_matches
    result["integrity"]["reused_parent_members_verified"] = reused_verified
    integrity_passed = all(result["integrity"].values())
    result["integrity_passed"] = integrity_passed

    electrons = result["populations"]["electrons"]
    interval = equivalence_interval(
        electrons["edupic"]["mean"],
        electrons["edupic"]["sample_standard_deviation"],
        electrons["aurorapic"]["mean"],
        electrons["aurorapic"]["sample_standard_deviation"])
    margin = 0.001 * float(
        rule["locked_inputs"]["initial_populations"]["electrons"])
    equivalent = (interval["lower_particles"] >= -margin and
                  interval["upper_particles"] <= margin)
    interval["equivalence_margin_particles"] = margin
    interval["equivalence_interval_particles"] = [-margin, margin]
    interval["passed"] = equivalent
    result["electron_balance_equivalence"] = interval

    ion_passed = bool(result["populations"]["ions"]["passed"])
    support = (integrity_passed and equivalent and ion_passed and
               result["collision_traffic_compatible"] and
               not result["wall_traffic_has_failure"] and
               result["field_compatible"])
    result.pop("four_period_support_gate_passed", None)
    result["five_member_support_gate_passed"] = support
    result["population_compatible"] = equivalent and ion_passed
    result["scope"] = "five_member_four_period_electron_balance_result"
    if not integrity_passed:
        result["formal_outcome"] = "inconclusive_integrity_failure"
    elif support:
        result["formal_outcome"] = (
            "four_period_electron_balance_equivalence_supported")
    else:
        result["formal_outcome"] = (
            "four_period_electron_balance_equivalence_not_supported")
    result["claim_boundary"] = rule["claim_boundary"]
    result["parent_execution_report_sha256"] = (
        rule["basis"]["parent_four_period_execution_sha256"])
    result["execution_report_sha256"] = sha256(report_path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    args = parser.parse_args()
    result = resolve(args.rule.resolve(), args.lock.resolve(),
                     args.report.resolve(), args.output_root.resolve())
    args.result.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
