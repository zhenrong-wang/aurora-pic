#!/usr/bin/env python3
"""Audit post-trend equilibration blocks and forecast the frozen horizon."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics

from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = "3e6e29425325e9e70557fa1a17545893fae816e49ff351bac14dc4dd82f37b27"
REPORT_HASHES = {
    120: "ac0ad23b00d83b23fe8d6e411834c308c79b0922604126e0921ac39c5b1f8c5b",
    124: "4b304e9cb97d76c7af4febf359ba5124471e0304918eec187799e77ca792bf61",
    128: "f6c54f82d00fccac6f081af0d8b2dd0ddfa1801cfa259d60be5327113827fd8d",
}


def linear_forecast(x: list[float], y: list[float], target: float) -> float:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("forecast requires equal vectors with two points")
    x_mean, y_mean = statistics.fmean(x), statistics.fmean(y)
    denominator = math.fsum((value - x_mean) ** 2 for value in x)
    if denominator == 0.0:
        raise ValueError("forecast coordinates are degenerate")
    slope = math.fsum((a - x_mean) * (b - y_mean)
                      for a, b in zip(x, y)) / denominator
    return y_mean + slope * (target - x_mean)


def analyze(rule_path: Path, report_paths: list[Path]) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("locked post-trend rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    reports = []
    for path in report_paths:
        report = json.loads(path.read_text(encoding="utf-8"))
        end_cycle = int(report["block"]["end_cycle"])
        if end_cycle not in REPORT_HASHES or sha256(path) != REPORT_HASHES[end_cycle]:
            raise ValueError("post-trend block is absent or differs")
        reports.append(report)
    reports.sort(key=lambda item: item["block"]["end_cycle"])
    if [item["block"]["end_cycle"] for item in reports] != [120, 124, 128]:
        raise ValueError("post-trend blocks are incomplete")
    previous_hash = rule["baseline"]["measurement_report_sha256"]
    for report in reports:
        if (report["inputs"]["prior_report_sha256"] != previous_hash or
                report["inputs"]["extension_rule_sha256"] != RULE_SHA256 or
                report["block"]["hard_safety_gates_passed"] is not True):
            raise ValueError("post-trend evidence chain is invalid")
        previous_hash = REPORT_HASHES[report["block"]["end_cycle"]]

    cycles = [float(item["block"]["end_cycle"]) for item in reports]
    metric_names = (
        "normalized_total_population_slope_per_cycle",
        "electron_source_loss_relative_imbalance",
        "ion_source_loss_relative_imbalance",
    )
    history = {}
    forecasts = {}
    for name in metric_names:
        values = [float(item["stationarity_screen"]["metrics"][name])
                  for item in reports]
        history[name] = values
        forecasts[name] = linear_forecast(cycles, values, 148.0)
    thresholds = reports[-1]["stationarity_screen"]["thresholds"]
    projected_gates = {
        "total_population_slope": abs(forecasts[
            "normalized_total_population_slope_per_cycle"]) <= float(
                thresholds[
                    "maximum_absolute_normalized_total_population_slope_per_cycle"]),
        "electron_source_loss_balance": abs(forecasts[
            "electron_source_loss_relative_imbalance"]) <= float(
                thresholds["maximum_electron_source_loss_relative_imbalance"]),
        "ion_source_loss_balance": abs(forecasts[
            "ion_source_loss_relative_imbalance"]) <= float(
                thresholds["maximum_ion_source_loss_relative_imbalance"]),
    }
    stages = [stage for report in reports for stage in report["stages"]]
    final = reports[-1]
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "post_trend_equilibration_horizon_sufficiency_audit",
        "evidence_sha256": {
            "rule": RULE_SHA256,
            "blocks": {str(cycle): digest
                       for cycle, digest in REPORT_HASHES.items()},
        },
        "completed": {
            "first_cycle": 117,
            "last_cycle": 128,
            "cycles": 12,
            "all_hard_safety_gates_passed": True,
            "maximum_peak_resident_set_kib": max(
                int(stage["resources"]["peak_resident_set_kib"])
                for stage in stages),
            "maximum_wall_seconds_per_cycle": max(
                float(stage["resources"]["wall_seconds"])
                for stage in stages),
            "final_checkpoint_sha256":
                final["stages"][-1]["output_checkpoint_sha256"],
        },
        "block_end_cycles": [int(value) for value in cycles],
        "metric_history": history,
        "linear_horizon_forecast_at_cycle148": forecasts,
        "projected_strict_gates": projected_gates,
        "all_projected_strict_gates_pass": all(projected_gates.values()),
        "decision": {
            "launch_remaining_declared_blocks_now": False,
            "campaign_complete": False,
            "comparison_measurement_eligible": False,
            "reason": (
                "All numerical safety gates pass, but recent block trends do "
                "not project population drift or source/loss balance near their "
                "frozen strict thresholds by cycle 148. Pause between authorized "
                "invocations and redesign equilibration rather than spend the "
                "remaining horizon without a plausible acceptance outcome."),
        },
        "forecast_boundary": (
            "The three-block linear forecast is a resource-planning diagnostic, "
            "not a physical asymptote, stopping proof, or validation result."),
        "claim_boundary": rule["claim_boundary"],
        "physics_claim": "none_resource_planning_and_stationarity_rejection_only",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("reports", type=Path, nargs=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.rule.resolve(),
                     [path.resolve() for path in args.reports])
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
