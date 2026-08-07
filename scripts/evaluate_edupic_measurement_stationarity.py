#!/usr/bin/env python3
"""Apply a predeclared eduPIC density-stationarity rule to block analysis."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

from run_edupic_stage import atomic_json, sha256


class StationarityEvaluationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise StationarityEvaluationError(message)


def finite_number(value: object, name: str) -> float:
    require(isinstance(value, (int, float)) and math.isfinite(float(value)),
            f"{name} is not a finite number")
    return float(value)


def read_json(path: Path, name: str) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise StationarityEvaluationError(f"cannot read {name}: {error}") from error
    require(isinstance(value, dict), f"{name} is not a JSON object")
    return value


def maximum_gate(threshold: float, value: float) -> dict:
    return {"threshold": threshold, "value": value,
            "passed": value <= threshold}


def minimum_gate(threshold: float, value: float) -> dict:
    return {"threshold": threshold, "value": value,
            "passed": value >= threshold}


def evaluate(analysis_path: Path, rule_path: Path) -> dict:
    analysis = read_json(analysis_path, "analysis")
    rule = read_json(rule_path, "stationarity rule")
    require(rule.get("scope") ==
            "predeclared_internal_native_measurement_stationarity_screen",
            "input rule has the wrong scope")
    require(analysis.get("scope") == "native_measurement_block_analysis",
            "input analysis has the wrong scope")
    require(rule.get("case_id") == analysis.get("case_id"),
            "analysis and rule case identities differ")

    contract = rule.get("analysis_contract")
    thresholds = rule.get("internal_density_stationarity_gates")
    campaign = analysis.get("campaign")
    profiles = analysis.get("density_profile")
    require(all(isinstance(value, dict) for value in
                (contract, thresholds, campaign, profiles)),
            "analysis or rule contract is incomplete")

    required_blocks = int(finite_number(
        contract.get("required_total_contiguous_blocks"),
        "required total blocks"))
    required_cycles = int(finite_number(
        contract.get("required_total_measurement_cycles"),
        "required total measurement cycles"))
    minimum_effective = finite_number(
        thresholds.get("minimum_ar1_effective_blocks_per_species"),
        "minimum effective blocks")
    maximum_drift = finite_number(
        thresholds.get("maximum_absolute_projected_fractional_drift_per_species"),
        "maximum projected drift")
    maximum_split = finite_number(
        thresholds.get("maximum_absolute_split_half_fractional_change_per_species"),
        "maximum split-half change")
    maximum_movement = finite_number(
        thresholds.get("maximum_adjacent_profile_relative_l2_per_species"),
        "maximum adjacent profile movement")

    blocks = int(finite_number(campaign.get("block_count"), "block count"))
    cycles = int(finite_number(
        campaign.get("completed_measurement_cycles"),
        "completed measurement cycles"))
    horizon = {
        "minimum_total_blocks": minimum_gate(required_blocks, blocks),
        "minimum_total_measurement_cycles": minimum_gate(required_cycles, cycles),
        "campaign_target_reached": {
            "required": True,
            "value": campaign.get("target_reached"),
            "passed": campaign.get("target_reached") is True,
        },
        "analysis_eligible": {
            "required": True,
            "value": analysis.get("analysis_eligible"),
            "passed": analysis.get("analysis_eligible") is True,
        },
    }
    horizon_complete = all(gate["passed"] for gate in horizon.values())

    species_gates: dict[str, dict] = {}
    for species in ("electron", "ion"):
        series = profiles.get(f"{species}_series")
        require(isinstance(series, dict), f"{species} density series is absent")
        effective = finite_number(series.get("ar1_effective_blocks"),
                                  f"{species} effective blocks")
        drift = abs(finite_number(
            series.get("projected_fractional_drift_across_series"),
            f"{species} projected drift"))
        split = abs(finite_number(series.get("split_half_fractional_change"),
                                  f"{species} split-half change"))
        movement = finite_number(
            series.get("maximum_adjacent_profile_relative_l2"),
            f"{species} maximum adjacent profile movement")
        species_gates[species] = {
            "minimum_ar1_effective_blocks": minimum_gate(
                minimum_effective, effective),
            "maximum_absolute_projected_fractional_drift": maximum_gate(
                maximum_drift, drift),
            "maximum_absolute_split_half_fractional_change": maximum_gate(
                maximum_split, split),
            "maximum_adjacent_profile_relative_l2": maximum_gate(
                maximum_movement, movement),
        }

    stationarity_passed = (
        horizon_complete and all(
            gate["passed"] for gates in species_gates.values()
            for gate in gates.values()))
    if not horizon_complete:
        classification = "stationarity_horizon_incomplete"
        decision: bool | None = None
    elif stationarity_passed:
        classification = "internal_density_stationarity_screen_passed"
        decision = True
    else:
        classification = "internal_density_stationarity_screen_failed"
        decision = False
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "predeclared_internal_native_measurement_stationarity_evaluation",
        "rule": {"path": str(rule_path.resolve()), "sha256": sha256(rule_path)},
        "analysis": {
            "path": str(analysis_path.resolve()),
            "sha256": sha256(analysis_path),
        },
        "horizon": horizon,
        "species_gates": species_gates,
        "horizon_complete": horizon_complete,
        "passed": decision,
        "classification": classification,
        "claim_boundary": (
            "This applies predeclared internal density sampling-readiness gates. "
            "It is not cross-code, experimental, or independent physical validation."),
        "physics_claim": "none_internal_sampling_readiness_only",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("analysis", type=Path)
    parser.add_argument("rule", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-pass", action="store_true")
    args = parser.parse_args()
    try:
        report = evaluate(args.analysis, args.rule)
    except (StationarityEvaluationError, OSError, json.JSONDecodeError) as error:
        print(f"eduPIC stationarity evaluation failed: {error}", file=sys.stderr)
        return 2
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        atomic_json(args.output, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 1 if args.require_pass and report["passed"] is not True else 0


if __name__ == "__main__":
    raise SystemExit(main())
