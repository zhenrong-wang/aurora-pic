#!/usr/bin/env python3
"""Analyze the preregistered upcoming-due held-density cadence control."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from analyze_edupic_common_state_divergence import (
    read_aurora, read_native, relative_rms, sha256)
from analyze_edupic_grid_field_sampling import spatial_mean_square


def analyze(rule_path: Path, lock_path: Path, parent_rule_path: Path,
            first_result_path: Path, first_report_path: Path,
            native_trace_path: Path, native_population_path: Path,
            report_path: Path, root: Path, runner_path: Path,
            binary_path: Path) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    parent = json.loads(parent_rule_path.read_text(encoding="utf-8"))
    first = json.loads(first_result_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    steps = parent["sampling_contract"]["edupic_pre_push_steps"]
    horizons = rule["sampling_contract"]["matching_horizons"]
    nodes = int(parent["physics_contract"]["nodes"])
    length = float(parent["physics_contract"]["length_m"])
    native = read_native(native_trace_path, steps, nodes, length)
    with native_population_path.open(newline="", encoding="utf-8") as stream:
        populations = {int(row["pre_push_step"]):
                       (int(row["electrons"]), int(row["ions"]))
                       for row in csv.DictReader(stream)}
    first_by_horizon = {int(item["aurorapic_horizon"]): item
                        for item in first["comparisons"]}
    members = {int(item["horizon"]): item for item in report["members"]}
    comparisons = []
    for step, horizon in zip(steps, horizons, strict=True):
        native_rho, native_field = native[step]
        field_path = root / f"horizon-{horizon:04d}" / "output" / f"fields_{horizon}.csv"
        rho, field = read_aurora(field_path, nodes, length)
        member = members[horizon]
        if sha256(field_path) != member["field_sha256"]:
            raise ValueError("cadence-control field hash differs from report")
        ratio = (
            spatial_mean_square(field, length, .2 * length, .4 * length) /
            spatial_mean_square(
                native_field, length, .2 * length, .4 * length))
        previous = first_by_horizon[horizon]
        previous_error = abs(math.log(previous["critical_field_energy_ratio"]))
        native_e, native_i = populations[step]
        comparisons.append({
            "edupic_pre_push_step": step,
            "aurorapic_horizon": horizon,
            "charge_relative_rms": relative_rms(rho, native_rho, length),
            "field_relative_rms": relative_rms(field, native_field, length),
            "critical_field_energy_ratio": ratio,
            "first_control_field_relative_rms": previous["field_relative_rms"],
            "critical_log_error_fraction_of_first_control":
                0.0 if horizon == 0 else abs(math.log(ratio)) /
                    max(previous_error, 1e-300),
            "electron_population_difference": member["electron_population"] - native_e,
            "ion_population_difference": member["ion_population"] - native_i,
        })
    integrity = {
        "cadence_rule_hash_matches_lock": sha256(rule_path) == lock["cadence_rule_sha256"],
        "binary_hash_matches_lock": sha256(binary_path) == lock["aurorapic_binary_sha256"] == report["binary_sha256"],
        "runner_hash_matches_lock": sha256(runner_path) == lock["runner_sha256"],
        "particle_state_hash_matches": report["particle_state_sha256"] == rule["locked_inputs"]["particle_state_sha256"],
        "first_result_hash_matches": sha256(first_result_path) == rule["locked_inputs"]["first_control_result_sha256"],
        "first_report_hash_matches": sha256(first_report_path) == rule["locked_inputs"]["first_control_report_sha256"],
        "native_trace_hash_matches": sha256(native_trace_path) == rule["locked_inputs"]["native_trace_sha256"],
        "native_population_hash_matches": sha256(native_population_path) == rule["locked_inputs"]["native_population_sha256"],
        "all_resource_gates_passed": report["all_resource_gates_passed"] is True,
        "all_horizons_complete": set(members) == set(horizons),
    }
    initial = comparisons[0]
    initial_parity = (initial["charge_relative_rms"] <= 2e-6 and
                      initial["field_relative_rms"] <= 2e-6 and
                      .995 <= initial["critical_field_energy_ratio"] <= 1.005)
    primary = [item for item in comparisons
               if item["aurorapic_horizon"] in {20, 50, 100}]
    cadence_supported = (all(
        item["critical_log_error_fraction_of_first_control"] <= .2 and
        item["field_relative_rms"] <= item["first_control_field_relative_rms"] and
        item["electron_population_difference"] == 0 and
        item["ion_population_difference"] == 0 for item in primary))
    material = [(item["field_relative_rms"] > .01 or
                 not .98 <= item["critical_field_energy_ratio"] <= 1.02)
                for item in comparisons]
    early_closure = cadence_supported and not any(
        flag for item, flag in zip(comparisons, material, strict=True)
        if item["aurorapic_horizon"] <= 100)
    full_closure = (not any(material) and all(
        item["electron_population_difference"] == 0 and
        item["ion_population_difference"] == 0 for item in comparisons))
    integrity_passed = all(integrity.values())
    if not integrity_passed or not initial_parity:
        outcome = "inconclusive_failed_gate"
    elif full_closure:
        outcome = "full_collision_free_common_state_trace_closure"
    elif cadence_supported:
        outcome = "held_density_refresh_cadence_supported"
    else:
        outcome = "held_density_refresh_cadence_not_supported"
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "upcoming_due_held_density_cadence_result",
        "provenance": {"cadence_rule_sha256": sha256(rule_path),
                       "execution_lock_sha256": sha256(lock_path),
                       "control_report_sha256": sha256(report_path)},
        "integrity": integrity,
        "integrity_gate_passed": integrity_passed,
        "initial_parity_gate_passed": initial_parity,
        "comparisons": comparisons,
        "material_field_divergence_flags": material,
        "cadence_support_gate_passed": cadence_supported,
        "early_trajectory_closure_through_horizon_100": early_closure,
        "full_trace_closure_gate_passed": full_closure,
        "formal_outcome": outcome,
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("rule", "lock", "parent_rule", "first_result",
                 "first_report", "native_trace", "native_population",
                 "report", "root", "runner", "binary", "output"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.rule, args.lock, args.parent_rule,
        args.first_result, args.first_report, args.native_trace,
        args.native_population, args.report, args.root, args.runner, args.binary)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "formal_outcome", "cadence_support_gate_passed",
        "early_trajectory_closure_through_horizon_100",
        "full_trace_closure_gate_passed")}, indent=2))


if __name__ == "__main__":
    main()
