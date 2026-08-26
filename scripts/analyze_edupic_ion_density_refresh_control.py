#!/usr/bin/env python3
"""Analyze the preregistered held-density control against native and baseline."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from analyze_edupic_common_state_divergence import (
    read_aurora, read_native, relative_rms, sha256)
from analyze_edupic_grid_field_sampling import spatial_mean_square


def analyze(control_rule_path: Path, execution_lock_path: Path,
            parent_rule_path: Path, baseline_result_path: Path,
            native_trace_path: Path, native_population_path: Path,
            control_report_path: Path, control_root: Path,
            runner_path: Path, binary_path: Path) -> dict[str, object]:
    rule = json.loads(control_rule_path.read_text(encoding="utf-8"))
    lock = json.loads(execution_lock_path.read_text(encoding="utf-8"))
    parent = json.loads(parent_rule_path.read_text(encoding="utf-8"))
    baseline_result = json.loads(baseline_result_path.read_text(encoding="utf-8"))
    report = json.loads(control_report_path.read_text(encoding="utf-8"))
    sampling = rule["sampling_contract"]
    steps = sampling["edupic_pre_push_steps"]
    horizons = sampling["matching_aurorapic_post_step_horizons"]
    nodes = int(parent["physics_contract"]["nodes"])
    length = float(parent["physics_contract"]["length_m"])
    native = read_native(native_trace_path, steps, nodes, length)
    with native_population_path.open(newline="", encoding="utf-8") as stream:
        populations = {int(row["pre_push_step"]):
                       (int(row["electrons"]), int(row["ions"]))
                       for row in csv.DictReader(stream)}
    baseline = {int(item["aurorapic_horizon"]): item
                for item in baseline_result["comparisons"]}
    members = {int(item["horizon"]): item for item in report["members"]}
    comparisons = []
    for step, horizon in zip(steps, horizons, strict=True):
        native_rho, native_field = native[step]
        field_path = control_root / f"horizon-{horizon:04d}" / "output" / f"fields_{horizon}.csv"
        rho, field = read_aurora(field_path, nodes, length)
        member = members[horizon]
        if sha256(field_path) != member["field_sha256"]:
            raise ValueError("control field hash differs from runner report")
        native_e2 = spatial_mean_square(
            native_field, length, .2 * length, .4 * length)
        ratio = spatial_mean_square(
            field, length, .2 * length, .4 * length) / native_e2
        baseline_ratio = float(baseline[horizon]["critical_field_energy_ratio"])
        baseline_log_error = abs(math.log(baseline_ratio))
        log_error = abs(math.log(ratio))
        native_e, native_i = populations[step]
        comparisons.append({
            "edupic_pre_push_step": step,
            "aurorapic_horizon": horizon,
            "charge_relative_rms": relative_rms(rho, native_rho, length),
            "field_relative_rms": relative_rms(field, native_field, length),
            "critical_field_energy_ratio": ratio,
            "baseline_field_relative_rms": baseline[horizon]["field_relative_rms"],
            "critical_log_error_fraction_of_baseline":
                0.0 if horizon == 0 else log_error / max(baseline_log_error, 1e-300),
            "electron_population_difference": member["electron_population"] - native_e,
            "ion_population_difference": member["ion_population"] - native_i,
        })
    integrity = {
        "control_rule_hash_matches_lock": sha256(control_rule_path) == lock["control_rule_sha256"],
        "binary_hash_matches_lock": sha256(binary_path) == lock["aurorapic_binary_sha256"] == report["binary_sha256"],
        "runner_hash_matches_lock": sha256(runner_path) == lock["runner_sha256"],
        "particle_state_hash_matches": report["particle_state_sha256"] == rule["locked_inputs"]["particle_state_sha256"],
        "parent_rule_hash_matches": sha256(parent_rule_path) == report["parent_rule_sha256"] == rule["basis"]["common_state_rule_sha256"],
        "baseline_result_hash_matches": sha256(baseline_result_path) == rule["locked_inputs"]["baseline_result_sha256"],
        "all_resource_gates_passed": report["all_resource_gates_passed"] is True,
        "all_horizons_complete": set(members) == set(horizons),
        "native_populations_complete": set(populations) == set(steps),
    }
    initial = comparisons[0]
    initial_parity = (initial["charge_relative_rms"] <= 2e-6 and
                      initial["field_relative_rms"] <= 2e-6 and
                      .995 <= initial["critical_field_energy_ratio"] <= 1.005)
    early = [item for item in comparisons if item["aurorapic_horizon"] in {1, 2, 5}]
    strong = (all(item["critical_log_error_fraction_of_baseline"] <= .2
                  for item in early) and
              all(item["field_relative_rms"] <= item["baseline_field_relative_rms"]
                  for item in early) and
              all(item["electron_population_difference"] == 0 and
                  item["ion_population_difference"] == 0 for item in early))
    material = [(item["field_relative_rms"] > .01 or not .98 <=
                 item["critical_field_energy_ratio"] <= 1.02)
                for item in comparisons]
    through_20 = [flag for item, flag in zip(comparisons, material, strict=True)
                  if item["aurorapic_horizon"] <= 20]
    early_explained = strong and not any(through_20)
    fractions = sorted(item["critical_log_error_fraction_of_baseline"] for item in early)
    partial = (fractions[1] <= .5 and all(
        item["field_relative_rms"] <= 1.2 * item["baseline_field_relative_rms"]
        for item in early))
    integrity_passed = all(integrity.values())
    if not integrity_passed or not initial_parity:
        outcome = "inconclusive_failed_gate"
    elif strong:
        outcome = "strong_ion_density_refresh_mechanism_support"
    elif partial:
        outcome = "partial_ion_density_refresh_mechanism_support"
    else:
        outcome = "ion_density_refresh_mechanism_not_supported"
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "held_density_control_result",
        "provenance": {
            "control_rule_sha256": sha256(control_rule_path),
            "execution_lock_sha256": sha256(execution_lock_path),
            "control_report_sha256": sha256(control_report_path),
            "native_trace_sha256": sha256(native_trace_path),
            "native_population_sha256": sha256(native_population_path),
        },
        "integrity": integrity,
        "integrity_gate_passed": integrity_passed,
        "initial_parity_gate_passed": initial_parity,
        "comparisons": comparisons,
        "material_field_divergence_flags": material,
        "strong_mechanism_support_gate_passed": strong,
        "early_divergence_explained_through_horizon_20": early_explained,
        "formal_outcome": outcome,
        "residual_observation": "The first remaining material flag is at horizon 20, the next ion-cache refresh boundary; its cadence alignment requires a separate prospective test.",
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("control_rule", "execution_lock", "parent_rule",
                 "baseline_result", "native_trace", "native_population",
                 "control_report", "control_root", "runner", "binary", "output"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.control_rule, args.execution_lock, args.parent_rule,
        args.baseline_result, args.native_trace, args.native_population,
        args.control_report, args.control_root, args.runner, args.binary)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "formal_outcome", "strong_mechanism_support_gate_passed",
        "early_divergence_explained_through_horizon_20")}, indent=2))


if __name__ == "__main__":
    main()
