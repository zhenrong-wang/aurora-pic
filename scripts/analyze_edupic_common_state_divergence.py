#!/usr/bin/env python3
"""Analyze the locked collision-free common-state divergence trace."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from analyze_edupic_grid_field_sampling import spatial_mean_square


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mean_square(values: list[float], length: float) -> float:
    return spatial_mean_square(values, length, 0.0, length)


def relative_rms(candidate: list[float], reference: list[float],
                 length: float) -> float:
    return math.sqrt(mean_square([a - b for a, b in zip(candidate, reference,
                                                        strict=True)], length) /
                     max(mean_square(reference, length), 1e-300))


def read_native(path: Path, steps: list[int], nodes: int,
                length: float) -> dict[int, tuple[list[float], list[float]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != len(steps) * nodes:
        raise ValueError("native trace shape differs")
    result = {}
    for index, step in enumerate(steps):
        block = rows[index * nodes:(index + 1) * nodes]
        rho, field = [], []
        for node, row in enumerate(block):
            if int(row["pre_push_step"]) != step or int(row["node"]) != node:
                raise ValueError("native trace ordering differs")
            if abs(float(row["x_m"]) - node * length / (nodes - 1)) > 1e-12:
                raise ValueError("native trace coordinates differ")
            rho.append(float(row["charge_density_C_m3"]))
            field.append(float(row["electric_field_V_m"]))
        result[step] = rho, field
    return result


def read_aurora(path: Path, nodes: int, length: float) -> tuple[list[float], list[float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != nodes:
        raise ValueError("AuroraPIC field shape differs")
    rho, field = [], []
    for node, row in enumerate(rows):
        if abs(float(row["x"]) - node * length / (nodes - 1)) > 1e-12:
            raise ValueError("AuroraPIC field coordinates differ")
        values = float(row["rho"]), float(row["E"])
        if not all(math.isfinite(value) for value in values):
            raise ValueError("AuroraPIC field contains non-finite values")
        rho.append(values[0]); field.append(values[1])
    return rho, field


def analyze(rule_path: Path, native_trace: Path, native_population: Path,
            native_primary_checkpoint: Path, native_extended_checkpoint: Path,
            native_extended_trace: Path, aurora_report: Path,
            aurora_root: Path) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    report = json.loads(aurora_report.read_text(encoding="utf-8"))
    locked = rule["locked_inputs"]; sampling = rule["sampling_contract"]
    nodes = int(rule["physics_contract"]["nodes"])
    length = float(rule["physics_contract"]["length_m"])
    steps = sampling["edupic_pre_push_steps"]
    horizons = sampling["matching_aurorapic_post_step_horizons"]
    native = read_native(native_trace, steps, nodes, length)
    with native_population.open(newline="", encoding="utf-8") as stream:
        populations = {int(row["pre_push_step"]):
                       (int(row["electrons"]), int(row["ions"]))
                       for row in csv.DictReader(stream)}
    report_members = {int(item["horizon"]): item for item in report["members"]}
    comparisons = []
    lower, upper = .2 * length, .4 * length
    for step, horizon in zip(steps, horizons, strict=True):
        native_rho, native_field = native[step]
        field_path = aurora_root / f"horizon-{horizon:04d}" / "output" / f"fields_{horizon}.csv"
        aurora_rho, aurora_field = read_aurora(field_path, nodes, length)
        member = report_members[horizon]
        if sha256(field_path) != member["field_sha256"]:
            raise ValueError("AuroraPIC field hash differs from runner report")
        native_e2 = spatial_mean_square(native_field, length, lower, upper)
        native_e, native_i = populations[step]
        comparisons.append({
            "edupic_pre_push_step": step, "aurorapic_horizon": horizon,
            "charge_relative_rms": relative_rms(aurora_rho, native_rho, length),
            "field_relative_rms": relative_rms(aurora_field, native_field, length),
            "critical_field_energy_ratio": spatial_mean_square(
                aurora_field, length, lower, upper) / native_e2,
            "electron_population_difference": member["electron_population"] - native_e,
            "ion_population_difference": member["ion_population"] - native_i,
        })
    integrity = {
        "rule_input_state_hash_matches": report["particle_state_sha256"] == locked["particle_state_sha256"],
        "rule_binary_hash_matches": report["binary_sha256"] == locked["aurorapic_binary_sha256"],
        "runner_rule_hash_matches": report["rule_sha256"] == sha256(rule_path),
        "all_runner_resource_gates_passed": report["all_resource_gates_passed"] is True,
        "population_rows_complete": set(populations) == set(steps),
        "passive_grid_trace_byte_identical": sha256(native_trace) == sha256(native_extended_trace),
        "passive_final_checkpoint_byte_identical": sha256(native_primary_checkpoint) == sha256(native_extended_checkpoint),
    }
    initial = comparisons[0]
    initial_parity = (initial["charge_relative_rms"] <= 2e-6 and
                      initial["field_relative_rms"] <= 2e-6 and
                      .995 <= initial["critical_field_energy_ratio"] <= 1.005)
    material = [(item["field_relative_rms"] > .01 or
                 not .98 <= item["critical_field_energy_ratio"] <= 1.02)
                for item in comparisons]
    sustained_index = next((index for index in range(len(material) - 2)
                            if all(material[index:index + 3])), None)
    earliest = None if sustained_index is None else horizons[sustained_index]
    integrity_passed = all(integrity.values())
    if not integrity_passed: outcome = "inconclusive_failed_integrity_gate"
    elif not initial_parity: outcome = "initial_state_or_deposition_mismatch"
    elif earliest in {1, 2}: outcome = "one_step_mover_or_boundary_mismatch"
    elif earliest is not None and earliest > 20: outcome = "later_nonlinear_divergence"
    elif earliest is None: outcome = "collision_free_trajectory_agreement"
    else: outcome = "mixed_or_transient_common_state_divergence_result"
    return {
        "schema_version": 1, "case_id": rule["case_id"],
        "scope": "collision_free_common_particle_state_divergence_result",
        "provenance": {"rule_sha256": sha256(rule_path),
                       "native_trace_sha256": sha256(native_trace),
                       "native_population_sha256": sha256(native_population),
                       "aurora_report_sha256": sha256(aurora_report)},
        "integrity": integrity, "integrity_gate_passed": integrity_passed,
        "initial_parity_gate_passed": initial_parity,
        "comparisons": comparisons,
        "material_field_divergence_flags": material,
        "earliest_sustained_material_divergence_horizon": earliest,
        "formal_outcome": outcome,
        "post_hoc_mechanism_candidate": {
            "name": "ion_density_refresh_staggering",
            "evidence": "eduPIC deposits ion density before each 20-step ion push and reuses that pre-push density between ion steps; AuroraPIC redeposits all species after the push. Initial profiles agree and the first localized field difference appears immediately after the first ion move, before any population difference.",
            "requires_prospective_control": True,
        },
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("rule", "native_trace", "native_population",
                 "native_primary_checkpoint", "native_extended_checkpoint",
                 "native_extended_trace", "aurora_report", "aurora_root", "output"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.rule, args.native_trace, args.native_population,
        args.native_primary_checkpoint, args.native_extended_checkpoint,
        args.native_extended_trace, args.aurora_report, args.aurora_root)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"formal_outcome": result["formal_outcome"],
                      "initial_parity_gate_passed": result["initial_parity_gate_passed"],
                      "earliest_sustained_material_divergence_horizon": result["earliest_sustained_material_divergence_horizon"]}, indent=2))


if __name__ == "__main__":
    main()
