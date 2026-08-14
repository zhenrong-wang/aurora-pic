#!/usr/bin/env python3
"""Audit phase/order-matched AuroraPIC electron heating against eduPIC."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_aurorapic_mesh_refinement import trapezoid
from compare_aurorapic_edupic_measurement_pilot import rows
from compare_edupic_phase_space import (
    REFERENCE_FILES, flatten_phase_major, read_matrix, spatial_phase_average,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


BRANCH_REPORT_SHA256 = (
    "4171678276ffb9fee0eecb843f9d9c0b44ecb5bbf13c65207618a9cd76820fc5")
RULE_SHA256 = (
    "3ae87d834ebbec16d8b59f964859a9fa398d6dca85f1d73dd204fd24a5d41842")
ELEMENTARY_CHARGE_C = 1.60217662e-19
LENGTH_M = 0.025
MAXIMUM_INTERNAL_POWER_DIFFERENCE = 0.02


def summarize(candidate_density: float, binned_power: float,
              exact_power: float, reference_density: float,
              reference_power: float) -> dict[str, float | bool]:
    values = (candidate_density, binned_power, exact_power,
              reference_density, reference_power)
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError("matched-heating inputs must be positive and finite")
    internal_difference = abs(binned_power / exact_power - 1.0)
    reference_per_particle = reference_power / reference_density
    binned_ratio = (binned_power / candidate_density) / reference_per_particle
    exact_ratio = (exact_power / candidate_density) / reference_per_particle
    return {
        "candidate_volume_phase_average_density_m-3": candidate_density,
        "candidate_phase_binned_j_dot_e_W_m-3": binned_power,
        "candidate_exact_electric_work_W_m-3": exact_power,
        "candidate_phase_binned_to_exact_power_ratio":
            binned_power / exact_power,
        "candidate_phase_binned_to_exact_relative_difference":
            internal_difference,
        "maximum_declared_internal_relative_difference":
            MAXIMUM_INTERNAL_POWER_DIFFERENCE,
        "internal_power_gate_passed":
            internal_difference <= MAXIMUM_INTERNAL_POWER_DIFFERENCE,
        "reference_volume_phase_average_density_m-3": reference_density,
        "reference_phase_binned_j_dot_e_W_m-3": reference_power,
        "candidate_to_reference_phase_binned_power_per_particle_ratio":
            binned_ratio,
        "candidate_exact_to_reference_binned_power_per_particle_ratio":
            exact_ratio,
        "matched_cross_code_power_per_particle_deficit_fraction":
            1.0 - binned_ratio,
    }


def phase_average(x: list[float], values: list[float], phases: int) -> float:
    nodes = len(x)
    if phases == 0 or nodes < 2 or len(values) != phases * nodes:
        raise ValueError("matched phase-space shape differs")
    length = x[-1] - x[0]
    if length <= 0.0:
        raise ValueError("candidate spatial coordinates are invalid")
    return math.fsum(
        trapezoid(x, values[phase * nodes:(phase + 1) * nodes]) / length
        for phase in range(phases)) / phases


def analyze(root: Path, reference: Path) -> dict[str, object]:
    report_path = root / "branch-report.json"
    if sha256(report_path) != BRANCH_REPORT_SHA256:
        raise ValueError("matched-heating branch report differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("rule_sha256") != RULE_SHA256 or not report.get(
            "all_gates_passed"):
        raise ValueError("matched-heating branch contract did not pass")
    output = root / "measurement" / "output"
    for name, expected in report["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"matched candidate output differs: {name}")
    metadata = json.loads(
        (output / "spatial_average_metadata.json").read_text(encoding="utf-8"))
    if (metadata.get("sampling_order") != "pre_collision" or
            metadata.get("phase_bins") != 200 or
            metadata.get("complete") is not True or
            set(metadata.get("phase_bin_samples", [])) != {80}):
        raise ValueError("candidate phase/order sampling contract differs")
    for key in ("electron_density", "electron_ohmic_power_density"):
        name, expected = REFERENCE_FILES[key]
        if sha256(reference / name) != expected:
            raise ValueError(f"locked eduPIC reference differs: {name}")

    fields = rows(output / "spatial_phase_fields.csv")
    moments = [row for row in rows(output / "spatial_phase_moments.csv")
               if row["species"] == "electrons"]
    phases = len({int(row["phase_bin"]) for row in fields})
    nodes = len(fields) // phases
    if len(fields) != len(moments) or phases * nodes != len(fields):
        raise ValueError("candidate field and moment shapes differ")
    x = [float(row["x_m"]) for row in fields[:nodes]]
    density = [float(row["number_density_mean_m-3"]) for row in moments]
    power = [
        -ELEMENTARY_CHARGE_C * number * float(moment["mean_velocity_x"]) *
        float(field["electric_field_mean_V_m"])
        for number, moment, field in zip(density, moments, fields)
    ]
    energy = json.loads((output / "energy-budget.json").read_text())
    exact_power = float(
        energy["electric_power_W_m-2"]["electric_work_electrons_J_m-2"]
    ) / LENGTH_M
    reference_density = spatial_phase_average(
        flatten_phase_major(read_matrix(
            reference / REFERENCE_FILES["electron_density"][0])), 200, 400)
    reference_power = spatial_phase_average(
        flatten_phase_major(read_matrix(
            reference / REFERENCE_FILES["electron_ohmic_power_density"][0])),
        200, 400)
    metrics = summarize(
        phase_average(x, density, phases), phase_average(x, power, phases),
        exact_power, reference_density, reference_power)
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "phase_and_order_matched_electron_heating_audit",
        "branch_report_sha256": BRANCH_REPORT_SHA256,
        "rule_sha256": RULE_SHA256,
        "candidate_phase_bins": phases,
        "candidate_sampling_order": "pre_collision",
        "reference_phase_bins": 200,
        "reference_sampling_order": "pre_collision",
        "metrics": metrics,
        "finding": (
            "Matched phase/order J.E closes against AuroraPIC's exact work "
            "ledger, while a material cross-code heating deficit remains."),
        "claim_boundary": (
            "This removes one diagnostic confounder for one stationary state; "
            "it is not predictive or experimental validation."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_root", type=Path)
    parser.add_argument("reference_raw_data", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.candidate_root.resolve(),
                     args.reference_raw_data.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
