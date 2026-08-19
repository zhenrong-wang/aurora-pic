#!/usr/bin/env python3
"""Evaluate the predeclared twelve-cycle matched-heating window."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_aurorapic_heating_localization import NODES, PHASES, reduced_profiles
from compare_aurorapic_edupic_measurement_pilot import rows
from compare_edupic_phase_space import (
    ELEMENTARY_CHARGE_C, REFERENCE_FILES, flatten_phase_major, metrics,
    read_matrix, spatial_phase_average,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "c85c7cbd9ede314ff5e744699fe630ac2ff19af74c0039c44cf2f5543fd0b2a0")
BRANCH_REPORT_SHA256 = (
    "1d87b9ebfbe3668513c844fc5314868e1bae625df259c7f659981f9ba03a9b4a")
ENSEMBLE_POWER_W_M3 = 1557.1549186535874
ENSEMBLE_DENSITY_M3 = 3339730898000002.5
ENSEMBLE_POWER_PER_ELECTRON_RATIO = 0.8725894597269717
LENGTH_M = 0.025


def gate_summary(power: float, exact_power: float, density: float,
                 power_per_electron_ratio: float) -> dict[str, object]:
    values = (power, exact_power, density, power_per_electron_ratio)
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError("long-window gate inputs must be positive and finite")
    differences = {
        "phase_binned_to_exact_power": abs(power / exact_power - 1.0),
        "power_from_four_cycle_ensemble_mean":
            abs(power / ENSEMBLE_POWER_W_M3 - 1.0),
        "density_from_four_cycle_ensemble_mean":
            abs(density / ENSEMBLE_DENSITY_M3 - 1.0),
        "power_per_electron_ratio_from_four_cycle_ensemble_mean":
            abs(power_per_electron_ratio /
                ENSEMBLE_POWER_PER_ELECTRON_RATIO - 1.0),
    }
    thresholds = {
        "phase_binned_to_exact_power": 0.02,
        "power_from_four_cycle_ensemble_mean": 0.03,
        "density_from_four_cycle_ensemble_mean": 0.01,
        "power_per_electron_ratio_from_four_cycle_ensemble_mean": 0.03,
    }
    gates = {name: differences[name] <= threshold
             for name, threshold in thresholds.items()}
    return {"relative_differences": differences, "thresholds": thresholds,
            "gates": gates, "all_gates_passed": all(gates.values())}


def analyze(root: Path, reference: Path) -> dict[str, object]:
    report_path = root / "branch-report.json"
    if sha256(report_path) != BRANCH_REPORT_SHA256:
        raise ValueError("long-window branch report differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if (report.get("rule_sha256") != RULE_SHA256 or
            report.get("all_gates_passed") is not True):
        raise ValueError("long-window execution contract did not pass")
    output = root / "measurement" / "output"
    for name, expected in report["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"long-window output differs: {name}")
    metadata = json.loads(
        (output / "spatial_average_metadata.json").read_text(encoding="utf-8"))
    if (metadata.get("sampling_order") != "pre_collision" or
            metadata.get("phase_bins") != PHASES or
            metadata.get("samples") != 48000 or
            set(metadata.get("phase_bin_samples", [])) != {240} or
            metadata.get("complete") is not True):
        raise ValueError("long-window sampling protocol differs")

    fields = rows(output / "spatial_phase_fields.csv")
    moments = [row for row in rows(output / "spatial_phase_moments.csv")
               if row["species"] == "electrons"]
    if len(fields) != PHASES * NODES or len(moments) != len(fields):
        raise ValueError("long-window phase-space shape differs")
    candidate = {
        "electron_density": [float(row["number_density_mean_m-3"])
                             for row in moments],
        "electric_field": [float(row["electric_field_mean_V_m"])
                           for row in fields],
    }
    candidate["electron_current_density"] = [
        -ELEMENTARY_CHARGE_C * density * float(row["mean_velocity_x"])
        for density, row in zip(candidate["electron_density"], moments)
    ]
    candidate["electron_ohmic_power_density"] = [
        current * field for current, field in zip(
            candidate["electron_current_density"], candidate["electric_field"])
    ]
    reference_values = {}
    for name in candidate:
        filename, expected = REFERENCE_FILES[name]
        path = reference / filename
        if sha256(path) != expected:
            raise ValueError(f"locked eduPIC reference differs: {filename}")
        reference_values[name] = flatten_phase_major(read_matrix(path))

    density = spatial_phase_average(
        candidate["electron_density"], PHASES, NODES)
    power = spatial_phase_average(
        candidate["electron_ohmic_power_density"], PHASES, NODES)
    reference_density = spatial_phase_average(
        reference_values["electron_density"], PHASES, NODES)
    reference_power = spatial_phase_average(
        reference_values["electron_ohmic_power_density"], PHASES, NODES)
    power_ratio = power / reference_power
    density_ratio = density / reference_density
    power_per_electron_ratio = power_ratio / density_ratio
    energy = json.loads((output / "energy-budget.json").read_text())
    exact_power = float(
        energy["electric_power_W_m-2"]["electric_work_electrons_J_m-2"]
    ) / LENGTH_M
    acceptance = gate_summary(
        power, exact_power, density, power_per_electron_ratio)
    if acceptance["all_gates_passed"] is not True:
        raise ValueError("long-window prospective acceptance failed")

    comparisons = {}
    reduced = {}
    for name, values in candidate.items():
        comparisons[name] = metrics(values, reference_values[name])
        candidate_phase, candidate_space = reduced_profiles(values)
        reference_phase, reference_space = reduced_profiles(
            reference_values[name])
        reduced[name] = {
            "spatially_integrated_phase_profile": metrics(
                candidate_phase, reference_phase),
            "cycle_average_spatial_profile": metrics(
                candidate_space, reference_space),
        }
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "long_window_matched_heating_validation",
        "rule_sha256": RULE_SHA256,
        "branch_report_sha256": BRANCH_REPORT_SHA256,
        "sampling": {
            "measurement_cycles": 12,
            "samples": 48000,
            "phase_bins": PHASES,
            "samples_per_phase_bin": 240,
            "sampling_order": "pre_collision",
        },
        "prospective_acceptance": acceptance,
        "volume_phase_averages": {
            "candidate_power_density_W_m-3": power,
            "candidate_exact_electric_work_W_m-3": exact_power,
            "candidate_electron_density_m-3": density,
            "reference_power_density_W_m-3": reference_power,
            "reference_electron_density_m-3": reference_density,
            "candidate_to_reference_power_density_ratio": power_ratio,
            "candidate_to_reference_electron_density_ratio": density_ratio,
            "candidate_to_reference_power_per_electron_ratio":
                power_per_electron_ratio,
        },
        "phase_space_comparisons": comparisons,
        "reduced_profile_comparisons": reduced,
        "resource_recovery_boundary": report["resources"],
        "finding": (
            "The twelve-cycle result passes every predeclared consistency "
            "gate and retains the density, power-per-electron, and spatial "
            "current mismatches seen in the four-cycle ensemble."),
        "claim_boundary": (
            "This tests window-length and continuation-seed sensitivity from "
            "one initial particle realization; it is not experimental or full "
            "initial-condition validation."),
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
