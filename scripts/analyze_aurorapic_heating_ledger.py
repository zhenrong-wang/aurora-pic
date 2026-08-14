#!/usr/bin/env python3
"""Compare AuroraPIC phase-binned J.E with exact electric work and eduPIC."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_aurorapic_mesh_refinement import phase_space_average
from compare_aurorapic_edupic_measurement_pilot import rows
from compare_edupic_phase_space import (
    REFERENCE_FILES, flatten_phase_major, read_matrix, spatial_phase_average,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


BRANCH_REPORT_SHA256 = (
    "c81b2590ff5b1998f0312f2d9b9979f94987c5b7c78aae5a159272a3b55f4fab")
ELEMENTARY_CHARGE_C = 1.60217662e-19
LENGTH_M = 0.025


def summarize(candidate_density: float, candidate_binned_power: float,
              candidate_exact_power: float, reference_density: float,
              reference_binned_power: float) -> dict[str, float | None]:
    values = (candidate_density, candidate_exact_power, reference_density,
              reference_binned_power)
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError("heating-ledger inputs must be positive and finite")
    if not math.isfinite(candidate_binned_power):
        raise ValueError("phase-binned candidate power must be finite")
    candidate_binned_per_particle = candidate_binned_power / candidate_density
    candidate_exact_per_particle = candidate_exact_power / candidate_density
    reference_per_particle = reference_binned_power / reference_density
    binned_ratio = candidate_binned_per_particle / reference_per_particle
    exact_ratio = candidate_exact_per_particle / reference_per_particle
    apparent_deficit = 1.0 - binned_ratio
    deficit_reduction = None
    if not math.isclose(apparent_deficit, 0.0, abs_tol=1.0e-15):
        deficit_reduction = (exact_ratio - binned_ratio) / apparent_deficit
    return {
        "candidate_phase_binned_j_dot_e_W_m-3": candidate_binned_power,
        "candidate_exact_electric_work_W_m-3": candidate_exact_power,
        "candidate_phase_binned_to_exact_power_ratio":
            candidate_binned_power / candidate_exact_power,
        "candidate_phase_binned_power_relative_deficit":
            1.0 - candidate_binned_power / candidate_exact_power,
        "candidate_volume_phase_average_density_m-3": candidate_density,
        "reference_phase_binned_j_dot_e_W_m-3": reference_binned_power,
        "reference_volume_phase_average_density_m-3": reference_density,
        "candidate_phase_binned_power_per_particle_W":
            candidate_binned_per_particle,
        "candidate_exact_power_per_particle_W": candidate_exact_per_particle,
        "reference_phase_binned_power_per_particle_W": reference_per_particle,
        "candidate_to_reference_phase_binned_power_per_particle_ratio":
            binned_ratio,
        "candidate_exact_to_reference_binned_power_per_particle_ratio":
            exact_ratio,
        "apparent_cross_code_deficit_fraction_phase_binned": apparent_deficit,
        "alternate_cross_code_deficit_fraction_exact_candidate": 1.0 - exact_ratio,
        "apparent_deficit_reduction_fraction": deficit_reduction,
    }


def analyze(output: Path, reference: Path) -> dict[str, object]:
    branch_path = output.parent.parent / "branch-report.json"
    if sha256(branch_path) != BRANCH_REPORT_SHA256:
        raise ValueError("region-matched candidate branch differs")
    branch = json.loads(branch_path.read_text(encoding="utf-8"))
    for name, expected in branch["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"candidate output differs: {name}")
    energy_path = output / "energy-budget.json"
    if sha256(energy_path) != branch["energy_analysis_sha256"]:
        raise ValueError("candidate energy analysis differs")
    energy = json.loads(energy_path.read_text(encoding="utf-8"))
    if energy.get("passes") is not True:
        raise ValueError("candidate energy ledger does not close")
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
        raise ValueError("candidate phase-space shapes differ")
    x = [float(row["x_m"]) for row in fields[:nodes]]
    density = [float(row["number_density_mean_m-3"]) for row in moments]
    current = [-ELEMENTARY_CHARGE_C * number * float(row["mean_velocity_x"])
               for number, row in zip(density, moments)]
    binned_power = [value * float(row["electric_field_mean_V_m"])
                    for value, row in zip(current, fields)]
    candidate_density = phase_space_average(x, density)
    candidate_binned_power = phase_space_average(x, binned_power)
    candidate_exact_power = float(
        energy["electric_power_W_m-2"]["electric_work_electrons_J_m-2"]
    ) / LENGTH_M

    reference_density = spatial_phase_average(
        flatten_phase_major(read_matrix(
            reference / REFERENCE_FILES["electron_density"][0])),
        phases=200, nodes=400)
    reference_power = spatial_phase_average(
        flatten_phase_major(read_matrix(
            reference / REFERENCE_FILES["electron_ohmic_power_density"][0])),
        phases=200, nodes=400)
    metrics = summarize(candidate_density, candidate_binned_power,
                        candidate_exact_power, reference_density,
                        reference_power)
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "post_measurement_electron_heating_ledger_localization",
        "candidate_branch_report_sha256": BRANCH_REPORT_SHA256,
        "candidate_energy_analysis_sha256": sha256(energy_path),
        "candidate_phase_bins": phases,
        "candidate_moment_sampling_order": "post_collision",
        "reference_phase_bins": 200,
        "reference_moment_sampling_order": "pre_collision",
        "metrics": metrics,
        "finding": (
            "The 16-bin post-collision candidate J.E estimator understates "
            "AuroraPIC's exact electron electric work. Correcting candidate "
            "power with the exact ledger reduces, but does not eliminate, the "
            "cross-code power-per-electron deficit."),
        "next_discriminator": (
            "Generate 200-bin pre-collision AuroraPIC moments to match eduPIC "
            "sampling order and phase resolution, and compare their integrated "
            "J.E against both the exact ledger and the reference matrix."),
        "claim_boundary": (
            "The exact AuroraPIC ledger and eduPIC's phase-binned J.E are "
            "different estimators. Their alternate ratio localizes diagnostic "
            "bias but is not a matched cross-code acceptance result."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_output", type=Path)
    parser.add_argument("reference_raw_data", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.candidate_output.resolve(),
                     args.reference_raw_data.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
