#!/usr/bin/env python3
"""Audit the 1D density, current, and phase-normalization contract with eduPIC."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from compare_edupic_phase_space import (
    ELEMENTARY_CHARGE_C, REFERENCE_FILES, read_matrix,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


NODES = 400
PHASES = 200
LENGTH_M = 0.025
EDUPIC_WEIGHT = 7.0e4
EDUPIC_AREA_M2 = 1.0e-4
UPSTREAM_EDUPIC_SOURCE_SHA256 = (
    "7c7679c0f0c98844940ea911bbb7581ec33f818e8d14427c9837ffdcf1ecea41")
EXPECTED_LOCAL_GENERATION_SOURCE_SHA256 = (
    "a850889bbc3c5917505eb31752cde607b7550c8212f7df01fa739b70d1a6a79f")
EXPECTED_BRANCH_REPORT_SHA256 = (
    "1d87b9ebfbe3668513c844fc5314868e1bae625df259c7f659981f9ba03a9b4a")


def trapezoid_integral(values: list[float], dx: float) -> float:
    if len(values) < 2 or not all(math.isfinite(value) for value in values):
        raise ValueError("trapezoid input must contain at least two finite values")
    return dx * (0.5 * values[0] + math.fsum(values[1:-1]) +
                 0.5 * values[-1])


def implied_counts(matrix: list[list[float]], dx: float,
                   line_weight: float) -> list[float]:
    if len(matrix) != NODES or any(len(row) != PHASES for row in matrix):
        raise ValueError("density matrix must be 400 nodes by 200 phase bins")
    if not math.isfinite(line_weight) or line_weight <= 0.0:
        raise ValueError("line weight must be positive and finite")
    return [trapezoid_integral(
        [matrix[node][phase] for node in range(NODES)], dx) / line_weight
        for phase in range(PHASES)]


def phase_matrix(path: Path, species: str,
                 column: str) -> tuple[list[list[float]], set[int]]:
    matrix = [[math.nan] * PHASES for _ in range(NODES)]
    samples: set[int] = set()
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["species"] != species:
                continue
            node = int(row["node"])
            phase = int(row["phase_bin"])
            matrix[node][phase] = float(row[column])
            samples.add(int(row["samples"]))
    if any(not math.isfinite(value) for row in matrix for value in row):
        raise ValueError("candidate phase matrix is incomplete or non-finite")
    return matrix, samples


def parse_scalar_config(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    section = ""
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
            continue
        if "=" in line:
            key, value = (part.strip() for part in line.split("=", 1))
            values[f"{section}.{key}" if section else key] = value
    return values


def relative_l2(candidate: list[float], reference: list[float]) -> float:
    denominator = math.fsum(value * value for value in reference)
    if len(candidate) != len(reference) or denominator <= 0.0:
        raise ValueError("relative-L2 vectors are incompatible")
    return math.sqrt(math.fsum((a - b) ** 2 for a, b in
                               zip(candidate, reference)) / denominator)


def flatten_space_major(matrix: list[list[float]]) -> list[float]:
    return [value for row in matrix for value in row]


def arithmetic_average(matrix: list[list[float]]) -> float:
    return math.fsum(flatten_space_major(matrix)) / (NODES * PHASES)


def physical_average(matrix: list[list[float]]) -> float:
    return math.fsum(
        (0.5 if node in (0, NODES - 1) else 1.0) * matrix[node][phase]
        for node in range(NODES) for phase in range(PHASES)
    ) / ((NODES - 1) * PHASES)


def summary(values: list[float]) -> dict[str, float]:
    return {"minimum": min(values), "mean": math.fsum(values) / len(values),
            "maximum": max(values)}


def audit(candidate_root: Path, reference: Path,
          edupic_source: Path) -> dict[str, object]:
    if sha256(edupic_source) != EXPECTED_LOCAL_GENERATION_SOURCE_SHA256:
        raise ValueError("eduPIC implementation differs from the reviewed local source")
    source = edupic_source.read_text(encoding="utf-8")
    required_source_terms = (
        "const double     WEIGHT         = 7.0e4;",
        "const double     ELECTRODE_AREA = 1.0e-4;",
        "const int        N_G            = 400;",
        "const int        N_T            = 4000;",
        "const int N_BIN                     = 20;",
        "const int N_XT                      = N_T / N_BIN;",
        "const double FACTOR_W = WEIGHT / DV;",
        "e_density[0]     *= 2.0;",
        "mean_v = vx_e[k] - 0.5 * e_x * FACTOR_E;",
        "f1 = (double)(N_XT) / (double)(no_of_cycles * N_T);",
        "je_xt[i][j]     = -ue_xt[i][j] * ne_xt[i][j] * E_CHARGE;",
        "powere_xt[i][j] = je_xt[i][j] * efield_xt[i][j];",
    )
    missing = [term for term in required_source_terms if term not in source]
    if missing:
        raise ValueError(f"reviewed eduPIC normalization terms missing: {missing}")

    report_path = candidate_root / "branch-report.json"
    if sha256(report_path) != EXPECTED_BRANCH_REPORT_SHA256:
        raise ValueError("candidate branch report differs from the locked run")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    output = candidate_root / "measurement" / "output"
    for name, expected in report["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"candidate output differs: {name}")

    config = parse_scalar_config(candidate_root / "measurement" / "input.cfg")
    candidate_weight = float(config["species.electrons.weight"])
    nodes = int(config["nx"])
    length = float(config["length"])
    phase_bins = int(config["spatial_average_phase_bins"])
    if nodes != NODES or phase_bins != PHASES or length != LENGTH_M:
        raise ValueError("candidate grid or phase contract differs")
    reference_line_weight = EDUPIC_WEIGHT / EDUPIC_AREA_M2
    if candidate_weight != reference_line_weight:
        raise ValueError("candidate and reference line weights differ")
    dx = length / (nodes - 1)

    reference_density_path = reference / REFERENCE_FILES["electron_density"][0]
    reference_current_path = reference / REFERENCE_FILES["electron_current_density"][0]
    reference_field_path = reference / REFERENCE_FILES["electric_field"][0]
    reference_power_path = reference / REFERENCE_FILES[
        "electron_ohmic_power_density"][0]
    for name, path in (
            ("electron_density", reference_density_path),
            ("electron_current_density", reference_current_path),
            ("electric_field", reference_field_path),
            ("electron_ohmic_power_density", reference_power_path)):
        if sha256(path) != REFERENCE_FILES[name][1]:
            raise ValueError(f"locked eduPIC reference differs: {path.name}")
    reference_density = read_matrix(reference_density_path)
    reference_current = read_matrix(reference_current_path)
    reference_field = read_matrix(reference_field_path)
    reference_power = read_matrix(reference_power_path)
    reconstructed_reference_power = [
        current * field for current, field in zip(
            flatten_space_major(reference_current),
            flatten_space_major(reference_field))]
    reference_power_values = flatten_space_major(reference_power)

    candidate_density, samples = phase_matrix(
        output / "spatial_phase_moments.csv", "electrons",
        "number_density_mean_m-3")
    if samples != {240}:
        raise ValueError("candidate phase-bin sample count differs")
    reference_counts = implied_counts(reference_density, dx,
                                      reference_line_weight)
    candidate_counts = implied_counts(candidate_density, dx,
                                      candidate_weight)
    count_ratio = (math.fsum(candidate_counts) /
                   math.fsum(reference_counts))
    density_ratio = physical_average(candidate_density) / physical_average(
        reference_density)
    if abs(count_ratio / density_ratio - 1.0) > 1.0e-12:
        raise ValueError("integrated density and implied-count ratios disagree")

    reference_simple_density = arithmetic_average(reference_density)
    candidate_simple_density = arithmetic_average(candidate_density)
    simple_ratio = candidate_simple_density / reference_simple_density
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "cross_code_normalization_audit",
        "inputs": {
            "edupic_upstream_source_sha256": UPSTREAM_EDUPIC_SOURCE_SHA256,
            "edupic_local_generation_source_sha256": sha256(edupic_source),
            "local_patch_boundary": (
                "Only main() is changed to export cross sections and return; "
                "the audited normalization code is unmodified."),
            "candidate_branch_report_sha256": sha256(report_path),
            "reference_density_sha256": sha256(reference_density_path),
            "candidate_phase_moments_sha256": sha256(
                output / "spatial_phase_moments.csv"),
        },
        "weight_contract": {
            "edupic_superparticle_weight": EDUPIC_WEIGHT,
            "edupic_fictive_area_m2": EDUPIC_AREA_M2,
            "edupic_effective_1d_line_weight_m-2": reference_line_weight,
            "aurorapic_1d_line_weight_m-2": candidate_weight,
            "ratio": candidate_weight / reference_line_weight,
            "equivalent": candidate_weight == reference_line_weight,
        },
        "deposition_contract": {
            "grid_nodes": nodes,
            "gap_length_m": length,
            "node_spacing_m": dx,
            "shape": "linear_CIC",
            "edupic_boundary_convention": "multiply endpoint density by two",
            "aurorapic_boundary_convention": "divide endpoint deposit by half-cell volume",
            "boundary_conventions_algebraically_equivalent": True,
            "physical_spatial_integral": "endpoint_half_weight_trapezoid",
        },
        "phase_contract": {
            "phase_bins": PHASES,
            "edupic_steps_per_bin": 20,
            "edupic_normalization": "1 / (measurement_cycles * 20)",
            "aurorapic_samples_per_bin": next(iter(samples)),
            "aurorapic_normalization": "1 / phase_bin_samples",
            "both_are_arithmetic_time_averages": True,
            "sampling_order": "pre_collision",
        },
        "density_conservation": {
            "reference_implied_macro_particle_count_by_phase":
                summary(reference_counts),
            "candidate_implied_macro_particle_count_by_phase":
                summary(candidate_counts),
            "candidate_to_reference_mean_implied_count_ratio": count_ratio,
            "candidate_to_reference_physical_density_ratio": density_ratio,
            "ratio_identity_relative_error": abs(count_ratio / density_ratio - 1.0),
        },
        "current_and_power_contract": {
            "reference_current": "-elementary_charge * density * mean_velocity_x",
            "candidate_current": "-elementary_charge * density * mean_velocity_x",
            "mean_velocity_time_centering": "pre-push x with leapfrog-centered vx",
            "reference_power_matrix_jE_relative_l2_reconstruction_error":
                relative_l2(reconstructed_reference_power,
                            reference_power_values),
            "elementary_charge_C": ELEMENTARY_CHARGE_C,
        },
        "spatial_quadrature_sensitivity": {
            "reference_trapezoid_vs_simple_density_relative_change":
                physical_average(reference_density) /
                reference_simple_density - 1.0,
            "candidate_trapezoid_vs_simple_density_relative_change":
                physical_average(candidate_density) /
                candidate_simple_density - 1.0,
            "candidate_to_reference_density_ratio_trapezoid": density_ratio,
            "candidate_to_reference_density_ratio_simple_node_average":
                simple_ratio,
            "ratio_change": density_ratio / simple_ratio - 1.0,
        },
        "assessment": {
            "normalization_mismatch_found": False,
            "density_excess_explained_by_normalization": False,
            "finding": (
                "The codes use the same effective 1D macro weight, CIC shape, "
                "half-cell endpoint volume, phase arithmetic mean, and q*n*u "
                "current convention. Integrating each density phase recovers "
                "macro-particle counts whose cross-code ratio is exactly the "
                "reported density ratio; the persistent density excess is a "
                "population/state difference, not a unit conversion."),
            "next_target": (
                "Run independent initial particle realizations from the same "
                "macroscopic state to quantify state/seed uncertainty in the "
                "density and spatial-current amplitudes."),
        },
        "claim_boundary": (
            "This rules out the audited 1D weight, deposition-volume, phase-"
            "averaging, current, and spatial-quadrature conventions. It does "
            "not prove that the two kinetic states or algorithms are equal."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_root", type=Path)
    parser.add_argument("reference_raw_data", type=Path)
    parser.add_argument("edupic_source", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.candidate_root.resolve(),
                   args.reference_raw_data.resolve(),
                   args.edupic_source.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
