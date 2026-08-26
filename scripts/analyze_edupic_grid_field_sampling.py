#!/usr/bin/env python3
"""Separate grid-field amplitude from near-threshold particle sampling."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from compare_edupic_phase_space import read_matrix


METRICS = (
    "grid_mean_squared_field_V2_m2",
    "particle_sampling_factor",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def boundary_value(values: list[float], x: float, length: float) -> float:
    scaled = x / length * (len(values) - 1)
    left = min(int(math.floor(scaled)), len(values) - 2)
    fraction = scaled - left
    return values[left] * (1.0 - fraction) + values[left + 1] * fraction


def spatial_mean_square(values: list[float], length: float,
                        lower: float, upper: float) -> float:
    if len(values) < 2 or not 0.0 <= lower < upper <= length:
        raise ValueError("invalid spatial mean-square contract")
    dx = length / (len(values) - 1)
    points = [(lower, boundary_value(values, lower, length))]
    points.extend((node * dx, value) for node, value in enumerate(values)
                  if lower < node * dx < upper)
    points.append((upper, boundary_value(values, upper, length)))
    integral = math.fsum(
        0.5 * (left_value ** 2 + right_value ** 2) *
        (right_x - left_x)
        for (left_x, left_value), (right_x, right_value)
        in zip(points, points[1:]))
    return integral / (upper - lower)


def phase_space_mean_square(matrix: list[list[float]], length: float,
                            lower_x: float, upper_x: float,
                            lower_phase: int, upper_phase: int) -> float:
    if (len(matrix) < 2 or lower_phase < 0 or lower_phase >= upper_phase or
            any(len(row) < upper_phase for row in matrix)):
        raise ValueError("field matrix is incomplete")
    profiles = [
        [matrix[node][phase] for node in range(len(matrix))]
        for phase in range(lower_phase, upper_phase)
    ]
    return math.fsum(
        spatial_mean_square(profile, length, lower_x, upper_x)
        for profile in profiles) / len(profiles)


def aurora_matrix(path: Path, nodes: int,
                  phases: int, length: float) -> list[list[float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != nodes * phases:
        raise ValueError("AuroraPIC phase-field table has the wrong shape")
    matrix = [[0.0] * phases for _ in range(nodes)]
    sample_counts = set()
    for index, row in enumerate(rows):
        phase, node = divmod(index, nodes)
        if int(row["phase_bin"]) != phase or int(row["node"]) != node:
            raise ValueError("AuroraPIC phase-field ordering differs")
        if abs(float(row["phase_fraction"]) -
               (phase + 0.5) / phases) > 1e-12:
            raise ValueError("AuroraPIC phase centers differ")
        if abs(float(row["x_m"]) - node * length / (nodes - 1)) > 1e-12:
            raise ValueError("AuroraPIC field coordinates differ")
        value = float(row["electric_field_mean_V_m"])
        if not math.isfinite(value):
            raise ValueError("AuroraPIC field is non-finite")
        matrix[node][phase] = value
        sample_counts.add(int(row["samples"]))
    if len(sample_counts) != 1 or min(sample_counts) <= 0:
        raise ValueError("AuroraPIC phase sample counts differ")
    return matrix


def relative_ranges(members: list[dict[str, object]]) -> dict[str, float]:
    result = {}
    for metric in METRICS:
        values = [float(member[metric]) for member in members]
        mean = sum(values) / len(values)
        result[metric] = (max(values) - min(values)) / max(abs(mean), 1e-300)
    return result


def means(members: list[dict[str, object]]) -> dict[str, float]:
    return {
        metric: sum(float(member[metric]) for member in members) / len(members)
        for metric in METRICS
    }


def analyze(rule_path: Path, native_result_path: Path,
            aurora_reports: list[Path], aurora_fields: list[Path],
            native_fields: list[Path]) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    native_result = json.loads(
        native_result_path.read_text(encoding="utf-8"))
    grid = rule["grid_contract"]
    locked = rule["locked_inputs"]
    length = float(grid["length_m"])
    nodes = int(grid["nodes"])
    phases = int(grid["phase_bins"])
    lower_x, upper_x = [float(value) * length
                        for value in grid["critical_x_over_L"]]
    lower_phase, upper_phase = map(int, grid["critical_phase_bins"])
    factor = (2.0 * float(grid["electron_mass_kg"]) /
              (float(grid["elementary_charge_C"]) *
               float(grid["electron_timestep_s"]) ** 2))
    reports = [json.loads(path.read_text(encoding="utf-8"))
               for path in aurora_reports]
    native_by_seed = {entry["seed"]: entry
                      for entry in native_result["members"]}
    aurora_members = []
    for report, field_path, expected in zip(
            reports, aurora_fields, locked["aurorapic_members"], strict=True):
        matrix = aurora_matrix(field_path, nodes, phases, length)
        grid_e2 = phase_space_mean_square(
            matrix, length, lower_x, upper_x, lower_phase, upper_phase)
        particle_e2 = factor * float(
            report["critical_scope"]["mean_quadratic_work_eV"])
        aurora_members.append({
            "id": expected["id"],
            "grid_mean_squared_field_V2_m2": grid_e2,
            "grid_rms_field_V_m": math.sqrt(grid_e2),
            "particle_sampled_mean_squared_field_V2_m2": particle_e2,
            "particle_sampling_factor": particle_e2 / grid_e2,
        })
    native_members = []
    for field_path, expected in zip(
            native_fields, locked["native_members"], strict=True):
        matrix = read_matrix(field_path, nodes, phases)
        grid_e2 = phase_space_mean_square(
            matrix, length, lower_x, upper_x, lower_phase, upper_phase)
        source = native_by_seed[expected["seed"]]
        particle_e2 = factor * float(source["mean_quadratic_work_eV"])
        native_members.append({
            "seed": expected["seed"],
            "grid_mean_squared_field_V2_m2": grid_e2,
            "grid_rms_field_V_m": math.sqrt(grid_e2),
            "particle_sampled_mean_squared_field_V2_m2": particle_e2,
            "particle_sampling_factor": particle_e2 / grid_e2,
        })
    aurora_ranges = relative_ranges(aurora_members)
    native_ranges = relative_ranges(native_members)
    aurora_means = means(aurora_members)
    native_means = means(native_members)
    ratios = [
        {metric: float(member[metric]) / native_means[metric]
         for metric in METRICS}
        for member in aurora_members
    ]
    mean_ratios = {
        metric: aurora_means[metric] / native_means[metric]
        for metric in METRICS
    }
    hashes_linked = (
        sha256(native_result_path) == locked["native_result_sha256"]
        and all(
            sha256(report_path) == expected["runner_report_sha256"]
            and sha256(field_path) == expected["spatial_phase_fields_sha256"]
            and report.get("all_gates_passed") is True
            and report.get("state_id") == expected["id"]
            for report_path, field_path, report, expected in zip(
                aurora_reports, aurora_fields, reports,
                locked["aurorapic_members"], strict=True))
        and all(
            sha256(path) == expected["efield_xt_sha256"]
            for path, expected in zip(
                native_fields, locked["native_members"], strict=True))
    )
    gates = {
        "locked_hashes_and_reports": hashes_linked,
        "member_counts": len(aurora_members) == 2 and
            len(native_members) == 3,
        "finite_positive_metrics": all(
            math.isfinite(float(member[metric])) and
            float(member[metric]) > 0.0
            for member in [*aurora_members, *native_members]
            for metric in METRICS),
        "aurorapic_repeatability": all(
            value <= 0.10 for value in aurora_ranges.values()),
        "native_repeatability": all(
            value <= 0.15 for value in native_ranges.values()),
    }
    all_gates = all(gates.values())
    grid_deficit = all(
        ratio["grid_mean_squared_field_V2_m2"] <= 0.90 for ratio in ratios)
    sampling_deficit = all(
        ratio["particle_sampling_factor"] <= 0.90 for ratio in ratios)
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "aurorapic_edupic_grid_field_vs_particle_sampling_result",
        "rule_sha256": sha256(rule_path),
        "gates": gates,
        "all_hash_shape_repeatability_and_metric_gates_passed": all_gates,
        "aurorapic_members": aurora_members,
        "native_edupic_members": native_members,
        "aurorapic_ensemble_means": aurora_means,
        "native_edupic_ensemble_means": native_means,
        "aurorapic_relative_ranges": aurora_ranges,
        "native_edupic_relative_ranges": native_ranges,
        "aurorapic_member_to_native_mean_ratios": ratios,
        "aurorapic_ensemble_mean_to_native_mean_ratios": mean_ratios,
        "prospective_decision_outcome": {
            "interpretation_allowed": all_gates,
            "grid_field_deficit_supported": all_gates and grid_deficit,
            "differential_particle_sampling_deficit_supported":
                all_gates and sampling_deficit,
            "result": (
                "grid_field_deficit_supported_particle_sampling_deficit_not_supported"
                if all_gates and grid_deficit and not sampling_deficit
                else "mixed_or_intermediate_grid_field_sampling_result"
            ),
        },
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "native_result_sha256": sha256(native_result_path),
            "aurorapic_field_sha256": [sha256(path) for path in aurora_fields],
            "native_field_sha256": [sha256(path) for path in native_fields],
        },
        "interpretation_note": (
            "Grid fields use squared phase-mean fields in both codes. The "
            "particle metric is an exact kick-derived mean E-squared. Their "
            "ratio localizes conditional spatial sampling but does not identify "
            "the first temporal source of any grid-field divergence."
        ),
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("native_result", type=Path)
    parser.add_argument("aurorapic_reports", nargs=2, type=Path)
    parser.add_argument("aurorapic_fields", nargs=2, type=Path)
    parser.add_argument("native_fields", nargs=3, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = analyze(args.rule, args.native_result, args.aurorapic_reports,
                     args.aurorapic_fields, args.native_fields)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result[
        "all_hash_shape_repeatability_and_metric_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
