#!/usr/bin/env python3
"""Descriptively compare the fresh AuroraPIC pilot with locked eduPIC tables."""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
from pathlib import Path

from run_aurorapic_edupic_pilot import atomic_json, sha256


REFERENCE_HASHES = {
    "density.csv": "4c860617882a551c502934cf38bcf22a6432565f2a8da277a35de216868400a9",
    "eepf.csv": "8ad8222107c7e61ea9d1a2d36b0550bcd69ded55dbf5fe0bcfd10810dfc56260",
    "ifed.csv": "98c7f40ee0f0ab12038c51f2dafc0b9f9e17781c7264d5399018f86d88dd2b3e",
}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        result = list(csv.DictReader(stream))
    if not result:
        raise ValueError(f"empty comparison table: {path}")
    return result


def interpolate(x: list[float], y: list[float], target: float) -> float:
    if target <= x[0]:
        return y[0]
    if target >= x[-1]:
        return y[-1]
    left = bisect.bisect_right(x, target) - 1
    fraction = (target - x[left]) / (x[left + 1] - x[left])
    return y[left] * (1.0 - fraction) + y[left + 1] * fraction


def relative_l2(candidate: list[float], reference: list[float]) -> float:
    return math.sqrt(
        math.fsum((a - b) ** 2 for a, b in zip(candidate, reference)) /
        math.fsum(value * value for value in reference))


def trapezoid(x: list[float], y: list[float]) -> float:
    return math.fsum(0.5 * (a + b) * (right - left)
                     for left, right, a, b in zip(x, x[1:], y, y[1:]))


def total_variation(a: list[tuple[float, float, float]],
                    b: list[tuple[float, float, float]]) -> float:
    edges = sorted({edge for low, high, _ in a + b for edge in (low, high)})

    def density(table: list[tuple[float, float, float]], value: float) -> float:
        for low, high, mass in table:
            if low <= value < high:
                return mass / (high - low)
        return 0.0

    return 0.5 * math.fsum(
        abs(density(a, 0.5 * (low + high)) -
            density(b, 0.5 * (low + high))) * (high - low)
        for low, high in zip(edges, edges[1:]))


def distribution_mean(table: list[tuple[float, float, float]]) -> float:
    normalization = math.fsum(mass for _, _, mass in table)
    return math.fsum(0.5 * (low + high) * mass
                     for low, high, mass in table) / normalization


def relative_difference(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b))


def analyze(candidate: Path, reference: Path) -> dict[str, object]:
    report_path = candidate.parent / "measurement-report.json"
    pilot = json.loads(report_path.read_text(encoding="utf-8"))
    if (pilot.get("scope") != "fresh_window_aurorapic_measurement_pilot" or
            pilot.get("all_gates_passed") is not True or
            pilot.get("window", {}).get("equilibration_statistics_excluded")
            is not True):
        raise ValueError("candidate is not a passing fresh-window pilot")
    for name, expected in REFERENCE_HASHES.items():
        if sha256(reference / name) != expected:
            raise ValueError(f"locked eduPIC {name} differs")

    density = rows(candidate / "spatial_average.csv")
    reference_density = rows(reference / "density.csv")
    species = {}
    for name in ("electrons", "ions"):
        selected = [row for row in density if row["species"] == name]
        species[name] = (
            [float(row["x_m"]) for row in selected],
            [float(row["number_density_mean_m-3"]) for row in selected],
        )
    x = [float(row["x_m"]) for row in reference_density]
    reference_profiles = {
        "electrons": [float(row["electron_density_duration_mean_m3"])
                      for row in reference_density],
        "ions": [float(row["ion_density_duration_mean_m3"])
                 for row in reference_density],
    }
    density_metrics = {}
    for name in ("electrons", "ions"):
        profile = [interpolate(*species[name], coordinate) for coordinate in x]
        candidate_integral = trapezoid(x, profile)
        reference_integral = trapezoid(x, reference_profiles[name])
        density_metrics[name] = {
            "profile_relative_l2": relative_l2(profile, reference_profiles[name]),
            "candidate_line_integrated_density_m-2": candidate_integral,
            "reference_line_integrated_density_m-2": reference_integral,
            "candidate_to_reference_line_integrated_ratio":
                candidate_integral / reference_integral,
        }

    wall = rows(candidate / "wall_impact_spectrum.csv")
    reference_ifed = rows(reference / "ifed.csv")

    def candidate_ifed(side: str) -> list[tuple[float, float, float]]:
        selected = [row for row in wall
                    if row["species"] == "ions" and row["electrode"] == side]
        width = (float(selected[1]["impact_energy_eV"]) -
                 float(selected[0]["impact_energy_eV"]))
        return [(float(row["impact_energy_eV"]) - 0.5 * width,
                 float(row["impact_energy_eV"]) + 0.5 * width,
                 float(row["probability_density"]) * width)
                for row in selected]

    def locked_ifed(field: str) -> list[tuple[float, float, float]]:
        return [(float(row["energy_ev"]) - 0.5,
                 float(row["energy_ev"]) + 0.5, float(row[field]))
                for row in reference_ifed]

    ifed_metrics = {}
    for side, electrode, field in (
        ("left", "powered", "equal_time_block_mixture_powered_ev-1"),
        ("right", "grounded", "equal_time_block_mixture_grounded_ev-1"),
    ):
        candidate_table = candidate_ifed(side)
        reference_table = locked_ifed(field)
        candidate_mean = distribution_mean(candidate_table)
        reference_mean = distribution_mean(reference_table)
        ifed_metrics[electrode] = {
            "total_variation": total_variation(candidate_table, reference_table),
            "candidate_mean_energy_eV": candidate_mean,
            "reference_mean_energy_eV": reference_mean,
            "mean_energy_relative_difference":
                relative_difference(candidate_mean, reference_mean),
        }

    histogram = rows(candidate / "phase_eedf.csv")
    moments = rows(candidate / "phase_eedf_moments.csv")
    total_observations = math.fsum(
        float(row["represented_observations"]) for row in moments)
    counts: dict[int, float] = {}
    centers: dict[int, float] = {}
    for row in histogram:
        energy_bin = int(row["energy_bin"])
        counts[energy_bin] = counts.get(energy_bin, 0.0) + float(
            row["represented_count"])
        centers[energy_bin] = float(row["energy_eV"])
    indices = sorted(counts)
    width = centers[indices[1]] - centers[indices[0]]
    candidate_eedf = [
        (centers[index] - 0.5 * width, centers[index] + 0.5 * width,
         counts[index] / total_observations) for index in indices]
    reference_eedf = [
        (float(row["energy_ev"]) - 0.025,
         float(row["energy_ev"]) + 0.025,
         float(row["equal_time_block_mixture_eepf_ev-1p5"]) *
         math.sqrt(float(row["energy_ev"])) * 0.05)
        for row in rows(reference / "eepf.csv")]
    candidate_mean = distribution_mean(candidate_eedf)
    reference_mean = distribution_mean(reference_eedf)
    return {
        "schema_version": 1, "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "descriptive_fresh_window_cross_code_pilot_comparison",
        "candidate_measurement_report_sha256": sha256(report_path),
        "reference_hashes": REFERENCE_HASHES,
        "density": density_metrics,
        "ion_impact_energy_distribution": ifed_metrics,
        "electron_energy_distribution": {
            "total_variation": total_variation(candidate_eedf, reference_eedf),
            "candidate_mean_energy_eV": candidate_mean,
            "reference_mean_energy_eV": reference_mean,
            "mean_energy_relative_difference":
                relative_difference(candidate_mean, reference_mean),
        },
        "acceptance": {"thresholds_declared": False, "passes": None},
        "reference_sampling_boundary": (
            "The 1024-cycle eduPIC reference passed density drift and shape "
            "gates but failed its predeclared minimum effective-block gate."),
        "claim_boundary": (
            "This four-cycle candidate pilot localizes discrepancies; it does "
            "not establish converged cross-code agreement or physical validation."),
        "physics_claim": "none_descriptive_cross_code_pilot_only",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_output", type=Path)
    parser.add_argument("reference_dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = analyze(args.candidate_output.resolve(), args.reference_dir.resolve())
    if args.output:
        atomic_json(args.output, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
