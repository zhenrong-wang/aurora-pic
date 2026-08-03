#!/usr/bin/env python3
"""Analyze restart-safe Turner CCP energy and sheath observables."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
from typing import Iterable


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_rows(path: Path) -> list[dict[str, str]]:
    require(path.is_file(), f"missing input: {path}")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    require(bool(rows), f"input has no data rows: {path}")
    return rows


def linear_value(x: list[float], y: list[float], point: float) -> float:
    require(len(x) == len(y) and len(x) >= 2, "invalid profile")
    require(x[0] <= point <= x[-1], "interpolation point is outside profile")
    for index in range(1, len(x)):
        if point <= x[index]:
            fraction = (point - x[index - 1]) / (x[index] - x[index - 1])
            return y[index - 1] + fraction * (y[index] - y[index - 1])
    return y[-1]


def integrate(x: list[float], y: list[float], left: float | None = None,
              right: float | None = None) -> float:
    require(len(x) == len(y) and len(x) >= 2, "invalid integration profile")
    a = x[0] if left is None else left
    b = x[-1] if right is None else right
    require(x[0] <= a < b <= x[-1], "invalid integration interval")
    points = [a] + [value for value in x if a < value < b] + [b]
    values = [linear_value(x, y, point) for point in points]
    return sum(
        0.5 * (values[index - 1] + values[index]) *
        (points[index] - points[index - 1])
        for index in range(1, len(points))
    )


def weighted_mean(x: list[float], values: list[float], weights: list[float],
                  left: float, right: float) -> float:
    denominator = integrate(x, weights, left, right)
    require(denominator > 0.0, "profile has zero regional weight")
    numerator = integrate(
        x, [value * weight for value, weight in zip(values, weights)],
        left, right)
    return numerator / denominator


def correlation(first: Iterable[float], second: Iterable[float]) -> float:
    a = list(first)
    b = list(second)
    require(len(a) == len(b) and len(a) >= 2, "invalid correlation vectors")
    mean_a = sum(a) / len(a)
    mean_b = sum(b) / len(b)
    numerator = sum((x - mean_a) * (y - mean_b) for x, y in zip(a, b))
    norm_a = sum((x - mean_a) ** 2 for x in a)
    norm_b = sum((y - mean_b) ** 2 for y in b)
    return numerator / math.sqrt(norm_a * norm_b) if norm_a > 0 and norm_b > 0 else 0.0


def parity_error(x: list[float], values: list[float], parity: int = 1) -> float:
    require(parity in (-1, 1), "profile parity must be -1 or 1")
    length = x[-1] - x[0]
    tolerance = max(1e-14, length * 1e-12)
    require(
        all(abs((x[index] + x[-1 - index]) - (x[0] + x[-1])) <= tolerance
            for index in range(len(x))),
        "symmetry metric requires a mirror-symmetric grid",
    )
    difference_squared = [
        (value - parity * values[-1 - index]) ** 2
        for index, value in enumerate(values)
    ]
    reference_squared = [value * value for value in values]
    denominator = integrate(x, reference_squared)
    return math.sqrt(integrate(x, difference_squared) / denominator) \
        if denominator > 0.0 else 0.0


def left_crossing(x: list[float], ratio: list[float], threshold: float) -> float | None:
    if ratio[0] >= threshold:
        return x[0]
    for index in range(1, len(x)):
        if ratio[index] >= threshold:
            previous = ratio[index - 1]
            if ratio[index] == previous:
                return x[index]
            fraction = (threshold - previous) / (ratio[index] - previous)
            return x[index - 1] + fraction * (x[index] - x[index - 1])
    return None


def profile_by_species(rows: list[dict[str, str]], species: str,
                       value_column: str) -> tuple[list[int], list[float], list[float]]:
    selected = [row for row in rows if row["species"] == species]
    require(bool(selected), f"species '{species}' is absent")
    selected.sort(key=lambda row: int(row["node"]))
    return (
        [int(row["node"]) for row in selected],
        [float(row["x_m"]) for row in selected],
        [float(row[value_column]) for row in selected],
    )


def analyze(args: argparse.Namespace) -> dict[str, object]:
    require(args.reported_midplane_ion_density > 0.0,
            "reported midplane ion density must be positive")
    require(args.reported_electron_temperature > 0.0,
            "reported electron temperature must be positive")
    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
    require(metadata.get("unit_system") == "si", "analysis requires SI output")
    require(metadata.get("complete") is True, "spatial average is incomplete")
    require(metadata.get("moments_complete") is True,
            "spatial moment history is incomplete")

    density_rows = read_rows(args.density)
    energy_rows = read_rows(args.kinetic_energy)
    field_rows = read_rows(args.field)
    electron_nodes, x, electron_density = profile_by_species(
        density_rows, args.electron_species, "number_density_mean_m-3")
    ion_nodes, ion_x, ion_density = profile_by_species(
        density_rows, args.ion_species, "number_density_mean_m-3")
    energy_nodes, energy_x, electron_energy = profile_by_species(
        energy_rows, args.electron_species, "mean_kinetic_energy_eV")
    _, _, electron_temperature = profile_by_species(
        energy_rows, args.electron_species, "effective_kinetic_temperature_eV")
    require(electron_nodes == ion_nodes == energy_nodes and x == ion_x == energy_x,
            "species density and energy grids differ")

    field_rows.sort(key=lambda row: int(row["node"]))
    field_nodes = [int(row["node"]) for row in field_rows]
    field_x = [float(row["x_m"]) for row in field_rows]
    require(field_nodes == electron_nodes and field_x == x,
            "density and field grids differ")
    potential = [float(row["potential_mean_V"]) for row in field_rows]
    electric_mean = [float(row["electric_field_mean_V_m"]) for row in field_rows]
    electric_rms = [float(row["electric_field_rms_V_m"]) for row in field_rows]

    origin = x[0]
    length = x[-1] - origin
    require(length > 0.0, "profile length must be positive")
    center = origin + 0.5 * length
    ratio = [ne / ni if ni > 0.0 else 0.0
             for ne, ni in zip(electron_density, ion_density)]
    reverse_x = [origin + x[-1] - value for value in reversed(x)]
    reverse_ratio = list(reversed(ratio))
    crossings: dict[str, object] = {}
    for threshold in (0.8, 0.9, 0.95):
        left_edge = left_crossing(x, ratio, threshold)
        right_mirrored = left_crossing(reverse_x, reverse_ratio, threshold)
        right_edge = None if right_mirrored is None else origin + x[-1] - right_mirrored
        crossings[f"ne_over_ni_{threshold:.2f}"] = {
            "left_edge_x_m": left_edge,
            "right_edge_x_m": right_edge,
            "left_width_m": None if left_edge is None else left_edge - origin,
            "right_width_m": None if right_edge is None else x[-1] - right_edge,
        }

    electric_squared = [value * value for value in electric_rms]
    field_total = integrate(x, electric_squared)
    edge_right_start = origin + 0.8 * length
    edge_left_end = origin + 0.2 * length
    edge_field = integrate(x, electric_squared, origin, edge_left_end) + \
        integrate(x, electric_squared, edge_right_start, x[-1])
    kinetic_density = [density * energy for density, energy
                       in zip(electron_density, electron_energy)]
    kinetic_total = integrate(x, kinetic_density)
    edge_kinetic = integrate(x, kinetic_density, origin, edge_left_end) + \
        integrate(x, kinetic_density, edge_right_start, x[-1])
    bulk_left = origin + 0.25 * length
    bulk_right = origin + 0.75 * length

    peak_index = max(range(len(x)), key=lambda index: electric_rms[index])
    left_indices = [index for index, coordinate in enumerate(x)
                    if coordinate <= center]
    right_indices = [index for index, coordinate in enumerate(x)
                     if coordinate >= center]
    left_peak_index = max(left_indices, key=lambda index: electric_rms[index])
    right_peak_index = max(right_indices, key=lambda index: electric_rms[index])
    peak_average = 0.5 * (
        electric_rms[right_peak_index] + electric_rms[left_peak_index])
    result: dict[str, object] = {
        "analysis_version": 1,
        "inputs": {
            "density": str(args.density),
            "kinetic_energy": str(args.kinetic_energy),
            "field": str(args.field),
            "metadata": str(args.metadata),
        },
        "contract": {
            "samples": metadata["samples"],
            "moment_samples": metadata["moment_samples"],
            "start_step": metadata["start_step"],
            "end_step": metadata["end_step"],
            "rf_cycles": metadata["rf_cycles"],
            "electron_species": args.electron_species,
            "ion_species": args.ion_species,
        },
        "geometry": {"left_x_m": origin, "right_x_m": x[-1],
                     "length_m": length, "nodes": len(x)},
        "density_structure": {
            "center_ne_over_ni": linear_value(x, ratio, center),
            "maximum_ne_over_ni": max(ratio),
            "electron_density_symmetry_relative_l2":
                parity_error(x, electron_density),
            "ion_density_symmetry_relative_l2":
                parity_error(x, ion_density),
        },
        "published_case1_scalar_context": {
            "acceptance_gate": "none_diagnostic_context_only",
            "midplane_ion_density_m-3":
                linear_value(x, ion_density, center),
            "turner_table_iii_midplane_ion_density_m-3":
                args.reported_midplane_ion_density,
            "midplane_ion_density_relative_difference":
                linear_value(x, ion_density, center) /
                args.reported_midplane_ion_density - 1.0,
            "midplane_effective_electron_temperature_eV":
                linear_value(x, electron_temperature, center),
            "turner_table_iii_electron_temperature_eV":
                args.reported_electron_temperature,
            "electron_temperature_relative_difference":
                linear_value(x, electron_temperature, center) /
                args.reported_electron_temperature - 1.0,
        },
        "density_ratio_sheath_proxies": crossings,
        "field_structure": {
            "peak_rms_V_m": electric_rms[peak_index],
            "peak_rms_x_m": x[peak_index],
            "left_peak_rms_V_m": electric_rms[left_peak_index],
            "left_peak_rms_x_m": x[left_peak_index],
            "right_peak_rms_V_m": electric_rms[right_peak_index],
            "right_peak_rms_x_m": x[right_peak_index],
            "peak_rms_left_right_relative_difference":
                (electric_rms[right_peak_index] - electric_rms[left_peak_index]) /
                peak_average if peak_average > 0.0 else 0.0,
            "outer_20_percent_field_squared_fraction":
                edge_field / field_total if field_total > 0.0 else 0.0,
            "left_mean_potential_drop_from_center_V":
                linear_value(x, potential, center) - potential[0],
            "right_mean_potential_drop_from_center_V":
                linear_value(x, potential, center) - potential[-1],
            "rms_field_symmetry_relative_l2": parity_error(x, electric_rms),
            "mean_field_antisymmetry_relative_l2":
                parity_error(x, electric_mean, -1),
        },
        "electron_energy_structure": {
            "bulk_density_weighted_mean_kinetic_energy_eV": weighted_mean(
                x, electron_energy, electron_density, bulk_left, bulk_right),
            "bulk_density_weighted_effective_temperature_eV": weighted_mean(
                x, electron_temperature, electron_density,
                bulk_left, bulk_right),
            "left_outer_20_percent_density_weighted_mean_kinetic_energy_eV":
                weighted_mean(x, electron_energy, electron_density,
                              origin, edge_left_end),
            "left_outer_20_percent_density_weighted_effective_temperature_eV":
                weighted_mean(x, electron_temperature, electron_density,
                              origin, edge_left_end),
            "right_outer_20_percent_density_weighted_mean_kinetic_energy_eV":
                weighted_mean(x, electron_energy, electron_density,
                              edge_right_start, x[-1]),
            "right_outer_20_percent_density_weighted_effective_temperature_eV":
                weighted_mean(x, electron_temperature, electron_density,
                              edge_right_start, x[-1]),
            "outer_20_percent_kinetic_energy_fraction":
                edge_kinetic / kinetic_total if kinetic_total > 0.0 else 0.0,
            "kinetic_energy_symmetry_relative_l2":
                parity_error(x, electron_energy),
            "effective_temperature_symmetry_relative_l2":
                parity_error(x, electron_temperature),
        },
        "coupling_indicators": {
            "node_correlation_ne_vs_rms_field":
                correlation(electron_density, electric_rms),
            "node_correlation_electron_energy_vs_rms_field":
                correlation(electron_energy, electric_rms),
        },
        "interpretation_limits": [
            "Density-ratio crossings are threshold-sensitive sheath proxies, not a unique sheath-edge definition.",
            "Effective kinetic temperature includes directed kinetic energy.",
            "Node correlations are descriptive and do not establish causality.",
        ],
    }
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--density", type=Path, required=True)
    parser.add_argument("--kinetic-energy", type=Path, required=True)
    parser.add_argument("--field", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--electron-species", default="electrons")
    parser.add_argument("--ion-species", default="ions")
    parser.add_argument(
        "--reported-midplane-ion-density", type=float, default=1.40e14)
    parser.add_argument(
        "--reported-electron-temperature", type=float, default=9.36)
    return parser.parse_args()


def main() -> int:
    try:
        args = parse_args()
        result = analyze(args)
        encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
        if args.output is None:
            sys.stdout.write(encoded)
        else:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(encoded, encoding="utf-8")
            print(args.output)
        return 0
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
