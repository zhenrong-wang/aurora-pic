#!/usr/bin/env python3
"""Analyze SI 1D spatial and RF-phase collision-energy diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        result = list(csv.DictReader(stream))
    if not result:
        raise ValueError(f"{path} contains no diagnostic rows")
    return result


def node_volumes(coordinates: dict[int, float], boundary: str) -> dict[int, float]:
    ordered = sorted(coordinates)
    if ordered != list(range(len(ordered))) or len(ordered) < 2:
        raise ValueError("spatial collision nodes must be contiguous")
    spacings = [coordinates[i + 1] - coordinates[i] for i in ordered[:-1]]
    dx = sum(spacings) / len(spacings)
    if not math.isfinite(dx) or dx <= 0.0 or any(
        abs(value - dx) > 1.0e-10 * dx for value in spacings
    ):
        raise ValueError("spatial collision coordinates are not uniform")
    volumes = {node: dx for node in ordered}
    if boundary == "dirichlet":
        volumes[ordered[0]] *= 0.5
        volumes[ordered[-1]] *= 0.5
    return volumes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--boundary", choices=("dirichlet", "periodic"),
                        default="dirichlet")
    parser.add_argument("--edge-fraction", type=float, default=0.1)
    parser.add_argument("--absolute-tolerance-j-m2", type=float, default=1e-15)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if not 0.0 < args.edge_fraction < 0.5:
        raise ValueError("edge fraction must lie between zero and one half")

    collision_rows = rows(args.output_dir / "collisions.csv")
    spatial_rows = rows(args.output_dir / "spatial_collision_power.csv")
    phase_rows = rows(args.output_dir / "spatial_phase_collision_power.csv")
    coordinates = {
        int(row["node"]): float(row["x_m"]) for row in spatial_rows
    }
    volumes = node_volumes(coordinates, args.boundary)
    length = (max(coordinates.values()) +
              (volumes[0] if args.boundary == "periodic" else 0.0))

    cumulative_prefix = "cumulative_tracked_kinetic_energy_change_"
    cumulative_suffix = "_J_m-2"
    global_energy: dict[str, float] = {}
    for column in collision_rows[0]:
        if column.startswith(cumulative_prefix) and column.endswith(
            cumulative_suffix
        ):
            channel = column[len(cumulative_prefix):-len(cumulative_suffix)]
            global_energy[channel] = (
                float(collision_rows[-1][column]) -
                float(collision_rows[0][column])
            )
    if not global_energy:
        raise ValueError("global SI collision-energy columns are missing")

    spatial_energy: defaultdict[str, float] = defaultdict(float)
    region_energy: defaultdict[str, dict[str, float]] = defaultdict(
        lambda: {"left_edge": 0.0, "bulk": 0.0, "right_edge": 0.0}
    )
    peak_loss: dict[str, tuple[float, float]] = {}
    durations = set()
    timestep_counts = set()
    for row in spatial_rows:
        channel = row["channel"]
        node = int(row["node"])
        x = float(row["x_m"])
        energy = float(row["energy_density_sum_J_m-3"]) * volumes[node]
        spatial_energy[channel] += energy
        region = (
            "left_edge" if x <= args.edge_fraction * length else
            "right_edge" if x >= (1.0 - args.edge_fraction) * length else
            "bulk"
        )
        region_energy[channel][region] += energy
        power_density = float(row["mean_power_density_W_m-3"])
        if channel not in peak_loss or power_density < peak_loss[channel][1]:
            peak_loss[channel] = (x, power_density)
        durations.add(float(row["duration_s"]))
        timestep_counts.add(int(row["timesteps"]))
    if len(durations) != 1 or len(timestep_counts) != 1:
        raise ValueError("spatial rows do not share one duration")
    duration = durations.pop()
    timesteps = timestep_counts.pop()

    phase_energy: defaultdict[str, float] = defaultdict(float)
    phase_channel_energy: defaultdict[str, defaultdict[int, float]] = (
        defaultdict(lambda: defaultdict(float))
    )
    phase_durations: dict[int, float] = {}
    phase_timesteps: dict[int, int] = {}
    for row in phase_rows:
        phase = int(row["phase_bin"])
        channel = row["channel"]
        node = int(row["node"])
        energy = float(row["energy_density_sum_J_m-3"]) * volumes[node]
        phase_energy[channel] += energy
        phase_channel_energy[channel][phase] += energy
        phase_durations[phase] = float(row["duration_s"])
        phase_timesteps[phase] = int(row["timesteps"])

    channels = sorted(global_energy)
    if sorted(spatial_energy) != channels or sorted(phase_energy) != channels:
        raise ValueError("global, spatial, and phase channel sets differ")
    spatial_residual = {
        channel: spatial_energy[channel] - global_energy[channel]
        for channel in channels
    }
    phase_residual = {
        channel: phase_energy[channel] - spatial_energy[channel]
        for channel in channels
    }
    max_spatial_residual = max(abs(value) for value in spatial_residual.values())
    max_phase_residual = max(abs(value) for value in phase_residual.values())

    summaries = {}
    for channel in channels:
        phase_power = {
            str(phase): (
                phase_channel_energy[channel][phase] / phase_durations[phase]
                if phase_durations[phase] > 0.0 else 0.0
            )
            for phase in sorted(phase_durations)
        }
        maximum_loss_phase = min(phase_power, key=phase_power.get)
        ordered_phase_power = [
            phase_power[str(phase)] for phase in sorted(phase_durations)
        ]
        losses = [-value for value in ordered_phase_power]
        mean_loss = sum(losses) / len(losses)
        half_cycle_relative_l2 = None
        if len(losses) % 2 == 0:
            half = len(losses) // 2
            numerator = sum(
                (losses[index] - losses[index + half]) ** 2
                for index in range(half)
            )
            denominator = sum(
                (0.5 * (losses[index] + losses[index + half])) ** 2
                for index in range(half)
            )
            half_cycle_relative_l2 = math.sqrt(numerator / denominator)
        total = spatial_energy[channel]
        summaries[channel] = {
            "mean_power_W_m-2": total / duration,
            "regional_energy_fraction": {
                region: (value / total if total != 0.0 else 0.0)
                for region, value in region_energy[channel].items()
            },
            "outer_edge_energy_fraction": (
                region_energy[channel]["left_edge"] +
                region_energy[channel]["right_edge"]
            ) / total if total != 0.0 else 0.0,
            "peak_loss_power_density_x_m": peak_loss[channel][0],
            "peak_loss_power_density_W_m-3": peak_loss[channel][1],
            "maximum_loss_phase_bin": int(maximum_loss_phase),
            "maximum_loss_phase_fraction": (
                int(maximum_loss_phase) + 0.5
            ) / len(phase_durations),
            "phase_conditional_loss_minimum_W_m-2": min(losses),
            "phase_conditional_loss_maximum_W_m-2": max(losses),
            "phase_conditional_loss_peak_to_trough_over_mean": (
                (max(losses) - min(losses)) / mean_loss
                if mean_loss != 0.0 else 0.0
            ),
            "half_cycle_loss_relative_l2": half_cycle_relative_l2,
            "phase_conditional_power_W_m-2": phase_power,
        }

    phase_timestep_uniform = len(set(phase_timesteps.values())) == 1
    passes = (
        max_spatial_residual <= args.absolute_tolerance_j_m2 and
        max_phase_residual <= args.absolute_tolerance_j_m2 and
        sum(phase_timesteps.values()) == timesteps and
        phase_timestep_uniform
    )
    report = {
        "schema_version": 1,
        "boundary": args.boundary,
        "edge_fraction": args.edge_fraction,
        "timesteps": timesteps,
        "duration_s": duration,
        "phase_timesteps": {
            str(key): phase_timesteps[key] for key in sorted(phase_timesteps)
        },
        "phase_timestep_uniform": phase_timestep_uniform,
        "channels": summaries,
        "closure": {
            "spatial_minus_global_J_m-2": spatial_residual,
            "phase_minus_spatial_J_m-2": phase_residual,
            "maximum_spatial_global_residual_J_m-2": max_spatial_residual,
            "maximum_phase_spatial_residual_J_m-2": max_phase_residual,
            "absolute_tolerance_J_m-2": args.absolute_tolerance_j_m2,
            "passes": passes,
        },
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if passes else 1


if __name__ == "__main__":
    raise SystemExit(main())
