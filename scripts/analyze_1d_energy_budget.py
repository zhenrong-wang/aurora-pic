#!/usr/bin/env python3
"""Close a 1D PIC tracked-particle kinetic-energy budget over a CSV window."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) < 2:
        raise ValueError(f"{path} must contain at least two samples")
    return rows


def delta(rows: list[dict[str, str]], column: str) -> float:
    return float(rows[-1][column]) - float(rows[0][column])


def selected_deltas(
    rows: list[dict[str, str]], predicate
) -> dict[str, float]:
    return {
        column: delta(rows, column)
        for column in rows[0]
        if predicate(column)
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument(
        "--absolute-tolerance-w-m2", type=float, default=1.0e-8
    )
    args = parser.parse_args()

    scalar_rows = read_rows(args.output_dir / "scalars.csv")
    collision_rows = read_rows(args.output_dir / "collisions.csv")
    power_rows = read_rows(args.output_dir / "power_transfer.csv")
    boundary_rows = read_rows(args.output_dir / "boundary_losses.csv")
    tables = (scalar_rows, collision_rows, power_rows, boundary_rows)
    endpoints = [
        (int(float(rows[0]["step"])), int(float(rows[-1]["step"])))
        for rows in tables
    ]
    if len(set(endpoints)) != 1:
        raise ValueError(f"diagnostic step windows differ: {endpoints}")
    duration = delta(scalar_rows, "time")
    if not math.isfinite(duration) or duration <= 0.0:
        raise ValueError("diagnostic duration must be positive and finite")

    collision_energy = selected_deltas(
        collision_rows,
        lambda name: name.startswith(
            "cumulative_tracked_kinetic_energy_change_"
        ),
    )
    electric_energy = selected_deltas(
        power_rows, lambda name: name.startswith("electric_work_")
    )
    wall_energy = selected_deltas(
        boundary_rows,
        lambda name: "_kinetic_energy_" in name and name.endswith("_J_m-2"),
    )
    if not collision_energy or not electric_energy or not wall_energy:
        raise ValueError("required SI energy-ledger columns are missing")

    collision_power = {
        name: value / duration for name, value in collision_energy.items()
    }
    electric_power = {
        name: value / duration for name, value in electric_energy.items()
    }
    wall_power = {name: value / duration for name, value in wall_energy.items()}
    kinetic_rate = delta(scalar_rows, "kinetic_energy") / duration
    predicted_rate = (
        sum(electric_power.values())
        + sum(collision_power.values())
        - sum(wall_power.values())
    )
    residual = kinetic_rate - predicted_rate
    scale = max(
        abs(kinetic_rate), abs(predicted_rate),
        abs(sum(electric_power.values())),
        abs(sum(collision_power.values())),
        abs(sum(wall_power.values())), 1.0,
    )
    report = {
        "schema_version": 1,
        "start_step": endpoints[0][0],
        "end_step": endpoints[0][1],
        "duration_s": duration,
        "electric_power_W_m-2": electric_power,
        "electric_power_total_W_m-2": sum(electric_power.values()),
        "collision_tracked_kinetic_power_W_m-2": collision_power,
        "collision_tracked_kinetic_power_total_W_m-2": sum(
            collision_power.values()
        ),
        "wall_kinetic_power_W_m-2": wall_power,
        "wall_kinetic_power_total_W_m-2": sum(wall_power.values()),
        "tracked_kinetic_energy_rate_W_m-2": kinetic_rate,
        "predicted_tracked_kinetic_energy_rate_W_m-2": predicted_rate,
        "closure_residual_W_m-2": residual,
        "relative_closure_residual": abs(residual) / scale,
        "absolute_tolerance_W_m-2": args.absolute_tolerance_w_m2,
        "passes": abs(residual) <= args.absolute_tolerance_w_m2,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if report["passes"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
