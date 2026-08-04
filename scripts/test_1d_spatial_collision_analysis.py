#!/usr/bin/env python3
"""Bounded regression for analyze_1d_spatial_collision.py."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def write_csv(path: Path, fields: list[str], data: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(data)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="aurorapic-spatial-collision-") as raw:
        output = Path(raw)
        energy_column = (
            "cumulative_tracked_kinetic_energy_change_test.loss_J_m-2"
        )
        write_csv(
            output / "collisions.csv",
            ["step", "time", energy_column],
            [[0, 0.0, 0.0], [2, 2.0, -3.0]],
        )
        spatial_fields = [
            "channel_id", "channel", "node", "x_m", "timesteps",
            "duration_s", "energy_density_sum_J_m-3",
            "mean_power_density_W_m-3",
        ]
        write_csv(
            output / "spatial_collision_power.csv",
            spatial_fields,
            [
                [0, "test.loss", 0, 0.0, 2, 2.0, -4.0, -2.0],
                [0, "test.loss", 1, 0.5, 2, 2.0, -2.0, -1.0],
                [0, "test.loss", 2, 1.0, 2, 2.0, -4.0, -2.0],
            ],
        )
        phase_fields = [
            "phase_bin", "phase_fraction", "timesteps", "duration_s",
            "channel_id", "channel", "node", "x_m",
            "energy_density_sum_J_m-3", "mean_power_density_W_m-3",
        ]
        phase_data = []
        for phase, fraction in ((0, 0.25), (1, 0.75)):
            for node, x, density in (
                (0, 0.0, -2.0), (1, 0.5, -1.0), (2, 1.0, -2.0)
            ):
                phase_data.append(
                    [phase, fraction, 1, 1.0, 0, "test.loss", node, x,
                     density, density]
                )
        write_csv(
            output / "spatial_phase_collision_power.csv",
            phase_fields, phase_data,
        )
        report_path = output / "report.json"
        analyzer = Path(__file__).with_name(
            "analyze_1d_spatial_collision.py"
        )
        completed = subprocess.run(
            [sys.executable, str(analyzer), str(output),
             "--json", str(report_path)],
            check=False, capture_output=True, text=True,
        )
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr or completed.stdout)
        report = json.loads(report_path.read_text(encoding="utf-8"))
        channel = report["channels"]["test.loss"]
        if not report["closure"]["passes"]:
            raise RuntimeError("synthetic spatial collision closure failed")
        if abs(channel["mean_power_W_m-2"] + 1.5) > 1.0e-14:
            raise RuntimeError("synthetic channel mean power is wrong")
        for fraction in channel["regional_energy_fraction"].values():
            if abs(fraction - 1.0 / 3.0) > 1.0e-14:
                raise RuntimeError("synthetic regional fraction is wrong")
        if report["phase_timesteps"] != {"0": 1, "1": 1}:
            raise RuntimeError("synthetic phase timestep counts are wrong")
    print("1D spatial collision analysis regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
