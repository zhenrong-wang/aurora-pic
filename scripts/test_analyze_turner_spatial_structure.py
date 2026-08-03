#!/usr/bin/env python3
"""Bounded regression for Turner spatial-structure analysis."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
ANALYZER = ROOT / "scripts" / "analyze_turner_spatial_structure.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_csv(path: Path, fields: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(rows)


def run(work: Path, metadata: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable, str(ANALYZER),
            "--density", str(work / "spatial_average.csv"),
            "--kinetic-energy", str(work / "spatial_kinetic_energy.csv"),
            "--field", str(work / "spatial_field_average.csv"),
            "--metadata", str(metadata),
            "--reported-midplane-ion-density", "10",
            "--reported-electron-temperature", str(4.0 / 3.0),
            "--output", str(work / "analysis.json"),
        ],
        cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_structure_", dir=ROOT / "tmp"
    ) as temporary:
        work = Path(temporary)
        x = [0.0, 0.25, 0.5, 0.75, 1.0]
        density_rows: list[list[object]] = []
        energy_rows: list[list[object]] = []
        for species_id, species, density, energy in (
            (0, "electrons", [1.0, 8.0, 10.0, 8.0, 1.0],
             [5.0, 3.0, 2.0, 3.0, 5.0]),
            (1, "ions", [10.0] * 5, [1.0] * 5),
        ):
            for node, coordinate in enumerate(x):
                density_rows.append(
                    [species_id, species, node, coordinate, density[node]])
                energy_rows.append(
                    [species_id, species, node, coordinate, energy[node],
                     2.0 * energy[node] / 3.0])
        write_csv(
            work / "spatial_average.csv",
            ["species_id", "species", "node", "x_m",
             "number_density_mean_m-3"], density_rows)
        write_csv(
            work / "spatial_kinetic_energy.csv",
            ["species_id", "species", "node", "x_m",
             "mean_kinetic_energy_eV", "effective_kinetic_temperature_eV"],
            energy_rows)
        write_csv(
            work / "spatial_field_average.csv",
            ["node", "x_m", "potential_mean_V", "electric_field_mean_V_m",
             "electric_field_rms_V_m"],
            [[node, coordinate, [0, 8, 10, 8, 0][node],
              [2, 1, 0, -1, -2][node],
              [10, 5, 1, 5, 10][node]]
             for node, coordinate in enumerate(x)])
        metadata = work / "metadata.json"
        metadata.write_text(json.dumps({
            "unit_system": "si", "complete": True,
            "moments_complete": True, "samples": 12800,
            "moment_samples": 12800, "start_step": 1, "end_step": 12800,
            "rf_cycles": 32,
        }), encoding="utf-8")

        completed = run(work, metadata)
        require(completed.returncode == 0,
                f"spatial-structure analysis failed: {completed.stderr}")
        result = json.loads((work / "analysis.json").read_text(encoding="utf-8"))
        proxy = result["density_ratio_sheath_proxies"]["ne_over_ni_0.90"]
        require(abs(proxy["left_width_m"] - 0.375) < 1e-14
                and abs(proxy["right_width_m"] - 0.375) < 1e-14,
                "density-ratio threshold interpolation is wrong")
        field = result["field_structure"]
        require(field["peak_rms_x_m"] == 0.0
                and field["peak_rms_V_m"] == 10.0
                and field["left_mean_potential_drop_from_center_V"] == 10.0
                and field["right_mean_potential_drop_from_center_V"] == 10.0
                and field["rms_field_symmetry_relative_l2"] == 0.0
                and field["mean_field_antisymmetry_relative_l2"] == 0.0,
                "symmetric field metrics are wrong")
        energy = result["electron_energy_structure"]
        require(energy["kinetic_energy_symmetry_relative_l2"] == 0.0
                and 0.0 < energy["outer_20_percent_kinetic_energy_fraction"] < 1.0,
                "electron-energy localization metrics are wrong")
        published = result["published_case1_scalar_context"]
        require(published["midplane_ion_density_relative_difference"] == 0.0
                and abs(published["electron_temperature_relative_difference"])
                < 1e-14,
                "published scalar context is wrong")

        incomplete = work / "incomplete.json"
        incomplete.write_text(metadata.read_text(encoding="utf-8").replace(
            '"moments_complete": true', '"moments_complete": false'),
            encoding="utf-8")
        rejected = run(work, incomplete)
        require(rejected.returncode != 0
                and "moment history is incomplete" in rejected.stderr,
                "analyzer accepted incomplete moment history")

    print("Turner spatial-structure analysis passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
