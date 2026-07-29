#!/usr/bin/env python3
"""Validate Hall convergence preparation and analysis without large runs."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "examples" / "hall_landmark_axial_azimuthal.case"
PREPARE = ROOT / "scripts" / "prepare_hall_convergence.py"
ANALYZE = ROOT / "scripts" / "analyze_hall_convergence.py"
ACK = "I_UNDERSTAND_THIS_IS_AN_OPT_IN_HALL_CONVERGENCE_PLAN"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    return subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
    )


def write_csv(
    path: Path, fields: list[str], rows: list[dict[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def synthetic_output(
    path: Path,
    multiplier: float,
    samples: int,
    nodes: int = 126,
    max_mode: int = 32,
) -> None:
    field_rows = []
    species_rows = []
    for index in range(nodes):
        coordinate = 0.025 * index / (nodes - 1)
        field_rows.append({
            "samples": samples,
            "profile_axis": "x",
            "coordinate": coordinate,
            "potential": multiplier * (10.0 + index),
            "electric_x": multiplier * (100.0 + index),
            "charge_density": multiplier * (0.1 + 0.01 * index),
        })
        for species_index, species in enumerate(("electrons", "ions"), 1):
            row: dict[str, object] = {
                "samples": samples,
                "profile_axis": "x",
                "coordinate": coordinate,
                "species": species,
            }
            for observable_index, observable in enumerate((
                "number_density", "temperature_ev",
                "mean_velocity_x", "mean_velocity_y",
                "current_density_x", "current_density_y",
            ), 1):
                row[observable] = (
                    multiplier
                    * (10 * species_index + observable_index + index)
                )
            species_rows.append(row)
    write_csv(
        path / "resolved_field_time_average.csv",
        [
            "samples", "profile_axis", "coordinate",
            "potential", "electric_x", "charge_density",
        ],
        field_rows,
    )
    write_csv(
        path / "resolved_species_time_average.csv",
        [
            "samples", "profile_axis", "coordinate", "species",
            "number_density", "temperature_ev",
            "mean_velocity_x", "mean_velocity_y",
            "current_density_x", "current_density_y",
        ],
        species_rows,
    )
    mode_rows = []
    for sample in range(samples):
        for mode in range(1, max_mode + 1):
            mode_rows.extend([
                {
                    "mode": mode,
                    "quantity": "electric_x",
                    "species": "",
                    "amplitude": multiplier * (mode + 1),
                },
                {
                    "mode": mode,
                    "quantity": "number_density",
                    "species": "electrons",
                    "amplitude": multiplier * (10 + mode),
                },
            ])
    write_csv(
        path / "resolved_modes.csv",
        ["mode", "quantity", "species", "amplitude"],
        mode_rows,
    )


def main() -> int:
    if len(sys.argv) != 2:
        raise RuntimeError("usage: test_hall_convergence.py <aurorapic_cli>")
    cli = Path(sys.argv[1]).resolve()
    require(cli.is_file(), f"missing AuroraPIC CLI: {cli}")
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_hall_convergence_"
    ) as temporary:
        work = Path(temporary)
        guarded = run([
            sys.executable, str(PREPARE), str(CASE),
            "--output-dir", str(work / "guarded"),
        ])
        require(
            guarded.returncode == 2 and ACK in guarded.stderr,
            "convergence planner bypassed its cost acknowledgement",
        )
        campaign_dir = work / "campaign"
        prepared = run([
            sys.executable, str(PREPARE), str(CASE),
            "--output-dir", str(campaign_dir),
            "--acknowledge-cost", ACK,
        ])
        require(
            prepared.returncode == 0,
            f"convergence preparation failed: {prepared.stderr}",
        )
        manifest_path = campaign_dir / "convergence.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        require(
            manifest["aggregate_initial_particle_updates"] == 7_680_000_000
            and len(manifest["runs"]) == 5
            and not manifest["launched"],
            "convergence manifest cost or launch contract is invalid",
        )
        multipliers = {
            "population_0p5": 1.2,
            "population_1": 1.0,
            "population_2": 1.05,
            "duration_0p5": 1.2,
            "duration_2": 1.05,
        }
        for item in manifest["runs"]:
            deck = campaign_dir / item["runtime_config"]
            validated = run([str(cli), "--validate-only", str(deck)])
            require(
                validated.returncode == 0,
                f"generated convergence deck is invalid: {validated.stderr}",
            )
            synthetic_output(
                Path(item["result_dir"]),
                multipliers[item["stage"]],
                item["diagnostic_samples"],
            )
        report_path = work / "report.json"
        analyzed = run([
            sys.executable, str(ANALYZE), str(manifest_path),
            "--report", str(report_path),
        ])
        report = json.loads(report_path.read_text(encoding="utf-8"))
        require(
            analyzed.returncode == 0
            and report["passed"]
            and report["comparisons"]["population"]["passed"]
            and report["comparisons"]["duration"]["passed"],
            f"convergent synthetic campaign failed: {analyzed.stderr}",
        )

        population_fine = next(
            item for item in manifest["runs"]
            if item["stage"] == "population_2"
        )
        synthetic_output(
            Path(population_fine["result_dir"]),
            1.8,
            population_fine["diagnostic_samples"],
        )
        failed_path = work / "failed.json"
        failed = run([
            sys.executable, str(ANALYZE), str(manifest_path),
            "--report", str(failed_path),
        ])
        failed_report = json.loads(
            failed_path.read_text(encoding="utf-8")
        )
        require(
            failed.returncode == 1
            and not failed_report["passed"]
            and not failed_report["comparisons"]["population"]["passed"],
            "convergence analyzer accepted a divergent population stage",
        )
    print("Hall convergence preparation and analysis passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
