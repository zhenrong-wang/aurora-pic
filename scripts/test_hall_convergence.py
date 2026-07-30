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
    particles_per_cell: int,
    steps: int,
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
    macro_charge = 1.6e-10 / particles_per_cell
    created = 10
    represented = macro_charge * created / 1.602176634e-19
    reverse_steps = int(round(0.04 * steps))
    cumulative_reverse = 0.02 * steps
    negative_charge = -1e-11 * steps
    positive_charge = 2e-12 * steps
    write_csv(
        path / "current_source.csv",
        [
            "macro_particles_created",
            "represented_particles_created",
            "control_updates",
            "reverse_diagnostics_start_step",
            "reverse_demand_steps",
            "reverse_demand_step_fraction",
            "cumulative_reverse_demand_macroparticles",
            "maximum_reverse_demand_macroparticles",
            "cumulative_monitored_negative_charge",
            "cumulative_monitored_positive_charge",
            "cumulative_processed_monitored_charge",
        ],
        [{
            "macro_particles_created": created,
            "represented_particles_created": represented,
            "control_updates": steps,
            "reverse_diagnostics_start_step": 0,
            "reverse_demand_steps": reverse_steps,
            "reverse_demand_step_fraction": reverse_steps / steps,
            "cumulative_reverse_demand_macroparticles":
                cumulative_reverse,
            "maximum_reverse_demand_macroparticles": 1,
            "cumulative_monitored_negative_charge": negative_charge,
            "cumulative_monitored_positive_charge": positive_charge,
            "cumulative_processed_monitored_charge":
                negative_charge + positive_charge,
        }],
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
            manifest["hall_convergence_version"] == 2
            and manifest["aggregate_initial_particle_updates"]
                == 7_680_000_000
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
                item["particles_per_cell_per_species"],
                item["steps"],
            )
        report_path = work / "report.json"
        analyzed = run([
            sys.executable, str(ANALYZE), str(manifest_path),
            "--report", str(report_path),
        ])
        require(
            report_path.is_file(),
            f"convergence analyzer did not write a report: "
            f"{analyzed.stderr}",
        )
        report = json.loads(report_path.read_text(encoding="utf-8"))
        require(
            analyzed.returncode == 0
            and report["schema_version"] == 2
            and report["passed"]
            and report["comparisons"]["population"]["passed"]
            and report["comparisons"]["duration"]["passed"]
            and report["comparisons"]["controller_population"]["passed"]
            and report["comparisons"]["controller_population"][
                "fine_to_baseline_reverse_charge_per_update_ratio"
            ] == 0.5
            and report["comparisons"]["controller_population"][
                "fine_to_baseline_reverse_impulse_ratio"
            ] == 0.5,
            f"convergent synthetic campaign failed: {analyzed.stderr}",
        )
        population_report_path = work / "population-report.json"
        population_analyzed = run([
            sys.executable, str(ANALYZE), str(manifest_path),
            "--axis", "population",
            "--report", str(population_report_path),
        ])
        population_report = json.loads(
            population_report_path.read_text(encoding="utf-8")
        )
        require(
            population_analyzed.returncode == 0
            and population_report["passed"]
            and population_report["analyzed_axes"] == ["population"]
            and len(population_report["runs"]) == 3
            and "duration" not in population_report["comparisons"],
            "population-only convergence analysis failed",
        )

        population_fine = next(
            item for item in manifest["runs"]
            if item["stage"] == "population_2"
        )
        synthetic_output(
            Path(population_fine["result_dir"]),
            1.8,
            population_fine["diagnostic_samples"],
            population_fine["particles_per_cell_per_species"],
            population_fine["steps"],
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

        synthetic_output(
            Path(population_fine["result_dir"]),
            1.05,
            population_fine["diagnostic_samples"],
            population_fine["particles_per_cell_per_species"],
            population_fine["steps"],
        )
        controller_path = (
            Path(population_fine["result_dir"]) / "current_source.csv"
        )
        with controller_path.open(
            newline="", encoding="utf-8"
        ) as stream:
            controller_rows = list(csv.DictReader(stream))
            controller_fields = list(
                controller_rows[0].keys()
            )
        controller_rows[-1][
            "cumulative_reverse_demand_macroparticles"
        ] = "400"
        controller_rows[-1][
            "maximum_reverse_demand_macroparticles"
        ] = "4"
        write_csv(
            controller_path, controller_fields, controller_rows
        )
        controller_failed_path = work / "controller-failed.json"
        controller_failed = run([
            sys.executable, str(ANALYZE), str(manifest_path),
            "--report", str(controller_failed_path),
        ])
        controller_failed_report = json.loads(
            controller_failed_path.read_text(encoding="utf-8")
        )
        require(
            controller_failed.returncode == 1
            and not controller_failed_report["passed"]
            and not controller_failed_report["comparisons"][
                "controller_population"
            ]["passed"]
            and controller_failed_report["comparisons"]["population"][
                "passed"
            ],
            "convergence analyzer accepted non-convergent controller demand",
        )
    print("Hall convergence preparation and analysis passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
