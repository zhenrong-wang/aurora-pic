#!/usr/bin/env python3
"""Regression tests for uncertainty-aware swarm reference comparison."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
COMPARATOR = ROOT / "scripts" / "compare_swarm.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_manifest(path: Path, data_file: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "[reference]",
                "swarm_reference_version = 2",
                f"data_file = {data_file.name}",
                "reference_id = aurorapic.synthetic.reference",
                "reference_version = 1",
                "gas = synthetic_swarm_gas",
                "population_model = fixed_population_no_avalanche",
                "coefficient_convention = flux_fixed_population",
                "provenance = AuroraPIC synthetic comparator fixture",
                "citation = AuroraPIC synthetic fixture",
                "retrieved = 2026-07-28",
                "license = Synthetic test data",
                "neutral_temperature_k = 300",
                "field_absolute_tolerance_td = 1e-12",
                "field_relative_tolerance = 1e-12",
                "",
                "[observable.drift]",
                "simulation_column = electron_drift_velocity_m_s",
                "reference_column = drift_velocity_m_s",
                "simulation_uncertainty_column = "
                "mean_velocity_x_standard_error_m_s",
                "reference_uncertainty_column = "
                "drift_velocity_standard_uncertainty_m_s",
                "relative_tolerance = 0.05",
                "absolute_tolerance = 0",
                "uncertainty_multiplier = 2",
                "",
                "[observable.mean_energy]",
                "simulation_column = mean_energy_ev",
                "reference_column = mean_energy_ev",
                "simulation_uncertainty_column = "
                "mean_energy_standard_error_ev",
                "reference_uncertainty_column = "
                "mean_energy_standard_uncertainty_ev",
                "relative_tolerance = 0.02",
                "absolute_tolerance = 0",
                "uncertainty_multiplier = 2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def run(
    simulation: Path,
    manifest: Path,
    output: Path,
    *,
    overwrite: bool = False,
) -> subprocess.CompletedProcess[str]:
    command = [
        sys.executable,
        str(COMPARATOR),
        str(simulation),
        str(manifest),
        "--output",
        str(output),
    ]
    if overwrite:
        command.append("--overwrite")
    return subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def main() -> int:
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_swarm_compare_"
    ) as temporary:
        work = Path(temporary)
        simulation = work / "simulation.csv"
        reference = work / "reference.csv"
        manifest = work / "reference.swarm-reference"
        report = work / "report.json"
        simulation.write_text(
            "dataset_id,dataset_version,gas,population_model,"
            "collision_model_signature,neutral_temperature_k,"
            "reduced_field_td,"
            "electron_drift_velocity_m_s,"
            "mean_velocity_x_standard_error_m_s,mean_energy_ev,"
            "mean_energy_standard_error_ev\n"
            "synthetic.dataset,1,synthetic_swarm_gas,"
            "fixed_population_no_avalanche,12345,300,1,100,2,1.0,0.02\n"
            "synthetic.dataset,1,synthetic_swarm_gas,"
            "fixed_population_no_avalanche,12345,300,5,200,3,2.0,0.03\n"
            "synthetic.dataset,1,synthetic_swarm_gas,"
            "fixed_population_no_avalanche,12345,300,10,300,4,3.0,0.04\n",
            encoding="utf-8",
        )
        reference.write_text(
            "reduced_field_td,drift_velocity_m_s,"
            "drift_velocity_standard_uncertainty_m_s,mean_energy_ev,"
            "mean_energy_standard_uncertainty_ev\n"
            "1,102,1,1.01,0.01\n"
            "5,195,2,1.98,0.02\n",
            encoding="utf-8",
        )
        write_manifest(manifest, reference)
        passed = run(simulation, manifest, report)
        require(
            passed.returncode == 0,
            f"valid swarm comparison failed: {passed.stderr}",
        )
        result = json.loads(report.read_text(encoding="utf-8"))
        require(
            result["passed"]
            and result["matched_reference_points"] == 2
            and result["extra_simulation_fields_td"] == [10.0]
            and len(result["simulation"]["sha256"]) == 64
            and len(result["reference"]["data_sha256"]) == 64,
            "passing swarm report is incomplete",
        )

        existing = run(simulation, manifest, report)
        require(
            existing.returncode == 2
            and "already exists" in existing.stderr,
            "swarm comparator overwrote a report without opt-in",
        )
        overwritten = run(
            simulation, manifest, report, overwrite=True
        )
        require(
            overwritten.returncode == 0,
            "swarm comparator rejected explicit report overwrite",
        )

        failed_reference = work / "failed_reference.csv"
        failed_manifest = work / "failed.swarm-reference"
        failed_report = work / "failed_report.json"
        failed_reference.write_text(
            "reduced_field_td,drift_velocity_m_s,"
            "drift_velocity_standard_uncertainty_m_s,mean_energy_ev,"
            "mean_energy_standard_uncertainty_ev\n"
            "1,102,1,1.01,0.01\n"
            "5,400,2,1.98,0.02\n",
            encoding="utf-8",
        )
        write_manifest(failed_manifest, failed_reference)
        failed = run(
            simulation, failed_manifest, failed_report
        )
        require(
            failed.returncode == 1
            and "did not meet acceptance criteria" in failed.stderr,
            "out-of-tolerance swarm comparison did not fail",
        )
        failed_result = json.loads(
            failed_report.read_text(encoding="utf-8")
        )
        require(
            not failed_result["passed"]
            and not failed_result["comparisons"][1]["observables"][0][
                "passed"
            ],
            "failed swarm report did not identify the residual",
        )

        mismatch_simulation = work / "mismatch.csv"
        mismatch_report = work / "mismatch_report.json"
        mismatch_simulation.write_text(
            simulation.read_text(encoding="utf-8").replace(
                "synthetic_swarm_gas", "different_gas"
            ),
            encoding="utf-8",
        )
        mismatch = run(
            mismatch_simulation, manifest, mismatch_report
        )
        require(
            mismatch.returncode == 2
            and "do not match reference gas" in mismatch.stderr
            and not mismatch_report.exists(),
            "swarm comparator accepted a gas identity mismatch",
        )

        temperature_simulation = work / "temperature.csv"
        temperature_report = work / "temperature_report.json"
        temperature_simulation.write_text(
            simulation.read_text(encoding="utf-8").replace(
                ",300,1,100", ",0,1,100"
            ),
            encoding="utf-8",
        )
        temperature_mismatch = run(
            temperature_simulation, manifest, temperature_report
        )
        require(
            temperature_mismatch.returncode == 2
            and "neutral_temperature_k values do not match"
            in temperature_mismatch.stderr
            and not temperature_report.exists(),
            "swarm comparator accepted a neutral-temperature mismatch",
        )

    print("swarm comparison validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
