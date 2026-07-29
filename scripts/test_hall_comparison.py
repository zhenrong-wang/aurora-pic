#!/usr/bin/env python3
"""Bounded regression for Hall comparison and production preflight tooling."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
COMPARATOR = ROOT / "scripts" / "compare_hall.py"
PREFLIGHT = ROOT / "scripts" / "preflight_hall.py"
PRODUCTION_CASE = (
    ROOT / "examples" / "hall_landmark_axial_azimuthal.case"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def write_fixture(work: Path) -> tuple[Path, Path, Path]:
    output = work / "simulation"
    output.mkdir()
    field = output / "resolved_field_time_average.csv"
    species = output / "resolved_species_time_average.csv"
    modes = output / "resolved_modes.csv"
    field.write_text(
        "start_time,end_time,duration,samples,profile_axis,"
        "coordinate,potential,electric_x,electric_y,charge_density\n"
        "0,0.5,0.5,3,x,0,200,100,0,0\n"
        "0,0.5,0.5,3,x,0.5,100,200,0,0\n"
        "0,0.5,0.5,3,x,1,0,300,0,0\n",
        encoding="utf-8",
    )
    species.write_text(
        "start_time,end_time,duration,samples,profile_axis,"
        "coordinate,species,macro_particle_equivalent,"
        "represented_number,number_density,mean_velocity_x,"
        "mean_velocity_y,mean_velocity_z,thermal_speed_x,"
        "thermal_speed_y,thermal_speed_z,temperature_ev,"
        "current_density_x,current_density_y,current_density_z\n"
        "0,0.5,0.5,3,x,0,electrons,1,1,10,0,0,0,1,1,1,10,0,0,0\n"
        "0,0.5,0.5,3,x,0.5,electrons,1,1,20,0,0,0,1,1,1,11,0,0,0\n"
        "0,0.5,0.5,3,x,1,electrons,1,1,30,0,0,0,1,1,1,12,0,0,0\n"
        "0,0.5,0.5,3,x,0,ions,1,1,9,0,0,0,1,1,1,1,0,0,0\n"
        "0,0.5,0.5,3,x,0.5,ions,1,1,19,0,0,0,1,1,1,1,0,0,0\n"
        "0,0.5,0.5,3,x,1,ions,1,1,29,0,0,0,1,1,1,1,0,0,0\n",
        encoding="utf-8",
    )
    modes.write_text(
        "step,time,mode_axis,mode,wavenumber,quantity,species,"
        "real,imaginary,amplitude\n"
        "0,0,y,1,3.141592653589793,electric_y,,1,0,2\n"
        "1,0.25,y,1,3.141592653589793,electric_y,,0,-1,2\n"
        "2,0.5,y,1,3.141592653589793,electric_y,,-1,0,2\n",
        encoding="utf-8",
    )
    profile_reference = work / "profile.csv"
    profile_reference.write_text(
        "coordinate_m,electric_x_v_m,electric_x_uncertainty_v_m,"
        "ion_density_m3,electron_temperature_ev\n"
        "0,101,1,9,10\n"
        "0.5,198,1,19,11\n"
        "1,302,1,29,12\n",
        encoding="utf-8",
    )
    mode_reference = work / "modes.csv"
    mode_reference.write_text(
        "mode,frequency_hz,frequency_uncertainty_hz\n"
        "1,1.01,0.01\n",
        encoding="utf-8",
    )
    case = work / "case.case"
    case.write_text(
        "case_manifest_version = 1\n"
        "case_id = landmark-synthetic\n"
        "status = reduced_integration_only\n"
        "\n[reference]\n"
        "production_nx = 500\n"
        "production_ny = 256\n"
        "production_steps = 4000000\n"
        "production_dt_s = 5e-12\n"
        "\n[reduced_contract]\n"
        "missing_physics = published_resolution\n",
        encoding="utf-8",
    )
    manifest = work / "reference.hall-reference"
    manifest.write_text(
        "[reference]\n"
        "hall_reference_version = 1\n"
        "case_id = landmark-synthetic\n"
        "case_variant = synthetic-case\n"
        f"case_manifest_sha256 = {digest(case)}\n"
        f"profile_data_file = {profile_reference.name}\n"
        f"profile_data_sha256 = {digest(profile_reference)}\n"
        f"mode_data_file = {mode_reference.name}\n"
        f"mode_data_sha256 = {digest(mode_reference)}\n"
        "profile_axis = x\n"
        "mode_axis = y\n"
        "coordinate_column = coordinate_m\n"
        "coordinate_absolute_tolerance = 1e-12\n"
        "provenance = AuroraPIC synthetic comparator fixture\n"
        "citation = AuroraPIC synthetic test\n"
        "retrieved = 2026-07-29\n"
        "license = synthetic test data\n"
        "\n[profile.axial_field]\n"
        "simulation_source = field\n"
        "simulation_column = electric_x\n"
        "reference_column = electric_x_v_m\n"
        "reference_uncertainty_column = electric_x_uncertainty_v_m\n"
        "relative_tolerance = 0.01\n"
        "uncertainty_multiplier = 2\n"
        "\n[profile.ion_density]\n"
        "simulation_source = species\n"
        "simulation_species = ions\n"
        "simulation_column = number_density\n"
        "reference_column = ion_density_m3\n"
        "absolute_tolerance = 1e-12\n"
        "\n[profile.electron_temperature]\n"
        "simulation_source = species\n"
        "simulation_species = electrons\n"
        "simulation_column = temperature_ev\n"
        "reference_column = electron_temperature_ev\n"
        "absolute_tolerance = 1e-12\n"
        "\n[mode.dominant_frequency]\n"
        "simulation_quantity = electric_y\n"
        "mode = 1\n"
        "metric = frequency_hz\n"
        "reference_column = frequency_hz\n"
        "reference_uncertainty_column = frequency_uncertainty_hz\n"
        "relative_tolerance = 0.01\n"
        "uncertainty_multiplier = 2\n",
        encoding="utf-8",
    )
    return output, manifest, case


def main() -> int:
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_hall_compare_"
    ) as temporary:
        work = Path(temporary)
        output, manifest, case = write_fixture(work)
        report = work / "comparison.json"
        command = [
            sys.executable,
            str(COMPARATOR),
            str(output),
            str(manifest),
            "--case-manifest",
            str(case),
            "--output",
            str(report),
        ]
        passed = run(command)
        require(
            passed.returncode == 0,
            f"valid Hall comparison failed: {passed.stderr}",
        )
        result = json.loads(report.read_text(encoding="utf-8"))
        mode_result = result["mode_comparisons"][0]
        require(
            result["passed"]
            and len(result["profile_comparisons"]) == 3
            and abs(mode_result["simulation"] - 1.0) < 1e-14
            and len(result["reference"]["profile_sha256"]) == 64
            and len(result["simulation"]["case_manifest_sha256"]) == 64,
            "passing Hall comparison report is incomplete",
        )
        existing = run(command)
        require(
            existing.returncode == 2
            and "already exists" in existing.stderr,
            "Hall comparator overwrote a report without opt-in",
        )

        failed_reference = work / "failed_profile.csv"
        failed_reference.write_text(
            (work / "profile.csv").read_text(encoding="utf-8").replace(
                "0.5,198", "0.5,500"
            ),
            encoding="utf-8",
        )
        failed_manifest = work / "failed.hall-reference"
        failed_manifest.write_text(
            manifest.read_text(encoding="utf-8")
            .replace(
                "profile_data_file = profile.csv",
                "profile_data_file = failed_profile.csv",
            )
            .replace(
                f"profile_data_sha256 = {digest(work / 'profile.csv')}",
                f"profile_data_sha256 = {digest(failed_reference)}",
            ),
            encoding="utf-8",
        )
        failed_report = work / "failed.json"
        failed = run(
            [
                sys.executable,
                str(COMPARATOR),
                str(output),
                str(failed_manifest),
                "--case-manifest",
                str(case),
                "--output",
                str(failed_report),
            ]
        )
        require(
            failed.returncode == 1
            and "did not meet acceptance criteria" in failed.stderr
            and not json.loads(
                failed_report.read_text(encoding="utf-8")
            )["passed"],
            "out-of-tolerance Hall profile comparison did not fail",
        )

        tampered = work / "tampered.hall-reference"
        tampered.write_text(
            manifest.read_text(encoding="utf-8").replace(
                "profile_data_sha256 = ",
                "profile_data_sha256 = " + "0" * 64 + "\n# original ",
                1,
            ),
            encoding="utf-8",
        )
        tampered_result = run(
            [
                sys.executable,
                str(COMPARATOR),
                str(output),
                str(tampered),
                "--case-manifest",
                str(case),
                "--output",
                str(work / "tampered.json"),
            ]
        )
        require(
            tampered_result.returncode == 2
            and "SHA-256 mismatch" in tampered_result.stderr,
            "Hall comparator accepted a bad reference hash",
        )

        preflight_report = work / "preflight.json"
        preflight = run(
            [
                sys.executable,
                str(PREFLIGHT),
                str(PRODUCTION_CASE),
                "--report",
                str(preflight_report),
                "--particle-push-rate",
                "1e8",
            ]
        )
        require(
            preflight.returncode == 0,
            f"valid Hall preflight failed: {preflight.stderr}",
        )
        estimate = json.loads(
            preflight_report.read_text(encoding="utf-8")
        )
        require(
            estimate["production_scale"]
            and estimate["within_declared_budgets"]
            and not estimate["launch_authorized"]
            and estimate["estimates"]["initial_macroparticles"]
                == 500 * 256 * 75 * 2
            and estimate["estimates"]["particle_updates_lower_bound"]
                == 500 * 256 * 75 * 2 * 4_000_000
            and estimate["estimates"]["estimated_wall_seconds"] > 0,
            "Hall preflight report is incomplete",
        )
        budget_failure = run(
            [
                sys.executable,
                str(PREFLIGHT),
                str(PRODUCTION_CASE),
                "--report",
                str(work / "budget_failure.json"),
                "--memory-budget-gib",
                "0.1",
            ]
        )
        require(
            budget_failure.returncode == 1
            and "exceeded declared resource budgets"
                in budget_failure.stderr,
            "Hall preflight ignored a memory budget failure",
        )

    print("Hall comparison and preflight validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
