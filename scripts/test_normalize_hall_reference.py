#!/usr/bin/env python3
"""Bounded end-to-end regression for Hall reference normalization."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "examples" / "hall_landmark_axial_azimuthal.case"
NORMALIZER = ROOT / "scripts" / "normalize_hall_reference.py"
COMPARATOR = ROOT / "scripts" / "compare_hall.py"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def write_simulation(output: Path) -> None:
    output.mkdir()
    (output / "resolved_field_time_average.csv").write_text(
        "start_time,end_time,duration,samples,profile_axis,"
        "coordinate,potential,electric_x,electric_y,charge_density\n"
        "1,2,1,3,x,0,200,100,0,0\n"
        "1,2,1,3,x,0.0125,100,200,0,0\n"
        "1,2,1,3,x,0.025,0,300,0,0\n",
        encoding="utf-8",
    )
    rows = [
        ("0", "electrons", "10", "10"),
        ("0.0125", "electrons", "20", "11"),
        ("0.025", "electrons", "30", "12"),
        ("0", "ions", "9", "1"),
        ("0.0125", "ions", "19", "1"),
        ("0.025", "ions", "29", "1"),
    ]
    with (output / "resolved_species_time_average.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow([
            "start_time", "end_time", "duration", "samples",
            "profile_axis", "coordinate", "species",
            "macro_particle_equivalent", "represented_number",
            "number_density", "mean_velocity_x", "mean_velocity_y",
            "mean_velocity_z", "thermal_speed_x", "thermal_speed_y",
            "thermal_speed_z", "temperature_ev", "current_density_x",
            "current_density_y", "current_density_z",
        ])
        for coordinate, species, density, temperature in rows:
            writer.writerow([
                1, 2, 1, 3, "x", coordinate, species, 1, 1, density,
                0, 0, 0, 1, 1, 1, temperature, 0, 0, 0,
            ])
    (output / "resolved_modes.csv").write_text(
        "step,time,mode_axis,mode,wavenumber,quantity,species,"
        "real,imaginary,amplitude\n"
        "0,0,y,1,490.8738521234052,electric_y,,1,0,1\n"
        "1,0.25,y,1,490.8738521234052,electric_y,,0,-1,1\n"
        "2,0.5,y,1,490.8738521234052,electric_y,,-1,0,1\n",
        encoding="utf-8",
    )


def main() -> int:
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_hall_normalize_"
    ) as temporary:
        work = Path(temporary)
        raw_profile = work / "published_profiles.csv"
        raw_profile.write_text(
            "x_cm,ex_a,ex_b,ni_a_cm3,ni_b_cm3,te_a,te_b\n"
            "0,99,101,8e-6,10e-6,9.5,10.5\n"
            "1.25,198,202,18e-6,20e-6,10.5,11.5\n"
            "2.5,297,303,28e-6,30e-6,11.5,12.5\n",
            encoding="utf-8",
        )
        raw_mode = work / "published_modes.csv"
        raw_mode.write_text(
            "azimuthal_mode,f_a_khz,f_b_khz\n"
            "1,0.00099,0.00101\n"
            "2,0.0019,0.0021\n",
            encoding="utf-8",
        )
        lock = work / "case2.hall-source"
        lock.write_text(
            "[source]\n"
            "hall_source_version = 1\n"
            "case_id = landmark-axial-azimuthal-2019\n"
            "case_variant = synthetic-multicode-envelope\n"
            f"case_manifest_sha256 = {digest(CASE)}\n"
            f"profile_file = {raw_profile.name}\n"
            f"profile_sha256 = {digest(raw_profile)}\n"
            f"mode_file = {raw_mode.name}\n"
            f"mode_sha256 = {digest(raw_mode)}\n"
            "source_url = https://example.invalid/hall-synthetic\n"
            "source_artifact_id = synthetic-fixture-v1\n"
            "provenance = AuroraPIC synthetic multi-code table\n"
            "citation = AuroraPIC normalization regression fixture\n"
            "retrieved = 2026-07-29\n"
            "license = synthetic test data\n"
            "\n[profile]\n"
            "coordinate_column = x_cm\n"
            "coordinate_scale_to_m = 0.01\n"
            "electric_field_columns = ex_a,ex_b\n"
            "electric_field_scale_to_v_m = 1\n"
            "ion_density_columns = ni_a_cm3,ni_b_cm3\n"
            "ion_density_scale_to_m3 = 1e6\n"
            "electron_temperature_columns = te_a,te_b\n"
            "electron_temperature_scale_to_ev = 1\n"
            "\n[mode]\n"
            "mode_column = azimuthal_mode\n"
            "frequency_columns = f_a_khz,f_b_khz\n"
            "frequency_scale_to_hz = 1000\n"
            "comparison_mode = 1\n"
            "\n[acceptance]\n"
            "coordinate_absolute_tolerance_m = 1e-12\n"
            "relative_tolerance = 0.01\n"
            "uncertainty_multiplier = 2\n",
            encoding="utf-8",
        )
        normalized = work / "normalized"
        result = run([
            sys.executable,
            str(NORMALIZER),
            str(lock),
            "--case-manifest",
            str(CASE),
            "--output-dir",
            str(normalized),
        ])
        require(
            result.returncode == 0,
            f"Hall reference normalization failed: {result.stderr}",
        )
        profile_rows = list(csv.DictReader(
            (normalized / "profiles.csv").open(encoding="utf-8")
        ))
        require(
            len(profile_rows) == 3
            and float(profile_rows[0]["coordinate_m"]) == 0.0
            and float(profile_rows[1]["coordinate_m"]) == 0.0125
            and float(profile_rows[0]["electric_x_v_m"]) == 100.0
            and float(profile_rows[0]["electric_x_uncertainty_v_m"]) == 1.0
            and float(profile_rows[0]["ion_density_m3"]) == 9.0
            and float(profile_rows[0]["ion_density_uncertainty_m3"]) == 1.0,
            "normalized profile units or envelope are incorrect",
        )
        audit = json.loads(
            (normalized / "normalization.json").read_text(encoding="utf-8")
        )
        require(
            audit["coordinate_policy"] == "native_no_interpolation"
            and audit["envelope_method"] == "midpoint_and_half_range"
            and audit["profile_rows"] == 3
            and audit["mode_rows"] == 2,
            "normalization audit is incomplete",
        )

        simulation = work / "simulation"
        write_simulation(simulation)
        comparison = work / "comparison.json"
        compared = run([
            sys.executable,
            str(COMPARATOR),
            str(simulation),
            str(normalized / "reference.hall-reference"),
            "--case-manifest",
            str(CASE),
            "--output",
            str(comparison),
        ])
        require(
            compared.returncode == 0
            and json.loads(comparison.read_text(encoding="utf-8"))["passed"],
            f"normalized reference was not comparator-ready: {compared.stderr}",
        )

        tampered_profile = raw_profile.read_text(encoding="utf-8").replace(
            "99,101", "90,110"
        )
        raw_profile.write_text(tampered_profile, encoding="utf-8")
        rejected_output = work / "rejected"
        rejected = run([
            sys.executable,
            str(NORMALIZER),
            str(lock),
            "--case-manifest",
            str(CASE),
            "--output-dir",
            str(rejected_output),
        ])
        require(
            rejected.returncode == 2
            and "profile SHA-256 mismatch" in rejected.stderr
            and not rejected_output.exists(),
            "normalizer accepted tampered raw reference data",
        )

    print("Hall reference normalization validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
