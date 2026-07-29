#!/usr/bin/env python3
"""Run AuroraPIC example configs as isolated CLI smoke tests.

The examples intentionally write diagnostics, fields, and optional particle samples.
This script copies each config to a temporary directory, rewrites only output_dir,
runs the CLI, and checks the generated files structurally so the smoke suite proves
more than process exit status while keeping repository output directories clean.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


class SmokeFailure(RuntimeError):
    """Raised when a smoke example does not produce expected diagnostics."""


def rewrite_output_dir(config_path: Path, temp_root: Path, output_name: str) -> tuple[Path, Path]:
    """Copy a config into temp_root with output_dir redirected to a temp child."""
    output_dir = temp_root / output_name
    copied_config = temp_root / config_path.name
    lines = config_path.read_text(encoding="utf-8").splitlines()
    rewritten: list[str] = []
    replaced = False
    for line in lines:
        if re.match(r"^\s*output_dir\s*=", line):
            rewritten.append(f"output_dir = {output_dir.as_posix()}")
            replaced = True
        else:
            rewritten.append(line)
    if not replaced:
        rewritten.append(f"output_dir = {output_dir.as_posix()}")
    copied_config.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
    for pattern in ("*.msh", "*.dat", "*.gas", "*.aps"):
        for support_path in EXAMPLES.glob(pattern):
            shutil.copy2(support_path, temp_root / support_path.name)
    return copied_config, output_dir


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SmokeFailure(message)


def require_file(path: Path) -> Path:
    require(path.is_file(), f"missing expected file: {path}")
    require(path.stat().st_size > 0, f"expected non-empty file: {path}")
    return path


def read_csv(path: Path) -> tuple[list[str], list[list[str]]]:
    require_file(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    require(len(rows) >= 2, f"expected header and at least one data row in {path}")
    require(rows[0], f"empty CSV header in {path}")
    for row_index, row in enumerate(rows[1:], start=2):
        require(len(row) == len(rows[0]), f"row {row_index} in {path} has {len(row)} columns; expected {len(rows[0])}")
    return rows[0], rows[1:]


def require_numeric_rows(header: Sequence[str], rows: Sequence[Sequence[str]], path: Path) -> None:
    for row_index, row in enumerate(rows, start=2):
        for column, value in zip(header, row):
            try:
                number = float(value)
            except ValueError as exc:
                raise SmokeFailure(f"non-numeric value in {path}:{row_index} column {column!r}: {value!r}") from exc
            require(math.isfinite(number), f"non-finite value in {path}:{row_index} column {column!r}: {value!r}")


def require_csv(path: Path, expected_header: Sequence[str] | None = None, min_rows: int = 1) -> tuple[list[str], list[list[str]]]:
    header, rows = read_csv(path)
    if expected_header is not None:
        require(header == list(expected_header), f"unexpected header in {path}: {header!r}")
    require(len(rows) >= min_rows, f"expected at least {min_rows} data rows in {path}; found {len(rows)}")
    require_numeric_rows(header, rows, path)
    return header, rows


def require_step(rows: Sequence[Sequence[str]], expected_step: int, path: Path) -> None:
    require(any(int(float(row[0])) == expected_step for row in rows), f"expected step {expected_step} in {path}")


def require_initialization_audits(output_dir: Path) -> None:
    require_file(output_dir / "initialization.csv")
    path = output_dir / "initialization_acceptance.csv"
    header, rows = read_csv(path)
    require(
        header == [
            "enabled", "overall_passed", "metric", "value", "scale",
            "relative_residual", "tolerance", "passed", "details",
        ],
        f"unexpected header in {path}: {header!r}",
    )
    require(len(rows) >= 1, f"expected an acceptance row in {path}")
    for row in rows:
        values = dict(zip(header, row))
        require_numeric_rows(
            ["enabled", "overall_passed", "value", "scale",
             "relative_residual", "tolerance", "passed"],
            [[
                values["enabled"], values["overall_passed"],
                values["value"], values["scale"],
                values["relative_residual"], values["tolerance"],
                values["passed"],
            ]],
            path,
        )
        require(
            values["overall_passed"] == "1" and
            values["passed"] == "1",
            "example initialization acceptance audit did not pass",
        )


def require_vtk(path: Path, dimensions: tuple[int, int, int]) -> None:
    require_file(path)
    text = path.read_text(encoding="utf-8")
    nx, ny, nz = dimensions
    points = nx * ny * nz
    required_fragments = [
        "# vtk DataFile Version 3.0",
        "ASCII",
        "DATASET STRUCTURED_GRID",
        f"DIMENSIONS {nx} {ny} {nz}",
        f"POINTS {points} double",
        f"POINT_DATA {points}",
        "SCALARS rho double 1",
        "SCALARS phi double 1",
        "VECTORS electric double",
    ]
    for fragment in required_fragments:
        require(fragment in text, f"missing VTK fragment {fragment!r} in {path}")


def require_vts(path: Path, dimensions: tuple[int, int, int]) -> None:
    require_file(path)
    text = path.read_text(encoding="utf-8")
    nx, ny, nz = dimensions
    required_fragments = [
        '<?xml version="1.0"?>',
        '<VTKFile type="StructuredGrid"',
        f'WholeExtent="0 {nx - 1} 0 {ny - 1} 0 {nz - 1}"',
        f'Extent="0 {nx - 1} 0 {ny - 1} 0 {nz - 1}"',
        '<PointData Scalars="rho" Vectors="electric">',
        'Name="rho" format="ascii"',
        'Name="phi" format="ascii"',
        'Name="electric" NumberOfComponents="3" format="ascii"',
        'Name="Points" NumberOfComponents="3" format="ascii"',
    ]
    for fragment in required_fragments:
        require(fragment in text, f"missing VTK XML fragment {fragment!r} in {path}")


def require_vtu(path: Path, points: int, cells: int) -> None:
    require_file(path)
    text = path.read_text(encoding="utf-8")
    required_fragments = [
        '<?xml version="1.0"?>',
        '<VTKFile type="UnstructuredGrid"',
        f'NumberOfPoints="{points}" NumberOfCells="{cells}"',
        '<PointData Scalars="rho" Vectors="electric">',
        'Name="rho" format="ascii"',
        'Name="phi" format="ascii"',
        'Name="electric" NumberOfComponents="3" format="ascii"',
        'Name="connectivity" format="ascii"',
        'Name="offsets" format="ascii"',
        'Name="types" format="ascii"',
    ]
    for fragment in required_fragments:
        require(fragment in text, f"missing VTK unstructured fragment {fragment!r} in {path}")


def require_particle_csv(path: Path, expected_header: Sequence[str], min_rows: int = 1) -> None:
    header, rows = read_csv(path)
    require(header == list(expected_header), f"unexpected header in {path}: {header!r}")
    require(len(rows) >= min_rows, f"expected at least {min_rows} data rows in {path}; found {len(rows)}")
    alive_index = header.index("alive")
    for row_index, row in enumerate(rows, start=2):
        for column_index, (column, value) in enumerate(zip(header, row)):
            if column == "species":
                require(value, f"empty species name in {path}:{row_index}")
                continue
            try:
                number = float(value)
            except ValueError as exc:
                raise SmokeFailure(f"non-numeric value in {path}:{row_index} column {column!r}: {value!r}") from exc
            require(math.isfinite(number), f"non-finite value in {path}:{row_index} column {column!r}: {value!r}")
            if column_index == alive_index:
                require(value in {"0", "1", "0.0", "1.0"}, f"invalid alive flag in {path}:{row_index}: {value!r}")


def run_example(cli: Path, config_name: str, output_name: str, temp_root: Path) -> Path:
    copied_config, output_dir = rewrite_output_dir(EXAMPLES / config_name, temp_root, output_name)
    print(f"[smoke] {config_name} -> {output_dir.relative_to(temp_root)}", flush=True)
    subprocess.run([str(cli), str(copied_config)], cwd=ROOT, check=True)
    require(output_dir.is_dir(), f"example did not create output directory: {output_dir}")
    unit_metadata = require_file(output_dir / "units.txt").read_text(
        encoding="utf-8"
    )
    config_text = copied_config.read_text(encoding="utf-8")
    unit_match = re.search(
        r"^\s*units\s*=\s*([A-Za-z]+)\s*$",
        config_text, re.MULTILINE,
    )
    expected_unit_system = (
        unit_match.group(1).lower() if unit_match else "normalized"
    )
    relative_match = re.search(
        r"^\s*relative_permittivity\s*=\s*([^\s#;]+)",
        config_text, re.MULTILINE,
    )
    expected_relative_permittivity = (
        float(relative_match.group(1)) if relative_match else 1.0
    )
    metadata_values = dict(
        line.split(maxsplit=1)
        for line in unit_metadata.splitlines()
        if len(line.split(maxsplit=1)) == 2
    )
    expected_permittivity = expected_relative_permittivity * (
        8.8541878128e-12
        if expected_unit_system == "si"
        else 1.0
    )
    require(
        metadata_values.get("unit_system") == expected_unit_system and
        math.isclose(
            float(metadata_values.get("relative_permittivity", "nan")),
            expected_relative_permittivity,
            rel_tol=1e-14,
        ) and
        math.isclose(
            float(metadata_values.get("permittivity", "nan")),
            expected_permittivity,
            rel_tol=1e-14,
        ),
        f"example unit metadata is missing or inconsistent: {output_dir}",
    )
    return output_dir


def check_two_stream(output_dir: Path) -> None:
    header, rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=["step", "time", "kinetic_energy", "field_energy", "total_energy", "charge_l1", "live_particles", "phi_left", "phi_right"],
        min_rows=2,
    )
    require_step(rows, 300, output_dir / "scalars.csv")
    require_csv(output_dir / "fields_0.csv", expected_header=["x", "rho", "phi", "E"], min_rows=128)
    require_csv(output_dir / "fields_300.csv", expected_header=["x", "rho", "phi", "E"], min_rows=128)


def check_external_state_1d(output_dir: Path) -> None:
    _, rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=[
            "step", "time", "kinetic_energy", "field_energy",
            "total_energy", "charge_l1", "live_particles",
            "phi_left", "phi_right",
        ],
        min_rows=3,
    )
    require_step(rows, 2, output_dir / "scalars.csv")
    header, initialization_rows = read_csv(
        output_dir / "initialization.csv")
    source_index = header.index("state_source")
    require(
        initialization_rows and
        all(row[source_index] == "external"
            for row in initialization_rows),
        "external-state example initialization source is not external",
    )
    state_metadata = require_file(
        output_dir / "initial_state_metadata.txt"
    ).read_text(encoding="utf-8")
    require(
        "format aurorapic_particle_state\n" in state_metadata and
        "version 1\n" in state_metadata and
        re.search(r"^signature [1-9][0-9]*$", state_metadata, re.MULTILINE) and
        "expected_signature 15689689587873937355\n" in state_metadata,
        "external-state example provenance metadata is incomplete",
    )


def check_sheath_steady(output_dir: Path) -> None:
    header, rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=["step", "time", "kinetic_energy", "field_energy", "total_energy", "charge_l1", "live_particles", "phi_left", "phi_right"],
        min_rows=2,
    )
    final_step = int(float(rows[-1][header.index("step")]))
    require(final_step > 0, f"steady sheath did not advance: final step {final_step}")
    require_file(output_dir / "fields_0.csv")
    require(any(path.name.startswith("fields_") and path.suffix == ".csv" for path in output_dir.iterdir()),
            f"steady sheath did not write field CSVs in {output_dir}")
    require_csv(output_dir / "collisions.csv", min_rows=2)


def check_mcc_relaxation(output_dir: Path) -> None:
    initialization_header, initialization_rows = read_csv(
        output_dir / "initialization.csv"
    )
    initialization = dict(zip(
        initialization_header, initialization_rows[0]
    ))
    require(
        float(initialization["thermal_velocity_y"]) > 0.05 and
        float(initialization["thermal_velocity_z"]) > 0.05,
        "1D3V MCC example did not initialize transverse velocity spread",
    )
    scalar_header, scalar_rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=[
            "step", "time", "kinetic_energy", "field_energy",
            "total_energy", "charge_l1", "live_particles",
            "phi_left", "phi_right",
        ],
        min_rows=6,
    )
    require_step(scalar_rows, 50, output_dir / "scalars.csv")
    initial_energy = float(
        scalar_rows[0][scalar_header.index("kinetic_energy")]
    )
    final_energy = float(
        scalar_rows[-1][scalar_header.index("kinetic_energy")]
    )
    require(
        final_energy < initial_energy,
        "MCC excitation case did not reduce kinetic energy",
    )
    collision_header, collision_rows = require_csv(
        output_dir / "collisions.csv", min_rows=6
    )
    require_step(collision_rows, 50, output_dir / "collisions.csv")
    final = {
        name: float(value)
        for name, value in zip(collision_header, collision_rows[-1])
    }
    require(
        final["cumulative_candidates"] > 0 and
        final["cumulative_collisions_elastic"] > 0 and
        final["cumulative_collisions_excitation"] > 0,
        "MCC example did not exercise every configured channel",
    )


def check_mcc_ionization_1d(output_dir: Path) -> None:
    scalar_header, scalar_rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=[
            "step", "time", "kinetic_energy", "field_energy",
            "total_energy", "charge_l1", "live_particles",
            "phi_left", "phi_right",
        ],
        min_rows=5,
    )
    require_step(scalar_rows, 4, output_dir / "scalars.csv")
    initial_energy = float(
        scalar_rows[0][scalar_header.index("total_energy")]
    )
    final_energy = float(
        scalar_rows[-1][scalar_header.index("total_energy")]
    )
    final_live = int(float(
        scalar_rows[-1][scalar_header.index("live_particles")]
    ))
    collision_header, collision_rows = require_csv(
        output_dir / "collisions.csv", min_rows=5
    )
    require_step(collision_rows, 4, output_dir / "collisions.csv")
    final = {
        name: int(float(value))
        for name, value in zip(collision_header, collision_rows[-1])
        if name != "time"
    }
    ionizations = final[
        "cumulative_collisions_electron_mcc.ionization"
    ]
    require(
        ionizations > 0 and
        final["cumulative_collisions_ion_mcc.elastic"] > 0,
        "1D multi-MCC example did not exercise both targets",
    )
    require(
        final_live == 128 + 2 * ionizations,
        "1D ionization did not create charge-paired products",
    )
    require(
        final_energy < initial_energy,
        "1D ionization case did not lose threshold energy",
    )
    require_file(output_dir / "checkpoint_2.apc")


def check_rf_electrode_1d(output_dir: Path) -> None:
    header, rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=[
            "step", "time", "kinetic_energy", "field_energy",
            "total_energy", "charge_l1", "live_particles",
            "phi_left", "phi_right",
        ],
        min_rows=5,
    )
    require_step(rows, 16, output_dir / "scalars.csv")
    by_step = {
        int(float(row[header.index("step")])): {
            name: float(value)
            for name, value in zip(header, row)
        }
        for row in rows
    }
    require(
        abs(by_step[0]["phi_right"]) < 1e-14,
        "RF electrode was not zero-phase at time zero",
    )
    require(
        abs(by_step[8]["phi_right"] - 2.0) < 1e-12,
        "RF electrode did not reach its configured quarter-cycle amplitude",
    )
    require(
        abs(by_step[16]["phi_right"]) < 1e-12,
        "RF electrode did not return to zero at half-cycle",
    )
    require_csv(
        output_dir / "fields_8.csv",
        expected_header=["x", "rho", "phi", "E"],
        min_rows=33,
    )


def check_plasma_2d(output_dir: Path) -> None:
    header, rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=[
            "step", "time", "kinetic_energy", "field_energy", "total_energy", "charge_l1", "live_particles",
            "absorbed_left", "absorbed_right", "absorbed_bottom", "absorbed_top",
            "live_particles_electrons", "live_particles_ions",
        ],
        min_rows=3,
    )
    require_step(rows, 20, output_dir / "scalars.csv")
    require_vtk(output_dir / "fields_0.vtk", (32, 32, 1))
    require_vtk(output_dir / "fields_20.vtk", (32, 32, 1))
    require_particle_csv(
        output_dir / "particles_20.csv",
        expected_header=["species_id", "species", "x", "y", "vx", "vy", "vz", "alive"],
        min_rows=1,
    )


def check_hall_field_profile_smoke(output_dir: Path) -> None:
    _, rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=[
            "step", "time", "kinetic_energy", "field_energy", "total_energy",
            "charge_l1", "live_particles", "absorbed_left",
            "absorbed_right", "absorbed_bottom", "absorbed_top",
            "live_particles_electrons", "live_particles_ions",
        ],
        min_rows=3,
    )
    require_step(rows, 4, output_dir / "scalars.csv")
    source_header, source_rows = read_csv(
        output_dir / "sources.csv")
    require(
        source_header == [
            "step", "time", "source", "macro_pairs_created",
            "represented_pairs_created",
            "fractional_macro_pair_remainder",
            "injected_kinetic_energy",
            "effective_profile_area", "out_of_plane_depth",
            "configured_peak_volumetric_pair_rate",
            "resolved_represented_pair_rate",
        ],
        f"unexpected header in {output_dir / 'sources.csv'}",
    )
    require(len(source_rows) >= 3, "Hall source diagnostics are incomplete")
    final_source = next(
        row for row in source_rows if int(float(row[0])) == 4)
    require(
        final_source[2] == "channel_pair_seed" and
        int(final_source[3]) == 4 and
        abs(float(final_source[5]) - 0.2) < 1e-12 and
        abs(float(final_source[7]) - 4.8e-5) < 1e-16 and
        abs(float(final_source[10]) - 2.5104e19) < 1e6,
        "Hall smoke did not create its configured source pairs",
    )
    current_path = output_dir / "current_source.csv"
    current_header, current_rows = read_csv(current_path)
    require(
        current_header == [
            "step", "time", "species", "monitor_boundary",
            "emission_boundary", "macro_particles_created",
            "represented_particles_created",
            "control_macro_remainder",
            "cumulative_processed_monitored_charge",
            "cumulative_emitted_charge",
            "charge_balance_residual",
            "injected_kinetic_energy",
        ],
        f"unexpected header in {current_path}",
    )
    require(
        len(current_rows) >= 3,
        "Hall current-source diagnostics are incomplete",
    )
    require_numeric_rows(
        [current_header[index] for index in (0, 1, 5, 6, 7, 8, 9, 10, 11)],
        [[row[index] for index in (0, 1, 5, 6, 7, 8, 9, 10, 11)]
         for row in current_rows],
        current_path,
    )
    final_current = next(
        row for row in current_rows if int(float(row[0])) == 4)
    require(
        final_current[2:5] == ["electrons", "left", "right"]
        and int(final_current[5]) == 0
        and abs(float(final_current[10])) < 1e-30,
        "Hall smoke current-control diagnostics are inconsistent",
    )
    potential_path = output_dir / "potential_reference.csv"
    potential_header, potential_rows = read_csv(potential_path)
    require(
        potential_header == [
            "step", "time", "axis", "coordinate", "target",
            "unshifted_line_mean", "applied_offset",
            "corrected_line_mean",
        ],
        f"unexpected header in {potential_path}",
    )
    require(
        len(potential_rows) >= 3,
        "Hall potential-reference diagnostics are incomplete",
    )
    require_numeric_rows(
        [potential_header[index] for index in (0, 1, 3, 4, 5, 6, 7)],
        [[row[index] for index in (0, 1, 3, 4, 5, 6, 7)]
         for row in potential_rows],
        potential_path,
    )
    final_potential = next(
        row for row in potential_rows
        if int(float(row[0])) == 4)
    require(
        final_potential[2] == "x"
        and abs(float(final_potential[3]) - 0.024) < 1e-15
        and abs(float(final_potential[7])) < 1e-12,
        "Hall smoke potential reference missed its target",
    )
    field_path = output_dir / "resolved_field_profiles.csv"
    field_header, field_rows = read_csv(field_path)
    require(
        field_header == [
            "step", "time", "profile_axis", "coordinate",
            "potential", "electric_x", "electric_y",
            "charge_density",
        ] and len(field_rows) == 5 * 16,
        "Hall resolved field profiles are incomplete",
    )
    require(
        all(row[2] == "x" for row in field_rows),
        "Hall field profile used the wrong axis",
    )
    require_numeric_rows(
        [field_header[index] for index in (0, 1, 3, 4, 5, 6, 7)],
        [[row[index] for index in (0, 1, 3, 4, 5, 6, 7)]
         for row in field_rows],
        field_path,
    )
    require_step(field_rows, 4, field_path)

    species_path = output_dir / "resolved_species_profiles.csv"
    species_header, species_rows = read_csv(species_path)
    require(
        species_header == [
            "step", "time", "profile_axis", "coordinate", "species",
            "macro_particle_equivalent", "represented_number",
            "number_density", "mean_velocity_x", "mean_velocity_y",
            "mean_velocity_z", "thermal_speed_x", "thermal_speed_y",
            "thermal_speed_z", "temperature_ev", "current_density_x",
            "current_density_y", "current_density_z",
        ] and len(species_rows) == 5 * 16 * 2,
        "Hall resolved species profiles are incomplete",
    )
    require(
        {row[4] for row in species_rows} == {"electrons", "ions"}
        and all(row[2] == "x" for row in species_rows),
        "Hall species profiles lost their axis or species identity",
    )
    require_numeric_rows(
        [species_header[index]
         for index in (0, 1, 3, *range(5, 18))],
        [[row[index]
          for index in (0, 1, 3, *range(5, 18))]
         for row in species_rows],
        species_path,
    )
    require_step(species_rows, 4, species_path)

    mode_path = output_dir / "resolved_modes.csv"
    mode_header, mode_rows = read_csv(mode_path)
    require(
        mode_header == [
            "step", "time", "mode_axis", "mode", "wavenumber",
            "quantity", "species", "real", "imaginary", "amplitude",
        ] and len(mode_rows) == 5 * 4 * 11,
        "Hall resolved mode history is incomplete",
    )
    require(
        all(row[2] == "y" for row in mode_rows)
        and {int(row[3]) for row in mode_rows} == {0, 1, 2, 3}
        and {"charge_density", "electric_x", "electric_y",
             "number_density", "current_x", "current_y", "current_z"}
            == {row[5] for row in mode_rows},
        "Hall mode history lost its axis, modes, or quantities",
    )
    require_numeric_rows(
        [mode_header[index] for index in (0, 1, 3, 4, 7, 8, 9)],
        [[row[index] for index in (0, 1, 3, 4, 7, 8, 9)]
         for row in mode_rows],
        mode_path,
    )
    require_step(mode_rows, 4, mode_path)

    for name, expected_rows in (
        ("resolved_field_time_average.csv", 16),
        ("resolved_species_time_average.csv", 32),
    ):
        average_path = output_dir / name
        average_header, average_rows = read_csv(average_path)
        require(
            len(average_rows) == expected_rows,
            f"Hall time average has the wrong row count: {average_path}",
        )
        require(
            all(
                abs(float(row[0])) < 1e-30
                and abs(float(row[1]) - 2e-11) < 1e-25
                and abs(float(row[2]) - 2e-11) < 1e-25
                and int(row[3]) == 5
                and row[4] == "x"
                for row in average_rows
            ),
            f"Hall time-average window metadata is wrong: {average_path}",
        )


def check_electrode_2d(output_dir: Path) -> None:
    header, rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=[
            "step", "time", "kinetic_energy", "field_energy", "total_energy", "charge_l1", "live_particles",
            "absorbed_left", "absorbed_right", "absorbed_bottom", "absorbed_top",
            "live_particles_electrons", "live_particles_ions",
        ],
        min_rows=3,
    )
    require_step(rows, 10, output_dir / "scalars.csv")
    require_vtk(output_dir / "fields_0.vtk", (32, 24, 1))
    require_vtk(output_dir / "fields_10.vtk", (32, 24, 1))


def check_plasma_3d(output_dir: Path) -> None:
    header, rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=[
            "step", "time", "kinetic_energy", "field_energy", "total_energy", "charge_l1", "live_particles",
            "absorbed_left", "absorbed_right", "absorbed_bottom", "absorbed_top", "absorbed_back", "absorbed_front",
            "live_particles_electrons", "live_particles_ions",
        ],
        min_rows=4,
    )
    require_step(rows, 3, output_dir / "scalars.csv")
    require_vtk(output_dir / "fields_0.vtk", (8, 8, 8))
    require_vtk(output_dir / "fields_3.vtk", (8, 8, 8))
    require_vts(output_dir / "fields_0.vts", (8, 8, 8))
    require_vts(output_dir / "fields_3.vts", (8, 8, 8))
    require_particle_csv(
        output_dir / "particles_3.csv",
        expected_header=["species_id", "species", "x", "y", "z", "vx", "vy", "vz", "alive"],
        min_rows=1,
    )


def check_imported_plasma_2d(output_dir: Path) -> None:
    header, rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=[
            "step", "time", "kinetic_energy", "field_energy", "total_energy",
            "charge_l1", "live_particles", "poisson_iterations",
            "poisson_initial_residual", "poisson_final_residual",
            "particle_seconds", "deposition_seconds", "field_solve_seconds",
            "location_cache_hits", "location_searches",
            "absorbed_electrode", "absorbed_inlet", "absorbed_outlet", "absorbed_wall",
            "injected_electron_inlet", "injected_ion_inlet",
            "emitted_outlet_secondary_electrons",
            "impact_macroparticles_electrons@outlet",
            "impact_physical_particles_electrons@outlet",
            "impact_charge_electrons@outlet",
            "impact_kinetic_energy_electrons@outlet",
            "impact_rate_electrons@outlet",
            "impact_flux_electrons@outlet",
            "impact_macroparticles_ions@outlet",
            "impact_physical_particles_ions@outlet",
            "impact_charge_ions@outlet",
            "impact_kinetic_energy_ions@outlet",
            "impact_rate_ions@outlet",
            "impact_flux_ions@outlet",
        ],
        min_rows=4,
    )
    require_step(rows, 3, output_dir / "scalars.csv")
    require_vtu(output_dir / "fields_0.vtu", 6, 3)
    require_vtu(output_dir / "fields_3.vtu", 6, 3)
    require_particle_csv(
        output_dir / "particles_3.csv",
        expected_header=["species_id", "species", "x", "y", "vx", "vy", "vz", "alive"],
        min_rows=1,
    )
    require_file(output_dir / "checkpoint_3.apc")


def check_imported_mcc_2d(output_dir: Path) -> None:
    _, scalar_rows = require_csv(
        output_dir / "scalars.csv", min_rows=7
    )
    require_step(scalar_rows, 6, output_dir / "scalars.csv")
    collision_header, collision_rows = require_csv(
        output_dir / "collisions.csv", min_rows=7
    )
    require_step(collision_rows, 6, output_dir / "collisions.csv")
    final = {
        name: float(value)
        for name, value in zip(collision_header, collision_rows[-1])
    }
    require(
        final["cumulative_candidates"] == 44 and
        final["cumulative_null_collisions"] == 20 and
        final["cumulative_synthetic_elastic"] == 24,
        "imported 2D3V MCC deterministic collision envelope changed",
    )
    require_file(output_dir / "checkpoint_3.apc")
    require_file(output_dir / "checkpoint_6.apc")


def check_imported_ionization_2d(output_dir: Path) -> None:
    scalar_header, scalar_rows = require_csv(
        output_dir / "scalars.csv", min_rows=5
    )
    require_step(scalar_rows, 4, output_dir / "scalars.csv")
    final_scalar = {
        name: float(value)
        for name, value in zip(scalar_header, scalar_rows[-1])
    }
    collision_header, collision_rows = require_csv(
        output_dir / "collisions.csv", min_rows=5
    )
    require_step(collision_rows, 4, output_dir / "collisions.csv")
    final_collision = {
        name: float(value)
        for name, value in zip(collision_header, collision_rows[-1])
    }
    require(
        final_collision["cumulative_candidates"] == 5 and
        final_collision["cumulative_null_collisions"] == 2 and
        final_collision["cumulative_synthetic_ionization"] == 3,
        "imported ionization deterministic collision envelope changed",
    )
    require(
        final_scalar["live_particles"] == 38 and
        abs(final_scalar["total_energy"] - 125.0) < 1e-11,
        "imported ionization product or energy envelope changed",
    )
    require_file(output_dir / "checkpoint_2.apc")
    require_file(output_dir / "checkpoint_4.apc")
    metadata = require_file(
        output_dir / "collision_data.txt"
    ).read_text(encoding="utf-8")
    require(
        'dataset_id "aurorapic.synthetic.ionization"' in metadata and
        'retrieved "2026-07-28"' in metadata and
        'channel "synthetic_ionization" "ionization"' in metadata,
        "imported ionization gas metadata output is incomplete",
    )


def check_imported_attachment_2d(output_dir: Path) -> None:
    scalar_header, scalar_rows = require_csv(
        output_dir / "scalars.csv", min_rows=5
    )
    require_step(scalar_rows, 4, output_dir / "scalars.csv")
    final_scalar = {
        name: float(value)
        for name, value in zip(scalar_header, scalar_rows[-1])
    }
    collision_header, collision_rows = require_csv(
        output_dir / "collisions.csv", min_rows=5
    )
    require_step(collision_rows, 4, output_dir / "collisions.csv")
    final_collision = {
        name: float(value)
        for name, value in zip(collision_header, collision_rows[-1])
    }
    require(
        final_collision["cumulative_candidates"] == 6 and
        final_collision["cumulative_null_collisions"] == 2 and
        final_collision["cumulative_synthetic_attachment"] == 4,
        "imported attachment deterministic collision envelope changed",
    )
    require(
        final_scalar["live_particles"] == 32 and
        abs(final_scalar["total_energy"] - 96.0) < 1e-11,
        "imported attachment product or energy envelope changed",
    )
    require_file(output_dir / "checkpoint_2.apc")
    require_file(output_dir / "checkpoint_4.apc")


def check_imported_charge_exchange_2d(output_dir: Path) -> None:
    scalar_header, scalar_rows = require_csv(
        output_dir / "scalars.csv", min_rows=7
    )
    require_step(scalar_rows, 6, output_dir / "scalars.csv")
    final_scalar = {
        name: float(value)
        for name, value in zip(scalar_header, scalar_rows[-1])
    }
    collision_header, collision_rows = require_csv(
        output_dir / "collisions.csv", min_rows=7
    )
    require_step(collision_rows, 6, output_dir / "collisions.csv")
    final_collision = {
        name: float(value)
        for name, value in zip(collision_header, collision_rows[-1])
    }
    require(
        final_collision["cumulative_candidates"] == 15 and
        final_collision["cumulative_null_collisions"] == 8 and
        final_collision["cumulative_synthetic_charge_exchange"] == 7,
        "imported charge-exchange deterministic envelope changed",
    )
    require(
        final_scalar["live_particles"] == 64 and
        abs(final_scalar["total_energy"] - 500.0) < 1e-10,
        "imported charge-exchange particle or energy envelope changed",
    )
    require_file(output_dir / "checkpoint_3.apc")
    require_file(output_dir / "checkpoint_6.apc")
    metadata = require_file(
        output_dir / "collision_data.txt"
    ).read_text(encoding="utf-8")
    require(
        'dataset_id "aurorapic.synthetic.charge_exchange"' in metadata and
        '"charge_exchange"' in metadata,
        "imported charge-exchange metadata output is incomplete",
    )


def check_biased_probe_2d(output_dir: Path) -> None:
    header, rows = require_csv(
        output_dir / "scalars.csv",
        min_rows=5,
    )
    required_columns = {
        "absorbed_outlet",
        "absorbed_probe",
        "absorbed_wall",
        "injected_electron_inlet",
        "injected_ion_inlet",
        "emitted_probe_secondary_electrons",
        "impact_flux_electrons@probe",
        "impact_flux_ions@probe",
    }
    require(
        required_columns.issubset(header),
        f"biased-probe diagnostics are missing columns: {sorted(required_columns - set(header))}",
    )
    require_step(rows, 20, output_dir / "scalars.csv")
    final = {name: float(value) for name, value in zip(header, rows[-1])}
    require(350 <= final["live_particles"] <= 450,
            f"biased-probe live-particle envelope changed: {final['live_particles']}")
    require(final["absorbed_outlet"] >= 250,
            "biased-probe case did not produce an outlet flux")
    require(final["absorbed_probe"] >= 80,
            "biased-probe case did not collect particles on the internal probe")
    require(final["absorbed_wall"] > 0,
            "biased-probe case did not exercise chamber-wall absorption")
    require(final["injected_electron_inlet"] == 40 and
            final["injected_ion_inlet"] == 40,
            "biased-probe inlet source schedule changed")
    require(final["emitted_probe_secondary_electrons"] > 0,
            "biased-probe case did not exercise secondary emission")
    require(final["poisson_final_residual"] < 1e-6,
            "biased-probe field solve exceeded its residual envelope")
    require_vtu(output_dir / "fields_0.vtu", 725, 1342)
    require_vtu(output_dir / "fields_20.vtu", 725, 1342)
    require_particle_csv(
        output_dir / "particles_20.csv",
        expected_header=["species_id", "species", "x", "y", "vx", "vy", "vz", "alive"],
        min_rows=1,
    )
    require_file(output_dir / "checkpoint_20.apc")


def run_smokes(cli: Path, temp_root: Path) -> None:
    checks = [
        ("external_state_1d.cfg", "external_state_1d", check_external_state_1d),
        ("two_stream.cfg", "two_stream", check_two_stream),
        ("sheath_steady.cfg", "sheath_steady", check_sheath_steady),
        ("mcc_relaxation.cfg", "mcc_relaxation", check_mcc_relaxation),
        ("mcc_ionization_1d.cfg", "mcc_ionization_1d", check_mcc_ionization_1d),
        ("rf_electrode_1d.cfg", "rf_electrode_1d", check_rf_electrode_1d),
        ("plasma_2d.cfg", "plasma_2d", check_plasma_2d),
        ("hall_field_profile_smoke.cfg", "hall_field_profile_smoke",
         check_hall_field_profile_smoke),
        ("electrode_2d.cfg", "electrode_2d", check_electrode_2d),
        ("imported_plasma_2d.cfg", "imported_plasma_2d", check_imported_plasma_2d),
        ("imported_mcc_2d.cfg", "imported_mcc_2d", check_imported_mcc_2d),
        ("imported_ionization_2d.cfg", "imported_ionization_2d", check_imported_ionization_2d),
        ("imported_attachment_2d.cfg", "imported_attachment_2d", check_imported_attachment_2d),
        ("imported_charge_exchange_2d.cfg", "imported_charge_exchange_2d", check_imported_charge_exchange_2d),
        ("biased_probe_2d.cfg", "biased_probe_2d", check_biased_probe_2d),
        ("plasma_3d.cfg", "plasma_3d", check_plasma_3d),
    ]
    for config_name, output_name, check in checks:
        output_dir = run_example(cli, config_name, output_name, temp_root)
        require_initialization_audits(output_dir)
        check(output_dir)


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cli", nargs="?", default=ROOT / "build" / "aurorapic_cli", type=Path,
                        help="Path to aurorapic_cli (default: build/aurorapic_cli)")
    parser.add_argument("--keep-output", action="store_true",
                        help="Keep the temporary smoke output directory after the run")
    return parser.parse_args(list(argv))


def main(argv: Iterable[str] = sys.argv[1:]) -> int:
    args = parse_args(argv)
    cli = args.cli if args.cli.is_absolute() else ROOT / args.cli
    require_file(cli)

    keep_output = args.keep_output or os.environ.get("KEEP_VERIFY_OUTPUTS") == "1"
    temp_parent = ROOT / "test_output_aurorapic_verify"
    temp_parent.mkdir(parents=True, exist_ok=True)
    temp_root = Path(tempfile.mkdtemp(prefix="aurorapic_verify_", dir=temp_parent))
    try:
        run_smokes(cli, temp_root)
        print(f"[smoke] example diagnostics verified under {temp_root}")
        return 0
    except (SmokeFailure, subprocess.CalledProcessError) as exc:
        print(f"[smoke] FAILED: {exc}", file=sys.stderr)
        print(f"[smoke] retained output for debugging: {temp_root}", file=sys.stderr)
        keep_output = True
        return 1
    finally:
        if keep_output:
            print(f"[smoke] kept output: {temp_root}")
        else:
            shutil.rmtree(temp_root, ignore_errors=True)
            parent = temp_root.parent
            try:
                parent.rmdir()
            except OSError:
                pass


if __name__ == "__main__":
    raise SystemExit(main())
