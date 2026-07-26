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
    return output_dir


def check_two_stream(output_dir: Path) -> None:
    header, rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=["step", "time", "kinetic_energy", "field_energy", "total_energy", "charge_l1", "live_particles"],
        min_rows=2,
    )
    require_step(rows, 300, output_dir / "scalars.csv")
    require_csv(output_dir / "fields_0.csv", expected_header=["x", "rho", "phi", "E"], min_rows=128)
    require_csv(output_dir / "fields_300.csv", expected_header=["x", "rho", "phi", "E"], min_rows=128)


def check_sheath_steady(output_dir: Path) -> None:
    header, rows = require_csv(
        output_dir / "scalars.csv",
        expected_header=["step", "time", "kinetic_energy", "field_energy", "total_energy", "charge_l1", "live_particles"],
        min_rows=2,
    )
    final_step = int(float(rows[-1][header.index("step")]))
    require(final_step > 0, f"steady sheath did not advance: final step {final_step}")
    require_file(output_dir / "fields_0.csv")
    require(any(path.name.startswith("fields_") and path.suffix == ".csv" for path in output_dir.iterdir()),
            f"steady sheath did not write field CSVs in {output_dir}")


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
        expected_header=["species_id", "species", "x", "y", "vx", "vy", "alive"],
        min_rows=1,
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
    require_particle_csv(
        output_dir / "particles_3.csv",
        expected_header=["species_id", "species", "x", "y", "z", "vx", "vy", "vz", "alive"],
        min_rows=1,
    )


def run_smokes(cli: Path, temp_root: Path) -> None:
    checks = [
        ("two_stream.cfg", "two_stream", check_two_stream),
        ("sheath_steady.cfg", "sheath_steady", check_sheath_steady),
        ("plasma_2d.cfg", "plasma_2d", check_plasma_2d),
        ("electrode_2d.cfg", "electrode_2d", check_electrode_2d),
        ("plasma_3d.cfg", "plasma_3d", check_plasma_3d),
    ]
    for config_name, output_name, check in checks:
        output_dir = run_example(cli, config_name, output_name, temp_root)
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
