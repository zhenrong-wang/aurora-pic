#!/usr/bin/env python3
"""Estimate and gate a published-scale Hall PIC campaign without launching it."""

from __future__ import annotations

import argparse
import configparser
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile


class HallPreflightError(RuntimeError):
    """Raised when a Hall production estimate is invalid."""


def positive(value: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise argparse.ArgumentTypeError("value must be positive and finite")
    return result


def nonnegative_integer(value: str) -> int:
    result = int(value)
    if result < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return result


def positive_integer(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return result


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise HallPreflightError(f"cannot hash {path}: {error}") from error


def load_case(path: Path) -> tuple[
    configparser.ConfigParser,
    configparser.SectionProxy,
    configparser.SectionProxy,
]:
    parser = configparser.ConfigParser(
        interpolation=None,
        strict=True,
        empty_lines_in_values=False,
    )
    try:
        parser.read_string(
            "[global]\n" + path.read_text(encoding="utf-8")
        )
    except (OSError, UnicodeError, configparser.Error) as error:
        raise HallPreflightError(
            f"cannot read Hall case manifest {path}: {error}"
        ) from error
    for section in ("global", "reference"):
        if section not in parser:
            raise HallPreflightError(
                f"Hall case manifest is missing [{section}]"
            )
    return parser, parser["global"], parser["reference"]


def required(
    section: configparser.SectionProxy,
    key: str,
) -> str:
    if key not in section or not section[key].strip():
        raise HallPreflightError(
            f"[{section.name}] requires {key!r}"
        )
    return section[key].strip()


def checked_int(
    section: configparser.SectionProxy,
    key: str,
) -> int:
    try:
        result = section.getint(key)
    except ValueError as error:
        raise HallPreflightError(
            f"[{section.name}] {key} must be an integer"
        ) from error
    if result <= 0:
        raise HallPreflightError(
            f"[{section.name}] {key} must be positive"
        )
    return result


def checked_float(
    section: configparser.SectionProxy,
    key: str,
) -> float:
    try:
        result = section.getfloat(key)
    except ValueError as error:
        raise HallPreflightError(
            f"[{section.name}] {key} must be numeric"
        ) from error
    if not math.isfinite(result) or result <= 0.0:
        raise HallPreflightError(
            f"[{section.name}] {key} must be positive and finite"
        )
    return result


def gibibytes(value: int | float) -> float:
    return float(value) / (1024.0 ** 3)


def write_json_atomic(
    path: Path,
    report: dict[str, object],
    overwrite: bool,
) -> None:
    if path.exists() and not overwrite:
        raise HallPreflightError(
            f"preflight report already exists: {path}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise


def estimate(args: argparse.Namespace) -> dict[str, object]:
    parser, global_section, reference = load_case(args.case_manifest)
    case_id = required(global_section, "case_id")
    nx = checked_int(reference, "production_nx")
    ny = checked_int(reference, "production_ny")
    steps = checked_int(reference, "production_steps")
    timestep = checked_float(reference, "production_dt_s")
    if args.max_mode > ny // 2:
        raise HallPreflightError(
            "max_mode exceeds the production azimuthal Nyquist limit"
        )
    if args.start_step > steps:
        raise HallPreflightError(
            "diagnostic start_step exceeds production steps"
        )
    if args.particles_per_cell > 10_000_000:
        raise HallPreflightError(
            "particles_per_cell exceeds the bounded preflight limit"
        )
    if args.species > 16:
        raise HallPreflightError(
            "species exceeds the bounded preflight limit"
        )

    cells = nx * ny
    initial_particles = (
        cells * args.particles_per_cell * args.species
    )
    particle_capacity = math.ceil(
        initial_particles * args.particle_capacity_factor
    )
    particle_memory = particle_capacity * args.particle_bytes
    mesh_points = cells
    field_memory = (
        mesh_points * args.field_arrays * 8
    )
    diagnostic_grid_memory = (
        (nx * (5 + 16 * args.species))
        + (ny * 4 * args.species)
        + ((args.max_mode + 1) * (3 + 4 * args.species) * 3)
    ) * 8
    estimated_memory = (
        particle_memory + field_memory + diagnostic_grid_memory
    )

    snapshots = (
        (steps - args.start_step) // args.diagnostic_interval + 1
    )
    field_profile_rows = snapshots * nx
    species_profile_rows = snapshots * nx * args.species
    mode_quantities = 3 + 4 * args.species
    mode_rows = (
        snapshots * (args.max_mode + 1) * mode_quantities
    )
    average_rows = nx * (1 + args.species)
    diagnostic_storage = (
        field_profile_rows * args.profile_row_bytes
        + species_profile_rows * args.species_row_bytes
        + mode_rows * args.mode_row_bytes
        + average_rows * args.average_row_bytes
    )
    checkpoint_bytes = particle_memory + field_memory
    checkpoint_storage = (
        checkpoint_bytes * args.retained_checkpoints
    )
    estimated_storage = diagnostic_storage + checkpoint_storage
    particle_updates_lower_bound = initial_particles * steps
    estimated_wall_seconds = (
        particle_updates_lower_bound / args.particle_push_rate
        if args.particle_push_rate is not None else None
    )

    memory_budget = args.memory_budget_gib * 1024.0 ** 3
    storage_budget = args.storage_budget_gib * 1024.0 ** 3
    checks = {
        "memory_budget": estimated_memory <= memory_budget,
        "storage_budget": estimated_storage <= storage_budget,
        "mode_nyquist": args.max_mode <= ny // 2,
        "bounded_inputs": True,
    }
    within_budgets = all(checks.values())
    reduced_status = global_section.get("status", "")
    missing_physics = (
        parser["reduced_contract"].get("missing_physics", "").split(",")
        if "reduced_contract" in parser else []
    )
    missing_physics = [
        item.strip() for item in missing_physics if item.strip()
    ]
    production_scale = (
        cells >= 100_000
        or initial_particles >= 10_000_000
        or steps >= 1_000_000
    )
    warnings = [
        "This tool estimates resources and never launches AuroraPIC.",
        "Particle-update work is a lower bound because pair and cathode "
        "sources can increase the live population.",
        "CSV estimates exclude VTK, particle dumps, logs, external "
        "reference data, filesystem replication, and temporary files.",
        "Checkpoint estimates include particle capacity and mesh fields "
        "but exclude metadata, compression effects, and temporary files.",
        "Measured throughput on the intended backend is required before "
        "using the optional wall-time estimate for scheduling.",
    ]
    blockers = []
    if reduced_status != "production_ready":
        blockers.append(
            f"case manifest status is {reduced_status!r}, not production_ready"
        )
    if missing_physics:
        blockers.append(
            "case manifest still lists: " + ",".join(missing_physics)
        )
    blockers.extend(
        [
            "no production input deck is generated by this preflight",
            "the mixed spectral-tridiagonal Poisson solver is single-rank "
            "and has not been production-grid qualified",
            "MPI decomposition is not yet available",
        ]
    )
    return {
        "schema_version": 1,
        "case_id": case_id,
        "case_manifest": str(args.case_manifest.resolve()),
        "case_manifest_sha256": sha256(args.case_manifest),
        "production_scale": production_scale,
        "within_declared_budgets": within_budgets,
        "launch_authorized": False,
        "checks": checks,
        "published_contract": {
            "nx": nx,
            "ny": ny,
            "cells": cells,
            "steps": steps,
            "timestep_s": timestep,
            "final_time_s": steps * timestep,
        },
        "assumptions": {
            "species": args.species,
            "particles_per_cell_per_species": args.particles_per_cell,
            "particle_capacity_factor": args.particle_capacity_factor,
            "particle_bytes": args.particle_bytes,
            "field_arrays": args.field_arrays,
            "diagnostic_interval_steps": args.diagnostic_interval,
            "diagnostic_start_step": args.start_step,
            "max_mode": args.max_mode,
            "retained_checkpoints": args.retained_checkpoints,
            "memory_budget_gib": args.memory_budget_gib,
            "storage_budget_gib": args.storage_budget_gib,
            "particle_push_rate_per_second": args.particle_push_rate,
        },
        "estimates": {
            "initial_macroparticles": initial_particles,
            "particle_capacity": particle_capacity,
            "particle_memory_bytes": particle_memory,
            "field_memory_bytes": field_memory,
            "diagnostic_working_memory_bytes": diagnostic_grid_memory,
            "total_memory_bytes": estimated_memory,
            "total_memory_gib": gibibytes(estimated_memory),
            "diagnostic_snapshots": snapshots,
            "field_profile_rows": field_profile_rows,
            "species_profile_rows": species_profile_rows,
            "mode_rows": mode_rows,
            "diagnostic_storage_bytes": diagnostic_storage,
            "checkpoint_bytes_each": checkpoint_bytes,
            "checkpoint_storage_bytes": checkpoint_storage,
            "total_storage_bytes": estimated_storage,
            "total_storage_gib": gibibytes(estimated_storage),
            "particle_updates_lower_bound": particle_updates_lower_bound,
            "estimated_wall_seconds": estimated_wall_seconds,
        },
        "warnings": warnings,
        "readiness_blockers": blockers,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate a published Hall campaign and gate declared memory/"
            "storage budgets without launching a simulation"
        )
    )
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--particles-per-cell", type=positive_integer, default=75
    )
    parser.add_argument("--species", type=positive_integer, default=2)
    parser.add_argument(
        "--particle-capacity-factor", type=positive, default=1.5
    )
    parser.add_argument(
        "--particle-bytes", type=positive_integer, default=96
    )
    parser.add_argument(
        "--field-arrays", type=positive_integer, default=16
    )
    parser.add_argument(
        "--diagnostic-interval", type=positive_integer, default=5000
    )
    parser.add_argument(
        "--start-step", type=nonnegative_integer, default=0
    )
    parser.add_argument(
        "--max-mode", type=nonnegative_integer, default=128
    )
    parser.add_argument(
        "--retained-checkpoints", type=nonnegative_integer, default=2
    )
    parser.add_argument(
        "--memory-budget-gib", type=positive, default=8.0
    )
    parser.add_argument(
        "--storage-budget-gib", type=positive, default=16.0
    )
    parser.add_argument(
        "--particle-push-rate",
        type=positive,
        default=None,
        help="Measured aggregate particle updates per second",
    )
    parser.add_argument(
        "--profile-row-bytes", type=positive_integer, default=180
    )
    parser.add_argument(
        "--species-row-bytes", type=positive_integer, default=360
    )
    parser.add_argument(
        "--mode-row-bytes", type=positive_integer, default=220
    )
    parser.add_argument(
        "--average-row-bytes", type=positive_integer, default=360
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = estimate(args)
        write_json_atomic(args.report, report, args.overwrite)
    except HallPreflightError as error:
        print(f"Hall preflight error: {error}", file=sys.stderr)
        return 2
    if not report["within_declared_budgets"]:
        print(
            "Hall preflight exceeded declared resource budgets",
            file=sys.stderr,
        )
        return 1
    print(
        "Hall preflight passed declared budgets; launch remains unauthorized: "
        f"{args.report}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
