#!/usr/bin/env python3
"""Audit a Turner diagnostic window's source, population, and wall balance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile


class BalanceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise BalanceError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_csv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            require(reader.fieldnames is not None, f"{path} has no header")
            rows = list(reader)
    except OSError as error:
        raise BalanceError(f"cannot read {path}: {error}") from error
    require(len(rows) >= 2, f"{path} must contain at least two samples")
    return rows


def integer(row: dict[str, str], name: str, source: str) -> int:
    require(name in row, f"{source} is missing column {name!r}")
    try:
        value = int(row[name])
    except ValueError as error:
        raise BalanceError(
            f"{source} column {name!r} is not an integer"
        ) from error
    require(value >= 0, f"{source} column {name!r} is negative")
    return value


def number(row: dict[str, str], name: str, source: str) -> float:
    require(name in row, f"{source} is missing column {name!r}")
    try:
        value = float(row[name])
    except ValueError as error:
        raise BalanceError(
            f"{source} column {name!r} is not numeric"
        ) from error
    require(math.isfinite(value), f"{source} column {name!r} is non-finite")
    return value


def write_report(path: Path, report: dict[str, object]) -> None:
    require(not path.exists(), f"refusing to overwrite report: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def analyze(
    scalars_path: Path,
    collisions_path: Path,
    boundary_path: Path,
    power_path: Path | None,
    electron_species: str,
    ion_species: str,
    ionization_channel: str,
    expected_steps: int,
    reported_ion_current: float,
    reported_electron_power: float,
    reported_ion_power: float,
) -> dict[str, object]:
    scalars = load_csv(scalars_path)
    collisions = load_csv(collisions_path)
    boundaries = load_csv(boundary_path)
    power = load_csv(power_path) if power_path is not None else None
    diagnostics = [
        (scalars, "scalars"),
        (collisions, "collisions"),
        (boundaries, "boundary losses"),
    ]
    if power is not None:
        diagnostics.append((power, "power transfer"))
    starts = [
        integer(rows[0], "step", label)
        for rows, label in diagnostics
    ]
    ends = [
        integer(rows[-1], "step", label)
        for rows, label in diagnostics
    ]
    require(len(set(starts)) == 1, "diagnostic start steps do not agree")
    require(len(set(ends)) == 1, "diagnostic end steps do not agree")
    start_step, end_step = starts[0], ends[0]
    require(end_step - start_step == expected_steps,
            "diagnostic window has the wrong step count")
    start_time = number(scalars[0], "time", "scalars")
    end_time = number(scalars[-1], "time", "scalars")
    duration = end_time - start_time
    require(duration > 0.0, "diagnostic duration must be positive")
    counter_origins = [
        integer(row, "counter_origin_step", "boundary losses")
        for row in (boundaries[0], boundaries[-1])
    ]
    require(
        counter_origins[0] == counter_origins[1]
        and counter_origins[0] <= start_step,
        "boundary counters do not cover the entire diagnostic window",
    )
    counter_origin = counter_origins[0]
    power_origin: int | None = None
    if power is not None:
        power_origins = [
            integer(row, "counter_origin_step", "power transfer")
            for row in (power[0], power[-1])
        ]
        require(
            power_origins[0] == power_origins[1]
            and power_origins[0] <= start_step,
            "power counters do not cover the entire diagnostic window",
        )
        power_origin = power_origins[0]

    cumulative_ionization = (
        "cumulative_collisions_" + ionization_channel
    )
    ionizations = (
        integer(collisions[-1], cumulative_ionization, "collisions")
        - integer(collisions[0], cumulative_ionization, "collisions")
    )
    require(ionizations >= 0, "cumulative ionization count decreased")

    species_report: dict[str, object] = {}
    electrode_report: dict[str, dict[str, object]] = {
        "left": {},
        "right": {},
    }
    for species in (electron_species, ion_species):
        live_column = "live_particles_" + species
        initial = integer(scalars[0], live_column, "scalars")
        final = integer(scalars[-1], live_column, "scalars")
        losses: dict[str, int] = {}
        for side in ("left", "right"):
            count_column = f"absorbed_{side}_count_{species}"
            initial_loss = integer(
                boundaries[0], count_column, "boundary losses"
            )
            final_loss = integer(
                boundaries[-1], count_column, "boundary losses"
            )
            losses[side] = final_loss - initial_loss
            require(losses[side] >= 0, "boundary loss count decreased")

            charge_column = (
                f"absorbed_{side}_charge_{species}_C_m-2"
            )
            energy_column = (
                f"absorbed_{side}_kinetic_energy_{species}_J_m-2"
            )
            charge = (
                number(
                    boundaries[-1], charge_column, "boundary losses"
                )
                - number(
                    boundaries[0], charge_column, "boundary losses"
                )
            )
            energy = (
                number(
                    boundaries[-1], energy_column, "boundary losses"
                )
                - number(
                    boundaries[0], energy_column, "boundary losses"
                )
            )
            require(energy >= 0.0, "absorbed kinetic energy decreased")
            electrode_report[side][species] = {
                "absorbed_macro_particles": losses[side],
                "absorbed_charge_C_m2": charge,
                "absorbed_charge_rate_A_m2": charge / duration,
                "absorbed_kinetic_energy_J_m2": energy,
                "wall_kinetic_power_W_m2": energy / duration,
            }

        expected_final = (
            initial + ionizations - losses["left"] - losses["right"]
        )
        residual = final - expected_final
        species_report[species] = {
            "initial_macro_particles": initial,
            "final_macro_particles": final,
            "change_macro_particles": final - initial,
            "ionization_products": ionizations,
            "absorbed_left": losses["left"],
            "absorbed_right": losses["right"],
            "balance_residual_macro_particles": residual,
            "balance_exact": residual == 0,
        }

    ion_currents = [
        abs(float(electrode_report[side][ion_species][
            "absorbed_charge_rate_A_m2"
        ]))
        for side in ("left", "right")
    ]
    mean_ion_current = math.fsum(ion_currents) / 2.0
    exact = all(
        bool(species_report[name]["balance_exact"])
        for name in (electron_species, ion_species)
    )
    report: dict[str, object] = {
        "turner_balance_version": 2,
        "scope": "post_benchmark_diagnostic_window",
        "window": {
            "start_step": start_step,
            "end_step": end_step,
            "steps": end_step - start_step,
            "start_time_s": start_time,
            "end_time_s": end_time,
            "duration_s": duration,
            "counter_origin_step": counter_origin,
        },
        "provenance": {
            "scalars": {
                "path": str(scalars_path.resolve()),
                "sha256": sha256(scalars_path),
            },
            "collisions": {
                "path": str(collisions_path.resolve()),
                "sha256": sha256(collisions_path),
            },
            "boundary_losses": {
                "path": str(boundary_path.resolve()),
                "sha256": sha256(boundary_path),
            },
        },
        "ionization_macro_events": ionizations,
        "species_balance": species_report,
        "electrodes": electrode_report,
        "ion_current_context": {
            "left_magnitude_A_m2": ion_currents[0],
            "right_magnitude_A_m2": ion_currents[1],
            "two_electrode_mean_magnitude_A_m2": mean_ion_current,
            "turner_table_iii_A_m2": reported_ion_current,
            "relative_difference":
                mean_ion_current / reported_ion_current - 1.0,
            "acceptance_gate": "none_diagnostic_context_only",
        },
        "balance_exact": exact,
        "wall_power_scope":
            "absorbed_particle_kinetic_energy_not_volume_power_transfer",
        "physics_claim":
            "none_post_benchmark_diagnostic_outside_published_duration",
    }
    if power is not None:
        electrical_power: dict[str, object] = {}
        for species, reported in (
            (electron_species, reported_electron_power),
            (ion_species, reported_ion_power),
        ):
            column = f"electric_work_{species}_J_m-2"
            work = (
                number(power[-1], column, "power transfer")
                - number(power[0], column, "power transfer")
            )
            measured = work / duration
            electrical_power[species] = {
                "electric_work_J_m2": work,
                "mean_electrical_power_W_m2": measured,
                "turner_table_iii_W_m2": reported,
                "relative_difference": measured / reported - 1.0,
            }
        report["window"]["power_counter_origin_step"] = power_origin
        report["provenance"]["power_transfer"] = {
            "path": str(power_path.resolve()),
            "sha256": sha256(power_path),
        }
        report["volume_electrical_power_context"] = {
            "definition":
                "discrete_species_kinetic_energy_change_due_to_electric_push",
            "species": electrical_power,
            "acceptance_gate": "none_diagnostic_context_only",
        }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scalars", type=Path, required=True)
    parser.add_argument("--collisions", type=Path, required=True)
    parser.add_argument("--boundary-losses", type=Path, required=True)
    parser.add_argument("--power-transfer", type=Path)
    parser.add_argument("--electron-species", default="electrons")
    parser.add_argument("--ion-species", default="ions")
    parser.add_argument(
        "--ionization-channel",
        default="electron_mcc.ionization",
    )
    parser.add_argument("--expected-steps", type=int, default=400)
    parser.add_argument(
        "--reported-ion-current", type=float, default=0.219
    )
    parser.add_argument(
        "--reported-electron-power", type=float, default=34.3
    )
    parser.add_argument(
        "--reported-ion-power", type=float, default=90.6
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        require(args.expected_steps > 0, "expected steps must be positive")
        require(
            math.isfinite(args.reported_ion_current)
            and args.reported_ion_current > 0.0,
            "reported ion current must be positive and finite",
        )
        require(
            math.isfinite(args.reported_electron_power)
            and args.reported_electron_power > 0.0
            and math.isfinite(args.reported_ion_power)
            and args.reported_ion_power > 0.0,
            "reported powers must be positive and finite",
        )
        report = analyze(
            args.scalars,
            args.collisions,
            args.boundary_losses,
            args.power_transfer,
            args.electron_species,
            args.ion_species,
            args.ionization_channel,
            args.expected_steps,
            args.reported_ion_current,
            args.reported_electron_power,
            args.reported_ion_power,
        )
        write_report(args.output, report)
    except (BalanceError, OSError) as error:
        print(f"Turner balance error: {error}", file=sys.stderr)
        return 2
    print(f"Turner balance report written to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
