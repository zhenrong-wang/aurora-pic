#!/usr/bin/env python3
"""Validate bounded Hall campaign diagnostics without making a physics claim."""

from __future__ import annotations

import argparse
import configparser
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile


class PilotError(RuntimeError):
    pass


def load_manifest(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    try:
        parser.read_string("[global]\n" + path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, configparser.Error) as error:
        raise PilotError(f"cannot read manifest {path}: {error}") from error
    return parser


def read_rows(path: Path, required: set[str]) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            fields = set(reader.fieldnames or ())
            if not required <= fields:
                missing = ", ".join(sorted(required - fields))
                raise PilotError(f"{path.name} is missing columns: {missing}")
            rows = list(reader)
    except OSError as error:
        raise PilotError(f"cannot read {path}: {error}") from error
    if not rows:
        raise PilotError(f"{path.name} has no data rows")
    return rows


def number(row: dict[str, str], key: str, source: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, ValueError) as error:
        raise PilotError(f"{source} contains invalid {key}") from error
    if not math.isfinite(value):
        raise PilotError(f"{source} contains non-finite {key}")
    return value


def integer(row: dict[str, str], key: str, source: str) -> int:
    value = number(row, key, source)
    result = int(value)
    if value != result:
        raise PilotError(f"{source} contains non-integral {key}")
    return result


def atomic_json(path: Path, report: dict[str, object]) -> None:
    if path.exists():
        raise PilotError(f"refusing to overwrite existing report: {path}")
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


def analyze(args: argparse.Namespace) -> dict[str, object]:
    manifest_path = args.case_manifest.resolve()
    manifest = load_manifest(manifest_path)
    tier_name = f"campaign.{args.tier}"
    if tier_name not in manifest:
        raise PilotError(f"manifest is missing [{tier_name}]")
    tier = manifest[tier_name]
    if tier["physics_claim"] != "none":
        raise PilotError(
            "pilot analyzer only accepts tiers with physics_claim = none"
        )
    steps = tier.getint("steps")
    interval = tier.getint("diagnostic_interval")
    start = tier.getint("diagnostic_start_step")
    cells_x = tier.getint("cells_x")
    nodes_x = cells_x + 1
    expected_scalar_samples = steps // interval + 1
    expected_resolved_samples = (steps - start) // interval + 1
    output = args.output_dir.resolve()

    scalars = read_rows(
        output / "scalars.csv",
        {"step", "time", "total_energy", "charge_l1", "live_particles"},
    )
    scalar_steps = [
        integer(row, "step", "scalars.csv") for row in scalars
    ]
    if (len(scalars) != expected_scalar_samples
            or scalar_steps[0] != 0 or scalar_steps[-1] != steps
            or scalar_steps != sorted(set(scalar_steps))):
        raise PilotError("scalar sampling cadence does not match the tier")
    maximum_live = 0
    for row in scalars:
        number(row, "time", "scalars.csv")
        number(row, "total_energy", "scalars.csv")
        number(row, "charge_l1", "scalars.csv")
        maximum_live = max(
            maximum_live, integer(row, "live_particles", "scalars.csv")
        )
    capacity = 2 * tier.getint("max_particles_per_species")
    if maximum_live > capacity:
        raise PilotError("reported live population exceeds configured capacity")

    sources = read_rows(
        output / "sources.csv",
        {
            "step", "resolved_represented_pair_rate",
            "fractional_macro_pair_remainder", "macro_pairs_created",
        },
    )
    expected_rate = manifest["pair_source"].getfloat(
        "derived_represented_pair_rate_s"
    )
    for row in sources:
        rate = number(row, "resolved_represented_pair_rate", "sources.csv")
        integer(row, "macro_pairs_created", "sources.csv")
        remainder = number(
            row, "fractional_macro_pair_remainder", "sources.csv"
        )
        if not math.isclose(rate, expected_rate, rel_tol=1e-14):
            raise PilotError("source represented rate drifted from the manifest")
        if remainder < 0.0 or remainder >= 1.0 + 1e-12:
            raise PilotError("source fractional remainder is outside [0,1)")

    current = read_rows(
        output / "current_source.csv",
        {
            "step", "charge_balance_residual",
            "cumulative_processed_monitored_charge",
            "cumulative_emitted_charge", "macro_particles_created",
            "control_macro_remainder",
        },
    )
    initial_particles_per_species = (
        tier.getint("cells_x")
        * tier.getint("cells_y")
        * tier.getint("particles_per_cell_per_species")
    )
    macro_weight = (
        manifest["reference"].getfloat("initial_density_m3")
        * manifest["reference"].getfloat("domain_x_m")
        * manifest["reference"].getfloat("domain_y_m")
        * manifest["pair_source"].getfloat("out_of_plane_depth_m")
        / initial_particles_per_species
    )
    macro_charge = 1.602176634e-19 * macro_weight
    maximum_charge_residual = 0.0
    maximum_controller_debt = 0.0
    controller_saturation_samples = 0
    for row in current:
        processed_charge = number(
            row, "cumulative_processed_monitored_charge",
            "current_source.csv"
        )
        emitted_charge = number(
            row, "cumulative_emitted_charge", "current_source.csv"
        )
        emitted_macroparticles = integer(
            row, "macro_particles_created", "current_source.csv"
        )
        remainder = number(
            row, "control_macro_remainder", "current_source.csv"
        )
        residual = number(
            row, "charge_balance_residual", "current_source.csv"
        )
        # The pinned emitted species is electrons, so its represented
        # macrocharge is negative. A negative remainder is valid actuator
        # debt: electron emission cannot retract particles already emitted.
        emitted_macrocharge = -macro_charge
        tolerance = macro_charge * 1e-9
        if not math.isclose(
            emitted_charge,
            emitted_macroparticles * emitted_macrocharge,
            rel_tol=1e-12,
            abs_tol=tolerance,
        ):
            raise PilotError(
                "cathode emitted-charge accounting is inconsistent"
            )
        if not math.isclose(
            residual, processed_charge - emitted_charge,
            rel_tol=1e-12, abs_tol=tolerance,
        ):
            raise PilotError(
                "cathode charge-balance residual is inconsistent"
            )
        if not math.isclose(
            residual, remainder * emitted_macrocharge,
            rel_tol=1e-12, abs_tol=tolerance,
        ):
            raise PilotError(
                "cathode control remainder is inconsistent with its residual"
            )
        if remainder >= 1.0 + 1e-12:
            raise PilotError(
                "cathode positive control remainder exceeds one macro charge"
            )
        if remainder < -1e-12:
            controller_saturation_samples += 1
            maximum_controller_debt = max(
                maximum_controller_debt, -remainder
            )
        maximum_charge_residual = max(
            maximum_charge_residual, abs(residual)
        )

    potential = read_rows(
        output / "potential_reference.csv",
        {"step", "target", "corrected_line_mean"},
    )
    maximum_reference_error = 0.0
    for row in potential:
        error = abs(
            number(row, "corrected_line_mean", "potential_reference.csv")
            - number(row, "target", "potential_reference.csv")
        )
        maximum_reference_error = max(maximum_reference_error, error)
    if maximum_reference_error > 1e-9:
        raise PilotError("potential-reference correction exceeded tolerance")

    field_average = read_rows(
        output / "resolved_field_time_average.csv",
        {"samples", "profile_axis", "coordinate", "electric_x"},
    )
    if len(field_average) != nodes_x or any(
        integer(row, "samples", "resolved_field_time_average.csv")
            != expected_resolved_samples
        for row in field_average
    ):
        raise PilotError("field time-average shape or sample count is invalid")
    for row in field_average:
        number(row, "coordinate", "resolved_field_time_average.csv")
        number(row, "electric_x", "resolved_field_time_average.csv")
    if ({row["profile_axis"] for row in field_average} != {"x"}
            or len({row["coordinate"] for row in field_average}) != nodes_x):
        raise PilotError("field time-average coordinate coverage is invalid")
    species_average = read_rows(
        output / "resolved_species_time_average.csv",
        {"samples", "profile_axis", "coordinate", "species", "number_density"},
    )
    if len(species_average) != 2 * nodes_x or any(
        integer(row, "samples", "resolved_species_time_average.csv")
            != expected_resolved_samples
        for row in species_average
    ):
        raise PilotError("species time-average shape or sample count is invalid")
    for row in species_average:
        number(row, "coordinate", "resolved_species_time_average.csv")
        number(row, "number_density", "resolved_species_time_average.csv")
    species_names = {row["species"] for row in species_average}
    if (species_names != {"electrons", "ions"}
            or any(
                sum(row["species"] == species for row in species_average)
                    != nodes_x
                for species in species_names
            )):
        raise PilotError("species time-average coverage is invalid")
    modes = read_rows(
        output / "resolved_modes.csv",
        {"step", "mode", "real", "imaginary", "amplitude"},
    )
    observed_modes = {
        integer(row, "mode", "resolved_modes.csv") for row in modes
    }
    if observed_modes != set(range(tier.getint("max_mode") + 1)):
        raise PilotError("resolved mode coverage does not match the tier")
    for row in modes:
        number(row, "real", "resolved_modes.csv")
        number(row, "imaginary", "resolved_modes.csv")
        number(row, "amplitude", "resolved_modes.csv")
    observed_mode_steps = {
        integer(row, "step", "resolved_modes.csv") for row in modes
    }
    if observed_mode_steps != set(range(start, steps + 1, interval)):
        raise PilotError("resolved mode sampling cadence is invalid")

    return {
        "schema_version": 1,
        "case_id": manifest["global"]["case_id"],
        "campaign_tier": args.tier,
        "physics_claim": "none",
        "passed": True,
        "case_manifest": str(manifest_path),
        "case_manifest_sha256": hashlib.sha256(
            manifest_path.read_bytes()
        ).hexdigest(),
        "output_dir": str(output),
        "metrics": {
            "steps": steps,
            "final_time_s": number(scalars[-1], "time", "scalars.csv"),
            "scalar_samples": len(scalars),
            "resolved_samples": expected_resolved_samples,
            "initial_macroparticles":
                2 * initial_particles_per_species,
            "final_live_particles": integer(
                scalars[-1], "live_particles", "scalars.csv"
            ),
            "maximum_live_particles": maximum_live,
            "macro_pairs_created": integer(
                sources[-1], "macro_pairs_created", "sources.csv"
            ),
            "cathode_macro_particles_created": integer(
                current[-1], "macro_particles_created", "current_source.csv"
            ),
            "maximum_charge_balance_residual_c": maximum_charge_residual,
            "controller_saturation_samples":
                controller_saturation_samples,
            "maximum_controller_debt_macroparticles":
                maximum_controller_debt,
            "maximum_potential_reference_error_v": maximum_reference_error,
            "resolved_modes": len(observed_modes),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check bounded Hall pilot outputs without a physics claim"
    )
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument(
        "--tier", choices=("micro", "workstation"), default="micro"
    )
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = analyze(args)
        atomic_json(args.report, report)
    except (PilotError, ValueError, configparser.Error) as error:
        print(f"Hall pilot analysis error: {error}", file=sys.stderr)
        return 2
    print(f"Hall {args.tier} pilot diagnostics passed: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
