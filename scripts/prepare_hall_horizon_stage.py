#!/usr/bin/env python3
"""Prepare one guarded, checkpoint-chained Hall time-horizon stage."""

from __future__ import annotations

import argparse
import configparser
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile


ACKNOWLEDGEMENT = (
    "I_UNDERSTAND_THIS_EXTENDS_A_HALL_RUN_FROM_A_PINNED_CHECKPOINT"
)


class HorizonError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise HorizonError(f"cannot hash {path}: {error}") from error


def load_ini(path: Path, add_global: bool) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(
        interpolation=None,
        strict=True,
        empty_lines_in_values=False,
    )
    try:
        text = path.read_text(encoding="utf-8")
        parser.read_string(("[global]\n" if add_global else "") + text)
    except (OSError, UnicodeError, configparser.Error) as error:
        raise HorizonError(f"cannot read {path}: {error}") from error
    return parser


def positive(value: float, context: str) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise HorizonError(f"{context} must be positive and finite")
    return value


def read_scalar_endpoints(
    path: Path,
) -> tuple[dict[str, str], dict[str, str]]:
    required = {
        "step", "time", "live_particles",
        "live_particles_electrons", "live_particles_ions",
    }
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            fields = set(reader.fieldnames or ())
            if not required <= fields:
                raise HorizonError(
                    f"{path} is missing columns {sorted(required - fields)}"
                )
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise HorizonError(f"cannot read {path}: {error}") from error
    if len(rows) < 2:
        raise HorizonError("prior scalar history requires at least two rows")
    return rows[0], rows[-1]


def integer(row: dict[str, str], key: str, context: str) -> int:
    try:
        numeric = float(row[key])
        result = int(numeric)
    except (KeyError, ValueError) as error:
        raise HorizonError(f"{context} {key} is not an integer") from error
    if numeric != result or result < 0:
        raise HorizonError(f"{context} {key} is not a non-negative integer")
    return result


def number(row: dict[str, str], key: str, context: str) -> float:
    try:
        return positive(float(row[key]), f"{context} {key}")
    except (KeyError, ValueError) as error:
        raise HorizonError(f"{context} {key} is not numeric") from error


def replace_global(text: str, key: str, value: object) -> str:
    pattern = re.compile(rf"(?m)^({re.escape(key)}\s*=).*$")
    if len(pattern.findall(text)) != 1:
        raise HorizonError(f"prior deck must contain one global {key}")
    return pattern.sub(rf"\g<1> {value}", text, count=1)


def insert_global(text: str, key: str, value: object) -> str:
    if re.search(rf"(?m)^{re.escape(key)}\s*=", text):
        return replace_global(text, key, value)
    marker = re.search(r"(?m)^initial_state_path\s*=.*$", text)
    if marker:
        raise HorizonError(
            "checkpoint extension cannot retain initial_state_path"
        )
    boundary = re.search(r"(?m)^\[", text)
    if boundary is None:
        raise HorizonError("prior deck has no sections")
    return text[:boundary.start()] + f"{key} = {value}\n\n" + text[boundary.start():]


def checkpoint_identity(path: Path) -> tuple[int, float]:
    try:
        with path.open(encoding="utf-8") as stream:
            magic = stream.readline().strip()
            dimension = stream.readline().split()
            units = stream.readline().split()
            step = stream.readline().split()
            time = stream.readline().split()
    except (OSError, UnicodeError) as error:
        raise HorizonError(f"cannot inspect checkpoint {path}: {error}") from error
    if (
        re.fullmatch(r"AuroraPIC-checkpoint-v[1-9][0-9]*", magic) is None
        or dimension != ["dimension", "2"]
        or not units
        or units[0] != "units"
        or len(step) != 2
        or step[0] != "step"
        or len(time) != 2
        or time[0] != "time"
    ):
        raise HorizonError("prior checkpoint header is invalid")
    try:
        return int(step[1]), float(time[1])
    except ValueError as error:
        raise HorizonError("prior checkpoint step/time is invalid") from error


def prepare(args: argparse.Namespace) -> Path:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise HorizonError(
            "horizon generation requires --acknowledge-cost "
            + ACKNOWLEDGEMENT
        )
    destination = args.output_dir.resolve()
    if destination.exists():
        raise HorizonError(f"refusing to overwrite output: {destination}")
    case_path = args.case_manifest.resolve()
    deck_path = args.prior_deck.resolve()
    prior_output = args.prior_output.resolve()
    case = load_ini(case_path, True)
    deck = load_ini(deck_path, True)
    for section in (
        "global", "campaign.workstation", "horizon.workstation",
    ):
        if section not in case:
            raise HorizonError(f"case manifest is missing [{section}]")
    contract = case["horizon.workstation"]
    if (
        contract.get("horizon_contract_version") != "1"
        or contract.get("base_tier") != "workstation"
        or contract.get("physics_claim") != "none"
    ):
        raise HorizonError("unsupported Hall horizon contract")
    global_deck = deck["global"]
    if (
        global_deck.get("dimension") != "2"
        or global_deck.get("mode") != "transient"
        or global_deck.get("runtime_backend") != "serial"
        or global_deck.getint("runtime_threads") != contract.getint("max_threads")
        or global_deck.getint("runtime_threads") != 1
    ):
        raise HorizonError(
            "prior deck must be a serial, single-thread transient 2D run"
        )
    configured_output = Path(global_deck["output_dir"]).resolve()
    if configured_output != prior_output:
        raise HorizonError("prior output does not match the prior deck")
    prior_steps = global_deck.getint("steps")
    dt = positive(global_deck.getfloat("dt"), "prior timestep")
    initial, final = read_scalar_endpoints(prior_output / "scalars.csv")
    if integer(final, "step", "final scalar") != prior_steps:
        raise HorizonError("prior run did not reach its configured final step")
    final_time = number(final, "time", "final scalar")
    if not math.isclose(final_time, prior_steps * dt, rel_tol=1e-10):
        raise HorizonError("prior final time does not match step times dt")
    checkpoint = prior_output / f"checkpoint_{prior_steps}.apc"
    checkpoint_step, checkpoint_time = checkpoint_identity(checkpoint)
    if (
        checkpoint_step != prior_steps
        or not math.isclose(checkpoint_time, final_time, rel_tol=1e-12)
    ):
        raise HorizonError("checkpoint does not match the prior final state")

    multiplier = contract.getint("step_multiplier")
    if multiplier < 2 or multiplier > 4:
        raise HorizonError("step_multiplier must be between 2 and 4")
    target_steps = prior_steps * multiplier
    production_steps = case["reference"].getint("production_steps")
    if target_steps > production_steps:
        target_steps = production_steps
    if target_steps <= prior_steps:
        raise HorizonError("prior run already reached the production horizon")
    averaging_fraction = contract.getfloat(
        "diagnostic_averaging_fraction"
    )
    samples = contract.getint("diagnostic_samples")
    if not 0.0 < averaging_fraction < 1.0 or samples < 2:
        raise HorizonError("invalid horizon diagnostic contract")
    diagnostic_start = int(round(
        target_steps * (1.0 - averaging_fraction)
    ))
    if diagnostic_start <= prior_steps:
        raise HorizonError(
            "new averaging window must begin after the prior checkpoint"
        )
    window = target_steps - diagnostic_start
    if window % (samples - 1):
        raise HorizonError(
            "target horizon cannot produce the requested diagnostic cadence"
        )
    interval = window // (samples - 1)

    growth_safety = positive(
        contract.getfloat("population_growth_safety_factor"),
        "population growth safety factor",
    )
    headroom = positive(
        contract.getfloat("capacity_headroom_factor"),
        "capacity headroom factor",
    )
    species = ("electrons", "ions")
    populations: dict[str, dict[str, object]] = {}
    projected_capacity = 0
    added_steps = target_steps - prior_steps
    for name in species:
        first = integer(
            initial, f"live_particles_{name}", "initial scalar"
        )
        last = integer(final, f"live_particles_{name}", "final scalar")
        net_growth_per_step = max(0.0, (last - first) / prior_steps)
        projected = math.ceil(
            last + growth_safety * net_growth_per_step * added_steps
        )
        required_capacity = math.ceil(max(projected, last) * headroom)
        projected_capacity = max(projected_capacity, required_capacity)
        populations[name] = {
            "initial": first,
            "final": last,
            "observed_net_growth_per_step": net_growth_per_step,
            "safety_projected_final": projected,
            "required_capacity_with_headroom": required_capacity,
        }
    prior_capacity = global_deck.getint("max_particles_per_species")
    capacity = max(prior_capacity, projected_capacity)
    final_total = integer(final, "live_particles", "final scalar")
    added_updates = final_total * added_steps
    maximum_updates = contract.getint("maximum_added_particle_updates")
    if added_updates > maximum_updates:
        raise HorizonError(
            "stage exceeds maximum_added_particle_updates"
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{destination.name}.", suffix=".tmp",
        dir=destination.parent,
    ))
    try:
        stage_output = destination / "output"
        text = deck_path.read_text(encoding="utf-8")
        for key, value in (
            ("steps", target_steps),
            ("output_interval", interval),
            ("output_dir", stage_output),
            ("resolved_diagnostic_interval", interval),
            ("resolved_diagnostic_start_step", diagnostic_start),
            ("checkpoint_interval", target_steps),
            ("max_particles_per_species", capacity),
            ("runtime_backend", "serial"),
            ("runtime_threads", 1),
        ):
            text = replace_global(text, key, value)
        text = insert_global(text, "restart_path", checkpoint)
        magnetic = Path(global_deck["magnetic_field_profile_file"])
        if not magnetic.is_absolute():
            magnetic = (deck_path.parent / magnetic).resolve()
        if not magnetic.is_file():
            raise HorizonError("prior magnetic-field asset is missing")
        text = replace_global(
            text, "magnetic_field_profile_file", magnetic
        )
        generated_deck = temporary / "horizon.cfg"
        generated_deck.write_text(text, encoding="utf-8")
        report = {
            "hall_horizon_stage_version": 1,
            "case_id": case["global"]["case_id"],
            "physics_claim": "none",
            "launched": False,
            "case_manifest": str(case_path),
            "case_manifest_sha256": sha256(case_path),
            "prior_deck": str(deck_path),
            "prior_deck_sha256": sha256(deck_path),
            "prior_output": str(prior_output),
            "prior_scalars_sha256": sha256(prior_output / "scalars.csv"),
            "restart_checkpoint": str(checkpoint),
            "restart_checkpoint_sha256": sha256(checkpoint),
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "prior_step": prior_steps,
            "prior_time_s": final_time,
            "target_step": target_steps,
            "target_time_s": target_steps * dt,
            "added_steps": added_steps,
            "diagnostic_start_step": diagnostic_start,
            "diagnostic_interval": interval,
            "diagnostic_samples": samples,
            "population_projection": populations,
            "max_particles_per_species": capacity,
            "estimated_added_particle_updates_lower_bound": added_updates,
            "maximum_added_particle_updates": maximum_updates,
            "runtime_backend": "serial",
            "runtime_threads": 1,
            "runtime_config": str(destination / generated_deck.name),
            "runtime_config_sha256": sha256(generated_deck),
            "result_dir": str(stage_output),
            "next_gate": (
                "Run analyze_hall_pilot.py-compatible integrity checks, "
                "compare_hall.py against the pinned screening reference, "
                "and inspect profile trend before preparing another stage."
            ),
            "warnings": [
                "This command prepared a deck and did not launch it.",
                "The CLI large-run acknowledgement remains required.",
                "Measured net growth is a projection, not a capacity proof.",
                "Published-profile agreement and time-window matching remain "
                "independent requirements.",
            ],
        }
        (temporary / "horizon.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, destination)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return destination


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare, but never launch, one checkpoint-chained Hall horizon stage"
        )
    )
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("--prior-deck", type=Path, required=True)
    parser.add_argument("--prior-output", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    try:
        output = prepare(parse_args())
    except (HorizonError, KeyError, ValueError, configparser.Error) as error:
        print(f"Hall horizon preparation error: {error}", file=sys.stderr)
        return 2
    print(f"Prepared Hall horizon stage without launching it: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
