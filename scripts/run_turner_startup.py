#!/usr/bin/env python3
"""Run and audit one checkpoint-split Turner Case 1 RF startup cycle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

from qualify_turner_runtime import (
    conservative_child_setup,
    set_global,
    sha256,
)


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_ONE_CYCLE_TURNER_STARTUP"
GENERATION_ACKNOWLEDGEMENT = (
    "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_RUN"
)
STEPS_PER_CYCLE = 400
SPLIT_STEP = 200
INITIAL_PARTICLES = 131_072
HARD_UPDATE_LIMIT = 60_000_000


class StartupError(RuntimeError):
    pass


def positive_integer(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return result


def insert_global(text: str, key: str, value: str) -> str:
    if re.search(rf"(?m)^\s*{re.escape(key)}\s*=", text):
        raise StartupError(f"generated deck already contains {key!r}")
    marker = text.find("\n[")
    if marker < 0:
        raise StartupError("generated deck has no configuration sections")
    return text[:marker] + f"\n{key} = {value}" + text[marker:]


def write_new(path: Path, text: str) -> None:
    if path.exists():
        raise StartupError(f"refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def rows(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            result = list(csv.DictReader(stream))
    except (OSError, UnicodeError, csv.Error) as error:
        raise StartupError(f"cannot read {path}: {error}") from error
    if not result:
        raise StartupError(f"diagnostic file is empty: {path}")
    return result


def number(row: dict[str, str], key: str, context: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, ValueError) as error:
        raise StartupError(f"{context} has invalid {key!r}") from error
    if not math.isfinite(value):
        raise StartupError(f"{context} has non-finite {key!r}")
    return value


def integer(row: dict[str, str], key: str, context: str) -> int:
    value = number(row, key, context)
    result = int(value)
    if value != result:
        raise StartupError(f"{context} has non-integer {key!r}")
    return result


def analyze(stage1: Path, stage2: Path, rf_frequency: float = 13.56e6) -> dict:
    scalars1 = rows(stage1 / "scalars.csv")
    scalars2 = rows(stage2 / "scalars.csv")
    collisions1 = rows(stage1 / "collisions.csv")
    collisions2 = rows(stage2 / "collisions.csv")
    if integer(scalars1[0], "step", "stage 1 scalars") != 0:
        raise StartupError("stage 1 did not start at step zero")
    if integer(scalars1[-1], "step", "stage 1 scalars") != SPLIT_STEP:
        raise StartupError("stage 1 did not end at the half-cycle")
    if integer(scalars2[0], "step", "stage 2 scalars") != SPLIT_STEP:
        raise StartupError("stage 2 did not restart at the half-cycle")
    if integer(scalars2[-1], "step", "stage 2 scalars") != STEPS_PER_CYCLE:
        raise StartupError("stage 2 did not complete one RF cycle")

    scalar_integer_keys = (
        "step", "live_particles",
        "live_particles_electrons", "live_particles_ions",
    )
    for key in scalar_integer_keys:
        if scalars1[-1].get(key) != scalars2[0].get(key):
            raise StartupError(f"checkpoint scalar discontinuity in {key}")
    scalar_float_keys = (
        "time", "kinetic_energy", "field_energy", "total_energy",
        "charge_l1", "phi_left", "phi_right",
    )
    maximum_restart_relative_difference = 0.0
    for key in scalar_float_keys:
        before = number(scalars1[-1], key, "stage 1 final scalars")
        after = number(scalars2[0], key, "stage 2 initial scalars")
        scale = max(abs(before), abs(after), 1e-300)
        relative = abs(before - after) / scale
        maximum_restart_relative_difference = max(
            maximum_restart_relative_difference, relative
        )
        if not math.isclose(before, after, rel_tol=1e-14, abs_tol=1e-18):
            raise StartupError(f"checkpoint scalar discontinuity in {key}")
    cumulative_keys = [
        key for key in collisions1[-1] if key.startswith("cumulative_")
    ]
    if not cumulative_keys:
        raise StartupError("collision diagnostics have no cumulative counters")
    for key in cumulative_keys:
        if collisions1[-1].get(key) != collisions2[0].get(key):
            raise StartupError(f"checkpoint collision discontinuity in {key}")

    combined_scalars = scalars1 + scalars2[1:]
    expected_steps = list(range(0, STEPS_PER_CYCLE + 1, 20))
    actual_steps = [integer(row, "step", "scalar history") for row in combined_scalars]
    if actual_steps != expected_steps:
        raise StartupError("scalar diagnostic cadence is incomplete")
    dt = number(combined_scalars[1], "time", "scalar history") / 20.0
    waveform_errors = []
    for row in combined_scalars:
        step = integer(row, "step", "scalar history")
        actual = number(row, "phi_right", "scalar history")
        expected = 450.0 * math.sin(
            2.0 * math.pi * rf_frequency * dt * step
        )
        waveform_errors.append(abs(actual - expected))
        for key in (
            "kinetic_energy", "field_energy", "total_energy", "charge_l1"
        ):
            number(row, key, "scalar history")
    maximum_waveform_error = max(waveform_errors)
    if maximum_waveform_error > 1e-9:
        raise StartupError("RF electrode waveform failed its exact-step check")

    initial_electrons = integer(
        combined_scalars[0], "live_particles_electrons", "initial scalars"
    )
    initial_ions = integer(
        combined_scalars[0], "live_particles_ions", "initial scalars"
    )
    final_electrons = integer(
        combined_scalars[-1], "live_particles_electrons", "final scalars"
    )
    final_ions = integer(
        combined_scalars[-1], "live_particles_ions", "final scalars"
    )
    final_collision = collisions2[-1]
    ionizations = integer(
        final_collision,
        "cumulative_collisions_electron_mcc.ionization",
        "final collisions",
    )
    electron_losses = initial_electrons + ionizations - final_electrons
    ion_losses = initial_ions + ionizations - final_ions
    if electron_losses < 0 or ion_losses < 0:
        raise StartupError("species balance implies negative electrode losses")
    if (
        final_electrons + final_ions
        != initial_electrons + initial_ions
        + 2 * ionizations - electron_losses - ion_losses
    ):
        raise StartupError("macro-particle balance is inconsistent")

    field_rows = rows(stage2 / f"fields_{STEPS_PER_CYCLE}.csv")
    electric = [abs(number(row, "E", "final field")) for row in field_rows]
    if len(electric) != 129:
        raise StartupError("final field does not contain the exact Case 1 grid")
    middle = electric[len(electric) // 4: 3 * len(electric) // 4]
    bulk_rms = math.sqrt(sum(value * value for value in middle) / len(middle))
    boundary_field = max(electric[0], electric[-1])
    maximum_field = max(electric)
    sheath_indicator = (
        boundary_field / bulk_rms if bulk_rms > 0 else math.inf
    )
    if not math.isfinite(sheath_indicator):
        raise StartupError("invalid sheath-field indicator")

    cumulative = {
        key.removeprefix("cumulative_"): integer(
            final_collision, key, "final collisions"
        )
        for key in cumulative_keys
    }
    return {
        "startup_checks_passed": True,
        "restart_continuity": {
            "split_step": SPLIT_STEP,
            "integer_fields_exact": True,
            "floating_fields_within_roundoff": True,
            "maximum_floating_relative_difference":
                maximum_restart_relative_difference,
            "collision_counters_exact": True,
        },
        "waveform": {
            "samples": len(combined_scalars),
            "maximum_absolute_error_v": maximum_waveform_error,
            "final_phi_right_v": number(
                combined_scalars[-1], "phi_right", "final scalars"
            ),
        },
        "particle_balance": {
            "initial_electrons": initial_electrons,
            "initial_ions": initial_ions,
            "ionization_pairs_created": ionizations,
            "inferred_electron_electrode_losses": electron_losses,
            "inferred_ion_electrode_losses": ion_losses,
            "final_electrons": final_electrons,
            "final_ions": final_ions,
        },
        "collisions": cumulative,
        "energy_and_charge": {
            "initial_total_energy_j": number(
                combined_scalars[0], "total_energy", "initial scalars"
            ),
            "final_total_energy_j": number(
                combined_scalars[-1], "total_energy", "final scalars"
            ),
            "maximum_charge_l1_c": max(
                number(row, "charge_l1", "scalar history")
                for row in combined_scalars
            ),
            "all_samples_finite": True,
        },
        "early_field_structure": {
            "final_boundary_field_max_v_m": boundary_field,
            "final_bulk_field_rms_v_m": bulk_rms,
            "final_domain_field_max_v_m": maximum_field,
            "boundary_to_bulk_rms_ratio": sheath_indicator,
            "interpretation": "startup indicator only; not a converged sheath",
        },
    }


def stage_deck(
    production_text: str,
    steps: int,
    output_dir: Path,
    checkpoint: Path,
    restart: Path | None = None,
) -> str:
    replacements = {
        "steps": str(steps),
        "output_interval": "20",
        "output_dir": str(output_dir),
        "checkpoint_output": "true",
        "checkpoint_interval": str(SPLIT_STEP),
        "spatial_average": "false",
        "spatial_average_interval": "1",
        "spatial_average_start_step": "0",
        "spatial_average_end_step": "0",
        "spatial_average_rf_frequency": "0",
        "spatial_average_rf_cycles": "0",
        "runtime_backend": "serial",
        "runtime_threads": "1",
    }
    result = production_text
    for key, value in replacements.items():
        result = set_global(result, key, value)
    result = insert_global(result, "checkpoint_path", str(checkpoint))
    if restart is not None:
        result = insert_global(result, "restart_path", str(restart))
    return result


def execute(args: argparse.Namespace) -> dict:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise StartupError(
            "startup run requires --acknowledge-cost " + ACKNOWLEDGEMENT
        )
    updates = INITIAL_PARTICLES * STEPS_PER_CYCLE
    if args.max_initial_updates > HARD_UPDATE_LIMIT:
        raise StartupError(
            "--max-initial-updates exceeds the built-in limit: "
            f"{args.max_initial_updates} > {HARD_UPDATE_LIMIT}"
        )
    if updates > args.max_initial_updates:
        raise StartupError(
            f"startup run exceeds --max-initial-updates: {updates} > "
            f"{args.max_initial_updates}"
        )
    executable = args.executable.resolve()
    work = args.work_dir.resolve()
    report_path = args.report.resolve()
    if not executable.is_file():
        raise StartupError(f"missing AuroraPIC CLI: {executable}")
    if work.exists():
        raise StartupError(f"refusing to overwrite work directory: {work}")
    if report_path.exists():
        raise StartupError(f"refusing to overwrite report: {report_path}")
    work.mkdir(parents=True)
    production = work / "production.cfg"
    preparer = Path(__file__).resolve().with_name("prepare_turner_case.py")
    generated = subprocess.run(
        [
            sys.executable, str(preparer), str(args.case_manifest.resolve()),
            str(args.normalized_dir.resolve()), "--output", str(production),
            "--acknowledge-cost", GENERATION_ACKNOWLEDGEMENT,
        ],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
    )
    if generated.returncode != 0:
        raise StartupError(
            "cannot generate exact campaign deck: " + generated.stderr.strip()
        )
    production_text = production.read_text(encoding="utf-8")
    stage1_output = work / "stage1-output"
    stage2_output = work / "stage2-output"
    checkpoint1 = work / "checkpoint-half-cycle.apc"
    checkpoint2 = work / "checkpoint-one-cycle.apc"
    stage1_config = work / "stage1.cfg"
    stage2_config = work / "stage2.cfg"
    write_new(
        stage1_config,
        stage_deck(
            production_text, SPLIT_STEP, stage1_output, checkpoint1
        ),
    )
    write_new(
        stage2_config,
        stage_deck(
            production_text, STEPS_PER_CYCLE, stage2_output, checkpoint2,
            checkpoint1,
        ),
    )
    environment = dict(os.environ)
    environment.update({
        "OMP_NUM_THREADS": "1",
        "OMP_DYNAMIC": "FALSE",
        "AURORA_OPENMP_THREADS": "1",
    })
    timings = []
    for label, config in (("half_cycle", stage1_config), ("full_cycle", stage2_config)):
        start = time.perf_counter()
        try:
            result = subprocess.run(
                [str(executable), str(config)],
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=args.timeout_seconds,
                preexec_fn=conservative_child_setup if os.name == "posix" else None,
            )
        except subprocess.TimeoutExpired as error:
            raise StartupError(
                f"{label} stage exceeded {args.timeout_seconds} seconds"
            ) from error
        elapsed = time.perf_counter() - start
        if result.returncode != 0:
            raise StartupError(
                f"{label} stage failed: {result.stdout.strip()} "
                f"{result.stderr.strip()}"
            )
        timings.append({"stage": label, "wall_seconds": elapsed})
    if not checkpoint1.is_file() or not checkpoint2.is_file():
        raise StartupError("startup stages did not write both checkpoints")
    diagnostics = analyze(stage1_output, stage2_output)
    return {
        "turner_startup_report_version": 1,
        "case_id": "turner-helium-ccp-2013-case-1",
        "scope": "one_checkpoint_split_rf_startup_cycle",
        "physics_claim": "none",
        "steady_state_claim": False,
        "production_launch_authorized": False,
        "work": {
            "steps": STEPS_PER_CYCLE,
            "rf_cycles": 1,
            "initial_particles": INITIAL_PARTICLES,
            "initial_particle_updates": updates,
            "stage_timings": timings,
            "total_wall_seconds": sum(item["wall_seconds"] for item in timings),
        },
        "diagnostics": diagnostics,
        "provenance": {
            "case_manifest_sha256": sha256(args.case_manifest.resolve()),
            "normalization_audit_sha256": sha256(
                args.normalized_dir.resolve() / "audit.json"
            ),
            "executable_sha256": sha256(executable),
            "production_config_sha256": sha256(production),
            "stage1_config_sha256": sha256(stage1_config),
            "stage2_config_sha256": sha256(stage2_config),
            "half_cycle_checkpoint_sha256": sha256(checkpoint1),
            "one_cycle_checkpoint_sha256": sha256(checkpoint2),
            "platform": platform.platform(),
        },
        "warnings": [
            "One RF cycle is far shorter than the 1280-cycle benchmark.",
            "Early field structure is not a converged sheath.",
            "The published chi-squared acceptance range is inapplicable.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("normalized_dir", type=Path)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument(
        "--max-initial-updates", type=positive_integer, default=60_000_000
    )
    parser.add_argument(
        "--timeout-seconds", type=positive_integer, default=120
    )
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = execute(args)
        write_new(
            args.report.resolve(),
            json.dumps(report, indent=2, sort_keys=True) + "\n",
        )
    except (StartupError, OSError, UnicodeError, json.JSONDecodeError) as error:
        print(f"Turner startup error: {error}", file=sys.stderr)
        return 2
    print(
        "Completed one bounded checkpoint-split RF startup cycle without "
        f"a physics claim: {args.report}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
