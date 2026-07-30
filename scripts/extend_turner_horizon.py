#!/usr/bin/env python3
"""Extend a Turner startup checkpoint through bounded whole-RF-cycle stages."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

from qualify_turner_runtime import conservative_child_setup, sha256
from run_turner_startup import (
    GENERATION_ACKNOWLEDGEMENT,
    INITIAL_PARTICLES,
    StartupError,
    integer,
    number,
    rows,
    stage_deck,
    write_new,
)


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_BOUNDED_TURNER_HORIZON"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"
STEPS_PER_CYCLE = 400
MAX_ADDITIONAL_CYCLES = 3
HARD_ADDITIONAL_UPDATE_LIMIT = 160_000_000


class HorizonError(RuntimeError):
    pass


def positive_integer(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return result


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise HorizonError(f"cannot read {path}: {error}") from error


def continuity(
    prior_scalar: dict[str, str],
    current_scalar: dict[str, str],
    prior_collision: dict[str, str],
    current_collision: dict[str, str],
) -> float:
    for key in (
        "step", "live_particles", "live_particles_electrons",
        "live_particles_ions",
    ):
        if (
            integer(prior_scalar, key, "prior final scalars")
            != integer(current_scalar, key, "current initial scalars")
        ):
            raise HorizonError(f"checkpoint scalar discontinuity in {key}")
    maximum_relative_difference = 0.0
    for key in (
        "time", "kinetic_energy", "field_energy", "total_energy",
        "charge_l1", "phi_left", "phi_right",
    ):
        before = number(prior_scalar, key, "prior final scalars")
        after = number(current_scalar, key, "current initial scalars")
        scale = max(abs(before), abs(after), 1e-300)
        relative = abs(before - after) / scale
        maximum_relative_difference = max(
            maximum_relative_difference, relative
        )
        if not math.isclose(before, after, rel_tol=1e-14, abs_tol=1e-18):
            raise HorizonError(f"checkpoint scalar discontinuity in {key}")
    cumulative = [
        key for key in prior_collision if key.startswith("cumulative_")
    ]
    if not cumulative:
        raise HorizonError("prior collision diagnostics lack cumulative counters")
    for key in cumulative:
        if (
            integer(prior_collision, key, "prior final collisions")
            != integer(current_collision, key, "current initial collisions")
        ):
            raise HorizonError(f"checkpoint collision discontinuity in {key}")
    return maximum_relative_difference


def field_metrics(path: Path) -> dict:
    field = rows(path)
    electric = [abs(number(row, "E", "phase-zero field")) for row in field]
    if len(electric) != 129:
        raise HorizonError("phase-zero field does not contain 129 nodes")
    middle = electric[len(electric) // 4: 3 * len(electric) // 4]
    bulk_rms = math.sqrt(sum(value * value for value in middle) / len(middle))
    boundary = max(electric[0], electric[-1])
    if bulk_rms <= 0 or not math.isfinite(bulk_rms):
        raise HorizonError("invalid bulk field RMS")
    return {
        "boundary_field_max_v_m": boundary,
        "bulk_field_rms_v_m": bulk_rms,
        "domain_field_max_v_m": max(electric),
        "boundary_to_bulk_rms_ratio": boundary / bulk_rms,
    }


def analyze_cycle(
    cycle: int,
    prior_scalar: dict[str, str],
    prior_collision: dict[str, str],
    output: Path,
) -> tuple[dict, dict[str, str], dict[str, str]]:
    scalar_rows = rows(output / "scalars.csv")
    collision_rows = rows(output / "collisions.csv")
    start_step = (cycle - 1) * STEPS_PER_CYCLE
    end_step = cycle * STEPS_PER_CYCLE
    if (
        integer(scalar_rows[0], "step", "cycle initial scalars") != start_step
        or integer(scalar_rows[-1], "step", "cycle final scalars") != end_step
    ):
        raise HorizonError(f"cycle {cycle} scalar endpoints are incorrect")
    expected = list(range(start_step, end_step + 1, 20))
    actual = [integer(row, "step", "cycle scalars") for row in scalar_rows]
    if actual != expected:
        raise HorizonError(f"cycle {cycle} scalar cadence is incomplete")
    restart_relative = continuity(
        prior_scalar, scalar_rows[0], prior_collision, collision_rows[0]
    )
    dt = number(scalar_rows[1], "time", "cycle scalars") - number(
        scalar_rows[0], "time", "cycle scalars"
    )
    dt /= 20.0
    waveform_error = 0.0
    for row in scalar_rows:
        step = integer(row, "step", "cycle scalars")
        actual_phi = number(row, "phi_right", "cycle scalars")
        expected_phi = 450.0 * math.sin(
            2.0 * math.pi * 13.56e6 * dt * step
        )
        waveform_error = max(waveform_error, abs(actual_phi - expected_phi))
        for key in (
            "kinetic_energy", "field_energy", "total_energy", "charge_l1"
        ):
            number(row, key, "cycle scalars")
    if waveform_error > 1e-9:
        raise HorizonError(f"cycle {cycle} electrode waveform is inconsistent")

    cumulative_keys = [
        key for key in collision_rows[-1] if key.startswith("cumulative_")
    ]
    collision_delta = {
        key.removeprefix("cumulative_"): (
            integer(collision_rows[-1], key, "cycle final collisions")
            - integer(collision_rows[0], key, "cycle initial collisions")
        )
        for key in cumulative_keys
    }
    if any(value < 0 for value in collision_delta.values()):
        raise HorizonError(f"cycle {cycle} collision counters decreased")
    ionizations = collision_delta.get(
        "collisions_electron_mcc.ionization"
    )
    if ionizations is None:
        raise HorizonError("ionization collision counter is absent")
    initial_electrons = integer(
        scalar_rows[0], "live_particles_electrons", "cycle initial scalars"
    )
    initial_ions = integer(
        scalar_rows[0], "live_particles_ions", "cycle initial scalars"
    )
    final_electrons = integer(
        scalar_rows[-1], "live_particles_electrons", "cycle final scalars"
    )
    final_ions = integer(
        scalar_rows[-1], "live_particles_ions", "cycle final scalars"
    )
    electron_losses = initial_electrons + ionizations - final_electrons
    ion_losses = initial_ions + ionizations - final_ions
    if electron_losses < 0 or ion_losses < 0:
        raise HorizonError(f"cycle {cycle} implies negative electrode losses")
    metrics = {
        "cycle": cycle,
        "start_step": start_step,
        "end_step": end_step,
        "restart_maximum_floating_relative_difference": restart_relative,
        "waveform_maximum_absolute_error_v": waveform_error,
        "population": {
            "initial_electrons": initial_electrons,
            "initial_ions": initial_ions,
            "ionization_pairs_created": ionizations,
            "inferred_electron_electrode_losses": electron_losses,
            "inferred_ion_electrode_losses": ion_losses,
            "final_electrons": final_electrons,
            "final_ions": final_ions,
            "electron_relative_change": (
                final_electrons - initial_electrons
            ) / initial_electrons,
            "ion_relative_change": (
                final_ions - initial_ions
            ) / initial_ions,
        },
        "collisions": collision_delta,
        "energy_and_charge": {
            "initial_total_energy_j": number(
                scalar_rows[0], "total_energy", "cycle initial scalars"
            ),
            "final_total_energy_j": number(
                scalar_rows[-1], "total_energy", "cycle final scalars"
            ),
            "maximum_charge_l1_c": max(
                number(row, "charge_l1", "cycle scalars")
                for row in scalar_rows
            ),
            "all_samples_finite": True,
        },
        "phase_zero_field": field_metrics(
            output / f"fields_{end_step}.csv"
        ),
    }
    return metrics, scalar_rows[-1], collision_rows[-1]


def execute(args: argparse.Namespace) -> dict:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise HorizonError(
            "horizon extension requires --acknowledge-cost " + ACKNOWLEDGEMENT
        )
    if args.additional_cycles > MAX_ADDITIONAL_CYCLES:
        raise HorizonError(
            f"--additional-cycles exceeds built-in limit "
            f"{MAX_ADDITIONAL_CYCLES}"
        )
    updates = (
        INITIAL_PARTICLES * STEPS_PER_CYCLE * args.additional_cycles
    )
    if args.max_initial_updates > HARD_ADDITIONAL_UPDATE_LIMIT:
        raise HorizonError(
            "--max-initial-updates exceeds built-in limit: "
            f"{args.max_initial_updates} > {HARD_ADDITIONAL_UPDATE_LIMIT}"
        )
    if updates > args.max_initial_updates:
        raise HorizonError(
            f"horizon exceeds --max-initial-updates: {updates} > "
            f"{args.max_initial_updates}"
        )
    executable = args.executable.resolve()
    prior_work = args.prior_work_dir.resolve()
    prior_report_path = args.prior_report.resolve()
    work = args.work_dir.resolve()
    report_path = args.report.resolve()
    if not executable.is_file():
        raise HorizonError(f"missing AuroraPIC CLI: {executable}")
    if work.exists():
        raise HorizonError(f"refusing to overwrite work directory: {work}")
    if report_path.exists():
        raise HorizonError(f"refusing to overwrite report: {report_path}")
    prior_report = read_json(prior_report_path)
    if (
        prior_report.get("turner_startup_report_version") != 1
        or prior_report.get("case_id") != "turner-helium-ccp-2013-case-1"
        or prior_report.get("work", {}).get("steps") != STEPS_PER_CYCLE
        or not prior_report.get("diagnostics", {}).get(
            "startup_checks_passed", False
        )
    ):
        raise HorizonError("prior report is not a completed one-cycle startup")
    checkpoint = prior_work / "checkpoint-one-cycle.apc"
    expected_checkpoint = prior_report.get("provenance", {}).get(
        "one_cycle_checkpoint_sha256"
    )
    if not checkpoint.is_file() or sha256(checkpoint) != expected_checkpoint:
        raise HorizonError("prior one-cycle checkpoint identity mismatch")
    prior_output = prior_work / "stage2-output"
    prior_scalars = rows(prior_output / "scalars.csv")
    prior_collisions = rows(prior_output / "collisions.csv")
    prior_scalar = prior_scalars[-1]
    prior_collision = prior_collisions[-1]
    if integer(prior_scalar, "step", "prior scalars") != STEPS_PER_CYCLE:
        raise HorizonError("prior diagnostics do not end at one RF cycle")

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
        raise HorizonError(
            "cannot generate exact campaign deck: " + generated.stderr.strip()
        )
    production_text = production.read_text(encoding="utf-8")
    environment = dict(os.environ)
    environment.update({
        "OMP_NUM_THREADS": "1",
        "OMP_DYNAMIC": "FALSE",
        "AURORA_OPENMP_THREADS": "1",
    })
    stages = []
    cycle_reports = []
    checkpoint_hashes = {}
    previous_checkpoint = checkpoint
    for cycle in range(2, 2 + args.additional_cycles):
        output = work / f"cycle-{cycle}-output"
        current_checkpoint = work / f"checkpoint-cycle-{cycle}.apc"
        config = work / f"cycle-{cycle}.cfg"
        write_new(
            config,
            stage_deck(
                production_text, cycle * STEPS_PER_CYCLE, output,
                current_checkpoint, previous_checkpoint,
            ),
        )
        start = time.perf_counter()
        try:
            result = subprocess.run(
                [
                    str(executable), "--allow-large-run",
                    CLI_ACKNOWLEDGEMENT, str(config),
                ],
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=args.timeout_seconds,
                preexec_fn=conservative_child_setup if os.name == "posix" else None,
            )
        except subprocess.TimeoutExpired as error:
            raise HorizonError(
                f"cycle {cycle} exceeded {args.timeout_seconds} seconds"
            ) from error
        elapsed = time.perf_counter() - start
        if result.returncode != 0:
            raise HorizonError(
                f"cycle {cycle} failed: {result.stdout.strip()} "
                f"{result.stderr.strip()}"
            )
        if not current_checkpoint.is_file():
            raise HorizonError(f"cycle {cycle} checkpoint is absent")
        metrics, prior_scalar, prior_collision = analyze_cycle(
            cycle, prior_scalar, prior_collision, output
        )
        cycle_reports.append(metrics)
        stages.append({
            "cycle": cycle,
            "wall_seconds": elapsed,
            "config_sha256": sha256(config),
        })
        checkpoint_hashes[str(cycle)] = sha256(current_checkpoint)
        previous_checkpoint = current_checkpoint

    electron_changes = [
        cycle["population"]["electron_relative_change"]
        for cycle in cycle_reports
    ]
    ion_changes = [
        cycle["population"]["ion_relative_change"]
        for cycle in cycle_reports
    ]
    electron_trend = (
        "decreasing" if all(value < 0 for value in electron_changes)
        else "increasing" if all(value > 0 for value in electron_changes)
        else "mixed"
    )
    ion_trend = (
        "decreasing" if all(value < 0 for value in ion_changes)
        else "increasing" if all(value > 0 for value in ion_changes)
        else "mixed"
    )
    return {
        "turner_horizon_report_version": 1,
        "case_id": "turner-helium-ccp-2013-case-1",
        "scope": (
            f"checkpointed_cycles_2_through_{1 + args.additional_cycles}"
        ),
        "physics_claim": "none",
        "steady_state_claim": False,
        "production_launch_authorized": False,
        "work": {
            "additional_cycles": args.additional_cycles,
            "additional_initial_particle_updates": updates,
            "stage_timings": stages,
            "total_wall_seconds": sum(
                stage["wall_seconds"] for stage in stages
            ),
        },
        "cycles": cycle_reports,
        "trend": {
            "electron_population": electron_trend,
            "ion_population": ion_trend,
            "final_electrons": cycle_reports[-1]["population"][
                "final_electrons"
            ],
            "final_ions": cycle_reports[-1]["population"]["final_ions"],
            "electron_fraction_of_original": (
                cycle_reports[-1]["population"]["final_electrons"] / 65536.0
            ),
            "maximum_last_cycle_population_relative_change": max(
                abs(electron_changes[-1]), abs(ion_changes[-1])
            ),
            "stationary": False,
            "interpretation":
                "short startup trend only; no stationarity acceptance applied",
        },
        "provenance": {
            "prior_report_sha256": sha256(prior_report_path),
            "prior_checkpoint_sha256": expected_checkpoint,
            "case_manifest_sha256": sha256(args.case_manifest.resolve()),
            "normalization_audit_sha256": sha256(
                args.normalized_dir.resolve() / "audit.json"
            ),
            "executable_sha256": sha256(executable),
            "production_config_sha256": sha256(production),
            "checkpoint_sha256_by_cycle": checkpoint_hashes,
            "platform": platform.platform(),
        },
        "warnings": [
            f"{1 + args.additional_cycles} RF cycles remain far shorter "
            "than the 1280-cycle benchmark.",
            "Population trends do not establish a stationary discharge.",
            "The published chi-squared acceptance range is inapplicable.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("normalized_dir", type=Path)
    parser.add_argument("--prior-work-dir", type=Path, required=True)
    parser.add_argument("--prior-report", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument(
        "--additional-cycles", type=positive_integer, default=3
    )
    parser.add_argument(
        "--max-initial-updates", type=positive_integer, default=160_000_000
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
    except (
        HorizonError, StartupError, OSError, UnicodeError,
        json.JSONDecodeError,
    ) as error:
        print(f"Turner horizon error: {error}", file=sys.stderr)
        return 2
    print(
        "Completed bounded checkpointed Turner horizon without a physics "
        f"claim: {args.report}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
