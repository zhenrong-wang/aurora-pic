#!/usr/bin/env python3
"""Advance eduPIC through adaptive, hash-locked, bounded stages."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

from run_edupic_stage import checkpoint_state, sha256, TIMESTEPS_PER_CYCLE


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_ADVANCES_BOUNDED_EDUPIC_EQUILIBRATION"
HARD_TOTAL_WALL_SECONDS = 600
HARD_STAGE_CYCLES = 16
HARD_STAGE_PARTICLE_STEPS = 1_000_000_000


class AdvanceError(RuntimeError):
    pass


def positive_integer(text: str) -> int:
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return value


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    os.replace(temporary, path)


def advance(args: argparse.Namespace) -> dict:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise AdvanceError("campaign requires --acknowledge-cost " + ACKNOWLEDGEMENT)
    if args.max_wall_seconds > HARD_TOTAL_WALL_SECONDS:
        raise AdvanceError("overall wall-time limit exceeds built-in ceiling")
    if args.max_stage_cycles > HARD_STAGE_CYCLES:
        raise AdvanceError("stage cycle limit exceeds built-in ceiling")
    if args.max_stage_initial_particle_steps > HARD_STAGE_PARTICLE_STEPS:
        raise AdvanceError("stage particle-step limit exceeds built-in ceiling")
    executable = args.executable.resolve()
    input_dir = args.input_state_dir.resolve()
    campaign_dir = args.campaign_dir.resolve()
    if campaign_dir.exists():
        raise AdvanceError(f"refusing to overwrite campaign directory: {campaign_dir}")
    if not executable.is_file() or sha256(executable) != args.expected_binary_sha256.lower():
        raise AdvanceError("external binary is missing or differs from locked SHA-256")
    initial = checkpoint_state(input_dir / "picdata.bin")
    if initial["sha256"] != args.expected_input_sha256.lower():
        raise AdvanceError("input checkpoint differs from locked SHA-256")
    if args.target_cycle <= initial["cycles"]:
        raise AdvanceError("target cycle must be after the input checkpoint")

    campaign_dir.mkdir(parents=True)
    report = {
        "schema_version": 1, "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "bounded_adaptive_equilibration_advance", "physics_claim": "none",
        "source_binary_sha256": args.expected_binary_sha256.lower(),
        "initial_state": initial, "target_cycle": args.target_cycle,
        "limits": {"maximum_total_wall_seconds": args.max_wall_seconds,
                   "maximum_stage_wall_seconds": args.stage_timeout_seconds,
                   "maximum_stage_cycles": args.max_stage_cycles,
                   "maximum_stage_initial_particle_steps":
                       args.max_stage_initial_particle_steps},
        "stages": [], "completed": False, "target_reached": False,
        "stop_reason": "campaign_started",
    }
    manifest = campaign_dir / "campaign-report.json"
    atomic_json(manifest, report)
    runner = Path(__file__).resolve().with_name("run_edupic_stage.py")
    current_dir = input_dir
    current = initial
    started = time.perf_counter()
    while current["cycles"] < args.target_cycle:
        elapsed = time.perf_counter() - started
        remaining_wall = args.max_wall_seconds - elapsed
        # Reserve five seconds for runner startup, validation, and teardown so
        # the outer watchdog never extends beyond the campaign wall ceiling.
        if remaining_wall < 6.0:
            report["stop_reason"] = "overall_wall_time_exhausted"
            break
        work_per_cycle = current["total_particles"] * TIMESTEPS_PER_CYCLE
        affordable_cycles = args.max_stage_initial_particle_steps // work_per_cycle
        cycles = min(args.max_stage_cycles, affordable_cycles,
                     args.target_cycle - current["cycles"])
        if cycles < 1:
            report["stop_reason"] = "stage_particle_step_budget_too_small"
            break
        expected_end = current["cycles"] + cycles
        stage_dir = campaign_dir / f"stage-{current['cycles']:06d}-{expected_end:06d}"
        predicted_seconds = None
        if report["stages"]:
            previous = report["stages"][-1]
            predicted_seconds = (previous["wall_seconds"] /
                                 previous["cycles"] * cycles * 1.5)
            if remaining_wall < predicted_seconds + 5.0:
                report["stop_reason"] = "insufficient_predicted_wall_time"
                report["next_stage_prediction"] = {
                    "planned_cycles": cycles,
                    "predicted_wall_seconds_with_safety_factor": predicted_seconds,
                    "remaining_campaign_wall_seconds": remaining_wall,
                }
                break
        timeout_seconds = min(args.stage_timeout_seconds,
                              max(1, int(math.floor(remaining_wall - 5.0))))
        if predicted_seconds is not None and timeout_seconds < predicted_seconds:
            report["stop_reason"] = "insufficient_stage_timeout"
            break
        command = [
            sys.executable, str(runner), str(executable), str(current_dir),
            str(stage_dir), "--cycles", str(cycles),
            "--expected-binary-sha256", args.expected_binary_sha256.lower(),
            "--expected-input-sha256", current["sha256"],
            "--timeout-seconds", str(timeout_seconds),
            "--max-initial-particle-steps",
            str(args.max_stage_initial_particle_steps),
            "--acknowledge-cost", "I_UNDERSTAND_THIS_IS_A_BOUNDED_EDUPIC_STAGE",
        ]
        stage = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               timeout=timeout_seconds + 3)
        if stage.returncode != 0:
            report["stop_reason"] = "stage_failed"
            report["failure"] = {"planned_start_cycle": current["cycles"],
                                 "planned_end_cycle": expected_end,
                                 "stderr": stage.stderr[-4000:]}
            atomic_json(manifest, report)
            raise AdvanceError(f"bounded stage failed; retained under {stage_dir}")
        stage_report_path = stage_dir / "stage-report.json"
        stage_report = json.loads(stage_report_path.read_text(encoding="utf-8"))
        final = stage_report["final_state"]
        if final["cycles"] != expected_end:
            raise AdvanceError("completed stage report has unexpected final cycle")
        report["stages"].append({
            "start_cycle": current["cycles"], "end_cycle": final["cycles"],
            "cycles": cycles, "wall_seconds": stage_report["stage"]["wall_seconds"],
            "predicted_wall_seconds_with_safety_factor": predicted_seconds,
            "initial_total_particles": current["total_particles"],
            "final_total_particles": final["total_particles"],
            "input_checkpoint_sha256": current["sha256"],
            "output_checkpoint_sha256": final["sha256"],
            "stage_report_sha256": sha256(stage_report_path),
        })
        current_dir = stage_dir
        current = final
        report["latest_state"] = current
        report["stop_reason"] = "stage_completed"
        atomic_json(manifest, report)
    report["total_wall_seconds"] = time.perf_counter() - started
    report["target_reached"] = current["cycles"] >= args.target_cycle
    report["completed"] = True
    if report["target_reached"]:
        report["stop_reason"] = "target_cycle_reached"
    atomic_json(manifest, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("input_state_dir", type=Path)
    parser.add_argument("campaign_dir", type=Path)
    parser.add_argument("--expected-binary-sha256", required=True)
    parser.add_argument("--expected-input-sha256", required=True)
    parser.add_argument("--target-cycle", type=positive_integer, required=True)
    parser.add_argument("--max-wall-seconds", type=positive_integer, default=60)
    parser.add_argument("--stage-timeout-seconds", type=positive_integer, default=30)
    parser.add_argument("--max-stage-cycles", type=positive_integer, default=4)
    parser.add_argument("--max-stage-initial-particle-steps", type=positive_integer,
                        default=250_000_000)
    parser.add_argument("--acknowledge-cost")
    args = parser.parse_args()
    try:
        report = advance(args)
    except (AdvanceError, OSError, subprocess.TimeoutExpired,
            json.JSONDecodeError) as error:
        print(f"eduPIC advance rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
