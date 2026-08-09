#!/usr/bin/env python3
"""Advance eduPIC through immutable, host-guarded native measurement blocks."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time

from advance_edupic_equilibration import (
    inspect_host_guard,
    nonnegative_integer,
    positive_float,
)
from run_edupic_measurement_stage import (
    HARD_CYCLE_LIMIT,
    HARD_INITIAL_PARTICLE_STEP_LIMIT,
    HARD_TIMEOUT_SECONDS,
)
from run_edupic_stage import (
    TIMESTEPS_PER_CYCLE,
    atomic_json,
    checkpoint_state,
    convergence_rows,
    sha256,
)


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_ADVANCES_BOUNDED_EDUPIC_MEASUREMENT"
HOST_GUARD_AMENDMENT_ACKNOWLEDGEMENT = (
    "I_UNDERSTAND_THIS_RELAXES_AN_EXISTING_CAMPAIGN_MEMORY_GUARD")
HARD_TOTAL_WALL_SECONDS = 600
HARD_STAGES_PER_INVOCATION = 16


class MeasurementAdvanceError(RuntimeError):
    pass


def apply_memory_guard_amendment(report: dict, limits: dict,
                                 args: argparse.Namespace) -> bool:
    """Audit an explicitly acknowledged reduction of only the RAM floor."""
    old_limits = report.get("limits")
    if not isinstance(old_limits, dict):
        return False
    old_guard = old_limits.get("host_health_guard")
    new_guard = limits.get("host_health_guard")
    if not isinstance(old_guard, dict) or not isinstance(new_guard, dict):
        return False
    old_memory = old_guard.get("minimum_available_memory_mib")
    new_memory = new_guard.get("minimum_available_memory_mib")
    if (not isinstance(old_memory, (int, float)) or
            not isinstance(new_memory, (int, float)) or
            new_memory >= old_memory):
        return False
    expected = json.loads(json.dumps(old_limits))
    expected["host_health_guard"]["minimum_available_memory_mib"] = new_memory
    if expected != limits:
        return False
    if (args.acknowledge_memory_guard_amendment !=
            HOST_GUARD_AMENDMENT_ACKNOWLEDGEMENT):
        raise MeasurementAdvanceError(
            "lowering an existing memory guard requires "
            "--acknowledge-memory-guard-amendment " +
            HOST_GUARD_AMENDMENT_ACKNOWLEDGEMENT)
    reason = args.memory_guard_amendment_reason
    if reason is None or not reason.strip():
        raise MeasurementAdvanceError(
            "lowering an existing memory guard requires a nonempty "
            "--memory-guard-amendment-reason")
    latest = report.get("latest_state", report.get("initial_state", {}))
    report.setdefault("operational_policy_amendments", []).append({
        "kind": "minimum_available_memory_mib_reduction",
        "applied_at_cycle": latest.get("cycles"),
        "previous_value_mib": old_memory,
        "new_value_mib": new_memory,
        "reason": reason.strip(),
        "scientific_analysis_contract_changed": False,
    })
    report["limits"] = limits
    return True


def positive_integer(text: str) -> int:
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return value


def diagnostic_hashes(stage_dir: Path, stage: dict) -> None:
    outputs = stage.get("outputs", {})
    simple = {
        "stdout.txt": outputs.get("stdout_sha256"),
        "stderr.txt": outputs.get("stderr_sha256"),
        "conv.dat": outputs.get("convergence_sha256"),
    }
    diagnostics = outputs.get("diagnostics", {})
    if not isinstance(diagnostics, dict) or not diagnostics:
        raise MeasurementAdvanceError(
            f"measurement diagnostics are absent in {stage_dir}")
    for name, metadata in diagnostics.items():
        if not isinstance(metadata, dict):
            raise MeasurementAdvanceError(
                f"invalid diagnostic metadata for {name} in {stage_dir}")
        simple[name] = metadata.get("sha256")
    for name, expected in simple.items():
        path = stage_dir / name
        if (not isinstance(expected, str) or len(expected) != 64 or
                not path.is_file() or sha256(path) != expected):
            raise MeasurementAdvanceError(
                f"measurement output hash differs for {name} in {stage_dir}")


def inspect_stage(stage_dir: Path, current: dict,
                  binary_sha256: str) -> tuple[dict, dict]:
    report_path = stage_dir / "stage-report.json"
    try:
        stage = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise MeasurementAdvanceError(
            f"cannot inspect completed stage {stage_dir}: {error}") from error
    if (stage.get("completed") is not True or
            stage.get("scope") != "bounded_external_native_measurement_stage" or
            stage.get("source_binary", {}).get("sha256") != binary_sha256 or
            stage.get("initial_state", {}).get("sha256") != current["sha256"] or
            stage.get("stage", {}).get("measurement_mode") is not True):
        raise MeasurementAdvanceError(
            f"completed measurement-stage contract differs in {stage_dir}")
    final = stage.get("final_state", {})
    on_disk = checkpoint_state(stage_dir / "picdata.bin")
    requested = stage.get("stage", {}).get("requested_cycles")
    if (final != on_disk or not isinstance(requested, int) or requested <= 0 or
            final.get("cycles") != current["cycles"] + requested):
        raise MeasurementAdvanceError(
            f"completed measurement checkpoint differs in {stage_dir}")
    history = convergence_rows(stage_dir / "conv.dat")
    if (history[-1]["cycle"] != final["cycles"] or
            history[-1]["electrons"] != final["electrons"] or
            history[-1]["ions"] != final["ions"]):
        raise MeasurementAdvanceError(
            f"measurement convergence state differs in {stage_dir}")
    diagnostic_hashes(stage_dir, stage)
    return stage, final


def stage_summary(stage_dir: Path, stage: dict, current: dict,
                  predicted_seconds: float | None,
                  recovered: bool = False) -> dict:
    final = stage["final_state"]
    result = {
        "start_cycle": current["cycles"],
        "end_cycle": final["cycles"],
        "measurement_cycles": final["cycles"] - current["cycles"],
        "wall_seconds": stage["stage"]["wall_seconds"],
        "predicted_wall_seconds_with_safety_factor": predicted_seconds,
        "initial_total_particles": current["total_particles"],
        "final_total_particles": final["total_particles"],
        "input_checkpoint_sha256": current["sha256"],
        "output_checkpoint_sha256": final["sha256"],
        "stage_report_sha256": sha256(stage_dir / "stage-report.json"),
        "diagnostic_output_bytes": sum(
            (stage_dir / name).stat().st_size
            for name in stage["outputs"]["diagnostics"]),
    }
    if recovered:
        result["recovered_after_coordinator_interruption"] = True
    return result


def expected_limits(args: argparse.Namespace) -> dict:
    host_policy = {
        "maximum_load_per_cpu": args.max_host_load_per_cpu,
        "minimum_available_memory_mib": args.min_available_memory_mib,
        "maximum_swap_io_pages_per_stage": args.max_swap_io_pages_per_stage,
        "minimum_free_disk_mib": args.min_free_disk_mib,
    }
    host_policy = {key: value for key, value in host_policy.items()
                   if value is not None}
    limits = {
        "maximum_total_wall_seconds": args.max_wall_seconds,
        "maximum_stage_wall_seconds": args.stage_timeout_seconds,
        "measurement_block_cycles": args.block_cycles,
        "maximum_stage_initial_particle_steps":
            args.max_stage_initial_particle_steps,
        "maximum_stages_per_invocation": args.max_stages_per_invocation,
        "qualified_seconds_per_cycle": args.qualified_seconds_per_cycle,
        "wall_prediction_safety_factor": args.wall_prediction_safety_factor,
    }
    if host_policy:
        limits["host_health_guard"] = host_policy
    return limits


def inspect_campaign_host_guard(
        report: dict, policy: dict, phase: str, cycle: int,
        previous_swap_total: int | None,
        campaign_dir: Path) -> tuple[list[str], int | None]:
    violations, swap_total = inspect_host_guard(
        report, policy, phase, cycle, previous_swap_total)
    sample = report["host_health_checks"][-1]
    minimum_disk = policy.get("minimum_free_disk_mib")
    if minimum_disk is not None:
        try:
            free_disk_mib = shutil.disk_usage(campaign_dir).free / (1024.0 ** 2)
        except OSError as error:
            raise MeasurementAdvanceError(
                f"cannot read campaign filesystem capacity: {error}") from error
        sample["free_disk_mib"] = free_disk_mib
        if free_disk_mib < minimum_disk:
            violations.append("host_free_disk_below_minimum")
            sample["violations"] = violations
    return violations, swap_total


def advance(args: argparse.Namespace) -> dict:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise MeasurementAdvanceError(
            "campaign requires --acknowledge-cost " + ACKNOWLEDGEMENT)
    if args.max_wall_seconds > HARD_TOTAL_WALL_SECONDS:
        raise MeasurementAdvanceError("overall wall-time limit exceeds built-in ceiling")
    if args.stage_timeout_seconds > HARD_TIMEOUT_SECONDS:
        raise MeasurementAdvanceError("stage timeout exceeds built-in ceiling")
    if args.block_cycles > HARD_CYCLE_LIMIT:
        raise MeasurementAdvanceError("measurement block exceeds built-in cycle limit")
    if args.max_stage_initial_particle_steps > HARD_INITIAL_PARTICLE_STEP_LIMIT:
        raise MeasurementAdvanceError("particle-step limit exceeds built-in ceiling")
    if args.max_stages_per_invocation > HARD_STAGES_PER_INVOCATION:
        raise MeasurementAdvanceError("stage-count limit exceeds built-in ceiling")

    executable = args.executable.resolve()
    input_dir = args.input_state_dir.resolve()
    campaign_dir = args.campaign_dir.resolve()
    binary_sha256 = args.expected_binary_sha256.lower()
    if not executable.is_file() or sha256(executable) != binary_sha256:
        raise MeasurementAdvanceError(
            "external binary is missing or differs from locked SHA-256")
    initial = checkpoint_state(input_dir / "picdata.bin")
    if initial["sha256"] != args.expected_input_sha256.lower():
        raise MeasurementAdvanceError("input checkpoint differs from locked SHA-256")
    target_cycle = initial["cycles"] + args.target_measurement_cycles
    limits = expected_limits(args)
    host_policy = limits.get("host_health_guard", {})
    manifest = campaign_dir / "campaign-report.json"
    new_report = {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "bounded_native_measurement_campaign",
        "physics_claim": "none",
        "claim_boundary": (
            "Native measurement blocks preserve recoverability and block evidence. "
            "They do not establish converged statistics or cross-code agreement."),
        "source_binary_sha256": binary_sha256,
        "initial_state": initial,
        "target_measurement_cycles": args.target_measurement_cycles,
        "target_cycle": target_cycle,
        "limits": limits,
        "stages": [],
        "completed": False,
        "target_reached": False,
        "stop_reason": "campaign_started",
    }
    current_dir = input_dir
    current = initial
    if campaign_dir.exists():
        if not args.resume_existing:
            raise MeasurementAdvanceError(
                f"refusing to overwrite campaign directory: {campaign_dir}")
        try:
            report = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise MeasurementAdvanceError(
                f"cannot resume campaign manifest: {error}") from error
        fixed_contract_differs = (
                report.get("scope") != new_report["scope"] or
                report.get("source_binary_sha256") != binary_sha256 or
                report.get("initial_state", {}).get("sha256") != initial["sha256"] or
                report.get("target_measurement_cycles") !=
                    args.target_measurement_cycles or
                report.get("target_cycle") != target_cycle)
        if fixed_contract_differs:
            raise MeasurementAdvanceError(
                "resume arguments differ from the measurement campaign contract")
        if report.get("limits") != limits:
            if not apply_memory_guard_amendment(report, limits, args):
                raise MeasurementAdvanceError(
                    "resume arguments differ from the measurement campaign contract")
            atomic_json(manifest, report)
        for recorded in report.get("stages", []):
            if recorded.get("start_cycle") != current["cycles"]:
                raise MeasurementAdvanceError(
                    "recorded measurement-stage chain is discontinuous")
            stage_dir = campaign_dir / (
                f"stage-{recorded['start_cycle']:06d}-{recorded['end_cycle']:06d}")
            stage, final = inspect_stage(stage_dir, current, binary_sha256)
            expected = stage_summary(
                stage_dir, stage, current,
                recorded.get("predicted_wall_seconds_with_safety_factor"),
                recovered=recorded.get(
                    "recovered_after_coordinator_interruption", False))
            if recorded != expected:
                raise MeasurementAdvanceError(
                    "recorded measurement-stage summary differs from its report")
            current_dir, current = stage_dir, final
        recovered = 0
        while current["cycles"] < target_cycle:
            candidates = sorted(
                path for path in campaign_dir.glob(
                    f"stage-{current['cycles']:06d}-*")
                if (path / "stage-report.json").is_file())
            if not candidates:
                break
            if len(candidates) != 1:
                raise MeasurementAdvanceError(
                    "multiple unrecorded measurement stages begin at the latest cycle")
            stage_dir = candidates[0]
            stage, final = inspect_stage(stage_dir, current, binary_sha256)
            if final["cycles"] > target_cycle:
                raise MeasurementAdvanceError(
                    "unrecorded measurement stage extends beyond campaign target")
            report["stages"].append(stage_summary(
                stage_dir, stage, current, None, recovered=True))
            current_dir, current = stage_dir, final
            report["latest_state"] = current
            recovered += 1
            atomic_json(manifest, report)
        report.pop("failure", None)
        report["recovered_unrecorded_stages"] = (
            report.get("recovered_unrecorded_stages", 0) + recovered)
        report["completed"] = False
        report["stop_reason"] = "resume_reconciled"
        atomic_json(manifest, report)
    else:
        if args.resume_existing:
            raise MeasurementAdvanceError(
                "cannot resume a campaign directory that does not exist")
        campaign_dir.mkdir(parents=True)
        report = new_report
        atomic_json(manifest, report)

    runner = Path(__file__).resolve().with_name(
        "run_edupic_measurement_stage.py")
    started = time.perf_counter()
    previous_swap_total = None
    if report.get("host_health_checks"):
        previous_swap_total = report["host_health_checks"][-1].get(
            "swap_io_pages_total")
    stages_this_invocation = 0
    while (current["cycles"] < target_cycle and
           stages_this_invocation < args.max_stages_per_invocation):
        elapsed = time.perf_counter() - started
        remaining_wall = args.max_wall_seconds - elapsed
        if remaining_wall < 6.0:
            report["stop_reason"] = "overall_wall_time_exhausted"
            break
        cycles = min(args.block_cycles, target_cycle - current["cycles"])
        work = current["total_particles"] * cycles * TIMESTEPS_PER_CYCLE
        if work > args.max_stage_initial_particle_steps:
            report["stop_reason"] = "stage_particle_step_budget_too_small"
            break
        observed_seconds_per_cycle = args.qualified_seconds_per_cycle
        if report["stages"]:
            previous = report["stages"][-1]
            observed_seconds_per_cycle = max(
                observed_seconds_per_cycle,
                previous["wall_seconds"] / previous["measurement_cycles"])
        predicted_seconds = (
            observed_seconds_per_cycle * cycles *
            args.wall_prediction_safety_factor)
        if remaining_wall < predicted_seconds + 5.0:
            report["stop_reason"] = "insufficient_predicted_wall_time"
            report["next_stage_prediction"] = {
                "planned_cycles": cycles,
                "observed_seconds_per_cycle_basis": observed_seconds_per_cycle,
                "predicted_wall_seconds_with_safety_factor": predicted_seconds,
                "remaining_campaign_wall_seconds": remaining_wall,
            }
            break
        if args.stage_timeout_seconds < math.ceil(predicted_seconds):
            report["stop_reason"] = "insufficient_stage_timeout"
            break
        if host_policy:
            violations, previous_swap_total = inspect_campaign_host_guard(
                report, host_policy, "before_stage", current["cycles"],
                previous_swap_total, campaign_dir)
            atomic_json(manifest, report)
            if violations:
                report["stop_reason"] = violations[0]
                break
        expected_end = current["cycles"] + cycles
        stage_dir = campaign_dir / (
            f"stage-{current['cycles']:06d}-{expected_end:06d}")
        command = [
            sys.executable, str(runner), str(executable), str(current_dir),
            str(stage_dir), "--cycles", str(cycles),
            "--expected-binary-sha256", binary_sha256,
            "--expected-input-sha256", current["sha256"],
            "--timeout-seconds", str(args.stage_timeout_seconds),
            "--max-initial-particle-steps",
            str(args.max_stage_initial_particle_steps),
            "--acknowledge-cost",
            "I_UNDERSTAND_THIS_IS_A_BOUNDED_EDUPIC_MEASUREMENT_STAGE",
        ]
        completed = subprocess.run(
            command, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=args.stage_timeout_seconds + 5)
        if completed.returncode != 0:
            report["stop_reason"] = "stage_failed"
            report["failure"] = {
                "planned_start_cycle": current["cycles"],
                "planned_end_cycle": expected_end,
                "stderr": completed.stderr[-4000:],
            }
            atomic_json(manifest, report)
            raise MeasurementAdvanceError(
                f"bounded measurement stage failed; retained under {stage_dir}")
        stage, final = inspect_stage(stage_dir, current, binary_sha256)
        report["stages"].append(stage_summary(
            stage_dir, stage, current, predicted_seconds))
        current_dir, current = stage_dir, final
        report["latest_state"] = current
        report["stop_reason"] = "stage_completed"
        stages_this_invocation += 1
        atomic_json(manifest, report)
        if host_policy:
            violations, previous_swap_total = inspect_campaign_host_guard(
                report, host_policy, "after_stage", current["cycles"],
                previous_swap_total, campaign_dir)
            atomic_json(manifest, report)
            if violations:
                report["stop_reason"] = violations[0]
                break
    if (current["cycles"] < target_cycle and
            stages_this_invocation >= args.max_stages_per_invocation and
            report["stop_reason"] == "stage_completed"):
        report["stop_reason"] = "stage_count_limit_reached"
    invocation_wall = time.perf_counter() - started
    report.setdefault("invocations", []).append({
        "resume_existing": args.resume_existing,
        "completed_stages": stages_this_invocation,
        "coordinator_wall_seconds": invocation_wall,
    })
    report["total_coordinator_wall_seconds_recorded"] = sum(
        value["coordinator_wall_seconds"] for value in report["invocations"])
    report["total_stage_wall_seconds"] = sum(
        value["wall_seconds"] for value in report["stages"])
    report["completed_measurement_cycles"] = (
        current["cycles"] - initial["cycles"])
    report["target_reached"] = current["cycles"] >= target_cycle
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
    parser.add_argument("--target-measurement-cycles", type=positive_integer,
                        required=True)
    parser.add_argument("--block-cycles", type=positive_integer, default=4)
    parser.add_argument("--qualified-seconds-per-cycle", type=positive_float,
                        required=True)
    parser.add_argument("--wall-prediction-safety-factor", type=positive_float,
                        default=1.5)
    parser.add_argument("--max-wall-seconds", type=positive_integer, default=600)
    parser.add_argument("--stage-timeout-seconds", type=positive_integer,
                        default=300)
    parser.add_argument("--max-stage-initial-particle-steps", type=positive_integer,
                        default=8_000_000_000)
    parser.add_argument("--max-stages-per-invocation", type=positive_integer,
                        default=1)
    parser.add_argument("--max-host-load-per-cpu", type=positive_float)
    parser.add_argument("--min-available-memory-mib", type=positive_float)
    parser.add_argument("--max-swap-io-pages-per-stage", type=nonnegative_integer)
    parser.add_argument("--min-free-disk-mib", type=positive_float)
    parser.add_argument("--acknowledge-memory-guard-amendment")
    parser.add_argument("--memory-guard-amendment-reason")
    parser.add_argument("--acknowledge-cost")
    parser.add_argument("--resume-existing", action="store_true")
    args = parser.parse_args()
    try:
        report = advance(args)
    except (MeasurementAdvanceError, OSError, subprocess.TimeoutExpired,
            json.JSONDecodeError) as error:
        print(f"eduPIC measurement advance rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
