#!/usr/bin/env python3
"""Run one locked electrode half-cell Gauss-correction control branch."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from run_aurorapic_edupic_pilot import (
    PilotError, atomic_json, atomic_text, available_memory_kib, run_process,
    set_global, sha256,
)
from run_aurorapic_ionizing_tail_block import analyze_output
from run_aurorapic_matched_half_step_thresholds import measurement_deck


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_AN_ENDPOINT_GAUSS_CONTROL"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"
APPROVED_RULE_SHA256S = {
    "ca162b2a1d0207ec6fc69a1b0a8e0bf2924538c6c27eb1e69f300fcc41852d70",
}


class EndpointGaussError(RuntimeError):
    pass


def relaxation_deck(base: str, rule: dict[str, object],
                    state: dict[str, object], output: Path,
                    checkpoint: Path) -> str:
    execution = rule["execution_contract"]
    diagnostic = rule["diagnostic_contract"]
    end = int(execution["relaxation_end_step"])
    result = measurement_deck(base, rule, state, output, checkpoint)
    for key, value in {
        "steps": end,
        "spatial_average_interval": diagnostic[
            "spatial_average_interval_steps"],
        "spatial_average_start_step": int(execution["initial_step"]) + 1,
        "spatial_average_end_step": end,
        "spatial_average_rf_cycles": execution["relaxation_cycles"],
        "spatial_average_phase_bins": diagnostic["phase_bins"],
        "checkpoint_interval": end,
    }.items():
        result = set_global(result, key, str(value))
    return result


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise EndpointGaussError("endpoint Gauss control was not acknowledged")
    rule_path = args.rule.resolve()
    rule_hash = sha256(rule_path)
    if rule_hash not in APPROVED_RULE_SHA256S:
        raise EndpointGaussError("rule is not an approved endpoint control")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    matches = [item for item in rule["locked_continuation_states"]
               if item["id"] == args.initial_state_id]
    if len(matches) != 1:
        raise EndpointGaussError("continuation-state id is not locked")
    state = matches[0]
    executable = args.executable.resolve()
    base_path = args.base_config.resolve()
    checkpoint = args.checkpoint.resolve()
    prior_report = args.prior_report.resolve()
    for path, expected, label in (
        (executable, state["solver_sha256"], "solver"),
        (base_path, state["base_config_sha256"], "base config"),
        (checkpoint, state["checkpoint_sha256"], "checkpoint"),
        (prior_report, state["prior_report_sha256"], "prior report"),
    ):
        if sha256(path) != expected:
            raise EndpointGaussError(f"locked {label} differs")
    previous = json.loads(prior_report.read_text(encoding="utf-8"))
    if (previous.get("all_gates_passed") is not True or
            previous.get("inputs", {}).get("initial_state_id") != state["id"] or
            previous.get("final_checkpoint_sha256") != state["checkpoint_sha256"]):
        raise EndpointGaussError("prior branch provenance differs")
    execution = rule["execution_contract"]
    if (int(execution["relaxation_end_step"]) -
            int(execution["initial_step"]) !=
            int(execution["relaxation_cycles"]) *
            int(execution["steps_per_cycle"]) or
            int(execution["end_step"]) - int(execution["start_step"]) !=
            int(execution["measurement_cycles"]) *
            int(execution["steps_per_cycle"])):
        raise EndpointGaussError("execution horizon is inconsistent")
    work = args.work_dir.resolve()
    if work.exists():
        raise EndpointGaussError(f"refusing to overwrite {work}")
    available_memory = available_memory_kib()
    available_disk = shutil.disk_usage(work.parent).free // 1024
    if available_memory < int(execution["minimum_available_memory_kib"]):
        raise EndpointGaussError("available memory is below the launch floor")
    if available_disk < int(execution["minimum_available_disk_kib"]):
        raise EndpointGaussError("available disk is below the launch floor")
    work.mkdir(parents=True)
    base = base_path.read_text(encoding="utf-8")

    relaxation = work / "relaxation"
    relaxation.mkdir()
    relaxation_output = relaxation / "output"
    relaxation_input = relaxation / "input.cfg"
    atomic_text(relaxation_input, relaxation_deck(
        base, rule, state, relaxation_output, checkpoint))
    relaxation_resources = run_process([
        str(executable), "--allow-large-run", CLI_ACKNOWLEDGEMENT,
        str(relaxation_input)], relaxation / "stdout.txt",
        relaxation / "stderr.txt", float(execution["timeout_seconds"]))
    relaxed_checkpoint = relaxation_output / (
        f"checkpoint_{int(execution['relaxation_end_step'])}.apc")
    if not relaxed_checkpoint.is_file():
        raise EndpointGaussError("relaxation checkpoint is missing")

    measurement = work / "measurement"
    measurement.mkdir()
    measurement_output = measurement / "output"
    measurement_input = measurement / "input.cfg"
    atomic_text(measurement_input, measurement_deck(
        base, rule, state, measurement_output, relaxed_checkpoint))
    measurement_resources = run_process([
        str(executable), "--allow-large-run", CLI_ACKNOWLEDGEMENT,
        str(measurement_input)], measurement / "stdout.txt",
        measurement / "stderr.txt", float(execution["timeout_seconds"]))
    peak = max(int(relaxation_resources["peak_resident_set_kib"]),
               int(measurement_resources["peak_resident_set_kib"]))
    result = analyze_output(measurement_output, rule, {
        **measurement_resources, "peak_resident_set_kib": peak})
    result.update({
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "aurorapic_electrode_half_cell_gauss_control_branch",
        "rule_sha256": rule_hash,
        "inputs": {
            "initial_state_id": state["id"],
            "solver_sha256": sha256(executable),
            "base_config_sha256": sha256(base_path),
            "input_checkpoint_sha256": sha256(checkpoint),
            "prior_report_sha256": sha256(prior_report),
            "relaxation_deck_sha256": sha256(relaxation_input),
            "measurement_deck_sha256": sha256(measurement_input),
        },
        "algorithm_contract": {
            "collision_velocity_sampling": "leapfrog_half_step",
            "random_stream_continued_from_checkpoint": True,
            "electrode_half_cell_gauss_correction": True,
            "rf_phase_rad": rule["locked_inputs"]["baseline_phase_rad"],
            "relaxation_cycles_after_solver_change":
                execution["relaxation_cycles"],
        },
        "resources": {
            "relaxation": relaxation_resources,
            "measurement": measurement_resources,
            "maximum_peak_resident_set_kib": peak,
            "available_memory_before_launch_kib": available_memory,
            "available_disk_before_launch_kib": available_disk,
        },
        "relaxation_checkpoint_sha256": sha256(relaxed_checkpoint),
        "claim_boundary": rule["claim_boundary"],
        "physics_claim": rule["physics_claim"],
    })
    atomic_json(work / "branch-report.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("executable", type=Path)
    parser.add_argument("base_config", type=Path)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("prior_report", type=Path)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("--initial-state-id", required=True)
    parser.add_argument("--acknowledge-cost", required=True)
    args = parser.parse_args()
    try:
        result = execute(args)
    except (EndpointGaussError, PilotError, OSError, ValueError,
            KeyError) as error:
        parser.error(str(error))
    print(json.dumps({
        "report": str(args.work_dir.resolve() / "branch-report.json"),
        "all_gates_passed": result["all_gates_passed"],
        "sampling": result["sampling"],
        "resources": result["resources"],
    }, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
