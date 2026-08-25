#!/usr/bin/env python3
"""Run one locked second block of matched-half-step threshold traffic."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from run_aurorapic_edupic_pilot import (
    PilotError, atomic_json, atomic_text, available_memory_kib, run_process,
    sha256,
)
from run_aurorapic_ionizing_tail_block import analyze_output
from run_aurorapic_matched_half_step_thresholds import measurement_deck


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_MATCHED_HALF_STEP_REPLICATION"
CLI_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LARGE_RUN"
APPROVED_RULE_SHA256S = {
    "d0af3ae0936353b4b904f3e77be6fc4ec214037f614832da18353756c020298d",
}


class ReplicationError(RuntimeError):
    pass


def execute(args: argparse.Namespace) -> dict[str, object]:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise ReplicationError("matched-half-step replication was not acknowledged")
    rule_path = args.rule.resolve()
    rule_hash = sha256(rule_path)
    if rule_hash not in APPROVED_RULE_SHA256S:
        raise ReplicationError("rule is not an approved matched-half-step replication")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    matches = [item for item in rule["locked_continuation_states"]
               if item["id"] == args.initial_state_id]
    if len(matches) != 1:
        raise ReplicationError("continuation-state id is not locked")
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
            raise ReplicationError(f"locked {label} differs")
    previous = json.loads(prior_report.read_text(encoding="utf-8"))
    if (previous.get("all_gates_passed") is not True or
            previous.get("inputs", {}).get("initial_state_id") != state["id"] or
            previous.get("final_checkpoint_sha256") != state["checkpoint_sha256"] or
            previous.get("algorithm_contract", {}).get(
                "collision_velocity_sampling") != "leapfrog_half_step"):
        raise ReplicationError("prior branch provenance differs")
    execution = rule["execution_contract"]
    if (int(execution["end_step"]) - int(execution["start_step"]) !=
            int(execution["measurement_cycles"]) *
            int(execution["steps_per_cycle"])):
        raise ReplicationError("replication horizon is inconsistent")
    work = args.work_dir.resolve()
    if work.exists():
        raise ReplicationError(f"refusing to overwrite {work}")
    available_memory = available_memory_kib()
    available_disk = shutil.disk_usage(work.parent).free // 1024
    if available_memory < int(execution["minimum_available_memory_kib"]):
        raise ReplicationError("available memory is below the launch floor")
    if available_disk < int(execution["minimum_available_disk_kib"]):
        raise ReplicationError("available disk is below the launch floor")
    work.mkdir(parents=True)
    output = work / "output"
    deck = work / "input.cfg"
    atomic_text(deck, measurement_deck(
        base_path.read_text(encoding="utf-8"), rule, state,
        output, checkpoint))
    resources = run_process([
        str(executable), "--allow-large-run", CLI_ACKNOWLEDGEMENT,
        str(deck)], work / "stdout.txt", work / "stderr.txt",
        float(execution["timeout_seconds"]))
    result = analyze_output(output, rule, resources)
    result.update({
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "aurorapic_matched_half_step_threshold_replication_block",
        "rule_sha256": rule_hash,
        "inputs": {
            "initial_state_id": state["id"],
            "solver_sha256": sha256(executable),
            "base_config_sha256": sha256(base_path),
            "input_checkpoint_sha256": sha256(checkpoint),
            "prior_report_sha256": sha256(prior_report),
            "deck_sha256": sha256(deck),
        },
        "algorithm_contract": {
            "collision_velocity_sampling": "leapfrog_half_step",
            "diagnostics_reset_on_restart": True,
            "random_stream_continued_from_checkpoint": True,
        },
        "resources": {
            **resources,
            "available_memory_before_launch_kib": available_memory,
            "available_disk_before_launch_kib": available_disk,
        },
        "claim_boundary": rule["claim_boundary"],
        "physics_claim": rule["physics_claim"],
    })
    atomic_json(work / "block-report.json", result)
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
    except (ReplicationError, PilotError, OSError, ValueError,
            KeyError) as error:
        parser.error(str(error))
    print(json.dumps({
        "report": str(args.work_dir.resolve() / "block-report.json"),
        "all_gates_passed": result["all_gates_passed"],
        "sampling": result["sampling"],
        "resources": result["resources"],
    }, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
