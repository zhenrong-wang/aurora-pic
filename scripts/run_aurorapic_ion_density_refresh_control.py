#!/usr/bin/env python3
"""Run preregistered held-density common-state horizons serially."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import resource
import subprocess
import time

from run_aurorapic_common_state_trace import (
    RunError, deck as baseline_deck, sha256, summarize_member)


def control_deck(rule: dict[str, object], state: Path, horizon: int,
                 output: Path) -> str:
    text = baseline_deck(rule, state, horizon, output)
    marker = "collision_velocity_sampling = leapfrog_half_step\n"
    if text.count(marker) != 1:
        raise RunError("baseline deck collision-staggering marker differs")
    return text.replace(
        marker,
        marker + "subcycle_charge_deposition = pre_push_held\n")


def execute(args: argparse.Namespace) -> dict[str, object]:
    control_rule_path = args.control_rule.resolve()
    parent_rule_path = args.parent_rule.resolve()
    state = args.state.resolve()
    binary = args.binary.resolve()
    root = args.output_root.resolve()
    control_rule = json.loads(control_rule_path.read_text(encoding="utf-8"))
    parent_rule = json.loads(parent_rule_path.read_text(encoding="utf-8"))
    locked = control_rule["locked_inputs"]
    contract = control_rule["execution_contract"]
    if sha256(state) != locked["particle_state_sha256"]:
        raise RunError("particle-state hash differs")
    if sha256(parent_rule_path) != control_rule["basis"][
            "common_state_rule_sha256"]:
        raise RunError("parent common-state rule hash differs")
    if sha256(binary) != args.binary_sha256:
        raise RunError("control binary hash differs from pre-execution lock")
    root.mkdir(parents=True, exist_ok=True)
    members = []
    horizons = control_rule["sampling_contract"][
        "matching_aurorapic_post_step_horizons"]
    for horizon in horizons:
        member = root / f"horizon-{horizon:04d}"
        output = member / "output"
        if member.exists():
            members.append(summarize_member(member, horizon))
            continue
        member.mkdir()
        config = member / "input.cfg"
        config.write_text(
            control_deck(parent_rule, state, horizon, output),
            encoding="utf-8")
        stdout = member / "stdout.txt"
        resources = member / "resources.txt"
        limit_bytes = int(contract["address_space_limit_kib"]) * 1024

        def limits() -> None:
            resource.setrlimit(resource.RLIMIT_AS, (limit_bytes, limit_bytes))
            os.nice(10)

        started = time.monotonic()
        with stdout.open("w", encoding="utf-8") as stream:
            completed = subprocess.run(
                ["/usr/bin/time", "-v", "-o", str(resources),
                 str(binary), "--allow-large-run",
                 "I_UNDERSTAND_THIS_IS_A_LARGE_RUN", str(config)],
                stdout=stream, stderr=subprocess.STDOUT,
                timeout=int(contract[
                    "timeout_seconds_each_aurorapic_horizon"]),
                preexec_fn=limits, check=False)
        if completed.returncode != (1 if horizon == 0 else 0):
            raise RunError(f"horizon {horizon} returned {completed.returncode}")
        summary = summarize_member(member, horizon)
        summary["orchestration_wall_seconds"] = time.monotonic() - started
        members.append(summary)
    return {
        "schema_version": 1,
        "scope": "aurorapic_held_density_common_state_horizons",
        "control_rule_sha256": sha256(control_rule_path),
        "parent_rule_sha256": sha256(parent_rule_path),
        "particle_state_sha256": sha256(state),
        "binary_sha256": sha256(binary),
        "subcycle_charge_deposition": "pre_push_held",
        "members": members,
        "all_resource_gates_passed": all(
            item["peak_resident_set_kib"] <=
            int(contract["maximum_peak_resident_set_kib"])
            for item in members),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control-rule", type=Path, required=True)
    parser.add_argument("--parent-rule", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = execute(args)
    except (RunError, OSError, ValueError, KeyError,
            subprocess.SubprocessError) as error:
        parser.error(str(error))
    args.report.write_text(json.dumps(result, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
