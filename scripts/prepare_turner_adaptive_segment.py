#!/usr/bin/env python3
"""Prepare, but never launch, a locked Turner adaptive segment."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_PREPARES_A_LARGE_LOCKED_RUN"


class PreparationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PreparationError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, label: str) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise PreparationError(f"cannot read {label}: {error}") from error
    require(isinstance(value, dict), f"{label} is not an object")
    return value


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(value, encoding="utf-8")
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def global_value(text: str, key: str) -> str:
    matches = re.findall(
        rf"(?m)^\s*{re.escape(key)}\s*=\s*(.*?)\s*$", text)
    require(len(matches) == 1, f"base deck must contain one {key!r}")
    return matches[0]


def set_or_insert_global(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(key)}\s*=\s*).*$")
    result, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    require(count <= 1, f"base deck contains duplicate {key!r}")
    if count:
        return result
    marker = text.find("\n[")
    require(marker >= 0, "base deck has no configuration sections")
    return text[:marker] + f"\n{key} = {value}" + text[marker:]


def checkpoint_step(path: Path) -> tuple[str, int]:
    with path.open("r", encoding="utf-8") as stream:
        magic = stream.readline().rstrip("\n")
        for line in stream:
            match = re.fullmatch(r"step\s+(\d+)\s*", line)
            if match:
                return magic, int(match.group(1))
    raise PreparationError("checkpoint has no step record")


def prepare(args: argparse.Namespace) -> dict:
    require(args.acknowledge_cost == ACKNOWLEDGEMENT,
            f"preparation requires --acknowledge-cost {ACKNOWLEDGEMENT}")
    rule_path = args.rule.resolve()
    lock_path = args.execution_lock.resolve()
    base = args.base_config.resolve()
    checkpoint = args.checkpoint.resolve()
    solver = args.solver.resolve()
    deck = args.output_config.resolve()
    report_path = args.report.resolve()
    output_dir = args.output_dir.resolve()
    prior_progress_arg = getattr(args, "prior_progress", None)
    prior_progress_path = (prior_progress_arg.resolve()
                           if prior_progress_arg is not None else None)
    for path, label in ((rule_path, "rule"), (lock_path, "execution lock"),
                        (base, "base config"), (checkpoint, "checkpoint"),
                        (solver, "solver")):
        require(path.is_file(), f"{label} does not exist: {path}")
    require(not deck.exists(), f"refusing to overwrite deck: {deck}")
    require(not report_path.exists(),
            f"refusing to overwrite report: {report_path}")
    require(not output_dir.exists(),
            f"refusing existing output directory: {output_dir}")

    rule = load_json(rule_path, "rule")
    lock = load_json(lock_path, "execution lock")
    require(lock.get("status") == "preregistered_not_launched",
            "execution lock is not in preregistered state")
    require(lock.get("rule", {}).get("sha256") == sha256(rule_path),
            "execution lock does not match the rule")
    require(lock.get("command_identity", {}).get("solver_binary_sha256") ==
            sha256(solver), "solver binary does not match execution lock")
    seed_order = lock.get("execution_order", {}).get("seeds", [])
    require(args.seed in seed_order, "seed is not in the locked execution order")
    initial = next((value for value in rule.get("locked_initial_states", [])
                    if value.get("seed") == args.seed), None)
    require(isinstance(initial, dict), "seed has no locked initial state")
    magic, source_step = checkpoint_step(checkpoint)

    rf = rule["rf_contract"]
    adaptive = rule["adaptive_equilibration"]
    steps_per_cycle = int(rf["steps_per_cycle"])
    cycles_per_block = int(rf["cycles_per_block"])
    segment_blocks = int(adaptive["execution_segment_blocks"])
    maximum_blocks = int(adaptive["maximum_blocks_per_seed"])
    initial_step = int(initial["checkpoint_step"])
    segment_index = 1
    completed_blocks = 0
    imports_prior_samples = False
    progress = None
    admitted = []
    if prior_progress_path is not None:
        require(prior_progress_path.is_file(),
                f"prior progress does not exist: {prior_progress_path}")
        progress = load_json(prior_progress_path, "prior progress")
        require(progress.get("case_id") == rule.get("case_id"),
                "prior progress case differs from the rule")
        require(progress.get("rule", {}).get("sha256") == sha256(rule_path),
                "prior progress does not match the rule")
        require(progress.get("execution_lock", {}).get("sha256") ==
                sha256(lock_path),
                "prior progress does not match the execution lock")
        all_admitted = progress.get("admitted_segments", [])
        require(isinstance(all_admitted, list),
                "prior progress admitted_segments is not a list")
        for earlier_seed in seed_order[:seed_order.index(args.seed)]:
            earlier = sorted(
                (value for value in all_admitted
                 if value.get("seed") == earlier_seed),
                key=lambda value: value.get("segment", 0))
            require(earlier and earlier[-1].get("classification") ==
                    "admitted_converged",
                    f"earlier seed {earlier_seed} has not converged")
        admitted = sorted(
            (value for value in all_admitted
             if value.get("seed") == args.seed),
            key=lambda value: value.get("segment", 0))
    require(seed_order.index(args.seed) == 0 or progress is not None,
            "later seeds require prior progress proving execution order")

    if not admitted:
        require(initial.get("base_config_sha256") == sha256(base),
                "base config hash differs from the locked state")
        require(initial.get("checkpoint_sha256") == sha256(checkpoint),
                "checkpoint hash differs from the locked state")
        require(source_step == initial_step,
                "checkpoint step differs from the locked state")
        require(magic.startswith("AuroraPIC-checkpoint-v") and magic !=
                "AuroraPIC-checkpoint-v25",
                "first segment requires the locked pre-v25 checkpoint")
    else:
        require([value.get("segment") for value in admitted] ==
                list(range(1, len(admitted) + 1)),
                "prior admitted segment sequence is not contiguous")
        require(all(value.get("integrity", {}).get("admitted") is True
                    for value in admitted),
                "prior progress contains an unadmitted segment")
        previous = admitted[-1]
        require(previous.get("classification") ==
                "admitted_horizon_incomplete",
                "prior segment does not authorize continuation")
        completed_blocks = sum(int(value.get("blocks", 0))
                               for value in admitted)
        require(0 < completed_blocks < maximum_blocks,
                "prior block count cannot be continued")
        require(previous.get("generated_deck_sha256") == sha256(base),
                "base config does not match the prior admitted deck")
        require(previous.get("output_checkpoint_sha256") ==
                sha256(checkpoint),
                "checkpoint does not match the prior admitted output")
        require(source_step == previous.get("end_step"),
                "checkpoint step differs from prior admitted endpoint")
        require(source_step == initial_step + completed_blocks *
                cycles_per_block * steps_per_cycle,
                "prior endpoint is inconsistent with its block count")
        require(magic == "AuroraPIC-checkpoint-v25",
                "continuation requires a v25 convergence checkpoint")
        segment_index = int(previous["segment"]) + 1
        imports_prior_samples = True

    require(source_step % steps_per_cycle == 0,
            "checkpoint is not at the locked RF phase")
    blocks_this_segment = min(segment_blocks, maximum_blocks - completed_blocks)
    target_step = (source_step + blocks_this_segment * cycles_per_block *
                   steps_per_cycle)
    maximum_step = (initial_step + maximum_blocks * cycles_per_block *
                    steps_per_cycle)

    text = base.read_text(encoding="utf-8")
    require(global_value(text, "config_version") == "1" and
            global_value(text, "units") == "si" and
            global_value(text, "dimension") == "1" and
            global_value(text, "runtime_backend") == "serial" and
            global_value(text, "runtime_threads") == "1",
            "base config differs from the locked serial 1D SI contract")
    require(re.findall(r"(?m)^\s*opportunity_sampling\s*=\s*(\S+)\s*$",
                       text) == ["single_bernoulli", "single_bernoulli"],
            "base config collision scheduling differs")
    values = {
        "steps": str(target_step),
        "mode": "steady_state",
        "max_steps": str(target_step),
        "output_interval": str(cycles_per_block * steps_per_cycle),
        "output_dir": str(output_dir),
        "restart_path": str(checkpoint),
        "checkpoint_output": "true",
        "checkpoint_interval": str(cycles_per_block * steps_per_cycle),
        "spatial_average": "true",
        "spatial_average_reset_on_restart": "true",
        "spatial_average_interval": "1",
        "spatial_average_start_step": str(target_step + 1),
        "spatial_average_end_step": str(
            target_step + cycles_per_block * steps_per_cycle),
        "spatial_average_rf_frequency": str(int(rf["frequency_hz"])),
        "spatial_average_rf_cycles": str(cycles_per_block),
        "periodic_convergence": "true",
        "periodic_convergence_reset_on_restart": (
            "false" if imports_prior_samples else "true"),
        "periodic_convergence_rf_frequency": str(int(rf["frequency_hz"])),
        "periodic_convergence_cycles_per_block": str(cycles_per_block),
        "periodic_convergence_minimum_blocks": str(
            adaptive["minimum_nominal_blocks"]),
        "periodic_convergence_minimum_effective_blocks": str(
            adaptive["minimum_ar1_effective_blocks_per_observable"]),
        "periodic_convergence_maximum_absolute_projected_fractional_drift": str(
            adaptive["maximum_absolute_projected_fractional_drift_per_observable"]),
        "periodic_convergence_maximum_absolute_split_half_fractional_change": str(
            adaptive["maximum_absolute_split_half_fractional_change_per_observable"]),
        "periodic_convergence_maximum_relative_standard_error": str(
            adaptive["maximum_relative_standard_error_per_observable"]),
    }
    for key, value in values.items():
        text = set_or_insert_global(text, key, value)
    atomic_text(deck, text)
    report = {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "locked_adaptive_convergence_segment_preflight",
        "rule": {"path": str(rule_path), "sha256": sha256(rule_path)},
        "execution_lock": {"path": str(lock_path),
                           "sha256": sha256(lock_path)},
        "seed": args.seed,
        "source": {"base_config": str(base),
                   "base_config_sha256": sha256(base),
                   "checkpoint": str(checkpoint),
                   "checkpoint_sha256": sha256(checkpoint),
                   "checkpoint_magic": magic, "source_step": source_step},
        "segment": {"index": segment_index, "blocks": blocks_this_segment,
                    "prior_admitted_blocks": completed_blocks,
                    "target_cumulative_blocks": (
                        completed_blocks + blocks_this_segment),
                    "start_step": source_step + 1,
                    "target_step": target_step,
                    "maximum_campaign_step": maximum_step},
        "generated_deck": str(deck),
        "generated_deck_sha256": sha256(deck),
        "output_dir": str(output_dir),
        "prior_progress": ({"path": str(prior_progress_path),
                            "sha256": sha256(prior_progress_path)}
                           if prior_progress_path is not None else None),
        "periodic_convergence_epoch_imports_prior_samples":
            imports_prior_samples,
        "launched": False,
        "physics_claim": "none_preparation_only"
    }
    atomic_text(report_path,
                json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--execution-lock", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--prior-progress", type=Path)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-config", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    try:
        report = prepare(parse_args())
    except (PreparationError, OSError, UnicodeError, ValueError) as error:
        print(f"Turner adaptive preparation error: {error}", file=sys.stderr)
        return 2
    print(json.dumps({"deck": report["generated_deck"],
                      "deck_sha256": report["generated_deck_sha256"],
                      "target_step": report["segment"]["target_step"],
                      "launched": False}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
