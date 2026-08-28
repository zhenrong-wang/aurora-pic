#!/usr/bin/env python3
"""Prepare one checksum-locked corrected Turner Case 1 density window."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys


STEPS_PER_CYCLE = 400
WINDOW_CYCLES = 32
WINDOW_STEPS = STEPS_PER_CYCLE * WINDOW_CYCLES
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_POST_BENCHMARK_DIAGNOSTIC"


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


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(value, encoding="utf-8")
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def set_or_insert_global(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(key)}\s*=\s*).*$")
    result, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    require(count <= 1, f"base deck contains duplicate {key!r}")
    if count == 1:
        return result
    marker = text.find("\n[")
    require(marker >= 0, "base deck has no configuration sections")
    return text[:marker] + f"\n{key} = {value}" + text[marker:]


def global_value(text: str, key: str) -> str:
    matches = re.findall(
        rf"(?m)^\s*{re.escape(key)}\s*=\s*(.*?)\s*$", text
    )
    require(len(matches) == 1, f"base deck must contain one {key!r}")
    return matches[0]


def build_deck(base: str, checkpoint: Path, output_dir: Path,
               source_step: int) -> str:
    require(global_value(base, "config_version") == "1",
            "base deck config_version must be 1")
    require(global_value(base, "units") == "si",
            "base deck must use SI units")
    require(global_value(base, "dimension") == "1",
            "base deck must be one-dimensional")
    require(global_value(base, "nx") == "129",
            "base deck must use Turner Case 1's 129 nodes")
    require(global_value(base, "runtime_backend") == "serial",
            "base deck runtime backend must be serial")
    require(global_value(base, "runtime_threads") == "1",
            "base deck runtime thread count must be one")
    opportunity_modes = re.findall(
        r"(?m)^\s*opportunity_sampling\s*=\s*(\S+)\s*$", base
    )
    require(
        opportunity_modes == ["single_bernoulli", "single_bernoulli"],
        "base deck must use single_bernoulli for both collision models",
    )
    require(
        Path(global_value(base, "output_dir")).resolve() != output_dir,
        "continuation output directory must differ from the base output",
    )

    end_step = source_step + WINDOW_STEPS
    values = {
        # A converged steady-state checkpoint may carry an already-satisfied
        # periodic controller.  Profile windows are fixed-duration
        # measurements, so they must not terminate from restored controller
        # state before collecting any spatial samples.
        "mode": "transient",
        "steps": str(end_step),
        "output_interval": str(STEPS_PER_CYCLE),
        "output_dir": str(output_dir),
        "restart_path": str(checkpoint),
        "checkpoint_output": "true",
        "checkpoint_interval": str(WINDOW_STEPS),
        "spatial_average": "true",
        "spatial_average_reset_on_restart": "true",
        "spatial_average_interval": "1",
        "spatial_average_start_step": str(source_step + 1),
        "spatial_average_end_step": str(end_step),
        "spatial_average_rf_frequency": "13560000",
        "spatial_average_rf_cycles": str(WINDOW_CYCLES),
        "periodic_convergence": "false",
        "periodic_convergence_reset_on_restart": "false",
        "periodic_convergence_rf_frequency": "0",
        "periodic_convergence_cycles_per_block": "0",
        "runtime_backend": "serial",
        "runtime_threads": "1",
    }
    result = base
    for key, value in values.items():
        result = set_or_insert_global(result, key, value)
    return result


def prepare(args: argparse.Namespace) -> dict[str, object]:
    require(args.acknowledge_cost == ACKNOWLEDGEMENT,
            f"preparation requires --acknowledge-cost {ACKNOWLEDGEMENT}")
    base = args.base_config.resolve()
    checkpoint = args.checkpoint.resolve()
    output_dir = args.output_dir.resolve()
    deck = args.output_config.resolve()
    report = args.report.resolve()
    require(base.is_file(), f"base config does not exist: {base}")
    require(checkpoint.is_file(), f"checkpoint does not exist: {checkpoint}")
    require(checkpoint.stat().st_size > 0, "checkpoint is empty")
    require(args.source_step >= 0, "source step must be non-negative")
    require(
        checkpoint.name == f"checkpoint_{args.source_step}.apc",
        "checkpoint filename does not match --source-step",
    )
    require(not output_dir.exists(),
            f"refusing existing output directory: {output_dir}")
    require(not deck.exists(), f"refusing to overwrite deck: {deck}")
    require(not report.exists(), f"refusing to overwrite report: {report}")
    base_text = base.read_text(encoding="utf-8")
    deck_text = build_deck(base_text, checkpoint, output_dir, args.source_step)
    atomic_text(deck, deck_text)
    end_step = args.source_step + WINDOW_STEPS
    value = {
        "schema_version": 1,
        "case_id": "turner-helium-ccp-2013-case-1",
        "scope": "corrected_schedule_post_benchmark_density_window_preflight",
        "source": {
            "base_config": str(base),
            "base_config_sha256": sha256(base),
            "checkpoint": str(checkpoint),
            "checkpoint_sha256": sha256(checkpoint),
            "source_step": args.source_step,
        },
        "window": {
            "start_step": args.source_step + 1,
            "end_step": end_step,
            "steps": WINDOW_STEPS,
            "rf_cycles": WINDOW_CYCLES,
            "averaging_samples": WINDOW_STEPS,
            "reset_on_restart": True,
        },
        "execution": {
            "mode": "transient",
            "runtime_backend": "serial",
            "runtime_threads": 1,
            "periodic_convergence": False,
            "collision_opportunity_sampling": "single_bernoulli",
            "output_dir": str(output_dir),
            "launched": False,
        },
        "generated_deck": str(deck),
        "generated_deck_sha256": sha256(deck),
        "physics_claim": "none_preparation_only",
    }
    atomic_text(report, json.dumps(value, indent=2, sort_keys=True) + "\n")
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base_config", type=Path)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--source-step", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-config", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    try:
        value = prepare(parse_args())
    except (PreparationError, OSError, UnicodeError, ValueError) as error:
        print(f"Turner density-window preparation error: {error}", file=sys.stderr)
        return 2
    print(json.dumps({
        "deck": value["generated_deck"],
        "deck_sha256": value["generated_deck_sha256"],
        "end_step": value["window"]["end_step"],
        "launched": False,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
