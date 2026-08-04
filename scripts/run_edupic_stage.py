#!/usr/bin/env python3
"""Run one conservative, checkpoint-preserving eduPIC equilibration stage."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import time


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_BOUNDED_EDUPIC_STAGE"
TIMESTEPS_PER_CYCLE = 4000
HARD_CYCLE_LIMIT = 16
HARD_INITIAL_PARTICLE_STEP_LIMIT = 1_000_000_000
HARD_TIMEOUT_SECONDS = 120
MAX_PARTICLES_PER_SPECIES = 1_000_000


class StageError(RuntimeError):
    pass


def positive_integer(text: str) -> int:
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checkpoint_state(path: Path) -> dict:
    """Read the small header/count contract without loading particle arrays."""
    try:
        size = path.stat().st_size
        with path.open("rb") as stream:
            def scalar() -> float:
                data = stream.read(8)
                if len(data) != 8:
                    raise StageError(f"truncated eduPIC checkpoint: {path}")
                return struct.unpack("=d", data)[0]

            simulation_time = scalar()
            cycles_value = scalar()
            electron_value = scalar()
            cycles = int(cycles_value)
            electrons = int(electron_value)
            if cycles_value != cycles or electron_value != electrons:
                raise StageError(f"non-integral eduPIC checkpoint header: {path}")
            if not 0 <= electrons <= MAX_PARTICLES_PER_SPECIES:
                raise StageError(f"invalid electron population in {path}")
            stream.seek(4 * electrons * 8, os.SEEK_CUR)
            ion_value = scalar()
            ions = int(ion_value)
            if ion_value != ions or not 0 <= ions <= MAX_PARTICLES_PER_SPECIES:
                raise StageError(f"invalid ion population in {path}")
    except OSError as error:
        raise StageError(f"cannot read eduPIC checkpoint {path}: {error}") from error
    expected_size = 32 + 32 * (electrons + ions)
    if size != expected_size:
        raise StageError(
            f"eduPIC checkpoint size mismatch: {size} != {expected_size}")
    if cycles < 0 or not math.isfinite(simulation_time) or simulation_time < 0.0:
        raise StageError(f"invalid time/cycle header in {path}")
    return {"cycles": cycles, "electrons": electrons, "ions": ions,
            "total_particles": electrons + ions,
            "simulation_time_seconds": simulation_time, "bytes": size,
            "sha256": sha256(path)}


def convergence_rows(path: Path) -> list[dict[str, int]]:
    rows = []
    try:
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 3:
                raise StageError(f"invalid conv.dat line {line_number}")
            cycle, electrons, ions = map(int, fields)
            if cycle <= 0 or electrons < 0 or ions < 0:
                raise StageError(f"invalid conv.dat values on line {line_number}")
            if rows and cycle != rows[-1]["cycle"] + 1:
                raise StageError("conv.dat cycles are not contiguous")
            rows.append({"cycle": cycle, "electrons": electrons, "ions": ions,
                         "total_particles": electrons + ions})
    except (OSError, UnicodeError, ValueError) as error:
        if isinstance(error, StageError):
            raise
        raise StageError(f"cannot parse {path}: {error}") from error
    if not rows:
        raise StageError(f"empty convergence history: {path}")
    return rows


def child_setup() -> None:
    try:
        os.nice(10)
    except OSError:
        pass
    if hasattr(os, "sched_getaffinity") and hasattr(os, "sched_setaffinity"):
        try:
            available = os.sched_getaffinity(0)
            if available:
                os.sched_setaffinity(0, {min(available)})
        except OSError:
            pass


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    os.replace(temporary, path)


def run(args: argparse.Namespace) -> dict:
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise StageError("stage requires --acknowledge-cost " + ACKNOWLEDGEMENT)
    if args.cycles > HARD_CYCLE_LIMIT:
        raise StageError(f"cycles exceed built-in limit {HARD_CYCLE_LIMIT}")
    if args.timeout_seconds > HARD_TIMEOUT_SECONDS:
        raise StageError(f"timeout exceeds built-in limit {HARD_TIMEOUT_SECONDS} seconds")
    executable = args.executable.resolve()
    input_dir = args.input_state_dir.resolve()
    stage_dir = args.stage_dir.resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise StageError(f"missing executable eduPIC binary: {executable}")
    executable_sha256 = sha256(executable)
    if executable_sha256 != args.expected_binary_sha256.lower():
        raise StageError("eduPIC binary SHA-256 differs from the locked value")
    if stage_dir.exists():
        raise StageError(f"refusing to overwrite stage directory: {stage_dir}")
    input_checkpoint = input_dir / "picdata.bin"
    input_convergence = input_dir / "conv.dat"
    initial = checkpoint_state(input_checkpoint)
    if initial["sha256"] != args.expected_input_sha256.lower():
        raise StageError("input checkpoint SHA-256 differs from the locked value")
    history = convergence_rows(input_convergence)
    if history[-1]["cycle"] != initial["cycles"]:
        raise StageError("checkpoint and convergence history cycles differ")
    if (history[-1]["electrons"] != initial["electrons"] or
            history[-1]["ions"] != initial["ions"]):
        raise StageError("checkpoint and convergence history populations differ")
    initial_particle_steps = (initial["total_particles"] * args.cycles *
                              TIMESTEPS_PER_CYCLE)
    if args.max_initial_particle_steps > HARD_INITIAL_PARTICLE_STEP_LIMIT:
        raise StageError("requested work limit exceeds built-in hard limit")
    if initial_particle_steps > args.max_initial_particle_steps:
        raise StageError(
            f"stage exceeds initial-particle-step limit: {initial_particle_steps} > "
            f"{args.max_initial_particle_steps}")

    stage_dir.mkdir(parents=True)
    shutil.copy2(input_checkpoint, stage_dir / "picdata.bin")
    shutil.copy2(input_convergence, stage_dir / "conv.dat")
    environment = dict(os.environ)
    environment.update({"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"})
    command = [str(executable), str(args.cycles)]
    start = time.perf_counter()
    try:
        completed = subprocess.run(
            command, cwd=stage_dir, env=environment, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=args.timeout_seconds,
            preexec_fn=child_setup if os.name == "posix" else None)
    except subprocess.TimeoutExpired as error:
        def timeout_text(value: str | bytes | None) -> str:
            if isinstance(value, bytes):
                return value.decode("utf-8", errors="replace")
            return value or ""
        (stage_dir / "stdout.txt").write_text(timeout_text(error.stdout), encoding="utf-8")
        (stage_dir / "stderr.txt").write_text(timeout_text(error.stderr), encoding="utf-8")
        raise StageError(
            f"eduPIC stage exceeded {args.timeout_seconds} second timeout") from error
    wall_seconds = time.perf_counter() - start
    (stage_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (stage_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise StageError(f"eduPIC stage failed with status {completed.returncode}")

    final = checkpoint_state(stage_dir / "picdata.bin")
    final_history = convergence_rows(stage_dir / "conv.dat")
    expected_cycle = initial["cycles"] + args.cycles
    if final["cycles"] != expected_cycle or final_history[-1]["cycle"] != expected_cycle:
        raise StageError("eduPIC stage did not reach its exact requested cycle")
    new_rows = final_history[len(history):]
    if len(new_rows) != args.cycles or new_rows[0]["cycle"] != initial["cycles"] + 1:
        raise StageError("eduPIC stage convergence history coverage is incomplete")
    if (new_rows[-1]["electrons"] != final["electrons"] or
            new_rows[-1]["ions"] != final["ions"]):
        raise StageError("final checkpoint and convergence population differ")
    report = {
        "schema_version": 1, "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "bounded_external_equilibration_stage", "physics_claim": "none",
        "source_binary": {"path": str(executable), "sha256": executable_sha256,
                          "locked_before_execution": True},
        "resource_policy": {"cpu_affinity_count": 1, "nice_increment": 10,
                            "timeout_seconds": args.timeout_seconds,
                            "initial_particle_step_limit": args.max_initial_particle_steps,
                            "initial_particle_steps": initial_particle_steps},
        "stage": {"requested_cycles": args.cycles,
                  "wall_seconds": wall_seconds,
                  "cycles_per_wall_second": args.cycles / wall_seconds,
                  "start_cycle": initial["cycles"], "end_cycle": final["cycles"]},
        "initial_state": initial, "final_state": final,
        "new_cycle_population": new_rows,
        "outputs": {"stdout_sha256": sha256(stage_dir / "stdout.txt"),
                    "stderr_sha256": sha256(stage_dir / "stderr.txt"),
                    "convergence_sha256": sha256(stage_dir / "conv.dat")},
        "completed": True,
    }
    atomic_json(stage_dir / "stage-report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("input_state_dir", type=Path)
    parser.add_argument("stage_dir", type=Path)
    parser.add_argument("--cycles", type=positive_integer, required=True)
    parser.add_argument("--expected-binary-sha256", required=True)
    parser.add_argument("--expected-input-sha256", required=True)
    parser.add_argument("--timeout-seconds", type=positive_integer, default=60)
    parser.add_argument("--max-initial-particle-steps", type=positive_integer,
                        default=50_000_000)
    parser.add_argument("--acknowledge-cost")
    args = parser.parse_args()
    try:
        report = run(args)
    except (StageError, OSError) as error:
        print(f"eduPIC stage rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
