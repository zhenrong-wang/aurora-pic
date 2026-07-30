#!/usr/bin/env python3
"""Run a strictly bounded Turner Case 1 slice and project campaign cost."""

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
import tempfile
import time


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_BOUNDED_TURNER_PROBE"
GENERATION_ACKNOWLEDGEMENT = (
    "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_RUN"
)
HARD_INITIAL_UPDATE_LIMIT = 2_000_000
PRODUCTION_STEPS = 512_000
PRODUCTION_PARTICLES = 131_072


class QualificationError(RuntimeError):
    pass


def positive_integer(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as error:
        raise QualificationError(f"cannot hash {path}: {error}") from error
    return digest.hexdigest()


def set_global(text: str, key: str, value: str) -> str:
    lines = text.splitlines()
    found = False
    for index, line in enumerate(lines):
        if line.lstrip().startswith("["):
            break
        match = re.match(rf"^(\s*{re.escape(key)}\s*=).*$", line)
        if match:
            if found:
                raise QualificationError(f"duplicate global key {key!r}")
            lines[index] = f"{match.group(1)} {value}"
            found = True
    if not found:
        raise QualificationError(f"generated campaign deck is missing {key!r}")
    return "\n".join(lines) + "\n"


def atomic_json(path: Path, report: dict) -> None:
    if path.exists():
        raise QualificationError(f"refusing to overwrite existing report: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def prepare_probe(args: argparse.Namespace, work: Path, output_dir: Path) -> tuple:
    preparer = Path(__file__).resolve().with_name("prepare_turner_case.py")
    if not preparer.is_file():
        raise QualificationError(f"missing campaign preparer: {preparer}")
    production = work / "production.cfg"
    preflight = production.with_suffix(".preflight.json")
    command = [
        sys.executable,
        str(preparer),
        str(args.case_manifest.resolve()),
        str(args.normalized_dir.resolve()),
        "--output",
        str(production),
        "--acknowledge-cost",
        GENERATION_ACKNOWLEDGEMENT,
    ]
    generated = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
    )
    if generated.returncode != 0:
        raise QualificationError(
            "cannot generate exact campaign deck: " + generated.stderr.strip()
        )
    try:
        production_report = json.loads(preflight.read_text(encoding="utf-8"))
        text = production.read_text(encoding="utf-8")
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise QualificationError(f"cannot inspect generated campaign: {error}") from error
    contract = production_report.get("contract", {})
    require_exact = (
        contract.get("steps") == PRODUCTION_STEPS
        and 2 * contract.get("particles_per_species", 0) == PRODUCTION_PARTICLES
        and contract.get("averaging_samples") == 12_800
        and not production_report.get("full_run_launched", True)
    )
    if not require_exact:
        raise QualificationError("generated deck is not the exact Case 1 contract")
    replacements = {
        "steps": str(args.steps),
        "output_interval": str(args.steps),
        "output_dir": str(output_dir),
        "checkpoint_output": "false",
        "checkpoint_interval": "0",
        "spatial_average": "false",
        "spatial_average_interval": "1",
        "spatial_average_start_step": "0",
        "spatial_average_end_step": "0",
        "spatial_average_rf_frequency": "0",
        "spatial_average_rf_cycles": "0",
        "runtime_backend": "serial",
        "runtime_threads": "1",
    }
    for key, value in replacements.items():
        text = set_global(text, key, value)
    probe = work / "probe.cfg"
    probe.write_text(text, encoding="utf-8")
    return probe, production, production_report


def conservative_child_setup() -> None:
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


def final_live_particles(output_dir: Path) -> int:
    path = output_dir / "scalars.csv"
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
    except (OSError, UnicodeError, csv.Error) as error:
        raise QualificationError(f"cannot read probe diagnostics: {error}") from error
    if not rows or "live_particles" not in rows[-1]:
        raise QualificationError("probe diagnostics lack a final live population")
    try:
        return int(rows[-1]["live_particles"])
    except ValueError as error:
        raise QualificationError("invalid final live population") from error


def qualify(args: argparse.Namespace) -> dict:
    executable = args.executable.resolve()
    case_manifest = args.case_manifest.resolve()
    normalized_dir = args.normalized_dir.resolve()
    if not executable.is_file():
        raise QualificationError(f"missing AuroraPIC CLI: {executable}")
    if args.acknowledge_cost != ACKNOWLEDGEMENT:
        raise QualificationError(
            "runtime qualification requires --acknowledge-cost "
            + ACKNOWLEDGEMENT
        )
    initial_updates = PRODUCTION_PARTICLES * args.steps
    if args.max_initial_updates > HARD_INITIAL_UPDATE_LIMIT:
        raise QualificationError(
            "--max-initial-updates exceeds the built-in conservative limit: "
            f"{args.max_initial_updates} > {HARD_INITIAL_UPDATE_LIMIT}"
        )
    if initial_updates > args.max_initial_updates:
        raise QualificationError(
            "probe exceeds --max-initial-updates: "
            f"{initial_updates} > {args.max_initial_updates}"
        )

    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    environment["OMP_DYNAMIC"] = "FALSE"
    environment["AURORA_OPENMP_THREADS"] = "1"
    parent_affinity = (
        sorted(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity") else None
    )
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_qualification_"
    ) as temporary:
        work = Path(temporary)
        output_dir = work / "output"
        probe, production, preflight = prepare_probe(args, work, output_dir)
        probe_digest = sha256(probe)
        production_digest = sha256(production)
        start = time.perf_counter()
        try:
            result = subprocess.run(
                [str(executable), str(probe)],
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=args.timeout_seconds,
                preexec_fn=conservative_child_setup if os.name == "posix" else None,
            )
        except subprocess.TimeoutExpired as error:
            raise QualificationError(
                f"bounded Turner probe exceeded {args.timeout_seconds} seconds"
            ) from error
        wall_seconds = time.perf_counter() - start
        if result.returncode != 0:
            raise QualificationError(
                "bounded Turner probe failed: "
                + result.stdout.strip() + " " + result.stderr.strip()
            )
        live_particles = final_live_particles(output_dir)

    update_rate = initial_updates / wall_seconds
    if not math.isfinite(update_rate) or update_rate <= 0:
        raise QualificationError("invalid measured particle-update rate")
    production_updates = PRODUCTION_PARTICLES * PRODUCTION_STEPS
    projected_seconds = production_updates / update_rate
    return {
        "turner_runtime_qualification_version": 1,
        "case_id": preflight["case_id"],
        "qualification_scope": "bounded_exact_population_runtime_probe",
        "physics_claim": "none",
        "production_launch_authorized": False,
        "probe": {
            "steps": args.steps,
            "initial_particles": PRODUCTION_PARTICLES,
            "initial_particle_updates": initial_updates,
            "final_live_particles": live_particles,
            "wall_seconds": wall_seconds,
            "initial_particle_updates_per_second": update_rate,
            "runtime_backend": "serial",
            "runtime_threads": 1,
            "nice_increment_requested": 10,
            "single_cpu_affinity_requested": os.name == "posix",
            "max_initial_updates": args.max_initial_updates,
            "timeout_seconds": args.timeout_seconds,
        },
        "projection": {
            "production_steps": PRODUCTION_STEPS,
            "initial_particles": PRODUCTION_PARTICLES,
            "initial_population_particle_updates": production_updates,
            "initial_population_only_wall_seconds": projected_seconds,
            "initial_population_only_wall_hours": projected_seconds / 3600.0,
            "initial_population_only_wall_days": projected_seconds / 86400.0,
        },
        "provenance": {
            "case_manifest": str(case_manifest),
            "case_manifest_sha256": sha256(case_manifest),
            "normalization_audit": str(normalized_dir / "audit.json"),
            "normalization_audit_sha256": sha256(normalized_dir / "audit.json"),
            "executable": str(executable),
            "executable_sha256": sha256(executable),
            "production_config_sha256": production_digest,
            "probe_config_sha256": probe_digest,
            "platform": platform.platform(),
            "python": platform.python_version(),
            "parent_cpu_affinity": parent_affinity,
        },
        "warnings": [
            "This shortened run cannot establish stationarity or physics agreement.",
            "The projection counts only the initial population; ionization and "
            "electrode losses change work during the full campaign.",
            "Startup and diagnostic overhead are included, while long-run cache, "
            "storage, and checkpoint behavior are not sampled.",
            "The measurement is specific to this binary, host, affinity, and data.",
            "A full launch still requires separate explicit authorization.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("normalized_dir", type=Path)
    parser.add_argument("--steps", type=positive_integer, default=4)
    parser.add_argument(
        "--max-initial-updates", type=positive_integer, default=1_000_000
    )
    parser.add_argument(
        "--timeout-seconds", type=positive_integer, default=60
    )
    parser.add_argument("--acknowledge-cost")
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.report.exists():
            raise QualificationError(
                f"refusing to overwrite existing report: {args.report}"
            )
        report = qualify(args)
        atomic_json(args.report, report)
    except (QualificationError, OSError, UnicodeError) as error:
        print(f"Turner runtime qualification error: {error}", file=sys.stderr)
        return 2
    print(
        "Bounded Turner runtime qualification passed without authorizing "
        f"production: {args.report}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
