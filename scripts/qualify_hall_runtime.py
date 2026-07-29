#!/usr/bin/env python3
"""Run a bounded Hall campaign slice and project published-case cost."""

from __future__ import annotations

import argparse
import configparser
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


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_BOUNDED_HALL_PROBE"
HARD_INITIAL_UPDATE_LIMIT = 25_000_000
WORKSTATION_GENERATION_ACKNOWLEDGEMENT = (
    "I_UNDERSTAND_THIS_IS_AN_OPT_IN_WORKSTATION_RUN"
)


class QualificationError(RuntimeError):
    pass


def positive_integer(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return result


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise QualificationError(f"cannot hash {path}: {error}") from error


def load_manifest(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    try:
        parser.read_string("[global]\n" + path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, configparser.Error) as error:
        raise QualificationError(f"cannot read {path}: {error}") from error
    for section in (
        "global",
        "reference",
        "campaign.micro",
        "campaign.workstation",
        "campaign.production",
    ):
        if section not in parser:
            raise QualificationError(f"manifest is missing [{section}]")
    return parser


def set_global(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(key)}\s*=).*$")
    if not pattern.search(text):
        raise QualificationError(
            f"generated campaign deck is missing {key!r}"
        )
    return pattern.sub(
        lambda match: f"{match.group(1)} {value}", text, count=1
    )


def atomic_json(path: Path, report: dict[str, object]) -> None:
    if path.exists():
        raise QualificationError(
            f"refusing to overwrite existing report: {path}"
        )
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


def tier_work(
    manifest: configparser.ConfigParser,
    tier_name: str,
    steps: int | None = None,
) -> dict[str, int]:
    tier = manifest[f"campaign.{tier_name}"]
    cells_x = tier.getint("cells_x")
    cells_y = tier.getint("cells_y")
    particles_per_cell = tier.getint("particles_per_cell_per_species")
    selected_steps = tier.getint("steps") if steps is None else steps
    initial_particles = 2 * cells_x * cells_y * particles_per_cell
    return {
        "cells_x": cells_x,
        "cells_y": cells_y,
        "particles_per_cell_per_species": particles_per_cell,
        "initial_particles": initial_particles,
        "steps": selected_steps,
        "initial_particle_updates": initial_particles * selected_steps,
    }


def prepare_probe(
    args: argparse.Namespace,
    work: Path,
    output: Path,
) -> Path:
    preparer = Path(__file__).resolve().with_name(
        "prepare_hall_campaign.py"
    )
    if not preparer.is_file():
        raise QualificationError(f"missing campaign preparer: {preparer}")
    generated = work / "generated.cfg"
    command = [
        sys.executable,
        str(preparer),
        str(args.case_manifest),
        "--tier",
        args.tier,
        "--output",
        str(generated),
        "--output-dir",
        str(output),
        "--seed",
        str(args.seed),
    ]
    if args.tier == "workstation":
        command.extend([
            "--acknowledge-cost",
            WORKSTATION_GENERATION_ACKNOWLEDGEMENT,
        ])
    result = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
    )
    if result.returncode != 0:
        raise QualificationError(
            "cannot generate qualification deck: " + result.stderr.strip()
        )
    text = generated.read_text(encoding="utf-8")
    replacements = {
        "steps": str(args.steps),
        "output_interval": str(args.steps),
        "output_dir": str(output),
        "resolved_diagnostics": "false",
        "resolved_diagnostic_interval": str(args.steps),
        "resolved_diagnostic_start_step": "0",
        "checkpoint_output": "false",
        "checkpoint_interval": "0",
        "vtk_output": "false",
        "particle_output": "false",
        "runtime_backend": "serial",
        "runtime_threads": "1",
    }
    for key, value in replacements.items():
        text = set_global(text, key, value)
    probe = work / "probe.cfg"
    probe.write_text(text, encoding="utf-8")
    return probe


def parse_final_population(stdout: str) -> int:
    matches = re.findall(r"\blive_particles=(\d+)\b", stdout)
    if not matches:
        raise QualificationError(
            "AuroraPIC output did not report a final live population"
        )
    return int(matches[-1])


def qualify(args: argparse.Namespace) -> dict[str, object]:
    executable = args.executable.resolve()
    case_manifest = args.case_manifest.resolve()
    if not executable.is_file():
        raise QualificationError(f"missing AuroraPIC CLI: {executable}")
    if args.tier == "workstation" and (
        args.acknowledge_cost != ACKNOWLEDGEMENT
    ):
        raise QualificationError(
            "workstation qualification requires --acknowledge-cost "
            + ACKNOWLEDGEMENT
        )
    manifest = load_manifest(case_manifest)
    probe_work = tier_work(manifest, args.tier, args.steps)
    tier_steps = manifest[f"campaign.{args.tier}"].getint("steps")
    if args.steps > tier_steps:
        raise QualificationError(
            f"probe steps exceed the {args.tier} tier: "
            f"{args.steps} > {tier_steps}"
        )
    if args.max_initial_updates > HARD_INITIAL_UPDATE_LIMIT:
        raise QualificationError(
            "--max-initial-updates exceeds the built-in conservative limit: "
            f"{args.max_initial_updates} > {HARD_INITIAL_UPDATE_LIMIT}"
        )
    if probe_work["initial_particle_updates"] > args.max_initial_updates:
        raise QualificationError(
            "probe exceeds --max-initial-updates: "
            f"{probe_work['initial_particle_updates']} > "
            f"{args.max_initial_updates}"
        )
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    environment["OMP_DYNAMIC"] = "FALSE"
    environment["AURORA_OPENMP_THREADS"] = "1"
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_hall_qualification_"
    ) as temporary:
        work = Path(temporary)
        output = work / "output"
        probe = prepare_probe(args, work, output)
        probe_digest = sha256(probe)
        start = time.perf_counter()
        try:
            result = subprocess.run(
                [str(executable), str(probe)],
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=args.timeout_seconds,
            )
        except subprocess.TimeoutExpired as error:
            raise QualificationError(
                f"bounded Hall probe exceeded {args.timeout_seconds} seconds"
            ) from error
        wall_seconds = time.perf_counter() - start
        if result.returncode != 0:
            raise QualificationError(
                "bounded Hall probe failed: "
                + result.stdout.strip() + " " + result.stderr.strip()
            )
        final_live_particles = parse_final_population(result.stdout)

    update_rate = (
        probe_work["initial_particle_updates"] / wall_seconds
    )
    if not math.isfinite(update_rate) or update_rate <= 0:
        raise QualificationError("invalid measured particle-update rate")
    projections: dict[str, object] = {}
    for name in ("workstation", "production"):
        work = tier_work(manifest, name)
        seconds = work["initial_particle_updates"] / update_rate
        projections[name] = {
            **work,
            "initial_population_only_wall_seconds": seconds,
            "initial_population_only_wall_hours": seconds / 3600.0,
            "initial_population_only_wall_days": seconds / 86400.0,
        }
    affinity = (
        sorted(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity") else None
    )
    return {
        "schema_version": 1,
        "case_id": manifest["global"]["case_id"],
        "qualification_scope": "bounded_runtime_probe",
        "physics_claim": "none",
        "production_launch_authorized": False,
        "probe": {
            "tier": args.tier,
            "seed": args.seed,
            **probe_work,
            "final_live_particles": final_live_particles,
            "wall_seconds": wall_seconds,
            "initial_particle_updates_per_second": update_rate,
            "runtime_backend": "serial",
            "runtime_threads": 1,
            "max_initial_updates": args.max_initial_updates,
            "timeout_seconds": args.timeout_seconds,
        },
        "projections": projections,
        "provenance": {
            "case_manifest": str(case_manifest),
            "case_manifest_sha256": sha256(case_manifest),
            "executable": str(executable),
            "executable_sha256": sha256(executable),
            "probe_config_sha256": probe_digest,
            "platform": platform.platform(),
            "python": platform.python_version(),
            "available_cpu_affinity": affinity,
        },
        "warnings": [
            "The probe is deliberately too short and/or coarse for a "
            "physics-agreement claim.",
            "Projection counts only the initial population; pair and cathode "
            "sources increase live particles during a full run.",
            "Measured wall time includes startup and output overhead and is "
            "specific to this binary, host, affinity, and storage.",
            "A production campaign still requires external reference data, "
            "convergence runs, multiple seeds, scheduler limits, and explicit "
            "execution authorization.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument(
        "case_manifest",
        type=Path,
        nargs="?",
        default=Path("examples/hall_landmark_axial_azimuthal.case"),
    )
    parser.add_argument(
        "--tier", choices=("micro", "workstation"), default="micro"
    )
    parser.add_argument("--steps", type=positive_integer, default=20)
    parser.add_argument(
        "--max-initial-updates", type=positive_integer, default=25_000_000
    )
    parser.add_argument(
        "--timeout-seconds", type=positive_integer, default=60
    )
    parser.add_argument("--seed", type=int, default=24680)
    parser.add_argument("--acknowledge-cost")
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.seed < 0 or args.seed > 0xFFFFFFFF:
        parser.error("--seed must be an unsigned 32-bit integer")
    return args


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
        print(f"Hall runtime qualification error: {error}", file=sys.stderr)
        return 2
    print(
        "Bounded Hall runtime qualification passed without authorizing "
        f"production: {args.report}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
