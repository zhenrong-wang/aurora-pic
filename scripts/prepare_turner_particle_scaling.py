#!/usr/bin/env python3
"""Prepare fresh-seed Turner 1x/2x macro-particle scaling arms."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from prepare_turner_case import PreparationError, prepare


ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_LONG_TURNER_PARTICLE_SCALING_CAMPAIGN"
TURNER_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_RUN"
STEPS_PER_RF_CYCLE = 400
MEASUREMENT_RF_CYCLES_PER_BLOCK = 32


class ParticleScalingPreparationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ParticleScalingPreparationError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_seeds(value: str) -> list[int]:
    try:
        seeds = [int(item.strip()) for item in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "seeds must be comma-separated integers"
        ) from error
    if (len(seeds) != 3 or len(seeds) != len(set(seeds))
            or any(seed < 0 or seed > 4_294_967_295 for seed in seeds)):
        raise argparse.ArgumentTypeError(
            "provide exactly three unique unsigned 32-bit seeds"
        )
    return seeds


def replace_value(text: str, section: str | None, key: str,
                  value: object) -> str:
    lines = text.splitlines()
    active: str | None = None
    matches: list[int] = []
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            active = stripped[1:-1]
            continue
        if active == section and stripped.startswith(f"{key} ="):
            matches.append(index)
    require(len(matches) == 1,
            f"expected one {key!r} value in section {section!r}")
    lines[matches[0]] = f"{key} = {value}"
    return "\n".join(lines) + "\n"


def derive_deck(base: str, level: dict[str, object], burn_in_steps: int,
                output_dir: Path) -> str:
    text = replace_value(base, None, "steps", burn_in_steps)
    text = replace_value(text, None, "output_dir", output_dir)
    text = replace_value(
        text, None, "spatial_average_start_step",
        burn_in_steps - STEPS_PER_RF_CYCLE * MEASUREMENT_RF_CYCLES_PER_BLOCK + 1,
    )
    text = replace_value(text, None, "spatial_average_end_step", burn_in_steps)
    text = replace_value(text, None, "max_particles_per_species",
                         level["max_particles_per_species"])
    for species in ("species.electrons", "species.ions"):
        text = replace_value(text, species, "particles",
                             level["particles_per_species"])
        text = replace_value(text, species, "weight", level["macro_weight"])
    return (
        "# Prospective particle-scaling arm; not the published numerical contract.\n"
        + text
    )


def prepare_campaign(args: argparse.Namespace) -> Path:
    require(args.acknowledge_cost == ACKNOWLEDGEMENT,
            f"preparation requires --acknowledge-cost {ACKNOWLEDGEMENT}")
    require(args.burn_in_rf_cycles >= 1280
            and args.burn_in_rf_cycles % MEASUREMENT_RF_CYCLES_PER_BLOCK == 0,
            "burn-in RF cycles must be a multiple of 32 and at least 1280")
    require(args.measurement_blocks >= 32,
            "at least 32 measurement blocks are required")
    executable = args.executable.resolve()
    require(executable.is_file(), f"solver executable does not exist: {executable}")
    destination = args.output_dir.resolve()
    require(not destination.exists(),
            f"refusing to overwrite particle-scaling directory: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    ))
    burn_in_steps = args.burn_in_rf_cycles * STEPS_PER_RF_CYCLE
    levels: list[dict[str, object]] = [
        {"id": "particles_1x", "factor": 1,
         "particles_per_species": 65536,
         "macro_weight": "261718750.00000003",
         "max_particles_per_species": 262144},
        {"id": "particles_2x", "factor": 2,
         "particles_per_species": 131072,
         "macro_weight": "130859375.00000001",
         "max_particles_per_species": 524288},
    ]
    records: list[dict[str, object]] = []
    try:
        for seed in args.seeds:
            source_dir = temporary / f"seed_{seed}" / "source"
            source_deck = source_dir / "turner_case1.cfg"
            source_report = source_dir / "turner_case1.preflight.json"
            prepare(argparse.Namespace(
                case_manifest=args.case_manifest,
                normalized_dir=args.normalized_dir,
                output=source_deck,
                report=source_report,
                output_dir=destination / "unused_source_output" / f"seed_{seed}",
                seed=seed,
                acknowledge_cost=TURNER_ACKNOWLEDGEMENT,
            ))
            base = source_deck.read_text(encoding="utf-8")
            for level in levels:
                identifier = str(level["id"])
                arm_dir = temporary / f"seed_{seed}" / identifier
                arm_dir.mkdir(parents=True)
                result_dir = (destination / "results" / f"seed_{seed}"
                              / identifier / "burn_in")
                deck = arm_dir / "burn_in.cfg"
                deck.write_text(derive_deck(base, level, burn_in_steps,
                                            result_dir), encoding="utf-8")
                validation = subprocess.run(
                    [str(executable), "--validate-only", str(deck)],
                    cwd=executable.parent.parent, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
                require(validation.returncode == 0,
                        f"arm seed={seed} level={identifier} failed validation: "
                        + validation.stdout + validation.stderr)
                represented = (int(level["particles_per_species"])
                               * float(str(level["macro_weight"])))
                records.append({
                    "seed": seed, "level": identifier,
                    "particle_count_factor": level["factor"],
                    "particles_per_species": level["particles_per_species"],
                    "macro_weight": float(str(level["macro_weight"])),
                    "represented_initial_particles_per_species": represented,
                    "burn_in_rf_cycles": args.burn_in_rf_cycles,
                    "burn_in_steps": burn_in_steps,
                    "measurement_blocks": args.measurement_blocks,
                    "measurement_rf_cycles": (args.measurement_blocks
                                              * MEASUREMENT_RF_CYCLES_PER_BLOCK),
                    "runtime_config": str(Path(f"seed_{seed}") / identifier
                                          / "burn_in.cfg"),
                    "runtime_config_sha256": sha256(deck),
                    "result_dir": str(result_dir),
                    "launched": False, "completed": False,
                })

        manifest = {
            "turner_particle_scaling_preparation_version": 1,
            "case_id": "turner-helium-ccp-2013-case-1",
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "case_manifest": str(args.case_manifest.resolve()),
            "case_manifest_sha256": sha256(args.case_manifest.resolve()),
            "normalization_audit": str(args.normalized_dir.resolve() / "audit.json"),
            "normalization_audit_sha256": sha256(
                args.normalized_dir.resolve() / "audit.json"),
            "solver": str(executable), "solver_sha256": sha256(executable),
            "seeds": args.seeds, "independent_seed_count": len(args.seeds),
            "levels": levels, "arms": records,
            "execution_policy": {
                "all_arms_prepared_before_launch": True,
                "maximum_concurrent_solver_processes": 1,
                "runtime_backend": "serial", "runtime_threads": 1,
                "early_stopping": False,
                "burn_in_measurements_are_acceptance_ineligible": True,
                "measurement_blocks_prepared_only_from_completed_burn_in_checkpoint": True,
            },
            "claim_boundary": (
                "Preparation does not launch a run or support a Turner validation "
                "claim. Particle scaling changes the published numerical contract."),
            "launched": False,
        }
        (temporary / "particle_scaling.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        os.replace(temporary, destination)
    except (PreparationError, OSError, ValueError, KeyError,
            subprocess.SubprocessError) as error:
        shutil.rmtree(temporary, ignore_errors=True)
        raise ParticleScalingPreparationError(str(error)) from error
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return destination / "particle_scaling.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("normalized_dir", type=Path)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seeds", type=parse_seeds, required=True)
    parser.add_argument("--burn-in-rf-cycles", type=int, default=4096)
    parser.add_argument("--measurement-blocks", type=int, default=64)
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    try:
        output = prepare_campaign(parse_args())
    except ParticleScalingPreparationError as error:
        print(f"Turner particle-scaling preparation error: {error}", file=sys.stderr)
        return 2
    print(f"Prepared Turner particle-scaling campaign without launch: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
