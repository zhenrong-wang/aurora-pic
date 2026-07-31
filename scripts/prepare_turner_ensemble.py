#!/usr/bin/env python3
"""Prepare checksum-locked full-duration Turner seed decks without launching."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

from prepare_turner_case import PreparationError, load_case, prepare


ENSEMBLE_ACKNOWLEDGEMENT = (
    "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_ENSEMBLE"
)
RUN_ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_RUN"


class EnsemblePreparationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EnsemblePreparationError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_seeds(value: str) -> list[int]:
    try:
        seeds = [int(item.strip()) for item in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "seeds must be comma-separated integers"
        ) from error
    if (
        len(seeds) < 3
        or len(seeds) > 16
        or len(seeds) != len(set(seeds))
        or any(seed < 0 or seed > 4_294_967_295 for seed in seeds)
    ):
        raise argparse.ArgumentTypeError(
            "provide 3-16 unique unsigned 32-bit seeds"
        )
    return seeds


def prepare_ensemble(args: argparse.Namespace) -> Path:
    require(
        args.acknowledge_cost == ENSEMBLE_ACKNOWLEDGEMENT,
        "ensemble generation requires --acknowledge-cost "
        f"{ENSEMBLE_ACKNOWLEDGEMENT}",
    )
    destination = args.output_dir.resolve()
    require(
        not destination.exists(),
        f"refusing to overwrite ensemble directory: {destination}",
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{destination.name}.", suffix=".tmp",
        dir=destination.parent,
    ))
    try:
        case = load_case(args.case_manifest.resolve())
        runs: list[dict[str, object]] = []
        aggregate_initial_updates = 0
        aggregate_capacity_particles = 0
        for seed in args.seeds:
            run_dir = temporary / f"seed_{seed}"
            deck = run_dir / "turner_case1.cfg"
            preflight = run_dir / "turner_case1.preflight.json"
            result_dir = destination / "results" / f"seed_{seed}"
            prepare(argparse.Namespace(
                case_manifest=args.case_manifest,
                normalized_dir=args.normalized_dir,
                output=deck,
                report=preflight,
                output_dir=result_dir,
                seed=seed,
                acknowledge_cost=RUN_ACKNOWLEDGEMENT,
            ))
            report = json.loads(preflight.read_text(encoding="utf-8"))
            contract = report["contract"]
            resources = report["resource_floor"]
            require(contract["seed"] == seed, "prepared deck seed drifted")
            aggregate_initial_updates += resources["initial_particle_updates"]
            aggregate_capacity_particles += resources["capacity_particles"]
            runs.append({
                "seed": seed,
                "runtime_config": f"seed_{seed}/turner_case1.cfg",
                "runtime_config_sha256": sha256(deck),
                "preflight_report":
                    f"seed_{seed}/turner_case1.preflight.json",
                "preflight_report_sha256": sha256(preflight),
                "result_dir": str(result_dir),
                "comparison_report": f"comparisons/seed_{seed}.json",
                "launched": False,
                "completed": False,
            })

        manifest = {
            "turner_ensemble_preparation_version": 1,
            "case_id": case["global"]["case_id"],
            "case_manifest": str(args.case_manifest.resolve()),
            "case_manifest_sha256": sha256(args.case_manifest),
            "normalization_audit": str(
                (args.normalized_dir.resolve() / "audit.json")
            ),
            "normalization_audit_sha256": sha256(
                args.normalized_dir.resolve() / "audit.json"
            ),
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "seed_count": len(args.seeds),
            "seeds": args.seeds,
            "independent_seed_contract": True,
            "launched": False,
            "launch_policy": "external_sequential_only",
            "runs": runs,
            "aggregate_resource_floor": {
                "initial_particle_updates": aggregate_initial_updates,
                "capacity_particles": aggregate_capacity_particles,
                "concurrent_runs_authorized": 1,
                "storage_bytes": "not_projected_by_preparer",
            },
            "claim_boundary": {
                "physics_claim": "none_before_all_runs_and_comparisons",
                "deck_generation_is_execution": False,
                "published_acceptance_requires_completed_independent_runs": True,
            },
            "warnings": [
                "No simulation was launched.",
                "Run at most one full Turner trajectory at a time on a "
                "shared workstation.",
                "Do not aggregate incomplete, restarted-with-changed-physics, "
                "or unverified comparison reports.",
            ],
        }
        (temporary / "ensemble.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, destination)
    except (PreparationError, OSError, ValueError, KeyError) as error:
        shutil.rmtree(temporary, ignore_errors=True)
        raise EnsemblePreparationError(str(error)) from error
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return destination / "ensemble.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("normalized_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seeds", type=parse_seeds, required=True)
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    try:
        manifest = prepare_ensemble(parse_args())
    except EnsemblePreparationError as error:
        print(f"Turner ensemble preparation error: {error}", file=sys.stderr)
        return 2
    print(f"Prepared Turner ensemble without launching it: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
