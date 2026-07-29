#!/usr/bin/env python3
"""Atomically prepare multiple seeded Hall decks without launching them."""

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

from prepare_hall_campaign import CampaignError, load_manifest, prepare


class EnsemblePreparationError(RuntimeError):
    pass


def parse_seeds(value: str) -> list[int]:
    try:
        seeds = [int(item.strip()) for item in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "seeds must be comma-separated integers"
        ) from error
    if (
        len(seeds) < 3
        or len(seeds) > 64
        or len(seeds) != len(set(seeds))
        or any(seed < 0 or seed > 4_294_967_295 for seed in seeds)
    ):
        raise argparse.ArgumentTypeError(
            "provide 3-64 unique unsigned 32-bit seeds"
        )
    return seeds


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare_ensemble(args: argparse.Namespace) -> Path:
    destination = args.output_dir.resolve()
    if destination.exists():
        raise EnsemblePreparationError(
            f"refusing to overwrite ensemble directory: {destination}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{destination.name}.", suffix=".tmp",
        dir=destination.parent,
    ))
    try:
        case = load_manifest(args.case_manifest)
        runs: list[dict[str, object]] = []
        for seed in args.seeds:
            deck_name = f"case_seed_{seed}.cfg"
            deck = temporary / deck_name
            result_dir = destination / "results" / f"seed_{seed}"
            prepare(argparse.Namespace(
                case_manifest=args.case_manifest,
                tier=args.tier,
                output=deck,
                output_dir=str(result_dir),
                seed=seed,
                acknowledge_cost=args.acknowledge_cost,
            ))
            runs.append({
                "seed": seed,
                "runtime_config": deck_name,
                "runtime_config_sha256": sha256(deck),
                "result_dir": str(result_dir),
                "comparison_report":
                    f"comparisons/seed_{seed}.json",
                "launched": False,
            })
        manifest_path = temporary / "ensemble.json"
        manifest = {
            "hall_ensemble_version": 1,
            "case_id": case["global"]["case_id"],
            "tier": args.tier,
            "case_manifest": str(args.case_manifest.resolve()),
            "case_manifest_sha256": sha256(args.case_manifest),
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "seed_count": len(args.seeds),
            "seeds": args.seeds,
            "launched": False,
            "runs": runs,
            "warnings": [
                "Deck generation did not launch a simulation.",
                "A physics claim requires independent completed runs and "
                "checksum-pinned comparison reports.",
            ],
        }
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, destination)
    except (CampaignError, OSError, ValueError) as error:
        shutil.rmtree(temporary, ignore_errors=True)
        raise EnsemblePreparationError(str(error)) from error
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return destination


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare three or more independent Hall seed decks"
    )
    parser.add_argument("case_manifest", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--tier",
        choices=("micro", "workstation", "production"),
        default="production",
    )
    parser.add_argument("--seeds", type=parse_seeds, required=True)
    parser.add_argument("--acknowledge-cost")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        output = prepare_ensemble(args)
    except EnsemblePreparationError as error:
        print(f"Hall ensemble preparation error: {error}", file=sys.stderr)
        return 2
    print(f"Prepared Hall ensemble without launching it: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
