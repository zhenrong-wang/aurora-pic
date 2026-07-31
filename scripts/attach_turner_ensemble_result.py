#!/usr/bin/env python3
"""Attach one verified full-duration result to a prepared Turner ensemble."""

from __future__ import annotations

import argparse
import configparser
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

from compare_turner import ComparisonError, compare


class AttachmentError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AttachmentError(message)


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise AttachmentError(f"cannot read {path}: {error}") from error


def load_json(path: Path, description: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise AttachmentError(f"cannot read {description}: {error}") from error
    require(isinstance(value, dict), f"{description} must be a JSON object")
    return value


def load_config(path: Path) -> dict[str, dict[str, str]]:
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    try:
        parser.read_string("[global]\n" + path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, configparser.Error) as error:
        raise AttachmentError(f"cannot parse config {path}: {error}") from error
    return {
        section: dict(parser[section].items())
        for section in parser.sections()
    }


def atomic_report(path: Path, value: dict[str, object]) -> None:
    require(not path.exists(), f"refusing to overwrite attachment: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
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


def attach(args: argparse.Namespace) -> dict[str, object]:
    manifest_path = args.ensemble_manifest.resolve()
    manifest = load_json(manifest_path, "ensemble manifest")
    require(
        manifest.get("turner_ensemble_preparation_version") == 1
        and manifest.get("case_id") == "turner-helium-ccp-2013-case-1"
        and manifest.get("launched") is False,
        "unsupported or already-mutated Turner ensemble manifest",
    )
    runs = manifest.get("runs")
    require(isinstance(runs, list), "ensemble manifest has no run list")
    matches = [
        run for run in runs
        if isinstance(run, dict) and run.get("seed") == args.seed
    ]
    require(len(matches) == 1, f"seed {args.seed} is not unique in ensemble")
    run = matches[0]
    require(
        run.get("launched") is False and run.get("completed") is False,
        "prepared run is not in its unattached state",
    )

    ensemble_dir = manifest_path.parent
    prepared_config = ensemble_dir / str(run["runtime_config"])
    preflight_path = ensemble_dir / str(run["preflight_report"])
    require(
        sha256(prepared_config) == run.get("runtime_config_sha256"),
        "prepared runtime config checksum mismatch",
    )
    require(
        sha256(preflight_path) == run.get("preflight_report_sha256"),
        "prepared preflight checksum mismatch",
    )
    preflight = load_json(preflight_path, "prepared preflight report")
    require(
        preflight.get("turner_case_preflight_version") == 1
        and preflight.get("full_run_launched") is False
        and preflight.get("production_launch_authorized") is False,
        "prepared preflight report has an invalid claim boundary",
    )
    contract = preflight.get("contract")
    provenance = preflight.get("provenance")
    require(
        isinstance(contract, dict)
        and contract.get("seed") == args.seed
        and contract.get("steps") == 512000
        and contract.get("averaging_samples") == 12800
        and isinstance(provenance, dict)
        and provenance.get("generated_deck_sha256") == sha256(prepared_config),
        "prepared preflight contract differs from the seed deck",
    )

    executed_config = args.executed_config.resolve()
    prepared_values = load_config(prepared_config)
    executed_values = load_config(executed_config)
    prepared_output = prepared_values.get("global", {}).pop("output_dir", None)
    executed_output = executed_values.get("global", {}).pop("output_dir", None)
    require(
        prepared_output is not None and executed_output is not None,
        "prepared or executed config has no output_dir",
    )
    require(
        prepared_values == executed_values,
        "executed config differs from the prepared seed contract beyond output_dir",
    )

    comparison_path = args.comparison_report.resolve()
    comparison = load_json(comparison_path, "comparison report")
    require(
        comparison.get("turner_comparison_version") == 1
        and comparison.get("case") == 1
        and comparison.get("comparison_scope")
            == "published_baseline_ion_density_statistic_only"
        and comparison.get("averaging_contract_verified") is True,
        "comparison is not a complete published-duration Turner Case 1 result",
    )
    candidate = comparison.get("candidate")
    reference = comparison.get("reference")
    require(
        isinstance(candidate, dict) and isinstance(reference, dict),
        "comparison provenance is incomplete",
    )
    profile_path = Path(str(candidate.get("path", "")))
    metadata_path = Path(str(candidate.get("averaging_metadata_path", "")))
    reference_path = Path(str(reference.get("path", "")))
    require(
        sha256(profile_path) == candidate.get("sha256")
        and sha256(metadata_path) == candidate.get("averaging_metadata_sha256"),
        "comparison candidate artifact checksum mismatch",
    )
    normalization_audit = Path(str(provenance.get("normalization_audit", "")))
    require(
        sha256(normalization_audit)
            == provenance.get("normalization_audit_sha256")
        and reference.get("normalization_audit_sha256")
            == provenance.get("normalization_audit_sha256"),
        "comparison and preflight normalization provenance differ",
    )
    recomputed = compare(
        1, reference_path, profile_path, normalization_audit,
        metadata_path, str(candidate.get("species", "ions")), False,
    )
    require(
        comparison == recomputed,
        "stored comparison differs from an independent recomputation",
    )
    require(
        Path(executed_output).resolve() == profile_path.parent.resolve(),
        "executed config output_dir differs from comparison candidate location",
    )

    checkpoint = args.final_checkpoint.resolve()
    try:
        checkpoint_bytes = checkpoint.stat().st_size
    except OSError as error:
        raise AttachmentError(f"cannot inspect final checkpoint: {error}") from error
    require(checkpoint_bytes > 0, "final checkpoint is empty")
    statistic = comparison.get("statistic")
    require(isinstance(statistic, dict), "comparison statistic is missing")
    accepted_99 = statistic.get("accepted_99_percent") is True
    return {
        "turner_ensemble_attachment_version": 1,
        "case_id": manifest["case_id"],
        "seed": args.seed,
        "ensemble": {
            "manifest": str(manifest_path),
            "manifest_sha256": sha256(manifest_path),
            "prepared_config_sha256": sha256(prepared_config),
            "prepared_preflight_sha256": sha256(preflight_path),
        },
        "executed_run": {
            "config": str(executed_config),
            "config_sha256": sha256(executed_config),
            "semantic_equivalence": "exact_except_output_dir",
            "prepared_output_dir": prepared_output,
            "executed_output_dir": executed_output,
            "final_checkpoint": str(checkpoint),
            "final_checkpoint_sha256": sha256(checkpoint),
            "final_checkpoint_bytes": checkpoint_bytes,
        },
        "comparison": {
            "report": str(comparison_path),
            "report_sha256": sha256(comparison_path),
            "independently_recomputed": True,
            "profile_sha256": candidate["sha256"],
            "averaging_metadata_sha256":
                candidate["averaging_metadata_sha256"],
            "x_squared": statistic.get("x_squared"),
            "accepted_95_percent": statistic.get("accepted_95_percent"),
            "accepted_99_percent": accepted_99,
        },
        "classification": (
            "single_seed_published_density_statistic_passed_99_percent"
            if accepted_99
            else "single_seed_published_density_statistic_failed_99_percent"
        ),
        "ensemble_member_attached": True,
        "ensemble_complete": False,
        "physics_claim": "single_seed_published_duration_density_comparison_only",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ensemble_manifest", type=Path)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--executed-config", type=Path, required=True)
    parser.add_argument("--comparison-report", type=Path, required=True)
    parser.add_argument("--final-checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = attach(args)
        atomic_report(args.output.resolve(), report)
    except (AttachmentError, ComparisonError, OSError) as error:
        print(f"Turner ensemble attachment error: {error}", file=sys.stderr)
        return 2
    print(f"Turner ensemble result attached: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
