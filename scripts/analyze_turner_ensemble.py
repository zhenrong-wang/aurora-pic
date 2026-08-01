#!/usr/bin/env python3
"""Audit and summarize a complete independent-seed Turner ensemble."""

from __future__ import annotations

import argparse
import configparser
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import sys
import tempfile

from compare_turner import ComparisonError, compare


class EnsembleAnalysisError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EnsembleAnalysisError(message)


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise EnsembleAnalysisError(f"cannot read {path}: {error}") from error


def load_json(path: Path, description: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise EnsembleAnalysisError(
            f"cannot read {description}: {error}"
        ) from error
    require(isinstance(value, dict), f"{description} must be a JSON object")
    return value


def atomic_report(path: Path, value: dict[str, object]) -> None:
    require(not path.exists(), f"refusing to overwrite ensemble analysis: {path}")
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


def load_config(path: Path) -> dict[str, dict[str, str]]:
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    try:
        parser.read_string("[global]\n" + path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, configparser.Error) as error:
        raise EnsembleAnalysisError(
            f"cannot parse config {path}: {error}"
        ) from error
    return {section: dict(parser[section].items())
            for section in parser.sections()}


def density_profile(path: Path, species: str) -> tuple[list[float], list[float]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            rows = [row for row in csv.DictReader(stream)
                    if row.get("species") == species]
        x = [float(row["x_m"]) for row in rows]
        density = [float(row["number_density_mean_m-3"]) for row in rows]
    except (OSError, UnicodeError, KeyError, ValueError) as error:
        raise EnsembleAnalysisError(
            f"cannot read candidate density profile {path}: {error}"
        ) from error
    require(len(x) >= 2 and len(x) == len(density),
            f"candidate density profile is incomplete: {path}")
    require(all(math.isfinite(value) for value in x + density),
            f"candidate density profile is non-finite: {path}")
    return x, density


def reference_profile(path: Path) -> tuple[list[float], list[float]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        x = [float(row["x_m"]) for row in rows]
        density = [float(row["ion_density_mean_m-3"]) for row in rows]
    except (OSError, UnicodeError, KeyError, ValueError) as error:
        raise EnsembleAnalysisError(
            f"cannot read reference density profile {path}: {error}"
        ) from error
    require(len(x) >= 2 and len(x) == len(density),
            f"reference density profile is incomplete: {path}")
    return x, density


def trapezoid(x: list[float], y: list[float]) -> float:
    return math.fsum(
        0.5 * (y[index] + y[index + 1]) * (x[index + 1] - x[index])
        for index in range(len(x) - 1)
    )


def binomial_upper_tail(trials: int, failures: int, probability: float) -> float:
    return math.fsum(
        math.comb(trials, count)
        * probability ** count
        * (1.0 - probability) ** (trials - count)
        for count in range(failures, trials + 1)
    )


def analyze(manifest_path: Path, attachment_paths: list[Path]) -> dict[str, object]:
    manifest_path = manifest_path.resolve()
    manifest = load_json(manifest_path, "ensemble manifest")
    require(
        manifest.get("turner_ensemble_preparation_version") == 1
        and manifest.get("case_id") == "turner-helium-ccp-2013-case-1",
        "unsupported Turner ensemble manifest",
    )
    runs = manifest.get("runs")
    require(isinstance(runs, list), "ensemble manifest has no run list")
    expected_seeds = {
        run.get("seed") for run in runs if isinstance(run, dict)
    }
    require(
        len(expected_seeds) == len(runs)
        and all(isinstance(seed, int) for seed in expected_seeds),
        "ensemble manifest seeds are incomplete or duplicated",
    )
    require(
        len(attachment_paths) == len(expected_seeds),
        "one attachment is required for every prepared seed",
    )

    manifest_hash = sha256(manifest_path)
    members: list[dict[str, object]] = []
    profiles: list[list[float]] = []
    reference_x: list[float] | None = None
    reference_density: list[float] | None = None
    seen: set[int] = set()

    for attachment_path in attachment_paths:
        attachment_path = attachment_path.resolve()
        attachment = load_json(attachment_path, "ensemble attachment")
        seed = attachment.get("seed")
        require(isinstance(seed, int) and seed in expected_seeds,
                f"attachment has an unexpected seed: {seed}")
        require(seed not in seen, f"duplicate attachment for seed {seed}")
        seen.add(seed)
        require(
            attachment.get("turner_ensemble_attachment_version") == 1
            and attachment.get("ensemble_member_attached") is True,
            f"seed {seed} is not a verified ensemble attachment",
        )
        ensemble = attachment.get("ensemble")
        executed = attachment.get("executed_run")
        summary = attachment.get("comparison")
        require(
            isinstance(ensemble, dict)
            and ensemble.get("manifest_sha256") == manifest_hash
            and Path(str(ensemble.get("manifest", ""))).resolve()
                == manifest_path,
            f"seed {seed} attachment belongs to a different ensemble",
        )
        require(isinstance(executed, dict) and isinstance(summary, dict),
                f"seed {seed} attachment provenance is incomplete")
        for key, hash_key in (
            ("config", "config_sha256"),
            ("final_checkpoint", "final_checkpoint_sha256"),
        ):
            artifact = Path(str(executed.get(key, "")))
            require(sha256(artifact) == executed.get(hash_key),
                    f"seed {seed} executed {key} checksum mismatch")

        ensemble_dir = manifest_path.parent
        run = next(run for run in runs if run.get("seed") == seed)
        prepared_config = ensemble_dir / str(run["runtime_config"])
        prepared_preflight = ensemble_dir / str(run["preflight_report"])
        require(
            sha256(prepared_config) == run.get("runtime_config_sha256")
            == ensemble.get("prepared_config_sha256")
            and sha256(prepared_preflight) == run.get("preflight_report_sha256")
            == ensemble.get("prepared_preflight_sha256"),
            f"seed {seed} prepared artifact checksum mismatch",
        )
        executed_config = Path(str(executed.get("config", "")))
        prepared_values = load_config(prepared_config)
        executed_values = load_config(executed_config)
        prepared_output = prepared_values.get("global", {}).pop("output_dir", None)
        executed_output = executed_values.get("global", {}).pop("output_dir", None)
        require(
            prepared_output is not None
            and executed_output is not None
            and prepared_values == executed_values
            and executed.get("semantic_equivalence") == "exact_except_output_dir",
            f"seed {seed} executed config differs from its prepared contract",
        )

        report_path = Path(str(summary.get("report", "")))
        require(sha256(report_path) == summary.get("report_sha256"),
                f"seed {seed} comparison report checksum mismatch")
        report = load_json(report_path, f"seed {seed} comparison report")
        candidate = report.get("candidate")
        reference = report.get("reference")
        statistic = report.get("statistic")
        require(
            report.get("turner_comparison_version") == 1
            and report.get("case") == 1
            and report.get("averaging_contract_verified") is True
            and isinstance(candidate, dict)
            and isinstance(reference, dict)
            and isinstance(statistic, dict),
            f"seed {seed} comparison contract is incomplete",
        )
        candidate_path = Path(str(candidate.get("path", "")))
        metadata_path = Path(str(candidate.get("averaging_metadata_path", "")))
        reference_path = Path(str(reference.get("path", "")))
        require(
            sha256(candidate_path) == candidate.get("sha256")
            and sha256(metadata_path) == candidate.get("averaging_metadata_sha256"),
            f"seed {seed} candidate checksum mismatch",
        )
        preflight = load_json(prepared_preflight, f"seed {seed} preflight")
        provenance = preflight.get("provenance")
        require(isinstance(provenance, dict),
                f"seed {seed} preflight provenance is incomplete")
        audit_path = Path(str(provenance.get("normalization_audit", "")))
        require(
            Path(executed_output).resolve() == candidate_path.parent.resolve(),
            f"seed {seed} executed output differs from comparison location",
        )
        recomputed = compare(
            1, reference_path, candidate_path, audit_path, metadata_path,
            str(candidate.get("species", "ions")), False,
        )
        require(recomputed == report,
                f"seed {seed} comparison does not recompute exactly")
        require(
            summary.get("independently_recomputed") is True
            and summary.get("x_squared") == statistic.get("x_squared")
            and summary.get("accepted_95_percent")
                == statistic.get("accepted_95_percent")
            and summary.get("accepted_99_percent")
                == statistic.get("accepted_99_percent"),
            f"seed {seed} attachment summary differs from comparison",
        )

        x, density = density_profile(
            candidate_path, str(candidate.get("species", "ions"))
        )
        if reference_x is None:
            reference_x, reference_density = reference_profile(reference_path)
        require(
            len(x) == len(reference_x)
            and max(abs(candidate_x - reference_value)
                    for candidate_x, reference_value in zip(x, reference_x))
                <= 0.00025 * (x[-1] - x[0]) / (len(x) - 1),
            f"seed {seed} coordinate grid differs",
        )
        require(reference_density is not None, "reference density is unavailable")
        bias = 100.0 * (
            trapezoid(x, density) / trapezoid(reference_x, reference_density) - 1.0
        )
        profiles.append(density)
        members.append({
            "seed": seed,
            "attachment": str(attachment_path),
            "attachment_sha256": sha256(attachment_path),
            "comparison_report_sha256": sha256(report_path),
            "x_squared": statistic["x_squared"],
            "accepted_95_percent": statistic["accepted_95_percent"],
            "accepted_99_percent": statistic["accepted_99_percent"],
            "relative_l2": report["secondary_metrics"]["relative_l2"],
            "integrated_ion_density_bias_percent": bias,
        })

    require(seen == expected_seeds, "ensemble attachment seed set is incomplete")
    members.sort(key=lambda member: int(member["seed"]))
    biases = [float(member["integrated_ion_density_bias_percent"])
              for member in members]
    x_squared = [float(member["x_squared"]) for member in members]
    pass_95 = sum(member["accepted_95_percent"] is True for member in members)
    pass_99 = sum(member["accepted_99_percent"] is True for member in members)
    failures_99 = len(members) - pass_99
    require(reference_x is not None and reference_density is not None,
            "ensemble has no density profiles")
    mean_profile = [math.fsum(values) / len(profiles)
                    for values in zip(*profiles)]
    mean_l2 = math.sqrt(
        math.fsum((candidate - reference) ** 2
                  for candidate, reference in zip(mean_profile, reference_density))
        / math.fsum(reference ** 2 for reference in reference_density)
    )

    return {
        "turner_ensemble_analysis_version": 1,
        "case_id": manifest["case_id"],
        "ensemble": {
            "manifest": str(manifest_path),
            "manifest_sha256": manifest_hash,
            "expected_members": len(expected_seeds),
            "verified_members": len(members),
            "complete": True,
        },
        "members": members,
        "published_individual_run_results": {
            "accepted_95_percent_count": pass_95,
            "accepted_99_percent_count": pass_99,
            "failed_99_percent_count": failures_99,
            "x_squared_minimum": min(x_squared),
            "x_squared_mean": statistics.fmean(x_squared),
            "x_squared_maximum": max(x_squared),
        },
        "density_amplitude": {
            "integrated_bias_percent_mean": statistics.fmean(biases),
            "integrated_bias_percent_sample_stddev": statistics.stdev(biases)
                if len(biases) > 1 else 0.0,
            "integrated_bias_percent_minimum": min(biases),
            "integrated_bias_percent_maximum": max(biases),
            "all_member_biases_same_sign": all(value > 0.0 for value in biases)
                or all(value < 0.0 for value in biases),
            "ensemble_mean_profile_relative_l2": mean_l2,
        },
        "nominal_independence_diagnostic": {
            "assumed_individual_99_percent_failure_probability": 0.01,
            "probability_of_at_least_observed_failures": binomial_upper_tail(
                len(members), failures_99, 0.01
            ),
            "interpretation": "diagnostic_only_not_an_ensemble_acceptance_rule",
        },
        "classification": (
            "all_members_passed_published_99_percent"
            if failures_99 == 0
            else "one_or_more_members_failed_published_99_percent"
        ),
        "formal_ensemble_acceptance_rule": "none_predeclared",
        "physics_claim": "independent_seed_descriptive_evidence_only",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ensemble_manifest", type=Path)
    parser.add_argument("--attachments", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = analyze(args.ensemble_manifest, args.attachments)
        atomic_report(args.output.resolve(), report)
    except (EnsembleAnalysisError, ComparisonError, OSError) as error:
        print(f"Turner ensemble analysis error: {error}", file=sys.stderr)
        return 2
    print(f"Turner ensemble analysis written to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
