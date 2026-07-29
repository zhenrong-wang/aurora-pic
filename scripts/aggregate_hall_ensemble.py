#!/usr/bin/env python3
"""Aggregate independent seeded Hall comparison reports conservatively."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import sys
import tempfile


class EnsembleError(RuntimeError):
    pass


T95 = {
    2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
    7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201,
    12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131, 16: 2.120,
    17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086, 21: 2.080,
    22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060, 26: 2.056,
    27: 2.052, 28: 2.048, 29: 2.045, 30: 2.042,
}


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise EnsembleError(f"cannot hash {path}: {error}") from error


def load_json(path: Path, label: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise EnsembleError(f"cannot read {label} {path}: {error}") from error
    if not isinstance(value, dict):
        raise EnsembleError(f"{label} must contain a JSON object")
    return value


def finite(value: object, context: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise EnsembleError(f"{context} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise EnsembleError(f"{context} must be finite")
    return result


def t95(seed_count: int) -> float:
    degrees = seed_count - 1
    if degrees <= 30:
        return T95[degrees]
    if degrees <= 40:
        return 2.04
    if degrees <= 60:
        return 2.02
    return 2.00


def summarize(
    values: list[float],
    reference: float,
    uncertainty: float,
    threshold: float,
    individual_passes: int,
    required_passes: int,
) -> dict[str, object]:
    mean = statistics.fmean(values)
    standard_deviation = statistics.stdev(values)
    standard_error = standard_deviation / math.sqrt(len(values))
    multiplier = t95(len(values))
    half_width = multiplier * standard_error
    bias = mean - reference
    passed = (
        abs(bias) + half_width <= threshold
        and individual_passes >= required_passes
    )
    return {
        "simulation_values": values,
        "ensemble_mean": mean,
        "sample_standard_deviation": standard_deviation,
        "standard_error": standard_error,
        "confidence_level": 0.95,
        "student_t_multiplier": multiplier,
        "confidence_half_width": half_width,
        "reference": reference,
        "reference_uncertainty": uncertainty,
        "ensemble_bias": bias,
        "acceptance_threshold": threshold,
        "conservative_error_bound": abs(bias) + half_width,
        "individual_passes": individual_passes,
        "required_individual_passes": required_passes,
        "passed": passed,
    }


def comparable_result(
    value: dict[str, object],
    context: str,
) -> tuple[float, float, float, float, bool]:
    simulation = finite(value.get("simulation"), f"{context} simulation")
    reference = finite(value.get("reference"), f"{context} reference")
    uncertainty = finite(
        value.get("reference_uncertainty"),
        f"{context} reference uncertainty",
    )
    threshold = finite(
        value.get("acceptance_threshold"),
        f"{context} acceptance threshold",
    )
    if uncertainty < 0.0 or threshold < 0.0:
        raise EnsembleError(f"{context} uncertainty/threshold must be non-negative")
    passed = value.get("passed")
    if not isinstance(passed, bool):
        raise EnsembleError(f"{context} passed must be boolean")
    return simulation, reference, uncertainty, threshold, passed


def same_float(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=0.0, abs_tol=0.0)


def aggregate(args: argparse.Namespace) -> dict[str, object]:
    manifest_path = args.ensemble_manifest.resolve()
    campaign_root = manifest_path.parent
    campaign = load_json(manifest_path, "Hall ensemble manifest")
    if campaign.get("hall_ensemble_version") != 1:
        raise EnsembleError("Hall ensemble manifest version must be 1")
    runs = campaign.get("runs")
    seeds = campaign.get("seeds")
    if not isinstance(runs, list) or not isinstance(seeds, list):
        raise EnsembleError("Hall ensemble manifest requires runs and seeds")
    if len(runs) < 3 or len(runs) != len(seeds):
        raise EnsembleError("Hall ensemble requires at least three aligned runs")
    if len(set(seeds)) != len(seeds):
        raise EnsembleError("Hall ensemble seeds must be unique")
    if campaign.get("seed_count") != len(seeds):
        raise EnsembleError("Hall ensemble seed_count is inconsistent")
    case_manifest = Path(str(campaign.get("case_manifest", "")))
    if sha256(case_manifest) != campaign.get("case_manifest_sha256"):
        raise EnsembleError("Hall ensemble case-manifest SHA-256 mismatch")

    required_passes = math.ceil(
        args.minimum_pass_fraction * len(runs) - 1e-15
    )
    reports: list[dict[str, object]] = []
    report_records: list[dict[str, object]] = []
    artifact_signatures: set[tuple[str, str, str]] = set()
    reference_identity: tuple[object, ...] | None = None
    averaging_window: dict[str, object] | None = None
    for index, run in enumerate(runs):
        if not isinstance(run, dict):
            raise EnsembleError(f"run {index} must be an object")
        seed = run.get("seed")
        if seed != seeds[index]:
            raise EnsembleError("run order/seed does not match the seed list")
        runtime_path = campaign_root / str(run.get("runtime_config", ""))
        runtime_hash = sha256(runtime_path)
        if runtime_hash != run.get("runtime_config_sha256"):
            raise EnsembleError(f"seed {seed} runtime-config SHA-256 mismatch")
        report_path = campaign_root / str(run.get("comparison_report", ""))
        report = load_json(report_path, f"seed {seed} comparison report")
        if report.get("schema_version") != 1:
            raise EnsembleError(f"seed {seed} comparison schema must be 1")
        if (
            report.get("case_id") != campaign.get("case_id")
            or not isinstance(report.get("passed"), bool)
        ):
            raise EnsembleError(f"seed {seed} report case/pass identity is invalid")
        simulation = report.get("simulation")
        reference = report.get("reference")
        if not isinstance(simulation, dict) or not isinstance(reference, dict):
            raise EnsembleError(f"seed {seed} report identity is incomplete")
        if (
            simulation.get("seed") != seed
            or simulation.get("runtime_config_sha256") != runtime_hash
            or simulation.get("case_manifest_sha256")
                != campaign.get("case_manifest_sha256")
        ):
            raise EnsembleError(f"seed {seed} report is not bound to its deck")
        if Path(str(simulation.get("runtime_config", ""))).resolve() != runtime_path.resolve():
            raise EnsembleError(f"seed {seed} runtime-config path mismatch")
        expected_result = Path(str(run.get("result_dir", ""))).resolve()
        if Path(str(simulation.get("output_dir", ""))).resolve() != expected_result:
            raise EnsembleError(f"seed {seed} result-directory mismatch")
        signature = tuple(
            str(simulation.get(key, ""))
            for key in (
                "field_average_sha256",
                "species_average_sha256",
                "mode_history_sha256",
            )
        )
        if any(len(item) != 64 for item in signature):
            raise EnsembleError(f"seed {seed} simulation hashes are incomplete")
        if signature in artifact_signatures:
            raise EnsembleError("independent seeds reuse identical simulation artifacts")
        artifact_signatures.add(signature)
        identity = (
            report.get("case_id"),
            report.get("case_variant"),
            reference.get("manifest_sha256"),
            reference.get("profile_sha256"),
            reference.get("mode_sha256"),
        )
        if reference_identity is None:
            reference_identity = identity
        elif identity != reference_identity:
            raise EnsembleError("seed reports do not share one reference identity")
        window = report.get("averaging_window")
        if not isinstance(window, dict):
            raise EnsembleError(f"seed {seed} averaging window is invalid")
        if averaging_window is None:
            averaging_window = window
        elif window != averaging_window:
            raise EnsembleError("seed reports use different averaging windows")
        reports.append(report)
        report_records.append({
            "seed": seed,
            "comparison_report": str(report_path.resolve()),
            "comparison_report_sha256": sha256(report_path),
            "individual_passed": report.get("passed"),
            "runtime_config_sha256": runtime_hash,
        })

    profile_shape: list[tuple[float, list[str]]] | None = None
    profile_values: dict[tuple[float, str], list[tuple[float, bool]]] = {}
    profile_contract: dict[tuple[float, str], tuple[float, float, float]] = {}
    for report_index, report in enumerate(reports):
        rows = report.get("profile_comparisons")
        if not isinstance(rows, list):
            raise EnsembleError("profile_comparisons must be a list")
        current_shape: list[tuple[float, list[str]]] = []
        for row_index, row in enumerate(rows):
            if not isinstance(row, dict) or not isinstance(row.get("observables"), list):
                raise EnsembleError("profile comparison row is invalid")
            coordinate = finite(
                row.get("reference_coordinate"),
                f"profile row {row_index} coordinate",
            )
            names: list[str] = []
            for observable in row["observables"]:
                if not isinstance(observable, dict) or not isinstance(
                    observable.get("name"), str
                ):
                    raise EnsembleError("profile observable is invalid")
                name = str(observable["name"])
                names.append(name)
                simulation, reference, uncertainty, threshold, passed = (
                    comparable_result(observable, f"profile {coordinate} {name}")
                )
                key = (coordinate, name)
                contract = (reference, uncertainty, threshold)
                if key in profile_contract and any(
                    not same_float(left, right)
                    for left, right in zip(profile_contract[key], contract)
                ):
                    raise EnsembleError("profile reference contract differs by seed")
                profile_contract[key] = contract
                profile_values.setdefault(key, []).append((simulation, passed))
            current_shape.append((coordinate, names))
        if profile_shape is None:
            profile_shape = current_shape
        elif current_shape != profile_shape:
            raise EnsembleError(
                f"profile comparison shape differs for report {report_index}"
            )

    mode_shape: list[tuple[str, int, str]] | None = None
    mode_values: dict[tuple[str, int, str], list[tuple[float, bool]]] = {}
    mode_contract: dict[
        tuple[str, int, str], tuple[float, float, float]
    ] = {}
    for report_index, report in enumerate(reports):
        rows = report.get("mode_comparisons")
        if not isinstance(rows, list):
            raise EnsembleError("mode_comparisons must be a list")
        current_shape: list[tuple[str, int, str]] = []
        for row in rows:
            if (
                not isinstance(row, dict)
                or not isinstance(row.get("name"), str)
                or not isinstance(row.get("mode"), int)
                or not isinstance(row.get("metric"), str)
            ):
                raise EnsembleError("mode comparison row is invalid")
            key = (str(row["name"]), int(row["mode"]), str(row["metric"]))
            current_shape.append(key)
            simulation, reference, uncertainty, threshold, passed = (
                comparable_result(row, f"mode {key}")
            )
            contract = (reference, uncertainty, threshold)
            if key in mode_contract and any(
                not same_float(left, right)
                for left, right in zip(mode_contract[key], contract)
            ):
                raise EnsembleError("mode reference contract differs by seed")
            mode_contract[key] = contract
            mode_values.setdefault(key, []).append((simulation, passed))
        if mode_shape is None:
            mode_shape = current_shape
        elif current_shape != mode_shape:
            raise EnsembleError(
                f"mode comparison shape differs for report {report_index}"
            )

    profile_results = []
    for key, values in profile_values.items():
        reference, uncertainty, threshold = profile_contract[key]
        summary = summarize(
            [item[0] for item in values],
            reference, uncertainty, threshold,
            sum(item[1] for item in values), required_passes,
        )
        profile_results.append({
            "reference_coordinate": key[0],
            "name": key[1],
            **summary,
        })
    mode_results = []
    for key, values in mode_values.items():
        reference, uncertainty, threshold = mode_contract[key]
        summary = summarize(
            [item[0] for item in values],
            reference, uncertainty, threshold,
            sum(item[1] for item in values), required_passes,
        )
        mode_results.append({
            "name": key[0], "mode": key[1], "metric": key[2], **summary,
        })
    if not profile_results or not mode_results:
        raise EnsembleError(
            "ensemble reports require profile and mode comparisons"
        )
    individual_reports_passed = sum(
        record["individual_passed"] is True for record in report_records
    )
    passed = (
        individual_reports_passed >= required_passes
        and all(item["passed"] for item in profile_results)
        and all(item["passed"] for item in mode_results)
    )
    return {
        "hall_ensemble_report_version": 1,
        "passed": passed,
        "physics_claim_eligible": passed and campaign.get("tier") == "production",
        "case_id": campaign.get("case_id"),
        "tier": campaign.get("tier"),
        "seed_count": len(seeds),
        "seeds": seeds,
        "minimum_pass_fraction": args.minimum_pass_fraction,
        "required_individual_passes": required_passes,
        "individual_reports_passed": individual_reports_passed,
        "confidence_policy":
            "two-sided 95% Student-t interval must fit inside acceptance band",
        "averaging_window": averaging_window,
        "reference_identity": {
            "case_id": reference_identity[0],
            "case_variant": reference_identity[1],
            "manifest_sha256": reference_identity[2],
            "profile_sha256": reference_identity[3],
            "mode_sha256": reference_identity[4],
        },
        "ensemble_manifest": str(manifest_path),
        "ensemble_manifest_sha256": sha256(manifest_path),
        "runs": report_records,
        "profile_ensemble": profile_results,
        "mode_ensemble": mode_results,
    }


def write_atomic(path: Path, report: dict[str, object]) -> None:
    if path.exists():
        raise EnsembleError(f"refusing to overwrite ensemble report: {path}")
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate three or more seeded Hall comparisons"
    )
    parser.add_argument("ensemble_manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--minimum-pass-fraction", type=float, default=2.0 / 3.0
    )
    args = parser.parse_args()
    if (
        not math.isfinite(args.minimum_pass_fraction)
        or args.minimum_pass_fraction <= 0.0
        or args.minimum_pass_fraction > 1.0
    ):
        parser.error("--minimum-pass-fraction must be in (0,1]")
    return args


def main() -> int:
    args = parse_args()
    try:
        report = aggregate(args)
        write_atomic(args.output, report)
    except EnsembleError as error:
        print(f"Hall ensemble input error: {error}", file=sys.stderr)
        return 2
    if not report["passed"]:
        print("Hall ensemble did not meet acceptance criteria", file=sys.stderr)
        return 1
    print(f"Hall ensemble comparison passed: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
