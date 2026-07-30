#!/usr/bin/env python3
"""Apply a predeclared convergence screen to adjacent Hall flux blocks."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile


METRICS = (
    "left_electron_charge_rate_a",
    "left_ion_charge_rate_a",
    "left_net_charge_rate_a",
    "reverse_step_fraction",
    "reverse_mean_demand_macroparticles",
)


class BlockComparisonError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise BlockComparisonError(f"cannot hash {path}: {error}") from error


def load_report(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise BlockComparisonError(f"cannot read {path}: {error}") from error
    if not isinstance(value, dict):
        raise BlockComparisonError(f"{path} does not contain a JSON object")
    return value


def finite(value: object, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise BlockComparisonError(f"{context} is not numeric")
    result = float(value)
    if not math.isfinite(result):
        raise BlockComparisonError(f"{context} is not finite")
    return result


def integer(value: object, context: str) -> int:
    numeric = finite(value, context)
    result = int(numeric)
    if numeric != result or result < 0:
        raise BlockComparisonError(
            f"{context} is not a non-negative integer"
        )
    return result


def report_contract(
    report: dict[str, object], path: Path
) -> tuple[int, int, int, dict[str, object]]:
    if report.get("schema_version") != 1:
        raise BlockComparisonError(f"{path} has unsupported schema_version")
    start = integer(report.get("start_step"), f"{path} start_step")
    end = integer(report.get("end_step"), f"{path} end_step")
    window = integer(report.get("window_steps"), f"{path} window_steps")
    if end <= start or window == 0 or (end - start) % window:
        raise BlockComparisonError(f"{path} has an invalid step range")
    summary = report.get("summary")
    if not isinstance(summary, dict):
        raise BlockComparisonError(f"{path} has no summary object")
    missing = set(METRICS) - set(summary)
    if missing:
        raise BlockComparisonError(
            f"{path} summary is missing metrics {sorted(missing)}"
        )
    return start, end, window, summary


def metric_value(
    summary: dict[str, object], metric: str, field: str, path: Path
) -> float:
    value = summary.get(metric)
    if not isinstance(value, dict):
        raise BlockComparisonError(f"{path} metric {metric} is not an object")
    return finite(value.get(field), f"{path} {metric} {field}")


def atomic_json(path: Path, value: dict[str, object]) -> None:
    if path.exists():
        raise BlockComparisonError(f"refusing to overwrite report: {path}")
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


def compare(args: argparse.Namespace) -> dict[str, object]:
    previous_path = args.previous_report.resolve()
    current_path = args.current_report.resolve()
    previous = load_report(previous_path)
    current = load_report(current_path)
    previous_start, previous_end, previous_window, previous_summary = (
        report_contract(previous, previous_path)
    )
    current_start, current_end, current_window, current_summary = (
        report_contract(current, current_path)
    )
    if previous_end != current_start:
        raise BlockComparisonError(
            "stationarity reports must describe adjacent blocks"
        )
    if previous_window != current_window:
        raise BlockComparisonError(
            "stationarity reports must use the same window_steps"
        )

    comparisons: dict[str, object] = {}
    passed = True
    for metric in METRICS:
        previous_mean = metric_value(
            previous_summary, metric, "mean", previous_path
        )
        current_mean = metric_value(
            current_summary, metric, "mean", current_path
        )
        current_cv = metric_value(
            current_summary, metric, "coefficient_of_variation", current_path
        )
        if current_cv < 0.0:
            raise BlockComparisonError(
                f"{current_path} {metric} coefficient_of_variation is negative"
            )
        if previous_mean == 0.0:
            raise BlockComparisonError(
                f"cannot normalize zero previous mean for {metric}"
            )
        relative_change = abs(current_mean - previous_mean) / abs(previous_mean)
        metric_passed = (
            relative_change <= args.max_mean_relative_change
            and current_cv <= args.max_window_cv
        )
        passed = passed and metric_passed
        comparisons[metric] = {
            "previous_mean": previous_mean,
            "current_mean": current_mean,
            "mean_absolute_relative_change": relative_change,
            "current_coefficient_of_variation": current_cv,
            "passed": metric_passed,
        }

    return {
        "schema_version": 1,
        "physics_claim": "none",
        "stationarity_claim": "none",
        "stationarity_screen_passed": passed,
        "criteria": {
            "max_mean_absolute_relative_change":
                args.max_mean_relative_change,
            "max_current_block_coefficient_of_variation":
                args.max_window_cv,
            "all_metrics_must_pass": True,
        },
        "previous_block": {
            "start_step": previous_start,
            "end_step": previous_end,
            "report": str(previous_path),
            "sha256": sha256(previous_path),
        },
        "current_block": {
            "start_step": current_start,
            "end_step": current_end,
            "report": str(current_path),
            "sha256": sha256(current_path),
        },
        "window_steps": current_window,
        "metrics": comparisons,
        "limitations": [
            "This is a same-seed adjacent-block convergence screen.",
            "Passing is necessary but not sufficient for a physics claim.",
            "Independent seeds, refinement, and reference agreement remain "
            "separate gates.",
        ],
    }


def fraction(value: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be numeric") from error
    if not math.isfinite(result) or result < 0.0:
        raise argparse.ArgumentTypeError("must be finite and non-negative")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Screen adjacent Hall flux blocks for convergence"
    )
    parser.add_argument("previous_report", type=Path)
    parser.add_argument("current_report", type=Path)
    parser.add_argument(
        "--max-mean-relative-change", type=fraction, required=True
    )
    parser.add_argument("--max-window-cv", type=fraction, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = compare(args)
        atomic_json(args.report, report)
    except BlockComparisonError as error:
        print(f"Hall flux block comparison error: {error}", file=sys.stderr)
        return 2
    if not report["stationarity_screen_passed"]:
        print("Hall flux blocks did not meet stationarity criteria")
        return 1
    print("Hall flux blocks met stationarity criteria")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
