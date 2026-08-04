#!/usr/bin/env python3
"""Apply a predeclared block/slope stationarity gate to eduPIC conv.dat."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys


class ConvergenceError(RuntimeError):
    pass


def positive_integer(text: str) -> int:
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return value


def nonnegative(text: str) -> float:
    value = float(text)
    if not math.isfinite(value) or value < 0.0:
        raise argparse.ArgumentTypeError("value must be finite and nonnegative")
    return value


def read_history(path: Path) -> list[dict[str, int]]:
    rows = []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise ConvergenceError(f"cannot read {path}: {error}") from error
    for line_number, line in enumerate(lines, 1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 3:
            raise ConvergenceError(f"line {line_number} does not have three fields")
        try:
            cycle, electrons, ions = map(int, fields)
        except ValueError as error:
            raise ConvergenceError(f"line {line_number} has a non-integer field") from error
        if cycle <= 0 or electrons < 0 or ions < 0:
            raise ConvergenceError(f"line {line_number} has an invalid value")
        if rows and cycle != rows[-1]["cycle"] + 1:
            raise ConvergenceError("cycles are not contiguous")
        rows.append({"cycle": cycle, "electrons": electrons, "ions": ions,
                     "total_particles": electrons + ions})
    if not rows:
        raise ConvergenceError("convergence history is empty")
    return rows


def metric(rows: list[dict[str, int]], field: str, block_cycles: int) -> dict:
    x = [float(row["cycle"]) for row in rows]
    y = [float(row[field]) for row in rows]
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    denominator = sum((value - mean_x) ** 2 for value in x)
    slope = (sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y)) /
             denominator if denominator > 0.0 else 0.0)
    relative_slope = slope / mean_y if mean_y > 0.0 else None
    block_means = []
    for start in range(0, len(rows), block_cycles):
        values = y[start:start + block_cycles]
        if len(values) == block_cycles:
            block_means.append(sum(values) / len(values))
    block_range = ((max(block_means) - min(block_means)) / mean_y
                   if mean_y > 0.0 and len(block_means) >= 2 else None)
    endpoint_change = ((y[-1] - y[0]) / mean_y if mean_y > 0.0 else None)
    return {"mean": mean_y, "linear_slope_particles_per_cycle": slope,
            "relative_slope_per_cycle": relative_slope,
            "endpoint_change_relative_to_mean": endpoint_change,
            "block_means": block_means,
            "block_mean_range_relative_to_window_mean": block_range}


def analyze(rows: list[dict[str, int]], window_cycles: int, block_cycles: int,
            max_relative_slope: float, max_block_range: float,
            minimum_cycle: int = 1500,
            provisional_recent_cycles: int = 25) -> dict:
    if window_cycles % block_cycles:
        raise ConvergenceError("window cycles must be an integer multiple of block cycles")
    required_blocks = window_cycles // block_cycles
    if required_blocks < 4:
        raise ConvergenceError("stationarity window must contain at least four blocks")
    if provisional_recent_cycles < 2:
        raise ConvergenceError("provisional recent window must contain at least two cycles")
    eligible = len(rows) >= window_cycles and rows[-1]["cycle"] >= minimum_cycle
    window = rows[-window_cycles:] if eligible else rows
    metrics = {field: metric(window, field, block_cycles)
               for field in ("electrons", "ions", "total_particles")}
    criteria = {
        field: {
            "relative_slope_passes": eligible and
                value["relative_slope_per_cycle"] is not None and
                abs(value["relative_slope_per_cycle"]) <= max_relative_slope,
            "block_range_passes": eligible and
                value["block_mean_range_relative_to_window_mean"] is not None and
                value["block_mean_range_relative_to_window_mean"] <= max_block_range,
        } for field, value in metrics.items()
    }
    stationary = eligible and all(
        item["relative_slope_passes"] and item["block_range_passes"]
        for item in criteria.values())
    recent = rows[-min(len(rows), provisional_recent_cycles):]
    recent_metrics = {field: metric(recent, field, len(recent))
                      for field in ("electrons", "ions", "total_particles")}
    return {
        "schema_version": 1, "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "population_stationarity_gate", "physics_claim": "none",
        "claim_boundary": "Population stationarity is necessary for measurement but does not validate plasma observables or cross-code agreement.",
        "history": {"available_cycles": len(rows),
                    "first_cycle": rows[0]["cycle"],
                    "last_cycle": rows[-1]["cycle"]},
        "contract": {"window_cycles": window_cycles,
                     "block_cycles": block_cycles,
                     "required_blocks": required_blocks,
                     "minimum_cycle": minimum_cycle,
                     "provisional_recent_cycles": provisional_recent_cycles,
                     "maximum_absolute_relative_slope_per_cycle": max_relative_slope,
                     "maximum_block_mean_range_relative_to_window_mean": max_block_range},
        "window": {"eligible": eligible, "samples": len(window),
                   "first_cycle": window[0]["cycle"],
                   "last_cycle": window[-1]["cycle"]},
        "metrics": metrics, "criteria": criteria, "stationary": stationary,
        "provisional_recent": {
            "claim": "descriptive_only_not_a_stationarity_gate",
            "samples": len(recent), "first_cycle": recent[0]["cycle"],
            "last_cycle": recent[-1]["cycle"], "metrics": recent_metrics,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("convergence", type=Path)
    parser.add_argument("--window-cycles", type=positive_integer, default=100)
    parser.add_argument("--block-cycles", type=positive_integer, default=25)
    parser.add_argument("--minimum-cycle", type=positive_integer, default=1500)
    parser.add_argument("--provisional-recent-cycles", type=positive_integer,
                        default=25)
    parser.add_argument("--max-relative-slope-per-cycle", type=nonnegative,
                        default=1e-4)
    parser.add_argument("--max-block-range-relative", type=nonnegative,
                        default=0.02)
    parser.add_argument("--require-stationary", action="store_true")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    try:
        report = analyze(read_history(args.convergence), args.window_cycles,
                         args.block_cycles, args.max_relative_slope_per_cycle,
                         args.max_block_range_relative, args.minimum_cycle,
                         args.provisional_recent_cycles)
    except ConvergenceError as error:
        print(f"eduPIC convergence analysis failed: {error}", file=sys.stderr)
        return 2
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 1 if args.require_stationary and not report["stationary"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
