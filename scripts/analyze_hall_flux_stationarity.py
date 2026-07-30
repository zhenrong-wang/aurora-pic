#!/usr/bin/env python3
"""Analyze adjacent Hall boundary-flux and controller windows."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import sys
import tempfile


class StationarityError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise StationarityError(f"cannot hash {path}: {error}") from error


def rows(path: Path, required: set[str]) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            fields = set(reader.fieldnames or ())
            if not required <= fields:
                raise StationarityError(
                    f"{path.name} is missing columns "
                    f"{sorted(required - fields)}"
                )
            result = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise StationarityError(f"cannot read {path}: {error}") from error
    if not result:
        raise StationarityError(f"{path.name} has no rows")
    return result


def finite(value: str, context: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise StationarityError(f"{context} is not numeric") from error
    if not math.isfinite(result):
        raise StationarityError(f"{context} is not finite")
    return result


def integer(value: str, context: str) -> int:
    numeric = finite(value, context)
    result = int(numeric)
    if numeric != result or result < 0:
        raise StationarityError(
            f"{context} is not a non-negative integer"
        )
    return result


def coefficient_of_variation(values: list[float]) -> float | None:
    mean = statistics.fmean(values)
    return statistics.pstdev(values) / abs(mean) if mean != 0.0 else None


def atomic_json(path: Path, value: dict[str, object]) -> None:
    if path.exists():
        raise StationarityError(f"refusing to overwrite report: {path}")
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


def analyze(args: argparse.Namespace) -> dict[str, object]:
    output = args.output_dir.resolve()
    if (
        args.start_step < 0
        or args.end_step <= args.start_step
        or args.window_steps <= 0
        or (args.end_step - args.start_step) % args.window_steps
    ):
        raise StationarityError(
            "step range must contain an integral positive number of windows"
        )
    flux_path = output / "boundary_flux.csv"
    flux_rows = rows(
        flux_path,
        {
            "step", "time", "window_start_step", "window_duration",
            "species", "boundary", "represented_charge", "charge_rate",
        },
    )
    current_path = output / "current_source.csv"
    current_rows = rows(
        current_path,
        {
            "step", "reverse_distribution_steps",
            "reverse_one_macro_steps", "reverse_two_macro_steps",
            "reverse_multi_macro_steps",
            "distributed_reverse_demand_macroparticles",
        },
    )
    current_by_step: dict[int, dict[str, str]] = {}
    for row in current_rows:
        step = integer(row["step"], "current-source step")
        if step in current_by_step:
            raise StationarityError("current-source steps are duplicated")
        current_by_step[step] = row

    windows: list[dict[str, object]] = []
    for start in range(
        args.start_step, args.end_step, args.window_steps
    ):
        end = start + args.window_steps
        if start not in current_by_step or end not in current_by_step:
            raise StationarityError(
                f"current-source rows do not bracket window {start}-{end}"
            )
        boundary_results: dict[str, object] = {}
        for boundary in ("left", "right"):
            species_results: dict[str, object] = {}
            species_rates: dict[str, float] = {}
            interval_rates: dict[str, list[float]] = {}
            for species in ("electrons", "ions"):
                selected = [
                    row for row in flux_rows
                    if start < integer(row["step"], "boundary-flux step")
                    <= end
                    and row["boundary"] == boundary
                    and row["species"] == species
                ]
                if not selected:
                    raise StationarityError(
                        f"missing {species}/{boundary} rows in "
                        f"window {start}-{end}"
                    )
                steps = [
                    integer(row["step"], "boundary-flux step")
                    for row in selected
                ]
                if (
                    steps != sorted(set(steps))
                    or steps[-1] != end
                    or integer(
                        selected[0]["window_start_step"],
                        "boundary-flux window start",
                    ) != start
                ):
                    raise StationarityError(
                        f"{species}/{boundary} flux cadence is incomplete"
                    )
                durations = [
                    finite(row["window_duration"], "flux duration")
                    for row in selected
                ]
                if any(value <= 0.0 for value in durations):
                    raise StationarityError(
                        "physical boundary-flux rows need positive duration"
                    )
                charges = [
                    finite(row["represented_charge"], "represented charge")
                    for row in selected
                ]
                rates = [
                    finite(row["charge_rate"], "charge rate")
                    for row in selected
                ]
                duration = sum(durations)
                rate = sum(charges) / duration
                species_rates[species] = rate
                interval_rates[species] = rates
                species_results[species] = {
                    "charge_rate_a": rate,
                    "interval_stddev_a": statistics.pstdev(rates),
                    "interval_coefficient_of_variation":
                        coefficient_of_variation(rates),
                    "intervals": len(rates),
                }
            net_rates = [
                electron + ion
                for electron, ion in zip(
                    interval_rates["electrons"],
                    interval_rates["ions"],
                )
            ]
            boundary_results[boundary] = {
                "species": species_results,
                "net_charge_rate_a":
                    species_rates["electrons"] + species_rates["ions"],
                "net_interval_stddev_a": statistics.pstdev(net_rates),
                "positive_net_interval_fraction":
                    sum(value > 0.0 for value in net_rates) / len(net_rates),
            }
        first = current_by_step[start]
        final = current_by_step[end]

        def delta(key: str) -> float:
            return finite(final[key], key) - finite(first[key], key)

        reverse_steps = integer(
            str(round(delta("reverse_distribution_steps"))),
            "window reverse steps",
        )
        one_steps = round(delta("reverse_one_macro_steps"))
        two_steps = round(delta("reverse_two_macro_steps"))
        multi_steps = round(delta("reverse_multi_macro_steps"))
        if (
            reverse_steps > args.window_steps
            or min(one_steps, two_steps, multi_steps) < 0
            or one_steps + two_steps + multi_steps != reverse_steps
        ):
            raise StationarityError(
                f"reverse distribution is inconsistent in {start}-{end}"
            )
        demand = delta(
            "distributed_reverse_demand_macroparticles"
        )
        windows.append({
            "start_step": start,
            "end_step": end,
            "boundaries": boundary_results,
            "reverse": {
                "steps": reverse_steps,
                "step_fraction": reverse_steps / args.window_steps,
                "one_macro_steps": one_steps,
                "two_macro_steps": two_steps,
                "multi_macro_steps": multi_steps,
                "demand_macroparticles": demand,
                "mean_demand_macroparticles":
                    demand / reverse_steps if reverse_steps else 0.0,
            },
        })

    summary_metrics: dict[str, object] = {}
    for name, values in {
        "left_electron_charge_rate_a": [
            window["boundaries"]["left"]["species"]["electrons"][
                "charge_rate_a"
            ]
            for window in windows
        ],
        "left_ion_charge_rate_a": [
            window["boundaries"]["left"]["species"]["ions"][
                "charge_rate_a"
            ]
            for window in windows
        ],
        "left_net_charge_rate_a": [
            window["boundaries"]["left"]["net_charge_rate_a"]
            for window in windows
        ],
        "reverse_step_fraction": [
            window["reverse"]["step_fraction"] for window in windows
        ],
        "reverse_mean_demand_macroparticles": [
            window["reverse"]["mean_demand_macroparticles"]
            for window in windows
        ],
    }.items():
        numeric = [float(value) for value in values]
        summary_metrics[name] = {
            "mean": statistics.fmean(numeric),
            "minimum": min(numeric),
            "maximum": max(numeric),
            "coefficient_of_variation":
                coefficient_of_variation(numeric),
        }
    return {
        "schema_version": 1,
        "physics_claim": "none",
        "stationarity_claim": "none",
        "output_dir": str(output),
        "start_step": args.start_step,
        "end_step": args.end_step,
        "window_steps": args.window_steps,
        "windows": windows,
        "summary": summary_metrics,
        "sha256": {
            "boundary_flux": sha256(flux_path),
            "current_source": sha256(current_path),
        },
        "limitations": [
            "This report describes adjacent windows from one random seed.",
            "No stationarity threshold is fitted or applied.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze adjacent Hall boundary-flux windows"
    )
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--start-step", type=int, required=True)
    parser.add_argument("--end-step", type=int, required=True)
    parser.add_argument("--window-steps", type=int, required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = analyze(args)
        atomic_json(args.report, report)
    except StationarityError as error:
        print(f"Hall flux stationarity error: {error}", file=sys.stderr)
        return 2
    print(f"Hall flux stationarity report written: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
