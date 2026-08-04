#!/usr/bin/env python3
"""Analyze phase-resolved regional EEDF histograms and exact moments."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{path} has no rows")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--threshold", action="append", default=[])
    parser.add_argument("--mirror", action="append", default=[])
    parser.add_argument("--max-overflow", type=float, default=1e-4)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    thresholds = {}
    for item in args.threshold:
        name, value = item.split(":", 1)
        thresholds[name] = float(value)
    mirrors = [tuple(item.split(":", 1)) for item in args.mirror]

    histogram_rows = read(args.output_dir / "phase_eedf.csv")
    moment_rows = read(args.output_dir / "phase_eedf_moments.csv")
    distributions: defaultdict[tuple[int, str], list[tuple[float, float, float]]] = (
        defaultdict(list)
    )
    for row in histogram_rows:
        key = (int(row["phase_bin"]), row["region"])
        distributions[key].append((
            float(row.get("energy_eV", row.get("energy_normalized", "0"))),
            float(row["represented_count"]),
            float(row["probability_density"]),
        ))
    moments = {
        (int(row["phase_bin"]), row["region"]): row
        for row in moment_rows
    }
    phases = sorted({key[0] for key in moments})
    regions = sorted({key[1] for key in moments})
    if not phases or phases != list(range(len(phases))):
        raise ValueError("phase EEDF bins must be contiguous")
    half = len(phases) // 2
    if len(phases) % 2:
        raise ValueError("phase EEDF analysis requires an even phase count")

    maximum_normalization_residual = 0.0
    maximum_overflow = 0.0
    region_reports = {}
    for region in regions:
        populated = []
        total_count = 0.0
        total_energy = 0.0
        aggregate_bins = None
        phase_thresholds = {}
        temperatures = []
        for phase in phases:
            row = moments[(phase, region)]
            count = float(row["represented_observations"])
            overflow = float(row["overflow_fraction"])
            maximum_overflow = max(maximum_overflow, overflow)
            bins = distributions[(phase, region)]
            if aggregate_bins is None:
                aggregate_bins = [0.0] * len(bins)
            histogram_count = sum(value[1] for value in bins)
            maximum_normalization_residual = max(
                maximum_normalization_residual,
                abs(histogram_count / count + overflow - 1.0)
                if count > 0.0 else abs(histogram_count),
            )
            if count <= 0.0:
                continue
            populated.append(phase)
            total_count += count
            total_energy += count * float(row["mean_energy"])
            temperatures.append(float(row["drift_separated_temperature"]))
            for index, (_, value, _) in enumerate(bins):
                aggregate_bins[index] += value
            phase_thresholds[str(phase)] = {
                name: sum(value for energy, value, _ in bins if energy >= level) /
                      count
                for name, level in thresholds.items()
            }
        aggregate_thresholds = {
            name: sum(
                aggregate_bins[index]
                for index, (energy, _, _) in enumerate(
                    distributions[(populated[0], region)] if populated else [])
                if energy >= level
            ) / total_count if total_count > 0.0 else 0.0
            for name, level in thresholds.items()
        }
        self_tv = []
        for phase in range(half):
            left = distributions[(phase, region)]
            right = distributions[(phase + half, region)]
            left_count = float(moments[(phase, region)]["represented_observations"])
            right_count = float(
                moments[(phase + half, region)]["represented_observations"])
            if left_count > 0.0 and right_count > 0.0:
                self_tv.append(0.5 * sum(
                    abs(a[1] / left_count - b[1] / right_count)
                    for a, b in zip(left, right)
                ))
        region_reports[region] = {
            "populated_phase_bins": populated,
            "empty_phase_bins": [p for p in phases if p not in populated],
            "represented_observations": total_count,
            "mean_energy": total_energy / total_count if total_count else 0.0,
            "drift_separated_temperature_min": min(temperatures) if temperatures else 0.0,
            "drift_separated_temperature_max": max(temperatures) if temperatures else 0.0,
            "aggregate_fraction_above_threshold": aggregate_thresholds,
            "phase_fraction_above_threshold": phase_thresholds,
            "half_cycle_histogram_total_variation_mean": (
                sum(self_tv) / len(self_tv) if self_tv else None
            ),
        }

    mirror_reports = {}
    for left_region, right_region in mirrors:
        values = []
        for phase in phases:
            opposite = (phase + half) % len(phases)
            left = distributions[(phase, left_region)]
            right = distributions[(opposite, right_region)]
            left_count = float(moments[(phase, left_region)]["represented_observations"])
            right_count = float(moments[(opposite, right_region)]["represented_observations"])
            if left_count > 0.0 and right_count > 0.0:
                values.append(0.5 * sum(
                    abs(a[1] / left_count - b[1] / right_count)
                    for a, b in zip(left, right)
                ))
        mirror_reports[f"{left_region}:{right_region}"] = {
            "populated_pairs": len(values),
            "half_cycle_mirrored_histogram_total_variation_mean": (
                sum(values) / len(values) if values else None
            ),
        }
    passes = (maximum_overflow <= args.max_overflow and
              maximum_normalization_residual <= 1e-10)
    report = {
        "schema_version": 1,
        "phase_bins": len(phases),
        "thresholds": thresholds,
        "regions": region_reports,
        "mirrored_regions": mirror_reports,
        "maximum_overflow_fraction": maximum_overflow,
        "maximum_histogram_normalization_residual": maximum_normalization_residual,
        "maximum_allowed_overflow_fraction": args.max_overflow,
        "passes": passes,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if passes else 1


if __name__ == "__main__":
    raise SystemExit(main())
