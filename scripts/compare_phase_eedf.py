#!/usr/bin/env python3
"""Compare two solver-neutral phase/region EEDF interchanges."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def load(directory: Path) -> tuple[dict, dict, dict]:
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    if (manifest.get("schema") != "aurorapic.phase-eedf-interchange" or
            manifest.get("schema_version") != 1):
        raise ValueError(f"unsupported EEDF interchange in {directory}")

    def rows(name: str) -> list[dict[str, str]]:
        path = directory / manifest["files"][name]
        with path.open(newline="", encoding="utf-8") as stream:
            return list(csv.DictReader(stream))

    distributions: defaultdict[tuple[int, str], list[tuple[float, float, float]]] = (
        defaultdict(list))
    for row in rows("distributions"):
        distributions[(int(row["phase_bin"]), row["region"])].append((
            float(row["energy_lower_eV"]), float(row["energy_upper_eV"]),
            float(row["probability_mass"])))
    moments = {(int(row["phase_bin"]), row["region"]): row
               for row in rows("moments")}
    if not distributions or set(distributions) != set(moments):
        raise ValueError(f"distribution/moment coverage is empty or inconsistent in {directory}")
    for key, values in distributions.items():
        values.sort()
        mass_sum = 0.0
        previous_high = None
        for low, high, mass in values:
            if not all(math.isfinite(item) for item in (low, high, mass)):
                raise ValueError(f"non-finite distribution value for {key}")
            if high <= low or mass < 0.0:
                raise ValueError(f"invalid energy bin or probability mass for {key}")
            if previous_high is not None and low < previous_high - 1e-12:
                raise ValueError(f"overlapping energy bins for {key}")
            previous_high = high
            mass_sum += mass
        overflow = float(moments[key]["overflow_fraction"])
        represented = float(moments[key]["represented_observations"])
        expected_total = 1.0 if represented > 0.0 else 0.0
        if (not math.isfinite(overflow) or not math.isfinite(represented) or
                represented < 0.0 or not 0.0 <= overflow <= 1.0 or
                abs(mass_sum + overflow - expected_total) > 1e-8):
            raise ValueError(f"distribution is not normalized for {key}")
    return manifest, dict(distributions), moments


def relative(a: float, b: float) -> float:
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale else 0.0


def tv_distance(a: list[tuple[float, float, float]],
                b: list[tuple[float, float, float]]) -> float:
    """Exact TV for piecewise-uniform histograms on their union grid."""
    edges = sorted({edge for low, high, _ in a + b for edge in (low, high)})

    def density(values: list[tuple[float, float, float]], x: float) -> float:
        for low, high, mass in values:
            if low <= x < high:
                return mass / (high - low)
        return 0.0

    integral = 0.0
    for low, high in zip(edges, edges[1:]):
        midpoint = 0.5 * (low + high)
        integral += abs(density(a, midpoint) - density(b, midpoint)) * (high - low)
    return 0.5 * integral


def tail(values: list[tuple[float, float, float]], threshold: float) -> float:
    total = 0.0
    for low, high, mass in values:
        overlap = max(0.0, high - max(low, threshold))
        total += mass * overlap / (high - low)
    return total


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--tail-eV", action="append", type=float, default=[])
    parser.add_argument("--max-tv", type=float)
    parser.add_argument("--max-mean-energy-relative", type=float)
    parser.add_argument("--max-temperature-relative", type=float)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    ref_manifest, ref_dist, ref_moments = load(args.reference)
    cand_manifest, cand_dist, cand_moments = load(args.candidate)
    if ref_manifest["case_id"] != cand_manifest["case_id"]:
        raise ValueError("case_id differs; refusing an unmatched comparison")
    if ref_manifest["species"] != cand_manifest["species"]:
        raise ValueError("species differs; refusing an unmatched comparison")
    if set(ref_dist) != set(cand_dist) or set(ref_moments) != set(cand_moments):
        raise ValueError("phase/region coverage differs")

    pairs = []
    for key in sorted(ref_dist):
        ref_moment = ref_moments[key]
        cand_moment = cand_moments[key]
        for field in ("phase_fraction", "x_min_m", "x_max_m"):
            if not math.isclose(float(ref_moment[field]), float(cand_moment[field]),
                                rel_tol=1e-12, abs_tol=1e-14):
                raise ValueError(f"{field} differs for phase/region {key}")
        ref_mean = float(ref_moment["mean_energy_eV"])
        cand_mean = float(cand_moment["mean_energy_eV"])
        ref_temperature = float(ref_moment["drift_separated_temperature_eV"])
        cand_temperature = float(cand_moment["drift_separated_temperature_eV"])
        tails = {}
        for threshold in args.tail_eV:
            ref_tail = tail(ref_dist[key], threshold)
            cand_tail = tail(cand_dist[key], threshold)
            tails[str(threshold)] = {
                "reference": ref_tail, "candidate": cand_tail,
                "absolute_difference": abs(ref_tail - cand_tail),
            }
        pairs.append({
            "phase_bin": key[0], "region": key[1],
            "total_variation": tv_distance(ref_dist[key], cand_dist[key]),
            "mean_energy_relative_difference": relative(ref_mean, cand_mean),
            "temperature_relative_difference": relative(ref_temperature,
                                                        cand_temperature),
            "overflow_absolute_difference": abs(
                float(ref_moment["overflow_fraction"]) -
                float(cand_moment["overflow_fraction"])),
            "tail_fraction": tails,
        })
    maximum = {
        "total_variation": max(pair["total_variation"] for pair in pairs),
        "mean_energy_relative_difference": max(
            pair["mean_energy_relative_difference"] for pair in pairs),
        "temperature_relative_difference": max(
            pair["temperature_relative_difference"] for pair in pairs),
        "overflow_absolute_difference": max(
            pair["overflow_absolute_difference"] for pair in pairs),
    }
    thresholds = {
        "maximum_total_variation": args.max_tv,
        "maximum_mean_energy_relative_difference": args.max_mean_energy_relative,
        "maximum_temperature_relative_difference": args.max_temperature_relative,
    }
    declared = all(value is not None for value in thresholds.values())
    passes = None
    if declared:
        passes = (maximum["total_variation"] <= args.max_tv and
                  maximum["mean_energy_relative_difference"] <=
                  args.max_mean_energy_relative and
                  maximum["temperature_relative_difference"] <=
                  args.max_temperature_relative)
    report = {
        "schema_version": 1, "case_id": ref_manifest["case_id"],
        "species": ref_manifest["species"],
        "reference_code": ref_manifest["code"],
        "candidate_code": cand_manifest["code"],
        "phase_region_pairs": len(pairs),
        "metrics": {"total_variation": "piecewise-uniform histogram TV",
                    "relative_difference_denominator":
                        "max(abs(reference), abs(candidate))"},
        "maximum": maximum,
        "acceptance": {"declared": declared, "thresholds": thresholds,
                       "passes": passes},
        "pairs": pairs,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 1 if passes is False else 0


if __name__ == "__main__":
    raise SystemExit(main())
