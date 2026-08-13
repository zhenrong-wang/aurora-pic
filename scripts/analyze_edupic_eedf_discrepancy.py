#!/usr/bin/env python3
"""Localize a fresh AuroraPIC/eduPIC EEDF discrepancy by energy band."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from compare_aurorapic_edupic_measurement_pilot import rows
from run_aurorapic_edupic_pilot import atomic_json, sha256


REFERENCE_EEPF_SHA256 = (
    "8ad8222107c7e61ea9d1a2d36b0550bcd69ded55dbf5fe0bcfd10810dfc56260")
BANDS_EV = ((0.0, 2.0), (2.0, 5.0), (5.0, 15.76),
            (15.76, 50.0), (50.0, 500.0))


def normalize(table: list[tuple[float, float, float]]) -> list[tuple[float, float, float]]:
    total = math.fsum(mass for _, _, mass in table)
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("distribution mass is not positive and finite")
    return [(low, high, mass / total) for low, high, mass in table]


def integrate(table: list[tuple[float, float, float]], low: float,
              high: float, moment: int = 0) -> float:
    total = 0.0
    for bin_low, bin_high, mass in table:
        overlap_low = max(low, bin_low)
        overlap_high = min(high, bin_high)
        if overlap_high <= overlap_low:
            continue
        density = mass / (bin_high - bin_low)
        total += (density * (overlap_high - overlap_low) if moment == 0 else
                  0.5 * density *
                  (overlap_high * overlap_high - overlap_low * overlap_low))
    return total


def candidate_distribution(output: Path) -> list[tuple[float, float, float]]:
    histogram = rows(output / "phase_eedf.csv")
    moments = rows(output / "phase_eedf_moments.csv")
    observations = math.fsum(
        float(row["represented_observations"]) for row in moments)
    counts: dict[int, float] = {}
    centers: dict[int, float] = {}
    for row in histogram:
        index = int(row["energy_bin"])
        counts[index] = counts.get(index, 0.0) + float(
            row["represented_count"])
        centers[index] = float(row["energy_eV"])
    indices = sorted(counts)
    if len(indices) < 2 or observations <= 0.0:
        raise ValueError("candidate EEDF is incomplete")
    width = centers[indices[1]] - centers[indices[0]]
    return normalize([
        (centers[index] - 0.5 * width, centers[index] + 0.5 * width,
         counts[index] / observations) for index in indices])


def reference_distribution(path: Path) -> list[tuple[float, float, float]]:
    if sha256(path) != REFERENCE_EEPF_SHA256:
        raise ValueError("locked eduPIC EEPF differs")
    return normalize([
        (float(row["energy_ev"]) - 0.025,
         float(row["energy_ev"]) + 0.025,
         float(row["equal_time_block_mixture_eepf_ev-1p5"]) *
         math.sqrt(float(row["energy_ev"])) * 0.05)
        for row in rows(path)])


def analyze(output: Path, reference: Path) -> dict[str, object]:
    candidate = candidate_distribution(output)
    locked = reference_distribution(reference)
    candidate_mean = integrate(candidate, 0.0, 500.0, 1)
    reference_mean = integrate(locked, 0.0, 500.0, 1)
    bands = []
    for low, high in BANDS_EV:
        candidate_mass = integrate(candidate, low, high)
        reference_mass = integrate(locked, low, high)
        candidate_energy = integrate(candidate, low, high, 1)
        reference_energy = integrate(locked, low, high, 1)
        bands.append({
            "low_eV": low,
            "high_eV": high,
            "candidate_probability_mass": candidate_mass,
            "reference_probability_mass": reference_mass,
            "probability_mass_difference": candidate_mass - reference_mass,
            "candidate_mean_energy_contribution_eV": candidate_energy,
            "reference_mean_energy_contribution_eV": reference_energy,
            "mean_energy_excess_contribution_eV":
                candidate_energy - reference_energy,
        })
    excess = candidate_mean - reference_mean
    dominant = max(bands, key=lambda item: abs(
        float(item["mean_energy_excess_contribution_eV"])))
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "post_measurement_eedf_energy_band_localization",
        "candidate_phase_eedf_sha256": sha256(output / "phase_eedf.csv"),
        "reference_eepf_sha256": REFERENCE_EEPF_SHA256,
        "candidate_mean_energy_eV": candidate_mean,
        "reference_mean_energy_eV": reference_mean,
        "mean_energy_excess_eV": excess,
        "energy_bands": bands,
        "dominant_absolute_excess_band_eV": [
            dominant["low_eV"], dominant["high_eV"]],
        "sub_ionization_threshold_excess_fraction": math.fsum(
            float(item["mean_energy_excess_contribution_eV"])
            for item in bands if float(item["high_eV"]) <= 15.76) / excess,
        "interpretation": (
            "The band decomposition distinguishes a bulk-shape discrepancy "
            "from an ionizing-tail-only discrepancy; it does not identify a "
            "specific collision or transport implementation defect."),
        "claim_boundary": (
            "This is a post-measurement diagnostic with no prospective "
            "acceptance threshold."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_output", type=Path)
    parser.add_argument("reference_eepf", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = analyze(args.candidate_output.resolve(),
                     args.reference_eepf.resolve())
    if args.output:
        atomic_json(args.output, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
