#!/usr/bin/env python3
"""Evaluate prospective repeatability gates for two fresh PIC windows."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import random

from compare_aurorapic_edupic_measurement_pilot import (
    distribution_mean, relative_difference, relative_l2, rows, total_variation,
    trapezoid,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


APPROVED_RULE_SHA256 = (
    "f9a5e33683986432f3c2050515ad6e7de02316b14dd35c14d3e6a05694e5a216")
ELEMENTARY_CHARGE_C = 1.60217662e-19


def passing_report(output: Path) -> tuple[Path, dict[str, object]]:
    path = output.parent / "measurement-report.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    if (report.get("scope") != "fresh_window_aurorapic_measurement_pilot" or
            report.get("all_gates_passed") is not True or
            report.get("window", {}).get("equilibration_statistics_excluded")
            is not True):
        raise ValueError(f"not a passing fresh measurement: {path}")
    for name, expected in report.get("output_hashes", {}).items():
        if sha256(output / name) != expected:
            raise ValueError(f"measurement output differs: {output / name}")
    return path, report


def density_integral(output: Path, species: str) -> float:
    selected = [row for row in rows(output / "spatial_average.csv")
                if row["species"] == species]
    return trapezoid(
        [float(row["x_m"]) for row in selected],
        [float(row["number_density_mean_m-3"]) for row in selected])


def eedf(output: Path) -> list[tuple[float, float, float]]:
    histogram = rows(output / "phase_eedf.csv")
    observations = math.fsum(float(row["represented_observations"])
                             for row in rows(output / "phase_eedf_moments.csv"))
    counts: dict[int, float] = {}
    centers: dict[int, float] = {}
    for row in histogram:
        index = int(row["energy_bin"])
        counts[index] = counts.get(index, 0.0) + float(row["represented_count"])
        centers[index] = float(row["energy_eV"])
    indices = sorted(counts)
    if len(indices) < 2 or observations <= 0.0:
        raise ValueError("incomplete phase EEDF")
    width = centers[indices[1]] - centers[indices[0]]
    return [(centers[index] - 0.5 * width,
             centers[index] + 0.5 * width,
             counts[index] / observations) for index in indices]


def ion_impact_distribution(
    output: Path, side: str,
) -> list[tuple[float, float, float]]:
    selected = [row for row in rows(output / "wall_impact_spectrum.csv")
                if row["species"] == "ions" and row["electrode"] == side]
    if len(selected) < 2:
        raise ValueError(f"incomplete ion impact distribution at {side}")
    width = (float(selected[1]["impact_energy_eV"]) -
             float(selected[0]["impact_energy_eV"]))
    return [(float(row["impact_energy_eV"]) - 0.5 * width,
             float(row["impact_energy_eV"]) + 0.5 * width,
             float(row["probability_density"]) * width)
            for row in selected]


def ion_impact_counts(output: Path, side: str) -> list[int]:
    return [int(row["macro_count"])
            for row in rows(output / "wall_impact_spectrum.csv")
            if row["species"] == "ions" and row["electrode"] == side]


def multinomial_tv_null(
    first: list[int], second: list[int], seed: int, trials: int = 3000,
) -> dict[str, object]:
    """Condition on pooled bins and estimate finite-count TV under one law."""
    if len(first) != len(second) or not first:
        raise ValueError("impact count tables differ")
    first_total, second_total = sum(first), sum(second)
    if first_total <= 0 or second_total <= 0:
        raise ValueError("impact count totals must be positive")
    weights = [a + b for a, b in zip(first, second)]
    population = list(range(len(weights)))
    observed = 0.5 * math.fsum(
        abs(a / first_total - b / second_total)
        for a, b in zip(first, second))
    rng = random.Random(seed)
    samples = []
    for _ in range(trials):
        a = [0] * len(weights)
        b = [0] * len(weights)
        for index in rng.choices(population, weights=weights, k=first_total):
            a[index] += 1
        for index in rng.choices(population, weights=weights, k=second_total):
            b[index] += 1
        samples.append(0.5 * math.fsum(
            abs(x / first_total - y / second_total) for x, y in zip(a, b)))
    samples.sort()
    quantile = lambda fraction: samples[int(fraction * (trials - 1))]
    return {
        "baseline_macro_impacts": first_total,
        "replication_macro_impacts": second_total,
        "observed_total_variation": observed,
        "null_mean_total_variation": math.fsum(samples) / trials,
        "null_95_percent_interval": [quantile(0.025), quantile(0.975)],
        "upper_tail_monte_carlo_probability":
            (1 + sum(value >= observed for value in samples)) / (trials + 1),
        "trials": trials,
        "seed": seed,
        "interpretation_boundary": (
            "This post-hoc conditional multinomial diagnostic assumes "
            "independent impacts. It can identify sparse-bin sampling as a "
            "plausible explanation but cannot override a prospective gate."),
    }


def ordered_values(path: Path, field: str, species: str | None = None,
                   current: bool = False) -> list[float]:
    selected = rows(path)
    if species is not None:
        selected = [row for row in selected if row["species"] == species]
    keys = [(int(row["phase_bin"]), int(row["node"])) for row in selected]
    if keys != sorted(keys) or len(keys) != len(set(keys)):
        raise ValueError(f"phase-space rows are not uniquely ordered: {path}")
    if current:
        charge = -ELEMENTARY_CHARGE_C if species == "electrons" else ELEMENTARY_CHARGE_C
        return [charge * float(row["number_density_mean_m-3"]) *
                float(row["mean_velocity_x"]) for row in selected]
    return [float(row[field]) for row in selected]


def average_ionization_rate(output: Path) -> float:
    selected = [row for row in rows(output / "spatial_collision_rate.csv")
                if row["channel"] == "electron_mcc.ionization"]
    x = [float(row["x_m"]) for row in selected]
    rate = [float(row["mean_event_rate_m-3_s-1"]) for row in selected]
    return trapezoid(x, rate) / (x[-1] - x[0])


def evaluate(metrics: dict[str, float], thresholds: dict[str, object]) -> dict[str, bool]:
    low, high = thresholds["allowed_average_ionization_rate_ratio"]
    return {
        "electron_density_integral": metrics["electron_density_integral_relative_change"] <=
            thresholds["maximum_electron_density_integral_relative_change"],
        "ion_density_integral": metrics["ion_density_integral_relative_change"] <=
            thresholds["maximum_ion_density_integral_relative_change"],
        "electron_energy_distribution": metrics["electron_energy_distribution_total_variation"] <=
            thresholds["maximum_electron_energy_distribution_total_variation"],
        "electron_mean_energy": metrics["electron_mean_energy_relative_change"] <=
            thresholds["maximum_electron_mean_energy_relative_change"],
        "powered_ion_energy_distribution": metrics["powered_ion_energy_distribution_total_variation"] <=
            thresholds["maximum_powered_ion_energy_distribution_total_variation"],
        "grounded_ion_energy_distribution": metrics["grounded_ion_energy_distribution_total_variation"] <=
            thresholds["maximum_grounded_ion_energy_distribution_total_variation"],
        "electrode_mean_ion_energy": metrics["maximum_electrode_mean_ion_energy_relative_change"] <=
            thresholds["maximum_electrode_mean_ion_energy_relative_change"],
        "electric_field_phase_space": metrics["electric_field_phase_space_relative_l2"] <=
            thresholds["maximum_electric_field_phase_space_relative_l2"],
        "electron_current_phase_space": metrics["electron_current_phase_space_relative_l2"] <=
            thresholds["maximum_electron_current_phase_space_relative_l2"],
        "average_ionization_rate": low <= metrics["average_ionization_rate_ratio"] <= high,
    }


def analyze(baseline: Path, replication: Path, rule_path: Path) -> dict[str, object]:
    if sha256(rule_path) != APPROVED_RULE_SHA256:
        raise ValueError("replication rule is not the prospectively approved rule")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    baseline_report_path, baseline_report = passing_report(baseline)
    replication_report_path, replication_report = passing_report(replication)
    expected_windows = rule["prospective_replication_acceptance"]
    if ([baseline_report["window"]["start_cycle"] + 1,
         baseline_report["window"]["end_cycle"]] !=
            expected_windows["reference_window_cycles"] or
            [replication_report["window"]["start_cycle"] + 1,
             replication_report["window"]["end_cycle"]] !=
            expected_windows["replication_window_cycles"]):
        raise ValueError("measurement windows differ from the replication rule")

    metrics: dict[str, float] = {}
    for species in ("electrons", "ions"):
        first = density_integral(baseline, species)
        second = density_integral(replication, species)
        metrics[f"{species[:-1] if species.endswith('s') else species}_density_integral_relative_change"] = relative_difference(second, first)

    first_eedf, second_eedf = eedf(baseline), eedf(replication)
    metrics["electron_energy_distribution_total_variation"] = total_variation(
        second_eedf, first_eedf)
    metrics["electron_mean_energy_relative_change"] = relative_difference(
        distribution_mean(second_eedf), distribution_mean(first_eedf))

    ion_mean_changes = []
    impact_sampling = {}
    for side, electrode, seed in (("left", "powered", 20260813),
                                  ("right", "grounded", 20260814)):
        first_ifed = ion_impact_distribution(baseline, side)
        second_ifed = ion_impact_distribution(replication, side)
        metrics[f"{electrode}_ion_energy_distribution_total_variation"] = (
            total_variation(second_ifed, first_ifed))
        ion_mean_changes.append(relative_difference(
            distribution_mean(second_ifed), distribution_mean(first_ifed)))
        impact_sampling[electrode] = multinomial_tv_null(
            ion_impact_counts(baseline, side),
            ion_impact_counts(replication, side), seed)
    metrics["maximum_electrode_mean_ion_energy_relative_change"] = max(ion_mean_changes)

    metrics["electric_field_phase_space_relative_l2"] = relative_l2(
        ordered_values(replication / "spatial_phase_fields.csv",
                       "electric_field_mean_V_m"),
        ordered_values(baseline / "spatial_phase_fields.csv",
                       "electric_field_mean_V_m"))
    metrics["electron_current_phase_space_relative_l2"] = relative_l2(
        ordered_values(replication / "spatial_phase_moments.csv", "",
                       "electrons", True),
        ordered_values(baseline / "spatial_phase_moments.csv", "",
                       "electrons", True))
    metrics["average_ionization_rate_ratio"] = (
        average_ionization_rate(replication) / average_ionization_rate(baseline))

    thresholds = dict(expected_windows)
    thresholds.pop("reference_window_cycles")
    thresholds.pop("replication_window_cycles")
    thresholds.pop("all_gates_required")
    gates = evaluate(metrics, thresholds)
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "prospective_stationary_measurement_repeatability",
        "rule_sha256": APPROVED_RULE_SHA256,
        "baseline_measurement_report_sha256": sha256(baseline_report_path),
        "replication_measurement_report_sha256": sha256(replication_report_path),
        "metrics": metrics,
        "thresholds": thresholds,
        "gates": gates,
        "all_gates_passed": all(gates.values()),
        "post_hoc_ion_impact_sampling_diagnostic": impact_sampling,
        "claim_boundary": (
            "Passing establishes repeatability of two adjacent four-cycle "
            "candidate windows under the declared gates. It does not establish "
            "independent-seed convergence, cross-code agreement, or physical validation."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_output", type=Path)
    parser.add_argument("replication_output", type=Path)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = analyze(args.baseline_output.resolve(), args.replication_output.resolve(),
                     args.rule.resolve())
    if args.output:
        atomic_json(args.output, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["all_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
