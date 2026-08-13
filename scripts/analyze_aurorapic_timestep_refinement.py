#!/usr/bin/env python3
"""Evaluate the locked common-state 2:1 timestep refinement."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from compare_aurorapic_measurement_windows import (
    average_ionization_rate, density_integral, eedf, ion_impact_distribution,
    ordered_values,
)
from compare_aurorapic_edupic_measurement_pilot import (
    distribution_mean, relative_difference, relative_l2, rows, total_variation,
)
from compare_edupic_phase_space import spatial_phase_average
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "d02b420adb3f457b5f6f471db6260fc127fda2435a83ddf7e307d22ab6640174")
ELEMENTARY_CHARGE_C = 1.60217662e-19


def report(output: Path, expected_branch: str) -> tuple[Path, dict[str, object]]:
    work = output.parent.parent
    path = work / "branch-report.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    if (value.get("scope") != "common_state_timestep_refinement_branch" or
            value.get("branch") != expected_branch or
            value.get("rule_sha256") != RULE_SHA256 or
            value.get("all_gates_passed") is not True):
        raise ValueError(f"invalid branch report: {path}")
    for name, expected in value["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"branch output differs: {output / name}")
    return path, value


def electron_power_per_particle(output: Path) -> float:
    fields = rows(output / "spatial_phase_fields.csv")
    electrons = [row for row in rows(output / "spatial_phase_moments.csv")
                 if row["species"] == "electrons"]
    if len(fields) != len(electrons):
        raise ValueError("field and electron phase-space shapes differ")
    density = [float(row["number_density_mean_m-3"]) for row in electrons]
    current = [-ELEMENTARY_CHARGE_C * number * float(row["mean_velocity_x"])
               for number, row in zip(density, electrons)]
    power = [value * float(row["electric_field_mean_V_m"])
             for value, row in zip(current, fields)]
    return spatial_phase_average(power) / spatial_phase_average(density)


def effective_ionization_frequency(output: Path) -> float:
    rate = average_ionization_rate(output)
    density = density_integral(output, "electrons") / 0.025
    return rate / density


def evaluate(metrics: dict[str, float], limits: dict[str, object]) -> dict[str, bool]:
    low, high = limits["allowed_average_ionization_rate_ratio"]
    return {
        "electron_density_integral": metrics["electron_density_integral_relative_change"] <=
            limits["maximum_electron_density_integral_relative_change"],
        "ion_density_integral": metrics["ion_density_integral_relative_change"] <=
            limits["maximum_ion_density_integral_relative_change"],
        "electron_energy_distribution": metrics["electron_energy_distribution_total_variation"] <=
            limits["maximum_electron_energy_distribution_total_variation"],
        "electron_mean_energy": metrics["electron_mean_energy_relative_change"] <=
            limits["maximum_electron_mean_energy_relative_change"],
        "electric_field_phase_space": metrics["electric_field_phase_space_relative_l2"] <=
            limits["maximum_electric_field_phase_space_relative_l2"],
        "electron_current_phase_space": metrics["electron_current_phase_space_relative_l2"] <=
            limits["maximum_electron_current_phase_space_relative_l2"],
        "average_ionization_rate": low <= metrics["average_ionization_rate_ratio"] <= high,
        "powered_ion_energy_distribution": metrics["powered_ion_energy_distribution_total_variation"] <=
            limits["maximum_powered_ion_energy_distribution_total_variation"],
        "grounded_ion_energy_distribution": metrics["grounded_ion_energy_distribution_total_variation"] <=
            limits["maximum_grounded_ion_energy_distribution_total_variation"],
        "electrode_mean_ion_energy": metrics["maximum_electrode_mean_ion_energy_relative_change"] <=
            limits["maximum_electrode_mean_ion_energy_relative_change"],
        "electron_power_per_particle": metrics["electron_power_per_particle_relative_change"] <=
            limits["maximum_electron_power_per_particle_relative_change"],
        "ionization_frequency": metrics["ionization_frequency_relative_change"] <=
            limits["maximum_ionization_frequency_relative_change"],
    }


def analyze(baseline: Path, refined: Path, rule_path: Path) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("timestep refinement rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    baseline_path, baseline_report = report(baseline, "baseline_dt")
    refined_path, refined_report = report(refined, "half_dt")
    if (float(baseline_report["numerics"]["timestep_s"]) /
            float(refined_report["numerics"]["timestep_s"]) != 2.0):
        raise ValueError("branch timestep ratio is not exactly 2:1")
    metrics: dict[str, float] = {}
    for species in ("electrons", "ions"):
        metrics[f"{species[:-1]}_density_integral_relative_change"] = (
            relative_difference(density_integral(refined, species),
                                density_integral(baseline, species)))
    baseline_eedf, refined_eedf = eedf(baseline), eedf(refined)
    metrics["electron_energy_distribution_total_variation"] = total_variation(
        refined_eedf, baseline_eedf)
    metrics["electron_mean_energy_relative_change"] = relative_difference(
        distribution_mean(refined_eedf), distribution_mean(baseline_eedf))
    metrics["electric_field_phase_space_relative_l2"] = relative_l2(
        ordered_values(refined / "spatial_phase_fields.csv", "electric_field_mean_V_m"),
        ordered_values(baseline / "spatial_phase_fields.csv", "electric_field_mean_V_m"))
    metrics["electron_current_phase_space_relative_l2"] = relative_l2(
        ordered_values(refined / "spatial_phase_moments.csv", "", "electrons", True),
        ordered_values(baseline / "spatial_phase_moments.csv", "", "electrons", True))
    metrics["average_ionization_rate_ratio"] = (
        average_ionization_rate(refined) / average_ionization_rate(baseline))
    rebin = int(rule["fresh_measurement_contract"]["ifed_comparison_rebin_factor"])
    ion_means = []
    for side, electrode in (("left", "powered"), ("right", "grounded")):
        first = ion_impact_distribution(baseline, side, rebin)
        second = ion_impact_distribution(refined, side, rebin)
        metrics[f"{electrode}_ion_energy_distribution_total_variation"] = (
            total_variation(second, first))
        ion_means.append(relative_difference(
            distribution_mean(second), distribution_mean(first)))
    metrics["maximum_electrode_mean_ion_energy_relative_change"] = max(ion_means)
    metrics["electron_power_per_particle_relative_change"] = relative_difference(
        electron_power_per_particle(refined), electron_power_per_particle(baseline))
    metrics["ionization_frequency_relative_change"] = relative_difference(
        effective_ionization_frequency(refined),
        effective_ionization_frequency(baseline))
    limits = dict(rule["prospective_acceptance"])
    limits.pop("all_gates_required")
    gates = evaluate(metrics, limits)
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "prospective_common_state_timestep_refinement_result",
        "rule_sha256": RULE_SHA256,
        "baseline_branch_report_sha256": sha256(baseline_path),
        "refined_branch_report_sha256": sha256(refined_path),
        "timestep_refinement_ratio": 2.0,
        "metrics": metrics,
        "thresholds": limits,
        "gates": gates,
        "all_gates_passed": all(gates.values()),
        "interpretation": (rule["interpretation"]["pass"] if all(gates.values())
                           else rule["interpretation"]["fail"]),
        "claim_boundary": (
            "A passing 2:1 paired refinement argues against material ordinary "
            "timestep error at the declared tolerances; it does not prove an "
            "asymptotic temporal order or full numerical convergence."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_output", type=Path)
    parser.add_argument("refined_output", type=Path)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.baseline_output.resolve(), args.refined_output.resolve(),
                     args.rule.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
