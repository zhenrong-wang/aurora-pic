#!/usr/bin/env python3
"""Evaluate the locked common-state 2:1 fixed-particle mesh refinement."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from compare_aurorapic_measurement_windows import (
    average_ionization_rate, density_integral, eedf, ion_impact_distribution,
)
from compare_aurorapic_edupic_measurement_pilot import (
    distribution_mean, relative_difference, relative_l2, rows, total_variation,
    trapezoid,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "c5397fc5a8c129dadaaa91d1134716a20caf53a115fb9cbf138ad7f84480bc57")
BASELINE_REPORT_SHA256 = (
    "d2024ff5b4aa82960bd24a5a38c2ee117f9c01694c32b28ecbe73dee489774c6")
ELEMENTARY_CHARGE_C = 1.60217662e-19
PHASES = 16


def branch_report(output: Path, refined: bool) -> tuple[Path, dict[str, object]]:
    path = output.parent.parent / "branch-report.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("all_gates_passed") is not True:
        raise ValueError(f"branch did not pass solver gates: {path}")
    if refined:
        valid = (value.get("scope") == "common_state_numerical_refinement_branch" and
                 value.get("axis") == "mesh_2x_fixed_particles" and
                 value.get("branch") == "refined_grid" and
                 value.get("rule_sha256") == RULE_SHA256)
    else:
        valid = sha256(path) == BASELINE_REPORT_SHA256
    if not valid:
        raise ValueError(f"branch identity differs: {path}")
    for name, expected in value["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"branch output differs: {output / name}")
    return path, value


def selected_species(output: Path, species: str) -> list[dict[str, str]]:
    return [row for row in rows(output / "spatial_phase_moments.csv")
            if row["species"] == species]


def coincident(refined: list[float], baseline_nodes: int,
               refined_nodes: int) -> list[float]:
    if refined_nodes != 2 * baseline_nodes - 1:
        raise ValueError("refined grid is not an exact 2:1 nodal refinement")
    if len(refined) != PHASES * refined_nodes:
        raise ValueError("refined phase-space vector has unexpected size")
    return [refined[phase * refined_nodes + 2 * node]
            for phase in range(PHASES) for node in range(baseline_nodes)]


def spatial_profile(output: Path, species: str) -> tuple[list[float], list[float]]:
    selected = [row for row in rows(output / "spatial_average.csv")
                if row["species"] == species]
    return ([float(row["x_m"]) for row in selected],
            [float(row["number_density_mean_m-3"]) for row in selected])


def phase_space_average(x: list[float], values: list[float]) -> float:
    nodes = len(x)
    if len(values) != PHASES * nodes:
        raise ValueError("phase-space average shape differs")
    return math.fsum(
        trapezoid(x, values[phase * nodes:(phase + 1) * nodes]) / (x[-1] - x[0])
        for phase in range(PHASES)) / PHASES


def electron_power_per_particle(output: Path) -> float:
    fields = rows(output / "spatial_phase_fields.csv")
    electrons = selected_species(output, "electrons")
    if len(fields) != len(electrons):
        raise ValueError("field and electron phase-space shapes differ")
    x = [float(row["x_m"]) for row in fields[:len(fields) // PHASES]]
    density = [float(row["number_density_mean_m-3"]) for row in electrons]
    current = [-ELEMENTARY_CHARGE_C * number * float(row["mean_velocity_x"])
               for number, row in zip(density, electrons)]
    power = [value * float(row["electric_field_mean_V_m"])
             for value, row in zip(current, fields)]
    return phase_space_average(x, power) / phase_space_average(x, density)


def effective_ionization_frequency(output: Path) -> float:
    return average_ionization_rate(output) / (
        density_integral(output, "electrons") / 0.025)


def evaluate(metrics: dict[str, float], limits: dict[str, object]) -> dict[str, bool]:
    low, high = limits["allowed_average_ionization_rate_ratio"]
    mapping = {
        "electron_density_integral": ("electron_density_integral_relative_change",
            "maximum_electron_density_integral_relative_change"),
        "ion_density_integral": ("ion_density_integral_relative_change",
            "maximum_ion_density_integral_relative_change"),
        "electron_density_profile": ("electron_density_profile_relative_l2",
            "maximum_electron_density_profile_relative_l2"),
        "ion_density_profile": ("ion_density_profile_relative_l2",
            "maximum_ion_density_profile_relative_l2"),
        "electron_energy_distribution": ("electron_energy_distribution_total_variation",
            "maximum_electron_energy_distribution_total_variation"),
        "electron_mean_energy": ("electron_mean_energy_relative_change",
            "maximum_electron_mean_energy_relative_change"),
        "electric_field_phase_space": ("electric_field_phase_space_relative_l2",
            "maximum_electric_field_phase_space_relative_l2"),
        "electron_current_phase_space": ("electron_current_phase_space_relative_l2",
            "maximum_electron_current_phase_space_relative_l2"),
        "powered_ion_energy_distribution": ("powered_ion_energy_distribution_total_variation",
            "maximum_powered_ion_energy_distribution_total_variation"),
        "grounded_ion_energy_distribution": ("grounded_ion_energy_distribution_total_variation",
            "maximum_grounded_ion_energy_distribution_total_variation"),
        "electrode_mean_ion_energy": ("maximum_electrode_mean_ion_energy_relative_change",
            "maximum_electrode_mean_ion_energy_relative_change"),
        "electron_power_per_particle": ("electron_power_per_particle_relative_change",
            "maximum_electron_power_per_particle_relative_change"),
        "ionization_frequency": ("ionization_frequency_relative_change",
            "maximum_ionization_frequency_relative_change"),
    }
    result = {gate: metrics[metric] <= limits[limit]
              for gate, (metric, limit) in mapping.items()}
    result["average_ionization_rate"] = (
        low <= metrics["average_ionization_rate_ratio"] <= high)
    return result


def analyze(baseline: Path, refined: Path, rule_path: Path) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("mesh refinement rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    baseline_path, _ = branch_report(baseline, False)
    refined_path, _ = branch_report(refined, True)
    baseline_nodes = int(rule["paired_baseline"]["nodes"])
    refined_nodes = int(rule["branches"]["refined_grid"]["nodes"])
    metrics: dict[str, float] = {}
    for species in ("electrons", "ions"):
        first_x, first_profile = spatial_profile(baseline, species)
        second_x, second_profile = spatial_profile(refined, species)
        if any(abs(second_x[2 * i] - x) > 1e-15 for i, x in enumerate(first_x)):
            raise ValueError("refined and baseline nodes do not coincide")
        metrics[f"{species[:-1]}_density_integral_relative_change"] = (
            relative_difference(trapezoid(second_x, second_profile),
                                trapezoid(first_x, first_profile)))
        metrics[f"{species[:-1]}_density_profile_relative_l2"] = relative_l2(
            second_profile[::2], first_profile)
    first_eedf, second_eedf = eedf(baseline), eedf(refined)
    metrics["electron_energy_distribution_total_variation"] = total_variation(
        second_eedf, first_eedf)
    metrics["electron_mean_energy_relative_change"] = relative_difference(
        distribution_mean(second_eedf), distribution_mean(first_eedf))
    first_fields = [float(row["electric_field_mean_V_m"])
                    for row in rows(baseline / "spatial_phase_fields.csv")]
    second_fields = [float(row["electric_field_mean_V_m"])
                     for row in rows(refined / "spatial_phase_fields.csv")]
    metrics["electric_field_phase_space_relative_l2"] = relative_l2(
        coincident(second_fields, baseline_nodes, refined_nodes), first_fields)
    first_e = selected_species(baseline, "electrons")
    second_e = selected_species(refined, "electrons")
    first_current = [-ELEMENTARY_CHARGE_C * float(row["number_density_mean_m-3"]) *
                     float(row["mean_velocity_x"]) for row in first_e]
    second_current = [-ELEMENTARY_CHARGE_C * float(row["number_density_mean_m-3"]) *
                      float(row["mean_velocity_x"]) for row in second_e]
    metrics["electron_current_phase_space_relative_l2"] = relative_l2(
        coincident(second_current, baseline_nodes, refined_nodes), first_current)
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
    passed = all(gates.values())
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "prospective_common_state_mesh_refinement_result",
        "rule_sha256": RULE_SHA256,
        "baseline_branch_report_sha256": sha256(baseline_path),
        "refined_branch_report_sha256": sha256(refined_path),
        "grid_refinement_ratio": 2,
        "fixed_total_particles": True,
        "metrics": metrics,
        "thresholds": limits,
        "gates": gates,
        "all_gates_passed": passed,
        "interpretation": rule["interpretation"]["pass" if passed else "fail"],
        "claim_boundary": (
            "This isolates a 2:1 field-grid change at fixed total particles. "
            "Passing does not establish same-particles-per-cell or asymptotic "
            "mesh convergence; failure may combine grid and particle-noise effects."),
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
