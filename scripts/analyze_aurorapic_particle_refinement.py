#!/usr/bin/env python3
"""Evaluate the locked refined-grid 2x macro-particle sensitivity pair."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from analyze_aurorapic_mesh_refinement import (
    ELEMENTARY_CHARGE_C, effective_ionization_frequency,
    electron_power_per_particle, evaluate, selected_species, spatial_profile,
)
from compare_aurorapic_measurement_windows import (
    average_ionization_rate, eedf, ion_impact_distribution,
)
from compare_aurorapic_edupic_measurement_pilot import (
    distribution_mean, relative_difference, relative_l2, rows, total_variation,
    trapezoid,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "b99ab4d96a1f3740bf6df6a084c5f42e45e7de0a58c9b8ae0a003dae7b18a306")
BASELINE_REPORT_SHA256 = (
    "d3348cb694b8c21e5ebad6c8fb4ffc0b1f6203396ede8553a30fddcbb8799d0e")


def branch_report(output: Path, doubled: bool) -> tuple[Path, dict[str, object]]:
    path = output.parent.parent / "branch-report.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("all_gates_passed") is not True:
        raise ValueError(f"branch did not pass solver gates: {path}")
    if doubled:
        valid = (value.get("scope") == "common_state_numerical_refinement_branch" and
                 value.get("axis") == "particles_2x_refined_grid" and
                 value.get("branch") == "double_particles" and
                 value.get("rule_sha256") == RULE_SHA256)
    else:
        valid = sha256(path) == BASELINE_REPORT_SHA256
    if not valid:
        raise ValueError(f"branch identity differs: {path}")
    for name, expected in value["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"branch output differs: {output / name}")
    return path, value


def require_same_grid(first_x: list[float], second_x: list[float]) -> None:
    if len(first_x) != len(second_x) or any(
            abs(first - second) > 1e-15
            for first, second in zip(first_x, second_x)):
        raise ValueError("particle-refinement outputs do not use the same grid")


def analyze(baseline: Path, doubled: Path, rule_path: Path) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("particle refinement rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    baseline_path, _ = branch_report(baseline, False)
    doubled_path, _ = branch_report(doubled, True)
    metrics: dict[str, float] = {}
    for species in ("electrons", "ions"):
        first_x, first_profile = spatial_profile(baseline, species)
        second_x, second_profile = spatial_profile(doubled, species)
        require_same_grid(first_x, second_x)
        metrics[f"{species[:-1]}_density_integral_relative_change"] = (
            relative_difference(trapezoid(second_x, second_profile),
                                trapezoid(first_x, first_profile)))
        metrics[f"{species[:-1]}_density_profile_relative_l2"] = relative_l2(
            second_profile, first_profile)

    first_eedf, second_eedf = eedf(baseline), eedf(doubled)
    metrics["electron_energy_distribution_total_variation"] = total_variation(
        second_eedf, first_eedf)
    metrics["electron_mean_energy_relative_change"] = relative_difference(
        distribution_mean(second_eedf), distribution_mean(first_eedf))

    first_fields = [float(row["electric_field_mean_V_m"])
                    for row in rows(baseline / "spatial_phase_fields.csv")]
    second_fields = [float(row["electric_field_mean_V_m"])
                     for row in rows(doubled / "spatial_phase_fields.csv")]
    metrics["electric_field_phase_space_relative_l2"] = relative_l2(
        second_fields, first_fields)
    first_e = selected_species(baseline, "electrons")
    second_e = selected_species(doubled, "electrons")
    first_current = [-ELEMENTARY_CHARGE_C * float(row["number_density_mean_m-3"]) *
                     float(row["mean_velocity_x"]) for row in first_e]
    second_current = [-ELEMENTARY_CHARGE_C * float(row["number_density_mean_m-3"]) *
                      float(row["mean_velocity_x"]) for row in second_e]
    metrics["electron_current_phase_space_relative_l2"] = relative_l2(
        second_current, first_current)
    metrics["average_ionization_rate_ratio"] = (
        average_ionization_rate(doubled) / average_ionization_rate(baseline))

    rebin = int(rule["fresh_measurement_contract"]["ifed_comparison_rebin_factor"])
    ion_means = []
    for side, electrode in (("left", "powered"), ("right", "grounded")):
        first = ion_impact_distribution(baseline, side, rebin)
        second = ion_impact_distribution(doubled, side, rebin)
        metrics[f"{electrode}_ion_energy_distribution_total_variation"] = (
            total_variation(second, first))
        ion_means.append(relative_difference(
            distribution_mean(second), distribution_mean(first)))
    metrics["maximum_electrode_mean_ion_energy_relative_change"] = max(ion_means)
    metrics["electron_power_per_particle_relative_change"] = relative_difference(
        electron_power_per_particle(doubled), electron_power_per_particle(baseline))
    metrics["ionization_frequency_relative_change"] = relative_difference(
        effective_ionization_frequency(doubled),
        effective_ionization_frequency(baseline))

    limits = dict(rule["prospective_acceptance"])
    limits.pop("all_gates_required")
    gates = evaluate(metrics, limits)
    passed = all(gates.values())
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "prospective_common_state_particle_refinement_result",
        "rule_sha256": RULE_SHA256,
        "baseline_branch_report_sha256": sha256(baseline_path),
        "doubled_branch_report_sha256": sha256(doubled_path),
        "particle_refinement_ratio": 2,
        "fixed_grid": True,
        "metrics": metrics,
        "thresholds": limits,
        "gates": gates,
        "all_gates_passed": passed,
        "interpretation": rule["interpretation"]["pass" if passed else "fail"],
        "claim_boundary": (
            "This isolates a 2x macro-particle count change at fixed 799-node grid "
            "and represented phase space. Initially coincident split children are "
            "collisionally decorrelated for two RF cycles; this is not an "
            "independent-seed ensemble or an asymptotic particle convergence proof."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_output", type=Path)
    parser.add_argument("doubled_output", type=Path)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.baseline_output.resolve(), args.doubled_output.resolve(),
                     args.rule.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
