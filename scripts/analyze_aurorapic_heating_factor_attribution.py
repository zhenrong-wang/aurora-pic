#!/usr/bin/env python3
"""Attribute matched heating differences exactly between current and field."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_aurorapic_heating_localization import NODES, PHASES, reduced_profiles
from analyze_aurorapic_heating_seed_ensemble import REPORT_HASHES, vector_mean
from compare_aurorapic_edupic_measurement_pilot import rows
from compare_edupic_phase_space import (
    ELEMENTARY_CHARGE_C, REFERENCE_FILES, flatten_phase_major, metrics,
    read_matrix, spatial_phase_average,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


SPATIAL_BOUNDARIES = [0, 40, 80, 160, 240, 320, 360, 400]


def shapley_current_field(candidate_current: list[float],
                          candidate_field: list[float],
                          reference_current: list[float],
                          reference_field: list[float]
                          ) -> tuple[list[float], list[float]]:
    sizes = {len(candidate_current), len(candidate_field),
             len(reference_current), len(reference_field)}
    if len(sizes) != 1 or not candidate_current:
        raise ValueError("current-field vectors must have one positive size")
    current = [
        0.5 * (candidate_e + reference_e) * (candidate_j - reference_j)
        for candidate_j, candidate_e, reference_j, reference_e in zip(
            candidate_current, candidate_field,
            reference_current, reference_field)
    ]
    field = [
        0.5 * (candidate_j + reference_j) * (candidate_e - reference_e)
        for candidate_j, candidate_e, reference_j, reference_e in zip(
            candidate_current, candidate_field,
            reference_current, reference_field)
    ]
    return current, field


def band_means(values: list[float]) -> list[dict[str, float | int]]:
    if len(values) != PHASES * NODES:
        raise ValueError("factor-attribution phase space has the wrong shape")
    spatial = reduced_profiles(values)[1]
    result = []
    for low, high in zip(SPATIAL_BOUNDARIES, SPATIAL_BOUNDARIES[1:]):
        result.append({
            "first_node": low,
            "last_node_inclusive": high - 1,
            "lower_gap_fraction": low / NODES,
            "upper_gap_fraction": high / NODES,
            "mean_contribution_W_m-3":
                math.fsum(spatial[low:high]) / (high - low),
        })
    return result


def load_member(seed: int, root: Path) -> dict[str, object]:
    report_path = root / "branch-report.json"
    if sha256(report_path) != REPORT_HASHES[seed]:
        raise ValueError(f"seed-{seed} branch report differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("all_gates_passed") is not True:
        raise ValueError(f"seed-{seed} branch did not pass")
    output = root / "measurement" / "output"
    for name, expected in report["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"seed-{seed} output differs: {name}")
    fields = rows(output / "spatial_phase_fields.csv")
    moments = [row for row in rows(output / "spatial_phase_moments.csv")
               if row["species"] == "electrons"]
    if len(fields) != PHASES * NODES or len(moments) != len(fields):
        raise ValueError(f"seed-{seed} phase-space shape differs")
    density = [float(row["number_density_mean_m-3"]) for row in moments]
    current = [
        -ELEMENTARY_CHARGE_C * number * float(row["mean_velocity_x"])
        for number, row in zip(density, moments)
    ]
    field = [float(row["electric_field_mean_V_m"]) for row in fields]
    return {"seed": seed, "density": density,
            "current": current, "field": field}


def analyze(roots: dict[int, Path], reference: Path) -> dict[str, object]:
    members = [load_member(seed, roots[seed]) for seed in sorted(roots)]
    reference_values = {}
    for name in ("electron_density", "electron_current_density",
                 "electric_field", "electron_ohmic_power_density"):
        filename, expected = REFERENCE_FILES[name]
        path = reference / filename
        if sha256(path) != expected:
            raise ValueError(f"locked eduPIC reference differs: {filename}")
        reference_values[name] = flatten_phase_major(read_matrix(path))
    reference_product = [
        current * field for current, field in zip(
            reference_values["electron_current_density"],
            reference_values["electric_field"])
    ]
    reference_factorization = metrics(
        reference_product,
        reference_values["electron_ohmic_power_density"])

    current_contributions = []
    field_contributions = []
    member_results = []
    for member in members:
        current, field = shapley_current_field(
            member["current"], member["field"],
            reference_values["electron_current_density"],
            reference_values["electric_field"])
        candidate_power = [
            value * electric for value, electric in zip(
                member["current"], member["field"])
        ]
        current_average = spatial_phase_average(current, PHASES, NODES)
        field_average = spatial_phase_average(field, PHASES, NODES)
        direct_difference = spatial_phase_average(
            candidate_power, PHASES, NODES) - spatial_phase_average(
                reference_product, PHASES, NODES)
        closure = current_average + field_average - direct_difference
        if abs(closure) > 1.0e-9:
            raise ValueError("current-field attribution does not close")
        current_contributions.append(current)
        field_contributions.append(field)
        member_results.append({
            "seed": member["seed"],
            "power_difference_W_m-3": direct_difference,
            "current_contribution_W_m-3": current_average,
            "field_contribution_W_m-3": field_average,
            "closure_residual_W_m-3": closure,
        })

    ensemble_current_contribution = vector_mean(current_contributions)
    ensemble_field_contribution = vector_mean(field_contributions)
    ensemble_density = vector_mean([member["density"] for member in members])
    ensemble_current = vector_mean([member["current"] for member in members])
    ensemble_current_average = spatial_phase_average(
        ensemble_current_contribution, PHASES, NODES)
    ensemble_field_average = spatial_phase_average(
        ensemble_field_contribution, PHASES, NODES)
    candidate_phase_current, candidate_spatial_current = reduced_profiles(
        ensemble_current)
    reference_phase_current, reference_spatial_current = reduced_profiles(
        reference_values["electron_current_density"])
    candidate_phase_density, candidate_spatial_density = reduced_profiles(
        ensemble_density)
    reference_phase_density, reference_spatial_density = reduced_profiles(
        reference_values["electron_density"])
    candidate_phase_current_per_electron = [
        current / density for current, density in zip(
            candidate_phase_current, candidate_phase_density)
    ]
    reference_phase_current_per_electron = [
        current / density for current, density in zip(
            reference_phase_current, reference_phase_density)
    ]

    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "matched_heating_current_field_factor_attribution",
        "method": {
            "power_identity": "P = J_e E",
            "attribution": "two-factor symmetric Shapley decomposition",
            "pointwise_closure": "delta_P = phi_current + phi_field",
            "density_drift_boundary": (
                "No exact n-u Shapley split is reported because reference "
                "velocity is undefined in empty sheath cells where n=J=0."),
            "acceptance_thresholds_declared": False,
        },
        "reference_power_factorization": reference_factorization,
        "members": member_results,
        "ensemble_volume_phase_attribution": {
            "current_contribution_W_m-3": ensemble_current_average,
            "field_contribution_W_m-3": ensemble_field_average,
            "total_power_difference_W_m-3":
                ensemble_current_average + ensemble_field_average,
            "interpretation": (
                "Current and field effects oppose one another; contribution "
                "fractions are not reported because cancellation makes them "
                "misleading."),
        },
        "ensemble_reduced_profiles": {
            "spatially_integrated_phase_current": metrics(
                candidate_phase_current, reference_phase_current),
            "spatially_integrated_phase_current_per_electron": metrics(
                candidate_phase_current_per_electron,
                reference_phase_current_per_electron),
            "spatially_integrated_phase_density": metrics(
                candidate_phase_density, reference_phase_density),
            "cycle_average_spatial_current": metrics(
                candidate_spatial_current, reference_spatial_current),
            "cycle_average_spatial_density": metrics(
                candidate_spatial_density, reference_spatial_density),
        },
        "spatial_band_attribution": {
            "current": band_means(ensemble_current_contribution),
            "field": band_means(ensemble_field_contribution),
        },
        "finding": (
            "A positive current contribution is over-cancelled by a negative "
            "field contribution. The small net power deficit cannot be assigned "
            "to one factor without losing this cancellation structure."),
        "next_discriminator": (
            "Test whether the persistent cycle-average spatial current and "
            "density amplitudes converge with independent initial particle "
            "realizations or a longer matched measurement window."),
        "claim_boundary": (
            "This exact algebraic attribution identifies correlated factors; "
            "it does not prove causal solver error or experimental accuracy."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_root", type=Path)
    parser.add_argument("seed_24601_root", type=Path)
    parser.add_argument("seed_35713_root", type=Path)
    parser.add_argument("reference_raw_data", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze({
        13507: args.baseline_root.resolve(),
        24601: args.seed_24601_root.resolve(),
        35713: args.seed_35713_root.resolve(),
    }, args.reference_raw_data.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
