#!/usr/bin/env python3
"""Analyze the predeclared constrained initial-microstate ensemble."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics

from analyze_aurorapic_heating_localization import NODES, PHASES, reduced_profiles
from analyze_aurorapic_heating_seed_ensemble import relative_range, vector_mean
from compare_aurorapic_edupic_measurement_pilot import rows
from compare_edupic_phase_space import (
    ELEMENTARY_CHARGE_C, REFERENCE_FILES, flatten_phase_major, metrics,
    read_matrix, spatial_phase_average,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "53f372d3ad10fb6c24c5265c962c559ca92d94cd4a3008c82c320b88912aabb4")
BASELINE_REPORT_SHA256 = (
    "4171678276ffb9fee0eecb843f9d9c0b44ecb5bbf13c65207618a9cd76820fc5")
LENGTH_M = 0.025


def load_member(name: str, root: Path, rule: dict[str, object],
                baseline: bool = False) -> dict[str, object]:
    report_path = root / "branch-report.json"
    report_hash = sha256(report_path)
    if baseline and report_hash != BASELINE_REPORT_SHA256:
        raise ValueError("locked source-microstate branch report differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("all_gates_passed") is not True:
        raise ValueError(f"{name} branch did not pass")
    if not baseline:
        if report.get("rule_sha256") != RULE_SHA256 or report.get("branch") != name:
            raise ValueError(f"{name} branch contract differs")
        expected = rule["particle_states"][name]
        if (report["inputs"].get("particle_state_name") != name or
                report["inputs"].get("particle_state_sha256") !=
                expected["particle_state_sha256"]):
            raise ValueError(f"{name} particle state differs")
    output = root / "measurement" / "output"
    for filename, expected in report["output_hashes"].items():
        if sha256(output / filename) != expected:
            raise ValueError(f"{name} output differs: {filename}")
    metadata_path = output / "spatial_average_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if (metadata.get("sampling_order") != "pre_collision" or
            metadata.get("phase_bins") != PHASES or
            metadata.get("complete") is not True or
            set(metadata.get("phase_bin_samples", [])) != {80}):
        raise ValueError(f"{name} sampling protocol differs")
    fields = rows(output / "spatial_phase_fields.csv")
    moments = [row for row in rows(output / "spatial_phase_moments.csv")
               if row["species"] == "electrons"]
    if len(fields) != PHASES * NODES or len(moments) != len(fields):
        raise ValueError(f"{name} phase-space shape differs")
    density = [float(row["number_density_mean_m-3"]) for row in moments]
    field = [float(row["electric_field_mean_V_m"]) for row in fields]
    current = [-ELEMENTARY_CHARGE_C * number * float(row["mean_velocity_x"])
               for number, row in zip(density, moments)]
    power = [j * electric for j, electric in zip(current, field)]
    for vector in (density, field, current, power):
        if any(not math.isfinite(value) for value in vector):
            raise ValueError(f"{name} contains a non-finite phase-space value")
    energy_path = output / "energy-budget.json"
    phase_eedf_path = output / "phase-eedf-analysis.json"
    if (sha256(energy_path) != report["energy_analysis_sha256"] or
            sha256(phase_eedf_path) != report["phase_eedf_analysis_sha256"]):
        raise ValueError(f"{name} analyzer output differs")
    energy = json.loads(energy_path.read_text())
    exact_power = float(
        energy["electric_power_W_m-2"]["electric_work_electrons_J_m-2"]
    ) / LENGTH_M
    return {
        "name": name,
        "branch_report_sha256": report_hash,
        "particle_state_sha256": report["inputs"]["particle_state_sha256"],
        "spatial_average_metadata_sha256": sha256(metadata_path),
        "density": density,
        "field": field,
        "current": current,
        "power": power,
        "average_density": spatial_phase_average(density, PHASES, NODES),
        "average_power": spatial_phase_average(power, PHASES, NODES),
        "exact_power": exact_power,
        "peak_resident_set_kib": report["resources"][
            "maximum_peak_resident_set_kib"],
    }


def profile_scatter(members: list[dict[str, object]], key: str) -> tuple[
        list[float], list[float], list[float]]:
    ensemble = vector_mean([member[key] for member in members])
    _ensemble_phase, ensemble_space = reduced_profiles(ensemble)
    scatter = []
    for member in members:
        _phase, space = reduced_profiles(member[key])
        scatter.append(metrics(space, ensemble_space)["relative_l2"])
    return ensemble, ensemble_space, scatter


def analyze(rule_path: Path, baseline_root: Path,
            branch_roots: dict[str, Path], reference: Path) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("microstate ensemble rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    expected_names = set(rule["particle_states"])
    if set(branch_roots) != expected_names:
        raise ValueError("microstate branch set differs")
    members = [load_member("locked_source_microstate", baseline_root, rule, True)]
    members.extend(load_member(name, branch_roots[name], rule)
                   for name in sorted(branch_roots))

    reference_vectors = {}
    for key, reference_name in (
            ("density", "electron_density"),
            ("field", "electric_field"),
            ("current", "electron_current_density"),
            ("power", "electron_ohmic_power_density")):
        filename, expected = REFERENCE_FILES[reference_name]
        path = reference / filename
        if sha256(path) != expected:
            raise ValueError(f"locked eduPIC reference differs: {filename}")
        reference_vectors[key] = flatten_phase_major(read_matrix(path))

    powers = [float(member["average_power"]) for member in members]
    densities = [float(member["average_density"]) for member in members]
    power_range = relative_range(powers)
    density_range = relative_range(densities)
    ensemble_vectors = {}
    ensemble_spaces = {}
    scatters = {}
    for key in ("density", "field", "current", "power"):
        ensemble_vectors[key], ensemble_spaces[key], scatters[key] = \
            profile_scatter(members, key)
    limits = rule["prospective_internal_repeatability"]
    gates = {
        "volume_phase_power_density_repeatability": power_range <= float(
            limits["maximum_relative_range_volume_phase_power_density"]),
        "volume_phase_electron_density_repeatability": density_range <= float(
            limits["maximum_relative_range_volume_phase_electron_density"]),
        "cycle_average_spatial_current_repeatability": max(
            scatters["current"]) <= float(limits[
                "maximum_member_to_ensemble_cycle_average_spatial_current_relative_l2"]),
    }
    reference_power = spatial_phase_average(
        reference_vectors["power"], PHASES, NODES)
    reference_density = spatial_phase_average(
        reference_vectors["density"], PHASES, NODES)
    member_results = []
    for index, member in enumerate(members):
        power_ratio = float(member["average_power"]) / reference_power
        density_ratio = float(member["average_density"]) / reference_density
        item = {
            "name": member["name"],
            "branch_report_sha256": member["branch_report_sha256"],
            "particle_state_sha256": member["particle_state_sha256"],
            "spatial_average_metadata_sha256":
                member["spatial_average_metadata_sha256"],
            "power_density_W_m-3": member["average_power"],
            "electron_density_m-3": member["average_density"],
            "exact_electric_work_W_m-3": member["exact_power"],
            "phase_binned_to_exact_power_relative_difference": abs(
                float(member["average_power"]) / float(member["exact_power"]) - 1.0),
            "candidate_to_reference_power_density_ratio": power_ratio,
            "candidate_to_reference_electron_density_ratio": density_ratio,
            "candidate_to_reference_power_per_electron_ratio":
                power_ratio / density_ratio,
            "peak_resident_set_kib": member["peak_resident_set_kib"],
            "member_to_ensemble_cycle_average_spatial_relative_l2": {
                key: scatters[key][index] for key in scatters},
        }
        member_results.append(item)

    cross_code = {}
    for key in ensemble_vectors:
        reference_phase, reference_space = reduced_profiles(
            reference_vectors[key])
        ensemble_phase, _space = reduced_profiles(ensemble_vectors[key])
        cross_code[key] = {
            "phase_space": metrics(
                ensemble_vectors[key], reference_vectors[key]),
            "spatially_integrated_phase": metrics(
                ensemble_phase, reference_phase),
            "cycle_average_spatial": metrics(
                ensemble_spaces[key], reference_space),
        }
    mean_power_per_electron_ratio = statistics.fmean(
        item["candidate_to_reference_power_per_electron_ratio"]
        for item in member_results)
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "constrained_initial_microstate_ensemble_analysis",
        "rule_sha256": RULE_SHA256,
        "members": member_results,
        "prospective_repeatability": {
            "power_density_relative_range": power_range,
            "maximum_power_density_relative_range": limits[
                "maximum_relative_range_volume_phase_power_density"],
            "electron_density_relative_range": density_range,
            "maximum_electron_density_relative_range": limits[
                "maximum_relative_range_volume_phase_electron_density"],
            "maximum_member_to_ensemble_cycle_average_spatial_current_relative_l2":
                max(scatters["current"]),
            "spatial_current_relative_l2_limit": limits[
                "maximum_member_to_ensemble_cycle_average_spatial_current_relative_l2"],
            "gates": gates,
            "all_gates_passed": all(gates.values()),
        },
        "ensemble_cross_code": {
            "mean_power_density_W_m-3": statistics.fmean(powers),
            "mean_electron_density_m-3": statistics.fmean(densities),
            "mean_power_per_electron_ratio": mean_power_per_electron_ratio,
            "comparisons": cross_code,
        },
        "microstate_scatter": {
            key: {
                "minimum_member_to_ensemble_cycle_average_spatial_relative_l2":
                    min(values),
                "maximum_member_to_ensemble_cycle_average_spatial_relative_l2":
                    max(values),
            } for key, values in scatters.items()
        },
        "assessment": {
            "constrained_microstate_repeatability_reached": all(gates.values()),
            "density_mismatch_exceeds_microstate_scatter":
                cross_code["density"]["cycle_average_spatial"]["relative_l2"] >
                max(scatters["density"]),
            "spatial_current_mismatch_exceeds_microstate_scatter":
                cross_code["current"]["cycle_average_spatial"]["relative_l2"] >
                max(scatters["current"]),
            "power_per_electron_deficit_fraction":
                1.0 - mean_power_per_electron_ratio,
            "finding": (
                "The density, power-per-electron, and cycle-average spatial-"
                "current discrepancies survive constrained independent "
                "microstates and materially exceed their observed scatter."),
            "next_target": (
                "Prospectively compare electron current per represented "
                "particle and source-loss balance against eduPIC to separate "
                "transport/heating differences from population regulation."),
        },
        "claim_boundary": (
            "This quantifies sensitivity to constrained conditional initial "
            "microstates that share cell populations, CIC first moments, and "
            "cellwise empirical velocity tuples. It is not experimental "
            "validation or unrestricted initial-condition uncertainty."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("baseline_root", type=Path)
    parser.add_argument("microstate_51949_root", type=Path)
    parser.add_argument("microstate_63059_root", type=Path)
    parser.add_argument("reference_raw_data", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(
        args.rule.resolve(), args.baseline_root.resolve(), {
            "microstate_51949": args.microstate_51949_root.resolve(),
            "microstate_63059": args.microstate_63059_root.resolve(),
        }, args.reference_raw_data.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
