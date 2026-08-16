#!/usr/bin/env python3
"""Analyze matched-heating continuation-seed repeatability."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics

from analyze_aurorapic_heating_localization import NODES, PHASES, reduced_profiles
from compare_aurorapic_edupic_measurement_pilot import rows
from compare_edupic_phase_space import (
    ELEMENTARY_CHARGE_C, REFERENCE_FILES, flatten_phase_major, metrics,
    read_matrix, spatial_phase_average,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "d44301f76436fcc832b626d9611fbb03d047c9ae6775ab6e94face1e4a01cd49")
REPORT_HASHES = {
    13507: "4171678276ffb9fee0eecb843f9d9c0b44ecb5bbf13c65207618a9cd76820fc5",
    24601: "9803394cd111e0ec2f15bbea9390f0d098beac510a217310b4c97270b697da8e",
    35713: "0ff115551bb527f39048254c5dc860c811cbe2c640ffbe23d2186837f948155f",
}
MAXIMUM_POWER_RELATIVE_RANGE = 0.08
MAXIMUM_DENSITY_RELATIVE_RANGE = 0.03
LENGTH_M = 0.025


def relative_range(values: list[float]) -> float:
    if (not values or any(not math.isfinite(value) or value <= 0.0
                          for value in values)):
        raise ValueError("relative-range values must be positive and finite")
    return (max(values) - min(values)) / statistics.fmean(values)


def vector_mean(vectors: list[list[float]]) -> list[float]:
    if not vectors or not vectors[0] or any(
            len(vector) != len(vectors[0]) for vector in vectors):
        raise ValueError("ensemble vectors must have one common positive size")
    return [math.fsum(vector[index] for vector in vectors) / len(vectors)
            for index in range(len(vectors[0]))]


def load_member(seed: int, root: Path) -> dict[str, object]:
    report_path = root / "branch-report.json"
    if sha256(report_path) != REPORT_HASHES[seed]:
        raise ValueError(f"seed-{seed} branch report differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("all_gates_passed") is not True:
        raise ValueError(f"seed-{seed} branch did not pass")
    if seed != 13507 and report.get("rule_sha256") != RULE_SHA256:
        raise ValueError(f"seed-{seed} ensemble rule differs")
    output = root / "measurement" / "output"
    for name, expected in report["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"seed-{seed} output differs: {name}")
    metadata = json.loads(
        (output / "spatial_average_metadata.json").read_text(encoding="utf-8"))
    if (metadata.get("sampling_order") != "pre_collision" or
            metadata.get("phase_bins") != PHASES or
            metadata.get("complete") is not True):
        raise ValueError(f"seed-{seed} sampling protocol differs")
    fields = rows(output / "spatial_phase_fields.csv")
    moments = [row for row in rows(output / "spatial_phase_moments.csv")
               if row["species"] == "electrons"]
    if len(fields) != PHASES * NODES or len(moments) != len(fields):
        raise ValueError(f"seed-{seed} phase-space shape differs")
    density = [float(row["number_density_mean_m-3"]) for row in moments]
    power = [
        -ELEMENTARY_CHARGE_C * number * float(moment["mean_velocity_x"]) *
        float(field["electric_field_mean_V_m"])
        for number, moment, field in zip(density, moments, fields)
    ]
    energy = json.loads((output / "energy-budget.json").read_text())
    exact_power = float(
        energy["electric_power_W_m-2"]["electric_work_electrons_J_m-2"]
    ) / LENGTH_M
    return {
        "seed": seed,
        "branch_report_sha256": REPORT_HASHES[seed],
        "density": density,
        "power": power,
        "average_density": spatial_phase_average(density, PHASES, NODES),
        "average_power": spatial_phase_average(power, PHASES, NODES),
        "exact_power": exact_power,
        "peak_resident_set_kib": report["resources"][
            "maximum_peak_resident_set_kib"],
    }


def analyze(roots: dict[int, Path], reference: Path) -> dict[str, object]:
    if set(roots) != set(REPORT_HASHES):
        raise ValueError("seed ensemble members differ")
    members = [load_member(seed, roots[seed]) for seed in sorted(roots)]
    reference_power_path = reference / REFERENCE_FILES[
        "electron_ohmic_power_density"][0]
    reference_density_path = reference / REFERENCE_FILES["electron_density"][0]
    for name, path in (("electron power", reference_power_path),
                       ("electron density", reference_density_path)):
        expected = REFERENCE_FILES[
            "electron_ohmic_power_density" if name == "electron power"
            else "electron_density"][1]
        if sha256(path) != expected:
            raise ValueError(f"locked eduPIC {name} reference differs")
    reference_power = flatten_phase_major(read_matrix(reference_power_path))
    reference_density = flatten_phase_major(read_matrix(reference_density_path))
    reference_average_power = spatial_phase_average(
        reference_power, PHASES, NODES)
    reference_average_density = spatial_phase_average(
        reference_density, PHASES, NODES)

    powers = [float(member["average_power"]) for member in members]
    densities = [float(member["average_density"]) for member in members]
    power_range = relative_range(powers)
    density_range = relative_range(densities)
    ensemble_power = vector_mean([member["power"] for member in members])
    ensemble_phase, ensemble_space = reduced_profiles(ensemble_power)
    reference_phase, reference_space = reduced_profiles(reference_power)
    member_results = []
    local_scatter = []
    phase_scatter = []
    spatial_scatter = []
    for member in members:
        member_power = member["power"]
        member_phase, member_space = reduced_profiles(member_power)
        local_difference = metrics(member_power, ensemble_power)["relative_l2"]
        phase_difference = metrics(member_phase, ensemble_phase)["relative_l2"]
        spatial_difference = metrics(member_space, ensemble_space)["relative_l2"]
        local_scatter.append(local_difference)
        phase_scatter.append(phase_difference)
        spatial_scatter.append(spatial_difference)
        power_ratio = float(member["average_power"]) / reference_average_power
        density_ratio = float(member["average_density"]) / reference_average_density
        member_results.append({
            "seed": member["seed"],
            "branch_report_sha256": member["branch_report_sha256"],
            "average_power_density_W_m-3": member["average_power"],
            "exact_electric_work_W_m-3": member["exact_power"],
            "phase_binned_to_exact_relative_difference": abs(
                float(member["average_power"]) /
                float(member["exact_power"]) - 1.0),
            "average_electron_density_m-3": member["average_density"],
            "candidate_to_reference_power_density_ratio": power_ratio,
            "candidate_to_reference_electron_density_ratio": density_ratio,
            "candidate_to_reference_power_per_electron_ratio":
                power_ratio / density_ratio,
            "member_to_ensemble_local_power_relative_l2": local_difference,
            "member_to_ensemble_phase_power_relative_l2": phase_difference,
            "member_to_ensemble_spatial_power_relative_l2": spatial_difference,
            "peak_resident_set_kib": member["peak_resident_set_kib"],
        })

    gates = {
        "volume_phase_power_density_repeatability":
            power_range <= MAXIMUM_POWER_RELATIVE_RANGE,
        "volume_phase_electron_density_repeatability":
            density_range <= MAXIMUM_DENSITY_RELATIVE_RANGE,
    }
    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "matched_heating_continuation_seed_ensemble",
        "rule_sha256": RULE_SHA256,
        "members": member_results,
        "prospective_repeatability": {
            "power_density_relative_range": power_range,
            "maximum_power_density_relative_range":
                MAXIMUM_POWER_RELATIVE_RANGE,
            "electron_density_relative_range": density_range,
            "maximum_electron_density_relative_range":
                MAXIMUM_DENSITY_RELATIVE_RANGE,
            "gates": gates,
            "all_gates_passed": all(gates.values()),
        },
        "ensemble_cross_code": {
            "mean_power_density_W_m-3": statistics.fmean(powers),
            "mean_electron_density_m-3": statistics.fmean(densities),
            "mean_power_per_electron_ratio": statistics.fmean(
                member["candidate_to_reference_power_per_electron_ratio"]
                for member in member_results),
            "minimum_power_per_electron_ratio": min(
                member["candidate_to_reference_power_per_electron_ratio"]
                for member in member_results),
            "maximum_power_per_electron_ratio": max(
                member["candidate_to_reference_power_per_electron_ratio"]
                for member in member_results),
            "local_power_phase_space": metrics(
                ensemble_power, reference_power),
            "spatially_integrated_phase_power": metrics(
                ensemble_phase, reference_phase),
            "cycle_average_spatial_power": metrics(
                ensemble_space, reference_space),
        },
        "stochastic_localization": {
            "member_to_ensemble_local_power_relative_l2_range":
                [min(local_scatter), max(local_scatter)],
            "member_to_ensemble_phase_power_relative_l2_range":
                [min(phase_scatter), max(phase_scatter)],
            "member_to_ensemble_spatial_power_relative_l2_range":
                [min(spatial_scatter), max(spatial_scatter)],
            "interpretation": (
                "Fine local and phase profiles retain material four-cycle "
                "noise; the cycle-average spatial mismatch exceeds member "
                "scatter and is the stronger persistent localization."),
        },
        "finding": (
            "Volume power, density, and their cross-code ratio are repeatable "
            "across continuation seeds. Fine local power is not converged."),
        "claim_boundary": (
            "These arms share one initial particle realization and vary only "
            "the continuation RNG; this is not full initial-condition uncertainty."),
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
