#!/usr/bin/env python3
"""Audit EEDF-folded argon rates and region-matched eduPIC agreement."""

from __future__ import annotations

import argparse
import bisect
import json
import math
from pathlib import Path

from analyze_edupic_eedf_discrepancy import reference_distribution
from compare_aurorapic_edupic_measurement_pilot import (
    distribution_mean, relative_difference, rows, total_variation,
)
from compare_edupic_phase_space import lower_bin_table, lower_bin_value
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "f497f445030a8d888ef3c9cf93b9ed1492e000eff1336eec9fbfcc578415f9f7")
REFERENCE_EEPF_SHA256 = (
    "8ad8222107c7e61ea9d1a2d36b0550bcd69ded55dbf5fe0bcfd10810dfc56260")
ELECTRON_MASS_KG = 9.10938356e-31
ELEMENTARY_CHARGE_C = 1.60217662e-19
REGIONS = {
    "full_gap": (0.0, 0.025),
    "edupic_center_10pct": (0.01125, 0.01375),
}
CHANNELS = {
    "elastic": "electron_mcc.elastic",
    "excitation": "electron_mcc.excitation",
    "ionization": "electron_mcc.ionization",
}


def candidate_distribution(
    output: Path, region: str,
) -> tuple[list[tuple[float, float, float]], float]:
    histogram = [row for row in rows(output / "phase_eedf.csv")
                 if row["region"] == region]
    moments = [row for row in rows(output / "phase_eedf_moments.csv")
               if row["region"] == region]
    if not histogram or not moments:
        raise ValueError(f"candidate EEDF region is missing: {region}")
    observations = math.fsum(float(row["represented_observations"])
                             for row in moments)
    counts: dict[int, float] = {}
    centers: dict[int, float] = {}
    for row in histogram:
        index = int(row["energy_bin"])
        counts[index] = counts.get(index, 0.0) + float(row["represented_count"])
        center = float(row["energy_eV"])
        if index in centers and centers[index] != center:
            raise ValueError("candidate EEDF energy centers differ by phase")
        centers[index] = center
    indices = sorted(counts)
    if len(indices) < 2 or observations <= 0.0:
        raise ValueError("candidate EEDF is incomplete")
    width = centers[indices[1]] - centers[indices[0]]
    table = [(centers[index] - 0.5 * width,
              centers[index] + 0.5 * width,
              counts[index] / observations) for index in indices]
    residuals = []
    for moment in moments:
        phase = int(moment["phase_bin"])
        represented = float(moment["represented_observations"])
        phase_count = math.fsum(
            float(row["represented_count"]) for row in histogram
            if int(row["phase_bin"]) == phase)
        overflow = float(moment["overflow_fraction"])
        residuals.append(abs(phase_count / represented + overflow - 1.0))
    return table, max(residuals)


def fold_frequency(
    distribution: list[tuple[float, float, float]], cross_section: Path,
    neutral_density: float,
) -> float:
    energies, values = lower_bin_table(cross_section)
    return neutral_density * math.fsum(
        mass * lower_bin_value(energies, values, 0.5 * (low + high)) *
        math.sqrt(2.0 * 0.5 * (low + high) * ELEMENTARY_CHARGE_C /
                  ELECTRON_MASS_KG)
        for low, high, mass in distribution)


def interpolate(x: list[float], values: list[float], target: float) -> float:
    if target < x[0] or target > x[-1]:
        raise ValueError("integration boundary lies outside the diagnostic grid")
    index = bisect.bisect_left(x, target)
    if index < len(x) and x[index] == target:
        return values[index]
    if index == 0 or index == len(x):
        raise ValueError("cannot interpolate diagnostic boundary")
    fraction = (target - x[index - 1]) / (x[index] - x[index - 1])
    return values[index - 1] + fraction * (values[index] - values[index - 1])


def integrate_region(x: list[float], values: list[float],
                     low: float, high: float) -> float:
    if len(x) != len(values) or len(x) < 2 or not low < high:
        raise ValueError("invalid regional integration vectors")
    selected_x = [low] + [value for value in x if low < value < high] + [high]
    selected_y = [interpolate(x, values, value) for value in selected_x]
    return math.fsum(
        0.5 * (selected_y[index] + selected_y[index + 1]) *
        (selected_x[index + 1] - selected_x[index])
        for index in range(len(selected_x) - 1))


def measured_frequency(output: Path, channel: str,
                       bounds: tuple[float, float]) -> float:
    rate_rows = [row for row in rows(output / "spatial_phase_collision_rate.csv")
                 if row["channel"] == channel]
    density_rows = [row for row in rows(output / "spatial_phase_moments.csv")
                    if row["species"] == "electrons"]
    phases = sorted({int(row["phase_bin"]) for row in density_rows})
    numerator = 0.0
    denominator = 0.0
    for phase in phases:
        selected_rate = [row for row in rate_rows
                         if int(row["phase_bin"]) == phase]
        selected_density = [row for row in density_rows
                            if int(row["phase_bin"]) == phase]
        rate_x = [float(row["x_m"]) for row in selected_rate]
        density_x = [float(row["x_m"]) for row in selected_density]
        if rate_x != density_x:
            raise ValueError("collision-rate and density grids differ")
        numerator += integrate_region(
            rate_x, [float(row["mean_event_rate_m-3_s-1"])
                     for row in selected_rate], *bounds)
        denominator += integrate_region(
            density_x, [float(row["number_density_mean_m-3"])
                        for row in selected_density], *bounds)
    if numerator <= 0.0 or denominator <= 0.0:
        raise ValueError("regional event rate or density is not positive")
    return numerator / denominator


def validated_branch(output: Path) -> tuple[Path, dict[str, object]]:
    path = output.parent.parent / "branch-report.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    if (report.get("all_gates_passed") is not True or
            report.get("axis") != "eedf_sampling_region" or
            report.get("branch") != "region_matched" or
            report.get("rule_sha256") != RULE_SHA256):
        raise ValueError("candidate branch identity differs")
    for name, expected in report["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"candidate output differs: {name}")
    return path, report


def analyze(output: Path, reference: Path, cross_sections: dict[str, Path],
            rule_path: Path) -> dict[str, object]:
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("region-matched audit rule differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    branch_path, _ = validated_branch(output)
    if sha256(reference) != REFERENCE_EEPF_SHA256:
        raise ValueError("locked eduPIC EEPF differs")
    expected_hashes = rule["cross_sections"]
    for name, path in cross_sections.items():
        if sha256(path) != expected_hashes[f"electron_{name}_sha256"]:
            raise ValueError(f"locked {name} cross section differs")
    neutral_density = float(rule["fixed_inputs"]["neutral_density_m-3"])
    distributions = {}
    normalizations = {}
    predicted = {}
    measured = {}
    closure = {}
    for region, bounds in REGIONS.items():
        distributions[region], normalizations[region] = candidate_distribution(
            output, region)
        predicted[region] = {
            name: fold_frequency(distributions[region], path, neutral_density)
            for name, path in cross_sections.items()}
        measured[region] = {
            name: measured_frequency(output, channel, bounds)
            for name, channel in CHANNELS.items()}
        closure[region] = {
            name: relative_difference(measured[region][name], value)
            for name, value in predicted[region].items()}

    reference_table = reference_distribution(reference)
    reference_predicted = {
        name: fold_frequency(reference_table, path, neutral_density)
        for name, path in cross_sections.items()}
    candidate_center = distributions["edupic_center_10pct"]
    cross_code = {
        "candidate_center_mean_energy_eV": distribution_mean(candidate_center),
        "reference_center_mean_energy_eV": distribution_mean(reference_table),
        "center_mean_energy_ratio": (
            distribution_mean(candidate_center) / distribution_mean(reference_table)),
        "center_eedf_total_variation": total_variation(
            candidate_center, reference_table),
        "candidate_folded_frequency_s-1": predicted["edupic_center_10pct"],
        "reference_folded_frequency_s-1": reference_predicted,
        "candidate_to_reference_folded_frequency_ratio": {
            name: predicted["edupic_center_10pct"][name] / value
            for name, value in reference_predicted.items()},
    }
    limits = rule["prospective_internal_acceptance"]
    gates = {
        "eedf_probability_normalization": max(normalizations.values()) <=
            float(limits["maximum_eedf_probability_normalization_residual"]),
    }
    for region, limit_key in (
        ("full_gap",
         "maximum_full_gap_measured_to_folded_frequency_relative_difference_per_channel"),
        ("edupic_center_10pct",
         "maximum_center_region_measured_to_folded_frequency_relative_difference_per_channel"),
    ):
        limit = float(limits[limit_key])
        for name, value in closure[region].items():
            gates[f"{region}_{name}_rate_closure"] = value <= limit
    passed = all(gates.values())
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "prospective_region_matched_eedf_collision_audit_result",
        "rule_sha256": RULE_SHA256,
        "candidate_branch_report_sha256": sha256(branch_path),
        "candidate_hashes": {
            "phase_eedf.csv": sha256(output / "phase_eedf.csv"),
            "phase_eedf_moments.csv": sha256(output / "phase_eedf_moments.csv"),
            "spatial_phase_collision_rate.csv": sha256(
                output / "spatial_phase_collision_rate.csv"),
        },
        "reference_eepf_sha256": REFERENCE_EEPF_SHA256,
        "sampling_protocol_correction": (
            "Prior AuroraPIC comparisons used a full-gap EEDF against eduPIC's "
            "central-10-percent EEPF. This result uses matched [0.45L, 0.55L] "
            "sampling for cross-code EEDF and folded-rate observables."),
        "candidate_internal_ledger": {
            "probability_normalization_maximum_residual": normalizations,
            "folded_frequency_s-1": predicted,
            "measured_frequency_s-1": measured,
            "measured_to_folded_relative_difference": closure,
        },
        "internal_gates": gates,
        "all_internal_gates_passed": passed,
        "cross_code_center_region_descriptive": cross_code,
        "interpretation": rule["interpretation"][
            "internal_pass" if passed else "internal_fail"],
        "claim_boundary": rule["interpretation"]["cross_code_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_output", type=Path)
    parser.add_argument("reference_eepf", type=Path)
    parser.add_argument("--elastic", type=Path, required=True)
    parser.add_argument("--excitation", type=Path, required=True)
    parser.add_argument("--ionization", type=Path, required=True)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(
        args.candidate_output.resolve(), args.reference_eepf.resolve(),
        {"elastic": args.elastic.resolve(),
         "excitation": args.excitation.resolve(),
         "ionization": args.ionization.resolve()}, args.rule.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["all_internal_gates_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
