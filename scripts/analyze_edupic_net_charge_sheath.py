#!/usr/bin/env python3
"""Compare stable net-charge and sheath modes in a matched CCP window."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from analyze_edupic_grid_field_sampling import boundary_value, phase_space_mean_square
from analyze_edupic_poisson_source_attribution import aurora_data, component_fields
from compare_edupic_phase_space import read_matrix


THRESHOLDS = (0.8, 0.9, 0.95)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def left_crossing(values: list[float], threshold: float, length: float) -> float:
    if len(values) < 2 or values[0] >= threshold:
        return 0.0
    dx = length / (len(values) - 1)
    for index in range(1, len(values)):
        if values[index] >= threshold:
            previous = values[index - 1]
            if values[index] == previous:
                return index * dx
            fraction = (threshold - previous) / (values[index] - previous)
            return (index - 1 + fraction) * dx
    raise ValueError("density-ratio sheath crossing is absent")


def integrate_profile(values: list[float], length: float,
                      lower: float, upper: float) -> float:
    if len(values) < 2 or not 0.0 <= lower < upper <= length:
        raise ValueError("invalid profile integration contract")
    dx = length / (len(values) - 1)
    points = [(lower, boundary_value(values, lower, length))]
    points.extend((node * dx, value) for node, value in enumerate(values)
                  if lower < node * dx < upper)
    points.append((upper, boundary_value(values, upper, length)))
    return math.fsum(0.5 * (left[1] + right[1]) * (right[0] - left[0])
                     for left, right in zip(points, points[1:]))


def spatial_mean_product(first: list[float], second: list[float], length: float,
                         lower: float, upper: float) -> float:
    if len(first) != len(second):
        raise ValueError("product profiles differ")
    dx = length / (len(first) - 1)
    coordinates = [lower]
    coordinates.extend(node * dx for node in range(len(first))
                       if lower < node * dx < upper)
    coordinates.append(upper)
    products = [boundary_value(first, x, length) *
                boundary_value(second, x, length) for x in coordinates]
    return math.fsum(0.5 * (a + b) * (right - left)
                     for left, right, a, b in zip(
                         coordinates, coordinates[1:], products, products[1:])) / (upper - lower)


def phase_mean_product(first: list[list[float]], second: list[list[float]],
                       length: float, lower_x: float, upper_x: float,
                       lower_phase: int, upper_phase: int) -> float:
    return math.fsum(spatial_mean_product(
        [first[node][phase] for node in range(len(first))],
        [second[node][phase] for node in range(len(second))],
        length, lower_x, upper_x) for phase in range(lower_phase, upper_phase)
    ) / (upper_phase - lower_phase)


def member_metrics(raw: dict[str, list[list[float]]], length: float,
                   lower_x: float, upper_x: float,
                   lower_phase: int, upper_phase: int,
                   permittivity: float, charge: float) -> dict[str, float]:
    ne, ni = raw["ne"], raw["ni"]
    nodes, phases = len(ne), len(ne[0])
    separation = [[0.0] * phases for _ in range(nodes)]
    ratios = [[0.0] * nodes for _ in range(phases)]
    for node in range(nodes):
        for phase in range(phases):
            denominator = ni[node][phase] + ne[node][phase]
            if denominator <= 0.0 or ni[node][phase] <= 0.0:
                raise ValueError("non-positive density in sheath reduction")
            separation[node][phase] = (ni[node][phase] - ne[node][phase]) / denominator
            ratios[phase][node] = ne[node][phase] / ni[node][phase]
    components = component_fields(ne, ni, raw["potential"], length,
                                  permittivity, charge)
    drive = components["boundary_drive"]
    space_charge = [[components["ion_space_charge"][node][phase] +
                     components["electron_space_charge"][node][phase]
                     for phase in range(phases)] for node in range(nodes)]
    drive_e2 = phase_space_mean_square(
        drive, length, lower_x, upper_x, lower_phase, upper_phase)
    cross = phase_mean_product(drive, space_charge, length, lower_x, upper_x,
                               lower_phase, upper_phase)
    metrics = {
        "critical_charge_separation_rms": math.sqrt(phase_space_mean_square(
            separation, length, lower_x, upper_x, lower_phase, upper_phase)),
        "drive_space_charge_cancellation_fraction": -2.0 * cross / drive_e2,
    }
    for threshold in THRESHOLDS:
        widths = []
        charges = []
        for phase in range(lower_phase, upper_phase):
            width = left_crossing(ratios[phase], threshold, length)
            widths.append(width)
            rho_positive = [max(charge * (ni[node][phase] - ne[node][phase]), 0.0)
                            for node in range(nodes)]
            charges.append(integrate_profile(rho_positive, length, 0.0, width))
        label = f"{threshold:.2f}"
        metrics[f"left_sheath_width_ne_over_ni_{label}_m"] = math.fsum(widths) / len(widths)
        metrics[f"positive_sheath_charge_ne_over_ni_{label}_C_m2"] = math.fsum(charges) / len(charges)
    return metrics


def relative_ranges(members: list[dict[str, float]]) -> dict[str, float]:
    keys = members[0].keys()
    result = {}
    for key in keys:
        values = [member[key] for member in members]
        mean = math.fsum(values) / len(values)
        result[key] = (max(values) - min(values)) / max(abs(mean), 1e-300)
    return result


def analyze(rule_path: Path, native_result_path: Path,
            aurora_reports: list[Path], aurora_fields: list[Path],
            aurora_moments: list[Path], native_inputs: list[list[Path]]) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    locked = rule["locked_inputs"]
    contract = rule["reduction_contract"]
    length = float(contract["length_m"])
    nodes, phases = int(contract["nodes"]), int(contract["phase_bins"])
    lower_phase, upper_phase = map(int, contract["critical_phase_bins"])
    lower_x, upper_x = [float(value) * length for value in contract["critical_x_over_L"]]
    permittivity = float(contract["permittivity_F_m"])
    charge = float(contract["elementary_charge_C"])
    reports = [json.loads(path.read_text(encoding="utf-8")) for path in aurora_reports]
    aurora_raw = [aurora_data(fields, moments, nodes, phases, length)
                  for fields, moments in zip(aurora_fields, aurora_moments, strict=True)]
    native_raw = []
    for paths in native_inputs:
        ne, ni, potential, field = [read_matrix(path, nodes, phases) for path in paths]
        native_raw.append({"ne": ne, "ni": ni, "potential": potential, "field": field})
    aurora_metrics = [member_metrics(raw, length, lower_x, upper_x,
                                     lower_phase, upper_phase, permittivity, charge)
                      for raw in aurora_raw]
    native_metrics = [member_metrics(raw, length, lower_x, upper_x,
                                     lower_phase, upper_phase, permittivity, charge)
                      for raw in native_raw]
    keys = tuple(aurora_metrics[0].keys())
    native_means = {key: math.fsum(member[key] for member in native_metrics) /
                    len(native_metrics) for key in keys}
    ratios = [{key: member[key] / native_means[key] for key in keys}
              for member in aurora_metrics]
    aurora_ranges = relative_ranges(aurora_metrics)
    native_ranges = relative_ranges(native_metrics)
    hashes_linked = (
        sha256(native_result_path) == locked["native_result_sha256"] and
        all(sha256(report_path) == expected["runner_report_sha256"] and
            sha256(fields_path) == expected["spatial_phase_fields_sha256"] and
            sha256(moments_path) == expected["spatial_phase_moments_sha256"] and
            report.get("all_gates_passed") is True and report.get("state_id") == expected["id"]
            for report_path, fields_path, moments_path, report, expected in zip(
                aurora_reports, aurora_fields, aurora_moments, reports,
                locked["aurorapic_members"], strict=True)) and
        all(all(sha256(path) == expected[key] for path, key in zip(
                paths, ("ne_xt_sha256", "ni_xt_sha256", "pot_xt_sha256", "efield_xt_sha256"),
                strict=True)) for paths, expected in zip(
                    native_inputs, locked["native_members"], strict=True)))
    monotonic = all(
        member["left_sheath_width_ne_over_ni_0.80_m"] <=
        member["left_sheath_width_ne_over_ni_0.90_m"] <=
        member["left_sheath_width_ne_over_ni_0.95_m"]
        for member in (*aurora_metrics, *native_metrics))
    gates = {
        "locked_hashes_and_reports": hashes_linked,
        "member_counts": len(aurora_metrics) == 2 and len(native_metrics) == 3,
        "finite_positive_metrics": all(math.isfinite(value) and value > 0.0
                                       for member in (*aurora_metrics, *native_metrics)
                                       for value in member.values()),
        "ordered_sheath_threshold_crossings": monotonic,
        "aurorapic_repeatability": all(value <= 0.10 for value in aurora_ranges.values()),
        "native_repeatability": all(value <= 0.15 for value in native_ranges.values()),
    }
    all_gates = all(gates.values())
    cancellation = all(member["drive_space_charge_cancellation_fraction"] >= 1.10
                       for member in ratios)
    separation = all(member["critical_charge_separation_rms"] >= 1.10
                     for member in ratios)
    wider = all(all(member[f"left_sheath_width_ne_over_ni_{threshold:.2f}_m"] >= 1.10
                    for threshold in THRESHOLDS) for member in ratios)
    lower_charge = all(all(member[
        f"positive_sheath_charge_ne_over_ni_{threshold:.2f}_C_m2"] <= 0.90
        for threshold in THRESHOLDS) for member in ratios)
    associations = {
        "stronger_space_charge_cancellation_supported": all_gates and cancellation,
        "larger_charge_separation_supported": all_gates and separation,
        "wider_left_sheath_supported": all_gates and wider,
        "lower_positive_sheath_charge_supported": all_gates and lower_charge,
    }
    supported = [key.removesuffix("_supported") for key, value in associations.items()
                 if value]
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "net_charge_and_sheath_structure_result",
        "rule_sha256": sha256(rule_path),
        "gates": gates,
        "all_hash_shape_crossing_positive_and_repeatability_gates_passed": all_gates,
        "aurorapic_members": [{"id": expected["id"], **metrics}
                              for expected, metrics in zip(
                                  locked["aurorapic_members"], aurora_metrics, strict=True)],
        "native_members": [{"seed": expected["seed"], **metrics}
                           for expected, metrics in zip(
                               locked["native_members"], native_metrics, strict=True)],
        "native_means": native_means,
        "aurorapic_member_to_native_mean_ratios": ratios,
        "aurorapic_relative_ranges": aurora_ranges,
        "native_relative_ranges": native_ranges,
        "prospective_decision_outcome": {
            "interpretation_allowed": all_gates,
            **associations,
            "result": ("_and_".join(supported) if all_gates and supported
                       else "mixed_or_intermediate_net_charge_sheath_result"),
        },
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "native_result_sha256": sha256(native_result_path),
            "aurorapic_field_sha256": [sha256(path) for path in aurora_fields],
            "aurorapic_moment_sha256": [sha256(path) for path in aurora_moments],
            "native_input_sha256": [[sha256(path) for path in paths]
                                    for paths in native_inputs],
        },
        "interpretation_note": (
            "Density-ratio crossings are threshold-sensitive sheath proxies, not "
            "unique physical sheath edges. All source-field quantities use net "
            "charge and avoid independent electron/ion substitution."),
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("native_result", type=Path)
    parser.add_argument("aurorapic_reports", nargs=2, type=Path)
    parser.add_argument("aurorapic_fields", nargs=2, type=Path)
    parser.add_argument("aurorapic_moments", nargs=2, type=Path)
    parser.add_argument("native_inputs", nargs=12, type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = analyze(args.rule, args.native_result, args.aurorapic_reports,
                     args.aurorapic_fields, args.aurorapic_moments,
                     [args.native_inputs[index:index + 4]
                      for index in range(0, 12, 4)])
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result[
        "all_hash_shape_crossing_positive_and_repeatability_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
