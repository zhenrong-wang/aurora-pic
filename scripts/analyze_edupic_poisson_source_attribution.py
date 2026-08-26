#!/usr/bin/env python3
"""Attribute a matched CCP field gap through the common discrete Poisson operator."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

from analyze_edupic_grid_field_sampling import phase_space_mean_square
from compare_edupic_phase_space import read_matrix


FACTORS = ("boundary_drive", "ion_space_charge", "electron_space_charge")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def solve_dirichlet(rho: list[float], left: float, right: float,
                      dx: float, permittivity: float) -> list[float]:
    """Solve AuroraPIC/eduPIC's nodal Dirichlet Poisson system."""
    if len(rho) < 3 or dx <= 0.0 or permittivity <= 0.0:
        raise ValueError("invalid Dirichlet Poisson contract")
    interior = len(rho) - 2
    diagonal = [0.0] * interior
    solution = [0.0] * interior
    rhs = [rho[index + 1] * dx * dx / permittivity
           for index in range(interior)]
    rhs[0] += left
    rhs[-1] += right
    diagonal[0] = -0.5
    solution[0] = 0.5 * rhs[0]
    for index in range(1, interior):
        denominator = 2.0 + diagonal[index - 1]
        diagonal[index] = 0.0 if index + 1 == interior else -1.0 / denominator
        solution[index] = (rhs[index] + solution[index - 1]) / denominator
    potential = [0.0] * len(rho)
    potential[0] = left
    potential[-1] = right
    potential[-2] = solution[-1]
    for index in range(interior - 2, -1, -1):
        potential[index + 1] = (
            solution[index] - diagonal[index] * potential[index + 2])
    return potential


def recover_field(potential: list[float], rho: list[float], dx: float,
                  permittivity: float) -> list[float]:
    if len(potential) != len(rho) or len(rho) < 3:
        raise ValueError("invalid electric-field recovery contract")
    field = [0.0] * len(rho)
    for index in range(1, len(rho) - 1):
        field[index] = -(potential[index + 1] - potential[index - 1]) / (2.0 * dx)
    field[0] = (-(potential[1] - potential[0]) / dx -
                rho[0] * dx / (2.0 * permittivity))
    field[-1] = (-(potential[-1] - potential[-2]) / dx +
                 rho[-1] * dx / (2.0 * permittivity))
    return field


def component_fields(ne: list[list[float]], ni: list[list[float]],
                     potential: list[list[float]], length: float,
                     permittivity: float,
                     elementary_charge: float) -> dict[str, list[list[float]]]:
    nodes = len(ne)
    phases = len(ne[0])
    if (nodes < 3 or len(ni) != nodes or len(potential) != nodes or
            any(len(row) != phases for matrix in (ne, ni, potential)
                for row in matrix)):
        raise ValueError("phase-space matrices have different shapes")
    dx = length / (nodes - 1)
    result = {factor: [[0.0] * phases for _ in range(nodes)]
              for factor in FACTORS}
    zeros = [0.0] * nodes
    for phase in range(phases):
        drive_phi = solve_dirichlet(
            zeros, potential[0][phase], potential[-1][phase], dx, permittivity)
        ion_rho = [elementary_charge * ni[node][phase]
                   for node in range(nodes)]
        electron_rho = [-elementary_charge * ne[node][phase]
                        for node in range(nodes)]
        profiles = {
            "boundary_drive": recover_field(drive_phi, zeros, dx, permittivity),
            "ion_space_charge": recover_field(
                solve_dirichlet(ion_rho, 0.0, 0.0, dx, permittivity),
                ion_rho, dx, permittivity),
            "electron_space_charge": recover_field(
                solve_dirichlet(electron_rho, 0.0, 0.0, dx, permittivity),
                electron_rho, dx, permittivity),
        }
        for factor, profile in profiles.items():
            for node, value in enumerate(profile):
                result[factor][node][phase] = value
    return result


def add_components(components: dict[str, list[list[float]]]) -> list[list[float]]:
    first = components[FACTORS[0]]
    return [[math.fsum(components[factor][node][phase] for factor in FACTORS)
             for phase in range(len(first[0]))]
            for node in range(len(first))]


def matrix_mean(matrices: list[list[list[float]]]) -> list[list[float]]:
    if not matrices:
        raise ValueError("cannot average an empty matrix ensemble")
    nodes, phases = len(matrices[0]), len(matrices[0][0])
    if any(len(matrix) != nodes or any(len(row) != phases for row in matrix)
           for matrix in matrices):
        raise ValueError("matrix ensemble shapes differ")
    return [[math.fsum(matrix[node][phase] for matrix in matrices) / len(matrices)
             for phase in range(phases)] for node in range(nodes)]


def relative_rms_error(candidate: list[list[float]],
                       reference: list[list[float]]) -> float:
    if len(candidate) != len(reference) or not candidate:
        raise ValueError("RMS matrices differ")
    difference = []
    scale = []
    for left, right in zip(candidate, reference, strict=True):
        if len(left) != len(right):
            raise ValueError("RMS matrix rows differ")
        difference.extend((a - b) ** 2 for a, b in zip(left, right, strict=True))
        scale.extend(b * b for b in right)
    return math.sqrt(math.fsum(difference) / max(math.fsum(scale), 1e-300))


def aurora_data(fields_path: Path, moments_path: Path, nodes: int,
                phases: int, length: float) -> dict[str, list[list[float]]]:
    with fields_path.open(newline="", encoding="utf-8") as stream:
        field_rows = list(csv.DictReader(stream))
    if len(field_rows) != nodes * phases:
        raise ValueError("AuroraPIC phase-field table has the wrong shape")
    potential = [[0.0] * phases for _ in range(nodes)]
    field = [[0.0] * phases for _ in range(nodes)]
    field_samples = set()
    for index, row in enumerate(field_rows):
        phase, node = divmod(index, nodes)
        if (int(row["phase_bin"]) != phase or int(row["node"]) != node or
                abs(float(row["phase_fraction"]) - (phase + 0.5) / phases) > 1e-12 or
                abs(float(row["x_m"]) - node * length / (nodes - 1)) > 1e-12):
            raise ValueError("AuroraPIC phase-field coordinates differ")
        potential[node][phase] = float(row["potential_mean_V"])
        field[node][phase] = float(row["electric_field_mean_V_m"])
        field_samples.add(int(row["samples"]))
    with moments_path.open(newline="", encoding="utf-8") as stream:
        moment_rows = list(csv.DictReader(stream))
    if len(moment_rows) != 2 * nodes * phases:
        raise ValueError("AuroraPIC phase-moment table has the wrong shape")
    ne = [[0.0] * phases for _ in range(nodes)]
    ni = [[0.0] * phases for _ in range(nodes)]
    moment_samples = set()
    expected_species = ("electrons", "ions")
    for index, row in enumerate(moment_rows):
        phase, remainder = divmod(index, 2 * nodes)
        species_id, node = divmod(remainder, nodes)
        if (int(row["phase_bin"]) != phase or int(row["species_id"]) != species_id or
                row["species"] != expected_species[species_id] or
                int(row["node"]) != node or
                abs(float(row["phase_fraction"]) - (phase + 0.5) / phases) > 1e-12 or
                abs(float(row["x_m"]) - node * length / (nodes - 1)) > 1e-12):
            raise ValueError("AuroraPIC phase-moment coordinates differ")
        target = ne if species_id == 0 else ni
        target[node][phase] = float(row["number_density_mean_m-3"])
        moment_samples.add(int(row["samples"]))
    if (len(field_samples) != 1 or len(moment_samples) != 1 or
            field_samples != moment_samples or min(field_samples) <= 0):
        raise ValueError("AuroraPIC phase sample counts differ")
    for matrix in (potential, field, ne, ni):
        if any(not math.isfinite(value) for row in matrix for value in row):
            raise ValueError("AuroraPIC phase-space data are non-finite")
    return {"potential": potential, "field": field, "ne": ne, "ni": ni}


def response(components: dict[str, list[list[float]]], length: float,
             lower_x: float, upper_x: float,
             lower_phase: int, upper_phase: int) -> float:
    return phase_space_mean_square(add_components(components), length,
                                   lower_x, upper_x,
                                   lower_phase, upper_phase)


def shapley_attribution(aurora: dict[str, list[list[float]]],
                        native: dict[str, list[list[float]]],
                        metric) -> tuple[dict[str, float], float, float, float]:
    values: dict[frozenset[str], float] = {}
    for size in range(4):
        for selected in itertools.combinations(FACTORS, size):
            subset = frozenset(selected)
            values[subset] = metric({
                factor: native[factor] if factor in subset else aurora[factor]
                for factor in FACTORS})
    allocation = {}
    weights = {0: 1.0 / 3.0, 1: 1.0 / 6.0, 2: 1.0 / 3.0}
    for factor in FACTORS:
        others = [item for item in FACTORS if item != factor]
        allocation[factor] = math.fsum(
            weights[len(subset)] *
            (values[frozenset((*subset, factor))] - values[frozenset(subset)])
            for size in range(3)
            for subset in itertools.combinations(others, size))
    start = values[frozenset()]
    end = values[frozenset(FACTORS)]
    residual = math.fsum(allocation.values()) - (end - start)
    return allocation, start, end, residual


def analyze(rule_path: Path, native_result_path: Path,
            aurora_reports: list[Path], aurora_fields: list[Path],
            aurora_moments: list[Path], native_inputs: list[list[Path]]) -> dict[str, object]:
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    native_result = json.loads(native_result_path.read_text(encoding="utf-8"))
    contract = rule["discrete_operator_contract"]
    locked = rule["locked_inputs"]
    length = float(contract["length_m"])
    nodes = int(contract["nodes"])
    phases = int(contract["phase_bins"])
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

    aurora_components = [component_fields(item["ne"], item["ni"], item["potential"],
                                          length, permittivity, charge)
                         for item in aurora_raw]
    native_components_members = [component_fields(
        item["ne"], item["ni"], item["potential"], length, permittivity, charge)
        for item in native_raw]
    native_raw_mean = {key: matrix_mean([item[key] for item in native_raw])
                       for key in ("ne", "ni", "potential", "field")}
    native_components = component_fields(
        native_raw_mean["ne"], native_raw_mean["ni"], native_raw_mean["potential"],
        length, permittivity, charge)

    aurora_errors = [relative_rms_error(add_components(parts), raw["field"])
                     for parts, raw in zip(aurora_components, aurora_raw, strict=True)]
    native_errors = [relative_rms_error(add_components(parts), raw["field"])
                     for parts, raw in zip(native_components_members, native_raw, strict=True)]
    metric = lambda parts: response(parts, length, lower_x, upper_x,
                                    lower_phase, upper_phase)
    native_response = metric(native_components)
    native_drive = metric({factor: (native_components[factor] if factor == "boundary_drive"
                                    else [[0.0] * phases for _ in range(nodes)])
                           for factor in FACTORS})
    members = []
    closure_limit_passes = []
    for expected, parts, error in zip(
            locked["aurorapic_members"], aurora_components, aurora_errors, strict=True):
        allocation, start, end, residual = shapley_attribution(parts, native_components, metric)
        gap = end - start
        fractions = {factor: allocation[factor] / gap for factor in FACTORS}
        aurora_drive = metric({factor: (parts[factor] if factor == "boundary_drive"
                                       else [[0.0] * phases for _ in range(nodes)])
                              for factor in FACTORS})
        scale = max(abs(start), abs(end), abs(gap), 1e-300)
        cancellation_amplification = (
            math.fsum(abs(value) for value in allocation.values()) /
            max(abs(gap), 1e-300))
        closure_limit_passes.append(abs(residual) <= 1e-10 * scale)
        members.append({
            "id": expected["id"],
            "operator_reconstruction_relative_rms_error": error,
            "aurorapic_reconstructed_grid_mean_squared_field_V2_m2": start,
            "native_ensemble_reconstructed_grid_mean_squared_field_V2_m2": end,
            "native_minus_aurorapic_gap_V2_m2": gap,
            "boundary_drive_mean_squared_field_V2_m2": aurora_drive,
            "boundary_drive_to_native_ratio": aurora_drive / native_drive,
            "shapley_attribution_V2_m2": allocation,
            "normalized_shapley_fractions": fractions,
            "post_hoc_absolute_attribution_to_gap_ratio":
                cancellation_amplification,
            "shapley_closure_residual_V2_m2": residual,
        })

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
                strict=True))
            for paths, expected in zip(native_inputs, locked["native_members"], strict=True)))
    gates = {
        "locked_hashes_and_reports": hashes_linked,
        "member_counts": len(members) == 2 and len(native_raw) == 3,
        "finite_metrics": all(math.isfinite(value)
                              for member in members
                              for value in (
                                  member["operator_reconstruction_relative_rms_error"],
                                  member["native_minus_aurorapic_gap_V2_m2"],
                                  member["boundary_drive_to_native_ratio"],
                                  *member["normalized_shapley_fractions"].values())),
        "operator_reconstruction": all(error <= 0.005
                                       for error in (*aurora_errors, *native_errors)),
        "positive_field_gap": all(member["native_minus_aurorapic_gap_V2_m2"] > 0.0
                                  for member in members),
        "boundary_drive_parity": all(0.99 <= member["boundary_drive_to_native_ratio"] <= 1.01
                                     for member in members),
        "shapley_closure": all(closure_limit_passes),
    }
    all_gates = all(gates.values())
    dominant = {factor: all(member["normalized_shapley_fractions"][factor] >= 0.50
                            for member in members) for factor in FACTORS}
    supported = [factor for factor, value in dominant.items() if value]
    outcome = (f"{supported[0]}_dominant" if len(supported) == 1
               else "mixed_space_charge_attribution")
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "phase_mean_poisson_source_attribution_result",
        "rule_sha256": sha256(rule_path),
        "gates": gates,
        "all_hash_shape_operator_boundary_gap_and_closure_gates_passed": all_gates,
        "aurorapic_members": members,
        "aurorapic_operator_reconstruction_relative_rms_errors": aurora_errors,
        "native_operator_reconstruction_relative_rms_errors": native_errors,
        "native_ensemble_reconstructed_grid_mean_squared_field_V2_m2": native_response,
        "native_ensemble_boundary_drive_mean_squared_field_V2_m2": native_drive,
        "prospective_decision_outcome": {
            "interpretation_allowed": all_gates,
            "dominance_supported": dominant if all_gates else {factor: False for factor in FACTORS},
            "result": outcome if all_gates else "interpretation_forbidden",
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
            "The Shapley allocation is an exact algebraic decomposition of a mature "
            "phase-mean field-energy gap under the shared linear Poisson operator. "
            "Its component substitutions are not independently evolved states. The "
            "absolute-attribution/gap ratio was added after observing the locked "
            "outcome as a transparent conditioning diagnostic; it does not alter "
            "any preregistered gate or decision."),
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rule", type=Path)
    parser.add_argument("native_result", type=Path)
    parser.add_argument("aurorapic_reports", nargs=2, type=Path)
    parser.add_argument("aurorapic_fields", nargs=2, type=Path)
    parser.add_argument("aurorapic_moments", nargs=2, type=Path)
    parser.add_argument("native_inputs", nargs=12, type=Path,
                        help="per seed: ne_xt ni_xt pot_xt efield_xt")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    native_inputs = [args.native_inputs[index:index + 4]
                     for index in range(0, 12, 4)]
    result = analyze(args.rule, args.native_result, args.aurorapic_reports,
                     args.aurorapic_fields, args.aurorapic_moments, native_inputs)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result[
        "all_hash_shape_operator_boundary_gap_and_closure_gates_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
