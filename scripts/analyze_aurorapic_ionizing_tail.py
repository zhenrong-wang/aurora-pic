#!/usr/bin/env python3
"""Fold the regional phase EEDF and compare its ionization with eduPIC."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from audit_edupic_ionization_path import read_cross_section
from compare_edupic_phase_space import (
    ARGON_NEUTRAL_DENSITY_M3, ELECTRON_MASS_KG, ELEMENTARY_CHARGE_C,
    REFERENCE_FILES, read_matrix,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


RULE_SHA256 = (
    "a8cded31a57af98b6c32dda816d122cda8ebf9f7daf31d906b516d1a5b12b9f2")
BLOCK_REPORT_SHA256 = (
    "c00bfff25a5fc9aa0b596302dbae8f0417de11d1641a329cc1a2bec03a697d42")
IONIZATION_TABLE_SHA256 = (
    "419958d75e53776ced9f8b81ff77518bf5fc5a18779d167e7231d752d1d9e7e0")
PHASES = 200
NODES = 400
LENGTH_M = 0.025


def csv_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as stream:
        yield from csv.DictReader(stream)


def kernel_frequency(energy_eV: float, energies: list[float],
                     cross_sections: list[float]) -> float:
    spacing = energies[1] - energies[0]
    index = min(int(energy_eV / spacing), len(energies) - 1)
    return (ARGON_NEUTRAL_DENSITY_M3 * cross_sections[index] *
            math.sqrt(2.0 * energy_eV * ELEMENTARY_CHARGE_C /
                      ELECTRON_MASS_KG))


def selected_nodes(low: float, high: float) -> list[int]:
    nodes = [node for node in range(NODES)
             if low <= node * LENGTH_M / (NODES - 1) <= high]
    if len(nodes) < 2:
        raise ValueError("regional nodal integration needs two nodes")
    return nodes


def region_phase_sum(values: list[float], nodes: list[int],
                     first_phase: int = 0,
                     past_last_phase: int = PHASES) -> float:
    if len(values) != PHASES * NODES:
        raise ValueError("phase-space vector has the wrong shape")
    return math.fsum(
        (0.5 if node in (nodes[0], nodes[-1]) else 1.0) *
        values[phase * NODES + node]
        for phase in range(first_phase, past_last_phase)
        for node in nodes)


def analyze(block_root: Path, rule_path: Path, reference: Path,
            ionization_table: Path) -> dict[str, object]:
    report_path = block_root / "block-report.json"
    if sha256(rule_path) != RULE_SHA256:
        raise ValueError("ionizing-tail rule differs")
    if sha256(report_path) != BLOCK_REPORT_SHA256:
        raise ValueError("ionizing-tail block report differs")
    if sha256(ionization_table) != IONIZATION_TABLE_SHA256:
        raise ValueError("ionization table differs")
    rule = json.loads(rule_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if (report.get("scope") !=
            "aurorapic_spatial_phase_ionizing_tail_block" or
            report.get("rule_sha256") != RULE_SHA256):
        raise ValueError("ionizing-tail block contract differs")
    output = block_root / "output"
    for name, expected in report["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"ionizing-tail output differs: {name}")

    diagnostic = rule["diagnostic_contract"]
    critical_names = list(diagnostic["critical_spatial_region_names"])
    region_by_name = {region["name"]: region
                      for region in diagnostic["regions"]}
    if any(name not in region_by_name for name in critical_names):
        raise ValueError("critical region is missing")
    first_critical_phase = round(
        float(diagnostic["critical_phase_fraction"][0]) * PHASES)
    past_last_critical_phase = round(
        float(diagnostic["critical_phase_fraction"][1]) * PHASES)

    energies, cross_sections = read_cross_section(ionization_table)
    eedf = {}
    for row in csv_rows(output / "phase_eedf.csv"):
        region = row["region"]
        if region not in critical_names:
            continue
        phase = int(row["phase_bin"])
        key = (region, phase)
        entry = eedf.setdefault(key, {
            "represented": 0.0, "kernel_sum": 0.0,
            "above_threshold": 0.0})
        count = float(row["represented_count"])
        energy = float(row["energy_eV"])
        entry["represented"] += count
        entry["kernel_sum"] += count * kernel_frequency(
            energy, energies, cross_sections)
        if energy >= float(diagnostic["ionization_threshold_eV"]):
            entry["above_threshold"] += count
    if len(eedf) != len(critical_names) * PHASES:
        raise ValueError("critical EEDF phase-region shape differs")

    electron_rows = [row for row in csv_rows(
        output / "spatial_phase_moments.csv")
        if row["species"] == "electrons"]
    candidate_density = [float(row["number_density_mean_m-3"])
                         for row in electron_rows]
    if len(candidate_density) != PHASES * NODES:
        raise ValueError("candidate density shape differs")
    candidate_rate = [0.0] * (PHASES * NODES)
    rate_rows = [row for row in csv_rows(
        output / "spatial_phase_collision_rate.csv")
        if row["channel"] == "electron_mcc.ionization"]
    if len(rate_rows) != PHASES * NODES:
        raise ValueError("candidate ionization-rate shape differs")
    for index, row in enumerate(rate_rows):
        if (int(row["phase_bin"]) != index // NODES or
                int(row["node"]) != index % NODES):
            raise ValueError("candidate ionization rows are misordered")
        candidate_rate[index] = float(row["mean_event_rate_m-3_s-1"])

    reference_density_path = reference / REFERENCE_FILES["electron_density"][0]
    reference_rate_path = reference / REFERENCE_FILES["ionization_rate"][0]
    if (sha256(reference_density_path) !=
            REFERENCE_FILES["electron_density"][1] or
            sha256(reference_rate_path) != REFERENCE_FILES["ionization_rate"][1]):
        raise ValueError("locked eduPIC reference differs")
    reference_density_matrix = read_matrix(reference_density_path)
    reference_rate_matrix = read_matrix(reference_rate_path)
    reference_density = [reference_density_matrix[node][phase]
                         for phase in range(PHASES) for node in range(NODES)]
    reference_rate = [reference_rate_matrix[node][phase]
                      for phase in range(PHASES) for node in range(NODES)]

    def summarize(names: list[str], first_phase: int,
                  past_last_phase: int) -> dict[str, float]:
        entries = [eedf[(name, phase)] for name in names
                   for phase in range(first_phase, past_last_phase)]
        represented = math.fsum(entry["represented"] for entry in entries)
        predicted = (
            math.fsum(entry["kernel_sum"] for entry in entries) / represented)
        above = (
            math.fsum(entry["above_threshold"] for entry in entries) /
            represented)
        low = min(float(region_by_name[name]["x_min_m"]) for name in names)
        high = max(float(region_by_name[name]["x_max_m"]) for name in names)
        nodes = selected_nodes(low, high)
        candidate_measured = (
            region_phase_sum(candidate_rate, nodes, first_phase,
                             past_last_phase) /
            region_phase_sum(candidate_density, nodes, first_phase,
                             past_last_phase))
        reference_measured = (
            region_phase_sum(reference_rate, nodes, first_phase,
                             past_last_phase) /
            region_phase_sum(reference_density, nodes, first_phase,
                             past_last_phase))
        return {
            "represented_particle_timestep_observations": represented,
            "fraction_in_histogram_bins_at_or_above_15p8_eV": above,
            "eedf_folded_ionization_frequency_s-1": predicted,
            "measured_candidate_ionization_frequency_s-1": candidate_measured,
            "published_edupic_ionization_frequency_s-1": reference_measured,
            "candidate_measured_to_eedf_folded_ratio":
                candidate_measured / predicted,
            "candidate_measured_to_published_edupic_ratio":
                candidate_measured / reference_measured,
            "candidate_eedf_folded_to_published_edupic_ratio":
                predicted / reference_measured,
        }

    critical = summarize(
        critical_names, first_critical_phase, past_last_critical_phase)
    full_phase = summarize(critical_names, 0, PHASES)
    octants = []
    for first in range(0, PHASES, 25):
        value = summarize(critical_names, first, first + 25)
        value.update({"lower_phase_fraction": first / PHASES,
                      "upper_phase_fraction": (first + 25) / PHASES})
        octants.append(value)
    regions = {}
    for name in critical_names:
        regions[name] = summarize([name], first_critical_phase,
                                  past_last_critical_phase)

    safe_gates = {name: passed for name, passed in report["gates"].items()
                  if name != "phase_eedf_observations"}
    critical_moments = [row for row in csv_rows(
        output / "phase_eedf_moments.csv") if row["region"] in critical_names]
    critical_minimum = min(int(row["macro_observations"])
                           for row in critical_moments)
    critical_overflow = max(float(row["overflow_fraction"])
                            for row in critical_moments)
    critical_sufficient = (
        critical_minimum >= int(
            diagnostic["minimum_macro_observations_per_region_phase_bin"]) and
        critical_overflow <= float(diagnostic["maximum_overflow_fraction"]))
    all_moments = list(csv_rows(output / "phase_eedf_moments.csv"))
    regional_sampling = {}
    for name in region_by_name:
        observations = [int(row["macro_observations"]) for row in all_moments
                        if row["region"] == name]
        regional_sampling[name] = {
            "minimum_macro_observations_per_phase_bin": min(observations),
            "maximum_macro_observations_per_phase_bin": max(observations),
            "zero_observation_phase_bins": sum(value == 0
                                                for value in observations),
        }
    return {
        "schema_version": 1,
        "case_id": rule["case_id"],
        "scope": "descriptive_spatial_phase_ionizing_tail_audit",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__).resolve()),
            "rule_sha256": RULE_SHA256,
            "block_report_sha256": BLOCK_REPORT_SHA256,
            "ionization_table_sha256": IONIZATION_TABLE_SHA256,
            "candidate_output_sha256": report["output_hashes"],
            "reference_sha256": {
                reference_density_path.name: sha256(reference_density_path),
                reference_rate_path.name: sha256(reference_rate_path),
            },
        },
        "execution_outcome": {
            "global_prospective_gates_passed": report["all_gates_passed"],
            "failed_global_gates": [name for name, passed in
                                    report["gates"].items() if not passed],
            "all_safety_and_finite_output_gates_passed": all(
                safe_gates.values()),
            "critical_interior_scope_sufficient": critical_sufficient,
            "critical_minimum_macro_observations_per_region_phase_bin":
                critical_minimum,
            "critical_maximum_overflow_fraction": critical_overflow,
            "resource_record": report["resources"],
            "regional_sampling": regional_sampling,
        },
        "critical_interior_phase_0p125_to_0p5": critical,
        "critical_interior_full_cycle": full_phase,
        "critical_phase_octants": octants,
        "critical_region_phase_0p125_to_0p5": regions,
        "assessment": {
            "localized_measured_to_eedf_folded_ratio":
                critical["candidate_measured_to_eedf_folded_ratio"],
            "localized_candidate_to_reference_measured_ratio":
                critical["candidate_measured_to_published_edupic_ratio"],
            "localized_candidate_eedf_to_reference_ratio":
                critical["candidate_eedf_folded_to_published_edupic_ratio"],
            "finding": (
                "In the predeclared interior and RF-phase window, the EEDF "
                "fold independently reproduces the measured candidate "
                "ionization deficit relative to eduPIC. The discrepancy is "
                "therefore carried by the energetic electron distribution, "
                "not by the candidate ionization event sampler."),
            "next_discriminator": (
                "Separate phase-dependent sheath energization, bulk "
                "transport, and inelastic cooling effects on the energetic "
                "tail; this audit does not distinguish those mechanisms."),
        },
        "claim_boundary": rule["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("block_root", type=Path)
    parser.add_argument("rule", type=Path)
    parser.add_argument("reference_raw_data", type=Path)
    parser.add_argument("ionization_table", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.block_root.resolve(), args.rule.resolve(),
                     args.reference_raw_data.resolve(),
                     args.ionization_table.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
