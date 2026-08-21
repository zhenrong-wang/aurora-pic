#!/usr/bin/env python3
"""Localize the long-window AuroraPIC/eduPIC ionization discrepancy."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from compare_aurorapic_edupic_measurement_pilot import rows
from compare_edupic_phase_space import (
    ELEMENTARY_CHARGE_C, REFERENCE_FILES, flatten_phase_major, metrics,
    read_matrix,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


LONG_REPORT_SHA256 = (
    "1d87b9ebfbe3668513c844fc5314868e1bae625df259c7f659981f9ba03a9b4a")
IONIZATION_PATH_AUDIT_SHA256 = (
    "10a3ee8261a23486525363e7674c228e63055ea7e743c7c33bba81f6a263f3ee")
PHASES = 200
NODES = 400
BAND_BOUNDARIES = [0, 40, 80, 160, 240, 320, 360, 400]
PHASE_BOUNDARIES = list(range(0, PHASES + 1, 25))
LENGTH_M = 0.025


def validate_phase_grid(data: list[dict[str, str]], count_field: str,
                        species: str | None = None) -> None:
    if len(data) != PHASES * NODES:
        raise ValueError("long-window phase-space shape differs")
    counts = set()
    for index, row in enumerate(data):
        phase, node = divmod(index, NODES)
        if int(row["phase_bin"]) != phase or int(row["node"]) != node:
            raise ValueError("long-window rows are not phase-major/node-minor")
        if abs(float(row["phase_fraction"]) - (phase + 0.5) / PHASES) > 1e-12:
            raise ValueError("long-window phase centers differ")
        if abs(float(row["x_m"]) - node * LENGTH_M / (NODES - 1)) > 1e-12:
            raise ValueError("long-window spatial coordinates differ")
        if species is not None and row["species"] != species:
            raise ValueError("long-window species differs")
        count = int(row[count_field])
        if count <= 0:
            raise ValueError("long-window sample count is not positive")
        counts.add(count)
    if len(counts) != 1:
        raise ValueError("long-window sample counts are not uniform")


def weighted_sum(values: list[float], first_node: int = 0,
                 past_last_node: int | None = None, phases: int = PHASES,
                 nodes: int = NODES) -> float:
    """Sum a phase-major nodal field using whole-gap trapezoid weights."""
    if past_last_node is None:
        past_last_node = nodes
    if (len(values) != phases * nodes or first_node < 0 or
            past_last_node > nodes or first_node >= past_last_node):
        raise ValueError("invalid phase-space integration contract")
    return math.fsum(
        (0.5 if node in (0, nodes - 1) else 1.0) *
        values[phase * nodes + node]
        for phase in range(phases)
        for node in range(first_node, past_last_node))


def band_summary(candidate_rate: list[float], reference_rate: list[float],
                 candidate_density: list[float], reference_density: list[float],
                 candidate_power: list[float], reference_power: list[float],
                 boundaries: list[int] = BAND_BOUNDARIES,
                 phases: int = PHASES, nodes: int = NODES
                 ) -> list[dict[str, object]]:
    if (boundaries[0] != 0 or boundaries[-1] != nodes or
            any(right <= left for left, right in
                zip(boundaries, boundaries[1:]))):
        raise ValueError("invalid spatial-band boundaries")
    vectors = (candidate_rate, reference_rate, candidate_density,
               reference_density, candidate_power, reference_power)
    if any(len(values) != phases * nodes for values in vectors):
        raise ValueError("spatial-band vector has the wrong shape")
    candidate_total = weighted_sum(candidate_rate, phases=phases, nodes=nodes)
    reference_total = weighted_sum(reference_rate, phases=phases, nodes=nodes)
    total_gap = reference_total - candidate_total
    result = []
    for low, high in zip(boundaries, boundaries[1:]):
        sums = [weighted_sum(values, low, high, phases, nodes)
                for values in vectors]
        candidate_source, reference_source, candidate_number, reference_number, \
            candidate_heating, reference_heating = sums
        source_ratio = candidate_source / reference_source
        density_ratio = candidate_number / reference_number
        power_ratio = candidate_heating / reference_heating
        source_gap = reference_source - candidate_source
        result.append({
            "first_node": low,
            "past_last_node": high,
            "lower_gap_fraction": low / nodes,
            "upper_gap_fraction": high / nodes,
            "reference_ionization_fraction":
                reference_source / reference_total,
            "candidate_ionization_fraction":
                candidate_source / candidate_total,
            "candidate_to_reference_ionization_rate_ratio": source_ratio,
            "candidate_to_reference_electron_density_ratio": density_ratio,
            "candidate_to_reference_effective_ionization_frequency_ratio":
                source_ratio / density_ratio,
            "reference_minus_candidate_source_fraction_of_reference_total":
                source_gap / reference_total,
            "fraction_of_net_source_gap": (
                source_gap / total_gap if total_gap != 0.0 else None),
            "candidate_to_reference_signed_electron_power_ratio": power_ratio,
            "candidate_to_reference_signed_power_per_electron_ratio":
                power_ratio / density_ratio,
        })
    return result


def phase_profile(values: list[float], phases: int = PHASES,
                  nodes: int = NODES) -> list[float]:
    """Trapezoid-integrate each phase bin over the spatial nodes."""
    if len(values) != phases * nodes:
        raise ValueError("phase-profile vector has the wrong shape")
    weights = [0.5] + [1.0] * (nodes - 2) + [0.5]
    return [math.fsum(weights[node] * values[phase * nodes + node]
                      for node in range(nodes))
            for phase in range(phases)]


def phase_band_summary(candidate_rate: list[float], reference_rate: list[float],
                       candidate_density: list[float],
                       reference_density: list[float],
                       candidate_power: list[float],
                       reference_power: list[float],
                       boundaries: list[int] = PHASE_BOUNDARIES,
                       phases: int = PHASES, nodes: int = NODES
                       ) -> list[dict[str, object]]:
    if (boundaries[0] != 0 or boundaries[-1] != phases or
            any(right <= left for left, right in
                zip(boundaries, boundaries[1:]))):
        raise ValueError("invalid phase-band boundaries")
    profiles = [phase_profile(values, phases, nodes) for values in (
        candidate_rate, reference_rate, candidate_density, reference_density,
        candidate_power, reference_power)]
    candidate_total = math.fsum(profiles[0])
    reference_total = math.fsum(profiles[1])
    total_gap = reference_total - candidate_total
    result = []
    for low, high in zip(boundaries, boundaries[1:]):
        sums = [math.fsum(profile[low:high]) for profile in profiles]
        candidate_source, reference_source, candidate_number, reference_number, \
            candidate_heating, reference_heating = sums
        source_ratio = candidate_source / reference_source
        density_ratio = candidate_number / reference_number
        power_ratio = candidate_heating / reference_heating
        result.append({
            "first_phase_bin": low,
            "past_last_phase_bin": high,
            "lower_phase_fraction": low / phases,
            "upper_phase_fraction": high / phases,
            "reference_ionization_fraction":
                reference_source / reference_total,
            "candidate_to_reference_ionization_rate_ratio": source_ratio,
            "candidate_to_reference_effective_ionization_frequency_ratio":
                source_ratio / density_ratio,
            "fraction_of_net_source_gap": (
                (reference_source - candidate_source) / total_gap
                if total_gap != 0.0 else None),
            "candidate_to_reference_signed_electron_power_ratio": power_ratio,
            "candidate_to_reference_signed_power_per_electron_ratio":
                power_ratio / density_ratio,
        })
    return result


def analyze(root: Path, reference: Path,
            ionization_path_audit: Path) -> dict[str, object]:
    report_path = root / "branch-report.json"
    if sha256(report_path) != LONG_REPORT_SHA256:
        raise ValueError("long-window report differs")
    if sha256(ionization_path_audit) != IONIZATION_PATH_AUDIT_SHA256:
        raise ValueError("ionization-path audit differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("all_gates_passed") is not True:
        raise ValueError("long-window execution did not pass")
    output = root / "measurement" / "output"
    required = ("spatial_phase_fields.csv", "spatial_phase_moments.csv",
                "spatial_phase_collision_rate.csv")
    for name in required:
        if sha256(output / name) != report["output_hashes"][name]:
            raise ValueError(f"long-window output differs: {name}")

    metadata = json.loads(
        (output / "spatial_average_metadata.json").read_text(encoding="utf-8"))
    if (metadata.get("phase_bins") != PHASES or
            metadata.get("samples") != 48000 or
            set(metadata.get("phase_bin_samples", [])) != {240} or
            metadata.get("sampling_order") != "pre_collision" or
            metadata.get("complete") is not True):
        raise ValueError("long-window sampling protocol differs")

    fields = rows(output / "spatial_phase_fields.csv")
    moments = [row for row in rows(output / "spatial_phase_moments.csv")
               if row["species"] == "electrons"]
    collision_rows = [
        row for row in rows(output / "spatial_phase_collision_rate.csv")
        if row["channel"] == "ionization" or
        row["channel"].endswith(".ionization")]
    if (len(fields) != PHASES * NODES or len(moments) != len(fields) or
            len(collision_rows) != len(fields)):
        raise ValueError("long-window phase-space shape differs")
    if sorted({row["channel"] for row in collision_rows}) != [
            "electron_mcc.ionization"]:
        raise ValueError("long-window ionization channel differs")
    validate_phase_grid(fields, "samples")
    validate_phase_grid(moments, "samples", "electrons")
    validate_phase_grid(collision_rows, "timesteps")

    candidate_density = [float(row["number_density_mean_m-3"])
                         for row in moments]
    candidate_rate = [float(row["mean_event_rate_m-3_s-1"])
                      for row in collision_rows]
    candidate_power = [
        -ELEMENTARY_CHARGE_C * density * float(moment["mean_velocity_x"]) *
        float(field["electric_field_mean_V_m"])
        for density, moment, field in zip(candidate_density, moments, fields)]

    reference_values = {}
    for name in ("electron_density", "electron_ohmic_power_density",
                 "ionization_rate"):
        filename, expected = REFERENCE_FILES[name]
        path = reference / filename
        if sha256(path) != expected:
            raise ValueError(f"locked eduPIC reference differs: {filename}")
        reference_values[name] = flatten_phase_major(read_matrix(path))

    reference_density = reference_values["electron_density"]
    reference_rate = reference_values["ionization_rate"]
    reference_power = reference_values["electron_ohmic_power_density"]
    candidate_source = weighted_sum(candidate_rate)
    reference_source = weighted_sum(reference_rate)
    candidate_number = weighted_sum(candidate_density)
    reference_number = weighted_sum(reference_density)
    candidate_heating = weighted_sum(candidate_power)
    reference_heating = weighted_sum(reference_power)
    source_ratio = candidate_source / reference_source
    density_ratio = candidate_number / reference_number
    power_ratio = candidate_heating / reference_heating
    bands = band_summary(
        candidate_rate, reference_rate, candidate_density, reference_density,
        candidate_power, reference_power)
    phase_bands = phase_band_summary(
        candidate_rate, reference_rate, candidate_density, reference_density,
        candidate_power, reference_power)
    central_gap_fraction = math.fsum(
        float(band["fraction_of_net_source_gap"]) for band in bands[2:4])
    first_half_gap_fraction = math.fsum(
        float(band["fraction_of_net_source_gap"]) for band in phase_bands[:4])

    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "descriptive_spatial_ionization_localization",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__).resolve()),
            "long_window_report_sha256": LONG_REPORT_SHA256,
            "ionization_path_audit_sha256": IONIZATION_PATH_AUDIT_SHA256,
            "candidate_output_sha256": {
                name: report["output_hashes"][name] for name in required},
            "reference_sha256": {
                REFERENCE_FILES[name][0]: REFERENCE_FILES[name][1]
                for name in ("electron_density",
                             "electron_ohmic_power_density",
                             "ionization_rate")},
        },
        "comparison_contract": {
            "phase_bins": PHASES,
            "spatial_nodes": NODES,
            "measurement_cycles": 12,
            "samples_per_phase_bin": 240,
            "phase_alignment": "direct_no_fitted_shift",
            "orientation": "direct_powered_left_no_reflection",
            "spatial_quadrature": "whole_gap_trapezoid_nodal_weights",
            "band_boundaries_node_indices": BAND_BOUNDARIES,
            "phase_band_boundaries": PHASE_BOUNDARIES,
            "acceptance_thresholds_declared": False,
        },
        "whole_gap": {
            "candidate_to_reference_ionization_rate_ratio": source_ratio,
            "candidate_to_reference_electron_density_ratio": density_ratio,
            "candidate_to_reference_effective_ionization_frequency_ratio":
                source_ratio / density_ratio,
            "candidate_to_reference_electron_power_ratio": power_ratio,
            "candidate_to_reference_power_per_electron_ratio":
                power_ratio / density_ratio,
        },
        "phase_space_comparisons": {
            "ionization_rate": metrics(candidate_rate, reference_rate),
            "electron_density": metrics(candidate_density, reference_density),
            "electron_power": metrics(candidate_power, reference_power),
        },
        "spatial_bands": bands,
        "phase_octants": phase_bands,
        "localization": {
            "fraction_of_net_source_gap_in_0p2_to_0p6_gap":
                central_gap_fraction,
            "fraction_of_net_source_gap_in_first_half_rf_cycle":
                first_half_gap_fraction,
            "finding": (
                "The 0.2-0.6 gap region and first half of the RF cycle account "
                "for most of the small net ionization-rate deficit. Effective "
                "ionization per electron is lower through nearly the entire "
                "gap despite a higher electron density."),
            "interpretation": (
                "This favors a distributed energetic-tail/heating mismatch "
                "over an ionization sampler defect confined to an electrode."),
        },
        "claim_boundary": (
            "This is a checksum-bound descriptive comparison with published "
            "simulation output. It has no prospective spatial acceptance "
            "threshold and does not uniquely identify the kinetic cause or "
            "constitute experimental validation."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_root", type=Path)
    parser.add_argument("reference_raw_data", type=Path)
    parser.add_argument("ionization_path_audit", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.candidate_root.resolve(),
                     args.reference_raw_data.resolve(),
                     args.ionization_path_audit.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
