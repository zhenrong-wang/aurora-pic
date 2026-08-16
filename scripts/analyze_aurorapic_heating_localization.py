#!/usr/bin/env python3
"""Localize matched AuroraPIC/eduPIC heating differences in phase and space."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from compare_aurorapic_edupic_measurement_pilot import rows
from compare_edupic_phase_space import (
    ELEMENTARY_CHARGE_C, REFERENCE_FILES, flatten_phase_major, metrics,
    read_matrix, spatial_phase_average,
)
from run_aurorapic_edupic_pilot import atomic_json, sha256


BRANCH_REPORT_SHA256 = (
    "4171678276ffb9fee0eecb843f9d9c0b44ecb5bbf13c65207618a9cd76820fc5")
PHASES = 200
NODES = 400


def reduced_profiles(values: list[float], phases: int = PHASES,
                     nodes: int = NODES) -> tuple[list[float], list[float]]:
    if len(values) != phases * nodes:
        raise ValueError("phase-space vector has the wrong shape")
    weights = [0.5] + [1.0] * (nodes - 2) + [0.5]
    phase_profile = [
        math.fsum(weights[node] * values[phase * nodes + node]
                  for node in range(nodes)) / (nodes - 1)
        for phase in range(phases)
    ]
    spatial_profile = [
        math.fsum(values[phase * nodes + node] for phase in range(phases)) /
        phases for node in range(nodes)
    ]
    return phase_profile, spatial_profile


def grouped_means(candidate: list[float], reference: list[float],
                  boundaries: list[int], scale: int) -> list[dict[str, float | int]]:
    if (len(candidate) != len(reference) or not candidate or
            boundaries[0] != 0 or boundaries[-1] != len(candidate) or
            any(right <= left for left, right in zip(boundaries, boundaries[1:]))):
        raise ValueError("grouped-mean contract is invalid")
    result = []
    for low, high in zip(boundaries, boundaries[1:]):
        candidate_mean = math.fsum(candidate[low:high]) / (high - low)
        reference_mean = math.fsum(reference[low:high]) / (high - low)
        result.append({
            "first_index": low,
            "last_index_inclusive": high - 1,
            "lower_fraction": low / scale,
            "upper_fraction": high / scale,
            "candidate_mean_W_m-3": candidate_mean,
            "reference_mean_W_m-3": reference_mean,
            "candidate_minus_reference_W_m-3":
                candidate_mean - reference_mean,
        })
    return result


def analyze(root: Path, reference: Path) -> dict[str, object]:
    report_path = root / "branch-report.json"
    if sha256(report_path) != BRANCH_REPORT_SHA256:
        raise ValueError("matched-heating branch report differs")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("all_gates_passed") is not True:
        raise ValueError("matched-heating branch did not pass")
    output = root / "measurement" / "output"
    for name, expected in report["output_hashes"].items():
        if sha256(output / name) != expected:
            raise ValueError(f"matched candidate output differs: {name}")
    metadata = json.loads(
        (output / "spatial_average_metadata.json").read_text(encoding="utf-8"))
    if (metadata.get("sampling_order") != "pre_collision" or
            metadata.get("phase_bins") != PHASES or
            metadata.get("complete") is not True):
        raise ValueError("candidate sampling protocol differs")

    fields = rows(output / "spatial_phase_fields.csv")
    moments = [row for row in rows(output / "spatial_phase_moments.csv")
               if row["species"] == "electrons"]
    if len(fields) != PHASES * NODES or len(moments) != len(fields):
        raise ValueError("candidate phase-space shape differs")
    candidate = {
        "electron_density": [float(row["number_density_mean_m-3"])
                             for row in moments],
        "electric_field": [float(row["electric_field_mean_V_m"])
                           for row in fields],
    }
    candidate["electron_current_density"] = [
        -ELEMENTARY_CHARGE_C * density * float(row["mean_velocity_x"])
        for density, row in zip(candidate["electron_density"], moments)
    ]
    candidate["electron_ohmic_power_density"] = [
        current * field for current, field in zip(
            candidate["electron_current_density"], candidate["electric_field"])
    ]

    reference_values = {}
    for name in candidate:
        filename, expected = REFERENCE_FILES[name]
        path = reference / filename
        if sha256(path) != expected:
            raise ValueError(f"locked eduPIC reference differs: {filename}")
        reference_values[name] = flatten_phase_major(read_matrix(path))

    comparisons = {
        name: metrics(values, reference_values[name])
        for name, values in candidate.items()
    }
    candidate_power = candidate["electron_ohmic_power_density"]
    reference_power = reference_values["electron_ohmic_power_density"]
    candidate_density = spatial_phase_average(
        candidate["electron_density"], PHASES, NODES)
    reference_density = spatial_phase_average(
        reference_values["electron_density"], PHASES, NODES)
    candidate_power_average = spatial_phase_average(
        candidate_power, PHASES, NODES)
    reference_power_average = spatial_phase_average(
        reference_power, PHASES, NODES)
    density_ratio = candidate_density / reference_density
    power_ratio = candidate_power_average / reference_power_average
    candidate_phase, candidate_space = reduced_profiles(candidate_power)
    reference_phase, reference_space = reduced_profiles(reference_power)

    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "matched_electron_heating_phase_space_localization",
        "branch_report_sha256": BRANCH_REPORT_SHA256,
        "comparison_contract": {
            "phase_bins": PHASES,
            "spatial_nodes": NODES,
            "phase_alignment": "direct_no_fitted_shift",
            "orientation": "direct_powered_left_no_reflection",
            "sampling_order": "pre_collision",
            "acceptance_thresholds_declared": False,
        },
        "volume_phase_averages": {
            "candidate_power_density_W_m-3": candidate_power_average,
            "reference_power_density_W_m-3": reference_power_average,
            "candidate_to_reference_power_density_ratio": power_ratio,
            "candidate_electron_density_m-3": candidate_density,
            "reference_electron_density_m-3": reference_density,
            "candidate_to_reference_electron_density_ratio": density_ratio,
            "candidate_to_reference_power_per_electron_ratio":
                power_ratio / density_ratio,
        },
        "phase_space_comparisons": comparisons,
        "reduced_power_comparisons": {
            "spatially_integrated_phase_profile": metrics(
                candidate_phase, reference_phase),
            "cycle_average_spatial_profile": metrics(
                candidate_space, reference_space),
        },
        "phase_octants": grouped_means(
            candidate_phase, reference_phase,
            list(range(0, PHASES + 1, 25)), PHASES),
        "spatial_bands": grouped_means(
            candidate_space, reference_space,
            [0, 40, 80, 160, 240, 320, 360, 400], NODES),
        "finding": (
            "Cycle-volume power density is close while candidate electron "
            "density is higher; the per-electron deficit is their ratio. "
            "Field structure agrees more closely than current and local power."),
        "next_discriminator": (
            "Use independent candidate seeds to distinguish persistent "
            "current/power structure differences from four-cycle PIC noise."),
        "claim_boundary": (
            "These are descriptive comparisons to public simulation output, "
            "not experimental validation or prospectively accepted agreement."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate_root", type=Path)
    parser.add_argument("reference_raw_data", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.candidate_root.resolve(),
                     args.reference_raw_data.resolve())
    if args.output:
        atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
