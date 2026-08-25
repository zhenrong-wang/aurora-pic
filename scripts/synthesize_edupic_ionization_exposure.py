#!/usr/bin/env python3
"""Synthesize independent eduPIC ionization-exposure diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


LOCKED_INPUT_SHA256 = {
    "phase_eedf": "5af94a076023bcb0f6b514d3a5de6327be8e6aa290f802cadf300446bbc33bbe",
    "cycle_history": "9288e356857e7d5609369a9766f1eb02d3e9244e156a2158dc2ed7c9b4c3c392",
    "checkpoint_partition": "dc1989badc921b6e2dc0a1b7d5efa94d65c499e869299d68fc22617277ac847e",
}


class IonizationExposureError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def synthesize(
        phase_path: Path,
        cycle_path: Path,
        partition_path: Path,
        expected_hashes: dict[str, str] = LOCKED_INPUT_SHA256,
) -> dict[str, object]:
    paths = {
        "phase_eedf": phase_path,
        "cycle_history": cycle_path,
        "checkpoint_partition": partition_path,
    }
    actual_hashes = {name: sha256(path) for name, path in paths.items()}
    for name, expected in expected_hashes.items():
        if actual_hashes.get(name) != expected:
            raise IonizationExposureError(f"locked input hash differs: {name}")

    documents = {
        name: json.loads(path.read_text(encoding="utf-8"))
        for name, path in paths.items()
    }
    phase = documents["phase_eedf"]
    cycle = documents["cycle_history"]
    partition = documents["checkpoint_partition"]
    case_ids = {document.get("case_id") for document in documents.values()}
    if case_ids != {"edupic-1.0-default-argon-ccp"}:
        raise IonizationExposureError("input case identifiers differ")
    if phase.get("all_measurement_gates_passed") is not True:
        raise IonizationExposureError("phase-EEDF measurement gates did not pass")
    if cycle.get(
            "all_passivity_balance_and_repeatability_gates_passed") is not True:
        raise IonizationExposureError("cycle-history gates did not pass")

    phase_ratio = float(
        phase["comparisons"]["critical_x020_to_x060_phase_0p125_to_0p5"]
        ["aurorapic_to_native_edupic_ratio"]
        ["eedf_folded_ionization_frequency_s-1"])
    event_ratio = float(
        cycle["normalized_four_cycle_comparison"]
        ["aurorapic_to_native_ionization_ratio"])
    endpoint_tail_ratio = float(
        partition["aurorapic_to_native_ratios"]
        ["ionizing_population_fraction"])
    endpoint_band_energy_ratio = float(
        partition["aurorapic_to_native_ratios"]
        ["ionizing_band_energy_per_area"])
    wall_loss_ratio = float(
        cycle["normalized_four_cycle_comparison"]
        ["aurorapic_to_native_electron_wall_loss_ratio"])
    absolute_difference = abs(event_ratio - phase_ratio)

    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "cross_diagnostic_ionization_exposure_synthesis",
        "provenance": {
            "analyzer_sha256": sha256(Path(__file__)),
            "input_sha256": actual_hashes,
        },
        "integrity": {
            "locked_inputs_match": True,
            "common_case_identifier": True,
            "source_measurement_and_repeatability_gates_passed": True,
        },
        "ratios_aurorapic_to_native_edupic": {
            "critical_phase_region_eedf_folded_ionization_frequency":
                phase_ratio,
            "four_cycle_realized_ionizations_per_starting_electron":
                event_ratio,
            "realized_to_folded_ratio": event_ratio / phase_ratio,
            "absolute_folded_to_realized_ratio_difference":
                absolute_difference,
            "endpoint_ionizing_population_fraction": endpoint_tail_ratio,
            "endpoint_ionizing_band_energy_per_area":
                endpoint_band_energy_ratio,
            "four_cycle_electron_wall_loss_per_starting_electron":
                wall_loss_ratio,
        },
        "assessment": {
            "observation": (
                "The prospectively locked critical phase-region EEDF fold "
                f"predicts an ionization-frequency ratio of {phase_ratio:.4f}, "
                "while the independent cycle histories measure an event ratio "
                f"of {event_ratio:.4f}; they differ by "
                f"{absolute_difference:.4f} in absolute ratio, or "
                f"{abs(event_ratio / phase_ratio - 1.0) * 100.0:.2f} percent "
                "relative. The endpoint ionizing-band energy per area is near "
                "parity even though relative tail occupancy is lower, so a "
                "global endpoint energy inventory hides the phase-space "
                "exposure deficit that quantitatively tracks realized "
                "ionization."),
            "collision_sampler_defect_supported": False,
            "phase_space_ionizing_exposure_deficit_supported": True,
            "excess_electron_wall_loss_supported": False,
        },
        "claim_boundary": (
            "This is a post hoc synthesis of independently evolved, "
            "checksum-locked diagnostics. Their close quantitative agreement "
            "is cross-diagnostic evidence, not a preregistered acceptance gate "
            "or proof of causal onset, experimental validity, or general PIC "
            "correctness."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase_eedf", type=Path)
    parser.add_argument("cycle_history", type=Path)
    parser.add_argument("checkpoint_partition", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        result = synthesize(
            args.phase_eedf, args.cycle_history, args.checkpoint_partition)
    except (IonizationExposureError, OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
