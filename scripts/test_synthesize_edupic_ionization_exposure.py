#!/usr/bin/env python3
"""Focused tests for the eduPIC ionization-exposure synthesis."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile

from synthesize_edupic_ionization_exposure import (
    IonizationExposureError, sha256, synthesize,
)


def write(path: Path, value: dict[str, object]) -> None:
    path.write_text(json.dumps(value) + "\n", encoding="utf-8")


def main() -> int:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        phase = root / "phase.json"
        cycle = root / "cycle.json"
        partition = root / "partition.json"
        base = {"case_id": "edupic-1.0-default-argon-ccp"}
        write(phase, base | {
            "all_measurement_gates_passed": True,
            "comparisons": {
                "critical_x020_to_x060_phase_0p125_to_0p5": {
                    "aurorapic_to_native_edupic_ratio": {
                        "eedf_folded_ionization_frequency_s-1": 0.85}}},
        })
        write(cycle, base | {
            "all_passivity_balance_and_repeatability_gates_passed": True,
            "normalized_four_cycle_comparison": {
                "aurorapic_to_native_ionization_ratio": 0.87,
                "aurorapic_to_native_electron_wall_loss_ratio": 0.97},
        })
        write(partition, base | {
            "aurorapic_to_native_ratios": {
                "ionizing_population_fraction": 0.94,
                "ionizing_band_energy_per_area": 1.02},
        })
        paths = {
            "phase_eedf": phase,
            "cycle_history": cycle,
            "checkpoint_partition": partition,
        }
        expected = {name: sha256(path) for name, path in paths.items()}
        result = synthesize(phase, cycle, partition, expected)
        ratios = result["ratios_aurorapic_to_native_edupic"]
        assert abs(ratios[
            "absolute_folded_to_realized_ratio_difference"] - 0.02) < 1e-12
        assert result["assessment"][
            "phase_space_ionizing_exposure_deficit_supported"]

        write(cycle, {"case_id": "different-case"})
        try:
            synthesize(phase, cycle, partition, expected)
        except IonizationExposureError:
            pass
        else:
            raise AssertionError("modified locked input was accepted")
    print("ionization-exposure synthesis tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
