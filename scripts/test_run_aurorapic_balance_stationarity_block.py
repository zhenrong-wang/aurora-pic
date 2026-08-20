#!/usr/bin/env python3
"""Focused tests for the bounded source/loss-stationarity runner."""

import csv
import json
from pathlib import Path
import tempfile

from run_aurorapic_balance_stationarity_block import (
    analyze_output, build_deck, normalized_slope,
)


def write_csv(path: Path, fields: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(rows)


def test_deck_and_slope() -> None:
    base = """steps = 1
output_interval = 1
output_dir = old
spatial_average = true
phase_eedf = true
wall_impact_spectrum = true
checkpoint_interval = 1
runtime_backend = serial
runtime_threads = 1
restart_path = old.apc
[species.electrons]
particles = 1
"""
    deck = build_deck(base, Path("new-output"), Path("state.apc"), 120000, 400)
    for expected in ("steps = 120000", "output_interval = 400",
                     "spatial_average = false", "phase_eedf = false",
                     "wall_impact_spectrum = false",
                     "checkpoint_interval = 120000",
                     "restart_path = state.apc"):
        assert expected in deck
    assert abs(normalized_slope([100.0, 101.0, 102.0]) - 1.0 / 101.0) < 1e-15


def test_exact_ledger_analysis() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        output = Path(temporary)
        scalar_fields = ["step", "time", "kinetic_energy", "field_energy",
                         "total_energy", "live_particles",
                         "live_particles_electrons", "live_particles_ions"]
        write_csv(output / "scalars.csv", scalar_fields, [
            [0, 0, 2, 1, 3, 2000, 1000, 1000],
            [4000, 1, 2, 1, 3, 2002, 1002, 1000],
            [8000, 2, 2, 1, 3, 2004, 1004, 1000],
        ])
        collision_fields = [
            "step", "cumulative_collisions_electron_mcc.ionization"]
        write_csv(output / "collisions.csv", collision_fields, [
            [0, 0], [4000, 100], [8000, 200]])
        boundary_fields = [
            "step", "absorbed_left_count_electrons",
            "absorbed_right_count_electrons", "absorbed_left_count_ions",
            "absorbed_right_count_ions"]
        write_csv(output / "boundary_losses.csv", boundary_fields, [
            [0, 0, 0, 0, 0], [8000, 98, 98, 100, 100]])
        write_csv(output / "fields_0.csv", ["E"], [[-2], [3]])
        (output / "checkpoint_8000.apc").write_bytes(b"checkpoint")
        rule = {
            "prospective_block_stationarity": {
                "maximum_absolute_electron_source_loss_relative_imbalance": .03,
                "maximum_absolute_ion_source_loss_relative_imbalance": .03,
                "maximum_absolute_normalized_electron_population_slope_per_cycle": .01,
                "maximum_absolute_normalized_ion_population_slope_per_cycle": .01,
                "maximum_absolute_normalized_field_energy_slope_per_cycle": .01,
                "maximum_ionization_count_coefficient_of_variation": .10,
            },
            "execution_contract": {
                "maximum_total_particles": 500000,
                "maximum_absolute_field_V_m": 1e7,
                "maximum_peak_resident_set_kib": 307200,
            },
        }
        result = analyze_output(
            output, 0, 8000, 4000, rule, {"peak_resident_set_kib": 100})
        assert result["particle_ledger"]["exact_species_closure"] is True
        assert result["metrics"]["electron_source_loss_relative_imbalance"] == .02
        assert result["stationarity_block_passed"] is True
        assert result["all_hard_safety_gates_passed"] is True


if __name__ == "__main__":
    test_deck_and_slope()
    test_exact_ledger_analysis()
    print("balance-stationarity runner tests passed")
