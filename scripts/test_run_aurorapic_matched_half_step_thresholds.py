#!/usr/bin/env python3
"""Focused deck tests for the matched-half-step threshold runner."""

from pathlib import Path

from run_aurorapic_matched_half_step_thresholds import (
    equilibration_deck, measurement_deck,
)


def main() -> None:
    base = Path(
        "tmp/edupic-argon-cross-sections-20260810/"
        "matched-heating-microstate-51949-20260820/equilibration/input.cfg"
    ).read_text(encoding="utf-8")
    rule = {
        "execution_contract": {
            "equilibration_steps": 4000,
            "equilibration_cycles": 1,
            "equilibration_phase_bins": 16,
            "equilibration_sampling_interval_steps": 10,
            "start_step": 4000,
            "end_step": 20000,
            "measurement_cycles": 4,
            "output_interval_steps": 400,
        },
        "diagnostic_contract": {
            "spatial_average_interval_steps": 2,
            "phase_bins": 200,
            "sampling_order": "pre_collision",
            "phase_eedf_species": "electrons",
            "energy_bins": 320,
            "energy_max_eV": 80.0,
            "tail_threshold_eV": 15.8,
            "regions": [{"name": "left", "x_min_m": 0.0,
                         "x_max_m": 0.01}],
        },
    }
    state = {
        "electrons": 12, "ions": 13, "seed": 17,
        "particle_state_signature": 12345,
    }
    initial = equilibration_deck(
        base, rule, state, Path("equil-output"), Path("state.aps"))
    measured = measurement_deck(
        base, rule, state, Path("measurement-output"),
        Path("checkpoint_4000.apc"))
    for expected in (
        "collision_velocity_sampling = leapfrog_half_step",
        "steps = 4000", "particles = 12", "particles = 13",
        "initial_state_path = state.aps",
        "initial_state_signature = 12345",
    ):
        assert expected in initial, expected
    for expected in (
        "collision_velocity_sampling = leapfrog_half_step",
        "steps = 20000", "spatial_average_start_step = 4001",
        "spatial_average_end_step = 20000",
        "phase_eedf_history = true", "phase_eedf_tail_threshold = 15.8",
        "phase_eedf_regions = left:0.0:0.01",
        "restart_path = checkpoint_4000.apc",
    ):
        assert expected in measured, expected
    assert "initial_state_path" not in measured
    assert "initial_state_signature" not in measured
    assert initial.count("collision_velocity_sampling =") == 1
    assert measured.count("collision_velocity_sampling =") == 1


if __name__ == "__main__":
    main()
    print("matched-half-step threshold runner tests passed")
