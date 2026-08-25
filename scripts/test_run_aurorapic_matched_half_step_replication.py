#!/usr/bin/env python3
"""Import/deck regression for the matched-half-step replication runner."""

from pathlib import Path

from run_aurorapic_matched_half_step_thresholds import measurement_deck


def main() -> None:
    base = Path(
        "tmp/edupic-argon-cross-sections-20260810/"
        "matched-heating-microstate-51949-20260820/equilibration/input.cfg"
    ).read_text(encoding="utf-8")
    rule = {
        "execution_contract": {
            "start_step": 20000, "end_step": 36000,
            "measurement_cycles": 4, "output_interval_steps": 400,
        },
        "diagnostic_contract": {
            "spatial_average_interval_steps": 2, "phase_bins": 200,
            "sampling_order": "pre_collision",
            "phase_eedf_species": "electrons", "energy_bins": 320,
            "energy_max_eV": 80.0, "tail_threshold_eV": 15.8,
            "regions": [{"name": "middle", "x_min_m": 0.005,
                         "x_max_m": 0.015}],
        },
    }
    state = {"electrons": 12, "ions": 13, "seed": 17}
    deck = measurement_deck(
        base, rule, state, Path("output"), Path("checkpoint_20000.apc"))
    for expected in (
        "steps = 36000", "spatial_average_start_step = 20001",
        "spatial_average_end_step = 36000",
        "collision_velocity_sampling = leapfrog_half_step",
        "restart_path = checkpoint_20000.apc",
        "phase_eedf_history = true",
    ):
        assert expected in deck, expected
    assert "initial_state_path" not in deck


if __name__ == "__main__":
    main()
    print("matched-half-step replication runner tests passed")
