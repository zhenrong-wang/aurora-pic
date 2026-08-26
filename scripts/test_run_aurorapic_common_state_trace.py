#!/usr/bin/env python3
"""Focused tests for the common-state horizon runner."""

from run_aurorapic_common_state_trace import deck


def main() -> None:
    rule = {
        "physics_contract": {"nodes": 400, "length_m": .025,
            "electron_timestep_s": 1e-11, "electrode_voltage_amplitude_V": 250,
            "rf_frequency_hz": 1e7, "aurorapic_initial_phase_rad": 1.5,
            "macro_weight_aurorapic_1d": 7e8, "ion_timestep_multiplier": 20},
        "locked_inputs": {"particle_state_signature": 123,
                          "source_populations": {"electrons": 10, "ions": 11}}}
    text = deck(rule, __import__("pathlib").Path("state.aps"), 20,
                __import__("pathlib").Path("output"))
    assert "steps = 20" in text and "output_interval = 20" in text
    assert "initial_state_signature = 123" in text
    assert "[collisions" not in text
    print("common-state horizon runner tests passed")


if __name__ == "__main__":
    main()
