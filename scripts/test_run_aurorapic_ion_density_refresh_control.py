#!/usr/bin/env python3
"""Focused tests for the held-density control runner."""

from pathlib import Path

from run_aurorapic_ion_density_refresh_control import control_deck


def main() -> None:
    rule = {
        "physics_contract": {"nodes": 400, "length_m": .025,
            "electron_timestep_s": 1e-11, "electrode_voltage_amplitude_V": 250,
            "rf_frequency_hz": 1e7, "aurorapic_initial_phase_rad": 1.5,
            "macro_weight_aurorapic_1d": 7e8, "ion_timestep_multiplier": 20},
        "locked_inputs": {"particle_state_signature": 123,
                          "source_populations": {"electrons": 10, "ions": 11}}}
    text = control_deck(rule, Path("state.aps"), 20, Path("output"))
    assert text.count("subcycle_charge_deposition = pre_push_held") == 1
    assert "timestep_multiplier = 20" in text
    print("held-density control runner tests passed")


if __name__ == "__main__":
    main()
