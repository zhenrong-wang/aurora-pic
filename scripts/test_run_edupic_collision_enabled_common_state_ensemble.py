#!/usr/bin/env python3
"""Focused deck test for the collision-enabled common-state runner."""

from pathlib import Path

from run_edupic_collision_enabled_common_state_ensemble import deck


def main() -> None:
    rule = {
        "physics_contract": {
            "nodes": 400, "length_m": 0.025,
            "electron_timestep_s": 1.0e-11,
            "electrode_voltage_amplitude_V": 250.0,
            "rf_frequency_hz": 13.56e6,
            "aurorapic_initial_phase_rad": 1.57,
            "neutral_density_m3": 2.0e21, "neutral_temperature_K": 350.0,
            "macro_weight_aurorapic_1d": 7.0e8,
            "ion_timestep_multiplier": 20,
        },
        "locked_inputs": {
            "aurorapic_particle_state_signature": 42,
            "initial_populations": {"electrons": 10, "ions": 11},
        },
        "ensemble_contract": {"electron_pushes": 3999},
    }
    value = deck(rule, Path("state.aps"), 13507, Path("e.gas"),
                 Path("i.gas"), Path("output"))
    assert "steps = 3999" in value
    assert "seed = 13507" in value
    assert "subcycle_charge_deposition = pre_push_held" in value
    assert value.count("model = null_collision") == 2
    assert "gas_data_file = e.gas" in value
    assert "gas_data_file = i.gas" in value
    print("collision-enabled common-state runner tests passed")


if __name__ == "__main__":
    main()
