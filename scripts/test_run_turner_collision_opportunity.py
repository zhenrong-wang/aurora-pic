#!/usr/bin/env python3
"""Focused deck test for the Turner collision-opportunity runner."""

from pathlib import Path

from run_turner_collision_opportunity import deck


def main() -> None:
    rule = {
        "physics_contract": {
            "nodes": 129, "length_m": .067, "timestep_s": 1e-10,
            "voltage_amplitude_V": 450, "rf_frequency_hz": 1e7,
            "neutral_density_m3": 1e20, "neutral_temperature_K": 300,
            "electron_max_frequency_s": 1e8, "ion_max_frequency_s": 1e7,
            "macro_weight_m2": 2e8},
        "ensemble_contract": {"steps": 400, "measurement_cycles": 1},
        "locked_inputs": {"particle_state_signature": 123,
                          "source_populations": {"electrons": 10, "ions": 11}}}
    text = deck(rule, Path("state.aps"), Path("e.gas"), Path("i.gas"),
                Path("output"), 13507, "single_bernoulli")
    assert text.count("opportunity_sampling = single_bernoulli") == 2
    assert "initial_state_signature = 123" in text
    assert "particles = 10" in text and "particles = 11" in text
    print("Turner collision-opportunity runner tests passed")


if __name__ == "__main__":
    main()
