#!/usr/bin/env python3
"""Small regression for the prospective heating-trend classifier."""

from analyze_edupic_heating_trend import classify, relative_change


def main() -> int:
    assert relative_change(2.0, 3.0) == 0.5
    favorable = {
        "electron_density_reference_ratio": 0.2,
        "electron_rf_power_per_particle_reference_ratio": -0.1,
        "electron_mean_energy_relative_l2": -0.1,
        "effective_ionization_frequency_reference_ratio": -0.1,
    }
    assert classify(favorable).startswith("all_predeclared")
    favorable["electron_mean_energy_relative_l2"] = 0.1
    assert classify(favorable) == "mixed_predeclared_directional_outcome"
    print("eduPIC heating-trend analysis regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
