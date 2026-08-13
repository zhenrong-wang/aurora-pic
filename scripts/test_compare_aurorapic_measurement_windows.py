#!/usr/bin/env python3
"""Small regression tests for prospective window acceptance."""

from compare_aurorapic_measurement_windows import evaluate


def main() -> None:
    thresholds = {
        "allowed_average_ionization_rate_ratio": [0.8, 1.2],
        "maximum_electron_density_integral_relative_change": 0.05,
        "maximum_ion_density_integral_relative_change": 0.05,
        "maximum_electron_energy_distribution_total_variation": 0.05,
        "maximum_electron_mean_energy_relative_change": 0.10,
        "maximum_powered_ion_energy_distribution_total_variation": 0.08,
        "maximum_grounded_ion_energy_distribution_total_variation": 0.08,
        "maximum_electrode_mean_ion_energy_relative_change": 0.10,
        "maximum_electric_field_phase_space_relative_l2": 0.05,
        "maximum_electron_current_phase_space_relative_l2": 0.10,
    }
    metrics = {
        "electron_density_integral_relative_change": 0.05,
        "ion_density_integral_relative_change": 0.05,
        "electron_energy_distribution_total_variation": 0.05,
        "electron_mean_energy_relative_change": 0.10,
        "powered_ion_energy_distribution_total_variation": 0.08,
        "grounded_ion_energy_distribution_total_variation": 0.08,
        "maximum_electrode_mean_ion_energy_relative_change": 0.10,
        "electric_field_phase_space_relative_l2": 0.05,
        "electron_current_phase_space_relative_l2": 0.10,
        "average_ionization_rate_ratio": 0.8,
    }
    assert all(evaluate(metrics, thresholds).values())
    metrics["average_ionization_rate_ratio"] = 1.2000001
    assert evaluate(metrics, thresholds)["average_ionization_rate"] is False


if __name__ == "__main__":
    main()
