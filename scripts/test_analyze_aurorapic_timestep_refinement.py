#!/usr/bin/env python3
"""Regression for timestep-refinement acceptance boundaries."""

from analyze_aurorapic_timestep_refinement import evaluate


def main() -> None:
    limits = {
        "maximum_electron_density_integral_relative_change": .03,
        "maximum_ion_density_integral_relative_change": .03,
        "maximum_electron_energy_distribution_total_variation": .04,
        "maximum_electron_mean_energy_relative_change": .05,
        "maximum_electric_field_phase_space_relative_l2": .03,
        "maximum_electron_current_phase_space_relative_l2": .08,
        "allowed_average_ionization_rate_ratio": [.9, 1.1],
        "maximum_powered_ion_energy_distribution_total_variation": .08,
        "maximum_grounded_ion_energy_distribution_total_variation": .08,
        "maximum_electrode_mean_ion_energy_relative_change": .1,
        "maximum_electron_power_per_particle_relative_change": .08,
        "maximum_ionization_frequency_relative_change": .1,
    }
    metrics = {
        "electron_density_integral_relative_change": .03,
        "ion_density_integral_relative_change": .03,
        "electron_energy_distribution_total_variation": .04,
        "electron_mean_energy_relative_change": .05,
        "electric_field_phase_space_relative_l2": .03,
        "electron_current_phase_space_relative_l2": .08,
        "average_ionization_rate_ratio": .9,
        "powered_ion_energy_distribution_total_variation": .08,
        "grounded_ion_energy_distribution_total_variation": .08,
        "maximum_electrode_mean_ion_energy_relative_change": .1,
        "electron_power_per_particle_relative_change": .08,
        "ionization_frequency_relative_change": .1,
    }
    assert all(evaluate(metrics, limits).values())
    metrics["electron_power_per_particle_relative_change"] = .080001
    assert evaluate(metrics, limits)["electron_power_per_particle"] is False


if __name__ == "__main__":
    main()
