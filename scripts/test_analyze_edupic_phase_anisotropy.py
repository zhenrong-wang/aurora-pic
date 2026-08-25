#!/usr/bin/env python3
"""Focused tests for phase-anisotropy aggregation."""

import unittest

from analyze_edupic_phase_anisotropy import aggregate, comparison


def row(phase: float, tail: float, longitudinal: float,
        imbalance: float) -> dict[str, str]:
    return {
        "phase_fraction": str(phase), "region": "critical",
        "macro_observations": "100", "represented_observations": "200",
        "tail_represented_observations": str(tail),
        "temperature_x": "2", "temperature_y": "1",
        "temperature_z": "1",
        "tail_longitudinal_energy_fraction": str(longitudinal),
        "tail_directional_population_imbalance": str(imbalance),
        "tail_mean_velocity_x": "3",
    }


class PhaseAnisotropyAnalysisTests(unittest.TestCase):
    def test_aggregate_uses_population_weights(self) -> None:
        result = aggregate([
            row(0.2, 20.0, 0.25, -0.5),
            row(0.3, 60.0, 0.75, 0.5),
        ], {"critical"}, 0.125, 0.5)
        self.assertAlmostEqual(result["tail_fraction"], 0.2)
        self.assertAlmostEqual(
            result["tail_longitudinal_energy_fraction"], 0.625)
        self.assertAlmostEqual(
            result["tail_directional_population_imbalance"], 0.25)
        self.assertAlmostEqual(result["temperature_x_fraction"], 0.5)
        self.assertEqual(result["tail_macro_observations"], 40)

    def test_comparison_reports_signed_difference_and_range(self) -> None:
        candidate = aggregate([row(0.2, 20.0, 0.5, 0.1)],
                              {"critical"}, 0.125, 0.5)
        first = aggregate([row(0.2, 40.0, 0.4, 0.2)],
                          {"critical"}, 0.125, 0.5)
        second = aggregate([row(0.2, 60.0, 0.6, 0.4)],
                           {"critical"}, 0.125, 0.5)
        result = comparison(candidate, [first, second])
        self.assertAlmostEqual(
            result["aurorapic_minus_native"]
                  ["tail_directional_population_imbalance"], -0.2)
        self.assertGreater(
            result["native_edupic_relative_range"]["tail_fraction"], 0.0)


if __name__ == "__main__":
    unittest.main()
