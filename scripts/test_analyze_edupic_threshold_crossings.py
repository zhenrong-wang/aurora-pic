#!/usr/bin/env python3
"""Unit checks for unconditional threshold-crossing aggregation."""

import unittest

from analyze_edupic_threshold_crossings import COUNT_COLUMNS, aggregate


def row(phase: int, observations: int, promotions: int) -> dict[str, str]:
    result = {column: "0" for column in COUNT_COLUMNS}
    result.update({
        "phase_fraction": str(0.2 + 0.1 * phase),
        "region": "critical",
        "electron_time_macro_observations": str(observations),
        "energetic_time_macro_observations": str(observations // 4),
        "interstep_promotions": str(promotions),
        "interstep_demotions": str(promotions - 1),
        "excitation_collision_demotions": "2",
    })
    return result


class ThresholdCrossingAnalysisTests(unittest.TestCase):
    def test_counts_are_summed_before_rates_are_formed(self) -> None:
        value = aggregate(
            [row(0, 100, 4), row(1, 300, 8)],
            {"critical"}, 0.0, 1.0)
        self.assertEqual(value["electron_time_macro_observations"], 400)
        self.assertEqual(value["interstep_promotions"], 12)
        self.assertAlmostEqual(
            value["interstep_promotions_per_million_electron_steps"],
            30000.0)
        self.assertAlmostEqual(value["energetic_fraction"], 0.25)

    def test_empty_or_zero_time_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "selection is empty"):
            aggregate([row(0, 100, 4)], {"elsewhere"}, 0.0, 1.0)
        with self.assertRaisesRegex(ValueError, "no electron time"):
            aggregate([row(0, 0, 1)], {"critical"}, 0.0, 1.0)


if __name__ == "__main__":
    unittest.main()
