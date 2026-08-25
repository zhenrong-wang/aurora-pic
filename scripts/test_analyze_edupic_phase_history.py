#!/usr/bin/env python3
"""Unit checks for energetic-particle history aggregation."""

import unittest

from analyze_edupic_phase_history import aggregate


def row(phase: int, tail: float, age: float,
        tail_macro: int | None = None) -> dict[str, str]:
    result = {
        "phase_bin": str(phase),
        "phase_fraction": str(0.2 + 0.1 * phase),
        "region": "x020_040",
        "macro_observations": "10",
        "represented_observations": "20",
        "tail_represented_observations": str(tail),
        "overflow_fraction": "0",
        "tail_mean_age_steps": str(age),
        "tail_mean_energetic_steps": "4",
        "tail_mean_energetic_duty_fraction": "0.4",
        "tail_mean_consecutive_energetic_steps": "2",
        "tail_mean_entries": "1",
        "tail_mean_elastic_collisions": "3",
        "tail_mean_excitation_collisions": "0.2",
        "tail_mean_ionization_collisions": "0.1",
        "tail_mean_charge_exchange_collisions": "0",
        "tail_mean_bgk_collisions": "0",
        "tail_born_during_window_fraction": "0.05",
    }
    if tail_macro is not None:
        result["tail_macro_observations"] = str(tail_macro)
    return result


class PhaseHistoryAnalysisTests(unittest.TestCase):
    def test_represented_and_macro_weighting_are_equivalent(self) -> None:
        candidate = aggregate(
            [row(0, 8, 10), row(1, 6, 20)], {"x020_040"}, 0.0, 1.0)
        native = aggregate(
            [row(0, 8, 10, 4), row(1, 6, 20, 3)],
            {"x020_040"}, 0.0, 1.0)
        self.assertEqual(candidate["tail_macro_observations"], 7)
        self.assertEqual(native["tail_macro_observations"], 7)
        self.assertAlmostEqual(candidate["age_steps"], native["age_steps"])
        self.assertAlmostEqual(candidate["age_steps"], 100.0 / 7.0)
        self.assertAlmostEqual(
            candidate["elastic_collisions_per_1000_age_steps"], 210.0)

    def test_empty_or_zero_tail_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "selection is empty"):
            aggregate([row(0, 8, 10)], {"elsewhere"}, 0.0, 1.0)
        with self.assertRaisesRegex(ValueError, "no energetic tail"):
            aggregate([row(0, 0, 10)], {"x020_040"}, 0.0, 1.0)


if __name__ == "__main__":
    unittest.main()
