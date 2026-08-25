#!/usr/bin/env python3
"""Contract tests for native eduPIC cycle-history instrumentation."""

from pathlib import Path
import unittest

from instrument_edupic_cycle_history import instrument


class InstrumentEdupicCycleHistoryTests(unittest.TestCase):
    def test_transform_is_passive_and_cycle_resolved(self) -> None:
        path = Path("tmp/edupic-upstream-review-20260804/C/eduPIC.cc")
        if not path.is_file():
            self.skipTest("local GPL eduPIC source is unavailable")
        source = path.read_text(encoding="utf-8").replace(
            "    test_cross_sections(); return 1;",
            "    //test_cross_sections(); return 1;",
        )
        transformed = instrument(source)
        self.assertEqual(source.count("R01(MTgen)"),
                         transformed.count("R01(MTgen)"))
        self.assertEqual(source.count("RMB(MTgen)"),
                         transformed.count("RMB(MTgen)"))
        self.assertIn("edupic_cycle_history.csv", transformed)
        self.assertIn("if (measurement_mode) save_cycle_history();",
                      transformed)
        self.assertIn("energetic_electrons++", transformed)
        self.assertIn("cycle_electron_ionization_collisions++", transformed)
        self.assertIn("save_threshold_crossings();", transformed)
        self.assertNotIn("std::random_device rd", transformed)

    def test_ambiguous_anchor_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "expected one source anchor"):
            instrument("not eduPIC")


if __name__ == "__main__":
    unittest.main()
