#!/usr/bin/env python3
"""Contract checks for native eduPIC threshold-crossing instrumentation."""

from pathlib import Path
import unittest

from instrument_edupic_threshold_crossings import instrument


class InstrumentEdupicThresholdCrossingTests(unittest.TestCase):
    def test_transform_is_passive_and_unconditional(self) -> None:
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
        self.assertIn("threshold_electron_time[phase][region]++", transformed)
        self.assertIn("record_threshold_collision", transformed)
        self.assertIn("record_threshold_birth", transformed)
        self.assertIn("edupic_phase_eedf_threshold_crossings.csv", transformed)

    def test_ambiguous_anchor_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "expected one source anchor"):
            instrument("not eduPIC")


if __name__ == "__main__":
    unittest.main()
