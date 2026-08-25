#!/usr/bin/env python3
"""Contract checks for native eduPIC field-push instrumentation."""

from pathlib import Path
import unittest

from instrument_edupic_field_push_thresholds import instrument


class InstrumentEdupicFieldPushThresholdTests(unittest.TestCase):
    def test_transform_is_passive_and_attributes_the_mover(self) -> None:
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
        self.assertIn("field_push_energetic_before", transformed)
        self.assertIn("record_field_push_threshold", transformed)
        self.assertIn("x_e[k] >= 0.0 && x_e[k] <= L", transformed)
        self.assertIn(
            "edupic_phase_eedf_field_push_thresholds.csv", transformed)

    def test_ambiguous_anchor_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "expected one source anchor"):
            instrument("not eduPIC")


if __name__ == "__main__":
    unittest.main()
