#!/usr/bin/env python3
"""Contract checks for native eduPIC promotion-band work instrumentation."""

from pathlib import Path
import unittest

from instrument_edupic_promotion_band_work import instrument


class InstrumentEdupicPromotionBandWorkTests(unittest.TestCase):
    def test_transform_is_passive_and_attributes_exact_work(self) -> None:
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
        self.assertIn("energy_before_eV < 11.5", transformed)
        self.assertIn("energy_before_eV >= 15.8", transformed)
        self.assertIn("work_eV = energy_after_eV - energy_before_eV",
                      transformed)
        self.assertIn("record_promotion_band_work", transformed)
        self.assertIn("edupic_phase_eedf_promotion_band_work.csv",
                      transformed)

    def test_ambiguous_anchor_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "expected one source anchor"):
            instrument("not eduPIC")


if __name__ == "__main__":
    unittest.main()
