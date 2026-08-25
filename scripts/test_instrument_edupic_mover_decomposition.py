#!/usr/bin/env python3
"""Contract checks for native eduPIC mover-decomposition instrumentation."""

from pathlib import Path
import unittest

from instrument_edupic_mover_decomposition import instrument


class InstrumentEdupicMoverDecompositionTests(unittest.TestCase):
    def test_transform_is_passive_and_exactly_decomposes_work(self) -> None:
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
        self.assertIn("E_MASS * origin_vx * delta_vx / EV_TO_J",
                      transformed)
        self.assertIn("0.5 * E_MASS * delta_vx * delta_vx / EV_TO_J",
                      transformed)
        self.assertIn("vx_e[k] - field_push_vx_before", transformed)
        self.assertIn("edupic_phase_eedf_mover_decomposition.csv",
                      transformed)

    def test_ambiguous_anchor_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "expected one source anchor"):
            instrument("not eduPIC")


if __name__ == "__main__":
    unittest.main()
