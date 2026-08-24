#!/usr/bin/env python3
"""Contract checks for the native eduPIC regional phase-EEDF transform."""

from pathlib import Path
import unittest

from instrument_edupic_phase_eedf import instrument


class InstrumentEdupicPhaseEedfTests(unittest.TestCase):
    def test_transform_is_passive_and_exactly_anchored(self) -> None:
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
        self.assertNotIn("std::random_device rd", transformed)
        self.assertEqual(transformed.count("void save_phase_eedf(void)"), 1)
        self.assertEqual(transformed.count("measure_phase_eedf("), 2)
        self.assertIn("if ((t % 2) == 0)", transformed)
        self.assertIn("particle_state_mutation_added", Path(
            "scripts/instrument_edupic_phase_eedf.py").read_text())

    def test_ambiguous_anchor_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "expected one source anchor"):
            instrument("not eduPIC")


if __name__ == "__main__":
    unittest.main()
