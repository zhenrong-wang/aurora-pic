#!/usr/bin/env python3
"""Contract checks for the native eduPIC velocity-anisotropy transform."""

from pathlib import Path
import unittest

from instrument_edupic_phase_anisotropy import instrument


class InstrumentEdupicPhaseAnisotropyTests(unittest.TestCase):
    def test_transform_is_passive_and_uses_synchronized_velocity(self) -> None:
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
        self.assertIn("measure_phase_eedf(x_e[k], mean_v, vy_e[k]",
                      transformed)
        self.assertIn("tail_longitudinal_energy_fraction", transformed)
        self.assertEqual(transformed.count("phase_eedf_tail_macro"), 3)

    def test_ambiguous_anchor_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "expected one source anchor"):
            instrument("not eduPIC")


if __name__ == "__main__":
    unittest.main()
