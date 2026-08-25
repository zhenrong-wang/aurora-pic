#!/usr/bin/env python3
"""Contract checks for native eduPIC energetic-history instrumentation."""

from pathlib import Path
import unittest

from instrument_edupic_phase_history import instrument


class InstrumentEdupicPhaseHistoryTests(unittest.TestCase):
    def test_transform_preserves_rng_calls_and_tracks_particle_identity(self) -> None:
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
        self.assertIn(
            "electron_history[k] = electron_history[N_e-1]", transformed)
        self.assertIn(
            "electron_history[N_e].born_during_window = true", transformed)
        self.assertIn("return collision_type", transformed)
        self.assertIn("edupic_phase_eedf_history.csv", transformed)

    def test_ambiguous_anchor_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "expected one source anchor"):
            instrument("not eduPIC")


if __name__ == "__main__":
    unittest.main()
