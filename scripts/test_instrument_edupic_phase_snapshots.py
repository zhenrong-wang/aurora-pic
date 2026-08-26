#!/usr/bin/env python3
"""Contract tests for passive native phase-grid snapshots."""

from pathlib import Path
import unittest

from instrument_edupic_phase_snapshots import instrument


class InstrumentEdupicPhaseSnapshotsTests(unittest.TestCase):
    def test_transform_is_passive_and_phase_aligned(self) -> None:
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
        self.assertIn("((t + 1) % 400) == 0", transformed)
        self.assertIn("save_phase_snapshot(t + 1, rho)", transformed)
        self.assertIn("edupic_phase_snapshots.csv", transformed)
        self.assertIn("save_mover_decomposition();", transformed)

    def test_ambiguous_anchor_is_rejected(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "expected one source anchor"):
            instrument("not eduPIC")


if __name__ == "__main__":
    unittest.main()
