#!/usr/bin/env python3
"""Synthetic contract checks for the eduPIC surface-flux instrumenter."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


SCRIPT = Path(__file__).with_name("instrument_edupic_surface_flux.py")
SPEC = importlib.util.spec_from_file_location("instrument_edupic", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class InstrumentEdupicSurfaceFluxTests(unittest.TestCase):
    def test_pinned_source_transforms_once_without_rng_changes(self) -> None:
        source_path = Path("tmp/edupic-upstream-review-20260804/C/eduPIC.cc")
        if not source_path.is_file():
            self.skipTest("local GPL eduPIC source is unavailable")
        source = source_path.read_text(encoding="utf-8").replace(
            "    test_cross_sections(); return 1;",
            "    //test_cross_sections(); return 1;",
        )
        transformed = MODULE.instrument(source)
        self.assertEqual(source.count("R01(MTgen)"), transformed.count("R01(MTgen)"))
        self.assertEqual(source.count("RMB(MTgen)"), transformed.count("RMB(MTgen)"))
        self.assertNotIn("std::random_device rd", transformed)
        self.assertIn("EDUPIC_DIAGNOSTIC_SEED", transformed)
        self.assertEqual(transformed.count("void save_surface_flux(void)"), 1)
        self.assertEqual(
            transformed.count("measure_electron_surface_crossings("), 2
        )
        self.assertIn("represented_per_macro = WEIGHT / ELECTRODE_AREA", transformed)
        self.assertIn("old_x, x_e[k], vx_e[k]", transformed)
        self.assertIn("eduPIC C implementation", MODULE.__doc__)

    def test_rejects_ambiguous_anchor(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "expected one source anchor"):
            MODULE.replace_once("x x", "x", "y", "test")


if __name__ == "__main__":
    unittest.main()
