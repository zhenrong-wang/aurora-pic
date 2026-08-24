#!/usr/bin/env python3
"""Conservative unit tests for the local Dutton reference preparer."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("prepare_lxcat_dutton_argon_reference.py")
SPEC = importlib.util.spec_from_file_location("prepare_dutton", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PrepareDuttonReferenceTests(unittest.TestCase):
    def test_selection_and_mobility_conversion(self) -> None:
        blocks = []
        fixtures = {
            "Brambring 1964.": [(52.8, 7.5758e23), (99.5, 6.3317e23)],
            "Jager et al 1962.": [(53.6, 8.0224e23), (90.8, 6.6079e23)],
            "Wagner 1964.": [(73.6, 7.6087e23), (99.7, 7.4223e23)],
        }
        for source, rows in fixtures.items():
            data = "\n".join(f" {field:.6e} {mobility:.6e}" for field, mobility in rows)
            blocks.append(
                "SPECIES: e / Ar\n"
                "PROCESS: Mobility x gas density (&mu;N)\n"
                f"COMMENT: {source}\n"
                "COLUMNS: Reduced electric field (Td) | Mobility x gas density ((m.V.s)-1)\n"
                "-----------------------------\n"
                f"{data}\n"
                "-----------------------------\n"
            )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_path = root / "source.txt"
            source_path.write_text("\n".join(blocks), encoding="utf-8")
            selected = MODULE.select_points(MODULE.parse_argon_mobility(source_path))
            self.assertEqual([row[0] for row in selected], [52.8, 53.6, 90.8, 99.5, 99.7])
            MODULE.write_reference(root / "out", selected)
            with (root / "out/dutton-argon-selected-reference.csv").open(
                newline="", encoding="utf-8"
            ) as stream:
                rows = list(csv.DictReader(stream))
            self.assertAlmostEqual(float(rows[0]["drift_velocity_m_s"]), 40000.224)
            self.assertEqual(rows[0]["source"], "Brambring 1964")


if __name__ == "__main__":
    unittest.main()
