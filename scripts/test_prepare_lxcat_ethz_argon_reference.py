#!/usr/bin/env python3
"""Conservative unit tests for the local ETHZ reference preparer."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("prepare_lxcat_ethz_argon_reference.py")
SPEC = importlib.util.spec_from_file_location("prepare_ethz", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PrepareEthzReferenceTests(unittest.TestCase):
    def test_interpolation_consensus_uncertainty_and_conversion(self) -> None:
        blocks = ["DATABASE: ETHZ (ETH Zurich, High Voltage Laboratory)"]
        for pressure in MODULE.PRESSURES_KPA:
            offset = pressure * 1.0e21
            blocks.append(
                "SPECIES: e / Ar\n"
                "PROCESS: Mobility x gas density (&mu;N)\n"
                f"COMMENT: p={pressure}kPa.\n"
                "COLUMNS: Reduced electric field (Td) | Mobility x gas density ((m.V.s)-1)\n"
                "-----------------------------\n"
                f" 5.0 {9.0e23 + offset}\n"
                f" 15.0 {10.0e23 + offset}\n"
                f" 25.0 {11.0e23 + offset}\n"
                f" 35.0 {12.0e23 + offset}\n"
                "-----------------------------\n"
            )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_path = root / "source.txt"
            source_path.write_text("\n".join(blocks), encoding="utf-8")
            datasets = MODULE.parse_argon_mobility(source_path)
            rows = MODULE.consensus_rows(datasets)
            self.assertEqual([row[0] for row in rows], [10.0, 20.0, 30.0])
            self.assertAlmostEqual(rows[0][1], 9562.0)
            self.assertAlmostEqual(rows[0][2], 28.635642126552707)
            self.assertAlmostEqual(rows[0][3], 9530.0)
            self.assertAlmostEqual(rows[0][4], 9600.0)

            MODULE.write_reference(root / "out", rows)
            with (root / "out/ethz-argon-five-pressure-consensus.csv").open(
                newline="", encoding="utf-8"
            ) as stream:
                output = list(csv.DictReader(stream))
            self.assertEqual(output[0]["pressure_series_count"], "5")
            self.assertEqual(output[0]["reduced_field_td"], "10")
            manifest = (
                root / "out/ethz-argon-five-pressure-consensus.swarm-reference"
            ).read_text(encoding="utf-8")
            self.assertIn("reference_uncertainty_column", manifest)
            self.assertIn("relative_tolerance = 0.10", manifest)

    def test_rejects_missing_pressure_series_and_extrapolation(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "outside a pressure series"):
            MODULE.interpolate([(10.0, 1.0), (20.0, 2.0)], 30.0)


if __name__ == "__main__":
    unittest.main()
