#!/usr/bin/env python3
"""Focused tests for the published-fit argon package generator."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "generate_maiorov_2024_argon.py"
SPEC = importlib.util.spec_from_file_location("maiorov_generator", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def main() -> int:
    assert math.isclose(MODULE.elastic_angstrom2(0.0), 7.66)
    assert MODULE.elastic_angstrom2(400.0) > 0.0
    assert MODULE.threshold_fit_angstrom2(
        MODULE.EXCITATION_EV,
        MODULE.EXCITATION_EV,
        0.802,
        0.229,
        1.55,
        0.702,
    ) == 0.0
    with tempfile.TemporaryDirectory() as temporary:
        output = Path(temporary) / "package"
        subprocess.run(
            [sys.executable, str(SCRIPT), str(output), "--points", "101"],
            check=True,
            capture_output=True,
            text=True,
        )
        manifest = (output / "maiorov_2024_argon.gas").read_text()
        assert "dataset_id = maiorov.2024.argon.analytic-fits" in manifest
        for name in (
            "electron_elastic.dat",
            "electron_excitation.dat",
            "electron_ionization.dat",
        ):
            rows = [
                line
                for line in (output / name).read_text().splitlines()
                if line and not line.startswith("#")
            ]
            assert len(rows) >= 104
            assert rows[0].split()[0] == "0"
    print("Maiorov argon generator tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
