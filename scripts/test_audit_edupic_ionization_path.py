#!/usr/bin/env python3
"""Focused tests for the eduPIC ionization-path audit."""

from pathlib import Path
import tempfile

from audit_edupic_ionization_path import cumulative_delta, read_cross_section


def test_cross_section_grid() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        path = Path(temporary) / "ionization.dat"
        path.write_text(
            "# energy cross section\n0 0\n0.001 1e-20\n0.002 2e-20\n",
            encoding="utf-8")
        energies, values = read_cross_section(path)
        assert energies == [0.0, .001, .002]
        assert values == [0.0, 1e-20, 2e-20]


def test_counter_delta() -> None:
    data = [{"ionization": "41"}, {"ionization": "58"}]
    assert cumulative_delta(data, "ionization") == 17


if __name__ == "__main__":
    test_cross_section_grid()
    test_counter_delta()
    print("ionization-path audit tests passed")
