#!/usr/bin/env python3
"""Regression tests for EEDF energy-band integration."""

import math

from analyze_edupic_eedf_discrepancy import integrate, normalize


def main() -> None:
    table = normalize([(0.0, 1.0, 1.0), (1.0, 3.0, 1.0)])
    assert math.isclose(integrate(table, 0.0, 3.0), 1.0)
    assert math.isclose(integrate(table, 0.0, 3.0, 1), 1.25)
    assert math.isclose(integrate(table, 0.5, 2.0), 0.5)
    assert math.isclose(integrate(table, 0.5, 2.0, 1), 0.5625)


if __name__ == "__main__":
    main()
