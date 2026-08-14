#!/usr/bin/env python3
"""Bounded regression tests for region-matched collision-audit math."""

import math

from analyze_aurorapic_region_matched_collisions import (
    integrate_region, interpolate,
)


def main() -> None:
    x = [0.0, 1.0, 2.0]
    values = [0.0, 2.0, 4.0]
    assert interpolate(x, values, 0.5) == 1.0
    assert math.isclose(integrate_region(x, values, 0.5, 1.5), 2.0)
    assert math.isclose(integrate_region(x, [3.0, 3.0, 3.0], 0.2, 1.8), 4.8)
    try:
        integrate_region(x, values, -0.1, 1.0)
    except ValueError:
        pass
    else:
        raise AssertionError("out-of-grid regional bound was accepted")


if __name__ == "__main__":
    main()
