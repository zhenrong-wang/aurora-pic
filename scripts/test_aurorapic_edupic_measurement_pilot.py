#!/usr/bin/env python3
"""Small regressions for AuroraPIC fresh-window comparison math."""

from compare_aurorapic_edupic_measurement_pilot import (
    interpolate, relative_l2, total_variation,
)


def main() -> int:
    assert interpolate([0.0, 1.0], [2.0, 4.0], 0.25) == 2.5
    assert relative_l2([1.0, 2.0], [1.0, 2.0]) == 0.0
    a = [(0.0, 1.0, 1.0)]
    b = [(0.0, 0.5, 1.0)]
    assert abs(total_variation(a, a)) < 1e-15
    assert abs(total_variation(a, b) - 0.5) < 1e-15
    print("AuroraPIC eduPIC measurement-pilot comparison regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
