#!/usr/bin/env python3
"""Unit tests for the eduPIC normalization audit helpers."""

import math

from audit_edupic_normalization import implied_counts, trapezoid_integral


def main() -> None:
    dx = 0.025 / 399
    density = 2.8e15
    assert math.isclose(
        trapezoid_integral([density] * 400, dx), density * 0.025,
        rel_tol=1.0e-15)

    line_weight = 7.0e8
    matrix = [[density] * 200 for _ in range(400)]
    counts = implied_counts(matrix, dx, line_weight)
    expected = density * 0.025 / line_weight
    assert len(counts) == 200
    assert all(math.isclose(value, expected, rel_tol=1.0e-15)
               for value in counts)

    try:
        trapezoid_integral([1.0], dx)
    except ValueError:
        pass
    else:
        raise AssertionError("short trapezoid input must fail")


if __name__ == "__main__":
    main()
