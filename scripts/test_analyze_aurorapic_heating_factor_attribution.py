#!/usr/bin/env python3
"""Regression tests for exact current-field heating attribution."""

import math

from analyze_aurorapic_heating_factor_attribution import shapley_current_field


def main() -> None:
    current, field = shapley_current_field(
        [3.0], [5.0], [2.0], [7.0])
    direct = 3.0 * 5.0 - 2.0 * 7.0
    assert math.isclose(current[0] + field[0], direct)
    assert math.isclose(current[0], 6.0)
    assert math.isclose(field[0], -5.0)
    try:
        shapley_current_field([1.0], [], [1.0], [1.0])
    except ValueError:
        pass
    else:
        raise AssertionError("mismatched factor vectors were accepted")


if __name__ == "__main__":
    main()
