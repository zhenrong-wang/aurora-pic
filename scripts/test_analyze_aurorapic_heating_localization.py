#!/usr/bin/env python3
"""Regression tests for matched heating localization helpers."""

import math

from analyze_aurorapic_heating_localization import grouped_means, reduced_profiles


def main() -> None:
    phase, space = reduced_profiles(
        [1.0, 3.0, 5.0, 2.0, 4.0, 6.0], phases=2, nodes=3)
    assert phase == [3.0, 4.0]
    assert space == [1.5, 3.5, 5.5]
    groups = grouped_means([1.0, 3.0], [2.0, 2.0], [0, 1, 2], 2)
    assert groups[0]["candidate_minus_reference_W_m-3"] == -1.0
    assert math.isclose(groups[1]["lower_fraction"], 0.5)
    try:
        grouped_means([1.0], [1.0], [1, 0], 1)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid grouped-mean boundaries were accepted")


if __name__ == "__main__":
    main()
