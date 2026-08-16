#!/usr/bin/env python3
"""Regression tests for matched-heating seed ensemble helpers."""

import math

from analyze_aurorapic_heating_seed_ensemble import relative_range, vector_mean


def main() -> None:
    assert math.isclose(relative_range([9.0, 10.0, 11.0]), 0.2)
    assert vector_mean([[1.0, 3.0], [3.0, 5.0]]) == [2.0, 4.0]
    try:
        vector_mean([[1.0], [1.0, 2.0]])
    except ValueError:
        pass
    else:
        raise AssertionError("unequal ensemble vectors were accepted")


if __name__ == "__main__":
    main()
