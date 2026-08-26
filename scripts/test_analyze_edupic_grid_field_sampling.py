#!/usr/bin/env python3
"""Focused tests for grid-field versus particle-sampling reductions."""

from analyze_edupic_grid_field_sampling import (
    boundary_value, spatial_mean_square,
)


def main() -> None:
    assert boundary_value([0.0, 2.0], 0.25, 1.0) == 0.5
    assert abs(spatial_mean_square(
        [3.0, 3.0], 1.0, 0.2, 0.8) - 9.0) < 1e-12
    value = spatial_mean_square([0.0, 1.0, 2.0], 1.0, 0.25, 0.75)
    assert abs(value - 1.125) < 1e-15
    print("grid-field sampling analyzer tests passed")


if __name__ == "__main__":
    main()
