#!/usr/bin/env python3
"""Focused tests for net-charge and sheath reductions."""

from analyze_edupic_net_charge_sheath import (
    integrate_profile, left_crossing, spatial_mean_product,
)


def main() -> None:
    assert abs(left_crossing([0.0, 0.5, 1.0], 0.75, 1.0) - 0.75) < 1e-15
    assert abs(integrate_profile([2.0, 2.0], 1.0, 0.2, 0.8) - 1.2) < 1e-12
    assert abs(spatial_mean_product(
        [2.0, 2.0], [3.0, 3.0], 1.0, 0.2, 0.8) - 6.0) < 1e-12
    print("net-charge sheath analyzer tests passed")


if __name__ == "__main__":
    main()
