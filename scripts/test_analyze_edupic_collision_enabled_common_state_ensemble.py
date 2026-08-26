#!/usr/bin/env python3
"""Focused arithmetic tests for the collision-enabled ensemble analyzer."""

from analyze_edupic_collision_enabled_common_state_ensemble import (
    summary, symmetric_relative)


def main() -> None:
    assert symmetric_relative(100.0, 100.0) == 0.0
    assert abs(symmetric_relative(95.0, 105.0) - 0.1) < 1.0e-15
    value = summary([1.0, 2.0, 3.0, 4.0, 5.0])
    assert value["mean"] == 3.0
    assert value["minimum"] == 1.0
    assert value["maximum"] == 5.0
    assert abs(value["sample_standard_deviation"] - 2.5 ** 0.5) < 1.0e-15
    print("collision-enabled common-state analyzer tests passed")


if __name__ == "__main__":
    main()
