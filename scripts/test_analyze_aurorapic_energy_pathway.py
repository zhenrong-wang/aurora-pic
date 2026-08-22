#!/usr/bin/env python3
"""Focused tests for the phase-resolved energy-pathway audit."""

import math

from analyze_aurorapic_energy_pathway import NODES, PHASES, periodic_derivative


def main() -> None:
    constant = periodic_derivative([3.0] * (PHASES * NODES))
    assert all(value == 0.0 for value in constant)
    phase_values = []
    for phase in range(PHASES):
        phase_values.extend([math.sin(2.0 * math.pi * phase / PHASES)] * NODES)
    derivative = periodic_derivative(phase_values)
    assert derivative[0] > 0.0
    assert abs(derivative[50 * NODES]) < 1e-6
    try:
        periodic_derivative([1.0])
    except ValueError:
        pass
    else:
        raise AssertionError("invalid energy-density shape was accepted")


if __name__ == "__main__":
    main()
    print("energy-pathway analyzer tests passed")
