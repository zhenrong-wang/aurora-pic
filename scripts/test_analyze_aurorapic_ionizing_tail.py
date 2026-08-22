#!/usr/bin/env python3
"""Focused tests for ionizing-tail analysis helpers."""

import math

from analyze_aurorapic_ionizing_tail import (
    kernel_frequency, region_phase_sum, selected_nodes,
)


def main() -> None:
    assert kernel_frequency(0.0, [0.0, .001], [0.0, 1.0]) == 0.0
    assert selected_nodes(0.005, 0.01)[0] > 0
    values = [1.0] * (200 * 400)
    nodes = [0, 1, 2]
    assert region_phase_sum(values, nodes, 0, 1) == 2.0
    try:
        region_phase_sum([1.0], nodes)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid phase-space shape was accepted")
    assert math.isfinite(kernel_frequency(
        .001, [0.0, .001], [0.0, 1e-20]))


if __name__ == "__main__":
    main()
    print("ionizing-tail analyzer tests passed")
