#!/usr/bin/env python3
"""Focused tests for spatial ionization localization helpers."""

import math

from analyze_edupic_spatial_ionization import (
    band_summary, phase_band_summary, phase_profile, validate_phase_grid,
    weighted_sum,
)


def test_weighted_sum() -> None:
    # Two phases, three nodes: endpoint weights give 2 spatial intervals.
    assert weighted_sum([1.0] * 6, phases=2, nodes=3) == 4.0
    assert weighted_sum([1.0] * 6, 0, 1, phases=2, nodes=3) == 1.0


def test_additive_band_attribution() -> None:
    reference_rate = [10.0, 10.0, 10.0]
    candidate_rate = [8.0, 10.0, 9.0]
    density = [2.0, 2.0, 2.0]
    power = [4.0, 4.0, 4.0]
    bands = band_summary(
        candidate_rate, reference_rate, density, density, power, power,
        boundaries=[0, 1, 3], phases=1, nodes=3)
    assert math.isclose(
        sum(float(band["fraction_of_net_source_gap"]) for band in bands),
        1.0)
    assert bands[0]["candidate_to_reference_effective_ionization_frequency_ratio"] == .8
    assert bands[1]["candidate_to_reference_electron_density_ratio"] == 1.0


def test_phase_attribution() -> None:
    # phase-major layout: [phase0 nodes..., phase1 nodes...]
    assert phase_profile([1.0] * 6, phases=2, nodes=3) == [2.0, 2.0]
    reference_rate = [10.0] * 6
    candidate_rate = [8.0] * 3 + [9.0] * 3
    density = [2.0] * 6
    power = [4.0] * 6
    bands = phase_band_summary(
        candidate_rate, reference_rate, density, density, power, power,
        boundaries=[0, 1, 2], phases=2, nodes=3)
    assert math.isclose(
        sum(float(band["fraction_of_net_source_gap"]) for band in bands),
        1.0)
    assert bands[0]["candidate_to_reference_ionization_rate_ratio"] == .8


def test_phase_grid_rejects_wrong_order() -> None:
    try:
        validate_phase_grid([], "samples")
    except ValueError:
        pass
    else:
        raise AssertionError("empty phase grid was accepted")


if __name__ == "__main__":
    test_weighted_sum()
    test_additive_band_attribution()
    test_phase_attribution()
    test_phase_grid_rejects_wrong_order()
    print("spatial-ionization analyzer tests passed")
