#!/usr/bin/env python3
"""Bounded regression tests for phase-space comparison math."""

import math

from compare_edupic_phase_space import (
    flatten_phase_major, lower_bin_value, metrics, periodic_overlap_average,
    phase_effective_frequency, phase_space_metrics, resample_matrix,
    spatial_phase_average,
)


def main() -> int:
    source = [float(index) for index in range(200)]
    reduced = periodic_overlap_average(source, 16)
    assert len(reduced) == 16
    assert abs(sum(reduced) / 16 - sum(source) / 200) < 1e-12
    assert abs(reduced[0] - 5.76) < 1e-12
    assert abs(reduced[1] - 18.24) < 1e-12

    exact = periodic_overlap_average([1.0, 1.0, 3.0, 3.0], 2)
    assert exact == [1.0, 3.0]
    matrix = resample_matrix([[1.0, 1.0, 3.0, 3.0],
                              [2.0, 2.0, 4.0, 4.0]], 2)
    assert flatten_phase_major(matrix) == [1.0, 2.0, 3.0, 4.0]

    same = metrics([1.0, -2.0, 3.0], [1.0, -2.0, 3.0])
    assert same["relative_l2"] == 0.0
    assert same["mean_bias_over_reference_rms"] == 0.0
    assert abs(same["pearson_correlation"] - 1.0) < 1e-15
    assert same["candidate_to_reference_peak_absolute_ratio"] == 1.0
    scaled = metrics([2.0, 4.0], [1.0, 2.0])
    assert abs(scaled["relative_l2"] - 1.0) < 1e-15
    assert math.isclose(scaled["pearson_correlation"], 1.0)
    phase = phase_space_metrics([1.0, 2.0, 2.0, 4.0],
                                [1.0, 2.0, 1.0, 2.0], phases=2, nodes=2)
    assert phase["phase_profile_relative_l2"] == [0.0, 1.0]
    assert phase["maximum_phase_profile_relative_l2_bin"] == 1
    assert math.isclose(phase["cycle_average_spatial_profile_relative_l2"], 0.5)
    assert spatial_phase_average(
        [1.0, 2.0, 3.0, 1.0, 2.0, 3.0], phases=2, nodes=3) == 2.0
    assert phase_effective_frequency(
        [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
        [1.0, 1.0, 1.0, 2.0, 2.0, 2.0], phases=2, nodes=3) == [2.0, 2.0]
    assert lower_bin_value([0.0, 1.0, 2.0], [0.0, 4.0, 9.0], 1.5) == 4.0
    print("eduPIC phase-space comparison regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
