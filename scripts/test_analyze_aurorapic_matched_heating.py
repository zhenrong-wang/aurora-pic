#!/usr/bin/env python3
"""Regression tests for matched electron-heating analysis."""

import math

from analyze_aurorapic_matched_heating import phase_average, summarize


def main() -> None:
    result = summarize(2.0, 10.1, 10.0, 4.0, 20.0)
    assert result["internal_power_gate_passed"] is True
    assert math.isclose(
        result["candidate_to_reference_phase_binned_power_per_particle_ratio"],
        1.01)
    assert math.isclose(phase_average([0.0, 1.0], [2.0, 2.0, 4.0, 4.0], 2),
                        3.0)
    try:
        phase_average([0.0, 1.0], [1.0], 2)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid phase-space shape was accepted")


if __name__ == "__main__":
    main()
