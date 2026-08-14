#!/usr/bin/env python3
"""Regression tests for electron-heating ledger localization."""

import math

from analyze_aurorapic_heating_ledger import summarize


def main() -> None:
    result = summarize(2.0, 8.0, 10.0, 4.0, 16.0)
    assert result["candidate_phase_binned_to_exact_power_ratio"] == 0.8
    assert result["candidate_to_reference_phase_binned_power_per_particle_ratio"] == 1.0
    assert math.isclose(
        result["candidate_exact_to_reference_binned_power_per_particle_ratio"],
        1.25)
    assert result["apparent_deficit_reduction_fraction"] is None

    deficit = summarize(2.0, 8.0, 10.0, 4.0, 20.0)
    assert math.isclose(
        deficit["apparent_deficit_reduction_fraction"], 1.0)
    try:
        summarize(0.0, 1.0, 1.0, 1.0, 1.0)
    except ValueError:
        pass
    else:
        raise AssertionError("zero candidate density was accepted")


if __name__ == "__main__":
    main()
