#!/usr/bin/env python3
"""Regression tests for long-window matched-heating gates."""

from analyze_aurorapic_long_window_heating import (
    ENSEMBLE_DENSITY_M3, ENSEMBLE_POWER_PER_ELECTRON_RATIO,
    ENSEMBLE_POWER_W_M3, gate_summary,
)


def main() -> None:
    passing = gate_summary(
        ENSEMBLE_POWER_W_M3, ENSEMBLE_POWER_W_M3,
        ENSEMBLE_DENSITY_M3, ENSEMBLE_POWER_PER_ELECTRON_RATIO)
    assert passing["all_gates_passed"] is True
    failing = gate_summary(
        1.1 * ENSEMBLE_POWER_W_M3, ENSEMBLE_POWER_W_M3,
        ENSEMBLE_DENSITY_M3, ENSEMBLE_POWER_PER_ELECTRON_RATIO)
    assert failing["all_gates_passed"] is False
    assert failing["gates"]["phase_binned_to_exact_power"] is False


if __name__ == "__main__":
    main()
