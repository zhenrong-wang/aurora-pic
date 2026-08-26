#!/usr/bin/env python3
"""Focused interval tests for the four-period balance analyzer."""

from analyze_edupic_four_period_electron_balance import equivalence_interval


def main() -> None:
    value = equivalence_interval(100.0, 10.0, 110.0, 20.0)
    assert value["aurorapic_minus_edupic_mean_particles"] == 10.0
    assert abs(value["standard_error_particles"] - 10.0) < 1.0e-15
    assert value["lower_particles"] < 0.0
    assert value["upper_particles"] > 30.0
    zero = equivalence_interval(100.0, 0.0, 100.0, 0.0)
    assert zero["lower_particles"] == 0.0
    assert zero["upper_particles"] == 0.0
    print("four-period electron-balance analyzer tests passed")


if __name__ == "__main__":
    main()
