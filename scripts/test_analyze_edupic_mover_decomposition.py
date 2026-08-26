#!/usr/bin/env python3
"""Focused tests for the cross-code mover-decomposition analyzer."""

from analyze_edupic_mover_decomposition import classify, means, relative_ranges


def main() -> None:
    members = [
        {"mean_origin_energy_eV": 0.9,
         "origin_longitudinal_energy_fraction": 1.8,
         "mean_positive_linear_work_eV": 2.7,
         "mean_quadratic_work_eV": 3.6},
        {"mean_origin_energy_eV": 1.1,
         "origin_longitudinal_energy_fraction": 2.2,
         "mean_positive_linear_work_eV": 3.3,
         "mean_quadratic_work_eV": 4.4},
    ]
    assert means(members)["mean_quadratic_work_eV"] == 4.0
    assert abs(relative_ranges(members)["mean_origin_energy_eV"] - 0.2) < 1e-15
    assert classify(0.94, (0.95, 1.05)) == "aurorapic_lower"
    assert classify(1.0, (0.95, 1.05)) == "parity"
    assert classify(1.06, (0.95, 1.05)) == "aurorapic_higher"
    print("cross-code mover-decomposition analyzer tests passed")


if __name__ == "__main__":
    main()
