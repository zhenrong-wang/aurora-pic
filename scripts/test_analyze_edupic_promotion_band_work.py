#!/usr/bin/env python3
"""Focused tests for promotion-band work comparison calculations."""

from analyze_edupic_promotion_band_work import (
    classify_ratio, means, relative_ranges,
)


def main() -> None:
    members = [
        {
            "band_supply_fraction": 0.9,
            "band_promotion_probability": 1.8,
            "mean_signed_work_eV": 2.7,
            "mean_positive_work_eV": 3.6,
            "mean_negative_work_eV": 0.9,
        },
        {
            "band_supply_fraction": 1.1,
            "band_promotion_probability": 2.2,
            "mean_signed_work_eV": 3.3,
            "mean_positive_work_eV": 4.4,
            "mean_negative_work_eV": 1.1,
        },
    ]
    assert means(members)["mean_positive_work_eV"] == 4.0
    assert abs(relative_ranges(members)["band_supply_fraction"] - 0.2) < 1e-15
    assert classify_ratio(0.94, (0.95, 1.05)) == "aurorapic_lower"
    assert classify_ratio(1.00, (0.95, 1.05)) == "parity"
    assert classify_ratio(1.06, (0.95, 1.05)) == "aurorapic_higher"


if __name__ == "__main__":
    main()
    print("promotion-band work analyzer tests passed")
