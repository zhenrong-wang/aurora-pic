#!/usr/bin/env python3
"""Focused ensemble arithmetic tests for corrected phase-EEDF analysis."""

from analyze_edupic_corrected_phase_eedf import compare_ensembles


def item(scale: float) -> dict[str, object]:
    names = (
        "histogram_mean_energy_eV", "fraction_11p5_to_15p8_eV",
        "fraction_15p8_to_30_eV", "fraction_at_or_above_30_eV",
        "fraction_at_or_above_15p8_eV",
        "eedf_folded_ionization_frequency_s-1")
    return {**{name: scale for name in names},
            "probability_by_energy_bin": [scale / 320.0] * 320}


def main() -> None:
    result = compare_ensembles([item(1.0), item(1.0)],
                               [item(1.0), item(1.0)])
    assert all(value == 1.0 for value in
               result["aurorapic_to_native_edupic_ratio"].values())
    assert result["total_variation_distance"] == 0.0
    shifted = compare_ensembles([item(0.9)], [item(1.0)])
    assert abs(shifted["aurorapic_to_native_edupic_ratio"][
        "eedf_folded_ionization_frequency_s-1"] - 0.9) < 1.0e-15
    print("corrected phase-EEDF analyzer tests passed")


if __name__ == "__main__":
    main()
