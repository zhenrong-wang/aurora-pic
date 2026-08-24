#!/usr/bin/env python3
"""Focused arithmetic tests for phase-EEDF comparison helpers."""

from analyze_edupic_phase_eedf_crosscode import compare, relative_range


def sample(mean: float, tail: float, folded: float, probabilities: list[float]) -> dict:
    return {
        "histogram_mean_energy_eV": mean,
        "fraction_11p5_to_15p8_eV": .03,
        "fraction_15p8_to_30_eV": tail,
        "fraction_at_or_above_30_eV": .01,
        "fraction_at_or_above_15p8_eV": tail + .01,
        "eedf_folded_ionization_frequency_s-1": folded,
        "probability_by_energy_bin": probabilities,
    }


def main() -> None:
    candidate = sample(4, .04, 80, [.5, .5] + [0] * 318)
    natives = [sample(5, .05, 100, [.4, .6] + [0] * 318),
               sample(5, .05, 100, [.4, .6] + [0] * 318)]
    result = compare(candidate, natives)
    assert result["aurorapic_to_native_edupic_ratio"][
        "eedf_folded_ionization_frequency_s-1"] == .8
    assert abs(result["total_variation_distance"] - .1) < 1e-15
    assert relative_range([9, 10, 11]) == .2


if __name__ == "__main__":
    main()
    print("native phase-EEDF cross-code analyzer tests passed")
