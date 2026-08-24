#!/usr/bin/env python3
"""Focused tests for the native eduPIC/AuroraPIC comparison arithmetic."""

from analyze_edupic_surface_flux_crosscode import evaluate


def member(critical: float, exceptional: float, tail: float) -> dict:
    return {
        "macro_crossings_by_surface": [20000, 21000],
        "maximum_overflow_fraction": 0.0,
        "critical_phase_0p125_to_0p5": {
            "direct_outward_energy_flux_divergence_W_m-2": critical,
        },
        "exceptional_phase_0p375_to_0p5": {
            "direct_outward_energy_flux_divergence_W_m-2": exceptional,
            "approximate_above_15p8_eV_outward_energy_flux_divergence_W_m-2":
                tail,
        },
    }


def main() -> None:
    rule = {
        "comparison": {"aurorapic_to_edupic_electron_density_ratio": 1.2},
        "prospective_acceptance": {
            "maximum_critical_direct_relative_range": .2,
            "maximum_exceptional_direct_relative_range": .15,
            "maximum_exceptional_tail_relative_range": .3,
            "minimum_macro_crossings_each_surface": 10000,
            "maximum_histogram_overflow_fraction": .001,
        },
    }
    aurora = {"metrics": {
        "mean_critical_phase_direct_flux_W_m-2": 24.0,
        "mean_exceptional_octant_direct_flux_W_m-2": 120.0,
        "mean_exceptional_octant_tail_flux_W_m-2": 18.0,
    }}
    metrics, gates = evaluate(rule, aurora, [
        member(20, 100, 15), member(21, 102, 14), member(19, 98, 16)])
    assert all(gates.values())
    assert metrics["edupic_ensemble_mean"]["critical_direct_W_m-2"] == 20
    assert metrics["density_normalized_aurorapic_to_edupic_ratio"][
        "exceptional_direct_W_m-2"] == 1.0
    _, failed = evaluate(rule, aurora, [
        member(20, 100, 15), member(30, 102, 14), member(19, 98, 16)])
    assert failed["critical_repeatability"] is False


if __name__ == "__main__":
    main()
    print("native eduPIC surface-flux cross-code analyzer tests passed")
