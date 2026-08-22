#!/usr/bin/env python3
"""Focused acceptance tests for surface-flux timestep refinement."""

from analyze_aurorapic_surface_flux_timestep import evaluate


def result(critical: float, exceptional: float, tail: float,
           critical_closure: float, exceptional_closure: float) -> dict:
    return {
        "critical_phase_0p125_to_0p5": {
            "direct_outward_energy_flux_divergence_W_m-2": critical,
            "relative_closure_error": critical_closure,
        },
        "exceptional_phase_0p375_to_0p5": {
            "direct_outward_energy_flux_divergence_W_m-2": exceptional,
            "approximate_above_15p8_eV_outward_energy_flux_divergence_W_m-2":
                tail,
            "relative_closure_error": exceptional_closure,
        },
    }


def main() -> None:
    limits = {
        "maximum_critical_phase_direct_flux_relative_change": 0.15,
        "maximum_exceptional_octant_direct_flux_relative_change": 0.10,
        "maximum_exceptional_octant_tail_flux_relative_change": 0.20,
        "maximum_critical_phase_closure_error_each_branch": 0.10,
        "maximum_exceptional_octant_closure_error_each_branch": 0.06,
    }
    metrics, gates = evaluate(
        result(27, 116, 20, .04, .01),
        result(26, 117, 19, .05, .02), limits)
    assert all(gates.values()) and metrics[
        "critical_phase_direct_flux_relative_change"] < .15
    _, failed = evaluate(
        result(27, 116, 20, .04, .01),
        result(20, 117, 19, .05, .02), limits)
    assert failed["critical_phase_direct_flux"] is False


if __name__ == "__main__":
    main()
    print("surface-flux timestep analyzer tests passed")
