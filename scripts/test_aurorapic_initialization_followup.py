#!/usr/bin/env python3
"""Lightweight contract tests for the warm-state transient follow-up."""

from __future__ import annotations

from run_aurorapic_initialization_followup import configured_base, evaluate


def stage(ionization: int, losses: int, growth: float,
          field: float, passes: bool = True) -> dict[str, object]:
    return {
        "passes": passes,
        "population": {
            "ionization_pairs": ionization,
            "electron_wall_losses": losses,
            "total_growth_factor": growth,
        },
        "maximum_sampled_absolute_electric_field_V_m": field,
    }


def main() -> int:
    base = (
        "[species.electrons]\nweight = 1\nparticles = 2\n\n"
        "[species.ions]\nweight = 1\nparticles = 2\n")
    rule = {
        "initial_state": {
            "electrons": 3, "ions": 4, "macro_weight": 2.5},
        "decision_rule": {"prospective_thresholds": {
            "maximum_first_cycle_ionization_pairs": 20,
            "minimum_first_cycle_electron_wall_losses": 5,
            "maximum_first_cycle_total_growth_factor": 1.02,
            "maximum_sampled_absolute_field_V_m": 100.0,
        }},
    }
    configured = configured_base(base, rule)
    if not all(token in configured for token in (
            "weight = 2.5", "particles = 3", "particles = 4")):
        raise RuntimeError("follow-up deck configuration changed")
    passing = evaluate([
        stage(20, 5, 1.02, 100.0), stage(1, 1, 1.0, 90.0)], rule)
    failing = evaluate([
        stage(21, 4, 1.021, 101.0, False),
        stage(1, 1, 1.0, 90.0)], rule)
    if not passing["passes"] or failing["passes"]:
        raise RuntimeError("prospective threshold evaluation changed")
    if any(failing["gates"].values()):
        raise RuntimeError("synthetic failure did not exercise every gate")
    print("initialization follow-up regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
