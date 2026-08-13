#!/usr/bin/env python3
"""Regression tests for the prospective density-bracket decision."""

from __future__ import annotations

from run_aurorapic_density_bracket import select_branch


def branch(electron: float, ion: float, slope: float,
           passes: bool = True) -> dict[str, object]:
    return {
        "metrics": {
            "electron_source_loss_relative_imbalance": electron,
            "ion_source_loss_relative_imbalance": ion,
            "normalized_total_population_slope_per_cycle": slope,
        },
        "safety_decision": {"passes": passes},
    }


def rule() -> dict[str, object]:
    return {
        "branches": {
            "control": {"added_pairs": 0},
            "plus25": {"added_pairs": 25},
            "plus50": {"added_pairs": 50},
        },
        "decision_rule": {
            "candidate_branches": ["plus25", "plus50"],
            "minimum_absolute_imbalance_improvement_vs_control": 0.05,
            "parsimony_tolerance": 0.02,
        },
    }


def main() -> None:
    result = select_branch({
        "control": branch(0.40, 0.42, 0.005),
        "plus25": branch(0.20, 0.22, 0.002),
        "plus50": branch(0.19, 0.21, 0.001),
    }, rule())
    assert result["selected_branch"] == "plus25"
    assert result["density_acceleration_supported"]

    rejected = select_branch({
        "control": branch(0.40, 0.42, 0.005),
        "plus25": branch(0.39, 0.40, 0.004),
        "plus50": branch(0.10, 0.12, 0.001, passes=False),
    }, rule())
    assert rejected["selected_branch"] is None
    assert not rejected["density_acceleration_supported"]

    best = select_branch({
        "control": branch(0.40, 0.42, 0.005),
        "plus25": branch(0.20, 0.25, 0.002),
        "plus50": branch(0.15, 0.18, 0.001),
    }, rule())
    assert best["selected_branch"] == "plus50"


if __name__ == "__main__":
    main()
