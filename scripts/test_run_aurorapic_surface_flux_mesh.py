#!/usr/bin/env python3
"""Focused deck-generation test for surface-flux mesh refinement."""

import json
from pathlib import Path

from run_aurorapic_surface_flux_mesh import build_deck


def main() -> None:
    rule = json.loads(Path(
        "benchmarks/ccp/edupic-argon-surface-flux-mesh-rule-20260822.json"
    ).read_text(encoding="utf-8"))
    base = Path(
        "tmp/edupic-argon-cross-sections-20260810/one-cycle-diagnostic.cfg"
    ).read_text(encoding="utf-8")
    deck = build_deck(
        base, Path("output"), Path("checkpoint.apc"), rule, "refined_grid")
    for expected in (
        "nx = 799", "steps = 32000",
        "spatial_average_start_step = 24001",
        "spatial_average_end_step = 32000",
        "spatial_average_phase_bins = 200",
        "phase_surface_flux_positions = 0.005,0.015",
        "phase_surface_flux_reset_on_restart = true",
        "wall_impact_reset_on_restart = true",
        "restart_path = checkpoint.apc",
    ):
        assert expected in deck, expected

    particle_rule = json.loads(Path(
        "benchmarks/ccp/"
        "edupic-argon-surface-flux-particle-rule-20260823.json"
    ).read_text(encoding="utf-8"))
    particle_deck = build_deck(
        base, Path("particle-output"), Path("particle-checkpoint.apc"),
        particle_rule, "double_particles")
    assert "nx = 799" in particle_deck
    assert particle_deck.count("weight = 350000000.0") == 2
    assert "particles = 238898" in particle_deck
    assert "particles = 249890" in particle_deck


if __name__ == "__main__":
    main()
    print("surface-flux mesh runner tests passed")
