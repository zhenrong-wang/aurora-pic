#!/usr/bin/env python3
"""Deck-generation regressions for the timestep refinement runner."""

import json
from pathlib import Path
import tempfile

from run_aurorapic_timestep_refinement import (
    branch_state, initial_deck, measurement_deck,
)


def main() -> None:
    rule = json.loads(Path(
        "benchmarks/ccp/edupic-argon-timestep-refinement-rule-20260813.json"
    ).read_text(encoding="utf-8"))
    base = Path(
        "tmp/edupic-argon-cross-sections-20260810/one-cycle-diagnostic.cfg"
    ).read_text(encoding="utf-8")
    with tempfile.TemporaryDirectory(prefix="aurorapic_dt_refinement_") as temporary:
        root = Path(temporary)
        first, equilibration, measurement = initial_deck(
            base, rule, "half_dt", root / "initial", root / "state.aps")
        assert equilibration == 16000 and measurement == 32000
        assert "dt = 9.218289085545723e-12" in first
        assert "particles = 119449" in first
        assert "initial_state_signature = 7816361585110270686" in first
        second = measurement_deck(
            base, rule, "half_dt", root / "measurement", root / "checkpoint.apc")
        assert "steps = 48000" in second
        assert "spatial_average_start_step = 16001" in second
        assert "spatial_average_phase_bins = 16" in second
        assert "spatial_average_sampling_order = post_collision" in second
        assert "wall_impact_reset_on_restart = true" in second
        mesh_rule = json.loads(Path(
            "benchmarks/ccp/edupic-argon-mesh-refinement-rule-20260814.json"
        ).read_text(encoding="utf-8"))
        mesh, _, _ = initial_deck(
            base, mesh_rule, "refined_grid", root / "mesh", root / "state.aps")
        assert "nx = 799" in mesh
        particle_rule = json.loads(Path(
            "benchmarks/ccp/edupic-argon-particle-refinement-rule-20260814.json"
        ).read_text(encoding="utf-8"))
        particle, _, _ = initial_deck(
            base, particle_rule, "double_particles",
            root / "particle", root / "split.aps")
        assert "particles = 238898" in particle
        assert particle.count("weight = 350000000.0") == 2
        region_rule = json.loads(Path(
            "benchmarks/ccp/edupic-argon-region-matched-eedf-rule-20260814.json"
        ).read_text(encoding="utf-8"))
        region = measurement_deck(
            base, region_rule, "region_matched", root / "region",
            root / "checkpoint.apc")
        assert ("phase_eedf_regions = full_gap:0:0.025,"
                "edupic_center_10pct:0.01125:0.01375") in region
        matched_rule = json.loads(json.dumps(region_rule))
        matched_rule["fresh_measurement_contract"]["spatial_phase_bins"] = 200
        matched_rule["fresh_measurement_contract"][
            "spatial_sampling_order"] = "pre_collision"
        matched = measurement_deck(
            base, matched_rule, "region_matched", root / "matched",
            root / "checkpoint.apc")
        assert "spatial_average_phase_bins = 200" in matched
        assert "spatial_average_sampling_order = pre_collision" in matched
        heating_rule = json.loads(Path(
            "benchmarks/ccp/edupic-argon-matched-heating-rule-20260814.json"
        ).read_text(encoding="utf-8"))
        heating = measurement_deck(
            base, heating_rule, "matched_heating", root / "heating",
            root / "checkpoint.apc")
        assert "spatial_average_phase_bins = 200" in heating
        assert "spatial_average_sampling_order = pre_collision" in heating
        seeded_rule = json.loads(json.dumps(heating_rule))
        seeded_rule["branches"]["matched_heating"]["seed"] = 24601
        seeded_rule["branches"]["matched_heating"]["output_interval"] = 400
        seeded, _, _ = initial_deck(
            base, seeded_rule, "matched_heating", root / "seeded",
            root / "state.aps")
        assert "seed = 24601" in seeded
        assert "output_interval = 400" in seeded
        multi_state = json.loads(json.dumps(seeded_rule))
        multi_state["particle_states"] = {
            "randomized": {
                "particle_state_sha256": "a" * 64,
                "particle_state_signature": 123456789,
            }}
        multi_state["branches"]["matched_heating"][
            "particle_state"] = "randomized"
        selected = branch_state(multi_state, "matched_heating")
        assert selected["particle_state_sha256"] == "a" * 64
        assert selected["particle_state_signature"] == 123456789
        assert selected["electrons"] == 119449
        randomized, _, _ = initial_deck(
            base, multi_state, "matched_heating", root / "randomized",
            root / "randomized.aps")
        assert "initial_state_signature = 123456789" in randomized
        long_rule = json.loads(Path(
            "benchmarks/ccp/edupic-argon-long-window-heating-rule-20260819.json"
        ).read_text(encoding="utf-8"))
        long_window, equilibration, measurement = initial_deck(
            base, long_rule, "long_window", root / "long", root / "state.aps")
        assert equilibration == 8000 and measurement == 48000
        assert "output_interval = 400" in long_window
        assert "seed = 46829" in long_window


if __name__ == "__main__":
    main()
