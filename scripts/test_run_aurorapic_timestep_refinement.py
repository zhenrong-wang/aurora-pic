#!/usr/bin/env python3
"""Deck-generation regressions for the timestep refinement runner."""

import json
from pathlib import Path
import tempfile

from run_aurorapic_timestep_refinement import initial_deck, measurement_deck


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


if __name__ == "__main__":
    main()
