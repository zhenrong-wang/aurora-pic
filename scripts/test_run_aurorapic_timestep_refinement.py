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
        assert "wall_impact_reset_on_restart = true" in second


if __name__ == "__main__":
    main()
