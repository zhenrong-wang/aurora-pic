#!/usr/bin/env python3
"""Focused deck test for the corrected phase-EEDF runner."""

import json
from pathlib import Path

from run_aurorapic_corrected_phase_eedf import deck


def main() -> None:
    rule = json.loads(Path(
        "benchmarks/ccp/edupic-corrected-cadence-phase-eedf-rule-20260826.json"
    ).read_text())
    value = deck(rule, Path("state.aps"), 13507, Path("e.gas"),
                 Path("i.gas"), Path("output"))
    assert "steps = 16000" in value
    assert "spatial_average_start_step = 1" in value
    assert "spatial_average_end_step = 16000" in value
    assert "subcycle_charge_deposition = pre_push_held" in value
    assert "phase_eedf_energy_bins = 320" in value
    assert "initial_state_signature = 16659304071935497571" in value
    assert value.count("model = null_collision") == 2
    print("corrected phase-EEDF runner tests passed")


if __name__ == "__main__":
    main()
