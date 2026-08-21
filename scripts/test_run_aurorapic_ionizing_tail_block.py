#!/usr/bin/env python3
"""Focused deck-generation tests for the ionizing-tail continuation."""

import json
from pathlib import Path

from run_aurorapic_ionizing_tail_block import build_deck


def main() -> None:
    rule = json.loads(Path(
        "benchmarks/ccp/edupic-argon-ionizing-tail-rule-20260821.json"
    ).read_text(encoding="utf-8"))
    base = Path(
        "tmp/edupic-argon-cross-sections-20260810/"
        "matched-heating-long-window-20260819/measurement/input.cfg"
    ).read_text(encoding="utf-8")
    deck = build_deck(base, Path("output"), Path("checkpoint.apc"), rule)
    required = (
        "steps = 72000", "spatial_average_interval = 2",
        "spatial_average_start_step = 56001",
        "spatial_average_phase_bins = 200",
        "phase_eedf_energy_bins = 320", "phase_eedf_energy_max = 80.0",
        "phase_eedf_regions = x000_010:0.0:0.0025",
        "runtime_backend = serial", "runtime_threads = 1",
        "restart_path = checkpoint.apc",
    )
    for text in required:
        assert text in deck, text
    assert deck.count("phase_eedf_regions =") == 1


if __name__ == "__main__":
    main()
    print("ionizing-tail runner tests passed")
