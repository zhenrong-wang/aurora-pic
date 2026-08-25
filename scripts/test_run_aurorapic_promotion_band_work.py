#!/usr/bin/env python3
"""Focused deck tests for the promotion-band work runner."""

import json
from pathlib import Path

from run_aurorapic_promotion_band_work import make_deck


def main() -> int:
    rule = json.loads(Path(
        "benchmarks/ccp/edupic-promotion-band-work-rule-20260825.json"
    ).read_text(encoding="utf-8"))
    prior = Path(
        "tmp/edupic-matched-half-step-replication-51949-20260825/input.cfg"
    ).read_text(encoding="utf-8")
    deck = make_deck(prior, Path("new-output"), Path("checkpoint.apc"), rule)
    for expected in (
        "steps = 52000", "spatial_average_start_step = 36001",
        "spatial_average_end_step = 52000", "checkpoint_interval = 52000",
        "restart_path = checkpoint.apc",
        "phase_eedf_promotion_band_min = 11.5",
        "collision_velocity_sampling = leapfrog_half_step",
    ):
        assert expected in deck, expected
    assert deck.count("phase_eedf_promotion_band_min") == 1
    print("promotion-band work runner tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
