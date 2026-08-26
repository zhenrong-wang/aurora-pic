#!/usr/bin/env python3
"""Focused tests for the guarded AuroraPIC mover-decomposition runner."""

import json
from pathlib import Path

from run_aurorapic_mover_decomposition import summarize
from run_aurorapic_promotion_band_work import make_deck


def main() -> None:
    row = {
        "field_push_macro_observations": "100",
        "field_push_promotion_band_observations": "10",
        "field_push_promotion_band_promotions": "2",
        "field_push_promotion_band_signed_macro_work_sum_eV": "5",
        "field_push_promotion_band_positive_macro_work_sum_eV": "8",
        "field_push_promotion_band_negative_macro_work_sum_eV": "3",
        "field_push_promotion_band_origin_macro_energy_sum_eV": "130",
        "field_push_promotion_band_origin_longitudinal_macro_energy_sum_eV": "52",
        "field_push_promotion_band_linear_macro_work_sum_eV": "4",
        "field_push_promotion_band_positive_linear_macro_work_sum_eV": "7",
        "field_push_promotion_band_negative_linear_macro_work_sum_eV": "3",
        "field_push_promotion_band_quadratic_macro_work_sum_eV": "1",
    }
    result = summarize([row])
    assert result["band_supply_fraction"] == 0.1
    assert result["band_promotion_probability"] == 0.2
    assert result["mean_origin_energy_eV"] == 13.0
    assert result["origin_longitudinal_energy_fraction"] == 0.4
    assert result["mean_positive_linear_work_eV"] == 0.7
    assert result["mean_quadratic_work_eV"] == 0.1
    assert result["linear_work_closure_residual_eV"] == 0.0
    assert result["total_work_decomposition_closure_residual_eV"] == 0.0
    rule = json.loads(Path(
        "benchmarks/ccp/edupic-mover-decomposition-rule-20260825.json"
    ).read_text(encoding="utf-8"))
    prior = Path(
        "tmp/edupic-matched-half-step-replication-51949-20260825/input.cfg"
    ).read_text(encoding="utf-8")
    deck = make_deck(prior, Path("new-output"), Path("checkpoint.apc"), rule)
    assert "phase_eedf_promotion_band_min = 11.5" in deck
    assert "steps = 52000" in deck
    print("AuroraPIC mover-decomposition runner tests passed")


if __name__ == "__main__":
    main()
