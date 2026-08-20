#!/usr/bin/env python3
"""Focused tests for balance-campaign aggregation."""

from analyze_aurorapic_balance_campaign import summarize


def report(index: int, created: int, electron_loss: int,
           ion_loss: int, electrons: int) -> dict[str, object]:
    return {
        "block_index": index,
        "window": {"cycles": 10},
        "particle_ledger": {
            "ionization_macro_events": created,
            "electron_wall_loss_macro_events": electron_loss,
            "ion_wall_loss_macro_events": ion_loss,
        },
        "cycle_endpoint_populations": {"electrons": [electrons] * 11},
        "metrics": {},
        "stationarity_gates": {},
        "stationarity_block_passed": False,
    }


def test_summary() -> None:
    result = summarize([
        report(1, 1000, 1050, 1040, 10000),
        report(2, 900, 1050, 1040, 9000),
    ])
    assert result["combined"]["cycles"] == 20
    assert result["combined"]["ionization_macro_events"] == 1900
    assert result["combined"]["electron_wall_loss_macro_events"] == 2100
    assert abs(result["first_to_last_change"][
        "ionization_macro_events_per_cycle_relative_change"] + .1) < 1e-15
    assert result["first_to_last_change"][
        "electron_wall_loss_macro_events_per_cycle_relative_change"] == 0.0
    # The lower raw source follows the lower population in this synthetic case.
    assert result["first_to_last_change"][
        "ionization_events_per_live_electron_per_cycle_relative_change"] == 0.0


if __name__ == "__main__":
    test_summary()
    print("balance-campaign analyzer tests passed")
