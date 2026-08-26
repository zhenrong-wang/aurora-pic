#!/usr/bin/env python3
"""Regression checks for the committed held-density control reduction."""

import json
from pathlib import Path


def main() -> None:
    result = json.loads(Path(
        "benchmarks/ccp/edupic-ion-density-refresh-control-result-20260826.json"
    ).read_text(encoding="utf-8"))
    assert result["integrity_gate_passed"]
    assert result["initial_parity_gate_passed"]
    assert result["strong_mechanism_support_gate_passed"]
    assert not result["early_divergence_explained_through_horizon_20"]
    early = {item["aurorapic_horizon"]: item
             for item in result["comparisons"]}
    assert all(early[h]["critical_log_error_fraction_of_baseline"] <= .2
               for h in (1, 2, 5))
    assert early[20]["critical_field_energy_ratio"] > 1.02
    print("held-density control analyzer regression passed")


if __name__ == "__main__":
    main()
