#!/usr/bin/env python3
"""Regression checks for the committed refresh-cadence result."""

import json
from pathlib import Path


def main() -> None:
    result = json.loads(Path(
        "benchmarks/ccp/edupic-ion-density-refresh-cadence-result-20260826.json"
    ).read_text(encoding="utf-8"))
    assert result["integrity_gate_passed"]
    assert result["initial_parity_gate_passed"]
    assert result["cadence_support_gate_passed"]
    assert result["early_trajectory_closure_through_horizon_100"]
    assert result["full_trace_closure_gate_passed"]
    assert not any(result["material_field_divergence_flags"])
    assert all(item["electron_population_difference"] == 0 and
               item["ion_population_difference"] == 0
               for item in result["comparisons"])
    print("held-density refresh-cadence analyzer regression passed")


if __name__ == "__main__":
    main()
