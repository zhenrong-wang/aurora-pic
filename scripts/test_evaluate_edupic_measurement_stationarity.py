#!/usr/bin/env python3
"""Regression for predeclared eduPIC stationarity-rule evaluation."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
EVALUATOR = ROOT / "scripts/evaluate_edupic_measurement_stationarity.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="aurorapic_edupic_stationarity_",
                                     dir=ROOT / "tmp") as temporary:
        work = Path(temporary)
        rule = {
            "case_id": "case", "scope":
                "predeclared_internal_native_measurement_stationarity_screen",
            "analysis_contract": {
                "required_total_contiguous_blocks": 16,
                "required_total_measurement_cycles": 64,
            },
            "internal_density_stationarity_gates": {
                "minimum_ar1_effective_blocks_per_species": 8.0,
                "maximum_absolute_projected_fractional_drift_per_species": 0.01,
                "maximum_absolute_split_half_fractional_change_per_species": 0.01,
                "maximum_adjacent_profile_relative_l2_per_species": 0.025,
            },
        }
        series = {
            "ar1_effective_blocks": 9.0,
            "projected_fractional_drift_across_series": -0.002,
            "split_half_fractional_change": 0.003,
            "maximum_adjacent_profile_relative_l2": 0.02,
        }
        analysis = {
            "case_id": "case", "scope": "native_measurement_block_analysis",
            "analysis_eligible": True,
            "campaign": {"block_count": 16,
                         "completed_measurement_cycles": 64,
                         "target_reached": True},
            "density_profile": {"electron_series": dict(series),
                                "ion_series": dict(series)},
        }
        rule_path, analysis_path = work / "rule.json", work / "analysis.json"
        rule_path.write_text(json.dumps(rule), encoding="utf-8")
        analysis_path.write_text(json.dumps(analysis), encoding="utf-8")
        passed = subprocess.run([
            sys.executable, str(EVALUATOR), str(analysis_path), str(rule_path),
            "--require-pass",
        ], text=True, capture_output=True)
        require(passed.returncode == 0 and
                json.loads(passed.stdout)["passed"] is True,
                "evaluator rejected a passing series")
        analysis["density_profile"]["electron_series"][
            "ar1_effective_blocks"] = 7.0
        analysis_path.write_text(json.dumps(analysis), encoding="utf-8")
        failed = subprocess.run([
            sys.executable, str(EVALUATOR), str(analysis_path), str(rule_path),
            "--require-pass",
        ], text=True, capture_output=True)
        report = json.loads(failed.stdout)
        require(failed.returncode == 1 and report["passed"] is False and
                report["classification"] ==
                "internal_density_stationarity_screen_failed" and
                report["species_gates"]["electron"]
                ["minimum_ar1_effective_blocks"]["passed"] is False,
                "evaluator did not preserve a failed effective-count gate")
    print("eduPIC stationarity evaluation regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
