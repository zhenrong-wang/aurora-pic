#!/usr/bin/env python3
"""Focused tests for common-state divergence reductions."""

import csv
import json
import tempfile
from pathlib import Path

from analyze_edupic_common_state_divergence import analyze, relative_rms


def main() -> None:
    assert relative_rms([2.0, 2.0], [2.0, 2.0], 1.0) == 0.0
    assert abs(relative_rms([2.2, 2.2], [2.0, 2.0], 1.0) - .1) < 1e-12
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        rule = {
            "case_id": "test", "claim_boundary": "test only",
            "locked_inputs": {"particle_state_sha256": "state",
                              "aurorapic_binary_sha256": "binary"},
            "physics_contract": {"nodes": 2, "length_m": 1.0},
            "sampling_contract": {
                "edupic_pre_push_steps": [1, 2, 3, 4],
                "matching_aurorapic_post_step_horizons": [0, 1, 2, 3]},
        }
        rule_path = root / "rule.json"
        rule_path.write_text(json.dumps(rule), encoding="utf-8")
        native_trace = root / "native.csv"
        with native_trace.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("pre_push_step", "node", "x_m",
                             "charge_density_C_m3", "electric_field_V_m"))
            for step in (1, 2, 3, 4):
                writer.writerows(((step, 0, 0, 1, 1),
                                  (step, 1, 1, 1, 1)))
        population = root / "population.csv"
        population.write_text(
            "pre_push_step,electrons,ions\n1,10,11\n2,10,11\n3,10,11\n4,10,11\n",
            encoding="utf-8")
        checkpoint = root / "checkpoint.bin"
        checkpoint.write_bytes(b"same")
        aurora_root = root / "aurora"
        members = []
        for horizon, field_value in ((0, 1.0), (1, .9), (2, .9), (3, .9)):
            output = aurora_root / f"horizon-{horizon:04d}" / "output"
            output.mkdir(parents=True)
            field = output / f"fields_{horizon}.csv"
            field.write_text(
                "x,rho,E\n0,1," + str(field_value) +
                "\n1,1," + str(field_value) + "\n", encoding="utf-8")
            members.append({"horizon": horizon,
                            "field_sha256": __import__("hashlib").sha256(
                                field.read_bytes()).hexdigest(),
                            "electron_population": 10,
                            "ion_population": 11})
        report = {"particle_state_sha256": "state",
                  "binary_sha256": "binary",
                  "rule_sha256": __import__("hashlib").sha256(
                      rule_path.read_bytes()).hexdigest(),
                  "all_resource_gates_passed": True, "members": members}
        report_path = root / "report.json"
        report_path.write_text(json.dumps(report), encoding="utf-8")
        result = analyze(rule_path, native_trace, population, checkpoint,
                         checkpoint, native_trace, report_path, aurora_root)
        assert result["integrity_gate_passed"]
        assert result["initial_parity_gate_passed"]
        assert result["earliest_sustained_material_divergence_horizon"] == 1
        assert result["formal_outcome"] == "one_step_mover_or_boundary_mismatch"
        assert result["post_hoc_mechanism_candidate"][
            "requires_prospective_control"] is True
    print("common-state divergence analyzer tests passed")


if __name__ == "__main__":
    main()
