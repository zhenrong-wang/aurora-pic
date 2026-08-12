#!/usr/bin/env python3
"""Synthetic regression for the fail-closed Turner credibility audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile

from audit_turner_credibility import AuditError, audit


def write(path: Path, value: dict[str, object]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        evidence_dir = root / "evidence"
        inputs = {
            "published_duration_ensemble": {
                "case_id": "turner-helium-ccp-2013-case-1",
                "verified_members": 3,
                "published_individual_run_results": {
                    "accepted_95_percent_count": 0,
                    "accepted_99_percent_count": 1,
                    "failed_99_percent_count": 2,
                },
                "density_amplitude": {
                    "all_member_biases_same_sign": True,
                    "integrated_bias_percent_mean": 2.2,
                },
            },
            "numerical_sensitivity": {
                "case_id": "turner-helium-ccp-2013-case-1",
                "interpretation": {
                    "conclusion":
                        "no_tested_refinement_materially_reduced_density_bias"},
                "global_observable_ranges": {
                    "exact_source_loss_balance_all_variants": True},
            },
            "spatial_collision": {
                "numerical_closure": {
                    "passes": True,
                    "particle_source_loss_balance_exact": True,
                },
                "power_balance_W_m-2": {"electric_work_total": 1.0},
                "turner_context": {
                    "midplane_ion_density_relative_difference": 0.03,
                    "midplane_electron_temperature_relative_difference": -0.06,
                    "electron_electrical_power_relative_difference": -0.002,
                    "ion_electrical_power_relative_difference": 0.005,
                    "mean_ion_current_relative_difference": 0.004,
                },
            },
            "phase_structure": {
                "scope": "post_benchmark_phase_resolved_diagnostic_window"},
            "phase_eedf": {
                "scope": "post_benchmark_32_cycle_phase_resolved_regional_eedf",
                "quality": {"passes": True}},
            "energy_ledger": {"acceptance": {"passes": True}},
        }
        locks = {}
        for name, value in inputs.items():
            relative = f"evidence/{name}.json"
            digest = write(root / relative, value)
            locks[name] = {"path": relative, "sha256": digest}
        contract = root / "benchmarks" / "ccp" / "contract.json"
        write(contract, {
            "schema_version": 1,
            "case_id": "turner-helium-ccp-2013-case-1",
            "scope": "turner_credibility_evidence_contract",
            "evidence": locks,
        })
        report = audit(contract)
        if not (
            report["classification"] ==
                "partial_code_to_code_verification_systematic_density_discrepancy"
            and report["numerical_integrity"]["passed"]
            and not report["published_duration_density"][
                "published_density_verification_passed"]
        ):
            raise RuntimeError("Turner credibility classification changed")
        (evidence_dir / "energy_ledger.json").write_text(
            "{}\n", encoding="utf-8")
        try:
            audit(contract)
        except AuditError as error:
            if "SHA-256" not in str(error):
                raise
        else:
            raise RuntimeError("Turner credibility audit accepted tampering")
    print("Turner credibility audit regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
