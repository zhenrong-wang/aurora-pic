#!/usr/bin/env python3
"""Produce a fail-closed Turner Case 1 scientific-credibility verdict."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import sys


class AuditError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise AuditError(f"cannot read evidence {path}: {error}") from error


def load_json(path: Path, label: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise AuditError(f"cannot read {label}: {error}") from error
    require(isinstance(value, dict), f"{label} must be a JSON object")
    return value


def finite(value: object, label: str) -> float:
    require(isinstance(value, (int, float)) and math.isfinite(float(value)),
            f"{label} must be finite")
    return float(value)


def atomic_json(path: Path, value: dict[str, object]) -> None:
    require(not path.exists(), f"refusing to overwrite audit: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def audit(contract_path: Path) -> dict[str, object]:
    contract_path = contract_path.resolve()
    contract = load_json(contract_path, "audit contract")
    require(contract.get("schema_version") == 1 and
            contract.get("case_id") == "turner-helium-ccp-2013-case-1" and
            contract.get("scope") == "turner_credibility_evidence_contract",
            "unsupported audit contract")
    root = contract_path.parents[2]
    locks = contract.get("evidence")
    require(isinstance(locks, dict), "audit contract has no evidence map")
    expected_names = {
        "published_duration_ensemble", "numerical_sensitivity",
        "spatial_collision", "phase_structure", "phase_eedf",
        "energy_ledger",
    }
    require(set(locks) == expected_names,
            "audit contract evidence set is incomplete or unexpected")
    evidence: dict[str, dict[str, object]] = {}
    provenance = {}
    for name in sorted(expected_names):
        lock = locks[name]
        require(isinstance(lock, dict), f"{name} lock must be an object")
        relative = lock.get("path")
        expected_hash = lock.get("sha256")
        require(isinstance(relative, str) and isinstance(expected_hash, str),
                f"{name} lock is incomplete")
        path = (root / relative).resolve()
        require(path.is_relative_to(root), f"{name} path escapes repository")
        actual_hash = sha256(path)
        require(actual_hash == expected_hash,
                f"{name} SHA-256 differs from the audit contract")
        evidence[name] = load_json(path, name)
        provenance[name] = {"path": relative, "sha256": actual_hash}

    ensemble = evidence["published_duration_ensemble"]
    sensitivity = evidence["numerical_sensitivity"]
    collision = evidence["spatial_collision"]
    phase = evidence["phase_structure"]
    eedf = evidence["phase_eedf"]
    ledger = evidence["energy_ledger"]
    require(ensemble.get("case_id") == contract["case_id"] and
            ensemble.get("verified_members") == 3,
            "published-duration ensemble identity is invalid")
    published = ensemble.get("published_individual_run_results")
    density = ensemble.get("density_amplitude")
    require(isinstance(published, dict) and isinstance(density, dict),
            "ensemble summary is incomplete")
    pass_95 = int(published.get("accepted_95_percent_count", -1))
    pass_99 = int(published.get("accepted_99_percent_count", -1))
    fail_99 = int(published.get("failed_99_percent_count", -1))
    require(pass_95 + (3 - pass_95) == 3 and pass_99 + fail_99 == 3 and
            0 <= pass_95 <= pass_99 <= 3,
            "ensemble acceptance counts are inconsistent")
    biases_same_sign = density.get("all_member_biases_same_sign") is True
    mean_bias = finite(density.get("integrated_bias_percent_mean"),
                       "ensemble mean density bias")
    require(mean_bias > 0.0 and biases_same_sign,
            "contract expects the observed positive systematic density bias")

    require(sensitivity.get("case_id") == contract["case_id"] and
            sensitivity.get("interpretation", {}).get("conclusion") ==
            "no_tested_refinement_materially_reduced_density_bias" and
            sensitivity.get("global_observable_ranges", {}).get(
                "exact_source_loss_balance_all_variants") is True,
            "numerical-sensitivity conclusion is missing or changed")
    numerical = collision.get("numerical_closure")
    power = collision.get("power_balance_W_m-2")
    context = collision.get("turner_context")
    require(isinstance(numerical, dict) and numerical.get("passes") is True and
            numerical.get("particle_source_loss_balance_exact") is True,
            "spatial-collision numerical closure did not pass")
    require(ledger.get("acceptance", {}).get("passes") is True,
            "independent energy ledger did not pass")
    require(isinstance(power, dict) and isinstance(context, dict),
            "spatial-collision observable context is incomplete")
    require(phase.get("scope") ==
            "post_benchmark_phase_resolved_diagnostic_window",
            "phase evidence is not the required diagnostic window")
    require(eedf.get("scope") ==
            "post_benchmark_32_cycle_phase_resolved_regional_eedf" and
            eedf.get("quality", {}).get("passes") is True,
            "phase EEDF evidence is invalid")

    scalar_differences = {
        "midplane_ion_density_percent": 100.0 * finite(
            context.get("midplane_ion_density_relative_difference"),
            "midplane ion-density difference"),
        "midplane_electron_temperature_percent": 100.0 * finite(
            context.get("midplane_electron_temperature_relative_difference"),
            "midplane electron-temperature difference"),
        "electron_power_percent": 100.0 * finite(
            context.get("electron_electrical_power_relative_difference"),
            "electron-power difference"),
        "ion_power_percent": 100.0 * finite(
            context.get("ion_electrical_power_relative_difference"),
            "ion-power difference"),
        "ion_current_percent": 100.0 * finite(
            context.get("mean_ion_current_relative_difference"),
            "ion-current difference"),
    }
    numerical_integrity = True
    published_density_verified = pass_95 == 3
    systematic_density_discrepancy = biases_same_sign and fail_99 >= 2
    classification = (
        "published_case_verified" if published_density_verified else
        "partial_code_to_code_verification_systematic_density_discrepancy"
        if numerical_integrity and systematic_density_discrepancy else
        "insufficient_or_inconsistent_evidence")
    return {
        "schema_version": 1,
        "case_id": contract["case_id"],
        "scope": "audited_turner_case1_scientific_credibility_status",
        "classification": classification,
        "published_duration_density": {
            "verified_members": 3,
            "accepted_95_percent": pass_95,
            "accepted_99_percent": pass_99,
            "failed_99_percent": fail_99,
            "ensemble_mean_integrated_bias_percent": mean_bias,
            "all_biases_same_sign": biases_same_sign,
            "published_density_verification_passed":
                published_density_verified,
        },
        "numerical_integrity": {
            "independent_energy_ledger_passed": True,
            "spatial_and_phase_collision_ledgers_passed": True,
            "particle_source_loss_balance_exact": True,
            "tested_resolution_refinements_remove_bias": False,
            "passed": numerical_integrity,
        },
        "published_scalar_context": {
            "relative_differences_percent": scalar_differences,
            "formal_acceptance_ranges_published": False,
            "interpretation": "descriptive_agreement_only",
        },
        "diagnostic_localization": {
            "phase_symmetry_supported": True,
            "regional_eedf_quality_passed": True,
            "independent_phase_matched_eedf_reference_available": False,
            "remaining_discrepancy_focus":
                "electron_sheath_transport_and_collision_model_conventions",
        },
        "claims": {
            "supported": [
                "The electrostatic PIC-MCC implementation executes this case with independently closed energy and particle ledgers.",
                "Turner-compatible diagnostic windows reproduce global current and power observables descriptively closely.",
                "Ordinary particle, timestep, and grid refinement did not remove the positive density bias in the tested matrix.",
            ],
            "not_supported": [
                "AuroraPIC passes the Turner Case 1 published ion-density benchmark across independent seeds.",
                "AuroraPIC is generally validated for CCP or arbitrary plasma devices.",
                "Post-benchmark diagnostic windows can replace the published-duration acceptance result.",
            ],
        },
        "next_credibility_gate": (
            "Run a prospectively matched independent implementation comparison "
            "of phase-resolved bulk/sheath electron energy distributions and "
            "collision-channel transport; retain the published-duration ion "
            "density discrepancy as unresolved."),
        "evidence_contract_sha256": sha256(contract_path),
        "audit_tool_sha256": sha256(Path(__file__).resolve()),
        "evidence": provenance,
        "claim_boundary": (
            "This audit consolidates existing code-to-code evidence. It does "
            "not create a new acceptance threshold, reinterpret failed Turner "
            "runs, or constitute experimental validation."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("contract", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = audit(args.contract)
        if args.output:
            atomic_json(args.output.resolve(), report)
    except (AuditError, KeyError, TypeError, ValueError) as error:
        print(f"Turner credibility audit rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
