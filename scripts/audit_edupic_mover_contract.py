#!/usr/bin/env python3
"""Audit the 1D electrostatic mover contract against pinned eduPIC 1.0."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


PINNED_EDUPIC_SHA256 = (
    "7c7679c0f0c98844940ea911bbb7581ec33f818e8d14427c9837ffdcf1ecea41"
)


class AuditError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(source: str, fragment: str, label: str) -> None:
    if fragment not in source:
        raise AuditError(f"missing {label} contract fragment")


def phase_metrics(steps_per_cycle: int, voltage_v: float) -> dict[str, float]:
    if steps_per_cycle <= 0:
        raise AuditError("steps_per_cycle must be positive")
    phase_advance = 2.0 * math.pi / steps_per_cycle
    normalized_maximum_difference = 2.0 * math.sin(phase_advance / 2.0)
    return {
        "phase_advance_rad": phase_advance,
        "rf_cycle_fraction": 1.0 / steps_per_cycle,
        "maximum_normalized_voltage_difference": normalized_maximum_difference,
        "maximum_voltage_difference_V":
            abs(voltage_v) * normalized_maximum_difference,
    }


def audit(native: Path, field_solver: Path, simulation: Path,
          steps_per_cycle: int, voltage_v: float) -> dict[str, object]:
    if sha256(native) != PINNED_EDUPIC_SHA256:
        raise AuditError("eduPIC source does not match the pinned source hash")
    native_text = native.read_text(encoding="utf-8")
    field_text = field_solver.read_text(encoding="utf-8")
    simulation_text = simulation.read_text(encoding="utf-8")

    require(native_text, "Time += DT_E;", "eduPIC pre-solve time advance")
    require(native_text, "solve_Poisson(rho,Time);", "eduPIC field solve")
    require(native_text,
            "e_x = c1 * efield[p] + c2 * efield[p+1];",
            "eduPIC CIC interpolation")
    require(native_text, "vx_e[k] -= e_x * FACTOR_E;", "eduPIC kick")
    require(native_text, "x_e[k]  += vx_e[k] * DT_E;", "eduPIC drift")
    require(native_text,
            "rho1[0]     * DX / (2.0 * EPSILON0)",
            "eduPIC left boundary half-cell correction")
    require(native_text,
            "rho1[N_G-1] * DX / (2.0 * EPSILON0)",
            "eduPIC right boundary half-cell correction")

    require(field_text,
            "E[i] * (1.0 - f) + E[i+1] * f",
            "AuroraPIC CIC interpolation")
    require(field_text,
            "E[i] = -(phi[i+1] - phi[i-1]) / (2.0 * grid.dx())",
            "AuroraPIC centered interior field")
    endpoint_half_cell_corrected = (
        "rho[0] * grid.dx() / (2.0 * permittivity_)" in field_text and
        "rho[n-1] * grid.dx() / (2.0 * permittivity_)" in field_text)
    if not endpoint_half_cell_corrected:
        require(field_text,
                "E[0] = -(phi[1] - phi[0]) / grid.dx();",
                "AuroraPIC left boundary field")
        require(field_text,
                "E[n-1] = -(phi[n-1] - phi[n-2]) / grid.dx();",
                "AuroraPIC right boundary field")
    require(simulation_text,
            "p, interpolate_electric(grid_, p.x),\n"
            "                        qm, timestep);\n"
            "                    drift_leapfrog(p, timestep);",
            "AuroraPIC kick-drift ordering")
    require(simulation_text, "deposit_and_solve(time_ + cfg_.dt);",
            "AuroraPIC end-step field solve")

    return {
        "schema_version": 1,
        "case_id": "edupic-1.0-default-argon-ccp",
        "scope": "static_1d_electrostatic_mover_contract_audit",
        "inputs": {
            "pinned_edupic_source_sha256": sha256(native),
            "aurorapic_field_solver_sha256": sha256(field_solver),
            "aurorapic_simulation_sha256": sha256(simulation),
            "steps_per_rf_cycle": steps_per_cycle,
            "electrode_voltage_amplitude_V": voltage_v,
        },
        "equivalent_contracts": {
            "interior_nodal_field":
                "centered potential difference at interior nodes",
            "particle_field_interpolation":
                "linear cloud-in-cell interpolation between adjacent nodes",
            "electron_kick": "v_half <- v_half + (q/m) E dt",
            "drift": "x <- x + v_half dt",
            "ordering": "charge deposit, Poisson solve, kick, drift",
        },
        "differences": {
            "rf_field_time": {
                "edupic": "advance time, solve E(t+dt), then kick",
                "aurorapic":
                    "kick with retained E(t), drift, then solve E(t+dt)",
                **phase_metrics(steps_per_cycle, voltage_v),
            },
            "electrode_node_field": {
                "edupic":
                    "one-sided potential gradient plus charged half-cell Gauss correction",
                "aurorapic": (
                    "one-sided potential gradient plus charged half-cell Gauss correction"
                    if endpoint_half_cell_corrected else
                    "one-sided potential gradient only"),
                "left_missing_term": "-rho[0] dx/(2 epsilon)",
                "right_missing_term": "+rho[n-1] dx/(2 epsilon)",
                "difference_present": not endpoint_half_cell_corrected,
                "direct_scope":
                    "the two wall-adjacent interpolation cells; indirect distribution effects can propagate inward",
            },
        },
        "interpretation": {
            "bulk_mover_mismatch_found": False,
            "rf_phase_offset_is_small": True,
            "phase_alignment_is_cleanly_testable_by_configuration": True,
            "boundary_field_difference_requires_a_prospective_solver_branch":
                not endpoint_half_cell_corrected,
            "boundary_half_cell_correction_present":
                endpoint_half_cell_corrected,
        },
        "claim_boundary":
            "This is a source-contract and scale audit, not a dynamic equivalence proof or validation result.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native_edupic", type=Path)
    parser.add_argument("aurorapic_field_solver", type=Path)
    parser.add_argument("aurorapic_simulation", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--steps-per-cycle", type=int, default=4000)
    parser.add_argument("--voltage-v", type=float, default=250.0)
    args = parser.parse_args()
    try:
        result = audit(args.native_edupic, args.aurorapic_field_solver,
                       args.aurorapic_simulation, args.steps_per_cycle,
                       args.voltage_v)
    except (AuditError, OSError, UnicodeError) as error:
        parser.error(str(error))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
