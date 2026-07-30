#!/usr/bin/env python3
"""Assess bounded Turner cycle reports for pre-benchmark stationarity."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile


class StationarityError(RuntimeError):
    pass


def positive_integer(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return result


def bounded_fraction(value: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise argparse.ArgumentTypeError("value must be positive and finite")
    return result


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise StationarityError(f"cannot hash {path}: {error}") from error


def load(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise StationarityError(f"cannot read {path}: {error}") from error


def relative_span(values: list[float], name: str) -> float:
    if not values or any(not math.isfinite(value) for value in values):
        raise StationarityError(f"{name} contains invalid values")
    mean = sum(values) / len(values)
    if mean == 0:
        return 0.0 if max(values) == min(values) else math.inf
    return (max(values) - min(values)) / abs(mean)


def startup_cycle(report: dict) -> dict:
    diagnostics = report.get("diagnostics", {})
    population = diagnostics.get("particle_balance", {})
    initial_electrons = population.get("initial_electrons")
    initial_ions = population.get("initial_ions")
    final_electrons = population.get("final_electrons")
    final_ions = population.get("final_ions")
    field = diagnostics.get("early_field_structure", {})
    energy = diagnostics.get("energy_and_charge", {})
    collisions = diagnostics.get("collisions", {})
    values = (
        initial_electrons, initial_ions, final_electrons, final_ions,
        field.get("final_boundary_field_max_v_m"),
        energy.get("final_total_energy_j"),
        collisions.get("collisions_electron_mcc.ionization"),
    )
    if any(not isinstance(value, (int, float)) for value in values):
        raise StationarityError("startup report lacks cycle-one diagnostics")
    return {
        "cycle": 1,
        "electron_relative_change":
            (final_electrons - initial_electrons) / initial_electrons,
        "ion_relative_change":
            (final_ions - initial_ions) / initial_ions,
        "final_electrons": final_electrons,
        "final_ions": final_ions,
        "ionization_pairs": collisions[
            "collisions_electron_mcc.ionization"
        ],
        "phase_zero_boundary_field_v_m":
            field["final_boundary_field_max_v_m"],
        "phase_zero_total_energy_j": energy["final_total_energy_j"],
    }


def horizon_cycles(report: dict) -> list[dict]:
    result = []
    for item in report.get("cycles", []):
        population = item.get("population", {})
        field = item.get("phase_zero_field", {})
        energy = item.get("energy_and_charge", {})
        values = (
            item.get("cycle"),
            population.get("electron_relative_change"),
            population.get("ion_relative_change"),
            population.get("final_electrons"),
            population.get("final_ions"),
            population.get("ionization_pairs_created"),
            field.get("boundary_field_max_v_m"),
            energy.get("final_total_energy_j"),
        )
        if any(not isinstance(value, (int, float)) for value in values):
            raise StationarityError("horizon report has incomplete cycle data")
        result.append({
            "cycle": item["cycle"],
            "electron_relative_change":
                population["electron_relative_change"],
            "ion_relative_change": population["ion_relative_change"],
            "final_electrons": population["final_electrons"],
            "final_ions": population["final_ions"],
            "ionization_pairs": population["ionization_pairs_created"],
            "phase_zero_boundary_field_v_m":
                field["boundary_field_max_v_m"],
            "phase_zero_total_energy_j": energy["final_total_energy_j"],
        })
    if not result:
        raise StationarityError("horizon report has no cycles")
    return result


def analyze(args: argparse.Namespace) -> dict:
    paths = [path.resolve() for path in args.reports]
    reports = [load(path) for path in paths]
    if reports[0].get("turner_startup_report_version") != 1:
        raise StationarityError("first report must be the one-cycle startup")
    case_id = reports[0].get("case_id")
    if case_id != "turner-helium-ccp-2013-case-1":
        raise StationarityError("unexpected Turner case identity")
    executable_hash = reports[0].get("provenance", {}).get(
        "executable_sha256"
    )
    audit_hash = reports[0].get("provenance", {}).get(
        "normalization_audit_sha256"
    )
    cycles = [startup_cycle(reports[0])]
    prior_digest = sha256(paths[0])
    for index, report in enumerate(reports[1:], 1):
        if (
            report.get("turner_horizon_report_version") != 1
            or report.get("case_id") != case_id
        ):
            raise StationarityError(f"report {index + 1} is not a Turner horizon")
        provenance = report.get("provenance", {})
        if provenance.get("prior_report_sha256") != prior_digest:
            raise StationarityError(
                f"report {index + 1} does not hash-chain to its predecessor"
            )
        if (
            provenance.get("executable_sha256") != executable_hash
            or provenance.get("normalization_audit_sha256") != audit_hash
        ):
            raise StationarityError(
                f"report {index + 1} changed executable or normalized data"
            )
        cycles.extend(horizon_cycles(report))
        prior_digest = sha256(paths[index])
    numbers = [cycle["cycle"] for cycle in cycles]
    if numbers != list(range(1, numbers[-1] + 1)):
        raise StationarityError("cycle reports are not contiguous from cycle one")
    if len(cycles) < args.window_cycles:
        raise StationarityError(
            f"need at least {args.window_cycles} cycles for the stationarity window"
        )
    window = cycles[-args.window_cycles:]
    population_change = max(
        max(
            abs(cycle["electron_relative_change"]),
            abs(cycle["ion_relative_change"]),
        )
        for cycle in window
    )
    ionization_span = relative_span(
        [cycle["ionization_pairs"] for cycle in window], "ionization"
    )
    field_span = relative_span(
        [
            cycle["phase_zero_boundary_field_v_m"]
            for cycle in window
        ],
        "phase-zero boundary field",
    )
    energy_span = relative_span(
        [cycle["phase_zero_total_energy_j"] for cycle in window],
        "phase-zero total energy",
    )
    gates = {
        "population_change": {
            "value": population_change,
            "maximum": args.max_population_change,
            "passed": population_change <= args.max_population_change,
        },
        "ionization_relative_span": {
            "value": ionization_span,
            "maximum": args.max_observable_span,
            "passed": ionization_span <= args.max_observable_span,
        },
        "boundary_field_relative_span": {
            "value": field_span,
            "maximum": args.max_observable_span,
            "passed": field_span <= args.max_observable_span,
        },
        "total_energy_relative_span": {
            "value": energy_span,
            "maximum": args.max_observable_span,
            "passed": energy_span <= args.max_observable_span,
        },
    }
    stationary = all(gate["passed"] for gate in gates.values())
    return {
        "turner_stationarity_screen_version": 1,
        "case_id": case_id,
        "scope": "pre_benchmark_stationarity_screen",
        "physics_claim": "none",
        "published_x2_applicable": False,
        "stationarity_screen_passed": stationary,
        "window": {
            "cycles": args.window_cycles,
            "start_cycle": window[0]["cycle"],
            "end_cycle": window[-1]["cycle"],
        },
        "gates": gates,
        "final_state": {
            "cycle": cycles[-1]["cycle"],
            "electrons": cycles[-1]["final_electrons"],
            "ions": cycles[-1]["final_ions"],
            "electron_fraction_of_initial":
                cycles[-1]["final_electrons"] / cycles[0]["final_electrons"]
                * (
                    1.0 + cycles[0]["electron_relative_change"]
                ),
        },
        "cycle_history": cycles,
        "provenance": {
            "reports": [
                {"path": str(path), "sha256": sha256(path)}
                for path in paths
            ],
            "executable_sha256": executable_hash,
            "normalization_audit_sha256": audit_hash,
        },
        "interpretation": (
            "screen passed; a longer confirmation window is still required"
            if stationary else
            "screen failed; continue checkpointed startup without X2 comparison"
        ),
    }


def atomic_json(path: Path, report: dict) -> None:
    if path.exists():
        raise StationarityError(f"refusing to overwrite output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="+", type=Path)
    parser.add_argument("--window-cycles", type=positive_integer, default=4)
    parser.add_argument(
        "--max-population-change", type=bounded_fraction, default=0.005
    )
    parser.add_argument(
        "--max-observable-span", type=bounded_fraction, default=0.05
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = analyze(args)
        atomic_json(args.output.resolve(), report)
    except StationarityError as error:
        print(f"Turner stationarity error: {error}", file=sys.stderr)
        return 2
    print(
        "Turner pre-benchmark stationarity screen "
        f"{'passed' if report['stationarity_screen_passed'] else 'failed'}: "
        f"{args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
