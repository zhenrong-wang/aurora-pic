#!/usr/bin/env python3
"""Split every record in a locked APS state at reciprocal macro weight."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

from export_checkpoint_particle_state import atomic_write, sha256, state_signature
from prepare_quasineutral_particle_state import (
    MAGIC, PreparationError, node_charge, read_state,
)


class SplitError(RuntimeError):
    pass


def split(species: dict[str, list[tuple[float, ...]]], factor: int) -> dict[
        str, list[tuple[float, ...]]]:
    if factor < 2:
        raise SplitError("particle split factor must be at least two")
    return {name: [record for record in records for _ in range(factor)]
            for name, records in species.items()}


def execute(args: argparse.Namespace) -> dict[str, object]:
    source = args.source.resolve()
    output = args.output.resolve()
    manifest = args.manifest.resolve()
    for target in (output, manifest):
        if target.exists():
            raise SplitError(f"refusing to overwrite: {target}")
    if sha256(source) != args.expected_source_sha256.lower():
        raise SplitError("source SHA-256 differs")
    if (args.factor < 2 or args.nodes < 3 or
            not math.isfinite(args.length) or args.length <= 0.0 or
            not math.isfinite(args.source_macro_weight) or
            args.source_macro_weight <= 0.0):
        raise SplitError("invalid particle-split contract")
    units, source_species = read_state(source)
    expected = {"electrons": args.expected_electrons,
                "ions": args.expected_ions}
    if {name: len(value) for name, value in source_species.items()} != expected:
        raise SplitError("source species counts differ")
    result_species = split(source_species, args.factor)
    total = sum(len(records) for records in result_species.values())
    lines = [
        MAGIC, "dimension 1", "velocity_dimensions 3", f"units {units}",
        "weighting species_constant", "velocity_staggering time_centered",
        f"particle_count {total}", "records",
    ]
    for name in sorted(result_species):
        for x, y, z, vx, vy, vz in result_species[name]:
            lines.append(
                f"particle {name} {x:.17g} {y:.17g} {z:.17g} "
                f"{vx:.17g} {vy:.17g} {vz:.17g}")
    lines.append("end")
    atomic_write(output, "\n".join(lines) + "\n")
    child_weight = args.source_macro_weight / args.factor
    source_charge = node_charge(
        source_species, args.source_macro_weight, args.length, args.nodes)
    result_charge = node_charge(
        result_species, child_weight, args.length, args.nodes)
    difference = [new - old for new, old in zip(result_charge, source_charge)]
    denominator = math.fsum(abs(value) for value in source_charge)
    signature_records = {
        name: [(x, vx, vy, vz) for x, _y, _z, vx, vy, vz in records]
        for name, records in result_species.items()}
    report = {
        "schema_version": 1,
        "scope": "exact_reciprocal_weight_particle_split",
        "source_sha256": args.expected_source_sha256.lower(),
        "output_sha256": sha256(output),
        "output_signature": state_signature(units, signature_records),
        "factor": args.factor,
        "source_macro_weight": args.source_macro_weight,
        "child_macro_weight": child_weight,
        "source_counts": expected,
        "output_counts": {
            name: len(records) for name, records in result_species.items()},
        "charge_preservation": {
            "signed_weighted_particle_error": math.fsum(difference),
            "node_charge_relative_l1_error": (
                math.fsum(abs(value) for value in difference) / denominator),
            "node_charge_maximum_absolute_weighted_particle_error": max(
                abs(value) for value in difference),
        },
        "claim_boundary": (
            "Exact record splitting with reciprocal macro weight preserves "
            "represented charge and phase space but creates initially "
            "coincident, velocity-identical children. Collisions can decorrelate "
            "them; this is a particle-count sensitivity initialization, not an "
            "independent statistical ensemble."),
    }
    atomic_write(manifest, json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-electrons", type=int, required=True)
    parser.add_argument("--expected-ions", type=int, required=True)
    parser.add_argument("--factor", type=int, required=True)
    parser.add_argument("--length", type=float, required=True)
    parser.add_argument("--nodes", type=int, required=True)
    parser.add_argument("--source-macro-weight", type=float, required=True)
    try:
        result = execute(parser.parse_args())
    except (SplitError, PreparationError, OSError, ValueError) as error:
        print(f"particle-state split rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
