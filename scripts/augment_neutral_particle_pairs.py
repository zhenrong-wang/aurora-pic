#!/usr/bin/env python3
"""Augment a locked 1D3V APS state with deterministic neutral pairs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

from export_checkpoint_particle_state import atomic_write, sha256, state_signature
from prepare_quasineutral_particle_state import (
    MAGIC, PreparationError, node_charge, read_state, stratified_indices,
)


class AugmentationError(RuntimeError):
    pass


def resampled_indices(size: int, count: int) -> list[int]:
    if size <= 0 or count < 0:
        raise AugmentationError("invalid resampling request")
    return [min(size - 1, ((2 * index + 1) * size) // (2 * count))
            for index in range(count)] if count else []


def paired_pool(species: dict[str, list[tuple[float, ...]]], bins: int,
                length: float) -> list[tuple[tuple[float, ...], tuple[float, ...]]]:
    grouped = {name: [[] for _ in range(bins)] for name in species}
    for name, records in species.items():
        for record in records:
            x = record[0]
            if x < 0.0 or x > length:
                raise AugmentationError(
                    "particle lies outside the declared domain")
            grouped[name][min(bins - 1, int(x / length * bins))].append(record)
    pool = []
    for bin_index in range(bins):
        electrons = sorted(grouped["electrons"][bin_index], key=lambda row: row[0])
        ions = sorted(grouped["ions"][bin_index], key=lambda row: row[0])
        count = min(len(electrons), len(ions))
        electron_indices = stratified_indices(len(electrons), count)
        ion_indices = stratified_indices(len(ions), count)
        pool.extend((electrons[electron], ions[ion])
                    for electron, ion in zip(electron_indices, ion_indices))
    if not pool:
        raise AugmentationError("source state has no pairable electron-ion bulk")
    return pool


def augment(species: dict[str, list[tuple[float, ...]]], bins: int,
            length: float, added_pairs: int) -> tuple[
                dict[str, list[tuple[float, ...]]], int]:
    pool = paired_pool(species, bins, length)
    result = {name: list(records) for name, records in species.items()}
    for index in resampled_indices(len(pool), added_pairs):
        electron, ion = pool[index]
        position = 0.5 * (electron[0] + ion[0])
        result["electrons"].append((position, *electron[1:]))
        result["ions"].append((position, *ion[1:]))
    return result, len(pool)


def execute(args: argparse.Namespace) -> dict[str, object]:
    source = args.source.resolve()
    output = args.output.resolve()
    manifest = args.manifest.resolve()
    for target in (output, manifest):
        if target.exists():
            raise AugmentationError(f"refusing to overwrite: {target}")
    if sha256(source) != args.expected_source_sha256.lower():
        raise AugmentationError("source SHA-256 does not match the locked value")
    if (not math.isfinite(args.length) or args.length <= 0.0 or
            args.bins < 2 or args.nodes < 3 or args.added_pairs <= 0 or
            not math.isfinite(args.macro_weight) or args.macro_weight <= 0.0):
        raise AugmentationError("invalid augmentation contract")
    units, source_species = read_state(source)
    expected = {
        "electrons": args.expected_electrons,
        "ions": args.expected_ions,
    }
    if {name: len(records) for name, records in source_species.items()} != expected:
        raise AugmentationError("source species counts do not match locked values")
    augmented, pool_size = augment(
        source_species, args.bins, args.length, args.added_pairs)
    total = sum(len(records) for records in augmented.values())
    lines = [
        MAGIC, "dimension 1", "velocity_dimensions 3", f"units {units}",
        "weighting species_constant", "velocity_staggering time_centered",
        f"particle_count {total}", "records",
    ]
    for name in sorted(augmented):
        for x, y, z, vx, vy, vz in augmented[name]:
            lines.append(
                f"particle {name} {x:.17g} {y:.17g} {z:.17g} "
                f"{vx:.17g} {vy:.17g} {vz:.17g}")
    lines.append("end")
    atomic_write(output, "\n".join(lines) + "\n")
    source_charge = node_charge(
        source_species, args.macro_weight, args.length, args.nodes)
    augmented_charge = node_charge(
        augmented, args.macro_weight, args.length, args.nodes)
    difference = [new - old for new, old in zip(augmented_charge, source_charge)]
    l1 = math.fsum(abs(value) for value in source_charge)
    signature_records = {
        name: [(x, vx, vy, vz) for x, _y, _z, vx, vy, vz in records]
        for name, records in augmented.items()}
    result = {
        "schema_version": 1,
        "scope": "exact_charge_neutral_pair_augmentation",
        "source_path": str(source),
        "source_sha256": args.expected_source_sha256.lower(),
        "output_path": str(output),
        "output_sha256": sha256(output),
        "output_signature": state_signature(units, signature_records),
        "spatial_bins": args.bins,
        "charge_diagnostic_nodes": args.nodes,
        "domain_length": args.length,
        "macro_weight": args.macro_weight,
        "source_counts": expected,
        "added_pairs": args.added_pairs,
        "pairing_pool_size": pool_size,
        "output_counts": {
            name: len(records) for name, records in augmented.items()},
        "charge_preservation": {
            "source_signed_charge_C_m-2": math.fsum(source_charge),
            "augmented_signed_charge_C_m-2": math.fsum(augmented_charge),
            "signed_charge_error_C_m-2": math.fsum(difference),
            "node_charge_relative_l1_error":
                math.fsum(abs(value) for value in difference) / l1,
            "node_charge_maximum_absolute_error_C_m-2": max(
                abs(value) for value in difference),
        },
        "claim_boundary": (
            "Every source record is retained and every added electron-ion "
            "pair is co-located at a common macro weight, so added grid "
            "charge cancels. Reused source velocities introduce deterministic "
            "particle correlations; this is an initialization hypothesis, "
            "not a kinetic equilibrium or validation result."),
    }
    atomic_write(manifest, json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-electrons", type=int, required=True)
    parser.add_argument("--expected-ions", type=int, required=True)
    parser.add_argument("--length", type=float, required=True)
    parser.add_argument("--bins", type=int, required=True)
    parser.add_argument("--nodes", type=int, required=True)
    parser.add_argument("--macro-weight", type=float, required=True)
    parser.add_argument("--added-pairs", type=int, required=True)
    try:
        report = execute(parser.parse_args())
    except (AugmentationError, PreparationError, OSError, ValueError) as error:
        print(f"neutral-pair augmentation rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
