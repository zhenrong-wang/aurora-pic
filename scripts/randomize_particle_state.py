#!/usr/bin/env python3
"""Create a constrained independent microstate from a locked 1D APS state."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import random
import sys

from export_checkpoint_particle_state import atomic_write, sha256, state_signature
from prepare_quasineutral_particle_state import (
    MAGIC, PreparationError, read_state,
)


class RandomizationError(RuntimeError):
    pass


def cell_and_fraction(x: float, length: float,
                      nodes: int) -> tuple[int, float]:
    if not math.isfinite(x) or x < 0.0 or x > length:
        raise RandomizationError("particle position lies outside the domain")
    dx = length / (nodes - 1)
    coordinate = x / dx
    cell = min(nodes - 2, int(math.floor(coordinate)))
    return cell, coordinate - cell


def nodal_number(records: list[tuple[float, ...]], length: float,
                 nodes: int) -> list[float]:
    result = [0.0] * nodes
    for record in records:
        cell, fraction = cell_and_fraction(record[0], length, nodes)
        result[cell] += 1.0 - fraction
        result[cell + 1] += fraction
    return result


def grouped_records(records: list[tuple[float, ...]], length: float,
                    nodes: int) -> list[list[tuple[float, ...]]]:
    grouped: list[list[tuple[float, ...]]] = [
        [] for _ in range(nodes - 1)]
    for record in records:
        cell, _fraction = cell_and_fraction(record[0], length, nodes)
        grouped[cell].append(record)
    return grouped


def randomize_species(records: list[tuple[float, ...]], length: float,
                      nodes: int, rng: random.Random) -> tuple[
                          list[tuple[float, ...]], dict[str, float | int]]:
    dx = length / (nodes - 1)
    output: list[tuple[float, ...]] = []
    squared_displacement = 0.0
    changed_positions = 0
    for cell, source in enumerate(grouped_records(records, length, nodes)):
        if not source:
            continue
        fractions = [cell_and_fraction(row[0], length, nodes)[1]
                     for row in source]
        order = list(range(len(source)))
        rng.shuffle(order)
        randomized = fractions.copy()
        for offset in range(0, len(order) - 1, 2):
            first, second = order[offset:offset + 2]
            pair_sum = fractions[first] + fractions[second]
            low = max(0.0, pair_sum - 1.0)
            high = min(1.0, pair_sum)
            value = rng.uniform(low, high)
            randomized[first] = value
            randomized[second] = pair_sum - value
        velocities = [row[3:] for row in source]
        rng.shuffle(velocities)
        for index, (record, velocity) in enumerate(zip(source, velocities)):
            x = (cell + randomized[index]) * dx
            if cell == nodes - 2:
                x = min(x, length)
            displacement = x - record[0]
            squared_displacement += displacement * displacement
            changed_positions += int(x != record[0])
            output.append((x, 0.0, 0.0, *velocity))
    if len(output) != len(records):
        raise RandomizationError("randomization changed the particle count")
    return output, {
        "particles": len(records),
        "changed_positions": changed_positions,
        "rms_position_displacement_m": math.sqrt(
            squared_displacement / len(records)) if records else 0.0,
    }


def maximum_relative_nodal_error(source: list[float],
                                 result: list[float]) -> tuple[float, float]:
    difference = [new - old for new, old in zip(result, source)]
    denominator = math.fsum(abs(value) for value in source)
    relative_l1 = (math.fsum(abs(value) for value in difference) / denominator
                   if denominator > 0.0 else 0.0)
    return relative_l1, max((abs(value) for value in difference), default=0.0)


def cell_velocity_fingerprint(records: list[tuple[float, ...]], length: float,
                              nodes: int) -> list[list[tuple[float, ...]]]:
    return [sorted(row[3:] for row in cell) for cell in
            grouped_records(records, length, nodes)]


def randomize_state(species: dict[str, list[tuple[float, ...]]],
                    length: float, nodes: int, seed: int) -> tuple[
                        dict[str, list[tuple[float, ...]]], dict[str, object]]:
    if nodes < 3 or not math.isfinite(length) or length <= 0.0:
        raise RandomizationError("invalid randomization geometry")
    if seed < 0 or seed > 0xFFFFFFFF:
        raise RandomizationError("seed must be an unsigned 32-bit integer")
    rng = random.Random(seed)
    result: dict[str, list[tuple[float, ...]]] = {}
    report: dict[str, object] = {}
    for name in sorted(species):
        source = species[name]
        randomized, item = randomize_species(source, length, nodes, rng)
        if (cell_velocity_fingerprint(source, length, nodes) !=
                cell_velocity_fingerprint(randomized, length, nodes)):
            raise RandomizationError("cellwise velocity tuples were not preserved")
        relative_l1, maximum = maximum_relative_nodal_error(
            nodal_number(source, length, nodes),
            nodal_number(randomized, length, nodes))
        item.update({
            "nodal_number_relative_l1_error": relative_l1,
            "nodal_number_maximum_absolute_error": maximum,
            "cellwise_velocity_tuple_multisets_preserved": True,
        })
        if relative_l1 > 1.0e-14 or maximum > 1.0e-10:
            raise RandomizationError("nodal number density was not preserved")
        result[name] = randomized
        report[name] = item
    return result, report


def execute(args: argparse.Namespace) -> dict[str, object]:
    source = args.source.resolve()
    output = args.output.resolve()
    manifest = args.manifest.resolve()
    for target in (output, manifest):
        if target.exists():
            raise RandomizationError(f"refusing to overwrite: {target}")
    if sha256(source) != args.expected_source_sha256.lower():
        raise RandomizationError("source SHA-256 differs")
    units, source_species = read_state(source)
    expected = {"electrons": args.expected_electrons,
                "ions": args.expected_ions}
    if {name: len(records) for name, records in source_species.items()} != expected:
        raise RandomizationError("source species counts differ")
    result_species, preservation = randomize_state(
        source_species, args.length, args.nodes, args.seed)
    lines = [
        MAGIC, "dimension 1", "velocity_dimensions 3", f"units {units}",
        "weighting species_constant", "velocity_staggering time_centered",
        f"particle_count {sum(map(len, result_species.values()))}", "records",
    ]
    for name in sorted(result_species):
        for x, y, z, vx, vy, vz in result_species[name]:
            lines.append(
                f"particle {name} {x:.17g} {y:.17g} {z:.17g} "
                f"{vx:.17g} {vy:.17g} {vz:.17g}")
    lines.append("end")
    atomic_write(output, "\n".join(lines) + "\n")
    signature_records = {
        name: [(x, vx, vy, vz) for x, _y, _z, vx, vy, vz in records]
        for name, records in result_species.items()}
    result = {
        "schema_version": 1,
        "scope": "constrained_independent_particle_microstate",
        "source_sha256": args.expected_source_sha256.lower(),
        "output_sha256": sha256(output),
        "output_signature": state_signature(units, signature_records),
        "seed": args.seed,
        "geometry": {"length_m": args.length, "nodes": args.nodes},
        "counts": expected,
        "preservation": preservation,
        "randomization_contract": {
            "cell_occupancy_preserved": True,
            "cell_fraction_sum_preserved": True,
            "nodal_CIC_number_density_preserved_to_roundoff": True,
            "cellwise_velocity_tuple_multisets_preserved": True,
            "subcell_positions_redrawn_under_pair_sum_constraints": True,
            "velocity_tuples_independently_re_paired_within_cell": True,
        },
        "claim_boundary": (
            "This is an independent constrained randomization conditional on "
            "the source cell populations, CIC first moments, and cellwise "
            "empirical velocity tuples. It varies microscopic subcell position "
            "and position-velocity pairing, but is not an independent draw from "
            "an unknown continuous distribution."),
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
    parser.add_argument("--nodes", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    try:
        result = execute(parser.parse_args())
    except (RandomizationError, PreparationError, OSError, ValueError) as error:
        print(f"particle-state randomization rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
