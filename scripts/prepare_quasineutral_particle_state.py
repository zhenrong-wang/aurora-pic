#!/usr/bin/env python3
"""Build a deterministic quasi-neutral-bulk APS v2 warm-start state."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

from export_checkpoint_particle_state import atomic_write, sha256, state_signature


MAGIC = "AuroraPIC-particle-state-v2"


class PreparationError(RuntimeError):
    pass


def read_state(path: Path) -> tuple[str, dict[str, list[tuple[float, ...]]]]:
    with path.open(encoding="utf-8") as stream:
        required = [
            MAGIC, "dimension 1", "velocity_dimensions 3",
        ]
        for expected in required:
            if stream.readline().rstrip("\n") != expected:
                raise PreparationError(
                    f"particle state does not match required header {expected!r}")
        units_fields = stream.readline().split()
        if len(units_fields) != 2 or units_fields[0] != "units":
            raise PreparationError("particle state is missing units")
        units = units_fields[1]
        if units not in {"si", "normalized"}:
            raise PreparationError("particle state has unsupported units")
        for expected in (
            "weighting species_constant",
            "velocity_staggering time_centered",
        ):
            if stream.readline().rstrip("\n") != expected:
                raise PreparationError(
                    f"particle state does not match required header {expected!r}")
        count_fields = stream.readline().split()
        if len(count_fields) != 2 or count_fields[0] != "particle_count":
            raise PreparationError("particle state is missing particle_count")
        declared_count = int(count_fields[1])
        if stream.readline().strip() != "records":
            raise PreparationError("particle state is missing records")
        species: dict[str, list[tuple[float, ...]]] = {
            "electrons": [], "ions": []}
        for index in range(declared_count):
            fields = stream.readline().split()
            if len(fields) != 8 or fields[0] != "particle":
                raise PreparationError(f"invalid particle record {index}")
            name = fields[1]
            if name not in species:
                raise PreparationError(f"unexpected species {name!r}")
            values = tuple(map(float, fields[2:]))
            if not all(math.isfinite(value) for value in values):
                raise PreparationError(f"non-finite particle record {index}")
            x, y, z, _vx, _vy, _vz, = values
            if y != 0.0 or z != 0.0:
                raise PreparationError("1D state has active transverse position")
            species[name].append(values)
        if stream.readline().strip() != "end" or stream.read().strip():
            raise PreparationError("particle state has missing end or trailing data")
    return units, species


def stratified_indices(size: int, count: int) -> list[int]:
    if count < 0 or count > size:
        raise PreparationError("invalid stratified selection size")
    if count == 0:
        return []
    result = [min(size - 1, ((2 * index + 1) * size) // (2 * count))
              for index in range(count)]
    if len(set(result)) != count:
        raise PreparationError("stratified selection produced duplicate indices")
    return result


def node_charge(records: dict[str, list[tuple[float, ...]]],
                weight: float, length: float, nodes: int) -> list[float]:
    result = [0.0] * nodes
    dx = length / (nodes - 1)
    for name, sign in (("electrons", -1.0), ("ions", 1.0)):
        for record in records[name]:
            coordinate = min(max(record[0], 0.0), length) / dx
            cell = min(nodes - 2, int(math.floor(coordinate)))
            fraction = coordinate - cell
            result[cell] += sign * weight * (1.0 - fraction)
            result[cell + 1] += sign * weight * fraction
    return result


def transform(species: dict[str, list[tuple[float, ...]]],
              bins: int, length: float,
              weight_factor: float) -> tuple[dict[str, list[tuple[float, ...]]], dict[str, object]]:
    grouped = {
        name: [[] for _ in range(bins)] for name in species}
    for name, records in species.items():
        for record in records:
            x = record[0]
            if x < 0.0 or x > length:
                raise PreparationError("particle lies outside the declared domain")
            grouped[name][min(bins - 1, int(x / length * bins))].append(record)
    output = {"electrons": [], "ions": []}
    bin_report = []
    residual_bins = {
        "electrons": [[] for _ in range(bins)],
        "ions": [[] for _ in range(bins)],
    }
    for bin_index in range(bins):
        electrons = sorted(grouped["electrons"][bin_index], key=lambda row: row[0])
        ions = sorted(grouped["ions"][bin_index], key=lambda row: row[0])
        pair_count = min(len(electrons), len(ions))
        selected = {}
        residual = {}
        for name, records in (("electrons", electrons), ("ions", ions)):
            paired_indices = set(stratified_indices(len(records), pair_count))
            selected[name] = [record for index, record in enumerate(records)
                              if index in paired_indices]
            residual[name] = [record for index, record in enumerate(records)
                              if index not in paired_indices]
        for electron, ion in zip(selected["electrons"], selected["ions"]):
            position = 0.5 * (electron[0] + ion[0])
            output["electrons"].append((position, *electron[1:]))
            output["ions"].append((position, *ion[1:]))
        for name in ("electrons", "ions"):
            residual_bins[name][bin_index] = residual[name]
        bin_report.append({
            "bin": bin_index,
            "source_electrons": len(electrons),
            "source_ions": len(ions),
            "paired_bulk_per_species": pair_count,
            "retained_residual_electrons": 0,
            "retained_residual_ions": 0,
        })
    for name in ("electrons", "ions"):
        cumulative_source = 0
        cumulative_retained = 0
        for bin_index, records in enumerate(residual_bins[name]):
            cumulative_source += len(records)
            target_cumulative = math.floor(
                cumulative_source / weight_factor + 0.5)
            target = target_cumulative - cumulative_retained
            indices = stratified_indices(len(records), target)
            retained = [records[index] for index in indices]
            output[name].extend(retained)
            bin_report[bin_index][
                f"retained_residual_{name}"] = len(retained)
            cumulative_retained = target_cumulative
    return output, {
        "bins": bin_report,
        "paired_bulk_per_species": sum(
            row["paired_bulk_per_species"] for row in bin_report),
        "source_residual_particles": sum(
            abs(row["source_electrons"] - row["source_ions"])
            for row in bin_report),
        "retained_residual_particles": sum(
            row["retained_residual_electrons"] +
            row["retained_residual_ions"] for row in bin_report),
    }


def execute(args: argparse.Namespace) -> dict[str, object]:
    source = args.source.resolve()
    output = args.output.resolve()
    manifest = args.manifest.resolve()
    for target in (output, manifest):
        if target.exists():
            raise PreparationError(f"refusing to overwrite: {target}")
    if sha256(source) != args.expected_source_sha256.lower():
        raise PreparationError("source SHA-256 does not match the locked value")
    if (not math.isfinite(args.length) or args.length <= 0.0 or
            args.bins < 2 or args.nodes < 3 or
            not math.isfinite(args.weight_factor) or args.weight_factor <= 1.0):
        raise PreparationError("invalid transform geometry or weight factor")
    units, source_species = read_state(source)
    expected = {
        "electrons": args.expected_electrons,
        "ions": args.expected_ions,
    }
    if {name: len(records) for name, records in source_species.items()} != expected:
        raise PreparationError("source species counts do not match locked values")
    transformed, decomposition = transform(
        source_species, args.bins, args.length, args.weight_factor)
    total = sum(len(records) for records in transformed.values())
    lines = [
        MAGIC, "dimension 1", "velocity_dimensions 3", f"units {units}",
        "weighting species_constant", "velocity_staggering time_centered",
        f"particle_count {total}", "records",
    ]
    for name in sorted(transformed):
        for x, y, z, vx, vy, vz in transformed[name]:
            lines.append(
                f"particle {name} {x:.17g} {y:.17g} {z:.17g} "
                f"{vx:.17g} {vy:.17g} {vz:.17g}")
    lines.append("end")
    atomic_write(output, "\n".join(lines) + "\n")
    source_charge = node_charge(source_species, 1.0, args.length, args.nodes)
    transformed_charge = node_charge(
        transformed, args.weight_factor, args.length, args.nodes)
    difference = [new - old for new, old in zip(transformed_charge, source_charge)]
    l1 = math.fsum(abs(value) for value in source_charge)
    charge_metrics = {
        "source_signed_macro_charge": math.fsum(source_charge),
        "transformed_signed_source_weight_equivalents":
            math.fsum(transformed_charge),
        "signed_charge_error": math.fsum(difference),
        "node_charge_relative_l1_error":
            math.fsum(abs(value) for value in difference) / l1,
        "node_charge_relative_l2_error": math.sqrt(
            math.fsum(value * value for value in difference) /
            math.fsum(value * value for value in source_charge)),
    }
    signature_records = {
        name: [(x, vx, vy, vz) for x, _y, _z, vx, vy, vz in records]
        for name, records in transformed.items()}
    result = {
        "schema_version": 1,
        "scope": "quasineutral_bulk_sheath_preserving_particle_transform",
        "source_path": str(source),
        "source_sha256": args.expected_source_sha256.lower(),
        "output_path": str(output),
        "output_sha256": sha256(output),
        "output_signature": state_signature(units, signature_records),
        "spatial_bins": args.bins,
        "charge_diagnostic_nodes": args.nodes,
        "domain_length": args.length,
        "source_macro_weight": args.source_weight,
        "transformed_macro_weight": args.source_weight * args.weight_factor,
        "weight_factor": args.weight_factor,
        "source_counts": expected,
        "output_counts": {
            name: len(records) for name, records in transformed.items()},
        "decomposition": decomposition,
        "charge_preservation": charge_metrics,
        "claim_boundary": (
            "The transform scales a binwise paired quasi-neutral bulk and "
            "approximately retains source residual charge by deterministic "
            "stratified sampling. It is an initialization hypothesis, not a "
            "self-consistent sheath solution or validation result."),
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
    parser.add_argument("--source-weight", type=float, required=True)
    parser.add_argument("--weight-factor", type=float, required=True)
    try:
        report = execute(parser.parse_args())
    except (PreparationError, OSError, ValueError) as error:
        print(f"quasi-neutral state preparation rejected: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
