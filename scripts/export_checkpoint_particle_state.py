#!/usr/bin/env python3
"""Export live particles from a locked 1D3V checkpoint to APS v2."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import tempfile


MAGIC = "AuroraPIC-particle-state-v2"
CHECKPOINT_MAGICS = {
    "AuroraPIC-checkpoint-v14", "AuroraPIC-checkpoint-v15",
    "AuroraPIC-checkpoint-v16"}
FNV_OFFSET = 14695981039346656037
FNV_PRIME = 1099511628211


class ExportError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_uint64(value: int, item: int) -> int:
    for byte in struct.pack("<Q", item):
        value = ((value ^ byte) * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
    return value


def hash_string(value: int, item: str) -> int:
    encoded = item.encode("utf-8")
    value = hash_uint64(value, len(encoded))
    for byte in encoded:
        value = ((value ^ byte) * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
    return value


def hash_double(value: int, item: float) -> int:
    return hash_uint64(value, struct.unpack("<Q", struct.pack("<d", item))[0])


def parse_expected_species(values: list[str]) -> dict[str, int]:
    result: dict[str, int] = {}
    for value in values:
        name, separator, count_text = value.partition("=")
        if not separator or not name or name in result:
            raise ExportError("--expected-species requires unique NAME=COUNT values")
        try:
            count = int(count_text)
        except ValueError as error:
            raise ExportError("expected species count is not an integer") from error
        if count <= 0:
            raise ExportError("expected species count must be positive")
        result[name] = count
    return result


def parse_checkpoint(path: Path, expected_step: int) -> dict[str, object]:
    with path.open(encoding="utf-8") as stream:
        checkpoint_magic = stream.readline().rstrip("\n")
        if checkpoint_magic not in CHECKPOINT_MAGICS:
            raise ExportError("only 1D checkpoint v14-v16 is supported")
        dimension = None
        units = None
        velocity_dimensions = None
        step = None
        for line in stream:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "dimension" and dimension is None:
                dimension = int(fields[1])
            elif fields[0] == "units" and units is None:
                units = fields[1]
            elif fields[0] == "velocity_dimensions" and velocity_dimensions is None:
                velocity_dimensions = int(fields[1])
            elif fields[0] == "step":
                step = int(fields[1])
                break
        if (dimension, velocity_dimensions) != (1, 3):
            raise ExportError("checkpoint is not 1D3V")
        if units not in {"normalized", "si"}:
            raise ExportError("checkpoint has unsupported units")
        if step != expected_step:
            raise ExportError(
                f"checkpoint step {step} does not match expected {expected_step}")
        time_fields = stream.readline().split()
        species_count_fields = stream.readline().split()
        rng_fields = stream.readline().split()
        if len(time_fields) != 2 or time_fields[0] != "time":
            raise ExportError("checkpoint is missing time after step")
        if (len(species_count_fields) != 2 or
                species_count_fields[0] != "species_count"):
            raise ExportError("checkpoint is missing species_count after time")
        if not rng_fields or rng_fields[0] != "rng":
            raise ExportError("checkpoint is missing RNG state")
        simulation_time = float(time_fields[1])
        species_count = int(species_count_fields[1])
        species: dict[str, list[tuple[float, float, float, float]]] = {}
        stored_counts: dict[str, int] = {}
        for expected_id in range(species_count):
            header = stream.readline().split()
            if len(header) != 4 or header[0] != "species":
                raise ExportError("checkpoint species header is malformed")
            species_id = int(header[1])
            name = header[2]
            stored_count = int(header[3])
            if species_id != expected_id or not name or name in species:
                raise ExportError("checkpoint species identity is invalid")
            records: list[tuple[float, float, float, float]] = []
            for _ in range(stored_count):
                fields = stream.readline().split()
                if len(fields) != 6:
                    raise ExportError("checkpoint particle record is malformed")
                x, vx, vy, vz, _v_half = map(float, fields[:5])
                alive = int(fields[5])
                if alive not in {0, 1} or not all(
                        math.isfinite(value) for value in (x, vx, vy, vz)):
                    raise ExportError("checkpoint particle record is invalid")
                if alive:
                    records.append((x, vx, vy, vz))
            if not records:
                raise ExportError(f"checkpoint species '{name}' has no live particles")
            species[name] = records
            stored_counts[name] = stored_count
        if stream.read().strip():
            raise ExportError("checkpoint contains trailing data")
    return {
        "checkpoint_format": checkpoint_magic,
        "step": step,
        "time": simulation_time,
        "units": units,
        "species": species,
        "stored_counts": stored_counts,
    }


def state_signature(units: str,
                    species: dict[str, list[tuple[float, float, float, float]]]) -> int:
    total = sum(len(records) for records in species.values())
    value = FNV_OFFSET
    value = hash_string(value, MAGIC)
    value = hash_uint64(value, 2)
    value = hash_uint64(value, 1)
    value = hash_uint64(value, 3)
    value = hash_string(value, units)
    value = hash_uint64(value, total)
    value = hash_uint64(value, len(species))
    for name in sorted(species):
        records = species[name]
        value = hash_string(value, name)
        value = hash_uint64(value, len(records))
        for x, vx, vy, vz in records:
            for item in (x, 0.0, 0.0, vx, vy, vz):
                value = hash_double(value, item)
    return value


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ExportError(f"refusing to overwrite: {path}")
    with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", dir=path.parent,
            prefix=path.name + ".", delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())
    try:
        os.link(temporary, path)
    except FileExistsError as error:
        raise ExportError(f"refusing to overwrite: {path}") from error
    finally:
        temporary.unlink(missing_ok=True)


def execute(args: argparse.Namespace) -> dict[str, object]:
    checkpoint = args.checkpoint.resolve()
    output = args.output.resolve()
    manifest = args.manifest.resolve()
    for target in (output, manifest):
        if target.exists():
            raise ExportError(f"refusing to overwrite: {target}")
    expected_hash = args.expected_checkpoint_sha256.lower()
    if sha256(checkpoint) != expected_hash:
        raise ExportError("checkpoint SHA-256 does not match the locked value")
    expected_species = parse_expected_species(args.expected_species)
    parsed = parse_checkpoint(checkpoint, args.expected_step)
    species = parsed["species"]
    assert isinstance(species, dict)
    realized_counts = {name: len(records) for name, records in species.items()}
    if realized_counts != expected_species:
        raise ExportError(
            f"live species counts {realized_counts} do not match locked {expected_species}")
    total = sum(realized_counts.values())
    lines = [
        MAGIC,
        "dimension 1",
        "velocity_dimensions 3",
        f"units {parsed['units']}",
        "weighting species_constant",
        "velocity_staggering time_centered",
        f"particle_count {total}",
        "records",
    ]
    for name in sorted(species):
        for x, vx, vy, vz in species[name]:
            lines.append(
                f"particle {name} {x:.17g} 0 0 "
                f"{vx:.17g} {vy:.17g} {vz:.17g}")
    lines.append("end")
    atomic_write(output, "\n".join(lines) + "\n")
    result = {
        "schema_version": 1,
        "scope": "checkpoint_live_particle_state_export",
        "source_checkpoint": str(checkpoint),
        "source_checkpoint_sha256": expected_hash,
        "source_checkpoint_format": parsed["checkpoint_format"],
        "source_step": parsed["step"],
        "source_time": parsed["time"],
        "spatial_dimension": 1,
        "velocity_dimensions": 3,
        "units": parsed["units"],
        "species": {
            name: {
                "stored_records": parsed["stored_counts"][name],
                "live_records": realized_counts[name],
                "discarded_inactive_records":
                    parsed["stored_counts"][name] - realized_counts[name],
            }
            for name in sorted(species)
        },
        "particle_count": total,
        "particle_state": str(output),
        "particle_state_sha256": sha256(output),
        "particle_state_signature": state_signature(str(parsed["units"]), species),
        "claim_boundary": (
            "The export preserves live time-centered particle coordinates and "
            "velocities; it does not preserve checkpoint time, RNG, fields, or diagnostics."),
    }
    atomic_write(manifest, json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--expected-checkpoint-sha256", required=True)
    parser.add_argument("--expected-step", type=int, required=True)
    parser.add_argument("--expected-species", action="append", required=True)
    try:
        result = execute(parser.parse_args())
    except (ExportError, OSError, ValueError) as error:
        print(f"checkpoint particle-state export rejected: {error}", file=os.sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
