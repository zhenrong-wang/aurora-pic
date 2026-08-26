#!/usr/bin/env python3
"""Export an eduPIC 1.0 binary checkpoint as half-step APS v3."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import tempfile

from export_checkpoint_particle_state import (
    FNV_OFFSET, MAGIC_V3, atomic_write, hash_double, hash_string, hash_uint64,
)


class ExportError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_nonnegative_integer(value: float, name: str) -> int:
    if not math.isfinite(value) or value < 0.0 or value != math.floor(value):
        raise ExportError(f"eduPIC checkpoint {name} is not an exact nonnegative integer")
    return int(value)


def parse_layout(data: bytes) -> dict[str, object]:
    if len(data) < 32 or len(data) % 8:
        raise ExportError("eduPIC checkpoint size is invalid")
    time, cycles_value, electrons_value = struct.unpack_from("<ddd", data)
    cycles = exact_nonnegative_integer(cycles_value, "cycle count")
    electrons = exact_nonnegative_integer(electrons_value, "electron count")
    ion_count_offset = 24 + 32 * electrons
    if ion_count_offset + 8 > len(data):
        raise ExportError("eduPIC checkpoint electron arrays are truncated")
    ions = exact_nonnegative_integer(
        struct.unpack_from("<d", data, ion_count_offset)[0], "ion count")
    expected_size = 8 * (4 + 4 * (electrons + ions))
    if len(data) != expected_size:
        raise ExportError("eduPIC checkpoint arrays or trailing bytes differ")
    if not math.isfinite(time) or time < 0.0:
        raise ExportError("eduPIC checkpoint time is invalid")
    return {
        "time_s": time, "cycles": cycles,
        "counts": {"electrons": electrons, "ions": ions},
        "offsets": {
            "electrons": (24, 24 + 8 * electrons, 24 + 16 * electrons,
                          24 + 24 * electrons),
            "ions": (ion_count_offset + 8,
                     ion_count_offset + 8 + 8 * ions,
                     ion_count_offset + 8 + 16 * ions,
                     ion_count_offset + 8 + 24 * ions),
        },
    }


def value_at(data: bytes, offset: int, index: int) -> float:
    value = struct.unpack_from("<d", data, offset + 8 * index)[0]
    if not math.isfinite(value):
        raise ExportError("eduPIC checkpoint contains a non-finite particle value")
    return value


def semantic_signature(data: bytes, layout: dict[str, object]) -> int:
    counts = layout["counts"]
    offsets = layout["offsets"]
    total = sum(counts.values())
    value = FNV_OFFSET
    value = hash_string(value, MAGIC_V3)
    value = hash_uint64(value, 3)
    value = hash_uint64(value, 1)
    value = hash_uint64(value, 3)
    value = hash_string(value, "si")
    value = hash_string(value, "leapfrog_half_step")
    value = hash_uint64(value, total)
    value = hash_uint64(value, 2)
    for species in ("electrons", "ions"):
        count = counts[species]
        value = hash_string(value, species)
        value = hash_uint64(value, count)
        for index in range(count):
            x, vx, vy, vz = (value_at(data, offset, index)
                              for offset in offsets[species])
            for item in (x, 0.0, 0.0, vx, vy, vz):
                value = hash_double(value, item)
    return value


def write_state(path: Path, data: bytes, layout: dict[str, object]) -> None:
    if path.exists():
        raise ExportError(f"refusing to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    counts = layout["counts"]
    offsets = layout["offsets"]
    with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", dir=path.parent,
            prefix=path.name + ".", delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(
            f"{MAGIC_V3}\n"
            "dimension 1\nvelocity_dimensions 3\nunits si\n"
            "weighting species_constant\n"
            "velocity_staggering leapfrog_half_step\n"
            f"particle_count {sum(counts.values())}\nrecords\n")
        for species in ("electrons", "ions"):
            for index in range(counts[species]):
                x, vx, vy, vz = (value_at(data, offset, index)
                                  for offset in offsets[species])
                stream.write(
                    f"particle {species} {x:.17g} 0 0 "
                    f"{vx:.17g} {vy:.17g} {vz:.17g}\n")
        stream.write("end\n")
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
    if output.exists() or manifest.exists():
        raise ExportError("refusing to overwrite an output or manifest")
    actual_hash = sha256(checkpoint)
    if actual_hash != args.expected_checkpoint_sha256.lower():
        raise ExportError("eduPIC checkpoint SHA-256 does not match")
    data = checkpoint.read_bytes()
    layout = parse_layout(data)
    if args.expected_cycles is not None and layout["cycles"] != args.expected_cycles:
        raise ExportError("eduPIC checkpoint cycle count does not match")
    if args.expected_time_s is not None and not math.isclose(
            layout["time_s"], args.expected_time_s, rel_tol=0.0,
            abs_tol=max(1e-30, abs(args.expected_time_s) * 1e-14)):
        raise ExportError("eduPIC checkpoint time does not match")
    if args.expected_electrons is not None and \
            layout["counts"]["electrons"] != args.expected_electrons:
        raise ExportError("eduPIC electron count does not match")
    if args.expected_ions is not None and \
            layout["counts"]["ions"] != args.expected_ions:
        raise ExportError("eduPIC ion count does not match")
    signature = semantic_signature(data, layout)
    write_state(output, data, layout)
    result = {
        "schema_version": 1,
        "scope": "edupic_checkpoint_half_step_particle_state_export",
        "source_checkpoint": str(checkpoint),
        "source_checkpoint_sha256": actual_hash,
        "source_format": "eduPIC-1.0-native-binary",
        "source_time_s": layout["time_s"],
        "source_cycles": layout["cycles"],
        "species_counts": layout["counts"],
        "particle_count": sum(layout["counts"].values()),
        "particle_state": str(output),
        "particle_state_sha256": sha256(output),
        "particle_state_signature": signature,
        "particle_state_version": 3,
        "velocity_staggering": "leapfrog_half_step",
        "claim_boundary": (
            "The export preserves eduPIC positions and all stored velocity "
            "components exactly as APS v3 leapfrog-half-step records. It does "
            "not preserve eduPIC RNG state, collision counters, field history, "
            "or diagnostic accumulators."),
    }
    atomic_write(manifest, json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--expected-checkpoint-sha256", required=True)
    parser.add_argument("--expected-cycles", type=int)
    parser.add_argument("--expected-time-s", type=float)
    parser.add_argument("--expected-electrons", type=int)
    parser.add_argument("--expected-ions", type=int)
    args = parser.parse_args()
    try:
        result = execute(args)
    except (ExportError, OSError, ValueError, struct.error) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
