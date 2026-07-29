#!/usr/bin/env python3
"""Plan or checksum-lock an externally acquired Hall reference artifact."""

from __future__ import annotations

import argparse
import configparser
from datetime import date
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile


SHA256 = re.compile(r"^[0-9a-f]{64}$")
LARGE_FILE_ACK = "I_UNDERSTAND_THIS_MAY_HASH_A_VERY_LARGE_FILE"


class SourceLockError(RuntimeError):
    pass


def load_registry(
    path: Path,
) -> tuple[configparser.ConfigParser, configparser.SectionProxy]:
    parser = configparser.ConfigParser(
        interpolation=None, strict=True, empty_lines_in_values=False
    )
    try:
        parser.read_string("[global]\n" + path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, configparser.Error) as error:
        raise SourceLockError(f"cannot read source registry {path}: {error}") from error
    if parser["global"].get("source_registry_version") != "1":
        raise SourceLockError("source registry version must be 1")
    return parser, parser["global"]


def atomic_json(path: Path, value: dict[str, object]) -> None:
    if path.exists():
        raise SourceLockError(f"refusing to overwrite existing lock: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
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


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            while block := stream.read(1024 * 1024):
                digest.update(block)
    except OSError as error:
        raise SourceLockError(f"cannot hash artifact {path}: {error}") from error
    return digest.hexdigest()


def lock(args: argparse.Namespace) -> dict[str, object]:
    registry_path = args.registry.resolve()
    parser, global_section = load_registry(registry_path)
    section_name = f"source.{args.source}"
    if section_name not in parser:
        available = ", ".join(
            name.removeprefix("source.")
            for name in parser.sections()
            if name.startswith("source.")
        )
        raise SourceLockError(
            f"unknown source {args.source!r}; available sources: {available}"
        )
    source = parser[section_name]
    required = (
        "variant", "kind", "cells_x", "cells_y", "landing_url", "doi",
        "license", "acquisition", "artifact_name", "display_size",
    )
    missing = [key for key in required if not source.get(key, "").strip()]
    if missing:
        raise SourceLockError(f"[{section_name}] is missing {missing}")
    try:
        cells_x = source.getint("cells_x")
        cells_y = source.getint("cells_y")
    except ValueError as error:
        raise SourceLockError("source cell counts must be integers") from error
    if cells_x <= 0 or cells_y <= 0:
        raise SourceLockError("source cell counts must be positive")

    registry_hash = hashlib.sha256(registry_path.read_bytes()).hexdigest()
    result: dict[str, object] = {
        "source_lock_version": 1,
        "case_id": global_section["case_id"],
        "source_id": args.source,
        "variant": source["variant"],
        "kind": source["kind"],
        "cells_x": cells_x,
        "cells_y": cells_y,
        "landing_url": source["landing_url"],
        "doi": source["doi"],
        "license": source["license"],
        "license_url": source.get("license_url", ""),
        "acquisition": source["acquisition"],
        "expected_artifact_name": source["artifact_name"],
        "display_size": source["display_size"],
        "file_set_id": source.get("file_set_id", ""),
        "registry": str(registry_path),
        "registry_sha256": registry_hash,
        "lock_date": date.today().isoformat(),
        "artifact_locked": False,
        "warnings": [
            "This tool never downloads reference data.",
            "A locally computed hash proves byte identity after locking, not "
            "that the first acquisition came from the declared repository.",
        ],
    }
    if args.artifact is None:
        result["status"] = "external_acquisition_required"
        return result

    artifact = args.artifact.resolve()
    try:
        size = artifact.stat().st_size
    except OSError as error:
        raise SourceLockError(f"cannot inspect artifact {artifact}: {error}") from error
    if not artifact.is_file():
        raise SourceLockError(f"artifact is not a regular file: {artifact}")
    expected_name = source["artifact_name"]
    if expected_name != "externally_supplied" and artifact.name != expected_name:
        raise SourceLockError(
            f"artifact must be named {expected_name!r}, got {artifact.name!r}"
        )
    if size > args.maximum_bytes and args.acknowledge_large_file_hash != LARGE_FILE_ACK:
        raise SourceLockError(
            f"artifact is {size} bytes, above the {args.maximum_bytes}-byte "
            "default hash limit; explicitly pass --acknowledge-large-file-hash "
            f"{LARGE_FILE_ACK}"
        )
    expected_hash = args.expected_sha256
    if expected_hash is not None and not SHA256.fullmatch(expected_hash):
        raise SourceLockError("--expected-sha256 must be lowercase SHA-256")
    digest = file_sha256(artifact)
    if expected_hash is not None and digest != expected_hash:
        raise SourceLockError(
            f"artifact SHA-256 mismatch: expected {expected_hash}, got {digest}"
        )
    result.update({
        "artifact_locked": True,
        "artifact_path": str(artifact),
        "artifact_name": artifact.name,
        "artifact_bytes": size,
        "artifact_sha256": digest,
        "expected_sha256": expected_hash,
        "repository_checksum_verified": expected_hash is not None,
        "status": (
            "checksum_verified"
            if expected_hash is not None
            else "locally_hashed_repository_checksum_unavailable"
        ),
    })
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plan or lock an externally acquired Hall source"
    )
    parser.add_argument("registry", type=Path)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--artifact", type=Path)
    parser.add_argument("--expected-sha256")
    parser.add_argument("--maximum-bytes", type=int, default=64 * 1024 * 1024)
    parser.add_argument("--acknowledge-large-file-hash")
    args = parser.parse_args()
    if args.maximum_bytes <= 0:
        parser.error("--maximum-bytes must be positive")
    return args


def main() -> int:
    args = parse_args()
    try:
        result = lock(args)
        atomic_json(args.output, result)
    except SourceLockError as error:
        print(f"Hall source lock error: {error}", file=sys.stderr)
        return 2
    print(f"Hall source lock written without downloading data: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
