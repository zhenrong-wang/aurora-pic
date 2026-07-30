#!/usr/bin/env python3
"""Verify a local Turner CCP publisher supplement without redistributing it."""

from __future__ import annotations

import argparse
import configparser
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import sys
import zipfile


SHA256 = re.compile(r"^[0-9a-f]{64}$")
PROCESS_TYPES = {"ELASTIC", "EXCITATION", "IONIZATION"}


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def pairs(value: str, context: str) -> dict[str, int]:
    result: dict[str, int] = {}
    try:
        for item in value.split(","):
            key, count = item.strip().split(":", 1)
            result[key] = int(count)
    except (ValueError, TypeError) as error:
        raise VerificationError(f"invalid {context}: {value!r}") from error
    require(all(count > 0 for count in result.values()),
            f"{context} counts must be positive")
    return result


def parse_registry(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(
        interpolation=None, strict=True, empty_lines_in_values=False
    )
    try:
        parser.read_string("[global]\n" + path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, configparser.Error) as error:
        raise VerificationError(f"cannot read registry {path}: {error}") from error
    require(parser["global"].get("source_registry_version") == "1",
            "source registry version must be 1")
    require("source.publisher_supplement" in parser,
            "registry has no publisher supplement")
    return parser


def verify_electron(text: str, expected: str) -> dict[str, int]:
    actual = {name: 0 for name in PROCESS_TYPES}
    for line in text.splitlines():
        value = line.strip()
        if value in actual:
            actual[value] += 1
    wanted = pairs(expected, "process_counts")
    require(actual == wanted,
            f"electron process counts differ: expected {wanted}, got {actual}")
    require("Biagi-v7.1" in text and "Generated on 28 Nov 2012." in text,
            "electron table lost the prescribed Biagi 7.1 identity")
    return actual


def numeric_rows(text: str, columns: int) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in text.splitlines():
        value = line.strip()
        if not value or value.startswith("#"):
            continue
        try:
            row = [float(item) for item in value.split()]
        except ValueError as error:
            raise VerificationError(f"invalid numeric row: {value!r}") from error
        require(len(row) == columns,
                f"expected {columns} columns, got {len(row)}")
        rows.append(row)
    return rows


def verify_ion(text: str, expected_rows: int) -> dict[str, object]:
    rows = numeric_rows(text, 3)
    require(len(rows) == expected_rows,
            f"ion row count differs: expected {expected_rows}, got {len(rows)}")
    require(rows[0][0] == 0.0 and rows[-1][0] == 10000.0,
            "ion center-of-mass energy range must be 0--10000 eV")
    require(all(a[0] < b[0] for a, b in zip(rows, rows[1:])),
            "ion energies must be strictly increasing")
    require(all(value >= 0.0 for row in rows for value in row),
            "ion table contains a negative value")
    require("Centre of mass energy" in text
            and "Isotropic scattering" in text
            and "Backward scattering" in text,
            "ion table lost its energy-frame or channel contract")
    return {"rows": len(rows), "energy_range_eV": [rows[0][0], rows[-1][0]]}


def verify_results(text: str, expected: str) -> dict[str, int]:
    cases: dict[str, list[list[float]]] = {}
    current: str | None = None
    for line in text.splitlines():
        match = re.fullmatch(r"# Case ([1-4])", line.strip())
        if match:
            current = match.group(1)
            require(current not in cases, f"duplicate Case {current}")
            cases[current] = []
            continue
        value = line.strip()
        if current is not None and value and not value.startswith("#"):
            try:
                row = [float(item) for item in value.split()]
            except ValueError as error:
                raise VerificationError(
                    f"invalid Case {current} result row: {value!r}"
                ) from error
            require(len(row) == 7, "Turner result rows must have seven columns")
            cases[current].append(row)
    actual = {case: len(rows) for case, rows in cases.items()}
    wanted = pairs(expected, "case_rows")
    require(actual == wanted,
            f"result case rows differ: expected {wanted}, got {actual}")
    for case, rows in cases.items():
        require(rows[0][0] == 0.0 and abs(rows[-1][0] - 0.067) < 5e-7,
                f"Case {case} does not span the published 6.7 cm gap")
        require(all(row[1] >= 0.0 and row[5] > 0.0 for row in rows),
                f"Case {case} has invalid density/statistical data")
    return actual


def verify(registry_path: Path, artifact_path: Path) -> dict[str, object]:
    registry = parse_registry(registry_path)
    source = registry["source.publisher_supplement"]
    artifact = artifact_path.resolve()
    require(artifact.is_file(), f"artifact is not a regular file: {artifact}")
    require(artifact.name == source["artifact_name"],
            f"artifact must be named {source['artifact_name']!r}")
    expected_bytes = source.getint("artifact_bytes")
    require(artifact.stat().st_size == expected_bytes,
            "publisher archive byte count differs")
    raw = artifact.read_bytes()
    expected_archive_hash = source["sha256"]
    require(SHA256.fullmatch(expected_archive_hash) is not None,
            "registry archive SHA-256 is invalid")
    require(digest(raw) == expected_archive_hash,
            "publisher archive SHA-256 differs")

    expected_names = [item.strip() for item in source["members"].split(",")]
    require(len(expected_names) == len(set(expected_names)),
            "registry contains duplicate member names")
    member_reports: dict[str, object] = {}
    try:
        with zipfile.ZipFile(artifact) as archive:
            infos = archive.infolist()
            names = [info.filename for info in infos]
            require(names == expected_names,
                    f"archive members differ: expected {expected_names}, got {names}")
            for info in infos:
                require(PurePosixPath(info.filename).name == info.filename,
                        f"archive member is not a flat safe filename: {info.filename}")
                require(not (info.flag_bits & 0x1),
                        f"encrypted archive member is unsupported: {info.filename}")
                section_name = f"member.{info.filename}"
                require(section_name in registry,
                        f"registry has no section [{section_name}]")
                section = registry[section_name]
                require(info.file_size == section.getint("bytes"),
                        f"{info.filename} byte count differs")
                data = archive.read(info)
                expected_hash = section["sha256"]
                require(SHA256.fullmatch(expected_hash) is not None,
                        f"{info.filename} registry SHA-256 is invalid")
                require(digest(data) == expected_hash,
                        f"{info.filename} SHA-256 differs")
                text = data.decode("utf-8")
                kind = section["kind"]
                if kind == "lxcat-electron-table":
                    semantics = verify_electron(text, section["process_counts"])
                elif kind == "turner-ion-table":
                    semantics = verify_ion(text, section.getint("rows"))
                elif kind == "turner-results":
                    semantics = verify_results(text, section["case_rows"])
                else:
                    raise VerificationError(f"unsupported member kind: {kind}")
                member_reports[info.filename] = {
                    "bytes": info.file_size,
                    "sha256": expected_hash,
                    "kind": kind,
                    "semantics": semantics,
                }
    except (OSError, zipfile.BadZipFile, UnicodeError) as error:
        raise VerificationError(f"cannot inspect publisher archive: {error}") from error

    return {
        "turner_source_verification_version": 1,
        "case_id": registry["global"]["case_id"],
        "doi": registry["global"]["doi"],
        "artifact_name": artifact.name,
        "artifact_bytes": len(raw),
        "artifact_sha256": expected_archive_hash,
        "license": source["license"],
        "redistribution": source["redistribution"],
        "verified": True,
        "members": member_reports,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify an externally acquired Turner CCP supplement"
    )
    parser.add_argument("registry", type=Path)
    parser.add_argument("artifact", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = verify(args.registry, args.artifact)
        rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.output is None:
            print(rendered, end="")
        else:
            require(not args.output.exists(),
                    f"refusing to overwrite existing report: {args.output}")
            args.output.write_text(rendered, encoding="utf-8")
    except (VerificationError, OSError) as error:
        print(f"Turner source verification error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
