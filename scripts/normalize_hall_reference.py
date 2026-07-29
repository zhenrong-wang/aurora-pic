#!/usr/bin/env python3
"""Normalize checksum-pinned Hall tables into the strict comparison contract."""

from __future__ import annotations

import argparse
import configparser
import csv
from datetime import date
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile


SHA256 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9_.:-]+$")


class NormalizeError(RuntimeError):
    pass


def required(
    section: configparser.SectionProxy,
    key: str,
) -> str:
    value = section.get(key, "").strip()
    if not value:
        raise NormalizeError(f"[{section.name}] requires {key!r}")
    return value


def finite(value: str, context: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise NormalizeError(f"{context}: invalid number {value!r}") from error
    if not math.isfinite(result):
        raise NormalizeError(f"{context}: value must be finite")
    return result


def positive(value: str, context: str) -> float:
    result = finite(value, context)
    if result <= 0.0:
        raise NormalizeError(f"{context}: value must be positive")
    return result


def columns(section: configparser.SectionProxy, key: str) -> tuple[str, ...]:
    result = tuple(
        item.strip() for item in required(section, key).split(",")
        if item.strip()
    )
    if not result or len(result) != len(set(result)):
        raise NormalizeError(f"[{section.name}] {key} has invalid columns")
    return result


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise NormalizeError(f"cannot hash {path}: {error}") from error


def load_lock(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(
        interpolation=None,
        strict=True,
        empty_lines_in_values=False,
    )
    try:
        parser.read(path, encoding="utf-8")
    except (OSError, UnicodeError, configparser.Error) as error:
        raise NormalizeError(f"cannot read source lock {path}: {error}") from error
    if parser.defaults():
        raise NormalizeError("[DEFAULT] values are not supported")
    for section in ("source", "profile", "mode", "acceptance"):
        if section not in parser:
            raise NormalizeError(f"source lock is missing [{section}]")
    unexpected_sections = sorted(
        set(parser.sections()) - {"source", "profile", "mode", "acceptance"}
    )
    if unexpected_sections:
        raise NormalizeError(
            f"source lock has unknown sections {unexpected_sections}"
        )
    source = parser["source"]
    allowed = {
        "source": {
            "hall_source_version", "case_id", "case_variant",
            "case_manifest_sha256", "profile_file", "profile_sha256",
            "mode_file", "mode_sha256", "source_url",
            "source_artifact_id", "provenance", "citation", "retrieved",
            "license",
        },
        "profile": {
            "coordinate_column", "coordinate_scale_to_m",
            "electric_field_columns", "electric_field_scale_to_v_m",
            "ion_density_columns", "ion_density_scale_to_m3",
            "electron_temperature_columns",
            "electron_temperature_scale_to_ev",
        },
        "mode": {
            "mode_column", "frequency_columns", "frequency_scale_to_hz",
            "comparison_mode",
        },
        "acceptance": {
            "coordinate_absolute_tolerance_m", "relative_tolerance",
            "uncertainty_multiplier",
        },
    }
    for section_name, permitted in allowed.items():
        unknown = sorted(set(parser[section_name]) - permitted)
        if unknown:
            raise NormalizeError(
                f"[{section_name}] contains unknown keys {unknown}"
            )
    if required(source, "hall_source_version") != "1":
        raise NormalizeError("only hall_source_version = 1 is supported")
    for key in ("case_id", "case_variant"):
        if not IDENTIFIER.fullmatch(required(source, key)):
            raise NormalizeError(f"[source] {key} is not a safe identifier")
    for key in (
        "source_url", "source_artifact_id", "provenance",
        "citation", "license",
    ):
        safe_manifest_value(required(source, key), f"[source] {key}")
    try:
        date.fromisoformat(required(source, "retrieved"))
    except ValueError as error:
        raise NormalizeError("[source] retrieved must be YYYY-MM-DD") from error
    for key in (
        "case_manifest_sha256", "profile_sha256", "mode_sha256"
    ):
        if not SHA256.fullmatch(required(source, key)):
            raise NormalizeError(f"[source] {key} must be lowercase SHA-256")
    return parser


def load_csv(
    path: Path,
    required_columns: set[str],
    label: str,
) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            fields = reader.fieldnames or []
            if (
                not fields
                or len(fields) != len(set(fields))
                or any(not item or item != item.strip() for item in fields)
            ):
                raise NormalizeError(f"{label} has an invalid header")
            missing = sorted(required_columns - set(fields))
            if missing:
                raise NormalizeError(f"{label} is missing columns {missing}")
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise NormalizeError(f"cannot read {label} {path}: {error}") from error
    if not rows:
        raise NormalizeError(f"{label} has no data rows")
    for line, row in enumerate(rows, start=2):
        if None in row or any(value is None for value in row.values()):
            raise NormalizeError(f"{label}:{line} has the wrong field count")
    return rows


def envelope(
    row: dict[str, str],
    selected: tuple[str, ...],
    scale: float,
    context: str,
) -> tuple[float, float]:
    values = [
        finite(row[column], f"{context} {column}") * scale
        for column in selected
    ]
    lower = min(values)
    upper = max(values)
    return 0.5 * (lower + upper), 0.5 * (upper - lower)


def normalized_profiles(
    rows: list[dict[str, str]],
    section: configparser.SectionProxy,
) -> list[list[object]]:
    coordinate_column = required(section, "coordinate_column")
    coordinate_scale = positive(
        required(section, "coordinate_scale_to_m"),
        "[profile] coordinate_scale_to_m",
    )
    quantities = (
        (
            columns(section, "electric_field_columns"),
            positive(required(section, "electric_field_scale_to_v_m"),
                     "[profile] electric_field_scale_to_v_m"),
        ),
        (
            columns(section, "ion_density_columns"),
            positive(required(section, "ion_density_scale_to_m3"),
                     "[profile] ion_density_scale_to_m3"),
        ),
        (
            columns(section, "electron_temperature_columns"),
            positive(required(section, "electron_temperature_scale_to_ev"),
                     "[profile] electron_temperature_scale_to_ev"),
        ),
    )
    result: list[list[object]] = []
    previous = -math.inf
    for line, row in enumerate(rows, start=2):
        coordinate = finite(
            row[coordinate_column], f"profile:{line} coordinate"
        ) * coordinate_scale
        if coordinate <= previous:
            raise NormalizeError(
                "profile coordinates must increase strictly in source order"
            )
        previous = coordinate
        values: list[object] = [coordinate]
        for selected, scale in quantities:
            center, uncertainty = envelope(
                row, selected, scale, f"profile:{line}"
            )
            values.extend((center, uncertainty))
        if values[3] < values[4] or values[5] < values[6]:
            raise NormalizeError(
                "profile density and temperature envelopes must be non-negative"
            )
        result.append(values)
    return result


def normalized_modes(
    rows: list[dict[str, str]],
    section: configparser.SectionProxy,
) -> list[list[object]]:
    mode_column = required(section, "mode_column")
    frequency_columns = columns(section, "frequency_columns")
    frequency_scale = positive(
        required(section, "frequency_scale_to_hz"),
        "[mode] frequency_scale_to_hz",
    )
    result: list[list[object]] = []
    observed: set[int] = set()
    for line, row in enumerate(rows, start=2):
        numeric_mode = finite(row[mode_column], f"mode:{line} mode")
        mode = int(numeric_mode)
        if numeric_mode != mode or mode <= 0 or mode in observed:
            raise NormalizeError(
                f"mode:{line} must contain a unique positive integer mode"
            )
        observed.add(mode)
        center, uncertainty = envelope(
            row, frequency_columns, frequency_scale, f"mode:{line}"
        )
        if center - uncertainty < 0.0:
            raise NormalizeError("mode frequency must be non-negative")
        result.append([mode, center, uncertainty])
    result.sort(key=lambda item: int(item[0]))
    return result


def csv_text(header: list[str], rows: list[list[object]]) -> str:
    import io
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return stream.getvalue()


def safe_manifest_value(value: str, context: str) -> str:
    if (
        not value
        or value != value.strip()
        or any(item in value for item in ("\n", "\r", "#", ";"))
    ):
        raise NormalizeError(f"{context} is unsafe for an INI manifest")
    return value


def normalize(args: argparse.Namespace) -> Path:
    lock_path = args.source_lock.resolve()
    lock = load_lock(lock_path)
    source = lock["source"]
    profile = lock["profile"]
    mode = lock["mode"]
    acceptance = lock["acceptance"]
    base = lock_path.parent
    profile_path = (base / required(source, "profile_file")).resolve()
    mode_path = (base / required(source, "mode_file")).resolve()
    case_path = args.case_manifest.resolve()
    for path, expected, label in (
        (profile_path, required(source, "profile_sha256"), "profile"),
        (mode_path, required(source, "mode_sha256"), "mode"),
        (
            case_path,
            required(source, "case_manifest_sha256"),
            "case manifest",
        ),
    ):
        if sha256(path) != expected:
            raise NormalizeError(f"{label} SHA-256 mismatch")
    case = configparser.ConfigParser(interpolation=None, strict=True)
    try:
        case.read_string(
            "[global]\n" + case_path.read_text(encoding="utf-8")
        )
    except (OSError, UnicodeError, configparser.Error) as error:
        raise NormalizeError(
            f"cannot read Hall case manifest {case_path}: {error}"
        ) from error
    if case["global"].get("case_id", "").strip() != source["case_id"]:
        raise NormalizeError("source lock case_id does not match case manifest")
    if "reference" not in case:
        raise NormalizeError("Hall case manifest is missing [reference]")

    profile_columns = {
        required(profile, "coordinate_column"),
        *columns(profile, "electric_field_columns"),
        *columns(profile, "ion_density_columns"),
        *columns(profile, "electron_temperature_columns"),
    }
    mode_columns = {
        required(mode, "mode_column"),
        *columns(mode, "frequency_columns"),
    }
    profile_rows = normalized_profiles(
        load_csv(profile_path, profile_columns, "raw profile"),
        profile,
    )
    mode_rows = normalized_modes(
        load_csv(mode_path, mode_columns, "raw mode"),
        mode,
    )
    comparison_mode = int(required(mode, "comparison_mode"))
    if comparison_mode not in {int(row[0]) for row in mode_rows}:
        raise NormalizeError("[mode] comparison_mode is absent from raw data")
    domain_x = positive(
        required(case["reference"], "domain_x_m"),
        "case [reference] domain_x_m",
    )
    if (
        float(profile_rows[0][0]) < 0.0
        or float(profile_rows[-1][0]) > domain_x
    ):
        raise NormalizeError("profile coordinates lie outside the axial domain")
    nyquist = int(required(case["reference"], "production_cells_y")) // 2
    if comparison_mode > nyquist:
        raise NormalizeError("comparison_mode exceeds the production Nyquist limit")

    output = args.output_dir.resolve()
    if output.exists():
        raise NormalizeError(f"refusing to overwrite existing output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent
    ))
    try:
        profile_output = temporary / "profiles.csv"
        mode_output = temporary / "modes.csv"
        profile_output.write_text(
            csv_text(
                [
                    "coordinate_m",
                    "electric_x_v_m",
                    "electric_x_uncertainty_v_m",
                    "ion_density_m3",
                    "ion_density_uncertainty_m3",
                    "electron_temperature_ev",
                    "electron_temperature_uncertainty_ev",
                ],
                profile_rows,
            ),
            encoding="utf-8",
        )
        mode_output.write_text(
            csv_text(
                [
                    "mode", "frequency_hz",
                    "frequency_uncertainty_hz",
                ],
                mode_rows,
            ),
            encoding="utf-8",
        )
        relative_tolerance = finite(
            required(acceptance, "relative_tolerance"),
            "[acceptance] relative_tolerance",
        )
        uncertainty_multiplier = finite(
            required(acceptance, "uncertainty_multiplier"),
            "[acceptance] uncertainty_multiplier",
        )
        if relative_tolerance < 0.0 or uncertainty_multiplier < 0.0:
            raise NormalizeError("acceptance values must be non-negative")
        coordinate_tolerance = finite(
            acceptance.get("coordinate_absolute_tolerance_m", "1e-12"),
            "[acceptance] coordinate_absolute_tolerance_m",
        )
        if coordinate_tolerance < 0.0:
            raise NormalizeError(
                "coordinate_absolute_tolerance_m must be non-negative"
            )
        metadata = {
            key: safe_manifest_value(required(source, key), f"[source] {key}")
            for key in (
                "source_url", "source_artifact_id", "provenance",
                "citation", "license",
            )
        }
        manifest = temporary / "reference.hall-reference"
        manifest.write_text(
            f"""[reference]
hall_reference_version = 1
case_id = {source['case_id']}
case_variant = {source['case_variant']}
case_manifest_sha256 = {source['case_manifest_sha256']}
profile_data_file = profiles.csv
profile_data_sha256 = {sha256(profile_output)}
mode_data_file = modes.csv
mode_data_sha256 = {sha256(mode_output)}
profile_axis = x
mode_axis = y
coordinate_column = coordinate_m
coordinate_absolute_tolerance = {coordinate_tolerance}
provenance = {metadata['provenance']}, artifact {metadata['source_artifact_id']} from {metadata['source_url']}, raw lock SHA-256 {sha256(lock_path)}
citation = {metadata['citation']}
retrieved = {source['retrieved']}
license = {metadata['license']}

[profile.axial_field]
simulation_source = field
simulation_column = electric_x
reference_column = electric_x_v_m
reference_uncertainty_column = electric_x_uncertainty_v_m
relative_tolerance = {relative_tolerance}
uncertainty_multiplier = {uncertainty_multiplier}

[profile.ion_density]
simulation_source = species
simulation_species = ions
simulation_column = number_density
reference_column = ion_density_m3
reference_uncertainty_column = ion_density_uncertainty_m3
relative_tolerance = {relative_tolerance}
uncertainty_multiplier = {uncertainty_multiplier}

[profile.electron_temperature]
simulation_source = species
simulation_species = electrons
simulation_column = temperature_ev
reference_column = electron_temperature_ev
reference_uncertainty_column = electron_temperature_uncertainty_ev
relative_tolerance = {relative_tolerance}
uncertainty_multiplier = {uncertainty_multiplier}

[mode.dominant_frequency]
simulation_quantity = electric_y
mode = {comparison_mode}
metric = frequency_hz
reference_column = frequency_hz
reference_uncertainty_column = frequency_uncertainty_hz
relative_tolerance = {relative_tolerance}
uncertainty_multiplier = {uncertainty_multiplier}
""",
            encoding="utf-8",
        )
        audit = temporary / "normalization.json"
        audit.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "case_id": source["case_id"],
                    "case_variant": source["case_variant"],
                    "source_lock": str(lock_path),
                    "source_lock_sha256": sha256(lock_path),
                    "source_url": source["source_url"],
                    "source_artifact_id": source["source_artifact_id"],
                    "raw_profile_sha256": source["profile_sha256"],
                    "raw_mode_sha256": source["mode_sha256"],
                    "profile_rows": len(profile_rows),
                    "mode_rows": len(mode_rows),
                    "envelope_method": "midpoint_and_half_range",
                    "coordinate_policy": "native_no_interpolation",
                    "profile_sha256": sha256(profile_output),
                    "mode_sha256": sha256(mode_output),
                    "reference_manifest_sha256": sha256(manifest),
                },
                indent=2,
                sort_keys=True,
            ) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, output)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Normalize checksum-pinned native Hall profiles and modes; "
            "never interpolate, extrapolate, download, or run a simulation"
        )
    )
    parser.add_argument("source_lock", type=Path)
    parser.add_argument("--case-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    try:
        output = normalize(parse_args())
    except (NormalizeError, ValueError, configparser.Error) as error:
        print(f"Hall reference normalization error: {error}", file=sys.stderr)
        return 2
    print(f"Normalized checksum-pinned Hall reference: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
