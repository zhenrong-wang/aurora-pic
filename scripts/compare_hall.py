#!/usr/bin/env python3
"""Compare AuroraPIC resolved 2D diagnostics with a pinned Hall reference."""

from __future__ import annotations

import argparse
import configparser
import csv
from dataclasses import dataclass
from datetime import date
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from typing import Iterable


IDENTIFIER = re.compile(r"^[A-Za-z0-9_.:-]+$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
PROFILE_SOURCES = {"field", "species"}
MODE_METRICS = {
    "frequency_hz",
    "signed_frequency_hz",
    "phase_velocity_m_s",
    "growth_rate_s",
    "rms_amplitude",
}


class HallComparisonInputError(RuntimeError):
    """Raised when a Hall comparison input is malformed or ambiguous."""


@dataclass(frozen=True)
class Acceptance:
    relative_tolerance: float
    absolute_tolerance: float
    uncertainty_multiplier: float


@dataclass(frozen=True)
class ProfileObservable:
    name: str
    simulation_source: str
    simulation_species: str | None
    simulation_column: str
    reference_column: str
    reference_uncertainty_column: str | None
    acceptance: Acceptance


@dataclass(frozen=True)
class ModeObservable:
    name: str
    simulation_quantity: str
    simulation_species: str
    mode: int
    metric: str
    reference_column: str
    reference_uncertainty_column: str | None
    acceptance: Acceptance


@dataclass(frozen=True)
class HallReference:
    path: Path
    case_id: str
    case_variant: str
    case_manifest_sha256: str
    profile_data: Path
    profile_sha256: str
    mode_data: Path
    mode_sha256: str
    profile_axis: str
    mode_axis: str
    coordinate_column: str
    coordinate_tolerance: float
    provenance: str
    citation: str
    retrieved: str
    license: str
    profile_observables: tuple[ProfileObservable, ...]
    mode_observables: tuple[ModeObservable, ...]


def required(
    section: configparser.SectionProxy,
    key: str,
    context: str,
) -> str:
    if key not in section or not section[key].strip():
        raise HallComparisonInputError(f"{context} requires {key!r}")
    return section[key].strip()


def finite_number(value: str, context: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise HallComparisonInputError(
            f"{context}: invalid number {value!r}"
        ) from error
    if not math.isfinite(result):
        raise HallComparisonInputError(f"{context}: number must be finite")
    return result


def nonnegative_number(value: str, context: str) -> float:
    result = finite_number(value, context)
    if result < 0.0:
        raise HallComparisonInputError(
            f"{context}: value must be non-negative"
        )
    return result


def optional_column(
    section: configparser.SectionProxy,
    key: str,
) -> str | None:
    value = section.get(key, "").strip()
    return value or None


def resolved_path(base: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def acceptance(
    section: configparser.SectionProxy,
    context: str,
) -> Acceptance:
    return Acceptance(
        relative_tolerance=nonnegative_number(
            section.get("relative_tolerance", "0"),
            f"{context} relative_tolerance",
        ),
        absolute_tolerance=nonnegative_number(
            section.get("absolute_tolerance", "0"),
            f"{context} absolute_tolerance",
        ),
        uncertainty_multiplier=nonnegative_number(
            section.get("uncertainty_multiplier", "0"),
            f"{context} uncertainty_multiplier",
        ),
    )


def load_reference(path: Path) -> HallReference:
    parser = configparser.ConfigParser(
        interpolation=None,
        strict=True,
        empty_lines_in_values=False,
    )
    parser.optionxform = str
    try:
        with path.open(encoding="utf-8") as stream:
            parser.read_file(stream)
    except (OSError, UnicodeError, configparser.Error) as error:
        raise HallComparisonInputError(
            f"cannot read Hall reference manifest {path}: {error}"
        ) from error
    if parser.defaults():
        raise HallComparisonInputError(
            f"{path}: [DEFAULT] values are not supported"
        )
    if "reference" not in parser:
        raise HallComparisonInputError(
            f"{path}: missing [reference] section"
        )
    reference = parser["reference"]
    allowed = {
        "hall_reference_version",
        "case_id",
        "case_variant",
        "case_manifest_sha256",
        "profile_data_file",
        "profile_data_sha256",
        "mode_data_file",
        "mode_data_sha256",
        "profile_axis",
        "mode_axis",
        "coordinate_column",
        "coordinate_absolute_tolerance",
        "provenance",
        "citation",
        "retrieved",
        "license",
    }
    unknown = sorted(set(reference) - allowed)
    if unknown:
        raise HallComparisonInputError(
            f"{path} [reference]: unknown keys {unknown}"
        )
    context = f"{path} [reference]"
    if required(reference, "hall_reference_version", context) != "1":
        raise HallComparisonInputError(
            f"{path} supports hall_reference_version = 1"
        )
    case_id = required(reference, "case_id", context)
    case_variant = required(reference, "case_variant", context)
    if (
        not IDENTIFIER.fullmatch(case_id)
        or not IDENTIFIER.fullmatch(case_variant)
    ):
        raise HallComparisonInputError(
            f"{context}: case identity contains unsupported characters"
        )
    retrieved = required(reference, "retrieved", context)
    try:
        date.fromisoformat(retrieved)
    except ValueError as error:
        raise HallComparisonInputError(
            f"{context}: retrieved must be YYYY-MM-DD"
        ) from error
    profile_hash = required(
        reference, "profile_data_sha256", context
    )
    mode_hash = required(reference, "mode_data_sha256", context)
    case_hash = required(reference, "case_manifest_sha256", context)
    if (
        not SHA256.fullmatch(profile_hash)
        or not SHA256.fullmatch(mode_hash)
        or not SHA256.fullmatch(case_hash)
    ):
        raise HallComparisonInputError(
            f"{context}: data hashes must be lowercase SHA-256"
        )
    profile_axis = required(reference, "profile_axis", context)
    mode_axis = required(reference, "mode_axis", context)
    if (
        profile_axis not in {"x", "y"}
        or mode_axis not in {"x", "y"}
        or profile_axis == mode_axis
    ):
        raise HallComparisonInputError(
            f"{context}: profile_axis and mode_axis must be distinct x/y"
        )

    common_allowed = {
        "reference_column",
        "reference_uncertainty_column",
        "relative_tolerance",
        "absolute_tolerance",
        "uncertainty_multiplier",
    }
    profile_allowed = common_allowed | {
        "simulation_source",
        "simulation_species",
        "simulation_column",
    }
    mode_allowed = common_allowed | {
        "simulation_quantity",
        "simulation_species",
        "mode",
        "metric",
    }
    profile_observables: list[ProfileObservable] = []
    mode_observables: list[ModeObservable] = []
    names: set[str] = set()
    for section_name in parser.sections():
        if section_name == "reference":
            continue
        if section_name.startswith("profile."):
            name = section_name.removeprefix("profile.").strip()
            section = parser[section_name]
            unknown = sorted(set(section) - profile_allowed)
            if unknown:
                raise HallComparisonInputError(
                    f"{path} [{section_name}]: unknown keys {unknown}"
                )
            if not name or not IDENTIFIER.fullmatch(name) or name in names:
                raise HallComparisonInputError(
                    f"{path}: invalid or duplicate [{section_name}]"
                )
            names.add(name)
            source = required(
                section, "simulation_source", f"{path} [{section_name}]"
            )
            if source not in PROFILE_SOURCES:
                raise HallComparisonInputError(
                    f"{path} [{section_name}]: simulation_source must "
                    "be field or species"
                )
            species = optional_column(section, "simulation_species")
            if source == "species" and species is None:
                raise HallComparisonInputError(
                    f"{path} [{section_name}] requires simulation_species"
                )
            if source == "field" and species is not None:
                raise HallComparisonInputError(
                    f"{path} [{section_name}] field source cannot "
                    "declare simulation_species"
                )
            profile_observables.append(
                ProfileObservable(
                    name=name,
                    simulation_source=source,
                    simulation_species=species,
                    simulation_column=required(
                        section,
                        "simulation_column",
                        f"{path} [{section_name}]",
                    ),
                    reference_column=required(
                        section,
                        "reference_column",
                        f"{path} [{section_name}]",
                    ),
                    reference_uncertainty_column=optional_column(
                        section, "reference_uncertainty_column"
                    ),
                    acceptance=acceptance(
                        section, f"{path} [{section_name}]"
                    ),
                )
            )
        elif section_name.startswith("mode."):
            name = section_name.removeprefix("mode.").strip()
            section = parser[section_name]
            unknown = sorted(set(section) - mode_allowed)
            if unknown:
                raise HallComparisonInputError(
                    f"{path} [{section_name}]: unknown keys {unknown}"
                )
            if not name or not IDENTIFIER.fullmatch(name) or name in names:
                raise HallComparisonInputError(
                    f"{path}: invalid or duplicate [{section_name}]"
                )
            names.add(name)
            metric = required(
                section, "metric", f"{path} [{section_name}]"
            )
            if metric not in MODE_METRICS:
                raise HallComparisonInputError(
                    f"{path} [{section_name}]: unsupported metric {metric!r}"
                )
            mode_text = required(
                section, "mode", f"{path} [{section_name}]"
            )
            try:
                mode = int(mode_text)
            except ValueError as error:
                raise HallComparisonInputError(
                    f"{path} [{section_name}]: mode must be an integer"
                ) from error
            if mode <= 0:
                raise HallComparisonInputError(
                    f"{path} [{section_name}]: mode must be positive"
                )
            mode_observables.append(
                ModeObservable(
                    name=name,
                    simulation_quantity=required(
                        section,
                        "simulation_quantity",
                        f"{path} [{section_name}]",
                    ),
                    simulation_species=section.get(
                        "simulation_species", ""
                    ).strip(),
                    mode=mode,
                    metric=metric,
                    reference_column=required(
                        section,
                        "reference_column",
                        f"{path} [{section_name}]",
                    ),
                    reference_uncertainty_column=optional_column(
                        section, "reference_uncertainty_column"
                    ),
                    acceptance=acceptance(
                        section, f"{path} [{section_name}]"
                    ),
                )
            )
        else:
            raise HallComparisonInputError(
                f"{path}: unknown section [{section_name}]"
            )
    if not profile_observables:
        raise HallComparisonInputError(
            f"{path}: at least one [profile.<name>] is required"
        )
    if not mode_observables:
        raise HallComparisonInputError(
            f"{path}: at least one [mode.<name>] is required"
        )
    return HallReference(
        path=path.resolve(),
        case_id=case_id,
        case_variant=case_variant,
        case_manifest_sha256=case_hash,
        profile_data=resolved_path(
            path.parent, required(reference, "profile_data_file", context)
        ),
        profile_sha256=profile_hash,
        mode_data=resolved_path(
            path.parent, required(reference, "mode_data_file", context)
        ),
        mode_sha256=mode_hash,
        profile_axis=profile_axis,
        mode_axis=mode_axis,
        coordinate_column=required(
            reference, "coordinate_column", context
        ),
        coordinate_tolerance=nonnegative_number(
            reference.get("coordinate_absolute_tolerance", "0"),
            f"{context} coordinate_absolute_tolerance",
        ),
        provenance=required(reference, "provenance", context),
        citation=required(reference, "citation", context),
        retrieved=retrieved,
        license=required(reference, "license", context),
        profile_observables=tuple(profile_observables),
        mode_observables=tuple(mode_observables),
    )


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise HallComparisonInputError(
            f"cannot hash {path}: {error}"
        ) from error


def load_csv(
    path: Path,
    label: str,
) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None:
                raise HallComparisonInputError(
                    f"{label} CSV has no header"
                )
            fields = reader.fieldnames
            if (
                any(not field or field != field.strip() for field in fields)
                or len(fields) != len(set(fields))
            ):
                raise HallComparisonInputError(
                    f"{label} CSV has invalid columns"
                )
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise HallComparisonInputError(
            f"cannot read {label} CSV {path}: {error}"
        ) from error
    if not rows:
        raise HallComparisonInputError(f"{label} CSV has no rows")
    for line, row in enumerate(rows, start=2):
        if None in row or any(value is None for value in row.values()):
            raise HallComparisonInputError(
                f"{label} CSV {path}:{line} has the wrong field count"
            )
    return fields, rows


def require_columns(
    fields: list[str],
    columns: Iterable[str | None],
    label: str,
) -> None:
    missing = sorted(
        column
        for column in columns
        if column is not None and column not in fields
    )
    if missing:
        raise HallComparisonInputError(
            f"{label} CSV is missing columns {missing}"
        )


def load_case_id(path: Path) -> str:
    parser = configparser.ConfigParser(interpolation=None)
    try:
        parser.read_string(
            "[global]\n" + path.read_text(encoding="utf-8")
        )
    except (OSError, UnicodeError, configparser.Error) as error:
        raise HallComparisonInputError(
            f"cannot read Hall case manifest {path}: {error}"
        ) from error
    if "case_id" not in parser["global"]:
        raise HallComparisonInputError(
            f"Hall case manifest {path} has no case_id"
        )
    return parser["global"]["case_id"].strip()


def averaging_identity(
    rows: list[dict[str, str]],
    label: str,
) -> dict[str, float | int | str]:
    required_fields = (
        "start_time",
        "end_time",
        "duration",
        "samples",
        "profile_axis",
    )
    identities: set[tuple[str, ...]] = set()
    for line, row in enumerate(rows, start=2):
        identity = tuple(row[field] for field in required_fields)
        identities.add(identity)
        for field in required_fields[:3]:
            finite_number(row[field], f"{label}:{line} {field}")
        try:
            samples = int(row["samples"])
        except ValueError as error:
            raise HallComparisonInputError(
                f"{label}:{line} samples must be an integer"
            ) from error
        if samples < 2:
            raise HallComparisonInputError(
                f"{label}:{line} requires at least two averaged samples"
            )
    if len(identities) != 1:
        raise HallComparisonInputError(
            f"{label} rows do not share one averaging window"
        )
    identity = next(iter(identities))
    start = finite_number(identity[0], f"{label} start_time")
    end = finite_number(identity[1], f"{label} end_time")
    duration = finite_number(identity[2], f"{label} duration")
    samples = int(identity[3])
    if start < 0.0 or not end > start or not duration > 0.0:
        raise HallComparisonInputError(
            f"{label} has an invalid averaging window"
        )
    if not math.isclose(end - start, duration, rel_tol=1e-12):
        raise HallComparisonInputError(
            f"{label} duration does not equal end_time - start_time"
        )
    return {
        "start_time": start,
        "end_time": end,
        "duration": duration,
        "samples": samples,
        "profile_axis": identity[4],
    }


def unique_coordinate_rows(
    rows: list[dict[str, str]],
    coordinate_column: str,
    label: str,
) -> list[tuple[float, dict[str, str]]]:
    indexed = sorted(
        (
            finite_number(
                row[coordinate_column],
                f"{label}:{line} {coordinate_column}",
            ),
            row,
        )
        for line, row in enumerate(rows, start=2)
    )
    for index in range(1, len(indexed)):
        if indexed[index][0] == indexed[index - 1][0]:
            raise HallComparisonInputError(
                f"{label} has duplicate coordinate {indexed[index][0]}"
            )
    return indexed


def match_coordinate(
    coordinate: float,
    simulation: list[tuple[float, dict[str, str]]],
    tolerance: float,
    label: str,
) -> tuple[float, dict[str, str]]:
    candidates = [
        item for item in simulation
        if abs(item[0] - coordinate) <= tolerance
    ]
    if len(candidates) != 1:
        raise HallComparisonInputError(
            f"{label} coordinate {coordinate} matched "
            f"{len(candidates)} simulation points"
        )
    return candidates[0]


def accepted_result(
    simulation: float,
    reference: float,
    reference_uncertainty: float,
    criteria: Acceptance,
) -> dict[str, float | bool]:
    residual = simulation - reference
    threshold = (
        criteria.absolute_tolerance
        + criteria.relative_tolerance * abs(reference)
        + criteria.uncertainty_multiplier * reference_uncertainty
    )
    return {
        "simulation": simulation,
        "reference": reference,
        "reference_uncertainty": reference_uncertainty,
        "residual": residual,
        "absolute_residual": abs(residual),
        "acceptance_threshold": threshold,
        "passed": abs(residual) <= threshold,
    }


def linear_slope(times: list[float], values: list[float]) -> float:
    mean_time = sum(times) / len(times)
    mean_value = sum(values) / len(values)
    denominator = sum((time - mean_time) ** 2 for time in times)
    if denominator == 0.0:
        raise HallComparisonInputError(
            "mode history has no time extent"
        )
    return sum(
        (time - mean_time) * (value - mean_value)
        for time, value in zip(times, values)
    ) / denominator


def mode_metrics(
    rows: list[dict[str, str]],
    observable: ModeObservable,
) -> dict[str, float]:
    selected: list[tuple[float, float, float, float, float]] = []
    for line, row in enumerate(rows, start=2):
        try:
            row_mode = int(row["mode"])
        except ValueError as error:
            raise HallComparisonInputError(
                f"mode history:{line} mode must be an integer"
            ) from error
        if (
            row_mode != observable.mode
            or row["quantity"] != observable.simulation_quantity
            or row["species"] != observable.simulation_species
        ):
            continue
        selected.append(
            (
                finite_number(row["time"], f"mode history:{line} time"),
                finite_number(row["wavenumber"], f"mode history:{line} wavenumber"),
                finite_number(row["real"], f"mode history:{line} real"),
                finite_number(row["imaginary"], f"mode history:{line} imaginary"),
                finite_number(row["amplitude"], f"mode history:{line} amplitude"),
            )
        )
    selected.sort()
    if len(selected) < 3:
        raise HallComparisonInputError(
            f"mode observable {observable.name!r} requires at least "
            "three history samples"
        )
    times = [item[0] for item in selected]
    if len(times) != len(set(times)):
        raise HallComparisonInputError(
            f"mode observable {observable.name!r} has duplicate times"
        )
    wavenumbers = {item[1] for item in selected}
    if len(wavenumbers) != 1:
        raise HallComparisonInputError(
            f"mode observable {observable.name!r} changed wavenumber"
        )
    amplitudes = [item[4] for item in selected]
    if any(not amplitude > 0.0 for amplitude in amplitudes):
        raise HallComparisonInputError(
            f"mode observable {observable.name!r} has zero amplitude"
        )
    phases = [math.atan2(item[3], item[2]) for item in selected]
    unwrapped = [phases[0]]
    for phase in phases[1:]:
        delta = phase - unwrapped[-1]
        while delta > math.pi:
            phase -= 2.0 * math.pi
            delta -= 2.0 * math.pi
        while delta < -math.pi:
            phase += 2.0 * math.pi
            delta += 2.0 * math.pi
        unwrapped.append(phase)
    phase_slope = linear_slope(times, unwrapped)
    signed_frequency = -phase_slope / (2.0 * math.pi)
    wavenumber = next(iter(wavenumbers))
    return {
        "frequency_hz": abs(signed_frequency),
        "signed_frequency_hz": signed_frequency,
        "phase_velocity_m_s": (
            abs(phase_slope) / abs(wavenumber)
            if wavenumber != 0.0 else math.nan
        ),
        "growth_rate_s": linear_slope(
            times, [math.log(value) for value in amplitudes]
        ),
        "rms_amplitude": math.sqrt(
            sum(value * value for value in amplitudes) /
            len(amplitudes)
        ),
    }


def compare(
    output_dir: Path,
    case_manifest: Path,
    reference: HallReference,
) -> dict[str, object]:
    case_id = load_case_id(case_manifest)
    if case_id != reference.case_id:
        raise HallComparisonInputError(
            f"case manifest id {case_id!r} does not match "
            f"reference case {reference.case_id!r}"
        )
    actual_case_hash = sha256(case_manifest)
    if actual_case_hash != reference.case_manifest_sha256:
        raise HallComparisonInputError(
            "Hall simulation case-manifest SHA-256 mismatch"
        )
    actual_profile_hash = sha256(reference.profile_data)
    actual_mode_hash = sha256(reference.mode_data)
    if actual_profile_hash != reference.profile_sha256:
        raise HallComparisonInputError(
            "Hall reference profile SHA-256 mismatch"
        )
    if actual_mode_hash != reference.mode_sha256:
        raise HallComparisonInputError(
            "Hall reference mode SHA-256 mismatch"
        )

    field_path = output_dir / "resolved_field_time_average.csv"
    species_path = output_dir / "resolved_species_time_average.csv"
    mode_path = output_dir / "resolved_modes.csv"
    field_fields, field_rows = load_csv(field_path, "field time average")
    species_fields, species_rows = load_csv(
        species_path, "species time average"
    )
    mode_fields, mode_rows = load_csv(mode_path, "mode history")
    require_columns(
        field_fields,
        {
            "start_time", "end_time", "duration", "samples",
            "profile_axis", "coordinate",
            *(item.simulation_column for item in reference.profile_observables
              if item.simulation_source == "field"),
        },
        "field time average",
    )
    require_columns(
        species_fields,
        {
            "start_time", "end_time", "duration", "samples",
            "profile_axis", "coordinate", "species",
            *(item.simulation_column for item in reference.profile_observables
              if item.simulation_source == "species"),
        },
        "species time average",
    )
    require_columns(
        mode_fields,
        {
            "time", "mode_axis", "mode", "wavenumber",
            "quantity", "species",
            "real", "imaginary", "amplitude",
        },
        "mode history",
    )
    field_window = averaging_identity(field_rows, "field time average")
    species_window = averaging_identity(
        species_rows, "species time average"
    )
    if field_window != species_window:
        raise HallComparisonInputError(
            "field and species averaging windows do not match"
        )
    if field_window["profile_axis"] != reference.profile_axis:
        raise HallComparisonInputError(
            "simulation profile axis does not match the reference contract"
        )
    if any(row["mode_axis"] != reference.mode_axis for row in mode_rows):
        raise HallComparisonInputError(
            "simulation mode axis does not match the reference contract"
        )

    profile_fields, profile_rows = load_csv(
        reference.profile_data, "Hall reference profile"
    )
    mode_reference_fields, mode_reference_rows = load_csv(
        reference.mode_data, "Hall reference mode"
    )
    require_columns(
        profile_fields,
        {
            reference.coordinate_column,
            *(item.reference_column
              for item in reference.profile_observables),
            *(item.reference_uncertainty_column
              for item in reference.profile_observables
              if item.reference_uncertainty_column is not None),
        },
        "Hall reference profile",
    )
    require_columns(
        mode_reference_fields,
        {
            "mode",
            *(item.reference_column
              for item in reference.mode_observables),
            *(item.reference_uncertainty_column
              for item in reference.mode_observables
              if item.reference_uncertainty_column is not None),
        },
        "Hall reference mode",
    )

    reference_profile = unique_coordinate_rows(
        profile_rows, reference.coordinate_column,
        "Hall reference profile",
    )
    simulation_profiles: dict[
        tuple[str, str | None],
        list[tuple[float, dict[str, str]]],
    ] = {
        ("field", None): unique_coordinate_rows(
            field_rows, "coordinate", "field time average"
        )
    }
    for observable in reference.profile_observables:
        key = (
            observable.simulation_source,
            observable.simulation_species,
        )
        if key in simulation_profiles or key[0] != "species":
            continue
        selected = [
            row for row in species_rows
            if row["species"] == key[1]
        ]
        if not selected:
            raise HallComparisonInputError(
                f"species time average has no species {key[1]!r}"
            )
        simulation_profiles[key] = unique_coordinate_rows(
            selected, "coordinate",
            f"species time average {key[1]}",
        )

    profile_results: list[dict[str, object]] = []
    passed = True
    for coordinate, reference_row in reference_profile:
        row_result: dict[str, object] = {
            "reference_coordinate": coordinate,
            "observables": [],
        }
        for observable in reference.profile_observables:
            key = (
                observable.simulation_source,
                observable.simulation_species,
            )
            matched_coordinate, simulation_row = match_coordinate(
                coordinate,
                simulation_profiles[key],
                reference.coordinate_tolerance,
                f"profile observable {observable.name!r}",
            )
            simulation_value = finite_number(
                simulation_row[observable.simulation_column],
                f"profile {observable.name} simulation",
            )
            reference_value = finite_number(
                reference_row[observable.reference_column],
                f"profile {observable.name} reference",
            )
            uncertainty = (
                nonnegative_number(
                    reference_row[
                        observable.reference_uncertainty_column
                    ],
                    f"profile {observable.name} uncertainty",
                )
                if observable.reference_uncertainty_column
                else 0.0
            )
            result = accepted_result(
                simulation_value,
                reference_value,
                uncertainty,
                observable.acceptance,
            )
            result.update(
                {
                    "name": observable.name,
                    "simulation_coordinate": matched_coordinate,
                }
            )
            passed = passed and bool(result["passed"])
            row_result["observables"].append(result)
        profile_results.append(row_result)

    reference_mode_by_number: dict[int, dict[str, str]] = {}
    for line, row in enumerate(mode_reference_rows, start=2):
        try:
            mode = int(row["mode"])
        except ValueError as error:
            raise HallComparisonInputError(
                f"Hall reference mode:{line} mode must be an integer"
            ) from error
        if mode in reference_mode_by_number:
            raise HallComparisonInputError(
                f"Hall reference mode has duplicate mode {mode}"
            )
        reference_mode_by_number[mode] = row
    mode_results: list[dict[str, object]] = []
    for observable in reference.mode_observables:
        if observable.mode not in reference_mode_by_number:
            raise HallComparisonInputError(
                f"Hall reference mode has no mode {observable.mode}"
            )
        metrics = mode_metrics(mode_rows, observable)
        simulation_value = metrics[observable.metric]
        if not math.isfinite(simulation_value):
            raise HallComparisonInputError(
                f"mode observable {observable.name!r} is not finite"
            )
        reference_row = reference_mode_by_number[observable.mode]
        reference_value = finite_number(
            reference_row[observable.reference_column],
            f"mode {observable.name} reference",
        )
        uncertainty = (
            nonnegative_number(
                reference_row[observable.reference_uncertainty_column],
                f"mode {observable.name} uncertainty",
            )
            if observable.reference_uncertainty_column else 0.0
        )
        result = accepted_result(
            simulation_value,
            reference_value,
            uncertainty,
            observable.acceptance,
        )
        result.update(
            {
                "name": observable.name,
                "mode": observable.mode,
                "metric": observable.metric,
                "simulation_quantity": observable.simulation_quantity,
                "simulation_species": observable.simulation_species,
                "derived_metrics": metrics,
            }
        )
        passed = passed and bool(result["passed"])
        mode_results.append(result)

    return {
        "schema_version": 1,
        "passed": passed,
        "case_id": reference.case_id,
        "case_variant": reference.case_variant,
        "profile_axis": reference.profile_axis,
        "mode_axis": reference.mode_axis,
        "averaging_window": field_window,
        "reference": {
            "manifest": str(reference.path),
            "manifest_sha256": sha256(reference.path),
            "profile_data": str(reference.profile_data),
            "profile_sha256": actual_profile_hash,
            "mode_data": str(reference.mode_data),
            "mode_sha256": actual_mode_hash,
            "provenance": reference.provenance,
            "citation": reference.citation,
            "retrieved": reference.retrieved,
            "license": reference.license,
        },
        "simulation": {
            "output_dir": str(output_dir.resolve()),
            "case_manifest": str(case_manifest.resolve()),
            "case_manifest_sha256": actual_case_hash,
            "field_average_sha256": sha256(field_path),
            "species_average_sha256": sha256(species_path),
            "mode_history_sha256": sha256(mode_path),
        },
        "profile_comparisons": profile_results,
        "mode_comparisons": mode_results,
    }


def write_json_atomic(
    path: Path,
    report: dict[str, object],
    overwrite: bool,
) -> None:
    if path.exists() and not overwrite:
        raise HallComparisonInputError(
            f"comparison report already exists: {path}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare AuroraPIC resolved profiles and modes with a "
            "checksum-pinned Hall reference"
        )
    )
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("reference_manifest", type=Path)
    parser.add_argument(
        "--case-manifest", type=Path, required=True
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        reference = load_reference(args.reference_manifest)
        report = compare(
            args.output_dir, args.case_manifest, reference
        )
        write_json_atomic(args.output, report, args.overwrite)
    except HallComparisonInputError as error:
        print(f"Hall comparison input error: {error}", file=sys.stderr)
        return 2
    if not report["passed"]:
        print(
            "Hall comparison did not meet acceptance criteria",
            file=sys.stderr,
        )
        return 1
    print(f"Hall comparison passed: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
