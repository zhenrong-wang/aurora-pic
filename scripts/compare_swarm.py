#!/usr/bin/env python3
"""Compare an AuroraPIC swarm CSV with traceable reference coefficients."""

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


class ComparisonInputError(RuntimeError):
    """Raised when a comparison input is ambiguous or malformed."""


@dataclass(frozen=True)
class Observable:
    name: str
    simulation_column: str
    reference_column: str
    simulation_uncertainty_column: str | None
    reference_uncertainty_column: str | None
    relative_tolerance: float
    absolute_tolerance: float
    uncertainty_multiplier: float


@dataclass(frozen=True)
class ReferenceContract:
    path: Path
    data_file: Path
    reference_id: str
    reference_version: str
    gas: str
    population_model: str
    coefficient_convention: str
    provenance: str
    citation: str
    retrieved: str
    license: str
    field_absolute_tolerance_td: float
    field_relative_tolerance: float
    observables: tuple[Observable, ...]


def finite_number(value: str, context: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise ComparisonInputError(
            f"{context}: invalid number {value!r}"
        ) from error
    if not math.isfinite(result):
        raise ComparisonInputError(f"{context}: number must be finite")
    return result


def nonnegative_number(value: str, context: str) -> float:
    result = finite_number(value, context)
    if result < 0.0:
        raise ComparisonInputError(f"{context}: value must be non-negative")
    return result


def required(
    section: configparser.SectionProxy,
    key: str,
    context: str,
) -> str:
    if key not in section or not section[key].strip():
        raise ComparisonInputError(f"{context} requires {key!r}")
    return section[key].strip()


def optional_column(
    section: configparser.SectionProxy,
    key: str,
) -> str | None:
    if key not in section:
        return None
    value = section[key].strip()
    return value or None


def load_contract(path: Path) -> ReferenceContract:
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
        raise ComparisonInputError(
            f"cannot read reference manifest {path}: {error}"
        ) from error
    if parser.defaults():
        raise ComparisonInputError(
            f"{path}: [DEFAULT] values are not supported"
        )

    if "reference" not in parser:
        raise ComparisonInputError(
            f"{path}: missing [reference] section"
        )
    reference = parser["reference"]
    reference_allowed = {
        "swarm_reference_version",
        "data_file",
        "reference_id",
        "reference_version",
        "gas",
        "population_model",
        "coefficient_convention",
        "provenance",
        "citation",
        "retrieved",
        "license",
        "field_absolute_tolerance_td",
        "field_relative_tolerance",
    }
    unknown = sorted(set(reference) - reference_allowed)
    if unknown:
        raise ComparisonInputError(
            f"{path} [reference]: unknown keys {unknown}"
        )
    if required(
        reference, "swarm_reference_version", f"{path} [reference]"
    ) != "1":
        raise ComparisonInputError(
            f"{path} supports swarm_reference_version = 1"
        )
    reference_id = required(
        reference, "reference_id", f"{path} [reference]"
    )
    if not IDENTIFIER.fullmatch(reference_id):
        raise ComparisonInputError(
            f"{path} [reference]: reference_id contains "
            "unsupported characters"
        )
    gas = required(reference, "gas", f"{path} [reference]")
    if not IDENTIFIER.fullmatch(gas):
        raise ComparisonInputError(
            f"{path} [reference]: gas contains unsupported characters"
        )
    population_model = required(
        reference, "population_model", f"{path} [reference]"
    )
    coefficient_convention = required(
        reference, "coefficient_convention", f"{path} [reference]"
    )
    for key, value in (
        ("population_model", population_model),
        ("coefficient_convention", coefficient_convention),
    ):
        if not IDENTIFIER.fullmatch(value):
            raise ComparisonInputError(
                f"{path} [reference]: {key} contains unsupported characters"
            )
    retrieved = required(
        reference, "retrieved", f"{path} [reference]"
    )
    try:
        date.fromisoformat(retrieved)
    except ValueError as error:
        raise ComparisonInputError(
            f"{path} [reference]: retrieved must be a valid YYYY-MM-DD date"
        ) from error
    absolute_field_tolerance = nonnegative_number(
        reference.get("field_absolute_tolerance_td", "1e-12"),
        f"{path} [reference] field_absolute_tolerance_td",
    )
    relative_field_tolerance = nonnegative_number(
        reference.get("field_relative_tolerance", "1e-12"),
        f"{path} [reference] field_relative_tolerance",
    )

    observable_allowed = {
        "simulation_column",
        "reference_column",
        "simulation_uncertainty_column",
        "reference_uncertainty_column",
        "relative_tolerance",
        "absolute_tolerance",
        "uncertainty_multiplier",
    }
    observables: list[Observable] = []
    names: set[str] = set()
    for section_name in parser.sections():
        if section_name == "reference":
            continue
        if not section_name.startswith("observable."):
            raise ComparisonInputError(
                f"{path}: unknown section [{section_name}]"
            )
        name = section_name.removeprefix("observable.").strip()
        if not name or not IDENTIFIER.fullmatch(name):
            raise ComparisonInputError(
                f"{path}: invalid observable section [{section_name}]"
            )
        if name in names:
            raise ComparisonInputError(
                f"{path}: duplicate observable name {name!r}"
            )
        names.add(name)
        section = parser[section_name]
        unknown = sorted(set(section) - observable_allowed)
        if unknown:
            raise ComparisonInputError(
                f"{path} [{section_name}]: unknown keys {unknown}"
            )
        context = f"{path} [{section_name}]"
        observables.append(
            Observable(
                name=name,
                simulation_column=required(
                    section, "simulation_column", context
                ),
                reference_column=required(
                    section, "reference_column", context
                ),
                simulation_uncertainty_column=optional_column(
                    section, "simulation_uncertainty_column"
                ),
                reference_uncertainty_column=optional_column(
                    section, "reference_uncertainty_column"
                ),
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
        )
    if not observables:
        raise ComparisonInputError(
            f"{path}: at least one [observable.<name>] is required"
        )

    data_file = Path(
        required(reference, "data_file", f"{path} [reference]")
    )
    if not data_file.is_absolute():
        data_file = (path.parent / data_file).resolve()
    return ReferenceContract(
        path=path.resolve(),
        data_file=data_file,
        reference_id=reference_id,
        reference_version=required(
            reference, "reference_version", f"{path} [reference]"
        ),
        gas=gas,
        population_model=population_model,
        coefficient_convention=coefficient_convention,
        provenance=required(
            reference, "provenance", f"{path} [reference]"
        ),
        citation=required(
            reference, "citation", f"{path} [reference]"
        ),
        retrieved=retrieved,
        license=required(
            reference, "license", f"{path} [reference]"
        ),
        field_absolute_tolerance_td=absolute_field_tolerance,
        field_relative_tolerance=relative_field_tolerance,
        observables=tuple(observables),
    )


def load_csv(path: Path, label: str) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None:
                raise ComparisonInputError(f"{label} CSV has no header")
            fields = [field.strip() for field in reader.fieldnames]
            if fields != reader.fieldnames:
                raise ComparisonInputError(
                    f"{label} CSV column names must not have "
                    "surrounding whitespace"
                )
            if any(not field for field in fields):
                raise ComparisonInputError(
                    f"{label} CSV has an empty column name"
                )
            if len(fields) != len(set(fields)):
                raise ComparisonInputError(
                    f"{label} CSV has duplicate column names"
                )
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise ComparisonInputError(
            f"cannot read {label} CSV {path}: {error}"
        ) from error
    if not rows:
        raise ComparisonInputError(f"{label} CSV has no data rows")
    for row_number, row in enumerate(rows, start=2):
        if None in row or any(value is None for value in row.values()):
            raise ComparisonInputError(
                f"{label} CSV {path}:{row_number} has the wrong field count"
            )
    return fields, rows


def require_columns(
    fields: list[str],
    columns: Iterable[str | None],
    label: str,
) -> None:
    missing = sorted(
        {
            column
            for column in columns
            if column is not None and column not in fields
        }
    )
    if missing:
        raise ComparisonInputError(
            f"{label} CSV is missing columns {missing}"
        )


def indexed_rows(
    rows: list[dict[str, str]],
    label: str,
) -> list[tuple[float, dict[str, str], int]]:
    result: list[tuple[float, dict[str, str], int]] = []
    for row_number, row in enumerate(rows, start=2):
        field = finite_number(
            row["reduced_field_td"],
            f"{label} CSV row {row_number} reduced_field_td",
        )
        if field < 0.0:
            raise ComparisonInputError(
                f"{label} CSV row {row_number}: "
                "reduced_field_td must be non-negative"
            )
        if any(existing[0] == field for existing in result):
            raise ComparisonInputError(
                f"{label} CSV has duplicate E/N value {field}"
            )
        result.append((field, row, row_number))
    return result


def uncertainty(
    row: dict[str, str],
    column: str | None,
    context: str,
) -> float:
    if column is None:
        return 0.0
    return nonnegative_number(row[column], f"{context} {column}")


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise ComparisonInputError(f"cannot hash {path}: {error}") from error


def compare(
    simulation_path: Path,
    contract: ReferenceContract,
) -> dict[str, object]:
    simulation_fields, simulation_rows = load_csv(
        simulation_path, "simulation"
    )
    reference_fields, reference_rows = load_csv(
        contract.data_file, "reference"
    )
    require_columns(
        simulation_fields,
        [
            "reduced_field_td",
            "gas",
            "dataset_id",
            "dataset_version",
            "population_model",
            "collision_model_signature",
            *(
                column
                for observable in contract.observables
                for column in (
                    observable.simulation_column,
                    observable.simulation_uncertainty_column,
                )
            ),
        ],
        "simulation",
    )
    require_columns(
        reference_fields,
        [
            "reduced_field_td",
            *(
                column
                for observable in contract.observables
                for column in (
                    observable.reference_column,
                    observable.reference_uncertainty_column,
                )
            ),
        ],
        "reference",
    )
    simulation_index = indexed_rows(
        simulation_rows, "simulation"
    )
    reference_index = indexed_rows(reference_rows, "reference")
    simulation_gases = {row["gas"] for _, row, _ in simulation_index}
    if simulation_gases != {contract.gas}:
        raise ComparisonInputError(
            "simulation gas values do not match reference gas "
            f"{contract.gas!r}: {sorted(simulation_gases)}"
        )
    simulation_dataset_ids = {
        row["dataset_id"] for _, row, _ in simulation_index
    }
    simulation_dataset_versions = {
        row["dataset_version"] for _, row, _ in simulation_index
    }
    simulation_population_models = {
        row["population_model"] for _, row, _ in simulation_index
    }
    simulation_collision_signatures = {
        row["collision_model_signature"]
        for _, row, _ in simulation_index
    }
    if len(simulation_dataset_ids) != 1 or "" in simulation_dataset_ids:
        raise ComparisonInputError(
            "simulation rows must share one non-empty dataset_id"
        )
    if (
        len(simulation_dataset_versions) != 1
        or "" in simulation_dataset_versions
    ):
        raise ComparisonInputError(
            "simulation rows must share one non-empty dataset_version"
        )
    if simulation_population_models != {contract.population_model}:
        raise ComparisonInputError(
            "simulation population_model values do not match reference "
            f"contract {contract.population_model!r}: "
            f"{sorted(simulation_population_models)}"
        )
    if (
        len(simulation_collision_signatures) != 1
        or "" in simulation_collision_signatures
    ):
        raise ComparisonInputError(
            "simulation rows must share one non-empty "
            "collision_model_signature"
        )

    comparisons: list[dict[str, object]] = []
    matched_simulation_fields: set[float] = set()
    overall_pass = True
    for reference_field, reference_row, reference_line in reference_index:
        tolerance = (
            contract.field_absolute_tolerance_td
            + contract.field_relative_tolerance
            * abs(reference_field)
        )
        matches = [
            item
            for item in simulation_index
            if abs(item[0] - reference_field) <= tolerance
        ]
        if len(matches) != 1:
            raise ComparisonInputError(
                f"reference E/N {reference_field} Td matches "
                f"{len(matches)} simulation rows; expected exactly one"
            )
        simulation_field, simulation_row, simulation_line = matches[0]
        if simulation_field in matched_simulation_fields:
            raise ComparisonInputError(
                f"simulation E/N {simulation_field} Td matched more than once"
            )
        matched_simulation_fields.add(simulation_field)
        observable_results: list[dict[str, object]] = []
        point_pass = True
        for observable in contract.observables:
            simulation_value = finite_number(
                simulation_row[observable.simulation_column],
                f"simulation CSV row {simulation_line} "
                f"{observable.simulation_column}",
            )
            reference_value = finite_number(
                reference_row[observable.reference_column],
                f"reference CSV row {reference_line} "
                f"{observable.reference_column}",
            )
            simulation_uncertainty = uncertainty(
                simulation_row,
                observable.simulation_uncertainty_column,
                f"simulation CSV row {simulation_line}",
            )
            reference_uncertainty = uncertainty(
                reference_row,
                observable.reference_uncertainty_column,
                f"reference CSV row {reference_line}",
            )
            combined_uncertainty = math.hypot(
                simulation_uncertainty,
                reference_uncertainty,
            )
            residual = simulation_value - reference_value
            absolute_residual = abs(residual)
            allowed = (
                observable.absolute_tolerance
                + observable.relative_tolerance
                * abs(reference_value)
                + observable.uncertainty_multiplier
                * combined_uncertainty
            )
            passed = absolute_residual <= allowed
            point_pass = point_pass and passed
            relative_residual = (
                None
                if reference_value == 0.0
                else residual / reference_value
            )
            normalized_residual = (
                (0.0 if absolute_residual == 0.0 else None)
                if allowed == 0.0
                else absolute_residual / allowed
            )
            observable_results.append(
                {
                    "name": observable.name,
                    "simulation_value": simulation_value,
                    "reference_value": reference_value,
                    "residual": residual,
                    "relative_residual": relative_residual,
                    "simulation_standard_uncertainty": (
                        simulation_uncertainty
                    ),
                    "reference_standard_uncertainty": (
                        reference_uncertainty
                    ),
                    "combined_standard_uncertainty": (
                        combined_uncertainty
                    ),
                    "allowed_absolute_residual": allowed,
                    "normalized_residual": normalized_residual,
                    "passed": passed,
                }
            )
        overall_pass = overall_pass and point_pass
        comparisons.append(
            {
                "reference_reduced_field_td": reference_field,
                "simulation_reduced_field_td": simulation_field,
                "passed": point_pass,
                "observables": observable_results,
            }
        )

    extra_fields = sorted(
        field
        for field, _, _ in simulation_index
        if field not in matched_simulation_fields
    )
    return {
        "report_version": 1,
        "passed": overall_pass,
        "simulation": {
            "file": str(simulation_path.resolve()),
            "sha256": sha256(simulation_path),
            "gas": contract.gas,
            "dataset_id": next(iter(simulation_dataset_ids)),
            "dataset_version": next(
                iter(simulation_dataset_versions)
            ),
            "population_model": contract.population_model,
            "collision_model_signature": next(
                iter(simulation_collision_signatures)
            ),
            "rows": len(simulation_rows),
        },
        "reference": {
            "manifest": str(contract.path),
            "manifest_sha256": sha256(contract.path),
            "data_file": str(contract.data_file),
            "data_sha256": sha256(contract.data_file),
            "reference_id": contract.reference_id,
            "reference_version": contract.reference_version,
            "gas": contract.gas,
            "population_model": contract.population_model,
            "coefficient_convention": contract.coefficient_convention,
            "provenance": contract.provenance,
            "citation": contract.citation,
            "retrieved": contract.retrieved,
            "license": contract.license,
            "rows": len(reference_rows),
        },
        "acceptance_rule": (
            "abs(simulation-reference) <= absolute_tolerance + "
            "relative_tolerance*abs(reference) + "
            "uncertainty_multiplier*combined_standard_uncertainty"
        ),
        "matched_reference_points": len(comparisons),
        "extra_simulation_fields_td": extra_fields,
        "comparisons": comparisons,
    }


def write_report(
    path: Path,
    report: dict[str, object],
    overwrite: bool,
) -> None:
    if path.exists() and not overwrite:
        raise ComparisonInputError(
            f"output already exists: {path}; use --overwrite explicitly"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            temporary_name = stream.name
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
        if overwrite:
            os.replace(temporary_name, path)
        else:
            try:
                os.link(temporary_name, path)
            except FileExistsError as error:
                raise ComparisonInputError(
                    f"output already exists: {path}; "
                    "use --overwrite explicitly"
                ) from error
            except OSError as error:
                raise ComparisonInputError(
                    f"cannot create comparison report {path}: {error}"
                ) from error
            Path(temporary_name).unlink()
        temporary_name = None
    except Exception:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
        raise


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("simulation_csv", type=Path)
    parser.add_argument("reference_manifest", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(list(argv))


def main(argv: Iterable[str] = sys.argv[1:]) -> int:
    args = parse_args(argv)
    try:
        contract = load_contract(args.reference_manifest)
        report = compare(args.simulation_csv, contract)
        write_report(args.output, report, args.overwrite)
    except ComparisonInputError as error:
        print(f"swarm comparison failed: {error}", file=sys.stderr)
        return 2
    if report["passed"]:
        print(
            "swarm comparison passed: "
            f"{report['matched_reference_points']} reference points"
        )
        return 0
    failures = sum(
        1
        for point in report["comparisons"]
        for observable in point["observables"]
        if not observable["passed"]
    )
    print(
        f"swarm comparison did not meet acceptance criteria: "
        f"{failures} observable values failed",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
