#!/usr/bin/env python3
"""Run and audit a serialized AuroraPIC swarm validation campaign."""

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
import subprocess
import sys
import tempfile
from typing import Iterable


IDENTIFIER = re.compile(r"^[A-Za-z0-9_.:-]+$")
MAX_CAMPAIGN_RUNS = 16
IDENTITY_COLUMNS = (
    "dataset_id",
    "dataset_version",
    "gas",
    "population_model",
    "collision_model_signature",
)


class CampaignInputError(RuntimeError):
    """Raised when a campaign input or generated artifact is invalid."""


@dataclass(frozen=True)
class CampaignRun:
    name: str
    config_file: Path
    result_file: Path


@dataclass(frozen=True)
class Observable:
    name: str
    simulation_column: str
    uncertainty_column: str | None
    relative_tolerance: float
    absolute_tolerance: float
    uncertainty_multiplier: float


@dataclass(frozen=True)
class Campaign:
    path: Path
    campaign_id: str
    campaign_version: str
    provenance: str
    retrieved: str
    reference_manifest: Path | None
    reference_run: str
    field_absolute_tolerance_td: float
    field_relative_tolerance: float
    runs: tuple[CampaignRun, ...]
    observables: tuple[Observable, ...]


def required(
    section: configparser.SectionProxy,
    key: str,
    context: str,
) -> str:
    if key not in section or not section[key].strip():
        raise CampaignInputError(f"{context} requires {key!r}")
    return section[key].strip()


def finite_number(value: str, context: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise CampaignInputError(
            f"{context}: invalid number {value!r}"
        ) from error
    if not math.isfinite(result):
        raise CampaignInputError(f"{context}: number must be finite")
    return result


def nonnegative_number(value: str, context: str) -> float:
    result = finite_number(value, context)
    if result < 0.0:
        raise CampaignInputError(f"{context}: value must be non-negative")
    return result


def resolved_path(base: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def load_campaign(path: Path) -> Campaign:
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
        raise CampaignInputError(
            f"cannot read campaign manifest {path}: {error}"
        ) from error
    if parser.defaults():
        raise CampaignInputError(
            f"{path}: [DEFAULT] values are not supported"
        )
    if "campaign" not in parser:
        raise CampaignInputError(f"{path}: missing [campaign] section")

    section = parser["campaign"]
    allowed = {
        "swarm_campaign_version",
        "campaign_id",
        "campaign_version",
        "provenance",
        "retrieved",
        "reference_manifest",
        "run_order",
        "reference_run",
        "field_absolute_tolerance_td",
        "field_relative_tolerance",
    }
    unknown = sorted(set(section) - allowed)
    if unknown:
        raise CampaignInputError(
            f"{path} [campaign]: unknown keys {unknown}"
        )
    context = f"{path} [campaign]"
    if required(section, "swarm_campaign_version", context) != "1":
        raise CampaignInputError(
            f"{path} supports swarm_campaign_version = 1"
        )
    campaign_id = required(section, "campaign_id", context)
    if not IDENTIFIER.fullmatch(campaign_id):
        raise CampaignInputError(
            f"{context}: campaign_id contains unsupported characters"
        )
    retrieved = required(section, "retrieved", context)
    try:
        date.fromisoformat(retrieved)
    except ValueError as error:
        raise CampaignInputError(
            f"{context}: retrieved must be a valid YYYY-MM-DD date"
        ) from error

    run_names = [
        item.strip()
        for item in required(section, "run_order", context).split(",")
    ]
    if (
        len(run_names) < 2
        or len(run_names) > MAX_CAMPAIGN_RUNS
        or any(not name or not IDENTIFIER.fullmatch(name) for name in run_names)
        or len(run_names) != len(set(run_names))
    ):
        raise CampaignInputError(
            f"{context}: run_order requires 2-{MAX_CAMPAIGN_RUNS} "
            "unique identifiers"
        )
    reference_run = required(section, "reference_run", context)
    if reference_run not in run_names:
        raise CampaignInputError(
            f"{context}: reference_run must occur in run_order"
        )

    runs: list[CampaignRun] = []
    for name in run_names:
        section_name = f"run.{name}"
        if section_name not in parser:
            raise CampaignInputError(
                f"{path}: missing [{section_name}] section"
            )
        run_section = parser[section_name]
        unknown = sorted(set(run_section) - {"config_file", "result_file"})
        if unknown:
            raise CampaignInputError(
                f"{path} [{section_name}]: unknown keys {unknown}"
            )
        run_context = f"{path} [{section_name}]"
        runs.append(
            CampaignRun(
                name=name,
                config_file=resolved_path(
                    path.parent,
                    required(run_section, "config_file", run_context),
                ),
                result_file=resolved_path(
                    path.parent,
                    required(run_section, "result_file", run_context),
                ),
            )
        )
    if len({run.result_file for run in runs}) != len(runs):
        raise CampaignInputError(
            f"{path}: every run must have a unique result_file"
        )

    observable_allowed = {
        "simulation_column",
        "uncertainty_column",
        "relative_tolerance",
        "absolute_tolerance",
        "uncertainty_multiplier",
    }
    observables: list[Observable] = []
    expected_sections = {"campaign", *(f"run.{name}" for name in run_names)}
    for section_name in parser.sections():
        if section_name in expected_sections:
            continue
        if not section_name.startswith("observable."):
            raise CampaignInputError(
                f"{path}: unknown section [{section_name}]"
            )
        name = section_name.removeprefix("observable.").strip()
        if not name or not IDENTIFIER.fullmatch(name):
            raise CampaignInputError(
                f"{path}: invalid observable section [{section_name}]"
            )
        observable_section = parser[section_name]
        unknown = sorted(set(observable_section) - observable_allowed)
        if unknown:
            raise CampaignInputError(
                f"{path} [{section_name}]: unknown keys {unknown}"
            )
        observable_context = f"{path} [{section_name}]"
        uncertainty_column = observable_section.get(
            "uncertainty_column", ""
        ).strip() or None
        observables.append(
            Observable(
                name=name,
                simulation_column=required(
                    observable_section,
                    "simulation_column",
                    observable_context,
                ),
                uncertainty_column=uncertainty_column,
                relative_tolerance=nonnegative_number(
                    observable_section.get("relative_tolerance", "0"),
                    f"{observable_context} relative_tolerance",
                ),
                absolute_tolerance=nonnegative_number(
                    observable_section.get("absolute_tolerance", "0"),
                    f"{observable_context} absolute_tolerance",
                ),
                uncertainty_multiplier=nonnegative_number(
                    observable_section.get("uncertainty_multiplier", "0"),
                    f"{observable_context} uncertainty_multiplier",
                ),
            )
        )
    if not observables:
        raise CampaignInputError(
            f"{path}: at least one [observable.<name>] is required"
        )
    if len({item.name for item in observables}) != len(observables):
        raise CampaignInputError(f"{path}: duplicate observable names")

    return Campaign(
        path=path.resolve(),
        campaign_id=campaign_id,
        campaign_version=required(
            section, "campaign_version", context
        ),
        provenance=required(section, "provenance", context),
        retrieved=retrieved,
        reference_manifest=(
            resolved_path(path.parent, section["reference_manifest"].strip())
            if section.get("reference_manifest", "").strip()
            else None
        ),
        reference_run=reference_run,
        field_absolute_tolerance_td=nonnegative_number(
            section.get("field_absolute_tolerance_td", "1e-12"),
            f"{context} field_absolute_tolerance_td",
        ),
        field_relative_tolerance=nonnegative_number(
            section.get("field_relative_tolerance", "1e-12"),
            f"{context} field_relative_tolerance",
        ),
        runs=tuple(runs),
        observables=tuple(observables),
    )


def sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise CampaignInputError(f"cannot hash {path}: {error}") from error


def reference_data_file(manifest: Path) -> Path:
    parser = configparser.ConfigParser(
        interpolation=None,
        strict=True,
        empty_lines_in_values=False,
    )
    parser.optionxform = str
    try:
        with manifest.open(encoding="utf-8") as stream:
            parser.read_file(stream)
    except (OSError, UnicodeError, configparser.Error) as error:
        raise CampaignInputError(
            f"cannot inspect reference manifest {manifest}: {error}"
        ) from error
    if "reference" not in parser:
        raise CampaignInputError(
            f"{manifest}: missing [reference] section"
        )
    value = required(
        parser["reference"],
        "data_file",
        f"{manifest} [reference]",
    )
    return resolved_path(manifest.parent, value)


def load_csv(
    path: Path,
    run_name: str,
    observables: tuple[Observable, ...],
) -> list[tuple[float, dict[str, str], int]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None:
                raise CampaignInputError(
                    f"run {run_name!r} result CSV has no header"
                )
            fields = reader.fieldnames
            if (
                any(not field or field != field.strip() for field in fields)
                or len(fields) != len(set(fields))
            ):
                raise CampaignInputError(
                    f"run {run_name!r} result CSV has invalid columns"
                )
            required_columns = {
                "reduced_field_td",
                *IDENTITY_COLUMNS,
                *(item.simulation_column for item in observables),
                *(
                    item.uncertainty_column
                    for item in observables
                    if item.uncertainty_column is not None
                ),
            }
            missing = sorted(required_columns - set(fields))
            if missing:
                raise CampaignInputError(
                    f"run {run_name!r} result CSV is missing columns {missing}"
                )
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise CampaignInputError(
            f"cannot read run {run_name!r} result CSV {path}: {error}"
        ) from error
    if not rows:
        raise CampaignInputError(
            f"run {run_name!r} result CSV has no rows"
        )
    indexed: list[tuple[float, dict[str, str], int]] = []
    for line, row in enumerate(rows, start=2):
        if None in row or any(value is None for value in row.values()):
            raise CampaignInputError(
                f"run {run_name!r} result CSV row {line} has wrong field count"
            )
        field = nonnegative_number(
            row["reduced_field_td"],
            f"run {run_name!r} row {line} reduced_field_td",
        )
        if any(existing[0] == field for existing in indexed):
            raise CampaignInputError(
                f"run {run_name!r} has duplicate E/N value {field}"
            )
        for column in IDENTITY_COLUMNS:
            if not row[column]:
                raise CampaignInputError(
                    f"run {run_name!r} row {line} has empty {column}"
                )
        indexed.append((field, row, line))
    return indexed


def common_identity(
    rows: list[tuple[float, dict[str, str], int]],
    run_name: str,
) -> dict[str, str]:
    result: dict[str, str] = {}
    for column in IDENTITY_COLUMNS:
        values = {row[column] for _, row, _ in rows}
        if len(values) != 1:
            raise CampaignInputError(
                f"run {run_name!r} rows do not share one {column}"
            )
        result[column] = next(iter(values))
    return result


def compare_resolution(
    campaign: Campaign,
    run: CampaignRun,
    run_rows: list[tuple[float, dict[str, str], int]],
    reference_rows: list[tuple[float, dict[str, str], int]],
) -> dict[str, object]:
    points: list[dict[str, object]] = []
    passed = True
    matched_fields: set[float] = set()
    for reference_field, reference_row, _ in reference_rows:
        field_tolerance = (
            campaign.field_absolute_tolerance_td
            + campaign.field_relative_tolerance * abs(reference_field)
        )
        matches = [
            item
            for item in run_rows
            if abs(item[0] - reference_field) <= field_tolerance
        ]
        if len(matches) != 1:
            raise CampaignInputError(
                f"reference run E/N {reference_field} Td matches "
                f"{len(matches)} rows in run {run.name!r}"
            )
        run_field, run_row, run_line = matches[0]
        if run_field in matched_fields:
            raise CampaignInputError(
                f"run {run.name!r} E/N {run_field} Td matched twice"
            )
        matched_fields.add(run_field)
        observable_results: list[dict[str, object]] = []
        point_passed = True
        for observable in campaign.observables:
            value = finite_number(
                run_row[observable.simulation_column],
                f"run {run.name!r} row {run_line} "
                f"{observable.simulation_column}",
            )
            reference_value = finite_number(
                reference_row[observable.simulation_column],
                f"reference run {campaign.reference_run!r} "
                f"{observable.simulation_column}",
            )
            uncertainty = (
                0.0
                if observable.uncertainty_column is None
                else nonnegative_number(
                    run_row[observable.uncertainty_column],
                    f"run {run.name!r} row {run_line} "
                    f"{observable.uncertainty_column}",
                )
            )
            reference_uncertainty = (
                0.0
                if observable.uncertainty_column is None
                else nonnegative_number(
                    reference_row[observable.uncertainty_column],
                    f"reference run {campaign.reference_run!r} "
                    f"{observable.uncertainty_column}",
                )
            )
            residual = value - reference_value
            allowed = (
                observable.absolute_tolerance
                + observable.relative_tolerance * abs(reference_value)
                + observable.uncertainty_multiplier
                * math.hypot(uncertainty, reference_uncertainty)
            )
            observable_passed = abs(residual) <= allowed
            point_passed = point_passed and observable_passed
            observable_results.append(
                {
                    "name": observable.name,
                    "value": value,
                    "reference_value": reference_value,
                    "residual": residual,
                    "combined_standard_uncertainty": math.hypot(
                        uncertainty, reference_uncertainty
                    ),
                    "allowed_absolute_residual": allowed,
                    "normalized_residual": (
                        0.0
                        if allowed == 0.0 and residual == 0.0
                        else None
                        if allowed == 0.0
                        else abs(residual) / allowed
                    ),
                    "passed": observable_passed,
                }
            )
        passed = passed and point_passed
        points.append(
            {
                "reference_reduced_field_td": reference_field,
                "run_reduced_field_td": run_field,
                "passed": point_passed,
                "observables": observable_results,
            }
        )
    if len(matched_fields) != len(run_rows):
        extras = sorted(
            field for field, _, _ in run_rows if field not in matched_fields
        )
        raise CampaignInputError(
            f"run {run.name!r} has E/N points absent from reference run: "
            f"{extras}"
        )
    return {
        "run": run.name,
        "reference_run": campaign.reference_run,
        "passed": passed,
        "points": points,
    }


def write_report(
    path: Path,
    report: dict[str, object],
    overwrite: bool,
) -> None:
    if path.exists() and not overwrite:
        raise CampaignInputError(
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
                raise CampaignInputError(
                    f"output already exists: {path}; "
                    "use --overwrite explicitly"
                ) from error
            Path(temporary_name).unlink()
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)


def run_campaign(
    campaign: Campaign,
    executable: Path,
    comparator: Path,
    output: Path,
    overwrite: bool,
) -> dict[str, object]:
    required_files = [(executable, "swarm executable")]
    if campaign.reference_manifest is not None:
        required_files.extend((
            (comparator, "swarm comparator"),
            (campaign.reference_manifest, "reference manifest"),
        ))
    for required_file, label in required_files:
        if not required_file.is_file():
            raise CampaignInputError(f"missing {label}: {required_file}")
    artifact_dir = output.parent / f"{output.stem}.artifacts"
    comparison_paths = (
        {
            run.name: artifact_dir / f"{run.name}.reference.json"
            for run in campaign.runs
        }
        if campaign.reference_manifest is not None
        else {}
    )
    protected_outputs = [
        output,
        *(run.result_file for run in campaign.runs),
        *comparison_paths.values(),
    ]
    if len(set(protected_outputs)) != len(protected_outputs):
        raise CampaignInputError(
            "aggregate report, result files, and per-run reference reports "
            "must all use unique paths"
        )
    reference_data = (
        reference_data_file(campaign.reference_manifest)
        if campaign.reference_manifest is not None
        else None
    )
    if reference_data is not None and not reference_data.is_file():
        raise CampaignInputError(f"missing reference data file: {reference_data}")
    protected_inputs = {
        campaign.path,
        executable,
        *(run.config_file for run in campaign.runs),
    }
    if campaign.reference_manifest is not None:
        protected_inputs.update((
            campaign.reference_manifest,
            reference_data,
            comparator,
        ))
    collisions = sorted(
        str(path) for path in protected_outputs if path in protected_inputs
    )
    if collisions:
        raise CampaignInputError(
            "campaign outputs must not overwrite campaign inputs: "
            + ", ".join(collisions)
        )
    if not overwrite:
        existing = [str(path) for path in protected_outputs if path.exists()]
        if existing:
            raise CampaignInputError(
                "campaign outputs already exist; use --overwrite explicitly: "
                + ", ".join(existing)
            )
    elif output.exists():
        output.unlink()
    for run in campaign.runs:
        if not run.config_file.is_file():
            raise CampaignInputError(
                f"missing config for run {run.name!r}: {run.config_file}"
            )

    if campaign.reference_manifest is not None:
        artifact_dir.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["OMP_DYNAMIC"] = "FALSE"
    environment["OMP_MAX_ACTIVE_LEVELS"] = "1"
    run_reports: list[dict[str, object]] = []
    reference_passed: bool | None = (
        True if campaign.reference_manifest is not None else None
    )
    for run in campaign.runs:
        if overwrite:
            run.result_file.unlink(missing_ok=True)
        try:
            completed = subprocess.run(
                [str(executable), str(run.config_file)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=environment,
            )
        except OSError as error:
            raise CampaignInputError(
                f"cannot execute swarm run {run.name!r}: {error}"
            ) from error
        if completed.returncode != 0:
            raise CampaignInputError(
                f"swarm run {run.name!r} failed with exit "
                f"{completed.returncode}: {completed.stderr.strip()}"
            )
        if not run.result_file.is_file():
            raise CampaignInputError(
                f"swarm run {run.name!r} did not create declared result "
                f"{run.result_file}"
            )
        run_report: dict[str, object] = {
            "name": run.name,
            "config_file": str(run.config_file),
            "config_sha256": sha256(run.config_file),
            "result_file": str(run.result_file),
            "result_sha256": sha256(run.result_file),
            "reference_report": None,
            "reference_report_sha256": None,
            "reference_passed": None,
            "swarm_stdout": completed.stdout,
            "swarm_stderr": completed.stderr,
            "reference_comparator_stdout": None,
            "reference_comparator_stderr": None,
        }
        if campaign.reference_manifest is not None:
            command = [
                sys.executable,
                str(comparator),
                str(run.result_file),
                str(campaign.reference_manifest),
                "--output",
                str(comparison_paths[run.name]),
            ]
            if overwrite:
                command.append("--overwrite")
            comparison = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if comparison.returncode not in (0, 1):
                raise CampaignInputError(
                    f"reference comparison for run {run.name!r} failed with "
                    f"exit {comparison.returncode}: {comparison.stderr.strip()}"
                )
            assert reference_passed is not None
            reference_passed = (
                reference_passed and comparison.returncode == 0
            )
            try:
                reference_report = json.loads(
                    comparison_paths[run.name].read_text(encoding="utf-8")
                )
            except (OSError, UnicodeError, json.JSONDecodeError) as error:
                raise CampaignInputError(
                    f"cannot read reference report for run {run.name!r}: {error}"
                ) from error
            report_passed = reference_report.get("passed")
            if (
                not isinstance(report_passed, bool)
                or report_passed != (comparison.returncode == 0)
            ):
                raise CampaignInputError(
                    f"reference report for run {run.name!r} has an "
                    "inconsistent pass status"
                )
            run_report.update({
                "reference_report": str(comparison_paths[run.name]),
                "reference_report_sha256": sha256(
                    comparison_paths[run.name]
                ),
                "reference_passed": report_passed,
                "reference_comparator_stdout": comparison.stdout,
                "reference_comparator_stderr": comparison.stderr,
            })
        run_reports.append(run_report)

    rows_by_run = {
        run.name: load_csv(run.result_file, run.name, campaign.observables)
        for run in campaign.runs
    }
    identities = {
        name: common_identity(rows, name)
        for name, rows in rows_by_run.items()
    }
    reference_identity = identities[campaign.reference_run]
    for name, identity in identities.items():
        if identity != reference_identity:
            raise CampaignInputError(
                f"run {name!r} physics identity differs from reference run "
                f"{campaign.reference_run!r}: {identity} != "
                f"{reference_identity}"
            )
    convergence = [
        compare_resolution(
            campaign,
            run,
            rows_by_run[run.name],
            rows_by_run[campaign.reference_run],
        )
        for run in campaign.runs
        if run.name != campaign.reference_run
    ]
    convergence_passed = all(item["passed"] for item in convergence)
    report: dict[str, object] = {
        "report_version": 1,
        "passed": convergence_passed and reference_passed is not False,
        "campaign": {
            "manifest": str(campaign.path),
            "manifest_sha256": sha256(campaign.path),
            "campaign_id": campaign.campaign_id,
            "campaign_version": campaign.campaign_version,
            "provenance": campaign.provenance,
            "retrieved": campaign.retrieved,
            "reference_manifest": (
                str(campaign.reference_manifest)
                if campaign.reference_manifest is not None
                else None
            ),
            "reference_manifest_sha256": (
                sha256(campaign.reference_manifest)
                if campaign.reference_manifest is not None
                else None
            ),
            "reference_run": campaign.reference_run,
            "execution_policy": {
                "run_order": "serial",
                "OMP_NUM_THREADS": "1",
                "OMP_DYNAMIC": "FALSE",
                "OMP_MAX_ACTIVE_LEVELS": "1",
            },
        },
        "tooling": {
            "swarm_executable": str(executable),
            "swarm_executable_sha256": sha256(executable),
            "reference_comparator": (
                str(comparator)
                if campaign.reference_manifest is not None
                else None
            ),
            "reference_comparator_sha256": (
                sha256(comparator)
                if campaign.reference_manifest is not None
                else None
            ),
            "python_version": sys.version,
        },
        "physics_identity": reference_identity,
        "reference_validation_passed": reference_passed,
        "external_reference_available": (
            campaign.reference_manifest is not None
        ),
        "convergence_validation_passed": convergence_passed,
        "convergence_acceptance_rule": (
            "abs(run-reference_run) <= absolute_tolerance + "
            "relative_tolerance*abs(reference_run) + "
            "uncertainty_multiplier*combined_standard_uncertainty"
        ),
        "runs": run_reports,
        "convergence": convergence,
        "claim_boundary": (
            "This reference-free campaign tests numerical convergence only; "
            "it does not validate the gas data or physical transport model."
            if campaign.reference_manifest is None
            else "Passing requires both the declared external-reference and "
                 "numerical-convergence criteria."
        ),
    }
    write_report(output, report, overwrite)
    return report


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign_manifest", type=Path)
    parser.add_argument(
        "--swarm-executable", required=True, type=Path
    )
    parser.add_argument(
        "--comparator",
        type=Path,
        default=Path(__file__).resolve().with_name("compare_swarm.py"),
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(list(argv))


def main(argv: Iterable[str] = sys.argv[1:]) -> int:
    args = parse_args(argv)
    try:
        campaign = load_campaign(args.campaign_manifest)
        report = run_campaign(
            campaign,
            args.swarm_executable.resolve(),
            args.comparator.resolve(),
            args.output.resolve(),
            args.overwrite,
        )
    except (CampaignInputError, OSError) as error:
        print(f"swarm campaign failed: {error}", file=sys.stderr)
        return 2
    if report["passed"]:
        print(
            f"swarm campaign passed: {len(campaign.runs)} serialized runs"
        )
        return 0
    print(
        "swarm campaign did not meet its declared criteria",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
