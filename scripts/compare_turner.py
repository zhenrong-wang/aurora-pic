#!/usr/bin/env python3
"""Compute the published Turner ion-density X-squared statistic."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile


EXPECTED_NODES = {1: 129, 2: 257, 3: 513, 4: 513}
RANGES_95 = {1: (55.0, 303.0), 2: (177.0, 435.0),
             3: (405.0, 693.0), 4: (417.0, 665.0)}
RANGES_99 = {1: (48.0, 405.0), 2: (160.0, 548.0),
             3: (382.0, 798.0), 4: (392.0, 730.0)}
REFERENCE_COLUMNS = (
    "x_m",
    "ion_density_mean_m-3",
    "ion_density_population_stddev_m-3",
)
CANDIDATE_COLUMNS = ("x_m", "ion_density_mean_m-3")


class ComparisonError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ComparisonError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_report(path: Path, report: dict[str, object]) -> None:
    require(not path.exists(),
            f"refusing to overwrite existing report: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
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


def load_columns(path: Path, columns: tuple[str, ...]) -> list[dict[str, float]]:
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            require(reader.fieldnames is not None, f"{path} has no CSV header")
            missing = [name for name in columns if name not in reader.fieldnames]
            require(not missing, f"{path} is missing columns {missing}")
            rows: list[dict[str, float]] = []
            for line, source in enumerate(reader, 2):
                try:
                    row = {name: float(source[name]) for name in columns}
                except (TypeError, ValueError) as error:
                    raise ComparisonError(
                        f"{path}:{line}: invalid numeric value"
                    ) from error
                require(all(math.isfinite(value) for value in row.values()),
                        f"{path}:{line}: non-finite numeric value")
                rows.append(row)
    except OSError as error:
        raise ComparisonError(f"cannot read {path}: {error}") from error
    return rows


def compare(case: int, reference: Path, candidate: Path,
            audit_path: Path) -> dict[str, object]:
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ComparisonError(f"cannot read normalization audit: {error}") from error
    require(audit.get("turner_normalization_version") == 1,
            "normalization audit version must be 1")
    files = audit.get("normalized_files")
    require(isinstance(files, dict) and reference.name in files,
            "reference is not present in the normalization audit")
    reference_lock = files[reference.name]
    require(isinstance(reference_lock, dict)
            and reference_lock.get("sha256") == sha256(reference),
            "reference SHA-256 differs from the normalization audit")
    expected_name = f"turner_case{case}_benchmark.csv"
    require(reference.name == expected_name,
            f"Case {case} statistical comparison requires {expected_name}")

    reference_rows = load_columns(reference, REFERENCE_COLUMNS)
    candidate_rows = load_columns(candidate, CANDIDATE_COLUMNS)
    expected_nodes = EXPECTED_NODES[case]
    require(len(reference_rows) == expected_nodes,
            f"Case {case} reference must have {expected_nodes} nodes")
    require(len(candidate_rows) == expected_nodes,
            f"Case {case} candidate must have {expected_nodes} nodes")

    terms: list[float] = []
    squared_error: list[float] = []
    squared_reference: list[float] = []
    relative_errors: list[float] = []
    for index, (ref, value) in enumerate(zip(reference_rows, candidate_rows)):
        x_ref = ref["x_m"]
        x_value = value["x_m"]
        require(abs(x_ref - x_value) <=
                1e-12 * max(1.0, abs(x_ref), abs(x_value)),
                f"candidate coordinate differs at node {index}")
        mean = ref["ion_density_mean_m-3"]
        sigma = ref["ion_density_population_stddev_m-3"]
        density = value["ion_density_mean_m-3"]
        require(mean > 0.0 and sigma > 0.0 and density >= 0.0,
                f"invalid density/statistical value at node {index}")
        delta = density - mean
        terms.append((delta / sigma) ** 2)
        squared_error.append(delta * delta)
        squared_reference.append(mean * mean)
        relative_errors.append(abs(delta) / mean)

    statistic = math.fsum(terms)
    interval_95 = RANGES_95[case]
    interval_99 = RANGES_99[case]
    return {
        "turner_comparison_version": 1,
        "case": case,
        "reference": {
            "path": str(reference.resolve()),
            "sha256": sha256(reference),
            "normalization_audit_sha256": sha256(audit_path),
        },
        "candidate": {
            "path": str(candidate.resolve()),
            "sha256": sha256(candidate),
        },
        "nodes": expected_nodes,
        "statistic": {
            "name": "Turner ion-density X-squared",
            "formula_variance": "population_standard_deviation_squared",
            "x_squared": statistic,
            "accepted_95_percent": interval_95[0] <= statistic <= interval_95[1],
            "accepted_99_percent": interval_99[0] <= statistic <= interval_99[1],
            "range_95_percent": list(interval_95),
            "range_99_percent": list(interval_99),
        },
        "secondary_metrics": {
            "relative_l2": math.sqrt(
                math.fsum(squared_error) / math.fsum(squared_reference)
            ),
            "maximum_pointwise_relative_error": max(relative_errors),
        },
        "comparison_scope": "published_baseline_ion_density_statistic_only",
        "physics_claim": "none_without_run_contract_and_final_32_period_average",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=int, choices=range(1, 5), required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--normalization-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = compare(
            args.case, args.reference, args.candidate,
            args.normalization_audit,
        )
        write_report(args.output, report)
    except (ComparisonError, OSError) as error:
        print(f"Turner comparison error: {error}", file=sys.stderr)
        return 2
    print(f"Turner Case {args.case} comparison written to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
