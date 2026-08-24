#!/usr/bin/env python3
"""Prepare a local, non-redistributed Dutton argon mobility reference.

The input is an LXCat Dutton ``Mobility x gas density`` export.  Raw data and
the generated reference remain local because the LXCat export terms prohibit
third-party redistribution.  This script is committed so the selection and
unit conversion are transparent and reproducible for an independently
downloaded export.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
import re


SOURCES = ("Brambring 1964.", "Jager et al 1962.", "Wagner 1964.")
TARGET_FIELDS_TD = (50.0, 100.0)
MAXIMUM_RELATIVE_FIELD_DISTANCE = 0.10
TD_TO_V_M2 = 1.0e-21


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_argon_mobility(path: Path) -> dict[str, list[tuple[float, float]]]:
    result = {source: [] for source in SOURCES}
    species_is_argon = False
    process_is_mobility = False
    active_source: str | None = None
    in_table = False

    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = raw_line.strip()
        if line.startswith("SPECIES:"):
            species_is_argon = bool(re.search(r"e\s*/\s*Ar\s*$", line))
            process_is_mobility = False
            active_source = None
            in_table = False
        elif line.startswith("PROCESS:"):
            process_is_mobility = (
                "Mobility x gas density" in line or "Mobility x gas density" in raw_line
            )
        elif line.startswith("COMMENT:"):
            comment = line.removeprefix("COMMENT:").strip()
            active_source = comment if comment in result else None
        elif line.startswith("COLUMNS:"):
            in_table = species_is_argon and process_is_mobility and active_source is not None
        elif in_table and line and not set(line) <= {"-"}:
            fields = line.split()
            if len(fields) != 2:
                raise RuntimeError(f"{path}:{line_number}: malformed data row")
            try:
                reduced_field_td, mobility_times_density = map(float, fields)
            except ValueError as error:
                raise RuntimeError(
                    f"{path}:{line_number}: malformed numeric data row"
                ) from error
            if not all(
                math.isfinite(value) and value > 0.0
                for value in (reduced_field_td, mobility_times_density)
            ):
                raise RuntimeError(f"{path}:{line_number}: invalid data value")
            assert active_source is not None
            result[active_source].append((reduced_field_td, mobility_times_density))

    missing = [source for source, rows in result.items() if not rows]
    if missing:
        raise RuntimeError(f"missing required argon mobility sources: {missing}")
    return result


def select_points(
    datasets: dict[str, list[tuple[float, float]]],
) -> list[tuple[float, float, str, float]]:
    selected: list[tuple[float, float, str, float]] = []
    for target in TARGET_FIELDS_TD:
        for source in SOURCES:
            field, mobility_n = min(
                datasets[source], key=lambda row: (abs(row[0] - target), row[0])
            )
            relative_distance = abs(field - target) / target
            if relative_distance <= MAXIMUM_RELATIVE_FIELD_DISTANCE:
                selected.append((field, mobility_n, source.rstrip("."), target))
    fields = [row[0] for row in selected]
    if len(fields) != len(set(fields)):
        raise RuntimeError("selection produced duplicate reduced fields")
    return sorted(selected)


def write_reference(output: Path, rows: list[tuple[float, float, str, float]]) -> None:
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "dutton-argon-selected-reference.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            [
                "reduced_field_td",
                "drift_velocity_m_s",
                "source",
                "selection_target_td",
            ]
        )
        for field, mobility_n, source, target in rows:
            writer.writerow(
                [
                    f"{field:.17g}",
                    f"{field * TD_TO_V_M2 * mobility_n:.17g}",
                    source,
                    f"{target:.17g}",
                ]
            )

    manifest_path = output / "dutton-argon-selected.swarm-reference"
    manifest_path.write_text(
        """[reference]
swarm_reference_version = 1
data_file = dutton-argon-selected-reference.csv
reference_id = dutton.1975.argon.mobility.selected-high-field
reference_version = lxcat-dutton-retrieved-2026-08-24
gas = argon
population_model = fixed_population_no_avalanche
coefficient_convention = historical_measured_mobility_unspecified_bulk_flux
provenance = Locally derived from the LXCat Dutton argon mobility export using the committed nearest-source selection rule; raw and derived tables are not redistributed
citation = Jack Dutton, Survey of Electron Swarm Data, J. Phys. Chem. Ref. Data 4, 577 (1975), doi:10.1063/1.555525; original source names retained in local CSV
retrieved = 2026-08-24
license = Local validation use under the LXCat terms; source data are not redistributed
field_absolute_tolerance_td = 1e-12
field_relative_tolerance = 1e-12

[observable.drift]
simulation_column = electron_drift_velocity_m_s
reference_column = drift_velocity_m_s
simulation_uncertainty_column = mean_velocity_x_standard_error_m_s
relative_tolerance = 0.20
absolute_tolerance = 0
uncertainty_multiplier = 2
""",
        encoding="utf-8",
        newline="\n",
    )
    print(f"reference_csv={csv_path}")
    print(f"reference_csv_sha256={sha256(csv_path)}")
    print(f"reference_manifest={manifest_path}")
    print(f"reference_manifest_sha256={sha256(manifest_path)}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lxcat_export", type=Path)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    if not args.lxcat_export.is_file():
        parser.error(f"LXCat export does not exist: {args.lxcat_export}")
    datasets = parse_argon_mobility(args.lxcat_export)
    selected = select_points(datasets)
    if len(selected) != 5:
        raise RuntimeError(f"expected five selected observations, got {len(selected)}")
    write_reference(args.output_directory.resolve(), selected)
    print(f"source_export_sha256={sha256(args.lxcat_export)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
