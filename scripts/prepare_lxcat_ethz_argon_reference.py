#!/usr/bin/env python3
"""Prepare a local ETHZ argon mobility consensus without redistributing data.

The input is an LXCat ETHZ ``Mobility x gas density`` export containing the
3, 4, 6, 8, and 10 kPa series from Haefliger and Franck (2018).  Raw data and
the generated reference remain local under the LXCat terms.  This committed
preparer makes the prospective interpolation, aggregation, uncertainty, and
unit-conversion rules reproducible.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
import re
import statistics


PRESSURES_KPA = (3, 4, 6, 8, 10)
TARGET_FIELDS_TD = (10.0, 20.0, 30.0)
TD_TO_V_M2 = 1.0e-21
PRESSURE_COMMENT = re.compile(r"p\s*=\s*(\d+)\s*kPa\.?$", re.IGNORECASE)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_argon_mobility(path: Path) -> dict[int, list[tuple[float, float]]]:
    result = {pressure: [] for pressure in PRESSURES_KPA}
    database_is_ethz = False
    species_is_argon = False
    process_is_mobility = False
    active_pressure: int | None = None
    in_table = False

    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = raw_line.strip()
        if line.startswith("DATABASE:"):
            database_is_ethz = "ETHZ" in line
        elif line.startswith("SPECIES:"):
            species_is_argon = bool(re.search(r"e\s*/\s*Ar\s*$", line))
            process_is_mobility = False
            active_pressure = None
            in_table = False
        elif line.startswith("PROCESS:"):
            process_is_mobility = "Mobility x gas density" in line
        elif line.startswith("COMMENT:"):
            match = PRESSURE_COMMENT.fullmatch(line.removeprefix("COMMENT:").strip())
            pressure = int(match.group(1)) if match else None
            active_pressure = pressure if pressure in result else None
        elif line.startswith("COLUMNS:"):
            in_table = (
                database_is_ethz
                and species_is_argon
                and process_is_mobility
                and active_pressure is not None
            )
        elif in_table and line and set(line) <= {"x", "X"}:
            in_table = False
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
            assert active_pressure is not None
            result[active_pressure].append((reduced_field_td, mobility_times_density))

    missing = [pressure for pressure, rows in result.items() if not rows]
    if missing:
        raise RuntimeError(f"missing required ETHZ argon pressure series: {missing}")
    for pressure, rows in result.items():
        fields = [row[0] for row in rows]
        if fields != sorted(fields) or len(fields) != len(set(fields)):
            raise RuntimeError(
                f"ETHZ argon {pressure} kPa fields must be unique and increasing"
            )
    return result


def interpolate(rows: list[tuple[float, float]], target: float) -> float:
    for field, value in rows:
        if field == target:
            return value
    for lower, upper in zip(rows, rows[1:]):
        if lower[0] < target < upper[0]:
            fraction = (target - lower[0]) / (upper[0] - lower[0])
            return lower[1] + fraction * (upper[1] - lower[1])
    raise RuntimeError(f"target {target:g} Td is outside a pressure series")


def consensus_rows(
    datasets: dict[int, list[tuple[float, float]]],
) -> list[tuple[float, float, float, float, float]]:
    result = []
    for target in TARGET_FIELDS_TD:
        mobility_values = [
            interpolate(datasets[pressure], target) for pressure in PRESSURES_KPA
        ]
        mean_mobility = statistics.fmean(mobility_values)
        mobility_standard_deviation = statistics.stdev(mobility_values)
        velocity_scale = target * TD_TO_V_M2
        result.append(
            (
                target,
                mean_mobility * velocity_scale,
                mobility_standard_deviation * velocity_scale,
                min(mobility_values) * velocity_scale,
                max(mobility_values) * velocity_scale,
            )
        )
    return result


def write_reference(
    output: Path, rows: list[tuple[float, float, float, float, float]]
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "ethz-argon-five-pressure-consensus.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            [
                "reduced_field_td",
                "drift_velocity_m_s",
                "drift_velocity_standard_uncertainty_m_s",
                "minimum_pressure_series_drift_velocity_m_s",
                "maximum_pressure_series_drift_velocity_m_s",
                "pressure_series_count",
            ]
        )
        for field, velocity, uncertainty, minimum, maximum in rows:
            writer.writerow(
                [
                    f"{field:.17g}",
                    f"{velocity:.17g}",
                    f"{uncertainty:.17g}",
                    f"{minimum:.17g}",
                    f"{maximum:.17g}",
                    len(PRESSURES_KPA),
                ]
            )

    manifest_path = output / "ethz-argon-five-pressure-consensus.swarm-reference"
    manifest_path.write_text(
        """[reference]
swarm_reference_version = 1
data_file = ethz-argon-five-pressure-consensus.csv
reference_id = ethz.haefliger-franck.2018.argon.mobility.five-pressure-consensus
reference_version = lxcat-ethz-updated-2018-03-12-retrieved-2026-08-24
gas = argon
population_model = fixed_population_no_avalanche
coefficient_convention = temporal_bulk_swarm_drift_equivalent_to_flux_without_nonconservative_collisions
provenance = Locally derived from five LXCat ETHZ pressure series by fixed-field linear interpolation and an unweighted cross-pressure mean; raw and derived tables are not redistributed
citation = P. Haefliger and C. M. Franck, Rev. Sci. Instrum. 89, 023114 (2018), doi:10.1063/1.5002762; ETHZ database, www.lxcat.net/ETHZ, retrieved 2026-08-24
retrieved = 2026-08-24
license = Local validation use under the LXCat terms; source data are not redistributed
field_absolute_tolerance_td = 1e-12
field_relative_tolerance = 1e-12

[observable.drift]
simulation_column = electron_drift_velocity_m_s
reference_column = drift_velocity_m_s
simulation_uncertainty_column = mean_velocity_x_standard_error_m_s
reference_uncertainty_column = drift_velocity_standard_uncertainty_m_s
relative_tolerance = 0.10
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
    rows = consensus_rows(datasets)
    if len(rows) != len(TARGET_FIELDS_TD):
        raise RuntimeError("unexpected ETHZ consensus row count")
    write_reference(args.output_directory.resolve(), rows)
    print(f"source_export_sha256={sha256(args.lxcat_export)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
