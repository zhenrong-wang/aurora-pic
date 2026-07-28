#!/usr/bin/env python3
"""Conservative end-to-end validation for the electron-swarm CLI."""

from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "synthetic_swarm.swarm"
GAS = ROOT / "examples" / "synthetic_swarm.gas"
COMPARATOR = ROOT / "scripts" / "compare_swarm.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main(argv: list[str]) -> int:
    executable = (
        Path(argv[0]).resolve()
        if argv
        else ROOT / "build" / "aurorapic_swarm"
    )
    require(executable.is_file(), f"missing swarm executable: {executable}")
    source = EXAMPLE.read_text(encoding="utf-8")
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_swarm_cli_"
    ) as temporary:
        work = Path(temporary)
        output = work / "swarm.csv"
        config = work / "swarm.swarm"
        source = source.replace(
            "gas_data_file = synthetic_swarm.gas",
            f"gas_data_file = {GAS}",
        ).replace(
            "output_file = synthetic_swarm.csv",
            f"output_file = {output}",
        )
        config.write_text(source, encoding="utf-8")
        subprocess.run([str(executable), str(config)], check=True)
        require(output.is_file(), "swarm CLI did not write its CSV")
        with output.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        require(len(rows) == 3, "swarm CLI did not emit all E/N rows")
        require(
            [float(row["reduced_field_td"]) for row in rows]
            == [1.0, 5.0, 10.0],
            "swarm CLI changed the configured E/N scan",
        )
        require(
            all(
                row["population_model"]
                == "fixed_population_no_avalanche"
                and float(row["neutral_temperature_k"]) == 300.0
                and 200.0
                < float(row["neutral_velocity_stddev_m_s"])
                < 300.0
                and float(row["neutral_speed_limit_sigma"]) == 8.0
                and float(row["maximum_observed_energy_ev"]) < 15.0
                for row in rows
            ),
            "swarm CLI omitted its model or energy-coverage contract",
        )
        reference = work / "reference.csv"
        with reference.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                ["reduced_field_td", "drift_velocity_m_s"]
            )
            for row in rows[:2]:
                writer.writerow(
                    [
                        row["reduced_field_td"],
                        row["electron_drift_velocity_m_s"],
                    ]
                )
        manifest = work / "reference.swarm-reference"
        manifest.write_text(
            "\n".join(
                [
                    "[reference]",
                    "swarm_reference_version = 2",
                    f"data_file = {reference.name}",
                    "reference_id = aurorapic.synthetic.cli",
                    "reference_version = 1",
                    "gas = synthetic_swarm_gas",
                    "population_model = fixed_population_no_avalanche",
                    "coefficient_convention = flux_fixed_population",
                    "provenance = AuroraPIC synthetic CLI output",
                    "citation = AuroraPIC synthetic fixture",
                    "retrieved = 2026-07-28",
                    "license = Synthetic test data",
                    "neutral_temperature_k = 300",
                    "",
                    "[observable.drift]",
                    "simulation_column = electron_drift_velocity_m_s",
                    "reference_column = drift_velocity_m_s",
                    "relative_tolerance = 0",
                    "absolute_tolerance = 0",
                    "uncertainty_multiplier = 0",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        comparison = work / "comparison.json"
        subprocess.run(
            [
                sys.executable,
                str(COMPARATOR),
                str(output),
                str(manifest),
                "--output",
                str(comparison),
            ],
            check=True,
        )
        require(
            comparison.is_file(),
            "swarm CLI output did not satisfy comparator schema",
        )

    print("swarm CLI validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
