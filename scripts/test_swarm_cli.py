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
                and float(row["maximum_observed_energy_ev"]) < 15.0
                for row in rows
            ),
            "swarm CLI omitted its model or energy-coverage contract",
        )

    print("swarm CLI validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
