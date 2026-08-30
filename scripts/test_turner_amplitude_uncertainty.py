#!/usr/bin/env python3
"""Bounded regression for the Turner amplitude uncertainty diagnostic."""

from __future__ import annotations

import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
ANALYZER = ROOT / "scripts" / "analyze_turner_amplitude_uncertainty.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fixture(path: Path, values: list[float]) -> None:
    path.write_text(json.dumps({
        "turner_density_block_analysis_version": 3,
        "case": 1,
        "species": "ions",
        "blocks": [
            {"line_integrated_density_m-2": value}
            for value in values
        ],
        "published_acceptance_applicable": False,
    }), encoding="utf-8")


def run(source: Path, output: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run([
        sys.executable, str(ANALYZER), str(source),
        "--replicates", "4000", "--random-seed", "123456",
        "--output", str(output),
    ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_uncertainty_", dir=ROOT / "tmp"
    ) as temporary:
        work = Path(temporary)
        generator = random.Random(98765)
        state = 0.0
        stationary: list[float] = []
        for _ in range(32):
            state = 0.65 * state + generator.gauss(0.0, 0.006)
            stationary.append(5.0e12 * (1.0 + state))
        stationary_source = work / "stationary.json"
        stationary_output = work / "stationary-output.json"
        fixture(stationary_source, stationary)
        stationary_run = run(stationary_source, stationary_output)
        stationary_value = json.loads(
            stationary_output.read_text(encoding="utf-8")
        )
        require(
            stationary_run.returncode == 0
            and stationary_value["classification"]
                == "exploratory_stationary_ar1_null_not_rejected"
            and stationary_value["parametric_stationary_null"][
                "stationary_null_rejected"
            ] is False
            and stationary_value["published_acceptance_applicable"] is False,
            "stationary AR(1) fixture was rejected or overstated",
        )

        drifting = [
            5.0e12 * (1.0 + 0.008 * index + 0.001 * (-1) ** index)
            for index in range(32)
        ]
        drifting_source = work / "drifting.json"
        drifting_output = work / "drifting-output.json"
        fixture(drifting_source, drifting)
        drifting_run = run(drifting_source, drifting_output)
        drifting_value = json.loads(
            drifting_output.read_text(encoding="utf-8")
        )
        require(
            drifting_run.returncode == 0
            and drifting_value["classification"]
                == "exploratory_stationary_ar1_null_rejected"
            and drifting_value["parametric_stationary_null"][
                "stationary_null_rejected"
            ] is True,
            "strong drifting fixture did not reject the stationary null",
        )

        repeated = work / "repeated.json"
        repeated_run = run(stationary_source, repeated)
        require(
            repeated_run.returncode == 0
            and repeated.read_bytes() == stationary_output.read_bytes(),
            "fixed-seed uncertainty report is not deterministic",
        )

        short_source = work / "short.json"
        fixture(short_source, stationary[:8])
        rejected = run(short_source, work / "short-output.json")
        require(
            rejected.returncode == 2
            and "at least 16" in rejected.stderr,
            "uncertainty diagnostic accepted a short series",
        )

    print("Turner amplitude uncertainty analysis passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
