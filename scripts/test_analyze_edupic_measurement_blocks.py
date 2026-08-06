#!/usr/bin/env python3
"""Synthetic regression for eduPIC native measurement-block analysis."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from test_edupic_measurement_stage import checkpoint, fake_source


ROOT = Path(__file__).resolve().parents[1]
COORDINATOR = ROOT / "scripts/advance_edupic_measurement.py"
ANALYZER = ROOT / "scripts/analyze_edupic_measurement_blocks.py"
ADVANCE_ACK = "I_UNDERSTAND_THIS_ADVANCES_BOUNDED_EDUPIC_MEASUREMENT"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
            prefix="aurorapic_edupic_block_analysis_",
            dir=ROOT / "tmp") as tmp:
        work = Path(tmp)
        source = work / "source"
        source.mkdir()
        checkpoint(source / "picdata.bin", 1, 3, 4)
        (source / "conv.dat").write_text("1 3 4\n", encoding="utf-8")
        fake = work / "fake.py"
        fake.write_text(fake_source(True), encoding="utf-8")
        fake.chmod(0o755)
        binary_hash = hashlib.sha256(fake.read_bytes()).hexdigest()
        input_hash = hashlib.sha256((source / "picdata.bin").read_bytes()).hexdigest()
        campaign = work / "campaign"
        advanced = subprocess.run([
            sys.executable, str(COORDINATOR), str(fake), str(source),
            str(campaign), "--expected-binary-sha256", binary_hash,
            "--expected-input-sha256", input_hash,
            "--target-measurement-cycles", "8", "--block-cycles", "2",
            "--qualified-seconds-per-cycle", "0.01",
            "--max-wall-seconds", "30", "--stage-timeout-seconds", "10",
            "--max-stage-initial-particle-steps", "1000000",
            "--max-stages-per-invocation", "4",
            "--acknowledge-cost", ADVANCE_ACK,
        ], text=True, capture_output=True)
        require(advanced.returncode == 0, advanced.stderr)
        density = work / "density.csv"
        eepf = work / "eepf.csv"
        ifed = work / "ifed.csv"
        output = work / "analysis.json"
        analyzed = subprocess.run([
            sys.executable, str(ANALYZER), str(campaign),
            "--minimum-blocks", "4", "--density-csv", str(density),
            "--eepf-csv", str(eepf), "--ifed-csv", str(ifed),
            "--output", str(output), "--require-eligible",
        ], text=True, capture_output=True)
        require(analyzed.returncode == 0, analyzed.stderr)
        report = json.loads(analyzed.stdout)
        require(report["analysis_eligible"] and
                report["campaign"]["block_count"] == 4 and
                report["aggregation_contract"]["exact_duration_weighted"]
                [0] == "density.dat" and
                abs(report["eepf"]["mixture_weighted_normalization"] - 1.0)
                < 2e-5 and
                max(report["eepf"]["block_total_variation_to_mixture"]) == 0.0,
                "block analyzer produced an invalid report")
        with density.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        require(len(rows) == 400 and
                float(rows[200]["electron_density_duration_mean_m3"]) == 1.0,
                "density aggregate is invalid")
        require(all(path.is_file() for path in (density, eepf, ifed, output)),
                "analyzer omitted an output")
    print("eduPIC measurement-block analysis regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
