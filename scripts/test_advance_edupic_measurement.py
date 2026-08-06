#!/usr/bin/env python3
"""Synthetic resume and integrity regression for eduPIC measurement campaigns."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from test_edupic_measurement_stage import checkpoint, fake_source


ROOT = Path(__file__).resolve().parents[1]
COORDINATOR = ROOT / "scripts/advance_edupic_measurement.py"
ACK = "I_UNDERSTAND_THIS_ADVANCES_BOUNDED_EDUPIC_MEASUREMENT"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
            prefix="aurorapic_edupic_measurement_advance_",
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
        common = [
            sys.executable, str(COORDINATOR), str(fake), str(source),
            str(campaign), "--expected-binary-sha256", binary_hash,
            "--expected-input-sha256", input_hash,
            "--target-measurement-cycles", "4", "--block-cycles", "2",
            "--qualified-seconds-per-cycle", "0.01",
            "--max-wall-seconds", "30", "--stage-timeout-seconds", "10",
            "--max-stage-initial-particle-steps", "1000000",
            "--max-stages-per-invocation", "1",
            "--acknowledge-cost", ACK,
        ]
        first = subprocess.run(common, text=True, capture_output=True)
        require(first.returncode == 0, first.stderr)
        first_report = json.loads(first.stdout)
        require(first_report["completed_measurement_cycles"] == 2 and
                not first_report["target_reached"] and
                first_report["stop_reason"] == "stage_count_limit_reached",
                "first bounded invocation did not stop after one block")
        second = subprocess.run(
            [*common, "--resume-existing"], text=True, capture_output=True)
        require(second.returncode == 0, second.stderr)
        report = json.loads(second.stdout)
        require(report["completed_measurement_cycles"] == 4 and
                report["target_reached"] and len(report["stages"]) == 2 and
                report["latest_state"]["cycles"] == 5,
                "resumed campaign did not reach its exact target")
        first_density = campaign / "stage-000001-000003" / "density.dat"
        first_density.write_text(
            first_density.read_text(encoding="utf-8") + "0 0 0\n",
            encoding="utf-8")
        rejected = subprocess.run(
            [*common, "--resume-existing"], text=True, capture_output=True)
        require(rejected.returncode == 2 and
                "output hash differs" in rejected.stderr,
                "resume accepted a changed native measurement output")
    print("eduPIC measurement advancement regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
