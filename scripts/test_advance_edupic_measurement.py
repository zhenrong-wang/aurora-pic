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
            "--min-free-disk-mib", "1",
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
        disk_campaign = work / "disk-campaign"
        disk_common = list(common)
        disk_common[4] = str(disk_campaign)
        minimum_index = disk_common.index("--min-free-disk-mib") + 1
        disk_common[minimum_index] = "1000000000000"
        disk_rejected = subprocess.run(
            disk_common, text=True, capture_output=True)
        disk_report = json.loads(disk_rejected.stdout)
        require(disk_rejected.returncode == 0 and
                disk_report["completed_measurement_cycles"] == 0 and
                disk_report["stop_reason"] == "host_free_disk_below_minimum" and
                disk_report["host_health_checks"][-1]["free_disk_mib"] > 0,
                "measurement campaign ignored its free-disk guard")

        amendment_campaign = work / "amendment-campaign"
        amendment_common = list(common)
        amendment_common[4] = str(amendment_campaign)
        amendment_common.extend(["--min-available-memory-mib", "1"])
        amendment_first = subprocess.run(
            amendment_common, text=True, capture_output=True)
        require(amendment_first.returncode == 0, amendment_first.stderr)
        lower_memory = list(amendment_common)
        lower_memory[lower_memory.index("--min-available-memory-mib") + 1] = "0.5"
        unacknowledged = subprocess.run(
            [*lower_memory, "--resume-existing"], text=True,
            capture_output=True)
        require(unacknowledged.returncode == 2 and
                "requires --acknowledge-memory-guard-amendment" in
                unacknowledged.stderr,
                "campaign silently accepted a lower memory guard")
        amended = subprocess.run([
            *lower_memory, "--resume-existing",
            "--acknowledge-memory-guard-amendment",
            "I_UNDERSTAND_THIS_RELAXES_AN_EXISTING_CAMPAIGN_MEMORY_GUARD",
            "--memory-guard-amendment-reason", "synthetic regression",
        ], text=True, capture_output=True)
        require(amended.returncode == 0, amended.stderr)
        amended_report = json.loads(amended.stdout)
        amendments = amended_report.get("operational_policy_amendments", [])
        require(amended_report["target_reached"] and len(amendments) == 1 and
                amendments[0]["previous_value_mib"] == 1.0 and
                amendments[0]["new_value_mib"] == 0.5 and
                amendments[0]["scientific_analysis_contract_changed"] is False,
                "acknowledged memory-guard amendment was not audited")
    print("eduPIC measurement advancement regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
