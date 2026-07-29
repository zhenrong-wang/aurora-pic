#!/usr/bin/env python3
"""Bounded regression for external Hall source planning and checksum locking."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "examples" / "hall_landmark_case2.sources"
LOCKER = ROOT / "scripts" / "lock_hall_source.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def main() -> int:
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_hall_source_lock_"
    ) as temporary:
        work = Path(temporary)
        plan = work / "plan.json"
        planned = run([
            sys.executable, str(LOCKER), str(REGISTRY),
            "--source", "warpx_deepblue", "--output", str(plan),
        ])
        plan_value = json.loads(plan.read_text(encoding="utf-8"))
        require(
            planned.returncode == 0
            and not plan_value["artifact_locked"]
            and plan_value["status"] == "external_acquisition_required"
            and plan_value["expected_artifact_name"] == "baseline_20us.tar"
            and plan_value["cells_x"] == 512
            and plan_value["cells_y"] == 256
            and plan_value["file_set_id"] == "m900nv362",
            "Hall external acquisition plan is incomplete",
        )

        artifact = work / "baseline_20us.tar"
        artifact.write_bytes(b"synthetic bounded AMReX archive fixture\n")
        artifact_hash = hashlib.sha256(artifact.read_bytes()).hexdigest()
        blocked_lock = work / "blocked.json"
        blocked = run([
            sys.executable, str(LOCKER), str(REGISTRY),
            "--source", "warpx_deepblue",
            "--artifact", str(artifact),
            "--maximum-bytes", "4",
            "--output", str(blocked_lock),
        ])
        require(
            blocked.returncode == 2 and not blocked_lock.exists()
            and "default hash limit" in blocked.stderr,
            "Hall source locker bypassed its file-size guard",
        )

        lock = work / "source-lock.json"
        locked = run([
            sys.executable, str(LOCKER), str(REGISTRY),
            "--source", "warpx_deepblue",
            "--artifact", str(artifact),
            "--expected-sha256", artifact_hash,
            "--maximum-bytes", "4",
            "--acknowledge-large-file-hash",
            "I_UNDERSTAND_THIS_MAY_HASH_A_VERY_LARGE_FILE",
            "--output", str(lock),
        ])
        lock_value = json.loads(lock.read_text(encoding="utf-8"))
        require(
            locked.returncode == 0
            and lock_value["artifact_locked"]
            and lock_value["repository_checksum_verified"]
            and lock_value["artifact_sha256"] == artifact_hash
            and lock_value["status"] == "checksum_verified",
            "Hall artifact checksum lock is incomplete",
        )

        bad_lock = work / "bad-lock.json"
        rejected = run([
            sys.executable, str(LOCKER), str(REGISTRY),
            "--source", "warpx_deepblue",
            "--artifact", str(artifact),
            "--expected-sha256", "0" * 64,
            "--maximum-bytes", "1024",
            "--output", str(bad_lock),
        ])
        require(
            rejected.returncode == 2 and not bad_lock.exists()
            and "SHA-256 mismatch" in rejected.stderr,
            "Hall source locker accepted a mismatched repository checksum",
        )

    print("Hall external source planning and checksum lock passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
