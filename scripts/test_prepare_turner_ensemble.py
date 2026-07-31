#!/usr/bin/env python3
"""Conservative regression for full-duration Turner ensemble preparation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from test_prepare_turner_case import create_fixture


ROOT = Path(__file__).resolve().parents[1]
PREPARER = ROOT / "scripts" / "prepare_turner_ensemble.py"
ACKNOWLEDGEMENT = (
    "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_ENSEMBLE"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command, cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60,
    )


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_ensemble_", dir=ROOT / "tmp"
    ) as temporary:
        work = Path(temporary)
        normalized, case_path, _ = create_fixture(work)
        destination = work / "ensemble"

        rejected = run([
            sys.executable, str(PREPARER), str(case_path), str(normalized),
            "--output-dir", str(destination),
            "--seeds", "13507,24680,97531",
        ])
        require(
            rejected.returncode == 2
            and ACKNOWLEDGEMENT in rejected.stderr
            and not destination.exists(),
            "Turner ensemble bypassed its cost acknowledgement",
        )

        completed = run([
            sys.executable, str(PREPARER), str(case_path), str(normalized),
            "--output-dir", str(destination),
            "--seeds", "13507,24680,97531",
            "--acknowledge-cost", ACKNOWLEDGEMENT,
        ])
        require(
            completed.returncode == 0,
            "Turner ensemble preparation failed: "
            + completed.stdout + completed.stderr,
        )
        manifest = json.loads(
            (destination / "ensemble.json").read_text(encoding="utf-8")
        )
        require(
            manifest["turner_ensemble_preparation_version"] == 1
            and manifest["seed_count"] == 3
            and manifest["seeds"] == [13507, 24680, 97531]
            and manifest["independent_seed_contract"] is True
            and manifest["launched"] is False
            and manifest["launch_policy"] == "external_sequential_only"
            and manifest["aggregate_resource_floor"][
                "concurrent_runs_authorized"
            ] == 1
            and manifest["claim_boundary"]["physics_claim"].startswith("none_"),
            "Turner ensemble manifest contract is incomplete",
        )
        deck_hashes = set()
        for run_record in manifest["runs"]:
            seed = run_record["seed"]
            deck = destination / run_record["runtime_config"]
            preflight = destination / run_record["preflight_report"]
            deck_text = deck.read_text(encoding="utf-8")
            report = json.loads(preflight.read_text(encoding="utf-8"))
            require(
                run_record["runtime_config_sha256"] == sha256(deck)
                and run_record["preflight_report_sha256"] == sha256(preflight)
                and f"seed = {seed}" in deck_text
                and report["contract"]["seed"] == seed
                and report["full_run_launched"] is False
                and str(destination / "results" / f"seed_{seed}")
                    in deck_text,
                f"Turner ensemble seed {seed} is not hash-consistent",
            )
            deck_hashes.add(sha256(deck))
        require(len(deck_hashes) == 3, "Turner seed decks are not distinct")

        overwrite = run([
            sys.executable, str(PREPARER), str(case_path), str(normalized),
            "--output-dir", str(destination),
            "--seeds", "13507,24680,97531",
            "--acknowledge-cost", ACKNOWLEDGEMENT,
        ])
        require(
            overwrite.returncode == 2 and "refusing to overwrite" in overwrite.stderr,
            "Turner ensemble preparer overwrote an existing campaign",
        )

        duplicate = run([
            sys.executable, str(PREPARER), str(case_path), str(normalized),
            "--output-dir", str(work / "duplicate"),
            "--seeds", "13507,13507,97531",
            "--acknowledge-cost", ACKNOWLEDGEMENT,
        ])
        require(
            duplicate.returncode == 2
            and "unique unsigned 32-bit" in duplicate.stderr,
            "Turner ensemble accepted duplicate seeds",
        )

    print("Turner ensemble preparation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
