#!/usr/bin/env python3
"""Bounded regression for seeded Hall ensemble preparation and aggregation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "examples" / "hall_landmark_axial_azimuthal.case"
PREPARE = ROOT / "scripts" / "prepare_hall_ensemble.py"
AGGREGATE = ROOT / "scripts" / "aggregate_hall_ensemble.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hash_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def comparison_report(
    campaign: dict[str, object],
    run_record: dict[str, object],
    campaign_root: Path,
    offset: float,
) -> dict[str, object]:
    seed = int(run_record["seed"])
    runtime = campaign_root / str(run_record["runtime_config"])
    return {
        "schema_version": 1,
        "passed": True,
        "case_id": campaign["case_id"],
        "case_variant": "synthetic-ensemble-reference",
        "profile_axis": "x",
        "mode_axis": "y",
        "averaging_window": {
            "start_time": 1.0,
            "end_time": 2.0,
            "duration": 1.0,
            "samples": 11,
            "profile_axis": "x",
        },
        "reference": {
            "manifest_sha256": "1" * 64,
            "profile_sha256": "2" * 64,
            "mode_sha256": "3" * 64,
        },
        "simulation": {
            "output_dir": run_record["result_dir"],
            "case_manifest_sha256": campaign["case_manifest_sha256"],
            "runtime_config": str(runtime.resolve()),
            "runtime_config_sha256": run_record["runtime_config_sha256"],
            "seed": seed,
            "field_average_sha256": hash_text(f"field-{seed}"),
            "species_average_sha256": hash_text(f"species-{seed}"),
            "mode_history_sha256": hash_text(f"mode-{seed}"),
        },
        "profile_comparisons": [{
            "reference_coordinate": 0.0125,
            "observables": [{
                "name": "axial_field",
                "simulation": 100.0 + offset,
                "reference": 100.0,
                "reference_uncertainty": 1.0,
                "residual": offset,
                "absolute_residual": abs(offset),
                "acceptance_threshold": 10.0,
                "passed": abs(offset) <= 10.0,
            }],
        }],
        "mode_comparisons": [{
            "name": "dominant_frequency",
            "mode": 16,
            "metric": "frequency_hz",
            "simulation": 10.0 + 0.1 * offset,
            "reference": 10.0,
            "reference_uncertainty": 0.1,
            "residual": 0.1 * offset,
            "absolute_residual": abs(0.1 * offset),
            "acceptance_threshold": 1.0,
            "passed": abs(0.1 * offset) <= 1.0,
        }],
    }


def write_reports(
    campaign_path: Path,
    offsets: list[float],
) -> tuple[dict[str, object], Path]:
    campaign = json.loads(campaign_path.read_text(encoding="utf-8"))
    root = campaign_path.parent
    for run_record, offset in zip(
        campaign["runs"], offsets, strict=True
    ):
        report_path = root / str(run_record["comparison_report"])
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(
            json.dumps(
                comparison_report(
                    campaign, run_record, root, offset
                ),
                indent=2,
                sort_keys=True,
            ) + "\n",
            encoding="utf-8",
        )
    return campaign, root


def main() -> int:
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_hall_ensemble_"
    ) as temporary:
        work = Path(temporary)
        guarded_dir = work / "guarded"
        guarded = run([
            sys.executable,
            str(PREPARE),
            str(CASE),
            "--tier", "workstation",
            "--seeds", "11,22,33",
            "--output-dir", str(guarded_dir),
        ])
        require(
            guarded.returncode == 2 and not guarded_dir.exists()
            and "requires --acknowledge-cost" in guarded.stderr,
            "Hall ensemble preparation bypassed tier cost authorization",
        )
        ensemble_dir = work / "ensemble"
        prepared = run([
            sys.executable,
            str(PREPARE),
            str(CASE),
            "--tier", "micro",
            "--seeds", "101,202,303",
            "--output-dir", str(ensemble_dir),
        ])
        require(
            prepared.returncode == 0,
            f"Hall ensemble preparation failed: {prepared.stderr}",
        )
        campaign_path = ensemble_dir / "ensemble.json"
        campaign, root = write_reports(
            campaign_path, [-2.0, 0.0, 2.0]
        )
        require(
            campaign["seed_count"] == 3
            and campaign["seeds"] == [101, 202, 303]
            and not campaign["launched"]
            and all(
                digest(root / str(run_record["runtime_config"]))
                    == run_record["runtime_config_sha256"]
                for run_record in campaign["runs"]
            ),
            "prepared Hall ensemble manifest is incomplete",
        )
        output = work / "aggregate.json"
        aggregated = run([
            sys.executable,
            str(AGGREGATE),
            str(campaign_path),
            "--output", str(output),
        ])
        result = json.loads(output.read_text(encoding="utf-8"))
        require(
            aggregated.returncode == 0
            and result["passed"]
            and not result["physics_claim_eligible"]
            and result["seed_count"] == 3
            and result["required_individual_passes"] == 2
            and result["individual_reports_passed"] == 3
            and result["profile_ensemble"][0]["ensemble_mean"] == 100.0
            and result["profile_ensemble"][0]["confidence_half_width"] < 5.0,
            f"Hall ensemble aggregation failed: {aggregated.stderr}",
        )

        failed_dir = work / "failed_ensemble"
        failed_prepared = run([
            sys.executable,
            str(PREPARE),
            str(CASE),
            "--tier", "micro",
            "--seeds", "404,505,606",
            "--output-dir", str(failed_dir),
        ])
        require(failed_prepared.returncode == 0, "failed fixture preparation failed")
        failed_campaign = failed_dir / "ensemble.json"
        write_reports(failed_campaign, [-20.0, 0.0, 20.0])
        failed_output = work / "failed.json"
        failed = run([
            sys.executable,
            str(AGGREGATE),
            str(failed_campaign),
            "--output", str(failed_output),
        ])
        require(
            failed.returncode == 1
            and not json.loads(
                failed_output.read_text(encoding="utf-8")
            )["passed"],
            "Hall ensemble accepted excessive seed uncertainty",
        )

        duplicate_dir = work / "duplicate_ensemble"
        duplicate_prepared = run([
            sys.executable,
            str(PREPARE),
            str(CASE),
            "--tier", "micro",
            "--seeds", "707,808,909",
            "--output-dir", str(duplicate_dir),
        ])
        require(
            duplicate_prepared.returncode == 0,
            "duplicate fixture preparation failed",
        )
        duplicate_campaign_path = duplicate_dir / "ensemble.json"
        duplicate_campaign, duplicate_root = write_reports(
            duplicate_campaign_path, [-1.0, 0.0, 1.0]
        )
        first_path = duplicate_root / str(
            duplicate_campaign["runs"][0]["comparison_report"]
        )
        second_path = duplicate_root / str(
            duplicate_campaign["runs"][1]["comparison_report"]
        )
        first_report = json.loads(first_path.read_text(encoding="utf-8"))
        second_report = json.loads(second_path.read_text(encoding="utf-8"))
        for key in (
            "field_average_sha256",
            "species_average_sha256",
            "mode_history_sha256",
        ):
            second_report["simulation"][key] = first_report["simulation"][key]
        second_path.write_text(
            json.dumps(second_report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        duplicate_output = work / "duplicate.json"
        duplicate = run([
            sys.executable,
            str(AGGREGATE),
            str(duplicate_campaign_path),
            "--output", str(duplicate_output),
        ])
        require(
            duplicate.returncode == 2
            and not duplicate_output.exists()
            and "reuse identical simulation artifacts" in duplicate.stderr,
            "Hall ensemble accepted duplicated seed artifacts",
        )

    print("Hall seeded ensemble preparation and aggregation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
