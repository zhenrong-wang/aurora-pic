#!/usr/bin/env python3
"""Bounded regression for the Turner density-block analyzer."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
ANALYZER = ROOT / "scripts" / "analyze_turner_density_blocks.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(work: Path, index: int, start: int) -> Path:
    profile = work / f"profile-{index}.csv"
    with profile.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "species_id", "species", "node", "x_m",
            "number_density_mean_m-3",
        ])
        for node in range(5):
            writer.writerow([
                0, "electrons", node, node * 0.01, 1.0e14
            ])
            writer.writerow([
                1, "ions", node, node * 0.01,
                (1.0 + 0.01 * index + 0.001 * node) * 1.0e14,
            ])
    metadata = work / f"metadata-{index}.json"
    metadata.write_text(json.dumps({
        "complete": True,
        "reset_on_restart": True,
        "start_step": start,
        "end_step": start + 12799,
        "final_step": start + 12799,
        "interval": 1,
        "samples": 12800,
        "expected_samples": 12800,
        "rf_cycles": 32,
    }), encoding="utf-8")
    report = work / f"comparison-{index}.json"
    report.write_text(json.dumps({
        "turner_comparison_version": 1,
        "case": 1,
        "comparison_scope": "post_benchmark_density_diagnostic_only",
        "candidate": {
            "path": str(profile.resolve()),
            "sha256": sha256(profile),
            "species": "ions",
            "averaging_metadata_path": str(metadata.resolve()),
            "averaging_metadata_sha256": sha256(metadata),
        },
        "statistic": {
            "x_squared": 100.0 + index,
            "published_acceptance_applicable": False,
        },
        "secondary_metrics": {"relative_l2": 0.01 * index},
    }), encoding="utf-8")
    return report


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_blocks_", dir=ROOT / "tmp"
    ) as temporary:
        work = Path(temporary)
        reports = [
            fixture(work, index, 525201 + index * 12800)
            for index in range(4)
        ]
        output = work / "analysis.json"
        completed = subprocess.run([
            sys.executable, str(ANALYZER), *map(str, reports),
            "--minimum-blocks", "4", "--output", str(output),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(
            completed.returncode == 0,
            f"density-block analyzer rejected valid fixtures: "
            f"{completed.stderr}",
        )
        value = json.loads(output.read_text(encoding="utf-8"))
        require(
            value["diagnostic_series_ready"] is True
            and value["classification"].startswith("diagnostic_series_ready")
            and len(value["blocks"]) == 4
            and len(value["adjacent_profile_relative_l2"]) == 3
            and value["series_metrics"][
                "projected_fractional_drift_across_series"
            ] > 0.0
            and value["series_metrics"]["ar1_effective_blocks"] is not None
            and value["published_acceptance_applicable"] is False
            and value["physics_claim"].startswith("none_"),
            "density-block analysis or claim boundary is incorrect",
        )

        insufficient = work / "insufficient.json"
        short = subprocess.run([
            sys.executable, str(ANALYZER), *map(str, reports[:2]),
            "--minimum-blocks", "4", "--output", str(insufficient),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        short_value = json.loads(insufficient.read_text(encoding="utf-8"))
        require(
            short.returncode == 0
            and short_value["classification"]
                == "insufficient_consecutive_blocks"
            and short_value["series_metrics"][
                "lag_one_integrated_density_correlation"
            ] is None,
            "short density series was not conservatively classified",
        )

        broken_metadata = json.loads(
            (work / "metadata-2.json").read_text(encoding="utf-8")
        )
        broken_metadata["start_step"] += 1
        broken_metadata["end_step"] += 1
        broken_metadata["final_step"] += 1
        (work / "metadata-2.json").write_text(
            json.dumps(broken_metadata), encoding="utf-8"
        )
        broken_report = json.loads(reports[2].read_text(encoding="utf-8"))
        broken_report["candidate"]["averaging_metadata_sha256"] = sha256(
            work / "metadata-2.json"
        )
        reports[2].write_text(json.dumps(broken_report), encoding="utf-8")
        rejected = subprocess.run([
            sys.executable, str(ANALYZER), *map(str, reports),
            "--minimum-blocks", "4", "--output", str(work / "rejected.json"),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(
            rejected.returncode == 2
            and "not contiguous" in rejected.stderr,
            "density-block analyzer accepted a discontinuous series",
        )

    print("Turner density-block analysis passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
