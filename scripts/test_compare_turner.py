#!/usr/bin/env python3
"""Bounded regression for the Turner ion-density comparator."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
COMPARATOR = ROOT / "scripts" / "compare_turner.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    project_tmp = ROOT / "tmp"
    project_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_compare_", dir=project_tmp
    ) as temporary:
        work = Path(temporary)
        reference = work / "turner_case1_benchmark.csv"
        candidate = work / "candidate.csv"
        reference_header = [
            "x_m", "electron_density_mean_m-3",
            "electron_density_mean_stddev_m-3",
            "electron_density_population_stddev_m-3",
            "ion_density_mean_m-3", "ion_density_mean_stddev_m-3",
            "ion_density_population_stddev_m-3",
        ]
        candidate_header = ["x_m", "ion_density_mean_m-3"]
        with reference.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(reference_header)
            for index in range(129):
                x = 0.067 * index / 128
                writer.writerow([x, 1, 1, 1, 1.0e14, 1, 1.0e11])
        with candidate.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(candidate_header)
            for index in range(129):
                x = 0.067 * index / 128
                writer.writerow([x, 1.001e14])
        audit = work / "audit.json"
        audit.write_text(json.dumps({
            "turner_normalization_version": 1,
            "normalized_files": {
                reference.name: {"sha256": sha256(reference)}
            },
        }), encoding="utf-8")
        report = work / "report.json"
        completed = subprocess.run([
            sys.executable, str(COMPARATOR), "--case", "1",
            "--reference", str(reference), "--candidate", str(candidate),
            "--normalization-audit", str(audit), "--output", str(report),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(completed.returncode == 0,
                f"Turner comparator rejected valid fixture: {completed.stderr}")
        value = json.loads(report.read_text(encoding="utf-8"))
        require(
            abs(value["statistic"]["x_squared"] - 129.0) < 1e-12
            and value["statistic"]["accepted_95_percent"]
            and value["statistic"]["accepted_99_percent"]
            and value["statistic"]["formula_variance"]
            == "population_standard_deviation_squared"
            and value["physics_claim"].startswith("none_"),
            "Turner comparator statistic or claim boundary is incorrect",
        )

        changed = reference.read_text(encoding="utf-8").replace(
            "100000000000000.0", "100000000000001.0", 1
        )
        reference.write_text(changed, encoding="utf-8")
        rejected = subprocess.run([
            sys.executable, str(COMPARATOR), "--case", "1",
            "--reference", str(reference), "--candidate", str(candidate),
            "--normalization-audit", str(audit),
            "--output", str(work / "bad.json"),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(rejected.returncode == 2 and "SHA-256 differs" in rejected.stderr,
                "Turner comparator accepted an unlocked reference")

    print("Turner ion-density comparison passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
