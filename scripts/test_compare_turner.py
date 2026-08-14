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
        candidate_header = [
            "species_id", "species", "node", "x_m",
            "number_density_mean_m-3",
        ]
        with reference.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(reference_header)
            for index in range(129):
                x = 0.067 * index / 128
                # Match the limited coordinate precision in the publisher
                # supplement rather than manufacturing exact binary equality.
                writer.writerow([
                    f"{x:.6g}", 1.0e14, 1, 1.0e11,
                    1.0e14, 1, 1.0e11
                ])
        with candidate.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(candidate_header)
            for index in range(129):
                x = 0.067 * index / 128
                writer.writerow([0, "electrons", index, x, 1.001e14])
                writer.writerow([1, "ions", index, x, 1.001e14])
        audit = work / "audit.json"
        audit.write_text(json.dumps({
            "turner_normalization_version": 1,
            "normalized_files": {
                reference.name: {"sha256": sha256(reference)}
            },
        }), encoding="utf-8")
        metadata = work / "spatial_average_metadata.json"
        metadata.write_text(json.dumps({
            "spatial_average_version": 6,
            "sampling_order": "post_collision",
            "moment_samples": 12800,
            "moments_complete": True,
            "unit_system": "si",
            "start_step": 499201,
            "end_step": 512000,
            "interval": 1,
            "samples": 12800,
            "expected_samples": 12800,
            "final_step": 512000,
            "dt": 1.0 / (13.56e6 * 400),
            "rf_frequency": 13.56e6,
            "rf_cycles": 32,
            "complete": True,
            "species": ["electrons", "ions"],
        }), encoding="utf-8")
        report = work / "report.json"
        completed = subprocess.run([
            sys.executable, str(COMPARATOR), "--case", "1",
            "--reference", str(reference), "--candidate", str(candidate),
            "--candidate-metadata", str(metadata),
            "--normalization-audit", str(audit), "--output", str(report),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(completed.returncode == 0,
                f"Turner comparator rejected valid fixture: {completed.stderr}")
        pre_collision = json.loads(metadata.read_text(encoding="utf-8"))
        pre_collision["sampling_order"] = "pre_collision"
        metadata.write_text(json.dumps(pre_collision), encoding="utf-8")
        rejected = subprocess.run([
            sys.executable, str(COMPARATOR), "--case", "1",
            "--reference", str(reference), "--candidate", str(candidate),
            "--candidate-metadata", str(metadata),
            "--normalization-audit", str(audit),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(rejected.returncode != 0,
                "Turner comparator accepted pre-collision diagnostics")
        pre_collision["sampling_order"] = "post_collision"
        metadata.write_text(json.dumps(pre_collision), encoding="utf-8")
        value = json.loads(report.read_text(encoding="utf-8"))
        require(
            abs(value["statistic"]["x_squared"] - 129.0) < 1e-12
            and value["statistic"]["accepted_95_percent"]
            and value["statistic"]["accepted_99_percent"]
            and value["statistic"]["formula_variance"]
            == "population_standard_deviation_squared"
            and value["candidate"]["species"] == "ions"
            and value["averaging_contract_verified"]
            and value["coordinate_contract"]["mapping"]
            == "ordered_prescribed_grid_no_interpolation"
            and value["coordinate_contract"]["maximum_reference_error_m"] > 0
            and value["physics_claim"].startswith("none_"),
            "Turner comparator statistic or claim boundary is incorrect",
        )

        electron_report = work / "electron-report.json"
        electron = subprocess.run([
            sys.executable, str(COMPARATOR), "--case", "1",
            "--reference", str(reference), "--candidate", str(candidate),
            "--candidate-metadata", str(metadata),
            "--normalization-audit", str(audit), "--species", "electrons",
            "--output", str(electron_report),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(electron.returncode == 0, electron.stderr)
        electron_result = json.loads(
            electron_report.read_text(encoding="utf-8")
        )
        require(
            abs(electron_result["statistic"]["x_squared"] - 129.0) < 1e-12
            and electron_result["statistic"]
                ["published_acceptance_applicable"] is False
            and "accepted_99_percent" not in electron_result["statistic"]
            and electron_result["comparison_scope"]
                == "published_baseline_electron_density_diagnostic_only"
            and electron_result["physics_claim"]
                == "none_published_electron_density_descriptive_only",
            "electron-density comparison overstated its claim",
        )

        diagnostic_metadata = work / "diagnostic_metadata.json"
        diagnostic_value = json.loads(
            metadata.read_text(encoding="utf-8")
        )
        diagnostic_value.update({
            "reset_on_restart": True,
            "start_step": 525201,
            "end_step": 538000,
            "final_step": 538000,
        })
        diagnostic_metadata.write_text(
            json.dumps(diagnostic_value), encoding="utf-8"
        )
        diagnostic_report = work / "diagnostic-report.json"
        diagnostic = subprocess.run([
            sys.executable, str(COMPARATOR), "--case", "1",
            "--reference", str(reference), "--candidate", str(candidate),
            "--candidate-metadata", str(diagnostic_metadata),
            "--normalization-audit", str(audit),
            "--post-benchmark-window",
            "--output", str(diagnostic_report),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(
            diagnostic.returncode == 0,
            f"Turner comparator rejected a valid post-benchmark "
            f"window: {diagnostic.stderr}",
        )
        diagnostic_result = json.loads(
            diagnostic_report.read_text(encoding="utf-8")
        )
        require(
            diagnostic_result["statistic"][
                "published_acceptance_applicable"
            ] is False
            and "accepted_99_percent"
                not in diagnostic_result["statistic"]
            and diagnostic_result["comparison_scope"]
                == "post_benchmark_density_diagnostic_only"
            and diagnostic_result["physics_claim"].startswith(
                "none_post_benchmark"
            ),
            "post-benchmark comparison overstated its claim",
        )

        refined_candidate = work / "refined-candidate.csv"
        with refined_candidate.open(
            "w", newline="", encoding="utf-8"
        ) as stream:
            writer = csv.writer(stream)
            writer.writerow(candidate_header)
            for index in range(257):
                x = 0.067 * index / 256
                writer.writerow([0, "electrons", index, x, 1.001e14])
                writer.writerow([1, "ions", index, x, 1.001e14])
        refined_metadata = work / "refined-metadata.json"
        refined_metadata.write_text(json.dumps({
            "spatial_average_version": 1,
            "unit_system": "si",
            "start_step": 998401,
            "end_step": 1024000,
            "interval": 1,
            "samples": 25600,
            "expected_samples": 25600,
            "final_step": 1024000,
            "dt": 1.0 / (13.56e6 * 800),
            "rf_frequency": 13.56e6,
            "rf_cycles": 32,
            "complete": True,
            "species": ["electrons", "ions"],
        }), encoding="utf-8")
        refined_report = work / "refined-report.json"
        refined = subprocess.run([
            sys.executable, str(COMPARATOR), "--case", "1",
            "--reference", str(reference),
            "--candidate", str(refined_candidate),
            "--candidate-metadata", str(refined_metadata),
            "--normalization-audit", str(audit),
            "--numerical-sensitivity", "--output", str(refined_report),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(
            refined.returncode == 0,
            f"Turner comparator rejected a valid numerical sensitivity: "
            f"{refined.stderr}",
        )
        refined_result = json.loads(
            refined_report.read_text(encoding="utf-8")
        )
        require(
            abs(refined_result["statistic"]["x_squared"] - 129.0) < 1e-12
            and refined_result["statistic"]
                ["published_acceptance_applicable"] is False
            and "accepted_99_percent" not in refined_result["statistic"]
            and refined_result["candidate_nodes"] == 257
            and refined_result["numerical_sensitivity_contract"]
                ["time_refinement_ratio"] == 2
            and refined_result["numerical_sensitivity_contract"]
                ["grid_refinement_ratio"] == 2
            and refined_result["comparison_scope"]
                == "numerical_sensitivity_density_diagnostic_only"
            and refined_result["physics_claim"]
                == "none_changed_published_numerical_contract",
            "numerical-sensitivity comparison overstated its claim",
        )

        changed = reference.read_text(encoding="utf-8").replace(
            "100000000000000.0", "100000000000001.0", 1
        )
        reference.write_text(changed, encoding="utf-8")
        rejected = subprocess.run([
            sys.executable, str(COMPARATOR), "--case", "1",
            "--reference", str(reference), "--candidate", str(candidate),
            "--candidate-metadata", str(metadata),
            "--normalization-audit", str(audit),
            "--output", str(work / "bad.json"),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(rejected.returncode == 2 and "SHA-256 differs" in rejected.stderr,
                "Turner comparator accepted an unlocked reference")

        invalid_reference = work / "turner_case1_benchmark.csv"
        rows = invalid_reference.read_text(encoding="utf-8").splitlines()
        columns = rows[1].split(",")
        columns[0] = "0.0001"
        rows[1] = ",".join(columns)
        invalid_reference.write_text(
            "\n".join(rows) + "\n", encoding="utf-8"
        )
        audit.write_text(json.dumps({
            "turner_normalization_version": 1,
            "normalized_files": {
                invalid_reference.name: {
                    "sha256": sha256(invalid_reference)
                }
            },
        }), encoding="utf-8")
        off_grid = subprocess.run([
            sys.executable, str(COMPARATOR), "--case", "1",
            "--reference", str(invalid_reference),
            "--candidate", str(candidate),
            "--candidate-metadata", str(metadata),
            "--normalization-audit", str(audit),
            "--output", str(work / "off-grid.json"),
        ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        require(
            off_grid.returncode == 2
            and "reference coordinate is off the prescribed grid"
            in off_grid.stderr,
            "Turner comparator accepted a reference on a different grid",
        )

    print("Turner ion-density comparison passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
