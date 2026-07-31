#!/usr/bin/env python3
"""Regression for checksum-verified Turner ensemble result attachment."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from test_prepare_turner_case import create_fixture, identity


ROOT = Path(__file__).resolve().parents[1]
ENSEMBLE_PREPARER = ROOT / "scripts" / "prepare_turner_ensemble.py"
COMPARATOR = ROOT / "scripts" / "compare_turner.py"
ATTACHER = ROOT / "scripts" / "attach_turner_ensemble_result.py"
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_ENSEMBLE"


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
        prefix="aurorapic_turner_attach_", dir=ROOT / "tmp"
    ) as temporary:
        work = Path(temporary)
        normalized, case_path, _ = create_fixture(work)
        reference = normalized / "turner_case1_benchmark.csv"
        with reference.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow([
                "x_m", "electron_density_mean_m-3",
                "electron_density_mean_stddev_m-3",
                "electron_density_population_stddev_m-3",
                "ion_density_mean_m-3", "ion_density_mean_stddev_m-3",
                "ion_density_population_stddev_m-3",
            ])
            for index in range(129):
                writer.writerow([
                    f"{0.067 * index / 128:.6g}", 1, 1, 1,
                    1.0e14, 1, 1.0e11,
                ])
        audit_path = normalized / "audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        audit["normalized_files"][reference.name] = identity(reference)
        audit_path.write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        case_text = case_path.read_text(encoding="utf-8")
        old_audit = next(
            line.split("=", 1)[1].strip()
            for line in case_text.splitlines()
            if line.startswith("normalized_audit_sha256")
        )
        case_path.write_text(
            case_text.replace(old_audit, sha256(audit_path)),
            encoding="utf-8",
        )

        ensemble = work / "ensemble"
        prepared = run([
            sys.executable, str(ENSEMBLE_PREPARER),
            str(case_path), str(normalized),
            "--output-dir", str(ensemble),
            "--seeds", "13507,24680,97531",
            "--acknowledge-cost", ACKNOWLEDGEMENT,
        ])
        require(prepared.returncode == 0, prepared.stderr)
        manifest = ensemble / "ensemble.json"
        manifest_hash = sha256(manifest)
        prepared_config = ensemble / "seed_13507" / "turner_case1.cfg"

        result_dir = work / "completed-seed-13507"
        result_dir.mkdir()
        executed_config = work / "executed.cfg"
        config_text = prepared_config.read_text(encoding="utf-8")
        old_output = next(
            line for line in config_text.splitlines()
            if line.startswith("output_dir = ")
        )
        executed_config.write_text(
            config_text.replace(old_output, f"output_dir = {result_dir}"),
            encoding="utf-8",
        )
        profile = result_dir / "spatial_average.csv"
        with profile.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow([
                "species_id", "species", "node", "x_m",
                "number_density_mean_m-3",
            ])
            for index in range(129):
                x = 0.067 * index / 128
                writer.writerow([0, "electrons", index, x, 1.0e13])
                writer.writerow([1, "ions", index, x, 1.001e14])
        metadata = result_dir / "spatial_average_metadata.json"
        metadata.write_text(json.dumps({
            "spatial_average_version": 1,
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
        comparison = work / "comparison.json"
        compared = run([
            sys.executable, str(COMPARATOR), "--case", "1",
            "--reference", str(reference), "--candidate", str(profile),
            "--candidate-metadata", str(metadata),
            "--normalization-audit", str(audit_path),
            "--output", str(comparison),
        ])
        require(compared.returncode == 0, compared.stderr)
        checkpoint = result_dir / "checkpoint_512000.apc"
        checkpoint.write_bytes(b"synthetic nonempty checkpoint")

        attachment = work / "attachment.json"
        attached = run([
            sys.executable, str(ATTACHER), str(manifest),
            "--seed", "13507", "--executed-config", str(executed_config),
            "--comparison-report", str(comparison),
            "--final-checkpoint", str(checkpoint),
            "--output", str(attachment),
        ])
        require(attached.returncode == 0, attached.stderr)
        value = json.loads(attachment.read_text(encoding="utf-8"))
        require(
            value["turner_ensemble_attachment_version"] == 1
            and value["seed"] == 13507
            and value["ensemble_member_attached"] is True
            and value["ensemble_complete"] is False
            and value["executed_run"]["semantic_equivalence"]
                == "exact_except_output_dir"
            and value["comparison"]["independently_recomputed"] is True
            and abs(value["comparison"]["x_squared"] - 129.0) < 1e-12
            and value["classification"].endswith("passed_99_percent")
            and value["physics_claim"].startswith("single_seed_")
            and sha256(manifest) == manifest_hash,
            "attached Turner result or claim boundary is incorrect",
        )

        changed_config = work / "changed.cfg"
        changed_config.write_text(
            executed_config.read_text(encoding="utf-8").replace(
                "phi_right_amplitude = 450", "phi_right_amplitude = 451"
            ),
            encoding="utf-8",
        )
        rejected = run([
            sys.executable, str(ATTACHER), str(manifest),
            "--seed", "13507", "--executed-config", str(changed_config),
            "--comparison-report", str(comparison),
            "--final-checkpoint", str(checkpoint),
            "--output", str(work / "rejected.json"),
        ])
        require(
            rejected.returncode == 2
            and "beyond output_dir" in rejected.stderr,
            "attachment accepted a changed executed physics contract",
        )

    print("Turner ensemble result attachment passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
