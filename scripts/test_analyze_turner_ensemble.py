#!/usr/bin/env python3
"""Regression for checksum-verified Turner ensemble aggregation."""

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
PREPARER = ROOT / "scripts" / "prepare_turner_ensemble.py"
COMPARATOR = ROOT / "scripts" / "compare_turner.py"
ATTACHER = ROOT / "scripts" / "attach_turner_ensemble_result.py"
ANALYZER = ROOT / "scripts" / "analyze_turner_ensemble.py"
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_ENSEMBLE"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command, cwd=ROOT, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=60,
    )


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_ensemble_analysis_", dir=ROOT / "tmp"
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
            case_text.replace(old_audit, sha256(audit_path)), encoding="utf-8"
        )

        ensemble = work / "ensemble"
        prepared = run([
            sys.executable, str(PREPARER), str(case_path), str(normalized),
            "--output-dir", str(ensemble), "--seeds", "13507,24680,97531",
            "--acknowledge-cost", ACKNOWLEDGEMENT,
        ])
        require(prepared.returncode == 0, prepared.stderr)
        manifest = ensemble / "ensemble.json"
        attachments: list[Path] = []
        checkpoints: list[Path] = []

        for seed, multiplier in ((13507, 1.001), (24680, 1.002),
                                 (97531, 0.999)):
            prepared_config = ensemble / f"seed_{seed}" / "turner_case1.cfg"
            result_dir = work / f"result-{seed}"
            result_dir.mkdir()
            executed_config = work / f"executed-{seed}.cfg"
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
                    writer.writerow([
                        1, "ions", index, x, multiplier * 1.0e14
                    ])
            metadata = result_dir / "spatial_average_metadata.json"
            metadata.write_text(json.dumps({
                "spatial_average_version": 1,
                "reset_on_restart": False,
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
            comparison = work / f"comparison-{seed}.json"
            compared = run([
                sys.executable, str(COMPARATOR), "--case", "1",
                "--reference", str(reference), "--candidate", str(profile),
                "--candidate-metadata", str(metadata),
                "--normalization-audit", str(audit_path),
                "--output", str(comparison),
            ])
            require(compared.returncode == 0, compared.stderr)
            checkpoint = result_dir / "checkpoint_512000.apc"
            checkpoint.write_bytes(f"checkpoint seed {seed}".encode())
            checkpoints.append(checkpoint)
            attachment = work / f"attachment-{seed}.json"
            attached = run([
                sys.executable, str(ATTACHER), str(manifest),
                "--seed", str(seed), "--executed-config", str(executed_config),
                "--comparison-report", str(comparison),
                "--final-checkpoint", str(checkpoint),
                "--output", str(attachment),
            ])
            require(attached.returncode == 0, attached.stderr)
            attachments.append(attachment)

        report_path = work / "ensemble-analysis.json"
        command = [
            sys.executable, str(ANALYZER), str(manifest), "--attachments",
            *(str(path) for path in attachments), "--output", str(report_path),
        ]
        analyzed = run(command)
        require(analyzed.returncode == 0, analyzed.stderr)
        report = json.loads(report_path.read_text(encoding="utf-8"))
        published = report["published_individual_run_results"]
        amplitude = report["density_amplitude"]
        require(
            report["turner_ensemble_analysis_version"] == 1
            and report["ensemble"]["complete"] is True
            and report["ensemble"]["verified_members"] == 3
            and published["accepted_95_percent_count"] == 2
            and published["accepted_99_percent_count"] == 2
            and published["failed_99_percent_count"] == 1
            and amplitude["all_member_biases_same_sign"] is False
            and abs(amplitude["integrated_bias_percent_mean"]
                    - 0.06666666666666667) < 1e-10
            and report["formal_ensemble_acceptance_rule"] == "none_predeclared"
            and report["physics_claim"].endswith("descriptive_evidence_only"),
            "Turner ensemble aggregate or claim boundary is incorrect",
        )

        checkpoints[0].write_bytes(b"tampered checkpoint")
        rejected = run(command[:-1] + [str(work / "rejected.json")])
        require(
            rejected.returncode == 2 and "checkpoint checksum mismatch"
            in rejected.stderr,
            "ensemble analysis accepted a modified completed-run artifact",
        )

    print("Turner ensemble analysis passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
