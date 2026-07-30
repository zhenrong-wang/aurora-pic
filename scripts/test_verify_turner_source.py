#!/usr/bin/env python3
"""Bounded regression for Turner publisher-supplement verification."""

from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
VERIFIER = ROOT / "scripts" / "verify_turner_source.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def result_table() -> bytes:
    lines = ["# synthetic Turner-format regression fixture"]
    for case in range(1, 5):
        lines.extend([
            f"# Case {case}",
            "0 1 1 1 1 1 1",
            "0.067 1 1 1 1 1 1",
        ])
    return ("\n".join(lines) + "\n").encode()


def main() -> int:
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_source_"
    ) as temporary:
        work = Path(temporary)
        members = {
            "electron.dat": (
                "Generated on 28 Nov 2012.\nBiagi-v7.1\n"
                "ELASTIC\nEXCITATION\nEXCITATION\nIONIZATION\n"
            ).encode(),
            "ion.dat": (
                "# Centre of mass energy\n"
                "# Isotropic scattering\n# Backward scattering\n"
                "0 1 1\n10000 1 1\n"
            ).encode(),
            "refined.dat": result_table(),
            "results.dat": result_table(),
        }
        artifact = work / "fixture.zip"
        with zipfile.ZipFile(artifact, "w", zipfile.ZIP_DEFLATED) as archive:
            for name, data in members.items():
                archive.writestr(name, data)

        sections = []
        kinds = {
            "electron.dat": (
                "lxcat-electron-table",
                "process_counts = ELASTIC:1,EXCITATION:2,IONIZATION:1",
            ),
            "ion.dat": ("turner-ion-table", "rows = 2"),
            "refined.dat": (
                "turner-results", "case_rows = 1:2,2:2,3:2,4:2"
            ),
            "results.dat": (
                "turner-results", "case_rows = 1:2,2:2,3:2,4:2"
            ),
        }
        for name, data in members.items():
            kind, semantic = kinds[name]
            sections.append(
                f"[member.{name}]\n"
                f"bytes = {len(data)}\n"
                f"sha256 = {sha256(data)}\n"
                f"kind = {kind}\n"
                f"{semantic}\n"
            )
        registry = work / "fixture.sources"
        registry.write_text(
            "source_registry_version = 1\n"
            "case_id = synthetic-turner-regression\n"
            "doi = 10.0/synthetic\n"
            "landing_url = https://example.invalid\n\n"
            "[source.publisher_supplement]\n"
            "artifact_name = fixture.zip\n"
            f"artifact_bytes = {artifact.stat().st_size}\n"
            f"sha256 = {sha256(artifact.read_bytes())}\n"
            "acquisition = synthetic\n"
            "license = synthetic\n"
            "redistribution = test-only\n"
            f"members = {', '.join(members)}\n\n"
            + "\n".join(sections),
            encoding="utf-8",
        )

        report = work / "report.json"
        valid = subprocess.run(
            [sys.executable, str(VERIFIER), str(registry), str(artifact),
             "--output", str(report)],
            cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        require(valid.returncode == 0 and report.exists(),
                f"valid Turner source fixture was rejected: {valid.stderr}")

        tampered_dir = work / "tampered"
        tampered_dir.mkdir()
        tampered = tampered_dir / "fixture.zip"
        tampered.write_bytes(artifact.read_bytes() + b"tampered")
        invalid = subprocess.run(
            [sys.executable, str(VERIFIER), str(registry), str(tampered)],
            cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        require(invalid.returncode == 2 and "byte count differs" in invalid.stderr,
                "tampered Turner source fixture was not rejected")

    print("Turner publisher-supplement verification passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
