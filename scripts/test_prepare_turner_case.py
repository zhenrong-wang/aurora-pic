#!/usr/bin/env python3
"""Conservative regression for Turner Case 1 campaign preparation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
PREPARER = ROOT / "scripts" / "prepare_turner_case.py"
CASE = ROOT / "examples" / "turner_helium_ccp_case1.case"
ACKNOWLEDGEMENT = "I_UNDERSTAND_THIS_IS_A_PRODUCTION_SCALE_TURNER_RUN"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def identity(path: Path) -> dict:
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
    )


def create_fixture(work: Path) -> tuple[Path, Path, str]:
    normalized = work / "normalized"
    normalized.mkdir()
    table = "# synthetic regression input\n0 1e-21\n10000 1e-23\n"
    (normalized / "electron.dat").write_text(table, encoding="utf-8")
    (normalized / "ion.dat").write_text(table, encoding="utf-8")
    electron_gas = """gas_data_version = 2
units = si
gas = He
neutral_mass = 6.67e-27
dataset_id = synthetic.turner.electron
dataset_version = test
data_provenance = synthetic regression
citation = none
retrieved = 2026-07-30
license = synthetic

[collision.elastic]
type = elastic
cross_section_file = electron.dat
energy_scale = 1.602176634e-19
angular_model = isotropic

[collision.ionization]
type = ionization
cross_section_file = electron.dat
energy_scale = 1.602176634e-19
threshold_energy = 3.9e-18
angular_model = isotropic
"""
    ion_gas = """gas_data_version = 2
units = si
gas = He
neutral_mass = 6.67e-27
dataset_id = synthetic.turner.ion
dataset_version = test
data_provenance = synthetic regression
citation = none
retrieved = 2026-07-30
license = synthetic

[collision.elastic]
type = elastic
cross_section_file = ion.dat
energy_scale = 1.602176634e-19
energy_frame = center_of_mass
angular_model = isotropic
"""
    (normalized / "turner_he_electron.gas").write_text(
        electron_gas, encoding="utf-8"
    )
    (normalized / "turner_he_ion.gas").write_text(
        ion_gas, encoding="utf-8"
    )
    (normalized / "turner_case1_benchmark.csv").write_text(
        "case,node,x_m,ion_mean_m3,ion_population_stddev_m3\n"
        "1,0,0,1,1\n",
        encoding="utf-8",
    )
    names = (
        "electron.dat", "ion.dat", "turner_he_electron.gas",
        "turner_he_ion.gas", "turner_case1_benchmark.csv",
    )
    audit = {
        "turner_normalization_version": 2,
        "case_id": "turner-helium-ccp-2013",
        "source_artifact": {
            "sha256":
                "a0a5fe93900d7d7b213157f1eab664e06aab6e718f2189910e65f23bd699d661"
        },
        "normalized_files": {
            name: identity(normalized / name) for name in names
        },
    }
    audit_path = normalized / "audit.json"
    audit_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    case_text = CASE.read_text(encoding="utf-8")
    pinned_audit = next(
        line.split("=", 1)[1].strip()
        for line in case_text.splitlines()
        if line.startswith("normalized_audit_sha256")
    )
    case_text = case_text.replace(
        pinned_audit, hashlib.sha256(audit_path.read_bytes()).hexdigest()
    )
    case_path = work / "case.case"
    case_path.write_text(case_text, encoding="utf-8")
    return normalized, case_path, table


def main() -> int:
    project_tmp = ROOT / "tmp"
    project_tmp.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_turner_prepare_", dir=project_tmp
    ) as temporary:
        work = Path(temporary)
        normalized, case_path, table = create_fixture(work)
        output = work / "campaign" / "turner.cfg"

        rejected = run([
            sys.executable, str(PREPARER), str(case_path), str(normalized),
            "--output", str(output),
        ])
        require(
            rejected.returncode == 2
            and ACKNOWLEDGEMENT in rejected.stderr
            and not output.exists(),
            "campaign generation bypassed its explicit acknowledgement",
        )
        completed = run([
            sys.executable, str(PREPARER), str(case_path), str(normalized),
            "--output", str(output),
            "--acknowledge-cost", ACKNOWLEDGEMENT,
        ])
        require(
            completed.returncode == 0,
            "campaign preparation failed: " + completed.stdout + completed.stderr,
        )
        report_path = output.with_suffix(".preflight.json")
        report = json.loads(report_path.read_text(encoding="utf-8"))
        contract = report["contract"]
        collision = report["collision_guard"]
        require(
            not report["full_run_launched"]
            and not report["production_launch_authorized"]
            and report["physics_claim"] == "none"
            and contract["steps"] == 512000
            and contract["rf_cycles"] == 1280
            and contract["particles_per_species"] == 65536
            and contract["averaging_samples"] == 12800
            and report["reported_case1_characteristics"][
                "total_macro_particles"
            ] == 31900
            and collision["electron"]["majorant_dt"] < 0.1
            and collision["ion"]["majorant_dt"] < 0.1
            and collision["electron"]["configured_majorant_s"]
            > collision["electron"]["sampled_peak_frequency_s"]
            and "runtime_backend = serial" in output.read_text(encoding="utf-8")
            and "spatial_average_start_step = 499201"
            in output.read_text(encoding="utf-8"),
            "campaign contract or preflight report is incomplete",
        )

        seeded_output = work / "seeded" / "turner.cfg"
        seeded = run([
            sys.executable, str(PREPARER), str(case_path), str(normalized),
            "--output", str(seeded_output), "--seed", "24680",
            "--acknowledge-cost", ACKNOWLEDGEMENT,
        ])
        seeded_report = json.loads(
            seeded_output.with_suffix(".preflight.json").read_text(
                encoding="utf-8"
            )
        )
        require(
            seeded.returncode == 0
            and seeded_report["contract"]["seed"] == 24680
            and "seed = 24680" in seeded_output.read_text(encoding="utf-8")
            and seeded_report["provenance"]["generated_deck_sha256"]
                == hashlib.sha256(seeded_output.read_bytes()).hexdigest(),
            "campaign seed override is not hash-consistent",
        )

        invalid_seed = run([
            sys.executable, str(PREPARER), str(case_path), str(normalized),
            "--output", str(work / "invalid-seed.cfg"),
            "--seed", "4294967296",
            "--acknowledge-cost", ACKNOWLEDGEMENT,
        ])
        require(
            invalid_seed.returncode == 2
            and "unsigned 32-bit" in invalid_seed.stderr,
            "campaign preparation accepted an invalid seed",
        )

        (normalized / "electron.dat").write_text(
            table + "20000 1e-22\n", encoding="utf-8"
        )
        tampered = run([
            sys.executable, str(PREPARER), str(case_path), str(normalized),
            "--output", str(work / "tampered.cfg"),
            "--acknowledge-cost", ACKNOWLEDGEMENT,
        ])
        require(
            tampered.returncode == 2
            and "mismatch" in tampered.stderr,
            "campaign preparation accepted a modified normalized table",
        )
    print("Turner campaign preparation regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
