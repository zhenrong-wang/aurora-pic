#!/usr/bin/env python3
"""Conservative end-to-end tests for serialized swarm campaigns."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "run_swarm_campaign.py"
COMPARATOR = ROOT / "scripts" / "compare_swarm.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_result(path: Path, drift: float, signature: str = "12345") -> None:
    path.write_text(
        "dataset_id,dataset_version,gas,population_model,"
        "collision_model_signature,reduced_field_td,"
        "electron_drift_velocity_m_s,"
        "mean_velocity_x_standard_error_m_s\n"
        "synthetic.dataset,1,synthetic_swarm_gas,"
        f"fixed_population_no_avalanche,{signature},1,{drift},0.5\n",
        encoding="utf-8",
    )


def run(
    manifest: Path,
    executable: Path,
    report: Path,
    *,
    overwrite: bool = False,
) -> subprocess.CompletedProcess[str]:
    command = [
        sys.executable,
        str(RUNNER),
        str(manifest),
        "--swarm-executable",
        str(executable),
        "--comparator",
        str(COMPARATOR),
        "--output",
        str(report),
    ]
    if overwrite:
        command.append("--overwrite")
    return subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def main() -> int:
    with tempfile.TemporaryDirectory(
        prefix="aurorapic_swarm_campaign_"
    ) as temporary:
        work = Path(temporary)
        executable = work / "fake_swarm.py"
        executable.write_text(
            "#!/usr/bin/env python3\n"
            "from pathlib import Path\n"
            "import sys\n"
            "values = {}\n"
            "for line in Path(sys.argv[1]).read_text().splitlines():\n"
            "    line = line.split('#', 1)[0].strip()\n"
            "    if line:\n"
            "        key, value = line.split('=', 1)\n"
            "        values[key.strip()] = value.strip()\n"
            "Path(values['output_file']).write_bytes("
            "Path(values['fixture_file']).read_bytes())\n",
            encoding="utf-8",
        )
        executable.chmod(0o755)

        coarse_fixture = work / "coarse_fixture.csv"
        fine_fixture = work / "fine_fixture.csv"
        write_result(coarse_fixture, 99.0)
        write_result(fine_fixture, 100.0)
        coarse_result = work / "coarse.csv"
        fine_result = work / "fine.csv"
        coarse_config = work / "coarse.swarm"
        fine_config = work / "fine.swarm"
        coarse_config.write_text(
            f"fixture_file = {coarse_fixture}\n"
            f"output_file = {coarse_result}\n",
            encoding="utf-8",
        )
        fine_config.write_text(
            f"fixture_file = {fine_fixture}\n"
            f"output_file = {fine_result}\n",
            encoding="utf-8",
        )

        reference_data = work / "reference.csv"
        reference_data.write_text(
            "reduced_field_td,drift_velocity_m_s\n"
            "1,100\n",
            encoding="utf-8",
        )
        reference_manifest = work / "reference.swarm-reference"
        reference_manifest.write_text(
            "[reference]\n"
            "swarm_reference_version = 1\n"
            f"data_file = {reference_data.name}\n"
            "reference_id = aurorapic.synthetic.campaign\n"
            "reference_version = 1\n"
            "gas = synthetic_swarm_gas\n"
            "population_model = fixed_population_no_avalanche\n"
            "coefficient_convention = flux_fixed_population\n"
            "provenance = AuroraPIC synthetic campaign fixture\n"
            "citation = AuroraPIC synthetic fixture\n"
            "retrieved = 2026-07-28\n"
            "license = Synthetic test data\n"
            "\n"
            "[observable.drift]\n"
            "simulation_column = electron_drift_velocity_m_s\n"
            "reference_column = drift_velocity_m_s\n"
            "relative_tolerance = 0.2\n"
            "absolute_tolerance = 0\n"
            "uncertainty_multiplier = 0\n",
            encoding="utf-8",
        )
        manifest = work / "campaign.swarm-campaign"
        manifest.write_text(
            "[campaign]\n"
            "swarm_campaign_version = 1\n"
            "campaign_id = aurorapic.synthetic.convergence\n"
            "campaign_version = 1\n"
            "provenance = AuroraPIC synthetic campaign test\n"
            "retrieved = 2026-07-28\n"
            f"reference_manifest = {reference_manifest.name}\n"
            "run_order = coarse, fine\n"
            "reference_run = fine\n"
            "\n"
            "[run.coarse]\n"
            f"config_file = {coarse_config.name}\n"
            f"result_file = {coarse_result.name}\n"
            "\n"
            "[run.fine]\n"
            f"config_file = {fine_config.name}\n"
            f"result_file = {fine_result.name}\n"
            "\n"
            "[observable.drift]\n"
            "simulation_column = electron_drift_velocity_m_s\n"
            "uncertainty_column = mean_velocity_x_standard_error_m_s\n"
            "relative_tolerance = 0.02\n"
            "absolute_tolerance = 0\n"
            "uncertainty_multiplier = 1\n",
            encoding="utf-8",
        )
        report = work / "campaign.json"
        passed = run(manifest, executable, report)
        require(
            passed.returncode == 0,
            f"valid campaign failed: {passed.stderr}",
        )
        result = json.loads(report.read_text(encoding="utf-8"))
        require(
            result["passed"]
            and result["reference_validation_passed"]
            and result["convergence_validation_passed"]
            and len(result["runs"]) == 2
            and result["campaign"]["execution_policy"]["run_order"]
            == "serial"
            and result["campaign"]["execution_policy"]["OMP_NUM_THREADS"]
            == "1"
            and len(result["runs"][0]["result_sha256"]) == 64,
            "passing campaign report is incomplete",
        )

        existing = run(manifest, executable, report)
        require(
            existing.returncode == 2
            and "already exist" in existing.stderr,
            "campaign overwrote artifacts without opt-in",
        )
        overwritten = run(
            manifest, executable, report, overwrite=True
        )
        require(
            overwritten.returncode == 0,
            "campaign rejected explicit overwrite",
        )

        write_result(coarse_fixture, 90.0)
        failed = run(manifest, executable, report, overwrite=True)
        require(
            failed.returncode == 1
            and "did not meet" in failed.stderr,
            "non-converged campaign did not fail",
        )
        failed_result = json.loads(report.read_text(encoding="utf-8"))
        require(
            not failed_result["passed"]
            and failed_result["reference_validation_passed"]
            and not failed_result["convergence_validation_passed"],
            "campaign did not distinguish reference and convergence failures",
        )

        write_result(coarse_fixture, 99.0, signature="different")
        invalid = run(manifest, executable, report, overwrite=True)
        require(
            invalid.returncode == 2
            and "physics identity differs" in invalid.stderr,
            "campaign accepted inconsistent collision models",
        )

    print("swarm campaign validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
