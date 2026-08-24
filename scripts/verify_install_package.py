#!/usr/bin/env python3
"""Smoke-test AuroraPIC install and package artifacts.

This script validates the release artifact path that static checks cannot cover:

1. install the already-built tree into a temporary prefix;
2. verify both simulation and swarm CLIs and run the simulation CLI;
3. configure, build, and run a downstream CMake consumer using
   ``find_package(AuroraPIC CONFIG REQUIRED)``;
4. build and inspect the CPack TGZ package for the same install surface.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]


class InstallSmokeError(RuntimeError):
    """Raised when an install/package smoke assertion fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise InstallSmokeError(message)


def run(command: list[str], *, cwd: Path | None = None) -> None:
    print("+ " + " ".join(command))
    subprocess.run(command, cwd=cwd, check=True)


def require_file(path: Path) -> Path:
    require(path.is_file(), f"missing expected file: {path}")
    require(path.stat().st_size > 0, f"expected non-empty file: {path}")
    return path


def require_dir(path: Path) -> Path:
    require(path.is_dir(), f"missing expected directory: {path}")
    return path


def write_cli_smoke_config(path: Path, output_dir: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "config_version = 1",
                "nx = 16",
                "length = 1.0",
                "dt = 0.001",
                "steps = 2",
                "output_interval = 1",
                "boundary = periodic",
                "mode = transient",
                "seed = 7",
                f"output_dir = {output_dir.as_posix()}",
                "",
                "[species]",
                "name = smoke_electrons",
                "charge = -1.0",
                "mass = 1.0",
                "weight = 0.01",
                "particles = 16",
                "drift_velocity = 0.0",
                "thermal_velocity = 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def smoke_installed_cli(prefix: Path, work: Path) -> None:
    cli = require_file(prefix / "bin" / "aurorapic_cli")
    require_file(prefix / "bin" / "aurorapic_swarm")
    output_dir = work / "cli-output"
    config = work / "install_smoke.cfg"
    write_cli_smoke_config(config, output_dir)
    run([str(cli), str(config)])
    scalars = require_file(output_dir / "scalars.csv")
    rows = scalars.read_text(encoding="utf-8").strip().splitlines()
    require(len(rows) >= 2, "installed CLI smoke should write scalar header plus data rows")
    require(rows[0].startswith("step,time,"), "installed CLI smoke scalar header is unexpected")


def write_consumer_project(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    (path / "CMakeLists.txt").write_text(
        "\n".join(
            [
                "cmake_minimum_required(VERSION 3.20)",
                "project(AuroraPICConsumer LANGUAGES CXX)",
                "find_package(AuroraPIC CONFIG REQUIRED)",
                "add_executable(consumer main.cpp)",
                "target_link_libraries(consumer PRIVATE AuroraPIC::aurorapic)",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (path / "main.cpp").write_text(
        r'''#include <cmath>
#include <iostream>
#include "pic/FieldSolver.hpp"
#include "pic/Grid.hpp"

int main() {
    pic::Grid grid(8, 1.0, pic::Boundary::Periodic);
    grid.rho()[1] = 0.5;
    grid.rho()[5] = -0.5;
    pic::FieldSolver solver;
    solver.solve(grid);
    const double sample = pic::interpolate_electric(grid, 0.25);
    if (!std::isfinite(sample)) {
        std::cerr << "non-finite interpolated field\n";
        return 1;
    }
    std::cout << "consumer smoke field=" << sample << "\n";
    return 0;
}
''',
        encoding="utf-8",
    )


def smoke_downstream_consumer(
    prefix: Path, work: Path, jobs: int
) -> None:
    consumer = work / "consumer"
    build = work / "consumer-build"
    write_consumer_project(consumer)
    run(
        [
            "cmake",
            "-S",
            str(consumer),
            "-B",
            str(build),
            f"-DCMAKE_PREFIX_PATH={prefix}",
            "-DCMAKE_BUILD_TYPE=Release",
        ]
    )
    run(["cmake", "--build", str(build), "--parallel", str(jobs)])
    run([str(build / "consumer")])


def smoke_install_tree(
    build_dir: Path, work: Path, jobs: int
) -> Path:
    prefix = work / "install-prefix"
    run(["cmake", "--install", str(build_dir), "--prefix", str(prefix)])
    require_file(prefix / "lib" / "cmake" / "AuroraPIC" / "AuroraPICConfig.cmake")
    require_file(prefix / "lib" / "cmake" / "AuroraPIC" / "AuroraPICConfigVersion.cmake")
    require_file(prefix / "lib" / "cmake" / "AuroraPIC" / "AuroraPICTargets.cmake")
    require_dir(prefix / "include" / "pic")
    require_dir(prefix / "share" / "aurorapic" / "examples")
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "compare_swarm.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "compare_hall.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "analyze_turner_balance.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "analyze_turner_density_blocks.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "prepare_turner_ensemble.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "attach_turner_ensemble_result.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "digitize_charoy_figure6.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "aggregate_hall_ensemble.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "analyze_hall_convergence.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "analyze_hall_flux_stationarity.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "compare_hall_flux_blocks.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "analyze_hall_pilot.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "lock_hall_source.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "normalize_hall_reference.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "examples" /
        "hall_landmark_case2.sources"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "preflight_hall.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "qualify_hall_runtime.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "prepare_hall_campaign.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "prepare_hall_ensemble.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "prepare_hall_convergence.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "prepare_hall_horizon_stage.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "analyze_hall_horizon_stage.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "run_swarm_campaign.py"
    )
    require_file(
        prefix / "share" / "aurorapic" / "scripts" /
        "validate_kinetic_benchmarks.py"
    )
    for script in (
        "export_phase_eedf.py", "compare_phase_eedf.py",
        "analyze_phase_eedf.py", "run_edupic_stage.py",
        "run_edupic_measurement_stage.py",
        "advance_edupic_measurement.py",
        "analyze_edupic_measurement_blocks.py",
        "evaluate_edupic_measurement_stationarity.py",
        "analyze_edupic_convergence.py",
        "advance_edupic_equilibration.py",
        "analyze_aurorapic_surface_flux.py",
        "analyze_aurorapic_surface_flux_timestep.py",
        "run_aurorapic_surface_flux_mesh.py",
        "analyze_aurorapic_surface_flux_mesh.py",
        "analyze_aurorapic_surface_flux_particle.py",
        "analyze_aurorapic_surface_flux_seed.py",
    ):
        require_file(prefix / "share" / "aurorapic" / "scripts" / script)
    require_file(
        prefix / "share" / "doc" / "AuroraPIC" /
        "phase-eedf-interchange.md"
    )
    require_file(
        prefix / "share" / "doc" / "AuroraPIC" /
        "phase-surface-flux.md"
    )
    require_file(
        prefix / "share" / "doc" / "AuroraPIC" /
        "kinetic-validation.md"
    )
    smoke_installed_cli(prefix, work)
    smoke_downstream_consumer(prefix, work, jobs)
    return prefix


def newest_package(build_dir: Path) -> Path:
    packages = sorted(build_dir.glob("AuroraPIC-*.tar.gz"), key=lambda p: p.stat().st_mtime)
    require(packages, f"no AuroraPIC TGZ package found under {build_dir}")
    return packages[-1]


def smoke_tgz_package(
    build_dir: Path, work: Path, jobs: int
) -> None:
    run([
        "cmake", "--build", str(build_dir), "--parallel", str(jobs),
        "--target", "package",
    ])
    package = newest_package(build_dir)
    require(tarfile.is_tarfile(package), f"package is not a valid tar archive: {package}")
    extract_dir = work / "package-extract"
    extract_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(package, "r:gz") as archive:
        archive.extractall(extract_dir)
    roots = [path for path in extract_dir.iterdir() if path.is_dir()]
    require(len(roots) == 1, f"expected one top-level package directory in {package}, found {len(roots)}")
    packaged_prefix = roots[0]
    require_file(packaged_prefix / "bin" / "aurorapic_cli")
    require_file(packaged_prefix / "bin" / "aurorapic_swarm")
    require_file(packaged_prefix / "lib" / "cmake" / "AuroraPIC" / "AuroraPICConfig.cmake")
    require_file(packaged_prefix / "lib" / "cmake" / "AuroraPIC" / "AuroraPICConfigVersion.cmake")
    require_file(packaged_prefix / "lib" / "cmake" / "AuroraPIC" / "AuroraPICTargets.cmake")
    require_dir(packaged_prefix / "include" / "pic")
    require_dir(packaged_prefix / "share" / "aurorapic" / "examples")
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "compare_swarm.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "compare_hall.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "digitize_charoy_figure6.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "aggregate_hall_ensemble.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "analyze_hall_convergence.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "analyze_hall_flux_stationarity.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "compare_hall_flux_blocks.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "analyze_hall_pilot.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "lock_hall_source.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "normalize_hall_reference.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "examples" /
        "hall_landmark_case2.sources"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "preflight_hall.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "qualify_hall_runtime.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "prepare_hall_campaign.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "prepare_hall_ensemble.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "prepare_hall_convergence.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "prepare_hall_horizon_stage.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "analyze_hall_horizon_stage.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "run_swarm_campaign.py"
    )
    require_file(
        packaged_prefix / "share" / "aurorapic" / "scripts" /
        "validate_kinetic_benchmarks.py"
    )
    require_file(
        packaged_prefix / "share" / "doc" / "AuroraPIC" /
        "kinetic-validation.md"
    )


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_dir", nargs="?", default="build", help="configured AuroraPIC build directory")
    parser.add_argument(
        "--jobs", type=int, default=1,
        help="maximum concurrent build jobs (default: 1)",
    )
    parser.add_argument("--keep-output", action="store_true", help="retain temporary install/package smoke outputs")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if args.jobs <= 0:
        print("install/package smoke failed: --jobs must be positive", file=sys.stderr)
        return 2
    build_dir = (ROOT / args.build_dir).resolve() if not Path(args.build_dir).is_absolute() else Path(args.build_dir)
    require_dir(build_dir)
    temp_root = Path(tempfile.mkdtemp(prefix="aurorapic_install_smoke_", dir=ROOT))
    try:
        smoke_install_tree(build_dir, temp_root, args.jobs)
        smoke_tgz_package(build_dir, temp_root, args.jobs)
    except (InstallSmokeError, subprocess.CalledProcessError) as exc:
        print(f"install/package smoke failed: {exc}", file=sys.stderr)
        print(f"retained smoke output: {temp_root}", file=sys.stderr)
        return 1
    if args.keep_output:
        print(f"install/package smoke passed; retained output: {temp_root}")
    else:
        shutil.rmtree(temp_root)
        print("install/package smoke passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
