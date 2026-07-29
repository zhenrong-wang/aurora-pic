#!/usr/bin/env python3
"""Validate AuroraPIC release-engineering artifacts.

This guard is intentionally static and lightweight so it can run in the smoke
suite before packaging. It catches drift in the CI matrix, CPack/install rules,
installable CMake package metadata, and the documented performance envelope
without depending on GitHub Actions being available locally.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CMAKE = ROOT / "CMakeLists.txt"
CMAKE_CONFIG_TEMPLATE = ROOT / "cmake" / "AuroraPICConfig.cmake.in"
CI = ROOT / ".github" / "workflows" / "ci.yml"
PERFORMANCE = ROOT / "docs" / "performance-envelope.md"
README = ROOT / "README.md"
ROADMAP = ROOT / "docs" / "multidimensional-roadmap.md"
VERIFY = ROOT / "scripts" / "verify.sh"
INSTALL_SMOKE = ROOT / "scripts" / "verify_install_package.py"
UNSTRUCTURED_BENCHMARK = ROOT / "scripts" / "benchmark_unstructured.py"


class ReleaseArtifactError(RuntimeError):
    """Raised when release-engineering artifacts drift."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ReleaseArtifactError(message)


def read(path: Path) -> str:
    require(path.is_file(), f"missing required file: {path.relative_to(ROOT)}")
    return path.read_text(encoding="utf-8")


def validate_cmake_packaging() -> None:
    cmake = read(CMAKE)
    for term in (
        "include(GNUInstallDirs)",
        "include(CMakePackageConfigHelpers)",
        "target_compile_features(aurorapic PUBLIC cxx_std_20)",
        "install(TARGETS aurorapic aurorapic_cli aurorapic_swarm",
        "install(DIRECTORY include/",
        "install(DIRECTORY examples/",
        "scripts/compare_swarm.py",
        "scripts/run_swarm_campaign.py",
        "scripts/validate_kinetic_benchmarks.py",
        "docs/kinetic-validation.md",
        "docs/performance-envelope.md",
        "docs/swarm-validation.md",
        "configure_package_config_file(",
        "AuroraPICConfigVersion.cmake",
        "set(CPACK_GENERATOR \"TGZ\")",
        "include(CPack)",
    ):
        require(term in cmake, f"CMake packaging must include {term!r}")

    config = read(CMAKE_CONFIG_TEMPLATE)
    for term in ("@PACKAGE_INIT@", "find_dependency(OpenMP)", "AuroraPICTargets.cmake"):
        require(term in config, f"CMake package config template must include {term!r}")


def validate_ci_matrix() -> None:
    ci = read(CI)
    for term in (
        "ubuntu-latest",
        "macos-latest",
        "compiler: gcc",
        "compiler: clang",
        "compiler: appleclang",
        "AURORA_ENABLE_OPENMP=${{ matrix.openmp }}",
        "python3 scripts/validate_release_artifacts.py",
        "python3 scripts/test_import_lxcat.py",
        "python3 scripts/test_swarm_cli.py build/aurorapic_swarm",
        "python3 scripts/test_compare_swarm.py",
        "python3 scripts/test_swarm_campaign.py",
        "python3 scripts/test_kinetic_benchmarks.py",
        "python3 scripts/validate_kinetic_benchmarks.py build/aurorapic_cli",
        "ctest --test-dir build --parallel 1 --output-on-failure",
        "python3 scripts/verify_examples.py build/aurorapic_cli",
        "python3 scripts/verify_install_package.py build --jobs 2",
        "actions/upload-artifact@v4",
    ):
        require(term in ci, f"CI workflow must include {term!r}")
    require("openmp: ON" in ci and "openmp: OFF" in ci, "CI matrix must cover OpenMP on and off")


def validate_install_smoke_script() -> None:
    script = read(INSTALL_SMOKE)
    for term in (
        '"cmake", "--install"',
        "find_package(AuroraPIC CONFIG REQUIRED)",
        "target_link_libraries(consumer PRIVATE AuroraPIC::aurorapic)",
        '"cmake", "--build"',
        '"--target",',
        '"package"',
        "AuroraPICConfig.cmake",
        "aurorapic_cli",
        "aurorapic_swarm",
        "compare_swarm.py",
        "run_swarm_campaign.py",
        "validate_kinetic_benchmarks.py",
        "kinetic-validation.md",
    ):
        require(term in script, f"install/package smoke script must include {term!r}")

    benchmark = read(UNSTRUCTURED_BENCHMARK)
    for term in (
        "time.perf_counter()",
        "statistics.median",
        "particle_seconds",
        "deposition_seconds",
        "field_solve_seconds",
        "location_cache_hits",
        "location_searches",
    ):
        require(term in benchmark, f"unstructured benchmark must include {term!r}")


def validate_performance_doc() -> None:
    doc = read(PERFORMANCE)
    for term in (
        "Verified smoke envelope",
        "Practical scaling expectations",
        "Before using larger runs",
        "Release-engineering envelope",
        "examples/two_stream.cfg",
        "examples/plasma_3d.cfg",
        "self-consistent electromagnetic field update",
        "not a general-purpose plasma production platform",
        "install/package smoke",
    ):
        require(term in doc, f"performance envelope must mention {term!r}")


def validate_cross_references() -> None:
    readme = read(README)
    roadmap = read(ROADMAP)
    verify = read(VERIFY)
    require("docs/performance-envelope.md" in readme, "README must link the performance envelope")
    require("docs/kinetic-validation.md" in readme, "README must link kinetic verification")
    require("CI workflow" in readme, "README must document CI workflow coverage")
    require("CPack" in readme, "README must document CPack packaging")
    require("find_package(AuroraPIC CONFIG REQUIRED)" in readme, "README must document downstream CMake package use")
    require("docs/performance-envelope.md" in roadmap, "roadmap must link the performance envelope")
    require("docs/kinetic-validation.md" in roadmap, "roadmap must link kinetic verification")
    require("CI matrix" in roadmap and "CPack" in roadmap, "roadmap must document CI matrix and CPack as M6 evidence")
    require("install/package smoke" in roadmap, "roadmap must document install/package smoke evidence")
    require(
        re.search(
            r"^python3\s+scripts/validate_release_artifacts\.py\s*$",
            verify, re.MULTILINE,
        ) is not None,
        "scripts/verify.sh must run scripts/validate_release_artifacts.py",
    )
    require(
        'python3 scripts/verify_install_package.py build --jobs "$BUILD_JOBS"'
        in verify,
        "scripts/verify.sh must run the resource-limited install/package smoke",
    )
    for term in (
        "AURORA_BUILD_JOBS",
        "AURORA_TEST_JOBS",
        "AURORA_OPENMP_THREADS",
        'cmake --build build --parallel "$BUILD_JOBS"',
        'ctest --test-dir build --parallel "$TEST_JOBS"',
    ):
        require(
            term in verify,
            f"scripts/verify.sh must enforce resource control {term!r}",
        )
    require(
        "python3 scripts/benchmark_unstructured.py build/aurorapic_cli --repeats 1"
        in verify,
        "scripts/verify.sh must run the unstructured benchmark smoke",
    )
    require(
        "python3 scripts/validate_kinetic_benchmarks.py build/aurorapic_cli"
        in verify,
        "scripts/verify.sh must run the quantitative kinetic benchmark",
    )


def main() -> int:
    try:
        validate_cmake_packaging()
        validate_ci_matrix()
        validate_install_smoke_script()
        validate_performance_doc()
        validate_cross_references()
    except ReleaseArtifactError as exc:
        print(f"release artifact validation failed: {exc}", file=sys.stderr)
        return 1
    print("release artifact validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
