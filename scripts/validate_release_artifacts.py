#!/usr/bin/env python3
"""Validate AuroraPIC release-engineering artifacts.

This guard is intentionally static and lightweight so it can run in the smoke
suite before packaging. It catches drift in the CI matrix, CPack install rules,
and the documented performance envelope without depending on GitHub Actions
being available locally.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CMAKE = ROOT / "CMakeLists.txt"
CI = ROOT / ".github" / "workflows" / "ci.yml"
PERFORMANCE = ROOT / "docs" / "performance-envelope.md"
README = ROOT / "README.md"
ROADMAP = ROOT / "docs" / "multidimensional-roadmap.md"
VERIFY = ROOT / "scripts" / "verify.sh"


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
        "install(TARGETS aurorapic aurorapic_cli",
        "install(DIRECTORY include/",
        "install(DIRECTORY examples/",
        "docs/performance-envelope.md",
        "set(CPACK_GENERATOR \"TGZ\")",
        "include(CPack)",
    ):
        require(term in cmake, f"CMake packaging must include {term!r}")


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
        "ctest --test-dir build --output-on-failure",
        "python3 scripts/verify_examples.py build/aurorapic_cli",
        "cmake --build build --target package",
        "actions/upload-artifact@v4",
    ):
        require(term in ci, f"CI workflow must include {term!r}")
    require("openmp: ON" in ci and "openmp: OFF" in ci, "CI matrix must cover OpenMP on and off")


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
    ):
        require(term in doc, f"performance envelope must mention {term!r}")


def validate_cross_references() -> None:
    readme = read(README)
    roadmap = read(ROADMAP)
    verify = read(VERIFY)
    require("docs/performance-envelope.md" in readme, "README must link the performance envelope")
    require("CI workflow" in readme, "README must document CI workflow coverage")
    require("CPack" in readme, "README must document CPack packaging")
    require("docs/performance-envelope.md" in roadmap, "roadmap must link the performance envelope")
    require("CI matrix" in roadmap and "CPack" in roadmap, "roadmap must document CI matrix and CPack as M6 evidence")
    pattern = re.compile(r"^python3\s+scripts/validate_release_artifacts\.py\s*$", re.MULTILINE)
    require(pattern.search(verify) is not None, "scripts/verify.sh must run scripts/validate_release_artifacts.py")


def main() -> int:
    try:
        validate_cmake_packaging()
        validate_ci_matrix()
        validate_performance_doc()
        validate_cross_references()
    except ReleaseArtifactError as exc:
        print(f"release artifact validation failed: {exc}", file=sys.stderr)
        return 1
    print("release artifact validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
