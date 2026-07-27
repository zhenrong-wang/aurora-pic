#!/usr/bin/env python3
"""Validate the documented AuroraPIC production milestone contract.

This is intentionally lightweight: the milestone ladder is a project-management
artifact, but it now gates the smoke suite so milestone IDs, statuses, README
linkage, and verify-suite integration do not silently drift.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ROADMAP = ROOT / "docs" / "multidimensional-roadmap.md"
README = ROOT / "README.md"
VERIFY = ROOT / "scripts" / "verify.sh"

EXPECTED_MILESTONES = {
    "M0": "Regression-preserving multidimensional PIC core",
    "M1": "Validation and benchmark suite",
    "M2": "Geometry and mesh import workflow",
    "M3": "Scalable data and restart formats",
    "M4": "Runtime scaling backend",
    "M5": "Higher-fidelity physics",
    "M6": "Release engineering and operability",
}
EXPECTED_STATUSES = {
    "M0": "Complete",
    "M1": "Complete",
    "M2": "Complete",
    "M3": "Complete",
    "M4": "Complete",
    "M5": "Complete",
    "M6": "Current baseline",
}
REQUIRED_EVIDENCE_TERMS = {
    "M0": ("CTest", "CLI examples", "pusher validation"),
    "M1": ("benchmarks", "documented tolerances"),
    "M2": ("Gmsh", "boundary labels"),
    "M3": ("VTK XML", "openPMD/HDF5", "compatibility tests"),
    "M4": ("OpenMP/MPI/GPU", "scaling smoke tests"),
    "M5": ("uniform-B", "Boris", "CTest", "CLI examples", "electromagnetic", "collision models"),
    "M6": ("config_version", "configuration compatibility", "clearer failure diagnostics", "CTest", "CI matrix"),
}


@dataclass(frozen=True)
class MilestoneRow:
    identifier: str
    status: str
    milestone: str
    evidence: str


class MilestoneValidationError(RuntimeError):
    """Raised when the milestone documentation contract is broken."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise MilestoneValidationError(message)


def section_text(markdown: str, heading: str) -> str:
    marker = f"## {heading}"
    start = markdown.find(marker)
    require(start >= 0, f"missing roadmap heading: {marker}")
    following = markdown.find("\n## ", start + len(marker))
    return markdown[start:] if following < 0 else markdown[start:following]


def parse_milestone_rows(section: str) -> list[MilestoneRow]:
    rows: list[MilestoneRow] = []
    for raw_line in section.splitlines():
        line = raw_line.strip()
        if not line.startswith("| M"):
            continue
        columns = [column.strip() for column in line.strip("|").split("|")]
        require(len(columns) == 4, f"milestone row should have 4 columns: {raw_line}")
        rows.append(MilestoneRow(*columns))
    return rows


def validate_roadmap() -> None:
    markdown = ROADMAP.read_text(encoding="utf-8")
    section = section_text(markdown, "Production-readiness milestone ladder")
    rows = parse_milestone_rows(section)
    seen = [row.identifier for row in rows]
    expected_ids = list(EXPECTED_MILESTONES)
    require(seen == expected_ids, f"milestone IDs changed: expected {expected_ids}, found {seen}")

    for row in rows:
        expected_title = EXPECTED_MILESTONES[row.identifier]
        expected_status = EXPECTED_STATUSES[row.identifier]
        require(row.milestone == expected_title, f"{row.identifier} title drifted: {row.milestone!r}")
        require(row.status == expected_status, f"{row.identifier} status drifted: {row.status!r}")
        require(row.evidence.endswith("."), f"{row.identifier} evidence should be a sentence ending with '.'")
        for term in REQUIRED_EVIDENCE_TERMS[row.identifier]:
            require(term in row.evidence, f"{row.identifier} evidence must mention {term!r}")
    require("### Immediate coding target" in section, "roadmap must identify the immediate coding target")
    require("scripts/validate_milestones.py" in section, "roadmap must reference this validation script")
    require("M2" in section and "Gmsh v2 ASCII" in section and "boundary labels" in section,
            "roadmap must keep the M2 importer baseline visible")
    require("M3" in section and "VTK XML" in section,
            "roadmap must keep the M3 output/restart target visible")
    require("M4" in section and "OpenMP/MPI/GPU" in section and "scaling smoke tests" in section,
            "roadmap must keep the M4 runtime-scaling baseline visible")
    require("M6" in section and "config_version" in section and "configuration compatibility" in section,
            "roadmap must keep the M6 configuration-compatibility baseline visible")


def validate_readme() -> None:
    readme = README.read_text(encoding="utf-8")
    require("## Production milestone baseline" in readme, "README missing production milestone baseline section")
    require("docs/multidimensional-roadmap.md#production-readiness-milestone-ladder" in readme,
            "README must link to the milestone ladder anchor")
    require("scripts/validate_milestones.py" in readme,
            "README must document that milestone validation is in the smoke suite")


def validate_verify_script() -> None:
    verify = VERIFY.read_text(encoding="utf-8")
    pattern = re.compile(r"^python3\s+scripts/validate_milestones\.py\s*$", re.MULTILINE)
    require(pattern.search(verify) is not None, "scripts/verify.sh must run scripts/validate_milestones.py")


def main() -> int:
    try:
        validate_roadmap()
        validate_readme()
        validate_verify_script()
    except MilestoneValidationError as exc:
        print(f"milestone validation failed: {exc}", file=sys.stderr)
        return 1
    print("milestone validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
