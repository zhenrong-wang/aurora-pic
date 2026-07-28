#!/usr/bin/env python3
"""Regenerate AuroraPIC's Gmsh-backed real-case integration mesh."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "examples" / "biased_probe_2d.geo"
MESH = ROOT / "examples" / "biased_probe_2d.msh"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gmsh",
        type=Path,
        help="Gmsh executable; defaults to the gmsh found on PATH",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="destination .msh path; defaults beside the geometry when writable",
    )
    args = parser.parse_args()

    executable = args.gmsh or (
        Path(found) if (found := shutil.which("gmsh")) else None
    )
    if executable is None:
        raise SystemExit(
            "gmsh was not found; install Gmsh 4.x or pass --gmsh /path/to/gmsh"
        )
    version = subprocess.run(
        [str(executable), "--version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    default_output = (
        MESH
        if os.access(GEOMETRY.parent, os.W_OK)
        else Path.cwd() / MESH.name
    )
    output = (args.output or default_output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(executable),
            "-2",
            "-format",
            "msh2",
            "-o",
            str(output),
            str(GEOMETRY),
        ],
        cwd=ROOT,
        check=True,
    )
    header = output.read_text(encoding="utf-8").splitlines()[:3]
    if header != ["$MeshFormat", "2.2 0 8", "$EndMeshFormat"]:
        raise SystemExit(f"Gmsh {version} did not produce a v2 ASCII mesh")
    print(f"generated {output} with Gmsh {version}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
