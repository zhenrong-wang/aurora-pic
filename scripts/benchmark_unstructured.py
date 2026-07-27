#!/usr/bin/env python3
"""Run a reproducible imported-geometry timing benchmark."""

import argparse
import csv
import pathlib
import re
import statistics
import subprocess
import tempfile
import time


def set_global(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(key)}\s*=).*$")
    if pattern.search(text):
        return pattern.sub(lambda match: f"{match.group(1)} {value}", text, count=1)
    section = re.search(r"(?m)^\s*\[", text)
    insertion = section.start() if section else len(text)
    return text[:insertion] + f"{key} = {value}\n" + text[insertion:]


def prepare_config(source: pathlib.Path, destination: pathlib.Path,
                   output_dir: pathlib.Path) -> None:
    text = source.read_text(encoding="utf-8")
    mesh_match = re.search(r"(?m)^\s*mesh_file\s*=\s*(.*?)\s*$", text)
    if not mesh_match:
        raise RuntimeError("benchmark config is missing mesh_file")
    mesh_path = pathlib.Path(mesh_match.group(1))
    if not mesh_path.is_absolute():
        mesh_path = (source.parent / mesh_path).resolve()
    text = set_global(text, "mesh_file", str(mesh_path))
    text = set_global(text, "output_dir", str(output_dir))
    text = set_global(text, "vtk_output", "false")
    text = set_global(text, "particle_output", "false")
    text = set_global(text, "checkpoint_output", "false")
    destination.write_text(text, encoding="utf-8")


def final_timing(csv_path: pathlib.Path) -> dict[str, float]:
    with csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError("benchmark produced no diagnostic samples")
    final = rows[-1]
    return {
        key: float(final[key])
        for key in (
            "particle_seconds",
            "deposition_seconds",
            "field_solve_seconds",
            "location_cache_hits",
            "location_searches",
        )
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark AuroraPIC's imported 2D runtime")
    parser.add_argument("executable", type=pathlib.Path)
    parser.add_argument(
        "config", type=pathlib.Path,
        nargs="?", default=pathlib.Path("examples/imported_plasma_2d.cfg"))
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats <= 0:
        parser.error("--repeats must be positive")

    executable = args.executable.resolve()
    config = args.config.resolve()
    wall_times: list[float] = []
    phase_times: dict[str, list[float]] = {
        "particle_seconds": [],
        "deposition_seconds": [],
        "field_solve_seconds": [],
        "location_cache_hits": [],
        "location_searches": [],
    }
    with tempfile.TemporaryDirectory(prefix="aurorapic_benchmark_") as temporary:
        root = pathlib.Path(temporary)
        for repeat in range(args.repeats):
            output_dir = root / f"output_{repeat}"
            benchmark_config = root / f"benchmark_{repeat}.cfg"
            prepare_config(config, benchmark_config, output_dir)
            start = time.perf_counter()
            subprocess.run(
                [str(executable), str(benchmark_config)],
                check=True, stdout=subprocess.DEVNULL)
            wall_times.append(time.perf_counter() - start)
            for key, value in final_timing(output_dir / "scalars.csv").items():
                phase_times[key].append(value)

    print(f"repeats={args.repeats}")
    print(f"wall_seconds_median={statistics.median(wall_times):.9g}")
    for key, values in phase_times.items():
        print(f"{key}_median={statistics.median(values):.9g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
