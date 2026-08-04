#!/usr/bin/env python3
"""Synthetic regression for the eduPIC population-stationarity gate."""

from __future__ import annotations

from pathlib import Path
import tempfile

from analyze_edupic_convergence import analyze, read_history


ROOT = Path(__file__).resolve().parents[1]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="aurorapic_edupic_convergence_",
                                     dir=ROOT / "tmp") as tmp:
        work = Path(tmp)
        stable = work / "stable.dat"
        stable.write_text("".join(
            f"{cycle} {50000 + cycle % 2} {60000 - cycle % 2}\n"
            for cycle in range(1401, 1501)), encoding="utf-8")
        stable_report = analyze(read_history(stable), 100, 25, 1e-4, 0.02)
        require(stable_report["stationary"] and stable_report["window"]["eligible"],
                "stable synthetic history failed stationarity")

        rising = work / "rising.dat"
        rising.write_text("".join(
            f"{cycle} {1000 + 20 * cycle} {2000 + 30 * cycle}\n"
            for cycle in range(1, 101)), encoding="utf-8")
        rising_report = analyze(read_history(rising), 100, 25, 1e-4, 0.02,
                                minimum_cycle=100)
        require(not rising_report["stationary"] and not rising_report["criteria"]
                ["total_particles"]["relative_slope_passes"] and
                rising_report["provisional_recent"]["samples"] == 25 and
                rising_report["provisional_recent"]["metrics"]
                ["total_particles"]["relative_slope_per_cycle"] > 0.0,
                "rising synthetic history passed stationarity")

        short_report = analyze(read_history(rising)[:8], 100, 25, 1e-4, 0.02)
        require(not short_report["stationary"] and not short_report["window"]["eligible"],
                "short history became stationarity-eligible")
    print("eduPIC convergence analysis regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
