#!/usr/bin/env python3
"""Dependency-free unit tests for the kinetic benchmark analyzer."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "validate_kinetic_benchmarks.py"
SPEC = importlib.util.spec_from_file_location(
    "validate_kinetic_benchmarks", MODULE_PATH
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load kinetic benchmark module")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def require_near(actual: float, expected: float, tolerance: float, message: str) -> None:
    if abs(actual - expected) > tolerance:
        raise RuntimeError(
            f"{message}: expected {expected}, got {actual}"
        )


def main() -> int:
    times = [0.025 * index for index in range(481)]
    damping = -0.1533
    frequency = 1.4156
    amplitudes = [
        0.02 * math.exp(damping * time) *
        abs(math.cos(frequency * time))
        for time in times
    ]
    fit = MODULE.analyze_damped_mode(times, amplitudes)
    require_near(
        float(fit["damping_rate"]), damping, 0.002,
        "synthetic damping-rate fit changed",
    )
    require_near(
        float(fit["angular_frequency"]), frequency, 0.01,
        "synthetic frequency fit changed",
    )
    growth = 0.2258
    growing_amplitudes = [
        0.001 * math.exp(growth * time)
        for time in times
    ]
    growth_fit = MODULE.analyze_exponential_growth(
        times, growing_amplitudes, 8.0, 12.0
    )
    require_near(
        float(growth_fit["growth_rate"]), growth, 1e-12,
        "synthetic growth-rate fit changed",
    )
    require_near(
        float(growth_fit["r_squared"]), 1.0, 1e-12,
        "synthetic growth fit quality changed",
    )
    print("kinetic benchmark analyzer tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
