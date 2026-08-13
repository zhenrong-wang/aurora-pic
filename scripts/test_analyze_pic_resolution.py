#!/usr/bin/env python3
"""Regression tests for explicit PIC resolution input handling."""

import math
import tempfile
from pathlib import Path

from analyze_pic_resolution import positive, sections


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="aurorapic_resolution_") as temporary:
        path = Path(temporary) / "input.cfg"
        path.write_text(
            "units = si\n# comment\nnx = 5\n[species.electrons]\n"
            "timestep_multiplier = 1\n", encoding="utf-8")
        parsed = sections(path)
        assert parsed["global"]["nx"] == "5"
        assert parsed["species.electrons"]["timestep_multiplier"] == "1"
    assert math.isclose(positive("1e-3", "test"), 1e-3)
    try:
        positive("0", "test")
    except ValueError:
        pass
    else:
        raise AssertionError("non-positive input was accepted")


if __name__ == "__main__":
    main()
