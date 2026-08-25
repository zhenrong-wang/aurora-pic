#!/usr/bin/env python3

import math
from pathlib import Path
import tempfile

from audit_edupic_frozen_collision_kernel import (
    opportunity_excess, scan_lookup_kernel,
)


def main() -> int:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        paths = []
        for channel, values in enumerate((
                (1.0, 2.0, 3.0, 4.0),
                (0.0, 1.0, 1.0, 2.0))):
            path = root / f"channel_{channel}.dat"
            path.write_text("# fixture\n" + "".join(
                f"{index * .001:.3f} {value}\n"
                for index, value in enumerate(values)), encoding="utf-8")
            paths.append(path)
        result = scan_lookup_kernel(
            paths, 1.0, 1.0, 1.0e-12, energy_limit_ev=.003)
        assert result["upper_half_bin_queries"] == 3
        assert result["maximum_absolute_event_probability_difference_per_step"] > 0
        assert result["maximum_channel_probability_total_variation"] > 0
    expected = (.02 - (1.0 - math.exp(-.02))) / .02
    assert math.isclose(opportunity_excess(2.0, .01), expected)
    print("frozen collision-kernel audit tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
