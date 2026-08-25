#!/usr/bin/env python3

from array import array
import math

from fold_edupic_frozen_collision_states import (
    accumulate_event, lookup_rates,
)


def main() -> int:
    tables = [array("d", [1.0, 2.0, 4.0]),
              array("d", [0.0, 1.0, 3.0])]
    native = lookup_rates(.00175, 2.0, tables, (0.0, 0.0), True)
    aurora = lookup_rates(.00175, 2.0, tables, (0.0, 0.0), False)
    assert native[0] == 2.0 * aurora[0]
    assert native[1] == 3.0 * aurora[1]
    channels = [0.0, 0.0]
    probability, mean = accumulate_event([2.0, 3.0], .1, channels)
    assert math.isclose(probability, 1.0 - math.exp(-.5))
    assert math.isclose(mean, .5)
    assert math.isclose(sum(channels), probability)
    print("frozen collision-state fold tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
