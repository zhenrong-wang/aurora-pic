#!/usr/bin/env python3

import math

from audit_edupic_mover_contract import AuditError, phase_metrics


def main() -> int:
    metrics = phase_metrics(4000, 250.0)
    assert math.isclose(metrics["phase_advance_rad"], math.pi / 2000.0)
    assert math.isclose(metrics["rf_cycle_fraction"], 0.00025)
    assert math.isclose(
        metrics["maximum_normalized_voltage_difference"],
        2.0 * math.sin(math.pi / 4000.0))
    assert metrics["maximum_voltage_difference_V"] < 0.4
    try:
        phase_metrics(0, 250.0)
    except AuditError:
        pass
    else:
        raise AssertionError("zero steps per cycle must be rejected")
    print("eduPIC mover-contract audit tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
