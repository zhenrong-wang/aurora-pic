#!/usr/bin/env python3

from analyze_edupic_endpoint_gauss_control import decision


def main() -> int:
    assert decision([0.80, 0.85], True)["deficit_persists"]
    assert decision([0.97, 1.02], True)[
        "endpoint_correction_explains_deficit"]
    assert decision([0.92, 0.94], True)["partial_endpoint_effect"]
    blocked = decision([0.80, 0.85], False)
    assert not blocked["interpretation_allowed"]
    assert blocked["directional_persistence_signal_without_interpretation"]
    print("endpoint Gauss-control analyzer tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
