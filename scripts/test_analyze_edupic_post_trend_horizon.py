#!/usr/bin/env python3
"""Regression for post-trend horizon forecasting math."""

from analyze_edupic_post_trend_horizon import linear_forecast


def main() -> int:
    assert linear_forecast([1.0, 2.0, 3.0], [3.0, 5.0, 7.0], 4.0) == 9.0
    assert linear_forecast([1.0, 3.0], [2.0, 2.0], 10.0) == 2.0
    print("eduPIC post-trend horizon audit regression passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
