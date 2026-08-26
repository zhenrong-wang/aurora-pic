#!/usr/bin/env python3
"""Focused test for passive common-state population instrumentation."""

from instrument_edupic_common_state_population import instrument


def main() -> None:
    source = """bool common_state_trace_step(int time_step) {
}
        if (common_state_trace_step(t + 1))
            save_common_state_trace(t + 1, rho);
"""
    result = instrument(source)
    assert "save_common_state_population(t + 1)" in result
    assert "particle" not in result.lower()
    print("common-state population instrumenter tests passed")


if __name__ == "__main__":
    main()
