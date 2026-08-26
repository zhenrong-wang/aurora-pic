#!/usr/bin/env python3
"""Focused tests for the Turner subcycle-policy invariance runner."""

from pathlib import Path

from run_turner_subcycle_policy_invariance import deck


def main() -> None:
    rule = {"paired_contract": {"seed": 13507, "start_step": 770401,
                                 "end_step": 770800}}
    current = deck(rule, Path("state.apc"), Path("electron.gas"),
                   Path("ion.gas"), Path("a"), "current_position")
    held = deck(rule, Path("state.apc"), Path("electron.gas"),
                Path("ion.gas"), Path("b"), "pre_push_held")
    assert "subcycle_charge_deposition = current_position" in current
    assert "subcycle_charge_deposition = pre_push_held" in held
    assert "timestep_multiplier" not in current
    normalized_current = current.replace("output_dir = a", "output_dir = X").replace(
        "subcycle_charge_deposition = current_position",
        "subcycle_charge_deposition = POLICY")
    normalized_held = held.replace("output_dir = b", "output_dir = X").replace(
        "subcycle_charge_deposition = pre_push_held",
        "subcycle_charge_deposition = POLICY")
    assert normalized_current == normalized_held
    print("Turner subcycle-policy invariance runner tests passed")


if __name__ == "__main__":
    main()
