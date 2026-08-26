#!/usr/bin/env python3
"""Create a collision-free, logarithmically sampled eduPIC grid trace."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from instrument_edupic_phase_eedf import EXPECTED_SOURCE_SHA256
from instrument_edupic_promotion_band_work import sha256
from instrument_edupic_surface_flux import replace_once


TRACE_STEPS = (1, 2, 3, 6, 11, 21, 51, 101, 201, 401,
               801, 1601, 2401, 3201, 4000)


def instrument(source: str) -> str:
    selected = " || ".join(f"time_step == {step}" for step in TRACE_STEPS)
    helper = rf'''
bool common_state_trace_step(int time_step) {{
    return {selected};
}}

void save_common_state_trace(int time_step, xvector rho) {{
    FILE *output = fopen(
        "edupic_common_state_trace.csv", time_step == 1 ? "w" : "a");
    if (time_step == 1) {{
        fprintf(output,
            "pre_push_step,node,x_m,charge_density_C_m3,potential_V,electric_field_V_m\n");
    }}
    for (int node=0; node<N_G; ++node) {{
        fprintf(output, "%d,%d,%.17g,%.17g,%.17g,%.17g\n",
            time_step, node, static_cast<double>(node) * DX,
            rho[node], pot[node], efield[node]);
    }}
    fclose(output);
}}

'''
    source = replace_once(
        source, "void do_one_cycle (void){\n",
        helper + "void do_one_cycle (void){\n", "trace helper")
    solve = (
        "        solve_Poisson(rho,Time);                                                 // compute potential and electric field\n")
    source = replace_once(
        source, solve,
        solve + "        if (common_state_trace_step(t + 1))\n"
        "            save_common_state_trace(t + 1, rho);\n",
        "pre-push trace hook")
    source = replace_once(
        source,
        "        for (k=0; k<N_e; k++){                              // checking for occurrence of a collision for all electrons in every time step\n",
        "        if (false) for (k=0; k<N_e; k++){                   // collision-free common-state trace\n",
        "electron collision disable")
    source = replace_once(
        source,
        "            for (k=0; k<N_i; k++){\n                vx_a = RMB(MTgen);",
        "            if (false) for (k=0; k<N_i; k++){\n                vx_a = RMB(MTgen);",
        "ion collision disable")
    return source


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    if sha256(args.source) != EXPECTED_SOURCE_SHA256:
        raise RuntimeError("eduPIC source differs from pinned C implementation")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(instrument(args.source.read_text(encoding="utf-8")),
                           encoding="utf-8", newline="\n")
    audit = {
        "schema_version": 1,
        "transform": "edupic_collision_free_common_state_grid_trace_v1",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "trace_pre_push_steps": list(TRACE_STEPS),
        "matching_aurorapic_post_step_horizons": [step - 1 for step in TRACE_STEPS],
        "snapshot_timing": "after deposition and Poisson solve, before particle push",
        "electron_collisions_disabled": True,
        "ion_collisions_disabled": True,
        "random_draws_added": 0,
    }
    if args.audit:
        args.audit.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n",
                              encoding="utf-8", newline="\n")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
