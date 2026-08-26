#!/usr/bin/env python3
"""Add passive ten-phase-per-cycle grid snapshots to pinned native eduPIC."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from instrument_edupic_mover_decomposition import instrument as instrument_mover
from instrument_edupic_phase_eedf import EXPECTED_SOURCE_SHA256
from instrument_edupic_promotion_band_work import sha256
from instrument_edupic_surface_flux import replace_once


def instrument(source: str) -> str:
    source = instrument_mover(source)
    helper = r'''
void save_phase_snapshot(int time_step, xvector rho){
    const bool first_snapshot =
        cycle == cycles_done + 1 && time_step == 400;
    FILE *output = fopen(
        "edupic_phase_snapshots.csv", first_snapshot ? "w" : "a");
    if (first_snapshot) {
        fprintf(output,
            "measurement_cycle,phase_index,phase_fraction,node,x_m,charge_density_C_m3,potential_V,electric_field_V_m\n");
    }
    const int measurement_cycle = cycle - cycles_done;
    const int phase_index = time_step / 400;
    const double phase_fraction =
        static_cast<double>(time_step) / static_cast<double>(N_T);
    for (int node=0; node<N_G; ++node) {
        fprintf(output, "%d,%d,%.17g,%d,%.17g,%.17g,%.17g,%.17g\n",
            measurement_cycle, phase_index, phase_fraction, node,
            static_cast<double>(node) * DX, rho[node], pot[node],
            efield[node]);
    }
    fclose(output);
}

'''
    source = replace_once(
        source,
        "void do_one_cycle (void){\n",
        helper + "void do_one_cycle (void){\n",
        "phase snapshot helper",
    )
    source = replace_once(
        source,
        "        solve_Poisson(rho,Time);                                                 // compute potential and electric field\n",
        "        solve_Poisson(rho,Time);                                                 // compute potential and electric field\n"
        "        if (measurement_mode && ((t + 1) % 400) == 0)\n"
        "            save_phase_snapshot(t + 1, rho);\n",
        "phase snapshot hook",
    )
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
        "transform": "passive_edupic_cycle_phase_grid_snapshots_v1",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "composes_mover_decomposition": True,
        "snapshots_per_cycle": 10,
        "snapshot_steps_within_cycle": list(range(400, 4001, 400)),
        "snapshot_phase_fractions": [value / 10 for value in range(1, 11)],
        "snapshot_timing": "after density deposition and Poisson solve, before particle push",
        "columns": ["charge_density_C_m3", "potential_V", "electric_field_V_m"],
        "random_draws_added": 0,
        "particle_state_mutation_added": False,
    }
    if args.audit:
        args.audit.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n",
                              encoding="utf-8", newline="\n")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
