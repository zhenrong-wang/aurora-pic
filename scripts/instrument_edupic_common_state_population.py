#!/usr/bin/env python3
"""Add passive population rows to the locked eduPIC common-state trace."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from instrument_edupic_promotion_band_work import sha256
from instrument_edupic_surface_flux import replace_once


EXPECTED_TRACE_SOURCE_SHA256 = (
    "b423fe256ab2d4757041886737215f56f45c9d49d61ebf9913b32e8226bd5c13")


def instrument(source: str) -> str:
    helper = r'''
void save_common_state_population(int time_step) {
    FILE *output = fopen(
        "edupic_common_state_population.csv", time_step == 1 ? "w" : "a");
    if (time_step == 1)
        fprintf(output, "pre_push_step,electrons,ions\n");
    fprintf(output, "%d,%d,%d\n", time_step, N_e, N_i);
    fclose(output);
}

'''
    source = replace_once(
        source, "bool common_state_trace_step(int time_step) {\n",
        helper + "bool common_state_trace_step(int time_step) {\n",
        "population helper")
    source = replace_once(
        source,
        "        if (common_state_trace_step(t + 1))\n"
        "            save_common_state_trace(t + 1, rho);\n",
        "        if (common_state_trace_step(t + 1)) {\n"
        "            save_common_state_trace(t + 1, rho);\n"
        "            save_common_state_population(t + 1);\n"
        "        }\n",
        "population hook")
    return source


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    if sha256(args.source) != EXPECTED_TRACE_SOURCE_SHA256:
        raise RuntimeError("locked common-state trace source differs")
    args.output.write_text(instrument(args.source.read_text(encoding="utf-8")),
                           encoding="utf-8", newline="\n")
    audit = {
        "schema_version": 1,
        "transform": "passive_common_state_population_extension_v1",
        "source_sha256": EXPECTED_TRACE_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "new_columns": ["pre_push_step", "electrons", "ions"],
        "random_draws_added": 0,
        "particle_state_mutation_added": False,
        "passivity_requirement": "grid trace and final checkpoint byte-identical to locked primary trace run"
    }
    if args.audit:
        args.audit.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n",
                              encoding="utf-8", newline="\n")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
