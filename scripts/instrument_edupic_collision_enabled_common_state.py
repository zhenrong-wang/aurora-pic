#!/usr/bin/env python3
"""Instrument pinned eduPIC for the collision-enabled common-state pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


EXPECTED_SOURCE_SHA256 = (
    "7c7679c0f0c98844940ea911bbb7581ec33f818e8d14427c9837ffdcf1ecea41")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_once(source: str, old: str, new: str, label: str) -> str:
    if source.count(old) != 1:
        raise RuntimeError(f"expected exactly one {label}, found {source.count(old)}")
    return source.replace(old, new, 1)


def instrument(source: str) -> str:
    source = replace_once(
        source,
        "Ullong   N_i_coll                   = 0;                     // counter for ion collisions\n",
        "Ullong   N_i_coll                   = 0;                     // counter for ion collisions\n"
        "Ullong   N_e_elastic                = 0;\n"
        "Ullong   N_e_excitation             = 0;\n"
        "Ullong   N_e_ionization             = 0;\n"
        "Ullong   N_i_isotropic              = 0;\n"
        "Ullong   N_i_backward               = 0;\n",
        "collision counter declaration")
    source = replace_once(
        source,
        "    if (rnd < (t0/t2)){                              // elastic scattering\n",
        "    if (rnd < (t0/t2)){                              // elastic scattering\n"
        "        N_e_elastic++;\n",
        "electron elastic branch")
    source = replace_once(
        source,
        "    } else if (rnd < (t1/t2)){                       // excitation\n",
        "    } else if (rnd < (t1/t2)){                       // excitation\n"
        "        N_e_excitation++;\n",
        "electron excitation branch")
    source = replace_once(
        source,
        "    } else {                                         // ionization\n",
        "    } else {                                         // ionization\n"
        "        N_e_ionization++;\n",
        "electron ionization branch")
    source = replace_once(
        source,
        "    if  (rnd < (t1 /t2)){                        // isotropic scattering\n",
        "    if  (rnd < (t1 /t2)){                        // isotropic scattering\n"
        "        N_i_isotropic++;\n",
        "ion isotropic branch")
    source = replace_once(
        source,
        "    } else {                                     // backward scattering\n",
        "    } else {                                     // backward scattering\n"
        "        N_i_backward++;\n",
        "ion backward branch")

    helper = r'''
void save_collision_enabled_common_state_endpoint(int pre_push_step) {
    FILE *field = fopen("edupic_collision_endpoint_field.csv", "w");
    fprintf(field, "pre_push_step,node,x_m,electric_field_V_m\n");
    for (int node=0; node<N_G; ++node) {
        fprintf(field, "%d,%d,%.17g,%.17g\n", pre_push_step, node,
                static_cast<double>(node) * DX, efield[node]);
    }
    fclose(field);
    FILE *metrics = fopen("edupic_collision_endpoint_metrics.csv", "w");
    fprintf(metrics,
        "pre_push_step,electrons,ions,electron_elastic,electron_excitation,"
        "electron_ionization,ion_isotropic,ion_backward,"
        "electron_absorbed_left,electron_absorbed_right,"
        "ion_absorbed_left,ion_absorbed_right\n");
    fprintf(metrics, "%d,%d,%d,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu\n",
        pre_push_step, N_e, N_i, N_e_elastic, N_e_excitation,
        N_e_ionization, N_i_isotropic, N_i_backward,
        N_e_abs_pow, N_e_abs_gnd, N_i_abs_pow, N_i_abs_gnd);
    fclose(metrics);
}

'''
    source = replace_once(
        source, "void do_one_cycle (void){\n",
        helper + "void do_one_cycle (void){\n", "cycle function")
    solve = (
        "        solve_Poisson(rho,Time);                                                 // compute potential and electric field\n")
    source = replace_once(
        source, solve,
        solve + "        if ((t + 1) == 4000) {\n"
        "            save_collision_enabled_common_state_endpoint(t + 1);\n"
        "            return;\n"
        "        }\n",
        "endpoint hook")
    load = (
        "        load_particle_data();                             // read previous configuration from file\n")
    source = replace_once(
        source, load,
        load + "        if (argc < 4) {\n"
        "            printf(\">> eduPIC: error = collision pilot needs RNG seed argv[3]\\n\");\n"
        "            return 2;\n"
        "        }\n"
        "        MTgen.seed(static_cast<unsigned long>(strtoul(argv[3], NULL, 10)));\n",
        "post-load seed hook")
    return source


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    if sha256(args.source) != EXPECTED_SOURCE_SHA256:
        raise RuntimeError("eduPIC source differs from locked input")
    transformed = instrument(args.source.read_text(encoding="utf-8"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(transformed, encoding="utf-8", newline="\n")
    audit = {
        "schema_version": 1,
        "transform": "edupic_collision_enabled_common_state_endpoint_v1",
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "instrumented_source_sha256": sha256(args.output),
        "endpoint_pre_push_step": 4000,
        "rng_seed_timing": "immediately_after_checkpoint_load",
        "channel_counters_added": [
            "electron_elastic", "electron_excitation", "electron_ionization",
            "ion_isotropic", "ion_backward"],
        "random_draws_added": 0,
    }
    if args.audit:
        args.audit.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n",
                              encoding="utf-8", newline="\n")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
