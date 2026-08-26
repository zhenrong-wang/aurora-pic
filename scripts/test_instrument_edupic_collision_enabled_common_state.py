#!/usr/bin/env python3
"""Focused tests for collision-enabled common-state instrumentation."""

from instrument_edupic_collision_enabled_common_state import instrument


def main() -> None:
    source = """Ullong   N_i_coll                   = 0;                     // counter for ion collisions
    if (rnd < (t0/t2)){                              // elastic scattering
    } else if (rnd < (t1/t2)){                       // excitation
    } else {                                         // ionization
    if  (rnd < (t1 /t2)){                        // isotropic scattering
    } else {                                     // backward scattering
void do_one_cycle (void){
        solve_Poisson(rho,Time);                                                 // compute potential and electric field
        load_particle_data();                             // read previous configuration from file
"""
    result = instrument(source)
    assert result.count("MTgen.seed") == 1
    assert result.count("save_collision_enabled_common_state_endpoint(") == 2
    assert result.count("N_e_elastic++") == 1
    assert result.count("N_e_excitation++") == 1
    assert result.count("N_e_ionization++") == 1
    assert result.count("N_i_isotropic++") == 1
    assert result.count("N_i_backward++") == 1
    assert "if ((t + 1) == N_T" in result
    assert "cycle == cycles_done + no_of_cycles" in result
    assert "global_pre_push_step" in result
    print("collision-enabled common-state instrumenter tests passed")


if __name__ == "__main__":
    main()
