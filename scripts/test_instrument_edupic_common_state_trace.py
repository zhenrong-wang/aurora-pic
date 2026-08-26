#!/usr/bin/env python3
"""Focused tests for the collision-free common-state trace transform."""

import unittest

from instrument_edupic_common_state_trace import TRACE_STEPS, instrument


class TraceInstrumenterTests(unittest.TestCase):
    def test_transform(self) -> None:
        source = """void do_one_cycle (void){
        solve_Poisson(rho,Time);                                                 // compute potential and electric field
        for (k=0; k<N_e; k++){                              // checking for occurrence of a collision for all electrons in every time step
        }
            for (k=0; k<N_i; k++){
                vx_a = RMB(MTgen);
            }
}
"""
        transformed = instrument(source)
        self.assertIn("save_common_state_trace(t + 1, rho)", transformed)
        self.assertIn("if (false) for (k=0; k<N_e; k++)", transformed)
        self.assertIn("if (false) for (k=0; k<N_i; k++)", transformed)
        self.assertEqual(TRACE_STEPS[0], 1)
        self.assertEqual(TRACE_STEPS[-1], 4000)


if __name__ == "__main__":
    unittest.main()
