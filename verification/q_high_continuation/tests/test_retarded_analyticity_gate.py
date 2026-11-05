"""Check q=.10 analytic Landau baseline and freeze high-q failure behavior."""
from __future__ import annotations
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"verification/q_high_continuation"))
from retarded_analyticity_gate import (
    REFERENCE_Q010,run,pack_catalog,derivative_analyticity_gate
)

class RetardedAnalyticityTests(unittest.TestCase):
    def test_published_q010_root_has_consistent_real_imag_frequency_derivatives(self):
        # Real archived 80-shell q0.10 source, not a mock spectrum.
        result=run(reference_q010=True,strict=True)
        self.assertEqual(result["q"],.1)
        self.assertEqual(result["spectral_radial_branches"],1440)
        self.assertEqual(result["spectral_carrier_branches"],80)
        self.assertTrue(result["global_QNM_residue_gate_pass"])
        self.assertLess(result["max_relative_CR_disagreement_at_1e5"],.03)
        self.assertFalse(result["qnm_pole_computed"])
        self.assertFalse(result["gravitational_residue_computed"])

    def test_scientific_gate_cannot_be_satisfied_by_numerical_flag(self):
        from retarded_analyticity_gate import PackedBranch
        import numpy as np
        x=np.array([0.37,0.375,0.37970,0.37971,0.382,0.385],float)
        weights=np.array([1.,1.2,.95,1.05,.8,1.])
        # This deliberately nonphysical test fixture creates a spectral
        # knot at the frequency of interest, not a fake Einstein-Vlasov mode.
        from hbh_phase_xxic8_landau_continuation import prepare_branch
        branch=prepare_branch(x,weights)
        p=PackedBranch([branch])
        self.assertGreater(p.knots_within(0.379701-0.02j,1e-5),0)
        self.assertEqual(len(p.x),1)

if __name__=="__main__":
    unittest.main()
