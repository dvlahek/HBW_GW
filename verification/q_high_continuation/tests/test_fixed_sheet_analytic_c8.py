"""Physical q0.10 and synthetic knot-crossing tests for fixed monodromy C8."""
from __future__ import annotations
import sys
import unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"verification/q_high_continuation"))
sys.path.insert(0,str(ROOT/"code"))
from fixed_sheet_analytic_c8 import (
    run,FixedBranch,FixedSheetResponse,evaluate_stress,ANCHOR_STRESS
)
from retarded_analyticity_gate import (
    pack_catalog,stable_chi,PackedBranch,REFERENCE_Q010
)
import hbh_phase_xxic8_landau_continuation as c8


class FrozenSheetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with np.load(
            ROOT/"verification/blind_multiple_shooting/blind_q010_spectral_nodes.npz",
            allow_pickle=False) as f:
            d={key:np.asarray(f[key]) for key in f.files}
        cls.radial,cls.carrier=c8.build_catalog(d)
        cls.packed=pack_catalog(cls.radial,cls.carrier)

    def test_known_q010_is_locally_holomorphic_and_physically_anchor_stable(self):
        result=evaluate_stress(self.packed,REFERENCE_Q010)
        self.assertTrue(result["analyticity_on_each_frozen_sheet_pass"])
        self.assertTrue(result["cross_anchor_physical_stability_pass"])
        self.assertLess(result["max_Cauchy_Riemann_directional_mismatch"],1e-6)
        self.assertLess(result["max_relative_derivative_anchor_variation"],.02)
        self.assertEqual(len(result["anchor_stress"]),len(ANCHOR_STRESS))
        self.assertFalse(result["global_QNM_claim_permitted"])
        self.assertFalse(result["physical_sheet_uniqueness_proven"])

    def test_one_anchor_matches_original_c8_at_its_crossing_frequency(self):
        z=REFERENCE_Q010
        new=FixedSheetResponse(self.packed,z.real).evaluate(z)
        old=stable_chi(z,self.packed)
        for field in ("chiJ","chiX"):
            self.assertLess(abs(new[field]-old[field])/abs(old[field]),1e-4)
        self.assertEqual(len(new["active_branch_counts"]),3)

    def test_real_frequency_crosses_knots_but_frozen_sheet_is_analytic(self):
        x=np.array([.376,.37965,.37970,.37975,.381,.383])
        vals=np.array([1e-4,2e-4,4e-4,1e-4,-1e-4,2e-4])
        br=c8.prepare_branch(x,vals)
        p=PackedBranch([br])
        frozen=FixedBranch(p,anchor=.379701)
        z=.379701-.024j
        h=1e-4
        vrp,_=frozen.value_and_derivative(z+h)
        vrm,_=frozen.value_and_derivative(z-h)
        vip,_=frozen.value_and_derivative(z+1j*h)
        vim,_=frozen.value_and_derivative(z-1j*h)
        d0=frozen.value_and_derivative(z)[1]
        dreal=(vrp-vrm)/(2*h)
        dimag=(vip-vim)/(2j*h)
        self.assertGreater(p.knots_within(z,h),0)
        self.assertLess(abs(dreal-dimag)/abs(d0),1e-5)
        self.assertLess(abs(dreal-d0)/abs(d0),1e-5)
        self.assertLess(abs(dimag-d0)/abs(d0),1e-5)

    def test_prevent_upper_sheet_and_same_input_from_generating_new_pole_claim(self):
        p=self.packed[0]
        frozen=FixedBranch(p,REFERENCE_Q010.real)
        with self.assertRaises(ValueError):
            frozen.value_and_derivative(REFERENCE_Q010.real+1e-5j)
        report=evaluate_stress(self.packed,REFERENCE_Q010)
        self.assertIs(report["exterior_gravitational_residue_computed"],False)


if __name__=="__main__":
    unittest.main()
