"""Verify exact endpoint claims and archived q-family provenance."""
from __future__ import annotations
import sys
import unittest
from pathlib import Path
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
from q_to_one_endpoint_audit import (
    thickness,transfer_identity_q1,archived_observations,
    source_audit,Q_CONFIRMED,SOURCE_Q,UNRESOLVED_Q,M_POLE,VAC_POLE
)

class QToOneTests(unittest.TestCase):
    def test_exact_q1_transfer_is_vacuum_identity(self):
        x=transfer_identity_q1()
        self.assertTrue(x["T_equals_identity_exactly"])
        self.assertEqual(x["thickness_M"],0)
        self.assertEqual(x["det_T"],1.)

    def test_shell_length_decreases_to_exact_vacuum(self):
        qs=[.05,.1,.2,.5,.9,.99,.999,1.]
        ls=[thickness(q) for q in qs]
        self.assertTrue(all(a>b for a,b in zip(ls,ls[1:])))
        self.assertLess(thickness(.99),.1)
        with self.assertRaises(ValueError):
            thickness(1.01)

    def test_reproduces_archived_poles_without_approximating_new_q(self):
        x=archived_observations()
        self.assertEqual(len(x),len(Q_CONFIRMED))
        self.assertEqual([r["q"] for r in x],list(Q_CONFIRMED))
        self.assertTrue(all(r["root_minus_im_Momega"]>0 for r in x))
        self.assertTrue(all(r["published_root_status"].startswith("archived") for r in x))
        self.assertTrue(all(q>.2 for q in UNRESOLVED_Q))

    def test_exterior_source_residue_ratios_are_not_injection_amplitudes(self):
        for file,n in (("results/source_data/Fig4_point_excitation_ratios.csv",6),
                       ("results/source_data/Fig4_smooth_excitation_ratios.csv",3)):
            audit=source_audit(file)
            self.assertEqual(len(audit),2)
            self.assertEqual([a["source_count"] for a in audit],[n,n])
            self.assertTrue(all(a["vacuum_normalized_residue_ratio_max"]<.05 for a in audit))
            self.assertGreater(audit[1]["vacuum_normalized_residue_ratio_median"],
                               audit[0]["vacuum_normalized_residue_ratio_median"])

    def test_q01_matter_pole_distinct_from_schwarzschild_gravitational(self):
        self.assertGreater(abs(M_POLE-VAC_POLE),.06)

if __name__=="__main__":
    unittest.main()
