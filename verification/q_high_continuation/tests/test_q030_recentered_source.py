"""Independent q=.30 recentered C8 source integrity/preflight tests.

This checks the numerical Taylor values, 641/801 provenance and hard
complex-radius guard. It does NOT certify the physical retarded sheet.
"""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from q030_taylor_q025_center import CENTER,RADIUS,SPECTRAL_SHA,cases,evaluate

class TestQ030Recentered(unittest.TestCase):
    def test_real_high_resolution_inputs_and_q025_center(self):
        d=cases()
        self.assertEqual(set(d),{"n641_a0","n801_a0","n801_am0005"})
        self.assertAlmostEqual(CENTER.real,.3810200143912758,places=13)
        self.assertAlmostEqual(CENTER.imag,-.025862906485263668,places=13)
        self.assertEqual(RADIUS,.004)
        for name,src in d.items():
            self.assertEqual(src["source_spectral_sha256"],SPECTRAL_SHA[name])
            self.assertEqual(len(SPECTRAL_SHA[name]),64)
            self.assertIs(src["physical_retarded_sheet_uniqueness_proven"],False)
            self.assertIs(src["global_qnm_confirmed"],False)
            self.assertEqual(src["q"],.3)
            self.assertEqual(src["nxi"],641 if name=="n641_a0" else 801)
            self.assertEqual(len(src["chiJ_coeff_scaled"]),13)
            self.assertEqual(len(src["chiX_coeff_scaled"]),13)
            with self.assertRaisesRegex(ValueError,"extrapolation"):
                evaluate(src,CENTER+.005)
    def test_independent_source_centers_and_analytic_locality(self):
        x=cases()
        values={name:evaluate(case,CENTER) for name,case in x.items()}
        self.assertLess(abs(values["n801_a0"]["chiJ"]-
                            complex(.06447048354304702,1.5536895727522444)),1e-12)
        self.assertLess(abs(values["n801_a0"]["chiX"]-
                            complex(1.9551187765346099,1.1404328976717615)),1e-12)
        for label in ("chiJ","chiX"):
            base=values["n801_a0"][label]
            self.assertLess(abs(values["n641_a0"][label]-base)/abs(base),.01)
            self.assertLess(abs(values["n801_am0005"][label]-base)/abs(base),.01)
            src=x["n801_a0"]
            h=1e-5
            re=(evaluate(src,CENTER+h)[label]-evaluate(src,CENTER-h)[label])/(2*h)
            im=(evaluate(src,CENTER+1j*h)[label]-evaluate(src,CENTER-1j*h)[label])/(2j*h)
            self.assertLess(abs(re-im)/abs(im),1e-6)

if __name__=="__main__":
    unittest.main()
