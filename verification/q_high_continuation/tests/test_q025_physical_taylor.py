import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from q025_taylor_regenerated import cases,evaluate

class TestQ025(unittest.TestCase):
    def test_payload_and_sources(self):
        d=cases()
        self.assertEqual(set(d),{'n641_a0','n801_a0','n801_am0005'})
        for x in d.values():
            self.assertEqual(x['q'],.25)
            self.assertTrue(x['physical_sheet_not_certified'])
            self.assertLess(x['worst_heldout_exact_cauchy_relative'],1e-9)
            z=complex(*x['center_Momega'])
            self.assertGreater(abs(evaluate(x,z)['chiJ']),0.)
            with self.assertRaises(ValueError):evaluate(x,z+.005)
    def test_stability(self):
        d=cases()
        for k in ('chiJ','chiX'):
            v=lambda name,i:complex(*d[name][k+'_coeff_scaled'][i])
            ref=v('n801_a0',0)
            self.assertLess(abs(v('n641_a0',0)-ref)/abs(ref),.01)
            self.assertLess(abs(v('n801_am0005',0)-ref)/abs(ref),.01)
            der=v('n801_a0',1)
            self.assertLess(abs(v('n641_a0',1)-der)/abs(der),.04)
            self.assertLess(abs(v('n801_am0005',1)-der)/abs(der),.05)

if __name__=='__main__':unittest.main()
