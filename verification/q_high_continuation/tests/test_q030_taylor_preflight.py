"""Executable CI gate for frozen q=.30 641/801 Taylor sources.

Unlike the previous CI, this test really calls cases() and inspects the
frozen data: a green CI means the checksum/schema checked, NOT a QNM proof.
"""
import sys
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
from q030_taylor_coefficients import cases

class FrozenTaylorPreflight(unittest.TestCase):
    def test_exact_payload_sha_and_all_three_high_resolution_sources(self):
        d=cases()
        self.assertEqual(set(d),{"n641_a0","n801_a0","n801_am0005"})
        for key,row in d.items():
            self.assertEqual(row["q"],.3)
            self.assertEqual(row["taylor_order"],12)
            self.assertEqual(row["strict_validity_radius"],.004)
            self.assertEqual(len(row["physical_spectral_sha256"]),64)
            self.assertEqual(row["nxi"],641 if key=="n641_a0" else 801)
            self.assertTrue(row["errors"])
        print("TAYLOR_Q030_PAYLOAD_SHA_AND_THREE_SOURCES_PASS")
        for key,row in d.items():
            print("Q030_TAYLOR_STRUCTURE",key,sorted(row))
            for subkey,val in row.items():
                if isinstance(val,dict):
                    print("Q030_TAYLOR_SUBKEYS",key,subkey,sorted(val))

if __name__=="__main__":
    unittest.main()
