"""No invented q>.20 spectral roots: physical input/provenance contract tests."""
from __future__ import annotations
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"verification/q_high_continuation"))
from run_new_q import (
    inspect_reference_q010,validate_manifest,load_spectral_arrays,
    inspect_new_q,seed_from_predecessor,sha256,
    EXPECTED_MODEL,ALLOWED_SOURCE_Q
)

class HighQIntakeTests(unittest.TestCase):
    def test_original_published_q010_spectral_input_is_verified_without_new_q(self):
        obj=inspect_reference_q010()
        self.assertEqual(obj["N"],80)
        self.assertEqual(obj["radial_branches"],1440)
        self.assertEqual(obj["carrier_branches"],80)
        self.assertEqual(obj["status"],"PUBLIC Q=.10 REFERENCE ONLY, NOT q>.20 data")
        self.assertEqual(obj["sha256"],
                         "04016df710591f61023338b20117b9b3fd046f4fc3ca886060f2b4a951d9114f")

    def test_new_q_requires_manifest_and_never_accepts_q1_as_new_matter_pole(self):
        self.assertIn(.3,ALLOWED_SOURCE_Q)
        self.assertNotIn(.2,ALLOWED_SOURCE_Q)
        self.assertNotIn(1.,ALLOWED_SOURCE_Q)
        with self.assertRaises((FileNotFoundError,ValueError)):
            inspect_new_q("/totally_missing_physical_q0p300_manifest.json")

    def test_original_published_seed_is_not_mistaken_for_a_new_q_result(self):
        prev,w,lams=seed_from_predecessor(.3,None)
        self.assertEqual(prev,.2)
        self.assertAlmostEqual(w.real,.3797013686231943,places=9)
        self.assertAlmostEqual(w.imag,-.0239854837986382,places=9)
        self.assertEqual(len(lams),2)
        with self.assertRaisesRegex(ValueError,"First new q"):
            seed_from_predecessor(.9,None)

    def test_rejects_fabricated_missing_upstream_generator_and_reused_q010(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp)
            spectral=tmp/"nodes.npz"
            np.savez(spectral, dummy=np.arange(4096))
            bg=tmp/"physical_bg.json"
            background={
                "q_integrated_from_background":.3,
                "relative_mass_balance_error":1e-8,
                "einstein_constraint_max":1e-8,
                "vlasov_stationarity_max":1e-8,
                "stress_energy_conservation_max":1e-8,
                "shell_count":20,
                "ocean_length_M":9*np.log(1/.3),
                "positive_phase_space_density":True,
                "physical_orbit_support":True,
                "stress_energy_from_same_distribution":True,
            }
            bg.write_text(json.dumps(background))
            manifest={
                "schema_version":1,
                "q":.3,
                "model_family":EXPECTED_MODEL,
                "new_self_consistent_equilibrium_constructed":False,
                "q010_spectral_reused":False,
                "spectral_convergence_verified":True,
                "spectral_convergence_max_relative":.01,
                "generator":{
                    "script_sha256":"a"*64,
                    "source_description":"actual upstream generator with documented equations",
                    "command":"python external_source.py --q .3"
                },
                "background":{"file":"physical_bg.json","sha256":sha256(bg)},
                "spectral_nodes":{"file":"nodes.npz","sha256":sha256(spectral)}
            }
            man=tmp/"manifest.json";man.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,"interpolation"):
                validate_manifest(man)
            manifest["new_self_consistent_equilibrium_constructed"]=True
            manifest["spectral_nodes"]["sha256"]="b"*64
            man.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,"checksum"):
                validate_manifest(man)
            manifest["spectral_nodes"]["sha256"]=sha256(spectral)
            manifest["spectral_convergence_max_relative"]=.25
            man.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,"convergence"):
                validate_manifest(man)
            manifest["spectral_convergence_max_relative"]=.01
            manifest["spectral_convergence_verified"]=False
            man.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,"convergence"):
                validate_manifest(man)
            manifest["spectral_convergence_verified"]=True
            manifest["q010_spectral_reused"]=True
            man.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,"Reused"):
                validate_manifest(man)
            manifest["q010_spectral_reused"]=False
            man.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,"schema|fields|node"):
                inspect_new_q(man)
            background["q_integrated_from_background"]=.1
            bg.write_text(json.dumps(background))
            manifest["background"]["sha256"]=sha256(bg)
            man.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,"mass ratio"):
                validate_manifest(man)

    def test_rejects_traversal_and_absent_source_physical_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp)
            from run_new_q import within_directory
            with self.assertRaises(ValueError):
                within_directory(tmp,"../fake_q0p3_data.npz")
            with self.assertRaises(ValueError):
                within_directory(tmp,"/tmp/fake_q0p3_data.npz")
            with self.assertRaises(FileNotFoundError):
                within_directory(tmp,"not_generated_q0p3_data.npz")

if __name__=="__main__":
    unittest.main()
