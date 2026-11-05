#!/usr/bin/env python3
"""Independent retarded susceptibility preflight for new physical high-q nodes.

This computes chiJ, chiX from the ORIGINAL C8 piecewise-linear Cauchy
transform on the retarded second sheet. It does NOT compute a global QNM
pole or a gravitational Green-function residue. The real purpose is to
audit the physical new-q input before the expensive global root search.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"code"))
import hbh_phase_xxic8_landau_continuation as c8

ROOT_Q020=0.3797013686231943-0.0239854837986382j
ROOT_Q010=0.3775148159474828-0.020871633919164j
UPPER_REFERENCE=0.370658573899+0.00590424942625j


def pack(z):
    return [float(z.real),float(z.imag)]


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()


def run(manifest,out):
    # Use same SHA, equilibrium and 61/81 spectral convergence checks as
    # full new-q intake; do not bypass on standalone diagnostics.
    from run_new_q import validate_manifest,load_spectral_arrays
    q,m,b,files=validate_manifest(manifest)
    nodes=load_spectral_arrays(files["spectral_nodes"],int(b["shell_count"]))
    radial,carrier=c8.build_catalog(nodes)
    if len(radial)!=1440 or len(carrier)!=80:
        raise RuntimeError("This preflight is only certified for the frozen 80-shell 18-harmonic experiment")
    probes={"published_q020_mode_seed":ROOT_Q020,
            "published_q010_mode_seed":ROOT_Q010,
            "upper_half_plane_reference":UPPER_REFERENCE}
    result={
        "q":q,
        "input_SHA256":sha256(files["spectral_nodes"]),
        "background_SHA256":sha256(files["background"]),
        "real_new_q_stationary_background_input":True,
        "published_seed_frequencies_are_not_new_q_roots":True,
        "q_high_complex_global_root_computed":False,
        "exterior_gravitational_residue_computed":False,
        "radial_branches":len(radial),"carrier_branches":len(carrier),
        "probes":{},
    }
    for label,w in probes.items():
        value=c8.evaluate_chi(w,radial,carrier,second_sheet=w.imag<0)
        if not all(np.isfinite(z) for z in (value["chiJ"],value["chiX"],value["PH"])):
            raise RuntimeError("Nonfinite retarded susceptibility")
        result["probes"][label]={
            "omega":[float(w.real),float(w.imag)],
            "chiJ":pack(value["chiJ"]),
            "chiX":pack(value["chiX"]),
            "carrier_PH":pack(value["PH"]),
            "on_retarded_second_sheet":bool(w.imag<0),
        }
    out=Path(out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n")
    print("Q_HIGH_RETARDED_SPECTRAL_SEED_PREFLIGHT_PASS",json.dumps({
        "q":q,
        "new_q_pole_not_yet_solved":True,
        "seed_chiJ":result["probes"]["published_q020_mode_seed"]["chiJ"],
        "seed_chiX":result["probes"]["published_q020_mode_seed"]["chiX"],
        "sha256":result["input_SHA256"],
    },separators=(",",":")),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",required=True,type=Path)
    p.add_argument("--output",type=Path,
                   default=ROOT/"results/q_high_continuation/retarded_seed_probe.json")
    o=p.parse_args()
    run(o.manifest,o.output)
