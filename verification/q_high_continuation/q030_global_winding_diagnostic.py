#!/usr/bin/env python3
"""Frozen local argument-winding cross-check for q=.30 global QNM search.

Numerical contour phase counts are NET zeros MINUS possible poles of
the Riccati/ECS representation, not an interval-analytic certificate.
Uses the independently reconstructed nxi=641/801 C8 Taylor sources
on EXACTLY their validated fixed local disk, not extrapolated data.
"""
import argparse
import cmath
import json
import math
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"code"))
import hbh_phase_xxic9_final_global_smooth_qnm_v4 as c9
from q030_taylor_regenerated import cases,evaluate
from q030_global_bounded_candidate import (
    OLD_Q020_ROOT,OLD_Q020_LAMBDA_PAIR,SOURCES
)

CHECKS=((.003,32),(.0035,48))


def run(path):
    data=cases()
    xivc=c9.load_xivc()
    phase_ix=xivc.load_phase_ix()
    Mfun,names=c9.build_generic_bulk_matrix(phase_ix)
    out={"q":.3,
         "numerical_winding_is_not_pole_free_analytic_proof":True,
         "root_search_not_replaced_by_contour":True,
         "physical_retarded_sheet_uniqueness_proven":False,
         "original_q020_seed":[OLD_Q020_ROOT.real,OLD_Q020_ROOT.imag],
         "results":{}}
    for name in SOURCES:
        case=data[name]
        def fixed(w,radial,carrier):
            return evaluate(case,w)
        c9.smooth_chi=fixed
        c9.ACTIVE_LAMBDA_SEED_PAIR=OLD_Q020_LAMBDA_PAIR
        det=c9.SmoothGlobalDet(xivc,Mfun,names,[],[],q=.3,tau=1.0,edge_t=1.0)
        entries=[]
        for radius,n in CHECKS:
            ws=[OLD_Q020_ROOT+radius*cmath.exp(2j*math.pi*i/n)
                for i in range(n+1)]
            try:
                vals=np.array([det.raw(w) for w in ws],dtype=complex)
                if not np.all(np.isfinite(vals)) or min(abs(vals))<1e-9:
                    raise RuntimeError("Nonfinite or near-zero determinant ON contour")
                phase_inc=np.angle(vals[1:]/vals[:-1])
                winding=float(np.sum(phase_inc)/(2*math.pi))
                nearest=int(round(winding))
                well_sampled=bool(
                    abs(winding-nearest)<.05 and
                    max(abs(phase_inc))<math.pi*.75)
                entries.append({
                    "radius":float(radius),"samples":n,
                    "net_winding":winding,
                    "nearest_integer":nearest,
                    "phase_sampling_gate_pass":well_sampled,
                    "min_abs_normalized_D_on_contour":float(min(abs(vals))),
                    "max_phase_increment":float(max(abs(phase_inc))),
                    "interpretation":"argument winding of numerically sampled Riccati determinant, zeros minus possible poles",
                })
            except Exception as e:
                entries.append({"radius":float(radius),"samples":n,
                                "status":"INCONCLUSIVE_CONTOUR_EVALUATION",
                                "exception":repr(e)})
        out["results"][name]=entries
        print("Q030_CONTOUR_CASE",json.dumps({
            "name":name,"contours":entries
        },separators=(",",":")),flush=True)
    rows=[r for es in out["results"].values() for r in es]
    valid=all(x.get("phase_sampling_gate_pass",False) for x in rows)
    stable=valid and len({x["nearest_integer"] for x in rows})==1
    out["sampling_and_anchor_winding_stable"]=stable
    out["status"]=("STABLE_DISCRETE_NET_WINDING" if stable
                  else "INCONCLUSIVE_NET_WINDING")
    out["physical_q030_QNM_confirmed"]=False
    output=Path(path);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2)+"\n")
    print("Q030_CONTOUR_FINAL",json.dumps({
        "status":out["status"],
        "integer":rows[0].get("nearest_integer") if stable else None,
        "physical_pole_confirmed":False
    }),flush=True)
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default=str(ROOT/"results/q030_global_local/q030_winding.json"))
    a=p.parse_args()
    run(a.output)
