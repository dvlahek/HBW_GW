#!/usr/bin/env python3
"""Bounded, NONPUBLICATION q=.30 smooth global-QNM *candidate* search.

This is the first real 7x4 C9/ECS determinant evaluation of independently
regenerated 641/801-node high-q Taylor susceptibilities. It compares
three frozen discrete retarded sheets with the EXACT same bulk operator,
old q=.20 complex seed, lambda seed pair and constrained search domain.

Do not call a root "physical" unless the real kinetic retarded sheet
is independently fixed; anchor sensitivity and nxi convergence are
necessary but not sufficient. Neither waveform nor detector strain
is present. No extrapolation past the .004 local Taylor validity disk.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"code"))
from q030_taylor_regenerated import cases, evaluate
import hbh_phase_xxic9_final_global_smooth_qnm_v4 as c9

OLD_Q020_ROOT=0.3797013686231943-0.0239854837986382j
OLD_Q020_LAMBDA_PAIR=(
    -0.0926645248365645+0.0917935565479283j,
    -0.1295576973856576-0.0917935565479283j,
)
LIMIT=.0024
MAX_NFEV=25
TARGET_ABS_NORMALIZED_D=2e-6
AGREEMENT_ROOT=2e-4
SOURCES=("n641_a0","n801_a0","n801_am0005")


def score_one(case,xivc,Mfun,names):
    def local_chi(w,radial,carrier):
        return evaluate(case,w)
    c9.smooth_chi=local_chi
    c9.ACTIVE_LAMBDA_SEED_PAIR=OLD_Q020_LAMBDA_PAIR
    determinant=c9.SmoothGlobalDet(
        xivc,Mfun,names,radial=[],carrier=[],
        q=.3,tau=1.0,edge_t=1.0,Rmax=80.)
    out={
      "q":.3,
      "nxi":case["nxi"],
      "anchor_real_Momega":case["anchor_real_Momega"],
      "source_spectrum_sha256":case["source_spectral_sha256"],
      "published_q020_seed":[OLD_Q020_ROOT.real,OLD_Q020_ROOT.imag],
      "max_allowed_root_displacement":float((2**.5)*LIMIT),
      "taylor_validity_radius":case["strict_validity_radius"],
      "numerical_global_root_candidate_only":True,
      "unique_physical_retarded_sheet_not_proved":True,
    }
    start=time.time()
    try:
        d0=determinant.raw(OLD_Q020_ROOT)
        out["normalized_global_det_at_old_q020_seed"]=[float(d0.real),float(d0.imag)]
        out["abs_seed_normalized_det"]=float(abs(d0))
    except Exception as e:
        out.update(status="SEED_EVALUATION_FAIL",failure=repr(e),
                   runtime_s=time.time()-start)
        return out
    def objective(x):
        w=OLD_Q020_ROOT+complex(float(x[0]),float(x[1]))
        # rectangle |dx|,|dy|<=.0024 lies inside disk radius .004
        return determinant.real2([w.real,w.imag])
    try:
        ans=least_squares(
            objective,np.zeros(2),
            bounds=([-LIMIT,-LIMIT],[LIMIT,LIMIT]),
            method="trf",max_nfev=MAX_NFEV,
            xtol=3e-10,ftol=3e-10,gtol=3e-10,
        )
        w=OLD_Q020_ROOT+complex(*map(float,ans.x))
        dd=determinant.raw(w)
        valid=bool(ans.success and abs(dd)<TARGET_ABS_NORMALIZED_D
                   and abs(w-OLD_Q020_ROOT)<.004
                   and np.all(np.abs(ans.x)<LIMIT*.995))
        out.update(
            status=("LOCAL_NUMERICAL_CANDIDATE" if valid
                    else "NO_CERTIFIED_ROOT_WITHIN_FROZEN_LOCAL_DISK"),
            candidate_Momega=[float(w.real),float(w.imag)],
            root_shift_from_seed_abs=float(abs(w-OLD_Q020_ROOT)),
            abs_normalized_global_det_at_candidate=float(abs(dd)),
            solver_success=bool(ans.success),
            solver_nfev=int(ans.nfev),
            solver_message=str(ans.message),
            passed_local_residual_and_interior_gate=valid,
        )
    except Exception as e:
        out.update(status="GLOBAL_SOLVER_EVALUATION_FAIL",
                   failure=repr(e),
                   determinant_last_exception=determinant.last_exception)
    out["runtime_s"]=float(time.time()-start)
    return out


def run(directory):
    output=Path(directory)
    output.mkdir(parents=True,exist_ok=True)
    data=cases()
    print("Q030_GLOBAL_BEGIN",json.dumps({
        "case_ids":list(SOURCES),"local_search_limit_each_component":LIMIT,
        "target_D_abs":TARGET_ABS_NORMALIZED_D,
        "physical_sheet_selection_unresolved":True,
    }),flush=True)
    xv=c9.load_xivc()
    phase_ix=xv.load_phase_ix()
    Mfun,names=c9.build_generic_bulk_matrix(phase_ix)
    results={}
    for name in SOURCES:
        out=score_one(data[name],xv,Mfun,names)
        results[name]=out
        (output/f"{name}.json").write_text(json.dumps(out,indent=2)+"\n")
        print("Q030_GLOBAL_CASE",json.dumps({
            "case":name,"status":out["status"],
            "candidate":out.get("candidate_Momega"),
            "abs_det":out.get("abs_normalized_global_det_at_candidate"),
            "sec":out["runtime_s"],
        }),flush=True)
    good=[name for name in SOURCES if results[name].get("passed_local_residual_and_interior_gate")]
    same_root=False
    root_spread=None
    if len(good)==len(SOURCES):
        ws=[complex(*results[name]["candidate_Momega"]) for name in SOURCES]
        root_spread=float(max(abs(z-ws[1]) for z in ws))
        same_root=root_spread<=AGREEMENT_ROOT
    verdict={
        "phase":"q030_nested_7x4_smooth_global_exploratory",
        "status":("NUMERICAL_MULTI_CASE_CANDIDATE" if same_root
                  else "INCONCLUSIVE_NO_RESOLUTION_AND_ANCHOR_STABLE_ROOT"),
        "global_model":"published smooth C9 effective ocean transfer + vacuum ECS, not the exact finite-shell radial perturbation problem",
        "physical_retarded_sheet_uniquely_selected":False,
        "physical_q030_QNM_confirmed":False,
        "gravitational_residue_computed":False,
        "astrophysical_detector_signal_computed":False,
        "max_cross_source_complex_root_spread":root_spread,
        "tolerance_max_root_spread":AGREEMENT_ROOT,
        "case_results":results,
        "next_gates":[
            "independent causal retarded action-angle contour/sheet selection",
            "continuation over physical q, N, orbital nxi and regulator path",
            "full-state analyticity and exterior Green-function residue",
        ],
    }
    (output/"q030_global_comparison.json").write_text(json.dumps(verdict,indent=2)+"\n")
    print("Q030_GLOBAL_FINAL",json.dumps({
        "status":verdict["status"],"candidate_count":len(good),
        "max_spread":root_spread,
        "physical_QNM_confirmed":False,
    }),flush=True)
    return verdict


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default=str(ROOT/"results/q030_global_local"))
    args=p.parse_args()
    run(args.output)
