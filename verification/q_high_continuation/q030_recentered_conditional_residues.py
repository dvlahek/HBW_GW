#!/usr/bin/env python3
"""Conditional q=.30 pole-excitation diagnostic on original effective C9 operator.

Inputs: independently calculated q=.30 80-shell C7 spectra, frozen
nxi641/801 C8 fixed-sheet Taylor continuations, C9 global root candidates.
Outputs: Wronksian derivative convergence and point/smooth exterior source
residues, with cross-resolution/anchor comparison. These are mathematical
residues of the EFFECTIVE-OCEAN model, NOT proven linearized metric
Green-function residues of the exact 80-shell background or merger strains.
"""
from __future__ import annotations
import argparse,json,math,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"code"))
import hbh_phase_xxic9_final_global_smooth_qnm_v4 as c9
import hbh_np_source_excitation_test_v2 as source
from q030_taylor_q025_center import cases,evaluate
from q030_recentered_global_bounded import OLD_Q020_LAMBDA_PAIR,SOURCES

TOL_DERIV=.025
TOL_SOURCE_SPREAD=.15
TOL_CROSS_SOURCE_ROOT=2e-4

def run(folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    source_report=folder/"q030_recentered_global_comparison.json"
    if not source_report.exists():
        raise FileNotFoundError("q030_recentered global candidates required before computing conditional residues")
    base=json.loads(source_report.read_text())
    if base["status"]!="NUMERICAL_MULTI_CASE_CANDIDATE":
        raise RuntimeError("No three-sheet stable global root: exterior residue calculation refused")
    d=cases()
    xivc=c9.load_xivc()
    phase_ix=xivc.load_phase_ix()
    Mfun,names=c9.build_generic_bulk_matrix(phase_ix)
    rows={}
    for label in SOURCES:
        case=d[label];candidate=base["case_results"][label]
        if not candidate.get("passed_local_residual_and_interior_gate"):
            raise RuntimeError("Source case has no interior global root")
        w=complex(*candidate["candidate_Momega"])
        def fixed(z,radial,carrier):
            return evaluate(case,z)
        c9.smooth_chi=fixed
        c9.ACTIVE_LAMBDA_SEED_PAIR=OLD_Q020_LAMBDA_PAIR
        try:
            summ,point,smooth=source.run_case(
                .30,[],[],w,OLD_Q020_LAMBDA_PAIR,
                xivc=xivc,Mfun=Mfun,row_names=names
            )
            drel=float(summ["dW_derivative_relative_difference"])
            if not np.isfinite(drel) or drel>TOL_DERIV:
                raise RuntimeError("Unnormalized W derivative unstable")
            def packed(x):
                return [float(x.real),float(x.imag)]
            out={
                "status":"EFFECTIVE_C9_SOURCE_RESIDUE_NUMERICAL_ONLY",
                "q":.30,"spectral_source":label,
                "physical_retarded_sheet_uniqueness_proven":False,
                "metric_perturbation_operator_on_exact_shells_verified":False,
                "merger_waveform_excitation_predicted":False,
                "detector_signal_analyzed":False,
                "mode_Momega":packed(w),
                "source_spectrum_SHA256":case["source_spectral_sha256"],
                "summary":summ,"point_source_residues":point,
                "smooth_source_residues":smooth,
            }
            rows[label]=out
            (folder/f"{label}_conditional_source_residues.json").write_text(json.dumps(out,indent=2)+"\n")
            print("Q030_RECENTERED_CONDITIONAL_SOURCE",json.dumps({
                "case":label,
                "dW_h2e5_h1e5_relative":drel,
                "smooth_residue_min_abs":summ["smooth_residue_min_abs"],
                "smooth_residue_max_abs":summ["smooth_residue_max_abs"],
                "point_residue_min_abs":summ["point_residue_min_abs"],
                "point_residue_max_abs":summ["point_residue_max_abs"]
            }),flush=True)
        except Exception as e:
            rows[label]={"status":"RESIDUE_NOT_NUMERICALLY_CERTIFIED",
                         "q":.30,"case":label,"error":repr(e)}
            print("Q030_RECENTERED_SOURCE_FAIL",json.dumps(rows[label]),flush=True)
    complete=all(rows[k]["status"]=="EFFECTIVE_C9_SOURCE_RESIDUE_NUMERICAL_ONLY" for k in SOURCES)
    cross={}
    if complete:
        for source_type in ("point_source_residues","smooth_source_residues"):
            keys=([("R_source",float)] if source_type=="point_source_residues"
                  else [("source_name",str)])
            ref=rows["n801_a0"][source_type]
            relative=[]
            for key in ("n641_a0","n801_am0005"):
                other=rows[key][source_type]
                if len(other)!=len(ref):
                    raise RuntimeError("Changed number of source probes")
                for a,b in zip(other,ref):
                    if source_type=="point_source_residues" and abs(a["R_source"]-b["R_source"])>1e-10:
                        raise RuntimeError("Changed point source locations")
                    if source_type=="smooth_source_residues" and a["source_name"]!=b["source_name"]:
                        raise RuntimeError("Changed smooth source family")
                    relative.append(float(abs(a["residue_abs"]-b["residue_abs"])/
                                          max(abs(b["residue_abs"]),1e-300)))
            cross[source_type+"_worst_relative_abs_residue_change"]=max(relative)
    stable=complete and bool(cross) and max(cross.values())<=TOL_SOURCE_SPREAD
    verdict={
        "phase":"q030_recentered_conditional_C9_exterior_residue",
        "status":("EXPLORATORY_RESIDUES_NUMERICALLY_STABLE" if stable
                  else "INCONCLUSIVE_EXTERIOR_RESIDUE"),
        "physical_metric_residue_confirmed":False,
        "astrophysical_merger_excitation_predicted":False,
        "detector_detectability_claimed":False,
        "physical_causal_sheet_uniqueness_proven":False,
        "source_comparison":cross,
        "max_allowed_relative_change":TOL_SOURCE_SPREAD,
        "cases":rows,
    }
    (folder/"q030_recentered_exterior_residue_comparison.json").write_text(json.dumps(verdict,indent=2)+"\n")
    print("Q030_RECENTERED_SOURCE_FINAL",json.dumps({
        "status":verdict["status"],"complete":complete,"cross":cross,
        "physical_residue_confirmed":False
    }),flush=True)
    return verdict

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default=str(ROOT/"results/q030_recentered_global_local"))
    a=p.parse_args()
    run(a.output)
