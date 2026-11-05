#!/usr/bin/env python3
"""Source-matched vacuum-220 normalizations AFTER exact-C1 length stress.

This probes effective C9 Green-function residues at r_obs=25 M with
identical exterior sources, normalized against the separately calculated
Schwarzschild vacuum ell=2 n=0 pole residues (at its OWN different
complex frequency). It is NOT the gravitational strain at infinity,
a merger-excitation amplitude, or a physical QNM proof.
"""
from __future__ import annotations
import argparse,csv,json,math,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"code"))
import hbh_phase_xxic9_final_global_smooth_qnm_v4 as c9
import hbh_np_source_excitation_test_v2 as source
from q025_taylor_regenerated import cases as q025_cases,evaluate as q025_chi
from q030_taylor_q025_center import cases as q030_cases,evaluate as q030_chi
from c1_metric_vs_c9_length_stress import OLD_LAMBDA_SEEDS,X_PHYS

def vacuum_source_benchmark():
    vals={}
    for kind,path,idkey in (
        ("point","results/source_data/Fig4_point_excitation_ratios.csv","R_source_over_M"),
        ("smooth","results/source_data/Fig4_smooth_excitation_ratios.csv","source_name"),
    ):
        with (ROOT/path).open(newline="") as f:
            rows=list(csv.DictReader(f))
        if not rows:raise RuntimeError("Missing published vacuum 220 source-reference table")
        for row in rows:
            key=(kind,str(float(row[idkey])) if kind=="point" else row[idkey])
            value=float(row["vacuum_residue_abs"])
            if value<=0 or not math.isfinite(value):
                raise ValueError("Invalid vacuum 220 reference")
            if key in vals and abs(vals[key]-value)>1e-10:
                raise ValueError("Vacuum same-source residue unexpectedly depends on old kinetic q case")
            vals[key]=value
    return vals

def run(length_audit,output_dir):
    source_path=Path(length_audit)
    report=json.loads(source_path.read_text())
    candidates=report["results"]
    refs=vacuum_source_benchmark()
    xivc=c9.load_xivc()
    Mfun,names=c9.build_generic_bulk_matrix(xivc.load_phase_ix())
    original_transfer=c9.ocean_transfer
    rows=[]
    summary=[]
    output_dir=Path(output_dir);output_dir.mkdir(parents=True,exist_ok=True)
    try:
        for q,load,chi in ((.25,q025_cases,q025_chi),(.30,q030_cases,q030_chi)):
            data=load()
            for candidate in [r for r in candidates if abs(r.get("q",-1)-q)<1e-12]:
                name=candidate["spectral_scheme"]
                if not candidate.get("numeric_candidate_pass",False):
                    raise RuntimeError("No accepted geometry-stressed local QNM candidate "+str((q,name)))
                qeff=math.exp(-X_PHYS[q]/9)
                def fix_length(xivc_,Mfun_,names_,w,radial,carrier,q_,tau,edge_t):
                    if abs(q_-q)>1e-10:
                        raise RuntimeError("Source physical q different from chosen background")
                    return original_transfer(xivc_,Mfun_,names_,w,radial,carrier,qeff,tau,edge_t)
                c9.ocean_transfer=fix_length
                c9.smooth_chi=lambda z,radial,carrier,case=data[name]:chi(case,z)
                c9.ACTIVE_LAMBDA_SEED_PAIR=OLD_LAMBDA_SEEDS
                w=complex(*candidate["geometry_length_corrected_root"])
                case_summary,points,smooth=source.run_case(q,[],[],w,OLD_LAMBDA_SEEDS,
                                            xivc=xivc,Mfun=Mfun,row_names=names)
                if case_summary["dW_derivative_relative_difference"]>.025:
                    raise RuntimeError("Geometry-stressed W derivative failed convergence")
                all_ratios={}
                for kind,items,field in (("point",points,"R_source"),
                                         ("smooth",smooth,"source_name")):
                    ratios=[]
                    for r in items:
                        src=(str(float(r[field])) if kind=="point" else r[field])
                        ref=refs[(kind,src)]
                        ratio=float(r["residue_abs"]/ref)
                        ratios.append(ratio)
                        rows.append({
                            "q":q,"case":name,"source_kind":kind,
                            "source_id":src,"r_obs_M":25.,
                            "matter_mode_omega_re":w.real,"matter_mode_omega_im":w.imag,
                            "vacuum_ref_mode":"Schwarzschild ell2 n0 at its own complex frequency",
                            "effective_C9_source_residue_abs":float(r["residue_abs"]),
                            "vacuum_220_source_residue_abs":ref,
                            "relative_pole_residue_abs_ratio":ratio,
                        })
                    all_ratios[kind]={"min":min(ratios),"median":float(np.median(ratios)),
                                      "max":max(ratios),"count":len(ratios)}
                summary.append({
                    "q":q,"case":name,
                    "geometry_stressed_omega":[w.real,w.imag],
                    "source_dW_step_relative_error":case_summary["dW_derivative_relative_difference"],
                    "ratio_statistics":all_ratios,
                    "metric_strain_excitation_confirmed":False,
                })
                print("C1_C9_LENGTH_STRESSED_SOURCE_BENCHMARK",json.dumps({
                    "q":q,"case":name,
                    "point_median":all_ratios["point"]["median"],
                    "smooth_median":all_ratios["smooth"]["median"],
                    "step_dW_rel":case_summary["dW_derivative_relative_difference"]
                },separators=(",",":")),flush=True)
    finally:
        c9.ocean_transfer=original_transfer

    with (output_dir/"effective_C9_length_stressed_same_source_ratios.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    report={
        "status":"EFFECTIVE_C9_LENGTH_STRESSED_VACUUM_220_SOURCE_BENCHMARK",
        "physical_metric_green_residue_confirmed":False,
        "physical_kinetic_sheet_uniquely_selected":False,
        "radiative_strain_at_infinity_computed":False,
        "merger_excitation_or_SNR_calculated":False,
        "vacuum_comparison":"pole-residue moduli of same sources at two different complex mode frequencies",
        "observer_R_over_M":25.,
        "source_reference_paths":[
            "results/source_data/Fig4_point_excitation_ratios.csv",
            "results/source_data/Fig4_smooth_excitation_ratios.csv"
        ],
        "results":summary,
    }
    (output_dir/"effective_C9_length_stressed_benchmark.json").write_text(json.dumps(report,indent=2)+"\n")
    print("C1_C9_LENGTH_STRESSED_BENCHMARK_FINAL",json.dumps({
        "source_case_count":len(summary),
        "physical_strain_claimed":False,
    },separators=(",",":")),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geometry-json",default=str(ROOT/"results/q_high_geometry_length_stress.json"))
    p.add_argument("--out",default=str(ROOT/"results/q_high_geometry_source"))
    args=p.parse_args()
    run(args.geometry_json,args.out)
