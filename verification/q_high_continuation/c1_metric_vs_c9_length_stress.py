#!/usr/bin/env python3
"""Physical-background tortoise-length stress test of EFFECTIVE C9 poles.

An independent exact-C1 metric audit integrated sqrt(g_rr)/alpha over
all 80 material shells AND intervening Schwarzschild vacuum gaps:
q=.25: X_phys=12.399467716175225; C9 length 12.476649250079015
q=.30: X_phys=10.768902085319294; C9 length 10.835755238933425
(0.62% length difference; independent grid change <2e-10).
Instead of silently claiming identical geometries, stress C9 by replacing
ONLY its ideal ocean length 9log(1/q) with measured X_phys. The C8
kinetic response, C9 7x4 operator and source normalization are UNCHANGED.

This does NOT derive the true finite-shell perturbation operator. Any
new root/residue is a controlled effective-model sensitivity result ONLY.
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"code"))
import hbh_phase_xxic9_final_global_smooth_qnm_v4 as c9
from q025_taylor_regenerated import cases as q025_cases, evaluate as q025_chi
from q030_taylor_q025_center import cases as q030_cases, evaluate as q030_chi

X_PHYS={.25:12.399467716175225,.30:10.768902085319294}
X_GEOM_COARSE={.25:12.399467733038104,.30:10.768902096500526}
ROOT_N801={
    .25:0.3810200143912758-0.025862906485263668j,
    .30:0.3832546292851808-0.027778360028813374j,
}
ROOT_OTHER={
    .25:{"n641_a0":0.38099074548973594-0.025873483312466364j,
          "n801_am0005":0.3809538387456802-0.025896337962632018j},
    .30:{"n641_a0":0.38326291395214707-0.02776839393639531j,
          "n801_am0005":0.3831737223916984-0.027820785440657844j},
}
OLD_LAMBDA_SEEDS=(
 -0.0926645248365645+0.0917935565479283j,
 -0.1295576973856576-0.0917935565479283j
)
MAX_LOCAL_COMPONENT=.00279
ABS_D_CUTOFF=2e-6

def good_geometry():
    for q,L in X_PHYS.items():
        L0=9*math.log(1/q)
        if not (0<L<L0 and abs(L-X_GEOM_COARSE[q])<2e-7):
            raise RuntimeError("Unexpected reconstructed C1 N80 geometry length")
    return True

def one(q,name,case,evaluate,Mfun,names,xivc):
    q_eff=math.exp(-X_PHYS[q]/9)
    w0=ROOT_N801[q] if name=="n801_a0" else ROOT_OTHER[q][name]
    center=complex(*case["center_Momega"])
    radius=case["strict_validity_radius"]
    if abs(w0-center)>=radius:
        raise ValueError("Existing candidate outside strict source analytic disk")

    original=c9.ocean_transfer
    c9.smooth_chi=lambda z,radial,carrier:evaluate(case,z)
    c9.ACTIVE_LAMBDA_SEED_PAIR=OLD_LAMBDA_SEEDS

    def run_det(length_adjusted):
        if length_adjusted:
            def changed(xivc_,Mfun_,names_,w,radial,carrier,q_,tau,edge_t):
                if abs(q_-q)>1e-12:
                    raise RuntimeError("Unexpected physical mass ratio in operator")
                return original(xivc_,Mfun_,names_,w,radial,carrier,
                                q_eff,tau,edge_t)
            c9.ocean_transfer=changed
        else:
            c9.ocean_transfer=original
        return c9.SmoothGlobalDet(xivc,Mfun,names,[],[],q=q,tau=1.,edge_t=1.,Rmax=80.)
    try:
        Dold=run_det(False)
        previous=Dold.raw(w0)
        Dnew=run_det(True)
        initial=Dnew.raw(w0)
        def F(x):
            w=center+complex(float(x[0]),float(x[1]))
            return Dnew.real2([w.real,w.imag])
        xinit=np.array([(w0-center).real,(w0-center).imag])
        xinit=np.clip(xinit,-MAX_LOCAL_COMPONENT+.000001,MAX_LOCAL_COMPONENT-.000001)
        result=least_squares(
            F,xinit,
            bounds=([-MAX_LOCAL_COMPONENT]*2,[MAX_LOCAL_COMPONENT]*2),
            method="trf",max_nfev=65,
            xtol=1e-11,ftol=1e-11,gtol=1e-11)
        w=center+complex(float(result.x[0]),float(result.x[1]))
        if abs(w-center)>=radius:
            raise RuntimeError("Attempted geometry-stressed root outside frozen analytic disk")
        D=Dnew.raw(w)
        accepted=(result.success and abs(D)<ABS_D_CUTOFF
                  and all(abs(result.x)<.995*MAX_LOCAL_COMPONENT))
        return {
            "q":q,"spectral_scheme":name,
            "physical_C1_tortoise_length":X_PHYS[q],
            "ideal_C9_tortoise_length":9*math.log(1/q),
            "model_q_parameter_only_for_length_override":q_eff,
            "physical_mass_ratio_remains_q":q,
            "same_C8_response_and_C9_bulk_operator":True,
            "is_exact_finite_shell_perturbation_validation":False,
            "uncorrected_complex_root":[w0.real,w0.imag],
            "uncorrected_D_at_own_root_abs":float(abs(previous)),
            "physical_length_D_at_uncorrected_root_abs":float(abs(initial)),
            "geometry_length_corrected_root":[w.real,w.imag],
            "absolute_complex_root_shift":float(abs(w-w0)),
            "geometry_length_corrected_D_abs":float(abs(D)),
            "solver_nfev":int(result.nfev),
            "numeric_candidate_pass":bool(accepted),
            "sheet_uniqueness_not_proven":True,
            "observational_detection_not_claimed":True,
        }
    finally:
        c9.ocean_transfer=original

def run(output):
    good_geometry()
    xivc=c9.load_xivc()
    phase=xivc.load_phase_ix()
    Mfun,names=c9.build_generic_bulk_matrix(phase)
    results=[]
    for q,get_data,eval_chi in ((.25,q025_cases,q025_chi),(.30,q030_cases,q030_chi)):
        src=get_data()
        for name in ("n641_a0","n801_a0","n801_am0005"):
            try:
                row=one(q,name,src[name],eval_chi,Mfun,names,xivc)
            except Exception as exc:
                row={"q":q,"spectral_scheme":name,
                     "status":"INCONCLUSIVE_NUMERICAL_GEOMETRY_STRESS",
                     "exception":repr(exc),
                     "physical_QNM_confirmed":False}
            print("C1_GEOMETRY_C9_LENGTH_STRESS",json.dumps(row,separators=(",",":")),flush=True)
            results.append(row)
    summary={
        "scope":"effective-C9 one-parameter ocean-length sensitivity, not exact 80-shell linearized Einstein-Vlasov dynamics",
        "physical_metric_tortoise_length_integrated_over_exact_N80_material_and_vacuum_gaps":X_PHYS,
        "source":"independently re-integrated local user C1-C7 80-shell physical backgrounds",
        "status":"ALL_NUMERICAL_STRESS_ROOTS_PASS" if all(r.get("numeric_candidate_pass") for r in results) else "INCONCLUSIVE_NUMERICAL_STRESS",
        "model_operator_physically_certified":False,
        "causal_kinetic_sheet_uniquely_selected":False,
        "radiative_strain_or_merger_excitation_computed":False,
        "results":results,
    }
    path=Path(output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(summary,indent=2)+"\n")
    print("C1_VS_C9_GEOMETRY_FINAL",json.dumps({
        "status":summary["status"],
        "max_root_shift":max([r.get("absolute_complex_root_shift",0) for r in results]),
        "physical_operator_certified":False
    },separators=(",",":")),flush=True)
    return summary

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default=str(ROOT/"results/q_high_geometry_length_stress.json"))
    args=p.parse_args()
    run(args.output)
