#!/usr/bin/env python3
"""A **genuine q=1 vacuum endpoint**, not fabricated q->1 continuation.

Valid results:
  1) q=1 has zero "ocean" transfer length and EXACT identity transfer.
  2) compute vacuum Zerilli/ECS determinant near the published q=.10
     matter-led pole; independently inspect its zero count on a fixed disk.
  3) audit genuinely computed q=.05..20 poles and q=.10/.20
     source-normalized vacuum/matter residues.
  4) mark ALL intermediate q>.2 as UNRESOLVED: missing newly built
     self-consistent spectral arrays. No interpolation/extrapolation.

Source data: this exact HBW_GW repository's publication-facing
q_family_results and Fig4 tables (those are archived numerical evidence,
not new QNM results). Input q=m_inner/M_outer; 1-q is the gravitational
mass difference fraction, NOT an observed material rest-mass fraction.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import root

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
import hbh_phase_xivc_global_polar_EV_ecs as vacuum

M_POLE=0.377514814460-0.020871654486j
VAC_POLE=0.373671684352-0.088962315642j
Q_CONFIRMED=(.05,.075,.1,.15,.2)
SOURCE_Q=(.1,.2)
UNRESOLVED_Q=(.3,.5,.7,.9,.95,.99,.999)
CONTOUR_R=.01
CONTOUR_N=32
RMAX=80.


def load_rows(name):
    with (ROOT/name).open(newline="") as file:
        return list(csv.DictReader(file))


def thickness(q):
    if not (0<q<=1):
        raise ValueError("Only physical mass ratio in (0,1] is allowed")
    return float(9.*math.log(1./q))


def transfer_identity_q1():
    """Check q=1 short-circuit of PUBLIC sharp-geometry solver.

    This branch deliberately requires no constitutive/spectral kernel.
    The smooth ocean solver has the same explicit q=1 identity
    condition; finite-q response CANNOT be obtained from this endpoint.
    """
    T,info=vacuum.ocean_transfer(None,M_POLE,1.)
    if np.shape(T)!=(2,2) or not np.allclose(T,np.eye(2),rtol=0,atol=0):
        raise RuntimeError("q=1 does not exactly collapse to the vacuum")
    return {"T_equals_identity_exactly":True,
            "thickness_M":thickness(1.),
            "det_T":float(np.linalg.det(T).real)}


def vacuum_D(w,Rmax=RMAX):
    """Normalized vacuum Zerilli determinant at the single q=1 endpoint."""
    zi=vacuum.vacuum_inner_state(w)
    zo=vacuum.vacuum_outer_state(w,Rmax=Rmax)
    det=zi[0]*zo[1]-zi[1]*zo[0]
    return complex(det/max(np.linalg.norm(zi)*np.linalg.norm(zo),1e-300))


def count_contour(center,radius=CONTOUR_R,n=CONTOUR_N,Rmax=RMAX):
    """Discrete winding of the vacuum determinant on a prescribed circle.

    The value is an empirical contour diagnostic, not an interval-
    arithmetic certificate. Vary contour density and Rmax to audit.
    """
    if n<16 or radius<=0:
        raise ValueError("Contour must be fixed sufficiently resolved circle")
    angle=np.linspace(0,2*math.pi,n+1)
    ww=center+radius*np.exp(1j*angle)
    dd=np.array([vacuum_D(z,Rmax) for z in ww])
    if not np.all(np.isfinite(dd)) or np.min(np.abs(dd))<1e-8:
        raise RuntimeError("q=1 vacuum determinant nearly vanishes on chosen contour")
    diff=np.angle(dd[1:]/dd[:-1])
    winding=float(np.sum(diff)/(2*math.pi))
    rounded=int(round(winding))
    if abs(winding-rounded)>0.05:
        raise RuntimeError("Poorly resolved determinant argument count")
    return {"center_re":float(center.real),"center_im":float(center.imag),
            "radius":float(radius),"N":n,
            "winding_unrounded":winding,"winding_number":rounded,
            "min_contour_abs_D":float(np.min(np.abs(dd))),
            "max_argument_step_rad":float(np.max(np.abs(diff)))}


def source_audit(name):
    rows=load_rows(name)
    result=[]
    for q in SOURCE_Q:
        cohort=[r for r in rows if abs(float(r["q"])-q)<1e-10]
        if not cohort:
            raise RuntimeError(f"No archived source residues for q={q}: {name}")
        errors=[abs(float(r["kinetic_to_vacuum_residue_ratio"])-
                    float(r["kinetic_residue_abs"])/float(r["vacuum_residue_abs"]))
                for r in cohort]
        if max(errors)>1e-11:
            raise RuntimeError("Residue ratios disagree with raw amplitudes")
        vals=[float(r["kinetic_to_vacuum_residue_ratio"]) for r in cohort]
        result.append({"q":q,"source_count":len(cohort),
                       "vacuum_normalized_residue_ratio_min":min(vals),
                       "vacuum_normalized_residue_ratio_max":max(vals),
                       "vacuum_normalized_residue_ratio_median":float(np.median(vals))})
    return result


def archived_observations():
    qrows=load_rows("results/q_family_results.csv")
    measured=sorted(float(x["q"]) for x in qrows)
    if len(qrows)!=5 or not np.allclose(measured,Q_CONFIRMED,atol=1e-12,rtol=0):
        raise RuntimeError("Published five-point q family changed unexpectedly")
    table=[]
    for row in qrows:
        q=float(row["q"])
        gamma=-float(row["root_im"])
        if gamma<=0 or not int(round(float(row["winding_number"])))==1:
            raise RuntimeError(f"Bad published q={q} root or argument count")
        table.append({
            "q":q,"medium_gravitational_mass_difference_fraction":1.-q,
            "ocean_transfer_length_M":thickness(q),
            "root_re_Momega":float(row["root_re"]),
            "root_minus_im_Momega":gamma,
            "root_quality_factor":float(row["quality_factor"]),
            "boundary_normalized_residue_abs":float(row["boundary_normalized_residue_abs"]),
            "published_root_status":"archived computed finite-q pole, not new continuation",
        })
    return table


def run(output,n=CONTOUR_N,Rmax=RMAX):
    transfer=transfer_identity_q1()
    table=archived_observations()
    p=source_audit("results/source_data/Fig4_point_excitation_ratios.csv")
    s=source_audit("results/source_data/Fig4_smooth_excitation_ratios.csv")
    measured_material=vacuum_D(M_POLE,Rmax)
    measured_vac=vacuum_D(VAC_POLE,Rmax)
    matter_contour=count_contour(M_POLE,n=n,Rmax=Rmax)
    if matter_contour["winding_number"]!=0:
        raise RuntimeError("Unexpected q=1 vacuum pole in q=.1 matter disk")
    if abs(measured_vac)>=abs(measured_material):
        raise RuntimeError("Vacuum known QNM determinant is not smaller at benchmark")
    source_ratios={}
    for key,item in (("point",p),("smooth",s)):
        source_ratios[key]={
            "q0p1":item[0],"q0p2":item[1],
            "q0p2_to_q0p1_median_ratio":float(
                item[1]["vacuum_normalized_residue_ratio_median"]/
                item[0]["vacuum_normalized_residue_ratio_median"])}
    data={
        "status":"q=1 exact Schwarzschild vacuum endpoint PASS; q>.2 self-consistent material poles NOT COMPUTED",
        "last_published_finite_q":max(Q_CONFIRMED),
        "endpoint_q":1.,
        "matter_pole_q0p1": [M_POLE.real,M_POLE.imag],
        "vacuum_grav220":[VAC_POLE.real,VAC_POLE.imag],
        "q1_transfer":transfer,
        "q1_normalized_vacuum_det_at_matter_pole":[measured_material.real,measured_material.imag],
        "q1_normalized_vacuum_det_at_vacuum_grav220":[measured_vac.real,measured_vac.imag],
        "q1_matter_pole_disk":matter_contour,
        "q1_exact_matter_branch_residue_at_fixed_old_pole":"NO distinct pole in tested q=1 vacuum disk",
        "published_q_family":table,
        "published_exterior_source_residue_ratios":source_ratios,
        "unresolved_mass_ratios":[{"q":q,"ocean_length_M":thickness(q),
                                   "mass_difference_fraction":1.-q,
                                   "pole_status":"NOT COMPUTED",
                                   "residue_status":"NOT COMPUTED"}
                                  for q in UNRESOLVED_Q],
        "physics_statement":{
            "conditional":"If T_q(omega)->I uniformly on compact complex neighborhoods away from vacuum poles, Hurwitz/Rouche implies no separate matter-only roots remain in such neighborhoods for q sufficiently near one.",
            "not_established":"Uniform bounded susceptibility/transfer generator and a self-consistent spectral family q>.20 have NOT been verified. Poles may move, merge, collide with a cut, lose observable residues or follow a singular limit.",
            "not_the_same_as":"q=1 vacuum endpoint does NOT numerically prove how finite-q matter-pole residue tends to zero or whether a hidden mode remains with vanishing gravitational coupling.",
        },
        "required_for_actual_q_to_1_continuation":[
            "Rebuild stationary Einstein-Vlasov atmosphere for each NEW q > .2; do not recycle the q=.1 spectral kernel",
            "Regenerate full radial and carrier action-angle spectral catalogs and validate stress-energy conservation",
            "Root-track the retarded global determinant with seed continuity, winding and step/multiprecision checks",
            "Compute dW/domega, its error across differencing scales, and source-normalized gravitational Green-function residue for unchanged exterior sources",
            "Measure continuum/cut stability and q-dependent source burst response; compare against Schwarzschild vacuum at q=1"
        ],
        "warnings":[
            "The mass-difference fraction 1-q is not necessarily a matter proper rest-mass fraction",
            "Original q=.05-.2 environment corresponds to 80-95% gravitational mass difference; no dilute-environment result exists",
            "Source residue ratios are not the initial amplitude of an astrophysical merger or detector waveform",
            "No detector strain is analyzed and no astrophysical detectability claim follows"
        ],
    }
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    (output/"q_to_one_endpoint_audit.json").write_text(json.dumps(data,indent=2)+"\n")
    outrows=table+data["unresolved_mass_ratios"]+[{"q":1.,"ocean_transfer_length_M":0.,
                         "pole_status":"EXACT VACUUM ENDPOINT",
                         "residue_status":"old matter pole absent in tested vacuum disk"}]
    with (output/"q_family_measured_vs_unresolved.csv").open("w",newline="") as f:
        fields=list(dict.fromkeys(k for r in outrows for k in r))
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader();w.writerows(outrows)
    print("Q_TO_ONE_VACUUM_GATE_PASS",json.dumps({
        "q1_vacuum_det_at_matter":abs(measured_material),
        "q1_vacuum_det_at_grav220":abs(measured_vac),
        "matter_disk_winding":matter_contour["winding_number"],
        "q_poles_measured":len(table),
        "q_poles_near_one_new":0,
        "point_q01_range":[p[0]["vacuum_normalized_residue_ratio_min"],p[0]["vacuum_normalized_residue_ratio_max"]],
        "smooth_q02_range":[s[1]["vacuum_normalized_residue_ratio_min"],s[1]["vacuum_normalized_residue_ratio_max"]]
    },separators=(",",":")),flush=True)
    return data


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default=str(ROOT/"results/q_to_one_endpoint"))
    parser.add_argument("--contour-n",type=int,default=CONTOUR_N)
    parser.add_argument("--Rmax",type=float,default=RMAX)
    args=parser.parse_args()
    run(args.output,n=args.contour_n,Rmax=args.Rmax)
