#!/usr/bin/env python3
"""Analytic-in-frequency full-state q=1 Zerilli Wronskian zero-count audit.

Unlike the Riccati logarithmic derivative, this propagates the two
LINEAR fundamental solution states (Z,dZ/dX), so their unnormalized
Wronskian is analytic in complex w on a fixed ECS integration ray
(within the convergence/regularity domain). No quotient Z'/Z is used.

This evaluates the q=1 SCHWARZSCHILD vacuum endpoint only.
It is not a q=.20->1 Einstein-Vlasov pole continuation and cannot
infer the changing kinetic residue. Numerical winding is not a
rigorous interval-arithmetic certificate.
"""
from __future__ import annotations

import argparse
import cmath
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"code"))
import hbh_phase_xivc_global_polar_EV_ecs as v

MATTER=0.377514814460-0.020871654486j
SCHWARZSCHILD_220=0.373671684352-0.088962315642j
EPS=1e-7
RMAX=80.
N=64
RADIUS=.01


def full_state_wronskian(w,theta,Rmax=RMAX):
    """Unnormalized analytic full-state W at match radius R=3.

    Fixed ECS ray angle theta is independent of the trial frequency
    across a contour, and boundary-series normalizations Z=1
    are also fixed analytically. Therefore no movable Riccati
    coordinate poles enter the Wronskian.
    """
    w=complex(w)
    ray=cmath.exp(1j*float(theta))
    r0=3.+float(Rmax)*ray
    def rhs_R(R,z):
        ff=1.-2./R
        return [z[1]/ff,(v.Vz(R)-w*w)*z[0]/ff]
    inner=solve_ivp(
        rhs_R,(2.+EPS,3.),[1.+0j,v.horizon_initial_y(w,EPS)],
        method="DOP853",rtol=2e-10,atol=2e-12,max_step=.05)
    if not inner.success:
        raise RuntimeError("Failed direct full-state ingoing integration: "+inner.message)

    def rhs_s(s,z):
        R=3.+s*ray
        ff=1.-2./R
        return [ray*z[1]/ff,ray*(v.Vz(R)-w*w)*z[0]/ff]
    outer=solve_ivp(
        rhs_s,(float(Rmax),0.),[1.+0j,v.infinity_initial_y_ecs(w,r0,order=8)],
        method="DOP853",rtol=2e-10,atol=2e-12,max_step=1.)
    if not outer.success:
        raise RuntimeError("Failed direct full-state outgoing integration: "+outer.message)
    zi=inner.y[:,-1]
    zo=outer.y[:,-1]
    W=zi[0]*zo[1]-zi[1]*zo[0]
    if not np.isfinite(W.real) or not np.isfinite(W.imag):
        raise RuntimeError("Nonfinite direct full-state Wronskian")
    scale=max(np.linalg.norm(zi)*np.linalg.norm(zo),1e-300)
    return complex(W),float(abs(W)/scale)


def analytic_zero_count(center=MATTER,radius=RADIUS,n=N,Rmax=RMAX):
    """Net contour winding of analytic full-state W at fixed ECS angle."""
    theta=float(-cmath.phase(center))
    ts=np.linspace(0.,2.*math.pi,n+1)
    seq=[full_state_wronskian(center+radius*np.exp(1j*t),
                              theta,Rmax)[0] for t in ts]
    values=np.asarray(seq)
    if not np.all(np.isfinite(values)) or np.min(np.abs(values))<1e-10:
        raise RuntimeError("Near-zero or invalid W on analytic contour")
    increments=np.angle(values[1:]/values[:-1])
    total=float(np.sum(increments)/(2.*math.pi))
    if abs(total-round(total))>.05 or np.max(np.abs(increments))>math.pi*.8:
        raise RuntimeError("Unresolved analytic W argument winding")
    return {
        "center_complex":[center.real,center.imag],
        "radius":float(radius),"N":int(n),
        "fixed_exterior_complex_scaling_theta_rad":theta,
        "Rmax":float(Rmax),
        "analytic_full_state_W_winding":int(round(total)),
        "unrounded":total,
        "min_abs_W_on_contour":float(np.min(np.abs(values))),
        "max_argument_step":float(np.max(np.abs(increments))),
        "boundary_normalization":"analytic fixed Z=1 at horizon and ECS outer radius",
        "riccati_pole_artifacts_possible":False,
    }


def run(output_dir):
    theta=-cmath.phase(MATTER)
    w_mat,norm_mat=full_state_wronskian(MATTER,theta)
    # Separate outgoing ray centered on the known q=1 vacuum QNM.
    w_vac,norm_vac=full_state_wronskian(
        SCHWARZSCHILD_220,-cmath.phase(SCHWARZSCHILD_220))
    if norm_mat<.05:
        raise RuntimeError("Old matter pole looks like q=1 vacuum root in full-state computation")
    if norm_vac>1e-6:
        raise RuntimeError("Known Schwarzschild fundamental failed full-state q=1 benchmark")
    c0=analytic_zero_count(n=32,Rmax=80.)
    c1=analytic_zero_count(n=64,Rmax=100.)
    if c0["analytic_full_state_W_winding"]!=0 or c1["analytic_full_state_W_winding"]!=0:
        raise RuntimeError("A q=1 vacuum zero was detected near the q=.1 matter frequency")
    result={
        "scientific_scope":"q=1 vacuum FULL-STATE analytic W contour; not finite-q kinetic continuation",
        "full_state_q1_W_at_old_q0p1_matter_pole":[w_mat.real,w_mat.imag],
        "full_state_q1_normalized_det_abs_at_matter":norm_mat,
        "full_state_q1_W_at_known_vacuum_qnm":[w_vac.real,w_vac.imag],
        "full_state_q1_normalized_det_abs_at_vacuum":norm_vac,
        "analytic_contours":[c0,c1],
        "numerical_endpoint_result":"NO q=1 VACUUM QNM inside fixed radius .01 disk around published q=.1 matter pole, within validated fixed-ray analytic full-state numerical method",
        "new_finite_q_between_0p2_and_1_residue_computed":False,
        "new_finite_q_between_0p2_and_1_pole_computed":False,
        "no_detector_strain_used":True,
        "remaining_physics_question":"track self-consistent q>.2 kinetic spectrum and exterior Green-function residues"
    }
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    (out/"q1_analytic_full_state_zero_count.json").write_text(json.dumps(result,indent=2)+"\n")
    print("Q1_ANALYTIC_FULL_STATE_NO_MATTER_POLE_PASS",json.dumps({
        "at_old_matter_abs_D":norm_mat,"at_vacuum_220_abs_D":norm_vac,
        "winding_N32_R80":c0["analytic_full_state_W_winding"],
        "winding_N64_R100":c1["analytic_full_state_W_winding"]
    },separators=(",",":")),flush=True)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default=str(ROOT/"results/q_to_one_endpoint"))
    args=parser.parse_args()
    run(args.output)
