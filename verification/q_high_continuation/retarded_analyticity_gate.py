#!/usr/bin/env python3
"""Check actual C8 second-sheet analyticity before high-q global root fitting.

The original C8 retarded branch picks a real-frequency interpolation segment
via searchsorted(Re omega). At complex omega this piecewise polynomial
continuation is NOT continuous when Re omega crosses a knot:
rho_left(omega) != rho_right(omega) for nonzero Im omega.
This matters for global root tracking and especially dW/domega residues.
The exact piecewise-linear branch integrals are evaluated with numerically
stable midpoint formulas, and real versus imaginary finite-difference
derivatives must agree before we allow a new high-q mode claim.

This is an *analytical/numerical diagnostic*, not a new physical continuation.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"code"))
import hbh_phase_xxic8_landau_continuation as c8

REFERENCE_Q020=0.3797013686231943-0.0239854837986382j
REFERENCE_Q010=0.3775148159474828-0.020871633919164j


class PackedBranch:
    def __init__(self,branches):
        self.a=np.concatenate([b["a"] for b in branches])
        self.b=np.concatenate([b["b"] for b in branches])
        self.m=np.concatenate([b["m"] for b in branches])
        self.c=np.concatenate([b["c"] for b in branches])
        self.x=[b["x"] for b in branches]
        self.xlo=np.array([z[0] for z in self.x])
        self.xhi=np.array([z[-1] for z in self.x])
        self.starts=np.r_[0,np.cumsum([len(b["a"]) for b in branches])]

    def at(self,z):
        z=complex(z)
        width=self.b-self.a
        midpoint=(self.a+self.b)*.5
        denom=z-midpoint
        u=width/(2*denom)
        if np.max(np.abs(u))>.45:
            raise RuntimeError("Unresolved spectral segment too near complex pole")
        u2=u*u
        summ=np.zeros_like(u,dtype=np.complex128)
        v=u2.copy()
        for n in range(1,21):
            summ+=v/(2*n+1)
            v*=u2
        L=width/denom*(1+summ)
        rho_mid=self.m*midpoint+self.c
        val=np.sum(rho_mid*L+self.m*width*summ)
        if z.imag<0:
            for i in np.flatnonzero((self.xlo<=z.real)&(z.real<=self.xhi)):
                grid=self.x[i]
                j=max(0,min(len(grid)-2,int(np.searchsorted(grid,z.real,side="right")-1)))
                k=self.starts[i]+j
                val-=2j*math.pi*(self.m[k]*z+self.c[k])
        return complex(val)

    def knots_within(self,z,h):
        xr=z.real
        return sum(int(np.count_nonzero(np.abs(grid-xr)<h))
                   for grid in self.x)


def pack_catalog(radial,carrier):
    return (PackedBranch([r["JJ"] for r in radial]),
            PackedBranch([r["JX"] for r in radial]),
            PackedBranch([r["HH"] for r in carrier]))


def stable_chi(z,packed):
    JJ,JX,H=(b.at(z) for b in packed)
    O=1./math.sqrt(27.)
    k=.5*((9*O/8)/(2*O-z))/H
    return {"chiJ":k*JJ,"chiX":k*JX,"PH":H}


def complex_pair(z):
    return [float(z.real),float(z.imag)]


def derivative_analyticity_gate(packed,z,large_h=1e-5,small_h=2e-7):
    """Compare valid within-cell holomorphic derivatives and old FD step.

    One-sided sheet selection is FROZEN by the complex observation point.
    If the first sheet changes as Re omega crosses nearby knots, real-axis
    derivatives are not complex derivatives; refuse residue extraction.
    """
    base=stable_chi(z,packed)
    result={"reference_complex_Momega":complex_pair(z),
            "radial_grid_knots_within_2e5":packed[0].knots_within(z,2e-5),
            "radial_grid_knots_within_1e5":packed[0].knots_within(z,1e-5),
            "per_field":{}}
    for name in ("chiJ","chiX"):
        vals={}
        for h in (large_h,small_h):
            vrp=stable_chi(z+h,packed)[name]
            vrm=stable_chi(z-h,packed)[name]
            vip=stable_chi(z+1j*h,packed)[name]
            vim=stable_chi(z-1j*h,packed)[name]
            dre=(vrp-vrm)/(2*h)
            dim=(vip-vim)/(2j*h)
            vals[str(h)]={
                "derivative_real_direction":complex_pair(dre),
                "derivative_imag_direction":complex_pair(dim),
                "relative_CR_disagreement":float(abs(dre-dim)/max(abs(dim),1e-300)),
            }
        result["per_field"][name]={"value":complex_pair(base[name]),
                                   "derivative_audits":vals}
    mismatch=max(result["per_field"][name]["derivative_audits"][str(large_h)]["relative_CR_disagreement"]
                 for name in ("chiJ","chiX"))
    # We require both numerical analyticity and a useful derivative
    # across the scale previously used for q0.10 source residues.
    result["max_relative_CR_disagreement_at_1e5"]=mismatch
    result["global_QNM_residue_gate_pass"]=bool(mismatch<=.03)
    result["status"]="PASS" if mismatch<=.03 else "BLOCKED_NONANALYTIC_C8_RETARDED_INTERPOLATION"
    result["qnm_pole_computed"]=False
    result["gravitational_residue_computed"]=False
    return result


def run(manifest=None,reference_q010=False,out=None,strict=False):
    if reference_q010:
        if manifest is not None:
            raise ValueError("Choose either a new-q manifest or q0.10 reference")
        p=ROOT/"verification/blind_multiple_shooting/blind_q010_spectral_nodes.npz"
        with np.load(p,allow_pickle=False) as f:
            arrays={k:np.array(f[k]) for k in f.files}
        z=REFERENCE_Q010
    else:
        if manifest is None:raise ValueError("New q requires a valid actual physical input manifest")
        from run_new_q import validate_manifest,load_spectral_arrays
        q,m,b,outputs=validate_manifest(manifest)
        arrays=load_spectral_arrays(outputs["spectral_nodes"],int(b["shell_count"]))
        z=REFERENCE_Q020
    radial,carrier=c8.build_catalog(arrays)
    packed=pack_catalog(radial,carrier)
    r=derivative_analyticity_gate(packed,z)
    r["q"]=.1 if reference_q010 else float(q)
    r["spectral_radial_branches"]=len(radial)
    r["spectral_carrier_branches"]=len(carrier)
    if out:
        path=Path(out)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(r,indent=2)+"\n")
    print("Q_HIGH_RETARDED_ANALYTICITY_GATE",json.dumps({
        "q":r["q"],"status":r["status"],
        "knots_near_reference":r["radial_grid_knots_within_2e5"],
        "CR_relative_disagreement":r["max_relative_CR_disagreement_at_1e5"]
    },separators=(",",":")),flush=True)
    if strict and not r["global_QNM_residue_gate_pass"]:
        raise RuntimeError("C8 nonholomorphic knot-dependent retarded continuation: NO reliable q>0.2 QNM/residue")
    return r


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path)
    p.add_argument("--reference-q010",action="store_true")
    p.add_argument("--output",type=Path)
    p.add_argument("--strict",action="store_true")
    a=p.parse_args()
    run(a.manifest,a.reference_q010,a.output,a.strict)
