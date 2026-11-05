#!/usr/bin/env python3
"""Fixed-path analytic continuation of finite-resolution C8 spectral branches.

The original piecewise-linear Cauchy transform IS analytic off the real axis,
but the C8 second-sheet -2*pi*i*rho(z) correction changes its interpolant
according to Re(z). That makes the legacy numerical retarded continuation
nonholomorphic on vertical spectral-knot lines.

For each physical spectral branch, choose ONE real cut-crossing x_anchor.
Freeze which linear density segment is selected by that path. Evaluate the
analytic complex-linear density extension from that selected segment at ALL
z in the tested lower-half-plane. This is a bona fide holomorphic sheet of
the DISCRETE piecewise-linear representation. It is NOT proof of a unique,
correct continuum kinetic retarded sheet; changing anchor can change it.

No new QNM root/physical exterior residue may be claimed unless the response
and its derivatives remain stable across predeclared anchors, q / quadrature
resolutions and independent physical retarded representations.
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
from retarded_analyticity_gate import (
    pack_catalog,REFERENCE_Q010,REFERENCE_Q020,complex_pair,
)

ANCHOR_STRESS=(-0.0005,0.,0.0005)
MAX_REL_DERIVATIVE_VARIATION=.05
MAX_REL_CHI_VARIATION=.05
DIFF_STEP=1e-5


class FixedBranch:
    def __init__(self,packed,anchor):
        self.p=packed
        self.anchor=float(anchor)
        self.anchor_selected_branches=0
        slope=0j
        intercept=0j
        active=np.flatnonzero((packed.xlo<=anchor)&(anchor<=packed.xhi))
        for i in active:
            x=packed.x[i]
            j=max(0,min(len(x)-2,int(np.searchsorted(
                x,anchor,side="right")-1)))
            k=int(packed.starts[i]+j)
            slope+=packed.m[k]
            intercept+=packed.c[k]
        self.anchor_selected_branches=int(len(active))
        self.slope=slope
        self.intercept=intercept

    def direct_and_derivative(self,z):
        p=self.p
        width=p.b-p.a
        midpoint=(p.a+p.b)*.5
        denom=z-midpoint
        u=width/(2*denom)
        if np.max(np.abs(u))>.45:
            raise RuntimeError("Frequencies approach spectral endpoints; stable midpoint expansion invalid")
        u2=u*u
        series=np.zeros_like(u,dtype=np.complex128)
        derivative_series=np.zeros_like(u,dtype=np.complex128)
        power=u2.copy()
        for n in range(1,22):
            series+=power/(2*n+1)
            derivative_series+=2*n*power/(2*n+1)
            power*=u2
        rho_mid=p.m*midpoint+p.c
        integral=np.sum(rho_mid*width/denom*(1.+series)+p.m*width*series)
        derivative=np.sum(
            -rho_mid*width/((z-p.a)*(z-p.b))
            -p.m*width/denom*derivative_series
        )
        return complex(integral),complex(derivative)

    def value_and_derivative(self,z):
        z=complex(z)
        if z.imag>=0:
            raise ValueError("Frozen second-sheet calculations require Im(omega)<0")
        value,der=self.direct_and_derivative(z)
        value-=2j*math.pi*(self.slope*z+self.intercept)
        der-=2j*math.pi*self.slope
        return value,der


class FixedSheetResponse:
    def __init__(self,packed,anchor):
        self.branches=[FixedBranch(p,anchor) for p in packed]
        self.anchor=float(anchor)

    def evaluate(self,z):
        z=complex(z)
        (J,dJ),(X,dX),(H,dH)=[
            br.value_and_derivative(z) for br in self.branches]
        O=1./math.sqrt(27.)
        kap=.5*((9.*O/8.)/(2.*O-z))/H
        dk_over_k=1./(2.*O-z)-dH/H
        return {"chiJ":kap*J,
                "chiX":kap*X,
                "dchiJ":kap*(dJ+J*dk_over_k),
                "dchiX":kap*(dX+X*dk_over_k),
                "PH":H,
                "active_branch_counts":[b.anchor_selected_branches
                                         for b in self.branches]}


def evaluate_stress(packed,z,offsets=ANCHOR_STRESS,h=DIFF_STEP):
    z=complex(z)
    if not z.imag<0 or not np.all(np.isfinite([z.real,z.imag])):
        raise ValueError("Retarded seed must be finite and below real axis")
    report=[]
    for offset in offsets:
        obj=FixedSheetResponse(packed,z.real+offset)
        v=obj.evaluate(z)
        rp=obj.evaluate(z+h)
        rm=obj.evaluate(z-h)
        ip=obj.evaluate(z+1j*h)
        im=obj.evaluate(z-1j*h)
        per={}
        for k in ("chiJ","chiX"):
            analytic=v["d"+k]
            dreal=(rp[k]-rm[k])/(2*h)
            dimag=(ip[k]-im[k])/(2j*h)
            per[k]={
                "value":complex_pair(v[k]),
                "analytic_derivative":complex_pair(analytic),
                "CR_direction_mismatch":float(abs(dreal-dimag)/max(abs(analytic),1e-300)),
                "real_vs_analytic_derivative":float(abs(dreal-analytic)/max(abs(analytic),1e-300)),
            }
        report.append({
            "anchor_offset_Momega":float(offset),
            "anchor_real_Momega":float(obj.anchor),
            "active_JJ_JX_HH_branches":v["active_branch_counts"],
            "fields":per,
        })
    nominal=min(report,key=lambda x:abs(x["anchor_offset_Momega"]))
    for row in report:
        for k in ("chiJ","chiX"):
            own=row["fields"][k]
            ref=nominal["fields"][k]
            own["anchor_relative_value_shift"]=float(
                abs(complex(*own["value"])-complex(*ref["value"]))/
                max(abs(complex(*ref["value"])),1e-300))
            own["anchor_relative_derivative_shift"]=float(
                abs(complex(*own["analytic_derivative"])-
                    complex(*ref["analytic_derivative"]))/
                max(abs(complex(*ref["analytic_derivative"])),1e-300))
    maxval=max(x["fields"][k]["anchor_relative_value_shift"]
               for x in report for k in ("chiJ","chiX"))
    maxder=max(x["fields"][k]["anchor_relative_derivative_shift"]
               for x in report for k in ("chiJ","chiX"))
    maxcr=max(x["fields"][k]["CR_direction_mismatch"]
               for x in report for k in ("chiJ","chiX"))
    return {
        "seed_Momega":complex_pair(z),
        "anchor_stress":report,
        "max_relative_value_anchor_variation":float(maxval),
        "max_relative_derivative_anchor_variation":float(maxder),
        "max_Cauchy_Riemann_directional_mismatch":float(maxcr),
        "analyticity_on_each_frozen_sheet_pass":bool(maxcr<1e-4),
        "cross_anchor_physical_stability_pass":bool(
            maxval<=MAX_REL_CHI_VARIATION and
            maxder<=MAX_REL_DERIVATIVE_VARIATION),
        "global_QNM_claim_permitted":False,
        "physical_sheet_uniqueness_proven":False,
        "exterior_gravitational_residue_computed":False,
        "interpretation":"analytic discrete sheet is necessary, but is NOT sufficient for a physical q>0.2 matter QNM",
    }


def run(*,reference_q010=False,manifest=None,output=None,strict=False):
    if reference_q010==bool(manifest):
        raise ValueError("Choose exactly one of published q0.1 or new-q manifest")
    if reference_q010:
        with np.load(ROOT/"verification/blind_multiple_shooting/blind_q010_spectral_nodes.npz",
                     allow_pickle=False) as f:
            arrays={k:np.asarray(f[k]) for k in f.files}
        z=REFERENCE_Q010
        q=.1
    else:
        from run_new_q import validate_manifest,load_spectral_arrays
        q,j,b,files=validate_manifest(manifest)
        arrays=load_spectral_arrays(files["spectral_nodes"],int(b["shell_count"]))
        z=REFERENCE_Q020
    radial,carrier=c8.build_catalog(arrays)
    packed=pack_catalog(radial,carrier)
    report=evaluate_stress(packed,z)
    report["q"]=q
    report["source"]=("published_q010_reference" if reference_q010
                      else "new_sha_locked_einstein_vlasov_background")
    report["radial_branches"]=len(radial)
    report["carrier_branches"]=len(carrier)
    if output:
        path=Path(output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(report,indent=2)+"\n")
    print("Q_HIGH_FIXED_SHEET_PATH_AUDIT",json.dumps({
        "q":q,
        "analyticity":report["analyticity_on_each_frozen_sheet_pass"],
        "max_value_anchor_variation":report["max_relative_value_anchor_variation"],
        "max_derivative_anchor_variation":report["max_relative_derivative_anchor_variation"],
        "physical_sheet_stable":report["cross_anchor_physical_stability_pass"],
        "root_and_residue_not_computed":True,
    },separators=(",",":")),flush=True)
    if strict and not report["cross_anchor_physical_stability_pass"]:
        raise RuntimeError("Fixed sheets are individually analytic but physically sheet/anchor ambiguous")
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--reference-q010",action="store_true")
    p.add_argument("--manifest",type=Path)
    p.add_argument("--output",type=Path)
    p.add_argument("--strict",action="store_true")
    a=p.parse_args()
    run(reference_q010=a.reference_q010,manifest=a.manifest,
        output=a.output,strict=a.strict)
