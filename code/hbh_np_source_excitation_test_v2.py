#!/usr/bin/env python3
"""Retarded Green-function source-coupling test for the matter-led QNM.

The source is placed entirely in the outer Schwarzschild vacuum.  The
calculation reports pole residues for point, Gaussian, and derivative-Gaussian
polar sources.  It is a linear source-coupling test, not a detector forecast.
"""

import argparse
import io
import math
import shutil
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

import hbh_phase_xxic8_landau_continuation as c8
import hbh_phase_xxic9_final_global_smooth_qnm_v4 as c9
from publication_inputs import load_publication_case as load_case_inputs

Q_CASES = [0.10, 0.20]
POINT_SOURCE_RADII = [3.5, 4.0, 5.0, 6.0, 8.0, 10.0]
R_SOURCE_MIN = 3.25
R_SOURCE_MAX = 12.0
R_OBS = 25.0

SMOOTH_SOURCES = [
    {"name":"gaussian_R5_sigmaX0p75","kind":"gaussian","R_center":5.0,"sigma_X":0.75},
    {"name":"gaussian_R8_sigmaX1p50","kind":"gaussian","R_center":8.0,"sigma_X":1.50},
    {"name":"derivative_gaussian_R6_sigmaX1p00","kind":"derivative_gaussian","R_center":6.0,"sigma_X":1.00},
]


def load_publication_case(family_zip, q):
    """Return the published root and kinetic nodes for a selected mass ratio."""
    radial, carrier, row = load_case_inputs(q, family_zip=family_zip)
    w = complex(float(row["root_re"]), float(row["root_im"]))
    lams = [
        complex(float(row["lambda1_re"]), float(row["lambda1_im"])),
        complex(float(row["lambda2_re"]), float(row["lambda2_im"])),
    ]
    return radial, carrier, w, lams


def tortoise_X(R):
    R = np.asarray(R, dtype=float)
    if np.any(R <= 2.0):
        raise ValueError("R must be larger than 2M.")
    return R + 2.0*np.log(R/2.0 - 1.0)


def propagate_vacuum_state(xivc, w, z3, R_end, *, dense=False):
    z3 = np.asarray(z3, dtype=np.complex128)

    def rhs(R, y):
        f = 1.0 - 2.0/R
        Z, ZX = y
        return np.asarray([ZX/f, (xivc.Vz(R)-w*w)*Z/f], dtype=np.complex128)

    sol = solve_ivp(
        rhs,
        (3.0, float(R_end)),
        z3,
        method="DOP853",
        rtol=2e-10,
        atol=2e-12,
        max_step=0.08,
        dense_output=dense,
    )
    if not sol.success:
        raise RuntimeError("Outer-vacuum propagation failed: " + sol.message)
    return sol


def edge_states(xivc, Mfun, row_names, q, w, lams, radial, carrier):
    old = c9.ACTIVE_LAMBDA_SEED_PAIR
    c9.ACTIVE_LAMBDA_SEED_PAIR = tuple(lams)
    try:
        zin = xivc.vacuum_inner_state(w, eps=1e-7, rtol=2e-10, atol=2e-12)
        T, info = c9.ocean_transfer(
            xivc, Mfun, row_names, w, radial, carrier, q, tau=1.0, edge_t=1.0
        )
        zL3 = T @ zin
        zR3 = xivc.vacuum_outer_state(
            w,
            Rmax=80.0,
            rtol=2e-10,
            atol=2e-12,
            asymptotic_order=8,
            theta_offset=0.0,
        )
        W = zL3[0]*zR3[1] - zL3[1]*zR3[0]
        return zL3, zR3, W, info
    finally:
        c9.ACTIVE_LAMBDA_SEED_PAIR = old


def wronskian(xivc, Mfun, row_names, q, w, lams, radial, carrier):
    return edge_states(xivc, Mfun, row_names, q, w, lams, radial, carrier)[2]


def wronskian_derivative(xivc, Mfun, row_names, q, w, lams, radial, carrier, h):
    Wp = wronskian(xivc, Mfun, row_names, q, w+h, lams, radial, carrier)
    Wm = wronskian(xivc, Mfun, row_names, q, w-h, lams, radial, carrier)
    return (Wp-Wm)/(2*h)


def smooth_source_profile(spec, Rgrid):
    X = tortoise_X(Rgrid)
    Xc = float(tortoise_X(spec["R_center"]))
    sig = float(spec["sigma_X"])
    u = (X-Xc)/sig

    if spec["kind"] == "gaussian":
        raw = np.exp(-0.5*u*u)
    elif spec["kind"] == "derivative_gaussian":
        raw = u*np.exp(-0.5*u*u)
    else:
        raise ValueError(spec["kind"])

    norm = math.sqrt(max(float(np.trapezoid(np.abs(raw)**2, X)), 1e-300))
    return X, raw/norm


def run_case(q, radial, carrier, w, lams, *, xivc, Mfun, row_names):
    zL3, zR3, W0, info = edge_states(xivc, Mfun, row_names, q, w, lams, radial, carrier)
    dW1 = wronskian_derivative(xivc, Mfun, row_names, q, w, lams, radial, carrier, 2e-5)
    dW2 = wronskian_derivative(xivc, Mfun, row_names, q, w, lams, radial, carrier, 1e-5)
    dW_rel = abs(dW2-dW1)/max(abs(dW2),1e-300)

    solL = propagate_vacuum_state(xivc, w, zL3, R_OBS, dense=True)
    solR = propagate_vacuum_state(xivc, w, zR3, R_OBS, dense=True)
    ZR_obs = complex(solR.sol(R_OBS)[0])

    point_rows = []
    for Rs in POINT_SOURCE_RADII:
        ZLs = complex(solL.sol(Rs)[0])
        res = ZLs*ZR_obs/dW2
        point_rows.append({
            "q":q,
            "R_source":Rs,
            "R_observer":R_OBS,
            "residue_re":res.real,
            "residue_im":res.imag,
            "residue_abs":abs(res),
            "derivative_relative_difference":dW_rel,
        })

    Rgrid = np.linspace(R_SOURCE_MIN, R_SOURCE_MAX, 700)
    Xgrid = tortoise_X(Rgrid)
    ZL_grid = np.asarray(solL.sol(Rgrid)[0], dtype=np.complex128)
    ZL_norm = math.sqrt(max(float(np.trapezoid(np.abs(ZL_grid)**2, Xgrid)),1e-300))

    smooth_rows = []
    for spec in SMOOTH_SOURCES:
        X,S = smooth_source_profile(spec, Rgrid)
        coupling = np.trapezoid(S*ZL_grid, X)
        overlap = abs(coupling)/max(ZL_norm,1e-300)
        signed_area = np.trapezoid(S, X)
        res = ZR_obs*coupling/dW2

        smooth_rows.append({
            "q":q,
            "source_name":spec["name"],
            "source_kind":spec["kind"],
            "R_center":spec["R_center"],
            "sigma_X":spec["sigma_X"],
            "normalized_overlap":overlap,
            "signed_source_area":signed_area,
            "residue_re":res.real,
            "residue_im":res.imag,
            "residue_abs":abs(res),
            "derivative_relative_difference":dW_rel,
        })

    summary = {
        "q":q,
        "root_re":w.real,
        "root_im":w.imag,
        "dW_abs":abs(dW2),
        "dW_derivative_relative_difference":dW_rel,
        "point_residue_min_abs":min(r["residue_abs"] for r in point_rows),
        "point_residue_max_abs":max(r["residue_abs"] for r in point_rows),
        "smooth_residue_min_abs":min(r["residue_abs"] for r in smooth_rows),
        "smooth_residue_max_abs":max(r["residue_abs"] for r in smooth_rows),
        "smooth_overlap_min":min(r["normalized_overlap"] for r in smooth_rows),
        "smooth_overlap_max":max(r["normalized_overlap"] for r in smooth_rows),
    }

    return summary, point_rows, smooth_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", default=None, help="Optional legacy full-family ZIP; without it only the public q=0.10 NPZ is available")
    ap.add_argument("--q", type=float, nargs="+", choices=Q_CASES, default=[0.10], help="Mass ratios to recompute; q=0.20 needs --family with its spectral nodes")
    ap.add_argument("--outdir", default="source_excitation_output")
    args = ap.parse_args()

    out = Path(args.outdir)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    xivc = c9.load_xivc()
    phase_ix = xivc.load_phase_ix()
    Mfun,row_names = c9.build_generic_bulk_matrix(phase_ix)

    summaries=[]
    point=[]
    smooth=[]

    for q in args.q:
        radial,carrier,w,lams = load_publication_case(args.family,q)
        s,p,g = run_case(q,radial,carrier,w,lams,xivc=xivc,Mfun=Mfun,row_names=row_names)
        summaries.append(s); point.extend(p); smooth.extend(g)

    pd.DataFrame(summaries).to_csv(out/"source_excitation_summary.csv",index=False)
    pd.DataFrame(point).to_csv(out/"point_source_residues.csv",index=False)
    pd.DataFrame(smooth).to_csv(out/"smooth_source_residues.csv",index=False)
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
