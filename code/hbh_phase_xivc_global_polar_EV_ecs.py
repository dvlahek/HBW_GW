#!/usr/bin/env python3
"""
HBH Phase XIVc — ECS-corrected global polar Einstein--Vlasov QNM solver
=======================================================================

Global Schwarzschild/Zerilli boundary solver used by the final smooth
Einstein--Vlasov continuation.  The q=1 Schwarzschild benchmark is treated as
a mandatory gate before target roots are interpreted.
"""

import argparse
import cmath
import json
import math
import shutil
import sys
import time
import zipfile
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import root

OMEGA0 = 1.0 / math.sqrt(27.0)
W_KIN = 2.0 * OMEGA0
W_POL = math.sqrt(5.0/27.0)


def mu2_exact(w):
    w2 = w*w
    num = -(
        (27*w2 - 5)
        *(
            1062882*w2*w2*w2
            - 295245*w2*w2
            + 20412*w2
            + 20
        )
    )
    den = 81*(
        354294*w2*w2*w2
        - 98415*w2*w2
        + 6264*w2
        + 100
    )
    return num/den


def lambdas_exact(w):
    mu = cmath.sqrt(mu2_exact(w))
    return (-1/9 + mu, -1/9 - mu)


def load_phase_ix():
    here = Path(__file__).resolve().parent
    path = here / "hbh_phase_ix_bulk_EV.py"
    if not path.exists():
        raise FileNotFoundError(
            "Put hbh_phase_ix_bulk_EV.py in the same folder."
        )

    spec = importlib.util.spec_from_file_location("phase_ix", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["phase_ix"] = mod
    spec.loader.exec_module(mod)
    return mod


def build_symbolic_ocean_matrix():
    mod = load_phase_ix()
    w, lam = sp.symbols("omega lambda")

    print("[setup] Building symbolic 7x4 Einstein--Vlasov ocean matrix...", flush=True)
    t0 = time.time()

    geom = mod.build_geometry(w, lam)
    stress = mod.build_stress(geom, w, lam)
    eq = mod.projected_equations(geom, stress)
    M = mod.coefficient_matrix(eq, geom["symbols"])

    print(f"[setup] symbolic matrix ready in {time.time()-t0:.1f} s", flush=True)

    Mfun = sp.lambdify((w, lam), M, "numpy")
    return Mfun


def nullvector_full_matrix(M):
    M = np.asarray(M, dtype=np.complex128)

    Mr = M.copy()
    for i in range(Mr.shape[0]):
        n = np.linalg.norm(Mr[i])
        if n > 0:
            Mr[i] /= n

    U, s, Vh = np.linalg.svd(Mr)
    v = Vh.conj().T[:, -1]

    k = int(np.argmax(np.abs(v)))
    if abs(v[k]) > 0:
        v = v / v[k]

    residual = np.linalg.norm(Mr @ v)
    return v, residual, float(s[-1]), float(s[0])


def orbital_edge_moments(w, lam, v):
    H0, H2, J, K = v
    S = H0 + K
    SX = lam*S

    wt = {-2: 3/8, 0: 1/4, 2: 3/8}

    xi = {}
    for n in (-2,0,2):
        dn = w - n*OMEGA0
        xi[n] = -SX/(2*dn*dn)

    xi_vec = -1j/6 * sum(n*wt[n]*xi[n] for n in wt)
    xi_tf = xi[0]/16 - (xi[2]+xi[-2])/32

    return xi, xi_vec, xi_tf


def interface_column(w, lam, v):
    H0, H2, J, K = v

    xi, xi_vec, xi_tf = orbital_edge_moments(w, lam, v)

    dVT_X = (4/81) * xi_vec
    dG_X = -(4/81) * xi_tf

    dK_gi = 3*dG_X
    dJ_gi = -3*math.sqrt(3)*dVT_X - 27j*w*dG_X

    K_eff = K + dK_gi
    J_eff = J + dJ_gi

    Z = K_eff - 1j*J_eff/(9*w)
    ZX = K_eff/27 + 26j*J_eff/(243*w)

    audit = {
        "H0": H0,
        "H2": H2,
        "J": J,
        "K": K,
        "S": H0+K,
        "SX": lam*(H0+K),
        "xi_vec": xi_vec,
        "xi_tf": xi_tf,
        "dVT_X": dVT_X,
        "dG_X": dG_X,
        "dK_gi": dK_gi,
        "dJ_gi": dJ_gi,
        "Z": Z,
        "ZX": ZX,
    }

    return np.array([Z,ZX], dtype=np.complex128), audit


def ocean_interface_basis(Mfun, w):
    lams = lambdas_exact(w)

    cols = []
    audits = []

    for branch, lam in zip(("plus","minus"), lams):
        M = Mfun(w, lam)
        v, res, smin, smax = nullvector_full_matrix(M)

        col, audit = interface_column(w, lam, v)

        audit.update({
            "branch": branch,
            "lambda": lam,
            "null_residual": res,
            "sv_min": smin,
            "sv_max": smax,
        })

        cols.append(col)
        audits.append(audit)

    C = np.column_stack(cols)

    return C, np.array(lams, dtype=np.complex128), audits


def ocean_transfer(Mfun, w, q):
    if q <= 0 or q > 1:
        raise ValueError("q must satisfy 0 < q <= 1")

    if abs(q-1.0) < 1e-15:
        return np.eye(2, dtype=np.complex128), {
            "C_cond": 1.0,
            "det_T": 1+0j,
            "expected_det": 1.0,
            "lambdas": [0j,0j],
            "audits": [],
        }

    C, lams, audits = ocean_interface_basis(Mfun, w)

    cond = np.linalg.cond(C)
    if not np.isfinite(cond) or cond > 1e13:
        raise FloatingPointError(
            f"Ill-conditioned interface basis C: cond={cond:.3e}"
        )

    L = 9*math.log(1/q)

    P = np.diag(np.exp(lams*L))
    T = C @ P @ np.linalg.inv(C)

    info = {
        "C_cond": float(cond),
        "det_T": np.linalg.det(T),
        "expected_det": q*q,
        "lambdas": lams,
        "audits": audits,
    }

    return T, info


def Vz(R):
    f = 1.0 - 2.0/R
    return (
        2*f
        *(
            12*R**3
            +12*R**2
            +18*R
            +9
        )
        /(R**3*(2*R+3)**2)
    )


def riccati_rhs(R, y, w):
    f = 1.0 - 2.0/R
    return [(Vz(R) - w*w - y[0]*y[0])/f]


def horizon_initial_y(w, eps):
    u = eps

    b1 = 27/(28*(1-4j*w))

    b2 = (
        -3j*(736*w*w + 305j*w - 91)
        /(196*(2*w+1j)*(4*w+1j)**2)
    )

    return -1j*w + b1*u + b2*u*u


def vacuum_inner_state(w, eps=1e-7, rtol=2e-10, atol=2e-12):
    R0 = 2.0 + eps
    y0 = horizon_initial_y(w, eps)

    sol = solve_ivp(
        lambda R, y: riccati_rhs(R, y, w),
        (R0, 3.0),
        [y0],
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=0.05,
    )

    if not sol.success:
        raise RuntimeError(
            "Inner Zerilli integration failed: " + sol.message
        )

    y3 = sol.y[0, -1]

    if not (
        np.isfinite(y3.real)
        and np.isfinite(y3.imag)
    ):
        raise FloatingPointError(
            "Inner Zerilli integration returned a non-finite state."
        )

    return np.array([1+0j, y3], dtype=np.complex128)


def infinity_initial_y_ecs(w, R, order=8):
    coeffs = {
        1: 0.0,
        2: -3j/w,
        3: 3*(4j*w-1)/w**2,
        4: 3*(-37j*w+32)/(4*w**2),
        5: 3*(42j*w**3-61*w**2-14j*w-6)/(2*w**4),
        6: 3*(-651j*w**4+1432*w**3+904j*w**2+296*w+144j)/(16*w**5),
        7: 3*(1197j*w**5-3633*w**4-3952j*w**3-172*w**2-1848j*w+288)/(16*w**6),
        8: 3*(-8505j*w**6+32382*w**5+53076j*w**4-18176*w**3+37736j*w**2-24816*w-1728j)/(64*w**7),
    }

    y = 1j*w
    for k in range(1, min(order, 8)+1):
        y += coeffs[k]/R**k
    return y


def vacuum_outer_state(
    w,
    Rmax=80.0,
    rtol=2e-10,
    atol=2e-12,
    asymptotic_order=8,
    theta_offset=0.0,
):
    theta = -cmath.phase(w) + theta_offset
    ray = cmath.exp(1j*theta)

    Smax = float(Rmax)
    R0 = 3.0 + Smax*ray
    y0 = infinity_initial_y_ecs(w, R0, order=asymptotic_order)

    def rhs(s, y):
        R = 3.0 + s*ray
        f = 1.0 - 2.0/R
        return [ray*(Vz(R) - w*w - y[0]*y[0])/f]

    sol = solve_ivp(
        rhs,
        (Smax, 0.0),
        [y0],
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=1.0,
    )

    if not sol.success:
        raise RuntimeError(
            "Outer ECS Zerilli integration failed: " + sol.message
        )

    y3 = sol.y[0,-1]
    return np.array([1+0j, y3], dtype=np.complex128)


class GlobalDeterminant:
    def __init__(
        self,
        Mfun,
        q,
        Rmax=80.0,
        horizon_eps=1e-7,
        rtol=2e-10,
        atol=2e-12,
    ):
        self.Mfun = Mfun
        self.q = q
        self.Rmax = Rmax
        self.horizon_eps = horizon_eps
        self.rtol = rtol
        self.atol = atol
        self.last_info = None

    def raw(self, w):
        zin = vacuum_inner_state(
            w,
            eps=self.horizon_eps,
            rtol=self.rtol,
            atol=self.atol,
        )

        zout = vacuum_outer_state(
            w,
            Rmax=self.Rmax,
            rtol=self.rtol,
            atol=self.atol,
        )

        T, info = ocean_transfer(self.Mfun, w, self.q)
        zpred = T @ zin

        D = zpred[0]*zout[1] - zpred[1]*zout[0]

        scale = max(
            np.linalg.norm(zpred)*np.linalg.norm(zout),
            1e-300,
        )
        Dn = D/scale

        info.update({
            "z_inner": zin,
            "z_outer": zout,
            "z_pred": zpred,
            "raw_det": D,
            "normalized_det": Dn,
        })

        self.last_info = info
        return Dn

    def real2(self, x):
        w = complex(float(x[0]), float(x[1]))

        try:
            D = self.raw(w)
            if not (np.isfinite(D.real) and np.isfinite(D.imag)):
                return np.array([1e6,1e6])
            return np.array([D.real,D.imag])
        except Exception as exc:
            self.last_exception = repr(exc)
            return np.array([1e6,1e6])


def solve_seed(detobj, seed, maxfev=120):
    ans = root(
        detobj.real2,
        [seed.real, seed.imag],
        method="hybr",
        options={"maxfev":maxfev, "xtol":1e-10},
    )

    w = complex(ans.x[0], ans.x[1])

    try:
        D = detobj.raw(w)
        absD = abs(D)
        info = detobj.last_info
    except Exception as exc:
        D = complex(np.nan,np.nan)
        absD = np.inf
        info = None
        final_exception = repr(exc)
    else:
        final_exception = ""

    return {
        "seed_re": seed.real,
        "seed_im": seed.imag,
        "root_re": w.real,
        "root_im": w.imag,
        "solver_success": bool(ans.success),
        "message": str(ans.message),
        "final_exception": final_exception,
        "abs_normalized_det": absD,
        "C_cond": info["C_cond"] if info is not None and "C_cond" in info else np.nan,
        "detT_re": info["det_T"].real if info is not None and "det_T" in info else np.nan,
        "detT_im": info["det_T"].imag if info is not None and "det_T" in info else np.nan,
    }


def dedupe_roots(df, tol=2e-6):
    if len(df)==0:
        return df

    df = df.sort_values("abs_normalized_det")
    kept=[]

    for _,r in df.iterrows():
        z=complex(r.root_re,r.root_im)

        if any(abs(z-complex(q["root_re"],q["root_im"]))<tol for q in kept):
            continue

        kept.append(r.to_dict())

    return pd.DataFrame(kept)


def stability_check(Mfun, q, w):
    variants = [
        ("base",       80.0, 1e-7, 8,  0.00),
        ("ecs_60",     60.0, 1e-7, 8,  0.00),
        ("ecs_100",   100.0, 1e-7, 8,  0.00),
        ("ecs_140",   140.0, 1e-7, 8,  0.00),
        ("series_6",   80.0, 1e-7, 6,  0.00),
        ("angle_p003", 80.0, 1e-7, 8, +0.03),
        ("angle_m003", 80.0, 1e-7, 8, -0.03),
        ("eps_1e-6",   80.0, 1e-6, 8,  0.00),
        ("eps_1e-8",   80.0, 1e-8, 8,  0.00),
    ]

    rows=[]

    for name,Smax,eps,order,angle in variants:
        zin = vacuum_inner_state(
            w,
            eps=eps,
            rtol=2e-10,
            atol=2e-12,
        )

        zout = vacuum_outer_state(
            w,
            Rmax=Smax,
            rtol=2e-10,
            atol=2e-12,
            asymptotic_order=order,
            theta_offset=angle,
        )

        T, info = ocean_transfer(Mfun, w, q)
        zpred = T @ zin

        D = zpred[0]*zout[1] - zpred[1]*zout[0]
        scale = max(
            np.linalg.norm(zpred)*np.linalg.norm(zout),
            1e-300,
        )
        Dn = D/scale

        rows.append({
            "variant":name,
            "ecs_length":Smax,
            "horizon_eps":eps,
            "asymptotic_order":order,
            "theta_offset":angle,
            "det_re":Dn.real,
            "det_im":Dn.imag,
            "abs_det":abs(Dn),
            "C_cond":info.get("C_cond",np.nan),
            "detT_error": abs(info.get("det_T",np.nan)-q*q) if q<1 else 0.0,
        })

    return pd.DataFrame(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--q",type=float,default=0.1)
    ap.add_argument("--Rmax",type=float,default=80.0)
    ap.add_argument("--horizon-eps",type=float,default=1e-7)
    args = ap.parse_args()

    Mfun = build_symbolic_ocean_matrix()
    schw_ref = 0.3736716844-0.0889623157j
    det = GlobalDeterminant(Mfun, args.q, Rmax=args.Rmax, horizon_eps=args.horizon_eps)
    result = solve_seed(det, schw_ref if abs(args.q-1.0)<1e-15 else 0.37-0.01j)
    print(result)
