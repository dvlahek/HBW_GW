#!/usr/bin/env python3
"""
HBH Phase IX local bulk Einstein--Vlasov derivation
===================================================

Purpose
-------
Build the full l=2 even-parity linearized Einstein equations directly in the
sharp HBH ocean, insert the conservation-closed collisionless source from
Phase VIII, and determine the allowed radial bulk exponents

    Y(X) ~ exp(lambda X - i omega T).

This is a BULK calculation.  It does NOT yet impose the inner/outer Zerilli
boundary conditions and therefore does not by itself produce a global QNM.

The script:
  1. constructs the HBH ocean metric
       ds^2 = F(X)[-dT^2+dX^2+27 dOmega^2],
       F=exp(2X/9)/3;
  2. constructs l=2 RW-gauge polar perturbations H0,H2,J,K;
  3. derives delta G_{mu nu} directly from the metric;
  4. inserts the Phase-VIII conservation-consistent Vlasov stress tensor;
  5. projects TT,TX,XX,Ttheta,Xtheta, angular-trace and angular-TF equations;
  6. builds the 7x4 algebraic coefficient matrix M(lambda,omega);
  7. obtains candidate lambda roots from independent 4x4 minors;
  8. validates every candidate against the FULL 7x4 system by smallest SVD
     singular value;
  9. writes CSV/TXT diagnostics and a ZIP.

Recommended first run
---------------------
py hbh_phase_ix_bulk_EV.py --omega-re 0.3849001795 --omega-im -0.005

Then try:
py hbh_phase_ix_bulk_EV.py --omega-re 0.3849001795 --omega-im -0.001
py hbh_phase_ix_bulk_EV.py --omega-re 0.4127 --omega-im -0.01

Requirements
------------
py -m pip install sympy numpy scipy pandas

This script is intentionally conservative: a candidate radial exponent is
accepted only if it nearly satisfies all seven projected Einstein equations.
"""

import argparse
import itertools
import json
import math
import shutil
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import sympy as sp


def log(msg):
    print(msg, flush=True)


def cplx(z):
    return complex(z)


def finite_complex(z):
    return np.isfinite(z.real) and np.isfinite(z.imag)


def numeric_matrix(Msym, lam_sym, lam_val, digits=30):
    out = np.empty(Msym.shape, dtype=np.complex128)
    for i in range(Msym.rows):
        for j in range(Msym.cols):
            z = Msym[i, j].subs(lam_sym, lam_val)
            out[i, j] = complex(sp.N(z, digits))
    return out


def smallest_sv(M):
    s = np.linalg.svd(M, compute_uv=False)
    return float(s[-1]), float(s[0])


def kinetic_coefficients(omega):
    """
    Exact l=2 collisionless coefficients.

    Definitions:
      S = H0 + K
      A = S_XX/2 + S_X/9
      C = J_X + 2J/9

      d = d_S S + d_A A + d_C C
      f = f_SX S_X
      e = e_S S + e_A A + e_C C
      p = p_SX S_X
      s = s_S S + s_A A + s_C C

    O^2=1/27.
    """
    O = sp.sqrt(sp.Rational(1, 27))
    w = omega

    dm = (2*O - w)
    dp = (2*O + w)

    return {
        "O": O,
        "d_S": 9*O**2 / (2*dm*dp),
        "d_A": (4*O**4 + O**2*w**2 + w**4) / (w**2 * dm**2 * dp**2),
        "d_C": sp.I*(O-w)*(O+w) / (w*dm*dp),
        "f_SX": sp.I*(O-w)*(O+w) / (2*w*dm*dp),
        "e_S": -3*sp.I*O*w / (4*dm*dp),
        "e_A": -sp.I*O*w / (dm**2*dp**2),
        "e_C": -O / (2*dm*dp),
        "p_SX": -O / (4*dm*dp),
        "s_S": -3*O**2 / (8*dm*dp),
        "s_A": O**2*(4*O**2 - 3*w**2) / (4*w**2*dm**2*dp**2),
        "s_C": sp.I*O**2 / (4*w*dm*dp),
    }


def build_geometry(omega, lam):
    """
    Returns:
      g0, h, ginv0, ginv1, Gamma0, Gamma1, R0, R1, G0, G1
    for the l=2,m=0 RW-gauge ansatz.

    We retain X,T,theta while differentiating, then evaluate at
      T=0, X=0, theta=pi/3.
    This removes almost all trigonometric algebra while retaining a generic
    point where Y and dY/dtheta are both nonzero.
    """

    T, X, th, ph = sp.symbols("T X theta phi", real=True)
    coords = [T, X, th, ph]

    H0, H2, J, K = sp.symbols("H0 H2 J K")

    F = sp.exp(sp.Rational(2, 9)*X) / 3
    Y = (3*sp.cos(th)**2 - 1) / 2
    E = sp.exp(lam*X - sp.I*omega*T)

    g0 = sp.diag(-F, F, 27*F, 27*F*sp.sin(th)**2)

    h = sp.zeros(4)
    h[0,0] = F*H0*E*Y
    h[0,1] = h[1,0] = F*J*E*Y
    h[1,1] = F*H2*E*Y
    h[2,2] = 27*F*K*E*Y
    h[3,3] = 27*F*sp.sin(th)**2*K*E*Y

    ginv0 = g0.inv()
    ginv1 = -ginv0*h*ginv0

    def der(expr, idx):
        return sp.diff(expr, coords[idx])

    C0 = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    C1 = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(4)]

    for d in range(4):
        for b in range(4):
            for c in range(4):
                C0[d][b][c] = der(g0[d,c], b) + der(g0[d,b], c) - der(g0[b,c], d)
                C1[d][b][c] = der(h[d,c], b) + der(h[d,b], c) - der(h[b,c], d)

    Gamma0 = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    Gamma1 = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(4)]

    for a in range(4):
        for b in range(4):
            for c in range(4):
                v0 = 0
                v1 = 0
                for d in range(4):
                    v0 += sp.Rational(1,2)*ginv0[a,d]*C0[d][b][c]
                    v1 += sp.Rational(1,2)*(ginv1[a,d]*C0[d][b][c] + ginv0[a,d]*C1[d][b][c])
                Gamma0[a][b][c] = sp.cancel(v0)
                Gamma1[a][b][c] = sp.expand(v1)

    R0 = sp.zeros(4)
    R1 = sp.zeros(4)

    log("  building Ricci tensor ...")
    for a in range(4):
        for b in range(4):
            v0 = 0
            v1 = 0

            for c in range(4):
                v0 += der(Gamma0[c][a][b], c) - der(Gamma0[c][a][c], b)
                v1 += der(Gamma1[c][a][b], c) - der(Gamma1[c][a][c], b)

                for d in range(4):
                    v0 += Gamma0[c][a][b]*Gamma0[d][c][d] - Gamma0[d][a][c]*Gamma0[c][b][d]
                    v1 += (
                        Gamma1[c][a][b]*Gamma0[d][c][d]
                        + Gamma0[c][a][b]*Gamma1[d][c][d]
                        - Gamma1[d][a][c]*Gamma0[c][b][d]
                        - Gamma0[d][a][c]*Gamma1[c][b][d]
                    )

            R0[a,b] = sp.cancel(v0)
            R1[a,b] = sp.expand(v1)

    eval_subs = {T: 0, X: 0, th: sp.pi/3, ph: 0}

    g0e = g0.subs(eval_subs)
    he = h.subs(eval_subs)
    ginv0e = ginv0.subs(eval_subs)
    ginv1e = ginv1.subs(eval_subs)
    R0e = R0.subs(eval_subs)
    R1e = R1.subs(eval_subs)

    Rscalar0 = sp.cancel(sum(ginv0e[a,b]*R0e[a,b] for a in range(4) for b in range(4)))

    Rscalar1 = sp.expand(sum(
        ginv1e[a,b]*R0e[a,b] + ginv0e[a,b]*R1e[a,b]
        for a in range(4) for b in range(4)
    ))

    G0 = sp.zeros(4)
    G1 = sp.zeros(4)

    for a in range(4):
        for b in range(4):
            G0[a,b] = sp.cancel(R0e[a,b] - sp.Rational(1,2)*g0e[a,b]*Rscalar0)
            G1[a,b] = sp.expand(
                R1e[a,b]
                - sp.Rational(1,2)*he[a,b]*Rscalar0
                - sp.Rational(1,2)*g0e[a,b]*Rscalar1
            )

    Y0 = sp.simplify(Y.subs(th, sp.pi/3))
    Yt = sp.simplify(sp.diff(Y, th).subs(th, sp.pi/3))
    Ytt = sp.simplify(sp.diff(Y, th, 2).subs(th, sp.pi/3))

    return {
        "symbols": (H0,H2,J,K),
        "g0": g0e,
        "h": he,
        "G0": G0,
        "G1": G1,
        "Y": Y0,
        "Ytheta": Yt,
        "Ythetatheta": Ytt,
        "theta": sp.pi/3,
    }


def build_stress(geom, omega, lam):
    H0,H2,J,K = geom["symbols"]
    theta = geom["theta"]
    Y = geom["Y"]
    Yt = geom["Ytheta"]
    Ytt = geom["Ythetatheta"]

    kin = kinetic_coefficients(omega)

    S = H0 + K
    SX = lam*S
    SXX = lam**2*S

    A = SXX/2 + SX/9
    C = lam*J + 2*J/9

    d = sp.expand(kin["d_S"]*S + kin["d_A"]*A + kin["d_C"]*C)
    f = sp.expand(kin["f_SX"]*SX)
    e = sp.expand(kin["e_S"]*S + kin["e_A"]*A + kin["e_C"]*C)
    p = sp.expand(kin["p_SX"]*SX)
    shear = sp.expand(kin["s_S"]*S + kin["s_A"]*A + kin["s_C"]*C)

    aamp = d + H0 - H2/2 - 3*K/2
    bamp = f - J
    Damp = d - H2/2 - 5*K/2

    pref = sp.Rational(1,36) / sp.pi

    s2 = sp.sin(theta)**2
    ct = sp.cos(theta)

    TFtt_cov = Ytt + 3*Y
    TFpp_cov = sp.sin(theta)*ct*Yt + 3*s2*Y

    TFtt_up = TFtt_cov
    TFpp_up = TFpp_cov / s2**2

    T1up = sp.zeros(4)

    T1up[0,0] = pref*aamp*Y
    T1up[0,1] = T1up[1,0] = pref*bamp*Y
    T1up[1,1] = 0

    O = kin["O"]
    T1up[0,2] = T1up[2,0] = pref*O*e*Yt
    T1up[1,2] = T1up[2,1] = pref*O*p*Yt

    T1up[2,2] = pref/sp.Integer(27) * (Damp*Y/2 + shear*TFtt_up)
    T1up[3,3] = pref/sp.Integer(27) * (Damp*Y/(2*s2) + shear*TFpp_up)

    T0up = sp.zeros(4)
    T0up[0,0] = pref
    T0up[2,2] = pref/sp.Integer(54)
    T0up[3,3] = pref/(sp.Integer(54)*s2)

    g0 = geom["g0"]
    h = geom["h"]

    T1cov = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            val = 0
            for c in range(4):
                for d in range(4):
                    val += g0[a,c]*g0[b,d]*T1up[c,d]
                    val += h[a,c]*g0[b,d]*T0up[c,d]
                    val += g0[a,c]*h[b,d]*T0up[c,d]
            T1cov[a,b] = sp.expand(val)

    return {
        "T1up": T1up,
        "T1cov": T1cov,
        "moments": {"d": d, "f": f, "e": e, "p": p, "s": shear, "A": A, "C": C}
    }


def projected_equations(geom, stress):
    G = geom["G1"]
    T = stress["T1cov"]
    Y = geom["Y"]
    Yt = geom["Ytheta"]
    th = geom["theta"]
    s2 = sp.sin(th)**2

    R = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            R[a,b] = sp.expand(G[a,b] - 8*sp.pi*T[a,b])

    eq = []
    eq.append(("TT", sp.cancel(R[0,0]/Y)))
    eq.append(("TX", sp.cancel(R[0,1]/Y)))
    eq.append(("XX", sp.cancel(R[1,1]/Y)))
    eq.append(("Ttheta", sp.cancel(R[0,2]/Yt)))
    eq.append(("Xtheta", sp.cancel(R[1,2]/Yt)))

    ang_trace = R[2,2] + R[3,3]/s2
    eq.append(("AngTrace", sp.cancel(ang_trace/Y)))

    ang_tf = R[2,2] - R[3,3]/s2
    TFtt = geom["Ythetatheta"] + 3*Y
    TFpp_over_s2 = (sp.sin(th)*sp.cos(th)*geom["Ytheta"] + 3*s2*Y) / s2
    TFshape = sp.simplify(TFtt - TFpp_over_s2)
    eq.append(("AngTF", sp.cancel(ang_tf/TFshape)))

    return eq


def coefficient_matrix(eq, symbols):
    rows = []
    for name, expr in eq:
        expr = sp.expand(expr)
        row = [sp.expand(sp.diff(expr, z)) for z in symbols]

        remainder = sp.expand(expr - sum(row[j]*symbols[j] for j in range(4)))
        remainder = sp.cancel(remainder)

        if remainder != 0:
            log(f"WARNING: nonhomogeneous remainder in {name}: {remainder}")

        rows.append(row)

    return sp.Matrix(rows)


def roots_from_minor(M, lam, row_ids, digits=30):
    sub = M[list(row_ids), :]
    det = sp.det(sub)

    num, den = sp.fraction(sp.together(det))

    try:
        poly = sp.Poly(sp.expand(num), lam)
    except Exception:
        return [], None

    coeffs = [complex(sp.N(c, digits)) for c in poly.all_coeffs()]

    while coeffs and abs(coeffs[0]) < 1e-13:
        coeffs.pop(0)

    if len(coeffs) <= 1:
        return [], poly

    rr = np.roots(np.array(coeffs, dtype=np.complex128))
    rr = [z for z in rr if finite_complex(z)]

    return rr, poly


def find_bulk_modes(M, lam, row_names):
    candidates = []

    for rows in itertools.combinations(range(M.rows), 4):
        log(f"  determinant rows {[row_names[i] for i in rows]}")

        try:
            roots, poly = roots_from_minor(M, lam, rows)
        except Exception as exc:
            log(f"    failed: {exc}")
            continue

        for z in roots:
            try:
                Mn = numeric_matrix(M, lam, z)
                smin, smax = smallest_sv(Mn)
                rel = smin/max(smax, 1e-300)

                candidates.append({
                    "lambda_re": z.real,
                    "lambda_im": z.imag,
                    "sv_min": smin,
                    "sv_max": smax,
                    "relative_sv": rel,
                    "minor_rows": ",".join(row_names[i] for i in rows),
                    "poly_degree": int(poly.degree()) if poly is not None else -1,
                })
            except Exception as exc:
                log(f"    validation failed for lambda={z}: {exc}")

    if not candidates:
        return pd.DataFrame()

    df = pd.DataFrame(candidates).sort_values("relative_sv")

    kept = []
    for _, r in df.iterrows():
        z = complex(r.lambda_re, r.lambda_im)

        duplicate = False
        for q in kept:
            zq = complex(q["lambda_re"], q["lambda_im"])
            if abs(z-zq) < 2e-6:
                duplicate = True
                break

        if not duplicate:
            kept.append(r.to_dict())

    return pd.DataFrame(kept).sort_values("relative_sv")


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--omega-re", type=float, default=2/math.sqrt(27))
    ap.add_argument("--omega-im", type=float, default=-0.005)
    ap.add_argument("--outdir", default="phase_ix_output")
    ap.add_argument(
        "--accept-sv",
        type=float,
        default=1e-7,
        help="relative smallest singular value used as a provisional good-mode threshold"
    )

    args = ap.parse_args()

    out = Path(args.outdir)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    omega_c = complex(args.omega_re, args.omega_im)
    omega = sp.Float(args.omega_re, 40) + sp.I*sp.Float(args.omega_im, 40)
    lam = sp.symbols("lambda")

    log("="*72)
    log("HBH Phase IX bulk Einstein--Vlasov")
    log(f"M omega = {omega_c}")
    log("="*72)

    t0 = time.time()

    log("[1/5] Building linearized Einstein tensor...")
    geom = build_geometry(omega, lam)
    log(f"      done in {time.time()-t0:.1f} s")

    log("[2/5] Building conservation-closed Vlasov source...")
    stress = build_stress(geom, omega, lam)

    log("[3/5] Projecting Einstein equations...")
    eq = projected_equations(geom, stress)

    row_names = [x[0] for x in eq]

    with open(out/"projected_equations.txt", "w", encoding="utf-8") as f:
        for name, expr in eq:
            f.write("\n" + "="*80 + "\n")
            f.write(name + "\n")
            f.write("="*80 + "\n")
            f.write(str(sp.factor(sp.together(expr))) + "\n")

    log("[4/5] Building 7x4 coefficient matrix...")
    M = coefficient_matrix(eq, geom["symbols"])

    with open(out/"coefficient_matrix.txt", "w", encoding="utf-8") as f:
        f.write("Columns: H0, H2, J, K\n")
        f.write("Rows: " + ", ".join(row_names) + "\n\n")
        for i, name in enumerate(row_names):
            f.write(f"[{name}]\n")
            for j, col in enumerate(["H0","H2","J","K"]):
                f.write(f"  {col}: {sp.factor(sp.together(M[i,j]))}\n")
            f.write("\n")

    log("[5/5] Searching radial lambda modes from all independent minors...")
    modes = find_bulk_modes(M, lam, row_names)

    if len(modes):
        modes["accepted"] = modes["relative_sv"] < args.accept_sv
        modes.to_csv(out/"bulk_lambda_modes.csv", index=False)

        log("\nBest bulk radial candidates:")
        print(modes[["lambda_re","lambda_im","relative_sv","minor_rows","accepted"]].head(20).to_string(index=False))
    else:
        log("No candidate roots extracted.")
        pd.DataFrame().to_csv(out/"bulk_lambda_modes.csv", index=False)

    with open(out/"kinetic_moments.txt", "w", encoding="utf-8") as f:
        for k,v in stress["moments"].items():
            f.write(f"{k} = {sp.factor(sp.together(v))}\n")

    metadata = {
        "omega_re": args.omega_re,
        "omega_im": args.omega_im,
        "omega_abs": abs(omega_c),
        "Omega0": 1/math.sqrt(27),
        "omega_kinetic": 2/math.sqrt(27),
        "polar_GW_threshold": math.sqrt(5/27),
        "elapsed_seconds": time.time()-t0,
        "rows": row_names,
        "columns": ["H0","H2","J","K"],
        "note": "Bulk ocean calculation only. Global HBH QNM requires matching to inner and outer vacuum Zerilli solutions."
    }

    (out/"metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    zip_path = Path(str(out) + ".zip")
    if zip_path.exists():
        zip_path.unlink()

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in out.rglob("*"):
            if p.is_file():
                zf.write(p, p.relative_to(out.parent))

    log("")
    log("="*72)
    log("DONE")
    log(f"Output folder: {out}")
    log(f"ZIP to send back: {zip_path}")
    log(f"Elapsed: {time.time()-t0:.1f} s")
    log("="*72)


if __name__ == "__main__":
    main()
