#!/usr/bin/env python3
"""
HBH q=0.1 microscopic spectral-profile robustness test.

Purpose
-------
Test whether the matter-led global QNM survives smooth, finite deformations of
how the collisionless response is distributed across the 80 publication
shells. This is a response-profile robustness test. It does not rebuild a new
self-consistent stationary Einstein--Vlasov background and must not be
presented as a second equilibrium family.

Run
---
  py profile_shape_robustness.py

Outputs
-------
  profile_shape_roots.csv
  profile_shape_stability.csv
  profile_shape_winding.csv
  profile_shape_summary.txt
  HBH_PROFILE_SHAPE_ROBUSTNESS_RESULTS.zip
"""
from pathlib import Path
import csv, math, zipfile
import numpy as np
from scipy.optimize import least_squares

import spectral
import independent_solver as solver

BASE_ROOT = 0.377514814460 - 0.020871654486j
LAMBDA_SEEDS = (-0.094 + 0.073j, -0.128 - 0.073j)
RE_WINDOW = (0.33, 0.43)
IM_WINDOW = (-0.065, -0.003)
RMAX_VALUES = (25.0, 40.0, 60.0)

PROFILE_CASES = [
    ("baseline", "baseline", 0.0),
    ("tilt_plus", "tilt", +0.35),
    ("tilt_minus", "tilt", -0.35),
    ("center_plus", "center", +0.45),
    ("center_minus", "center", -0.45),
    ("two_lobe_plus", "two_lobe", +0.35),
    ("two_lobe_minus", "two_lobe", -0.35),
]


def shell_weights(shell_ids, kind, amp):
    s = np.asarray(shell_ids, dtype=float)
    smin, smax = float(np.min(s)), float(np.max(s))
    x = 2.0 * (s - smin) / max(smax - smin, 1.0) - 1.0
    if kind == "baseline":
        raw = np.ones_like(x)
    elif kind == "tilt":
        raw = np.exp(amp * x)
    elif kind == "center":
        raw = np.exp(amp * (1.0 - x*x))
    elif kind == "two_lobe":
        raw = np.exp(amp * np.cos(2.0 * np.pi * x))
    else:
        raise ValueError(kind)
    uniq = np.unique(s)
    means = [float(np.mean(raw[s == u])) for u in uniq]
    raw = raw / np.mean(means)
    return raw


def deformed_catalog(npz, kind, amp):
    data = {k: np.array(npz[k], copy=True) for k in npz.files}
    wJ = shell_weights(data["shell_J"], kind, amp)
    wH = shell_weights(data["shell_H"], kind, amp)
    data["q_JJ"] *= wJ
    data["q_JX"] *= wJ
    data["q_HH"] *= wH
    return spectral.build_catalog(data), (float(wJ.min()), float(wJ.max()))


def set_catalog(radial, carrier):
    solver._radial = radial
    solver._carrier = carrier


def refine_root(seed, Rmax=40.0):
    def F(x):
        w = complex(x[0], x[1])
        try:
            D, lams, relsvs, chis = solver.global_condition(w, LAMBDA_SEEDS, Rmax=Rmax)
            if not (np.isfinite(D.real) and np.isfinite(D.imag)):
                return [1e3, 1e3]
            return [D.real, D.imag]
        except Exception:
            return [1e3, 1e3]
    ans = least_squares(F, [seed.real, seed.imag], xtol=3e-12, ftol=3e-12,
                        gtol=3e-12, max_nfev=180)
    w = complex(ans.x[0], ans.x[1])
    D, lams, relsvs, chis = solver.global_condition(w, LAMBDA_SEEDS, Rmax=Rmax)
    return w, D, lams, relsvs, chis, bool(ans.success)


def winding_number(root, half_re=0.0035, half_im=0.0035, nside=8):
    pts = []
    for t in np.linspace(0, 1, nside, endpoint=False):
        pts.append(root + complex(-half_re + 2*half_re*t, -half_im))
    for t in np.linspace(0, 1, nside, endpoint=False):
        pts.append(root + complex(+half_re, -half_im + 2*half_im*t))
    for t in np.linspace(0, 1, nside, endpoint=False):
        pts.append(root + complex(+half_re - 2*half_re*t, +half_im))
    for t in np.linspace(0, 1, nside, endpoint=False):
        pts.append(root + complex(-half_re, +half_im - 2*half_im*t))
    vals = []
    for w in pts:
        D, *_ = solver.global_condition(w, LAMBDA_SEEDS, Rmax=40.0)
        vals.append(D)
    vals = np.asarray(vals, dtype=np.complex128)
    phases = np.unwrap(np.angle(np.r_[vals, vals[0]]))
    winding = float((phases[-1] - phases[0])/(2*np.pi))
    return winding, pts, vals


def main():
    here = Path(__file__).resolve().parent
    npz_path = here / "blind_q010_spectral_nodes.npz"
    if not npz_path.exists():
        raise FileNotFoundError(f"Missing {npz_path.name} next to this script")
    npz = np.load(npz_path)

    roots_rows = []
    stab_rows = []
    winding_rows = []
    all_pass = True

    print("\nHBH microscopic spectral-profile robustness test")
    print("This test deforms shell-wise response weights and recomputes the full global pole.\n")

    for name, kind, amp in PROFILE_CASES:
        print(f"[{name}] building profile ...", flush=True)
        (radial, carrier), (wmin, wmax) = deformed_catalog(npz, kind, amp)
        set_catalog(radial, carrier)

        w40, D40, lams, relsvs, chis, ok = refine_root(BASE_ROOT, 40.0)
        root_shift = abs(w40 - BASE_ROOT)
        in_window = RE_WINDOW[0] < w40.real < RE_WINDOW[1] and IM_WINDOW[0] < w40.imag < IM_WINDOW[1]
        svmax = max(relsvs)

        rroots = []
        for Rmax in RMAX_VALUES:
            wr, Dr, lr, svr, chir, okr = refine_root(w40, Rmax)
            rroots.append(wr)
            stab_rows.append({
                "profile": name, "kind": kind, "amp": amp, "Rmax": Rmax,
                "root_re": wr.real, "root_im": wr.imag, "abs_D": abs(Dr),
                "relsv_max": max(svr), "success": okr,
            })
        rshift = max(abs(r - rroots[1]) for r in rroots)

        try:
            winding, pts, vals = winding_number(w40)
        except Exception as exc:
            winding = float("nan")
            pts, vals = [], []
            print(f"  winding failed: {exc!r}")

        for j, (p, v) in enumerate(zip(pts, vals)):
            winding_rows.append({
                "profile": name, "index": j, "omega_re": p.real, "omega_im": p.imag,
                "D_re": v.real, "D_im": v.imag, "abs_D": abs(v),
            })

        case_pass = (
            ok and in_window and abs(D40) < 5e-7 and svmax < 1e-8 and
            rshift < 5e-5 and np.isfinite(winding) and abs(winding - 1.0) < 0.08
        )
        all_pass &= case_pass
        roots_rows.append({
            "profile": name, "kind": kind, "amp": amp,
            "weight_min": wmin, "weight_max": wmax,
            "root_re": w40.real, "root_im": w40.imag,
            "delta_from_baseline_abs": root_shift,
            "delta_from_baseline_rel": root_shift/max(abs(BASE_ROOT),1e-300),
            "abs_D": abs(D40), "relsv_max": svmax,
            "lambda1_re": lams[0].real, "lambda1_im": lams[0].imag,
            "lambda2_re": lams[1].real, "lambda2_im": lams[1].imag,
            "chiJ_re": chis[0].real, "chiJ_im": chis[0].imag,
            "chiX_re": chis[1].real, "chiX_im": chis[1].imag,
            "Rmax_max_shift": rshift, "winding": winding,
            "pass": case_pass,
        })
        print(f"  omega = {w40.real:.12f} {w40.imag:+.12f} i")
        print(f"  weight range = [{wmin:.3f}, {wmax:.3f}]")
        print(f"  |D|={abs(D40):.3e}, relSV={svmax:.3e}, Rmax shift={rshift:.3e}, winding={winding:.6f}")
        print(f"  {'PASS' if case_pass else 'FAIL'}\n")

    def write_csv(path, rows):
        if not rows:
            return
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)

    write_csv(here/"profile_shape_roots.csv", roots_rows)
    write_csv(here/"profile_shape_stability.csv", stab_rows)
    write_csv(here/"profile_shape_winding.csv", winding_rows)

    verdict = "HBH_PROFILE_SHAPE_ROBUSTNESS_PASS" if all_pass else "HBH_PROFILE_SHAPE_ROBUSTNESS_REVIEW"
    lines = [
        verdict,
        "",
        "Interpretation:",
        "This is a controlled microscopic response-profile robustness test.",
        "It does not by itself constitute a second self-consistent stationary Einstein--Vlasov family.",
        "PASS means the same simple global matter-led pole survives all smooth shell-profile deformations,",
        "with unit winding, full-bulk residual control, and stable outer-vacuum matching.",
        "",
    ]
    for r in roots_rows:
        lines.append(f"{r['profile']}: omega={r['root_re']:.12f}{r['root_im']:+.12f}i, "
                     f"rel_shift={r['delta_from_baseline_rel']:.4e}, winding={r['winding']:.6f}, pass={r['pass']}")
    (here/"profile_shape_summary.txt").write_text("\n".join(lines)+"\n", encoding="utf-8")

    zout = here/"HBH_PROFILE_SHAPE_ROBUSTNESS_RESULTS.zip"
    with zipfile.ZipFile(zout, "w", zipfile.ZIP_DEFLATED) as z:
        for fn in ["profile_shape_roots.csv","profile_shape_stability.csv","profile_shape_winding.csv","profile_shape_summary.txt"]:
            z.write(here/fn, fn)
    print(verdict)
    print(f"Results: {zout}")

if __name__ == "__main__":
    main()
