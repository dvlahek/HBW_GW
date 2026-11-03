"""Physical-sheet upper-half-plane winding audit for CQG-117441.

The contour crosses Re(M omega)=0 but remains strictly above omega=0:
-0.05 <= Re(M omega) <= 1, 0.001 <= Im(M omega) <= 0.5.
For Im(omega)>0 spectral.retarded_transform uses the direct transform only,
so no second-sheet Plemelj correction enters this audit.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from independent_solver import global_condition

HERE = Path(__file__).resolve().parent
BASE_SEEDS = (-0.06 + 0.02j, -0.15 - 0.02j)
FALLBACKS = (
    (-0.04 + 0.05j, -0.18 - 0.05j),
    (-0.08 + 0.08j, -0.14 - 0.08j),
)
ANCHOR_RE = 0.377514814460


def contour(xmin, xmax, ymin, ymax, n):
    width, height = xmax - xmin, ymax - ymin
    per = 2 * (width + height)
    nb = max(2, round(n * width / per))
    nr = max(2, round(n * height / per))
    nt = max(2, round(n * width / per))
    nl = max(2, n - nb - nr - nt)
    pts = (
        [complex(x, ymin) for x in np.linspace(xmin, xmax, nb, endpoint=False)]
        + [complex(xmax, y) for y in np.linspace(ymin, ymax, nr, endpoint=False)]
        + [complex(x, ymax) for x in np.linspace(xmax, xmin, nt, endpoint=False)]
        + [complex(xmin, y) for y in np.linspace(ymax, ymin, nl, endpoint=False)]
    )
    i0 = min(range(len(pts)), key=lambda i: abs(pts[i].imag-ymin)+abs(pts[i].real-ANCHOR_RE))
    pts = pts[i0:] + pts[:i0]
    return pts + [pts[0]]


def solve_point(w, previous, rmax):
    """Use root continuation first; invoke fallback seeds only if the gate fails."""
    attempts = [previous]
    best = None
    errors = []
    for j, seeds in enumerate(attempts + list(FALLBACKS)):
        try:
            D, lams, relsvs, chis = global_condition(w, seeds, Rmax=rmax)
            sep = abs(lams[0] - lams[1])
            maxsv = max(float(relsvs[0]), float(relsvs[1]))
            valid = np.isfinite(D.real) and np.isfinite(D.imag) and sep > 1e-5
            score = maxsv + (0 if valid else 1e6)
            cand = (score, D, tuple(lams), tuple(relsvs), tuple(chis), sep)
            if best is None or cand[0] < best[0]:
                best = cand
            # Normal path: previous roots remain distinct and satisfy a generous
            # full-system residual gate, so do not spend extra seed searches.
            if j == 0 and valid and maxsv < 1e-5:
                break
        except Exception as exc:  # fallback is part of the audit robustness
            errors.append(f"{type(exc).__name__}: {exc}")
    if best is None or best[0] >= 1e6:
        raise RuntimeError(f"failed at w={w}: {' | '.join(errors)}")
    _, D, lams, relsvs, chis, sep = best
    return D, lams, relsvs, chis, sep


def run(label, n, bounds, rmax, outdir):
    pts = contour(*bounds, n)
    previous = BASE_SEEDS
    vals, rows = [], []
    for k, w in enumerate(pts):
        D, lams, relsvs, chis, sep = solve_point(w, previous, rmax)
        previous = lams
        vals.append(D)
        rows.append({
            "k": k, "w_re": w.real, "w_im": w.imag,
            "D_re": D.real, "D_im": D.imag, "abs_D": abs(D),
            "arg_D": float(np.angle(D)),
            "lambda1_re": lams[0].real, "lambda1_im": lams[0].imag,
            "lambda2_re": lams[1].real, "lambda2_im": lams[1].imag,
            "lambda_sep": sep,
            "rel_sv1": float(relsvs[0]), "rel_sv2": float(relsvs[1]),
            "chiJ_re": chis[0].real, "chiJ_im": chis[0].imag,
            "chiX_re": chis[1].real, "chiX_im": chis[1].imag,
        })
        if k % 25 == 0 or k == len(pts)-1:
            print(f"[{label}] {k}/{len(pts)-1} w={w.real:+.6f}{w.imag:+.6f}i |D|={abs(D):.6e} maxSV={max(relsvs):.3e} sep={sep:.3e}", flush=True)

    z = np.asarray(vals, complex)
    phase_steps = np.angle(z[1:] / z[:-1])
    unwrapped = np.unwrap(np.angle(z))
    winding = float(np.sum(phase_steps)/(2*np.pi))
    for row, ph in zip(rows, unwrapped):
        row["arg_D_unwrapped"] = float(ph)

    csv_path = outdir / f"upper_half_plane_{label}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader(); wr.writerows(rows)

    a = np.abs(z); i = int(np.argmin(a))
    summary = {
        "label": label,
        "physical_sheet": True,
        "bounds": {"re_min": bounds[0], "re_max": bounds[1], "im_min": bounds[2], "im_max": bounds[3]},
        "requested_contour_points": n,
        "evaluated_points_including_closure": len(pts),
        "rmax": rmax,
        "winding": winding,
        "winding_integer": int(round(winding)),
        "winding_unwrapped": float((unwrapped[-1]-unwrapped[0])/(2*np.pi)),
        "min_abs_D": float(a[i]),
        "min_abs_D_at": {"re": rows[i]["w_re"], "im": rows[i]["w_im"]},
        "max_abs_phase_step_rad": float(np.max(np.abs(phase_steps))),
        "phase_step_below_pi_over_2": bool(np.max(np.abs(phase_steps)) < np.pi/2),
        "closure_relative_D_mismatch": float(abs(z[-1]-z[0])/max(abs(z[0]),1e-300)),
        "max_full_bulk_relative_sv": float(max(max(r["rel_sv1"],r["rel_sv2"]) for r in rows)),
        "min_bulk_root_separation": float(min(r["lambda_sep"] for r in rows)),
        "csv": csv_path.name,
    }
    return summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-primary", type=int, default=384)
    p.add_argument("--n-dense", type=int, default=768)
    p.add_argument("--re-min", type=float, default=-0.05)
    p.add_argument("--re-max", type=float, default=1.0)
    p.add_argument("--im-min", type=float, default=0.001)
    p.add_argument("--im-max", type=float, default=0.5)
    p.add_argument("--rmax", type=float, default=40.0)
    p.add_argument("--outdir", type=Path, default=HERE/"upper_half_plane_audit_output")
    a = p.parse_args()
    if a.im_min <= 0 or not (a.re_min < 0 < a.re_max):
        raise ValueError("contour must cross Re(omega)=0 and stay strictly above Im(omega)=0")
    a.outdir.mkdir(parents=True, exist_ok=True)
    bounds = (a.re_min,a.re_max,a.im_min,a.im_max)
    ss = [run("primary",a.n_primary,bounds,a.rmax,a.outdir), run("dense",a.n_dense,bounds,a.rmax,a.outdir)]
    agree = ss[0]["winding_integer"] == ss[1]["winding_integer"]
    zero = all(s["winding_integer"] == 0 for s in ss)
    phase = all(s["phase_step_below_pi_over_2"] for s in ss)
    final = {
        "audit":"physical-sheet upper-half-plane growing-mode audit",
        "interpretation_scope":"finite-window linear spectral audit only",
        "contours":ss,
        "integer_windings_agree":agree,
        "both_windings_zero":zero,
        "phase_continuity_gate":phase,
        "verdict":"PASS_NO_GROWING_ROOTS_IN_TESTED_REGION" if agree and zero and phase else "REVIEW_REQUIRED",
    }
    (a.outdir/"upper_half_plane_audit_summary.json").write_text(json.dumps(final,indent=2),encoding="utf-8")
    lines=["PHYSICAL-SHEET UPPER-HALF-PLANE AUDIT","======================================",f"Region: {a.re_min} <= Re(Mw) <= {a.re_max}, {a.im_min} <= Im(Mw) <= {a.im_max}",""]
    for s in ss:
        lines += [f"{s['label']}:",f"  points: {s['requested_contour_points']}",f"  winding: {s['winding']:.12f} (integer {s['winding_integer']})",f"  min |D|: {s['min_abs_D']:.12e} at {s['min_abs_D_at']['re']:+.9f}{s['min_abs_D_at']['im']:+.9f}i",f"  max |phase step|: {s['max_abs_phase_step_rad']:.6f} rad",f"  closure relative D mismatch: {s['closure_relative_D_mismatch']:.3e}",f"  max full-bulk relative SV: {s['max_full_bulk_relative_sv']:.3e}",f"  min bulk-root separation: {s['min_bulk_root_separation']:.6e}",""]
    lines.append(f"VERDICT: {final['verdict']}")
    (a.outdir/"upper_half_plane_audit_summary.txt").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print("\n".join(lines),flush=True)
    if final["verdict"] != "PASS_NO_GROWING_ROOTS_IN_TESTED_REGION":
        raise SystemExit(2)

if __name__ == "__main__":
    main()
