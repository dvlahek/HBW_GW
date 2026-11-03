"""
Physical-sheet upper-half-plane argument-principle audit for the q=0.1
Einstein--Vlasov global determinant.

Purpose
-------
Test a rectangular region that crosses Re(M omega)=0 while staying strictly
above omega=0:

    -0.05 <= Re(M omega) <= 1.0
     0.001 <= Im(M omega) <= 0.5

For Im(omega)>0, spectral.retarded_transform() uses the direct Cauchy
transform only, so this is the physical sheet (no lower-half-plane/second-sheet
Plemelj correction).

The script evaluates two nested contour resolutions, follows the two physical
bulk roots continuously around each contour, writes the full contour CSVs and
summary JSON, and reports winding, min |D|, phase-step continuity, full-bulk
relative singular-value residuals, and bulk-root separation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

from independent_solver import global_condition


HERE = Path(__file__).resolve().parent
DEFAULT_SEEDS = (-0.06 + 0.02j, -0.15 - 0.02j)
ALT_SEED_PAIRS = (
    DEFAULT_SEEDS,
    (-0.04 + 0.05j, -0.18 - 0.05j),
    (-0.08 + 0.08j, -0.14 - 0.08j),
)
ANCHOR_RE = 0.377514814460


def rectangular_contour(
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
    n_total: int,
    anchor_re: float = ANCHOR_RE,
) -> list[complex]:
    """Return a CCW rectangular contour, rotated to start near anchor_re on bottom."""
    width = xmax - xmin
    height = ymax - ymin
    perimeter = 2.0 * (width + height)

    # Allocate points approximately uniformly in arclength, retaining at least
    # two points per edge. endpoint=False avoids corner duplication.
    nb = max(2, int(round(n_total * width / perimeter)))
    nr = max(2, int(round(n_total * height / perimeter)))
    nt = max(2, int(round(n_total * width / perimeter)))
    nl = max(2, n_total - nb - nr - nt)

    bottom = [complex(x, ymin) for x in np.linspace(xmin, xmax, nb, endpoint=False)]
    right = [complex(xmax, y) for y in np.linspace(ymin, ymax, nr, endpoint=False)]
    top = [complex(x, ymax) for x in np.linspace(xmax, xmin, nt, endpoint=False)]
    left = [complex(xmin, y) for y in np.linspace(ymax, ymin, nl, endpoint=False)]
    pts = bottom + right + top + left

    # Start on the lower edge close to the already-audited matter-mode real
    # frequency. This gives the bulk-root continuation a well-conditioned anchor
    # while leaving the mathematical contour unchanged.
    i0 = min(
        range(len(pts)),
        key=lambda i: abs(pts[i].imag - ymin) + abs(pts[i].real - anchor_re),
    )
    pts = pts[i0:] + pts[:i0]
    pts.append(pts[0])
    return pts


def _candidate_seed_pairs(previous: tuple[complex, complex] | None) -> Iterable[tuple[complex, complex]]:
    if previous is not None:
        yield previous
        # A small perturbation helps if one least-squares solve falls into the
        # same bulk root twice at an isolated contour point.
        yield (previous[0] + 0.005j, previous[1] - 0.005j)
    yield from ALT_SEED_PAIRS


def evaluate_point(
    w: complex,
    previous: tuple[complex, complex] | None,
    rmax: float,
) -> tuple[complex, tuple[complex, complex], tuple[float, float], tuple[complex, complex]]:
    """Evaluate global_condition robustly, preferring distinct low-residual bulk roots."""
    candidates = []
    errors = []
    seen = set()

    for seeds in _candidate_seed_pairs(previous):
        key = tuple((round(z.real, 10), round(z.imag, 10)) for z in seeds)
        if key in seen:
            continue
        seen.add(key)
        try:
            D, lams, relsvs, chis = global_condition(w, seeds, Rmax=rmax)
            if not (np.isfinite(D.real) and np.isfinite(D.imag)):
                raise FloatingPointError("non-finite determinant")
            sep = abs(lams[0] - lams[1])
            max_rel = max(float(relsvs[0]), float(relsvs[1]))
            # Penalize duplicate/near-duplicate bulk roots heavily. Otherwise
            # choose the candidate with the smallest full-system residual.
            penalty = 1e6 if sep < 1e-5 else 0.0
            score = penalty + max_rel
            candidates.append((score, -sep, D, tuple(lams), tuple(relsvs), tuple(chis)))
        except Exception as exc:  # noqa: BLE001 -- audit should try fallbacks
            errors.append(f"{type(exc).__name__}: {exc}")

    if not candidates:
        raise RuntimeError(f"all bulk-root attempts failed at w={w}: {' | '.join(errors)}")

    candidates.sort(key=lambda x: (x[0], x[1]))
    _, _, D, lams, relsvs, chis = candidates[0]
    if abs(lams[0] - lams[1]) < 1e-5:
        raise RuntimeError(f"bulk roots unresolved/duplicated at w={w}: {lams}")
    return D, lams, relsvs, chis


def run_contour(
    label: str,
    n_total: int,
    bounds: tuple[float, float, float, float],
    rmax: float,
    outdir: Path,
) -> dict:
    xmin, xmax, ymin, ymax = bounds
    points = rectangular_contour(xmin, xmax, ymin, ymax, n_total)

    rows = []
    vals: list[complex] = []
    previous = DEFAULT_SEEDS

    for k, w in enumerate(points):
        D, lams, relsvs, chis = evaluate_point(w, previous, rmax)
        previous = lams
        vals.append(D)
        rows.append(
            {
                "k": k,
                "w_re": w.real,
                "w_im": w.imag,
                "D_re": D.real,
                "D_im": D.imag,
                "abs_D": abs(D),
                "arg_D": float(np.angle(D)),
                "lambda1_re": lams[0].real,
                "lambda1_im": lams[0].imag,
                "lambda2_re": lams[1].real,
                "lambda2_im": lams[1].imag,
                "lambda_sep": abs(lams[0] - lams[1]),
                "rel_sv1": float(relsvs[0]),
                "rel_sv2": float(relsvs[1]),
                "chiJ_re": chis[0].real,
                "chiJ_im": chis[0].imag,
                "chiX_re": chis[1].real,
                "chiX_im": chis[1].imag,
            }
        )
        if k % 25 == 0 or k == len(points) - 1:
            print(
                f"[{label}] {k:4d}/{len(points)-1} "
                f"w={w.real:+.6f}{w.imag:+.6f}i |D|={abs(D):.6e} "
                f"maxSV={max(relsvs):.3e} sep={abs(lams[0]-lams[1]):.3e}",
                flush=True,
            )

    vals_arr = np.asarray(vals, dtype=np.complex128)
    # Incremental principal phase is the numerically robust argument-principle
    # sum. np.unwrap is also reported as an independent bookkeeping check.
    phase_steps = np.angle(vals_arr[1:] / vals_arr[:-1])
    winding_incremental = float(np.sum(phase_steps) / (2.0 * np.pi))
    unwrapped = np.unwrap(np.angle(vals_arr))
    winding_unwrapped = float((unwrapped[-1] - unwrapped[0]) / (2.0 * np.pi))

    for r, phu in zip(rows, unwrapped):
        r["arg_D_unwrapped"] = float(phu)

    csv_path = outdir / f"upper_half_plane_{label}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    absD = np.abs(vals_arr)
    min_idx = int(np.argmin(absD))
    max_rel = max(max(r["rel_sv1"], r["rel_sv2"]) for r in rows)
    min_sep = min(r["lambda_sep"] for r in rows)
    closure_rel = float(abs(vals_arr[-1] - vals_arr[0]) / max(abs(vals_arr[0]), 1e-300))
    max_phase_step = float(np.max(np.abs(phase_steps)))

    summary = {
        "label": label,
        "physical_sheet": True,
        "bounds": {
            "re_min": xmin,
            "re_max": xmax,
            "im_min": ymin,
            "im_max": ymax,
        },
        "requested_contour_points": n_total,
        "evaluated_points_including_closure": len(points),
        "rmax": rmax,
        "winding_incremental": winding_incremental,
        "winding_unwrapped": winding_unwrapped,
        "winding_integer": int(round(winding_incremental)),
        "min_abs_D": float(absD[min_idx]),
        "min_abs_D_at": {
            "re": rows[min_idx]["w_re"],
            "im": rows[min_idx]["w_im"],
        },
        "max_abs_phase_step_rad": max_phase_step,
        "phase_step_below_pi_over_2": bool(max_phase_step < np.pi / 2),
        "closure_relative_D_mismatch": closure_rel,
        "max_full_bulk_relative_sv": float(max_rel),
        "min_bulk_root_separation": float(min_sep),
        "csv": csv_path.name,
    }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-primary", type=int, default=384)
    parser.add_argument("--n-dense", type=int, default=768)
    parser.add_argument("--re-min", type=float, default=-0.05)
    parser.add_argument("--re-max", type=float, default=1.0)
    parser.add_argument("--im-min", type=float, default=0.001)
    parser.add_argument("--im-max", type=float, default=0.5)
    parser.add_argument("--rmax", type=float, default=40.0)
    parser.add_argument("--outdir", type=Path, default=HERE / "upper_half_plane_audit_output")
    args = parser.parse_args()

    if args.im_min <= 0:
        raise ValueError("im-min must be strictly positive so the contour never crosses omega=0")
    if not (args.re_min < 0 < args.re_max):
        raise ValueError("the audit rectangle must cross Re(omega)=0")

    args.outdir.mkdir(parents=True, exist_ok=True)
    bounds = (args.re_min, args.re_max, args.im_min, args.im_max)

    summaries = []
    for label, n in (("primary", args.n_primary), ("dense", args.n_dense)):
        summaries.append(run_contour(label, n, bounds, args.rmax, args.outdir))

    agree = summaries[0]["winding_integer"] == summaries[1]["winding_integer"]
    both_zero = all(s["winding_integer"] == 0 for s in summaries)
    phase_ok = all(s["phase_step_below_pi_over_2"] for s in summaries)

    final = {
        "audit": "physical-sheet upper-half-plane growing-mode audit",
        "interpretation_scope": "finite-window linear spectral audit only",
        "contours": summaries,
        "integer_windings_agree": agree,
        "both_windings_zero": both_zero,
        "phase_continuity_gate": phase_ok,
        "verdict": "PASS_NO_GROWING_ROOTS_IN_TESTED_REGION" if (agree and both_zero and phase_ok) else "REVIEW_REQUIRED",
    }

    summary_path = args.outdir / "upper_half_plane_audit_summary.json"
    summary_path.write_text(json.dumps(final, indent=2), encoding="utf-8")

    txt = [
        "PHYSICAL-SHEET UPPER-HALF-PLANE AUDIT",
        "======================================",
        f"Region: {args.re_min} <= Re(Mw) <= {args.re_max}, {args.im_min} <= Im(Mw) <= {args.im_max}",
        "Physical sheet: yes (Im w > 0; direct spectral transform only)",
        "",
    ]
    for s in summaries:
        txt.extend(
            [
                f"{s['label']}:",
                f"  requested points: {s['requested_contour_points']}",
                f"  winding: {s['winding_incremental']:.12f} (integer {s['winding_integer']})",
                f"  min |D|: {s['min_abs_D']:.12e} at {s['min_abs_D_at']['re']:+.9f}{s['min_abs_D_at']['im']:+.9f}i",
                f"  max |phase step|: {s['max_abs_phase_step_rad']:.6f} rad",
                f"  closure relative D mismatch: {s['closure_relative_D_mismatch']:.3e}",
                f"  max full-bulk relative SV: {s['max_full_bulk_relative_sv']:.3e}",
                f"  min bulk-root separation: {s['min_bulk_root_separation']:.6e}",
                "",
            ]
        )
    txt.append(f"VERDICT: {final['verdict']}")
    (args.outdir / "upper_half_plane_audit_summary.txt").write_text("\n".join(txt) + "\n", encoding="utf-8")

    print("\n".join(txt), flush=True)
    if final["verdict"] != "PASS_NO_GROWING_ROOTS_IN_TESTED_REGION":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
