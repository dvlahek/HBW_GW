"""
Blind root search over
0.34 < Re(Momega) < 0.42 and -0.05 < Im(Momega) < -0.005.

No target matter-mode frequency is used as a seed. The search uses a grid of
frequency seeds and two bulk-root seed pairs. Each converged candidate is
checked at three outer-vacuum radii.
"""
import csv
import time
import numpy as np
from scipy.optimize import least_squares

from independent_solver import global_condition

RE_LO, RE_HI = 0.34, 0.42
IM_LO, IM_HI = -0.05, -0.005

seed_grid = []
for re in np.linspace(RE_LO + 0.01, RE_HI - 0.01, 4):
    for im in np.linspace(IM_LO + 0.003, IM_HI - 0.001, 3):
        seed_grid.append(complex(re, im))

LAMBDA_SEED_PAIRS = [
    (-0.06 + 0.02j, -0.15 - 0.02j),
    (-0.04 + 0.05j, -0.18 - 0.05j),
]

RESOLUTIONS = [25.0, 40.0, 60.0]

all_seed_results = []
converged_roots = []

t_start = time.time()
for wseed in seed_grid:
    for lam_pair in LAMBDA_SEED_PAIRS:
        try:
            def real2(x):
                w = complex(x[0], x[1])
                try:
                    D, lams, relsvs, chis = global_condition(w, lam_pair, Rmax=40.0)
                    if not (np.isfinite(D.real) and np.isfinite(D.imag)):
                        return [1e3, 1e3]
                    return [D.real, D.imag]
                except Exception:
                    return [1e3, 1e3]

            ans = least_squares(
                real2,
                [wseed.real, wseed.imag],
                bounds=([RE_LO - 0.05, IM_LO - 0.05], [RE_HI + 0.05, IM_HI + 0.05]),
                xtol=2e-12,
                ftol=2e-12,
                gtol=2e-12,
                max_nfev=150,
            )
            w_root = complex(ans.x[0], ans.x[1])

            res_vals = []
            for Rmax in RESOLUTIONS:
                D, lams, relsvs, chis = global_condition(w_root, lam_pair, Rmax=Rmax)
                res_vals.append((Rmax, D, lams, relsvs))

            Dfinal = res_vals[-1][1]
            success = ans.success and abs(Dfinal) < 1e-4

            all_seed_results.append({
                "seed_re": wseed.real,
                "seed_im": wseed.imag,
                "lambda_seed_pair": str(lam_pair),
                "root_re": w_root.real,
                "root_im": w_root.imag,
                "abs_D_final": abs(Dfinal),
                "success": success,
                "in_window": (RE_LO < w_root.real < RE_HI) and (IM_LO < w_root.imag < IM_HI),
            })

            if success and (RE_LO < w_root.real < RE_HI) and (IM_LO < w_root.imag < IM_HI):
                converged_roots.append({
                    "w": w_root,
                    "seed": wseed,
                    "lambda_seed_pair": lam_pair,
                    "res_vals": res_vals,
                })
        except Exception as exc:
            all_seed_results.append({
                "seed_re": wseed.real,
                "seed_im": wseed.imag,
                "lambda_seed_pair": str(lam_pair),
                "root_re": np.nan,
                "root_im": np.nan,
                "abs_D_final": np.nan,
                "success": False,
                "in_window": False,
                "error": repr(exc),
            })

print(f"Total time: {time.time()-t_start:.1f}s")
print(f"Converged in window: {len(converged_roots)}")

deduped = []
for r in converged_roots:
    if not any(abs(r["w"] - d["w"]) < 1e-4 for d in deduped):
        deduped.append(r)

with open("all_frequency_seeds.csv", "w", newline="") as f:
    fields = ["seed_re", "seed_im", "lambda_seed_pair", "root_re", "root_im", "abs_D_final", "success", "in_window"]
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for r in all_seed_results:
        writer.writerow({k: r.get(k) for k in fields})

with open("all_converged_roots.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["root_re", "root_im", "seed_re", "seed_im", "lambda_seed_pair"])
    for d in deduped:
        writer.writerow([d["w"].real, d["w"].imag, d["seed"].real, d["seed"].imag, str(d["lambda_seed_pair"])])

with open("resolution_convergence.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["root_re", "root_im", "Rmax", "D_re", "D_im", "abs_D", "lambda1_re", "lambda1_im", "lambda2_re", "lambda2_im", "relsv1", "relsv2"])
    for d in deduped:
        for Rmax, D, lams, relsvs in d["res_vals"]:
            writer.writerow([
                d["w"].real, d["w"].imag, Rmax, D.real, D.imag, abs(D),
                lams[0].real, lams[0].imag, lams[1].real, lams[1].imag,
                relsvs[0], relsvs[1],
            ])
