"""
Argument-principle winding-number audit around the root obtained by the blind
search. This script does not rerun the search.
"""
import csv
import numpy as np
from independent_solver import global_condition

W_ROOT = 0.3775148114965411 - 0.020871653994053295j
LAMBDA_SEED_PAIR = (-0.06 + 0.02j, -0.15 - 0.02j)
N_POINTS = 24
RADIUS = 3e-4


def run_winding_audit(w_center=W_ROOT, lam_pair=LAMBDA_SEED_PAIR, n=N_POINTS, radius=RADIUS):
    rows = []
    vals = []
    for k in range(n + 1):
        th = 2 * np.pi * k / n
        w = w_center + radius * np.exp(1j * th)
        D, lams, relsvs, chis = global_condition(w, lam_pair, Rmax=40.0)
        vals.append(D)
        rows.append({
            "k": k,
            "theta": th,
            "w_re": w.real,
            "w_im": w.imag,
            "D_re": D.real,
            "D_im": D.imag,
            "abs_D": abs(D),
        })

    phases = np.unwrap(np.angle(vals))
    winding = (phases[-1] - phases[0]) / (2 * np.pi)
    return winding, rows


if __name__ == "__main__":
    winding, rows = run_winding_audit()
    print(f"Winding number around candidate root {W_ROOT}: {winding:.6f}")

    with open("winding_audit.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["k", "theta", "w_re", "w_im", "D_re", "D_im", "abs_D"])
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    with open("winding_audit_result.txt", "w") as f:
        f.write(f"Center (converged root): {W_ROOT.real:.10f} {W_ROOT.imag:+.10f}i\n")
        f.write(f"Contour radius: {RADIUS}\n")
        f.write(f"Number of contour points: {N_POINTS}\n")
        f.write(f"Winding number: {winding:.6f}\n")
        f.write("Interpretation: winding=1 indicates a simple isolated pole.\n")
