"""
Independent blind verification solver.
Method: real-axis multiple shooting (Riccati/log-derivative form) for both
vacuum regions, with the outer (far-field) boundary handled by starting an
8th-order analytic asymptotic series at a finite Rmax and integrating inward
-- validated stable in a window Rmax in [15,80] against the known vacuum
benchmark (see schwarzschild_benchmark.csv). NO exterior complex scaling and
NO Chebyshev/spectral collocation matrix is used anywhere in this file.

Matter bulk: independently-derived 7x4 Einstein-Vlasov system (bulk_blind.py)
with chi_J(omega), chi_X(omega) evaluated on the fly via the retarded
spectral continuation (spectral.py) from the supplied blind spectral-node
data. Matter transfer uses the two physical bulk eigenmodes propagated by
their exponential factors over the shell length L=9 ln(1/q), matched to the
Zerilli state via Z=K-iJ/(9w), ZX=K/27+26iJ/(243w).
"""
import math
import time
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from bulk_blind import build_geometry, build_stress, projected_equations, coefficient_matrix
from spectral import build_catalog, chi_JX

Q = 0.1
L_SHELL = 9.0 * math.log(1.0 / Q)
O_CONST = 1.0 / math.sqrt(27.0)


def Vz(R):
    f = 1.0 - 2.0 / R
    return 2 * f * (12 * R**3 + 12 * R**2 + 18 * R + 9) / (R**3 * (2 * R + 3) ** 2)


def horizon_y0(w, eps):
    u = eps
    b1 = 27 / (28 * (1 - 4j * w))
    b2 = -3j * (736 * w * w + 305j * w - 91) / (196 * (2 * w + 1j) * (4 * w + 1j) ** 2)
    return -1j * w + b1 * u + b2 * u * u


def inner_ratio(w, eps=1e-7, rtol=3e-12, atol=1e-14):
    R0 = 2.0 + eps
    y0 = horizon_y0(w, eps)
    sol = solve_ivp(lambda R, y: [(Vz(R) - w * w - y[0] * y[0]) / (1 - 2.0 / R)],
                    (R0, 3.0), [y0], method="DOP853", rtol=rtol, atol=atol, max_step=0.05)
    if not sol.success:
        raise RuntimeError("inner integration failed")
    return sol.y[0, -1]


def infinity_series_y(w, R, order=8):
    coeffs = {
        1: 0.0, 2: -3j / w, 3: 3 * (4j * w - 1) / w**2, 4: 3 * (-37j * w + 32) / (4 * w**2),
        5: 3 * (42j * w**3 - 61 * w**2 - 14j * w - 6) / (2 * w**4),
        6: 3 * (-651j * w**4 + 1432 * w**3 + 904j * w**2 + 296 * w + 144j) / (16 * w**5),
        7: 3 * (1197j * w**5 - 3633 * w**4 - 3952j * w**3 - 172 * w**2 - 1848j * w + 288) / (16 * w**6),
        8: 3 * (-8505j * w**6 + 32382 * w**5 + 53076j * w**4 - 18176 * w**3 + 37736j * w**2 - 24816 * w - 1728j) / (64 * w**7),
    }
    y = 1j * w
    for k in range(1, min(order, 8) + 1):
        y += coeffs[k] / R**k
    return y


def outer_ratio(w, Rmax=40.0, order=8, rtol=3e-12, atol=1e-14):
    y0 = infinity_series_y(w, Rmax, order=order)
    sol = solve_ivp(lambda R, y: [(Vz(R) - w * w - y[0] * y[0]) / (1 - 2.0 / R)],
                    (Rmax, 3.0), [y0], method="DOP853", rtol=rtol, atol=atol, max_step=0.5)
    if not sol.success:
        raise RuntimeError("outer integration failed")
    return sol.y[0, -1]


print("[setup] building independent 7x4 bulk matrix...", flush=True)
_t0 = time.time()
w_sym, lam_sym, chiJ_sym, chiX_sym = sp.symbols("omega lambda chi_J chi_X")
_geom = build_geometry(w_sym, lam_sym)
_stress = build_stress(_geom, w_sym, lam_sym, chiJ_sym, chiX_sym)
_eq = projected_equations(_geom, _stress)
_row_names = [n for n, _ in _eq]
_M = coefficient_matrix(_eq, _geom["symbols"])
Mfun = sp.lambdify((w_sym, lam_sym, chiJ_sym, chiX_sym), _M, "numpy")
print(f"[setup] done in {time.time()-_t0:.1f}s", flush=True)

print("[setup] loading spectral data...", flush=True)
_npz = np.load("blind_q010_spectral_nodes.npz")
_radial, _carrier = build_catalog(_npz)
print(f"[setup] {len(_radial)} radial branches, {len(_carrier)} carrier shells", flush=True)


def row_normalized(A):
    A = A.copy()
    for i in range(A.shape[0]):
        n = np.linalg.norm(A[i])
        if n > 0:
            A[i] /= n
    return A


def nullvector(Mn):
    A = row_normalized(Mn)
    U, s, Vh = np.linalg.svd(A)
    v = Vh.conj().T[:, -1]
    k = int(np.argmax(np.abs(v)))
    if abs(v[k]) > 0:
        v = v / v[k]
    return v, float(s[-1] / max(s[0], 1e-300))


def find_bulk_lambda(w, chiJ, chiX, seed):
    idx = [_row_names.index(n) for n in ["TX", "XX", "Xtheta", "AngTrace"]]

    def F(x):
        lam = complex(x[0], x[1])
        Mn = np.asarray(Mfun(w, lam, chiJ, chiX), dtype=np.complex128)
        A = Mn[idx, :].copy()
        for i in range(4):
            n = np.linalg.norm(A[i])
            if n > 0:
                A[i] /= n
        d = np.linalg.det(A)
        if not (np.isfinite(d.real) and np.isfinite(d.imag)):
            return [1e3, 1e3]
        return [d.real, d.imag]

    ans = least_squares(F, [seed.real, seed.imag], xtol=2e-13, ftol=2e-13, gtol=2e-13, max_nfev=400)
    lam = complex(ans.x[0], ans.x[1])
    Mn = np.asarray(Mfun(w, lam, chiJ, chiX), dtype=np.complex128)
    v, relsv = nullvector(Mn)
    return lam, v, relsv


def bulk_two_lambdas(w, seeds):
    chiJ, chiX = chi_JX(w, _radial, _carrier, O_CONST)
    lam1, v1, rel1 = find_bulk_lambda(w, chiJ, chiX, seeds[0])
    lam2, v2, rel2 = find_bulk_lambda(w, chiJ, chiX, seeds[1])
    return (lam1, v1, rel1), (lam2, v2, rel2), (chiJ, chiX)


def direct_column(w, v):
    H0, H2, J, K = v
    Z = K - 1j * J / (9 * w)
    ZX = K / 27 + 26j * J / (243 * w)
    return np.array([Z, ZX], dtype=np.complex128)


def matter_transfer(w, seeds):
    (lam1, v1, rel1), (lam2, v2, rel2), (chiJ, chiX) = bulk_two_lambdas(w, seeds)
    c1 = direct_column(w, v1)
    c2 = direct_column(w, v2)
    c1 = c1 / c1[0] if abs(c1[0]) > 1e-12 else c1 / c1[1]
    c2 = c2 / c2[0] if abs(c2[0]) > 1e-12 else c2 / c2[1]
    C = np.column_stack([c1, c2])
    P = np.diag(np.exp(np.array([lam1, lam2]) * L_SHELL))
    T = C @ P @ np.linalg.inv(C)
    return T, (lam1, lam2), (rel1, rel2), (chiJ, chiX)


def global_condition(w, seeds, Rmax=40.0):
    yin = inner_ratio(w)
    yout = outer_ratio(w, Rmax=Rmax)
    zin = np.array([1.0 + 0j, yin], dtype=np.complex128)
    T, lams, relsvs, chis = matter_transfer(w, seeds)
    zpred = T @ zin
    zoutvec = np.array([1.0 + 0j, yout], dtype=np.complex128)
    D = zpred[0] * zoutvec[1] - zpred[1] * zoutvec[0]
    scale = max(np.linalg.norm(zpred) * np.linalg.norm(zoutvec), 1e-300)
    return D / scale, lams, relsvs, chis
