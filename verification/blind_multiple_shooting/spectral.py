"""
Retarded spectral continuation for chi_J(z), chi_X(z) used in the blind
verification.
"""
import numpy as np
import math


def trap_weights(nu):
    nu = np.asarray(nu, dtype=float)
    w = np.empty_like(nu)
    w[0] = (nu[1] - nu[0]) / 2
    w[-1] = (nu[-1] - nu[-2]) / 2
    if len(nu) > 2:
        w[1:-1] = (nu[2:] - nu[:-2]) / 2
    return w


def make_branch(nu, q):
    idx = np.argsort(nu)
    x = np.asarray(nu, dtype=float)[idx]
    qq = np.asarray(q, dtype=np.complex128)[idx]
    w = trap_weights(x)
    rho = qq / w
    a, b = x[:-1], x[1:]
    ra, rb = rho[:-1], rho[1:]
    m = (rb - ra) / (b - a)
    c = ra - m * a
    return {"x": x, "a": a, "b": b, "m": m, "c": c}


def direct_transform(z, br):
    z = complex(z)
    a, b, m, c = br["a"], br["b"], br["m"], br["c"]
    L = np.log(z - a) - np.log(z - b)
    return np.sum((m * z + c) * L - m * (b - a))


def local_rho(z, br):
    x = br["x"]
    xr = z.real
    if xr < x[0] or xr > x[-1]:
        return 0.0 + 0.0j
    j = np.searchsorted(x, xr, side="right") - 1
    j = max(0, min(j, len(x) - 2))
    return br["m"][j] * z + br["c"][j]


def retarded_transform(z, br):
    val = direct_transform(z, br)
    if z.imag < 0 and br["x"][0] <= z.real <= br["x"][-1]:
        val -= 2j * math.pi * local_rho(z, br)
    return val


def build_catalog(npz):
    nuJ, qJJ, qJX = npz["nu_J"], npz["q_JJ"], npz["q_JX"]
    nr, nphi, shellJ = npz["n_r"], npz["n_phi"], npz["shell_J"]

    radial = []
    for s in np.unique(shellJ):
        for a in [-3, -2, -1, 1, 2, 3]:
            for b in [-2, 0, 2]:
                mask = (shellJ == s) & (nr == a) & (nphi == b)
                if not np.any(mask) or mask.sum() < 2:
                    continue
                radial.append({
                    "JJ": make_branch(nuJ[mask], qJJ[mask]),
                    "JX": make_branch(nuJ[mask], qJX[mask]),
                })

    nuH, qHH, shellH = npz["nu_H"], npz["q_HH"], npz["shell_H"]
    carrier = []
    for s in np.unique(shellH):
        mask = shellH == s
        if mask.sum() < 2:
            continue
        carrier.append(make_branch(nuH[mask], qHH[mask]))

    return radial, carrier


def P_H(z, carrier):
    return sum(direct_transform(z, br) for br in carrier)


def P_JJ_II(z, radial):
    return sum(retarded_transform(z, item["JJ"]) for item in radial)


def P_JX_II(z, radial):
    return sum(retarded_transform(z, item["JX"]) for item in radial)


def chi_JX(z, radial, carrier, O):
    z = complex(z)
    dSplus = (9.0 * O / 8.0) / (2.0 * O - z)
    ph = P_H(z, carrier)
    kappa = 0.5 * dSplus / ph
    chiJ = kappa * P_JJ_II(z, radial)
    chiX = kappa * P_JX_II(z, radial)
    return chiJ, chiX
