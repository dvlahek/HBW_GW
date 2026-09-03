#!/usr/bin/env python3
"""Retarded spectral continuation used by the smooth Einstein--Vlasov QNM solver."""

import io
import math
import zipfile
from pathlib import Path

import numpy as np


def load_c7(path):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)

    with zipfile.ZipFile(p, "r") as zf:
        matches = [
            n for n in zf.namelist()
            if n.endswith("spectral_nodes_N80.npz")
        ]
        if len(matches) != 1:
            raise RuntimeError("Could not uniquely locate spectral_nodes_N80.npz.")
        raw = zf.read(matches[0])

    with np.load(io.BytesIO(raw)) as npz:
        return {k: np.array(npz[k]) for k in npz.files}


def trap_weights(x):
    x = np.asarray(x, dtype=float)
    if len(x) < 2:
        raise ValueError("Need at least two spectral nodes.")

    w = np.empty_like(x)
    w[0] = 0.5*(x[1]-x[0])
    w[-1] = 0.5*(x[-1]-x[-2])
    if len(x) > 2:
        w[1:-1] = 0.5*(x[2:]-x[:-2])
    return w


def prepare_branch(nu, q):
    idx = np.argsort(np.asarray(nu))
    x = np.asarray(nu, dtype=float)[idx]
    q = np.asarray(q, dtype=np.complex128)[idx]

    if not np.all(np.diff(x) > 0):
        raise RuntimeError("Spectral branch is not strictly monotone.")

    rho = q/trap_weights(x)
    a = x[:-1]
    b = x[1:]
    ra = rho[:-1]
    rb = rho[1:]
    m = (rb-ra)/(b-a)
    c = ra-m*a

    return {"x": x, "rho": rho, "a": a, "b": b, "m": m, "c": c}


def direct_branch(z, br):
    z = complex(z)
    a = br["a"]
    b = br["b"]
    m = br["m"]
    c = br["c"]
    L = np.log(z-a) - np.log(z-b)
    return np.sum((m*z+c)*L - m*(b-a))


def local_rho(z, br):
    z = complex(z)
    x = br["x"]
    xr = z.real

    if xr < x[0] or xr > x[-1]:
        return 0.0+0.0j

    j = np.searchsorted(x, xr, side="right") - 1
    j = max(0, min(int(j), len(x)-2))
    return br["m"][j]*z + br["c"][j]


def continued_branch(z, br):
    """Retarded second-sheet continuation through a vertical path."""
    z = complex(z)
    val = direct_branch(z, br)

    if z.imag < 0 and br["x"][0] <= z.real <= br["x"][-1]:
        val -= 2j*math.pi*local_rho(z, br)

    return val


def build_catalog(spc):
    nuJ = spc["nu_J"]
    qJJ = spc["q_JJ"]
    qJX = spc["q_JX"]
    nr = spc["n_r"]
    nphi = spc["n_phi"]
    shell = spc["shell_J"]

    radial = []
    for s in np.unique(shell):
        for a in [-3, -2, -1, 1, 2, 3]:
            for b in [-2, 0, 2]:
                mask = (shell == s) & (nr == a) & (nphi == b)
                if not np.any(mask):
                    continue
                radial.append({
                    "shell": int(s),
                    "n_r": int(a),
                    "n_phi": int(b),
                    "JJ": prepare_branch(nuJ[mask], qJJ[mask]),
                    "JX": prepare_branch(nuJ[mask], qJX[mask]),
                })

    nuH = spc["nu_H"]
    qHH = spc["q_HH"]
    shellH = spc["shell_H"]

    carrier = []
    for s in np.unique(shellH):
        mask = shellH == s
        carrier.append({
            "shell": int(s),
            "HH": prepare_branch(nuH[mask], qHH[mask]),
        })

    return radial, carrier


def transform_radial(z, radial, which, *, second_sheet):
    total = 0.0+0.0j
    for item in radial:
        br = item[which]
        total += continued_branch(z, br) if second_sheet else direct_branch(z, br)
    return total
