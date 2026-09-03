#!/usr/bin/env python3
"""Final smooth Einstein--Vlasov global QNM machinery.

This publication-facing module keeps the numerical core used by the verified
C9-v4 calculation: the retarded smooth constitutive response, the
conservation-complete 7x4 bulk operator, the confluent ocean transfer, and the
global Zerilli determinant.  Development homotopy/checkpoint code is omitted.
"""

import cmath
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import least_squares

import hbh_phase_xxic8_landau_continuation as c8

Q = 0.1
W_SHARP = 0.370658573899 - 0.00590424942625j
C8_LAM_TARGETS = [
    -0.0363374333142 + 0.00491398482649j,
    -0.185884788908 - 0.00491398482649j,
]
ACTIVE_LAMBDA_SEED_PAIR = None


def load_xivc():
    here = Path(__file__).resolve().parent
    p = here / "hbh_phase_xivc_global_polar_EV_ecs.py"
    if not p.exists():
        raise FileNotFoundError("Put hbh_phase_xivc_global_polar_EV_ecs.py in the same folder.")

    spec = importlib.util.spec_from_file_location("xivc_c9", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["xivc_c9"] = mod
    spec.loader.exec_module(mod)
    return mod


def carrier_transform(z, carrier, *, second_sheet):
    total = 0.0+0.0j
    for item in carrier:
        br = item["HH"]
        total += c8.continued_branch(z, br) if second_sheet else c8.direct_branch(z, br)
    return total


def smooth_chi(z, radial, carrier):
    z = complex(z)
    second = z.imag < 0

    PH = carrier_transform(z, carrier, second_sheet=second)
    PJJ = c8.transform_radial(z, radial, "JJ", second_sheet=second)
    PJX = c8.transform_radial(z, radial, "JX", second_sheet=second)

    O = 1.0/math.sqrt(27.0)
    dSplus = (9.0*O/8.0)/(2.0*O-z)
    kappa = 0.5*dSplus/PH

    return {
        "PH": PH,
        "PJJ": PJJ,
        "PJX": PJX,
        "kappa": kappa,
        "chiJ": kappa*PJJ,
        "chiX": kappa*PJX,
    }


def build_generic_bulk_matrix(phase_ix):
    print("[setup] Building generic conservation-complete 7x4 matrix...", flush=True)

    w, lam = sp.symbols("omega lambda")
    chiX, chiJ, tau = sp.symbols("chi_X chi_J tau")

    geom = phase_ix.build_geometry(w, lam)
    base = phase_ix.build_stress(geom, w, lam)

    H0,H2,J,K = geom["symbols"]
    S = H0+K
    SX = lam*S

    kin = phase_ix.kinetic_coefficients(w)
    fSX = kin["f_SX"]

    df = tau*((chiX-fSX)*SX + chiJ*J)

    Oo = sp.sqrt(sp.Rational(1,27))
    dp = -sp.I*w*df/(6*Oo)
    de = (lam+sp.Rational(2,9))*df/(6*Oo)
    ds = -sp.Rational(9,2)*sp.I*w*(lam+sp.Rational(2,9))*df

    theta = geom["theta"]
    Y = geom["Y"]
    Yt = geom["Ytheta"]
    Ytt = geom["Ythetatheta"]
    s2 = sp.sin(theta)**2
    ct = sp.cos(theta)
    pref = sp.Rational(1,36)/sp.pi

    TFtt_up = Ytt + 3*Y
    TFpp_up = (sp.sin(theta)*ct*Yt + 3*s2*Y)/s2**2

    dTup = sp.zeros(4)
    dTup[0,0] = 0
    dTup[0,1] = dTup[1,0] = pref*df*Y
    dTup[1,1] = 0
    dTup[0,2] = dTup[2,0] = pref*Oo*de*Yt
    dTup[1,2] = dTup[2,1] = pref*Oo*dp*Yt
    dTup[2,2] = pref/sp.Integer(27)*ds*TFtt_up
    dTup[3,3] = pref/sp.Integer(27)*ds*TFpp_up

    g0 = geom["g0"]
    dTcov = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            val = 0
            for c in range(4):
                for d in range(4):
                    val += g0[a,c]*g0[b,d]*dTup[c,d]
            dTcov[a,b] = sp.expand(val)

    corrected = {
        "T1up": base["T1up"] + dTup,
        "T1cov": base["T1cov"] + dTcov,
        "moments": dict(base["moments"]),
    }

    eq = phase_ix.projected_equations(geom, corrected)
    names = [x[0] for x in eq]
    M = phase_ix.coefficient_matrix(eq, geom["symbols"])
    Mfun = sp.lambdify((w,lam,chiX,chiJ,tau), M, "numpy")

    print("[setup] generic matrix ready", flush=True)
    return Mfun, names


def row_normalized(M):
    A = np.asarray(M, dtype=np.complex128).copy()
    for i in range(A.shape[0]):
        n = np.linalg.norm(A[i])
        if n > 0:
            A[i] /= n
    return A


def nullvector_full(M):
    A = row_normalized(M)
    U,s,Vh = np.linalg.svd(A)
    v = Vh.conj().T[:,-1]
    k = int(np.argmax(np.abs(v)))
    if abs(v[k]) > 0:
        v /= v[k]
    return v, float(np.linalg.norm(A@v)), float(s[-1]), float(s[0])


def selected_minor_indices(row_names):
    preferred = ["TX", "XX", "Xtheta", "AngTrace"]
    return [row_names.index(name) for name in preferred]


def solve_lambda_one(Mfun, row_names, w, chiX, chiJ, tau, seed):
    idx = selected_minor_indices(row_names)
    re_bounds = (-0.80, 0.45)
    im_bounds = (-0.80, 0.80)

    def F(x):
        lam = complex(float(x[0]), float(x[1]))
        M = np.asarray(Mfun(w, lam, chiX, chiJ, tau), dtype=np.complex128)
        A = M[idx,:].copy()
        for i in range(4):
            n = np.linalg.norm(A[i])
            if n > 0:
                A[i] /= n
        d = np.linalg.det(A)
        if not (np.isfinite(d.real) and np.isfinite(d.imag)):
            return np.array([1e3,1e3])
        return np.array([d.real,d.imag])

    offsets = [
        0j, 0.005, -0.005, 0.005j, -0.005j,
        0.015, -0.015, 0.015j, -0.015j,
        0.025+0.025j, 0.025-0.025j, -0.025+0.025j, -0.025-0.025j,
    ]

    candidates = []
    for off in offsets:
        s0 = seed + off
        x0 = np.array([
            np.clip(s0.real, re_bounds[0]+1e-8, re_bounds[1]-1e-8),
            np.clip(s0.imag, im_bounds[0]+1e-8, im_bounds[1]-1e-8),
        ])

        ans = least_squares(
            F,
            x0,
            bounds=([re_bounds[0],im_bounds[0]],[re_bounds[1],im_bounds[1]]),
            xtol=2e-13,
            ftol=2e-13,
            gtol=2e-13,
            max_nfev=400,
            x_scale="jac",
        )

        lam = complex(ans.x[0], ans.x[1])
        if not (np.isfinite(lam.real) and np.isfinite(lam.imag)):
            continue

        M = np.asarray(Mfun(w, lam, chiX, chiJ, tau), dtype=np.complex128)
        v,res,smin,smax = nullvector_full(M)
        rel = smin/max(smax,1e-300)
        minor_res = float(np.linalg.norm(F([lam.real,lam.imag])))

        if ans.success and abs(lam) < 1.0 and rel < 3e-7 and minor_res < 3e-7:
            candidates.append((abs(lam-seed), rel, minor_res, lam, v, res, smin, smax))

    if not candidates:
        raise FloatingPointError(f"No bounded physical lambda root near seed={seed} at w={w}.")

    candidates.sort(key=lambda x:(x[0],x[1],x[2]))
    return candidates[0]


def validate_lambda_candidate(Mfun, w, chiX, chiJ, tau, lam):
    M = np.asarray(Mfun(w, lam, chiX, chiJ, tau), dtype=np.complex128)
    v,res,smin,smax = nullvector_full(M)
    return {
        "lambda": complex(lam),
        "vector": v,
        "residual": float(res),
        "sv_min": float(smin),
        "sv_max": float(smax),
        "relative_sv": float(smin/max(smax,1e-300)),
    }


def smooth_modes(xivc, Mfun, row_names, w, radial, carrier, tau):
    if tau == 0:
        chiX = 0j
        chiJ = 0j
    else:
        ch = smooth_chi(w, radial, carrier)
        chiX = ch["chiX"]
        chiJ = ch["chiJ"]

    sharp_lams = xivc.lambdas_exact(w)
    correction_seeds = [
        C8_LAM_TARGETS[0] - xivc.lambdas_exact(W_SHARP)[0],
        C8_LAM_TARGETS[1] - xivc.lambdas_exact(W_SHARP)[1],
    ]

    global ACTIVE_LAMBDA_SEED_PAIR
    if tau > 0 and ACTIVE_LAMBDA_SEED_PAIR is not None:
        seeds = [complex(ACTIVE_LAMBDA_SEED_PAIR[0]), complex(ACTIVE_LAMBDA_SEED_PAIR[1])]
    else:
        seeds = [sharp_lams[i] + tau*correction_seeds[i] for i in range(2)]

    r0 = solve_lambda_one(Mfun, row_names, w, chiX, chiJ, tau, seeds[0])
    cand0 = validate_lambda_candidate(Mfun, w, chiX, chiJ, tau, complex(r0[3]))
    partner = -2.0/9.0 - cand0["lambda"]
    cand1 = validate_lambda_candidate(Mfun, w, chiX, chiJ, tau, partner)

    if cand1["relative_sv"] > 3e-6:
        r1 = solve_lambda_one(Mfun, row_names, w, chiX, chiJ, tau, seeds[1])
        cand1 = validate_lambda_candidate(Mfun, w, chiX, chiJ, tau, complex(r1[3]))

    candidates = [cand0,cand1]
    direct = abs(candidates[0]["lambda"]-seeds[0]) + abs(candidates[1]["lambda"]-seeds[1])
    swap = abs(candidates[1]["lambda"]-seeds[0]) + abs(candidates[0]["lambda"]-seeds[1])
    if swap < direct:
        candidates = [candidates[1],candidates[0]]

    lams = np.asarray([candidates[0]["lambda"], candidates[1]["lambda"]], dtype=np.complex128)
    relsv = [candidates[0]["relative_sv"], candidates[1]["relative_sv"]]

    if max(relsv) > 3e-6:
        raise FloatingPointError(f"Poor common 7x4 lambda pair: {lams}")
    if abs(np.sum(lams)+2.0/9.0) > 5e-5:
        raise FloatingPointError("Lambda-pair identity lost.")

    return {
        "lambdas": lams,
        "vectors": [candidates[0]["vector"], candidates[1]["vector"]],
        "relative_sv": relsv,
        "chiX": chiX,
        "chiJ": chiJ,
        "lambda_separation": float(abs(lams[0]-lams[1])),
        "lambda_center": 0.5*np.sum(lams),
    }


def direct_smooth_column(w, v):
    H0,H2,J,K = v
    Z = K - 1j*J/(9*w)
    ZX = K/27 + 26j*J/(243*w)
    return np.asarray([Z,ZX], dtype=np.complex128)


def _normalize_state_column(z):
    z = np.asarray(z, dtype=np.complex128).copy()
    if abs(z[0]) > 1e-12:
        return z/z[0]
    if abs(z[1]) > 1e-12:
        return z/z[1]
    raise FloatingPointError("Vanishing interface state column.")


def ocean_basis(xivc, Mfun, row_names, w, radial, carrier, tau, edge_t):
    modes = smooth_modes(xivc, Mfun, row_names, w, radial, carrier, tau)
    cols = []
    for lam,v in zip(modes["lambdas"], modes["vectors"]):
        direct = direct_smooth_column(w, v)
        sharp_col, _ = xivc.interface_column(w, lam, v)
        col = (1.0-edge_t)*sharp_col + edge_t*direct
        cols.append(_normalize_state_column(col))
    C = np.column_stack(cols)
    return C, modes, float(np.linalg.cond(C))


def _sinh_over_delta(delta, L):
    x = delta*L
    if abs(x) < 1e-6:
        return L + delta**2*L**3/6.0 + delta**4*L**5/120.0
    return np.sinh(x)/delta


def _generalized_boundary_derivative(xivc, Mfun, row_names, w, chiX, chiJ, tau, edge_t, lamc):
    eps = 2e-5
    cols = []
    for lam in [lamc+eps, lamc-eps]:
        M = np.asarray(Mfun(w, lam, chiX, chiJ, tau), dtype=np.complex128)
        v,_,_,_ = nullvector_full(M)
        direct = direct_smooth_column(w, v)
        sharp_col, _ = xivc.interface_column(w, lam, v)
        cols.append(_normalize_state_column((1.0-edge_t)*sharp_col + edge_t*direct))
    return (cols[0]-cols[1])/(2.0*eps)


def ocean_transfer(xivc, Mfun, row_names, w, radial, carrier, q, tau, edge_t):
    if abs(q-1.0) < 1e-15:
        return np.eye(2,dtype=np.complex128), {
            "det_T":1+0j, "expected_det":1+0j, "C_cond":1.0,
            "confluent_basis_cond":1.0, "lambda_separation":0.0,
            "lambdas":np.array([0j,0j]), "relative_sv":[0.0,0.0],
            "chiX":0j, "chiJ":0j,
        }

    C,modes,ordinary_cond = ocean_basis(xivc, Mfun, row_names, w, radial, carrier, tau, edge_t)
    L = 9.0*math.log(1.0/q)
    lams = modes["lambdas"]
    lamc = 0.5*np.sum(lams)
    delta = 0.5*(lams[0]-lams[1])
    z1,z2 = C[:,0],C[:,1]
    zc = 0.5*(z1+z2)

    if abs(delta) > 2e-7:
        zd = (z1-z2)/(2.0*delta)
    else:
        zd = _generalized_boundary_derivative(
            xivc,Mfun,row_names,w,modes["chiX"],modes["chiJ"],tau,edge_t,lamc
        )

    B = np.column_stack([zc,zd])
    bcond = np.linalg.cond(B)
    if not np.isfinite(bcond) or bcond > 1e11:
        raise FloatingPointError(f"Ill-conditioned confluent ocean basis: {bcond:.3e}")

    c = np.cosh(delta*L)
    sod = _sinh_over_delta(delta,L)
    ds = delta*np.sinh(delta*L)
    Qc = np.exp(lamc*L)*np.asarray([[c,sod],[ds,c]],dtype=np.complex128)
    T = B@Qc@np.linalg.inv(B)
    expected = np.exp(np.sum(lams)*L)
    detT = np.linalg.det(T)

    rel_det_error = abs(detT-expected)/max(abs(expected),1e-300)
    if rel_det_error > 2e-6:
        raise FloatingPointError("Ocean transfer determinant identity failed.")

    return T, {
        "det_T":detT,
        "expected_det":expected,
        "C_cond":ordinary_cond,
        "confluent_basis_cond":float(bcond),
        "lambda_separation":float(abs(lams[0]-lams[1])),
        "lambdas":lams,
        "relative_sv":modes["relative_sv"],
        "chiX":modes["chiX"],
        "chiJ":modes["chiJ"],
    }


class SmoothGlobalDet:
    def __init__(self, xivc, Mfun, row_names, radial, carrier, *, q, tau=1.0, edge_t=1.0,
                 Rmax=80.0, horizon_eps=1e-7, theta_offset=0.0, asymptotic_order=8):
        self.xivc=xivc; self.Mfun=Mfun; self.row_names=row_names
        self.radial=radial; self.carrier=carrier
        self.q=float(q); self.tau=float(tau); self.edge_t=float(edge_t)
        self.Rmax=float(Rmax); self.horizon_eps=float(horizon_eps)
        self.theta_offset=float(theta_offset); self.asymptotic_order=int(asymptotic_order)
        self.last_info=None; self.last_exception=""

    def raw(self,w):
        w=complex(w)
        zin=self.xivc.vacuum_inner_state(w,eps=self.horizon_eps,rtol=2e-10,atol=2e-12)
        zout=self.xivc.vacuum_outer_state(
            w,Rmax=self.Rmax,rtol=2e-10,atol=2e-12,
            asymptotic_order=self.asymptotic_order,theta_offset=self.theta_offset
        )
        T,info=ocean_transfer(
            self.xivc,self.Mfun,self.row_names,w,self.radial,self.carrier,
            self.q,self.tau,self.edge_t
        )
        zpred=T@zin
        D=zpred[0]*zout[1]-zpred[1]*zout[0]
        scale=max(np.linalg.norm(zpred)*np.linalg.norm(zout),1e-300)
        Dn=D/scale
        info.update({"z_inner":zin,"z_outer":zout,"z_pred":zpred,"normalized_det":Dn})
        self.last_info=info
        return Dn

    def real2(self,x):
        w=complex(float(x[0]),float(x[1]))
        try:
            d=self.raw(w)
            if not (np.isfinite(d.real) and np.isfinite(d.imag)):
                return np.array([1e3,1e3])
            return np.array([d.real,d.imag])
        except Exception as exc:
            self.last_exception=repr(exc)
            return np.array([1e3,1e3])


def solve_global(det, seed, *, max_nfev=160):
    trial_seeds=[seed,seed+0.001,seed-0.001,seed+0.001j,seed-0.001j,seed+0.003-0.002j,seed-0.003-0.002j]
    cand=[]
    for s in trial_seeds:
        ans=least_squares(
            det.real2,[s.real,min(s.imag,-1e-5)],
            bounds=([0.20,-0.18],[0.55,-1e-6]),
            xtol=2e-11,ftol=2e-11,gtol=2e-11,max_nfev=max_nfev
        )
        w=complex(ans.x[0],ans.x[1])
        try:
            d=det.raw(w); ad=abs(d); info=det.last_info
        except Exception:
            ad=np.inf; info=None
        cand.append((ad,abs(w-seed),w,ans,info))
    cand.sort(key=lambda x:(x[0],x[1]))
    return cand[0]


def winding_number(det, center, radius=4e-4, n=40):
    vals=[]
    for k in range(n+1):
        th=2.0*math.pi*k/n
        vals.append(det.raw(center+radius*cmath.exp(1j*th)))
    phases=np.unwrap(np.angle(vals))
    return float((phases[-1]-phases[0])/(2.0*math.pi))


def finite_difference_derivative(det, w, h):
    return (det.raw(w+h)-det.raw(w-h))/(2*h)
