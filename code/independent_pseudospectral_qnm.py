#!/usr/bin/env python3
import argparse
import io
import json
import math
import shutil
import traceback
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import sympy as sp
from scipy.optimize import least_squares
import hbh_phase_ix_bulk_EV as phase_ix
import hbh_phase_xxic8_landau_continuation as landau
Q = 0.1
L_OCEAN = 9.0 * math.log(1.0 / Q)
REFERENCE_W = 0.37751481446 - 0.020871654486j
ROOT_SEEDS = [0.378 - 0.021j, 0.373 - 0.018j, 0.383 - 0.025j]
LAMBDA_SEEDS = [-0.094 + 0.0734j, -0.1282 - 0.0734j]

def zip_folder(out):
    zp = Path(str(out) + '.zip')
    if zp.exists():
        zp.unlink()
    with zipfile.ZipFile(zp, 'w', zipfile.ZIP_DEFLATED) as zf:
        for p in out.rglob('*'):
            if p.is_file():
                zf.write(p, p.relative_to(out.parent))
    return zp

def load_publication_spectral_cache(path):
    p = Path(path)
    with zipfile.ZipFile(p, 'r') as zf:
        result_names = [n for n in zf.namelist() if n.endswith('consistent_q_family_results.csv')]
        cache_names = [n for n in zf.namelist() if n.endswith('spectral_cache/q_0p100_N80.npz')]
        if len(result_names) != 1 or len(cache_names) != 1:
            raise RuntimeError('Could not locate the q=0.1 publication result/cache.')
        results = pd.read_csv(zf.open(result_names[0]))
        row = results[np.isclose(results.q, Q)]
        if len(row) != 1:
            raise RuntimeError('Expected one publication q=0.1 row.')
        stored = row.iloc[0].to_dict()
        raw = io.BytesIO(zf.read(cache_names[0]))
        with np.load(raw, allow_pickle=False) as npz:
            spc = {k: np.array(npz[k]) for k in npz.files}
    radial, carrier = landau.build_catalog(spc)
    return (radial, carrier, stored)

def row_normalized(A):
    A = np.asarray(A, dtype=np.complex128).copy()
    for i in range(A.shape[0]):
        n = np.linalg.norm(A[i])
        if n > 0:
            A[i] /= n
    return A

def nullvector_full(M):
    A = row_normalized(M)
    _, s, Vh = np.linalg.svd(A)
    v = Vh.conj().T[:, -1]
    k = int(np.argmax(np.abs(v)))
    if abs(v[k]) > 0:
        v /= v[k]
    return {'vector': v, 'relative_sv': float(s[-1] / max(s[0], 1e-300)), 'residual': float(np.linalg.norm(A @ v))}

def build_bulk_matrix():
    print('[setup] building independent generic 7x4 bulk matrix...', flush=True)
    omega, lam = sp.symbols('omega lambda')
    chiX, chiJ = sp.symbols('chi_X chi_J')
    geom = phase_ix.build_geometry(omega, lam)
    base = phase_ix.build_stress(geom, omega, lam)
    H0, H2, J, K = geom['symbols']
    S = H0 + K
    SX = lam * S
    sharp_kin = phase_ix.kinetic_coefficients(omega)
    fSX = sharp_kin['f_SX']
    df = (chiX - fSX) * SX + chiJ * J
    O0 = sp.sqrt(sp.Rational(1, 27))
    dp = -sp.I * omega * df / (6 * O0)
    de = (lam + sp.Rational(2, 9)) * df / (6 * O0)
    ds = -sp.Rational(9, 2) * sp.I * omega * (lam + sp.Rational(2, 9)) * df
    theta = geom['theta']
    Y = geom['Y']
    Yt = geom['Ytheta']
    Ytt = geom['Ythetatheta']
    s2 = sp.sin(theta) ** 2
    ct = sp.cos(theta)
    pref = sp.Rational(1, 36) / sp.pi
    TFtt_cov = Ytt + 3 * Y
    TFpp_cov = sp.sin(theta) * ct * Yt + 3 * s2 * Y
    dTup = sp.zeros(4)
    dTup[0, 1] = dTup[1, 0] = pref * df * Y
    dTup[0, 2] = dTup[2, 0] = pref * O0 * de * Yt
    dTup[1, 2] = dTup[2, 1] = pref * O0 * dp * Yt
    dTup[2, 2] = pref * ds * TFtt_cov / sp.Integer(27)
    dTup[3, 3] = pref * ds * (TFpp_cov / s2 ** 2) / sp.Integer(27)
    g0 = geom['g0']
    dTcov = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            val = 0
            for c in range(4):
                for d in range(4):
                    val += g0[a, c] * g0[b, d] * dTup[c, d]
            dTcov[a, b] = sp.expand(val)
    corrected = {'T1up': base['T1up'] + dTup, 'T1cov': base['T1cov'] + dTcov, 'moments': dict(base['moments'])}
    equations = phase_ix.projected_equations(geom, corrected)
    names = [e[0] for e in equations]
    M = phase_ix.coefficient_matrix(equations, geom['symbols'])
    Mfun = sp.lambdify((omega, lam, chiX, chiJ), M, 'numpy')
    return (Mfun, names)

def selected_minor_indices(row_names):
    preferred = ['TX', 'XX', 'Xtheta', 'AngTrace']
    return [row_names.index(name) for name in preferred]

def solve_one_lambda(Mfun, row_names, omega, chiX, chiJ, seed):
    idx = selected_minor_indices(row_names)
    def residual(x):
        lam = complex(float(x[0]), float(x[1]))
        M = np.asarray(Mfun(omega, lam, chiX, chiJ), dtype=np.complex128)
        A = row_normalized(M[idx, :])
        d = np.linalg.det(A)
        if not np.isfinite(d.real) or not np.isfinite(d.imag):
            return np.array([1000.0, 1000.0])
        return np.array([d.real, d.imag])
    seed_cloud = [seed, seed + 0.006, seed - 0.006, seed + 0.006j, seed - 0.006j, seed + 0.015 + 0.01j, seed - 0.015 - 0.01j]
    candidates = []
    for s in seed_cloud:
        ans = least_squares(residual, [s.real, s.imag], bounds=([-0.8, -0.8], [0.45, 0.8]), xtol=2e-12, ftol=2e-12, gtol=2e-12, max_nfev=250)
        lam = complex(ans.x[0], ans.x[1])
        M = np.asarray(Mfun(omega, lam, chiX, chiJ), dtype=np.complex128)
        nv = nullvector_full(M)
        minor_res = float(np.linalg.norm(residual([lam.real, lam.imag])))
        if ans.success and abs(lam) < 1.0 and (nv['relative_sv'] < 2e-06) and (minor_res < 2e-06):
            candidates.append((abs(lam - seed), nv['relative_sv'], minor_res, lam, nv['vector']))
    if not candidates:
        raise RuntimeError(f'No physical lambda root found near {seed} at omega={omega}.')
    candidates.sort(key=lambda x: (x[0], x[1], x[2]))
    _, relsv, minor_res, lam, vector = candidates[0]
    return {'lambda': lam, 'vector': vector, 'relative_sv': relsv, 'minor_residual': minor_res}

def bulk_modes(Mfun, row_names, omega, radial, carrier):
    ch = landau.evaluate_chi(omega, radial, carrier, second_sheet=True)
    m1 = solve_one_lambda(Mfun, row_names, omega, ch['chiX'], ch['chiJ'], LAMBDA_SEEDS[0])
    m2 = solve_one_lambda(Mfun, row_names, omega, ch['chiX'], ch['chiJ'], LAMBDA_SEEDS[1])
    if abs(m1['lambda'] - m2['lambda']) < 0.0001:
        raise RuntimeError('The two independently solved bulk roots collapsed.')
    return ([m1, m2], ch)

def interface_state(omega, vector):
    H0, H2, J, K = vector
    Z = K - 1j * J / (9 * omega)
    ZX = K / 27 + 26j * J / (243 * omega)
    return np.asarray([Z, ZX], dtype=np.complex128)

def cheb_lobatto(n, a, b):
    if n < 4:
        raise ValueError('n must be >= 4')
    j = np.arange(n + 1)
    z = np.cos(np.pi * j / n)
    c = np.ones(n + 1)
    c[0] = 2
    c[-1] = 2
    c = c * (-1.0) ** j
    Z = np.tile(z, (n + 1, 1))
    dZ = Z.T - Z
    D = np.outer(c, 1 / c) / (dZ + np.eye(n + 1))
    D = D - np.diag(np.sum(D, axis=1))
    x = a + (b - a) * (z + 1) / 2
    D = 2 * D / (b - a)
    p = np.arange(n, -1, -1)
    x = x[p]
    D = D[np.ix_(p, p)]
    return (x, D)

def solve_overdetermined(A, b):
    A = np.asarray(A, dtype=np.complex128)
    b = np.asarray(b, dtype=np.complex128)
    for i in range(A.shape[0]):
        n = np.linalg.norm(A[i])
        if n > 0:
            A[i] /= n
            b[i] /= n
    u, *_ = np.linalg.lstsq(A, b, rcond=None)
    rel = np.linalg.norm(A @ u - b) / max(np.linalg.norm(b), 1.0)
    return (u, float(rel))

def Vz(R):
    f = 1.0 - 2.0 / R
    return 2 * f * (12 * R ** 3 + 12 * R ** 2 + 18 * R + 9) / (R ** 3 * (2 * R + 3) ** 2)

def inner_robin(omega, n):
    y, D = cheb_lobatto(n, 0.0, 1.0)
    R = 2.0 + y
    D1 = D
    D2 = D1 @ D1
    f = 1.0 - 2.0 / R
    fp = 2.0 / R ** 2
    V = Vz(R)
    Arows = []
    brows = []
    row = np.zeros(n + 1, dtype=np.complex128)
    row[0] = 1.0
    Arows.append(row)
    brows.append(1.0 + 0j)
    alpha = 27.0 / 28.0 / (0.5 - 2j * omega)
    Arows.append(D1[0] - alpha * np.eye(n + 1)[0])
    brows.append(0.0 + 0j)
    for k in range(1, n + 1):
        row = f[k] ** 2 * D2[k] + f[k] * (fp[k] - 2j * omega) * D1[k] - V[k] * np.eye(n + 1)[k]
        Arows.append(row)
        brows.append(0.0 + 0j)
    u, rel = solve_overdetermined(np.vstack(Arows), np.asarray(brows))
    u3 = u[-1]
    uR3 = D1[-1] @ u
    if abs(u3) < 1e-13:
        raise RuntimeError('Inner factored solution vanished at R=3.')
    f3 = 1.0 / 3.0
    y3 = -1j * omega + f3 * uR3 / u3
    return (complex(y3), rel)

def outer_robin(omega, n):
    x, D = cheb_lobatto(n, 0.0, 1.0 / 3.0)
    D1 = D
    D2 = D1 @ D1
    R = np.empty_like(x)
    R[0] = np.inf
    R[1:] = 1.0 / x[1:]
    f = 1.0 - 2.0 * x
    V_over_x2 = np.empty_like(x)
    V_over_x2[0] = 6.0
    V_over_x2[1:] = Vz(R[1:]) / x[1:] ** 2
    Arows = []
    brows = []
    row = np.zeros(n + 1, dtype=np.complex128)
    row[0] = 1.0
    Arows.append(row)
    brows.append(1.0 + 0j)
    alpha = 3j / omega
    Arows.append(D1[0] - alpha * np.eye(n + 1)[0])
    brows.append(0.0 + 0j)
    for k in range(1, n + 1):
        row = x[k] ** 2 * f[k] ** 2 * D2[k] + (2 * x[k] * f[k] ** 2 - 2 * f[k] * (x[k] ** 2 + 1j * omega)) * D1[k] - V_over_x2[k] * np.eye(n + 1)[k]
        Arows.append(row)
        brows.append(0.0 + 0j)
    u, rel = solve_overdetermined(np.vstack(Arows), np.asarray(brows))
    ued = u[-1]
    ux = D1[-1] @ u
    if abs(ued) < 1e-13:
        raise RuntimeError('Outer factored solution vanished at R=3.')
    x3 = 1.0 / 3.0
    uR3 = -x3 ** 2 * ux
    f3 = 1.0 / 3.0
    y3 = 1j * omega + f3 * uR3 / ued
    return (complex(y3), rel)

class PseudospectralMatcher:
    def __init__(self, Mfun, row_names, radial, carrier, ncoll):
        self.Mfun = Mfun
        self.row_names = row_names
        self.radial = radial
        self.carrier = carrier
        self.ncoll = int(ncoll)
        self.last = None
    def value(self, omega):
        omega = complex(omega)
        yH, rin = inner_robin(omega, self.ncoll)
        yO, rout = outer_robin(omega, self.ncoll)
        modes, ch = bulk_modes(self.Mfun, self.row_names, omega, self.radial, self.carrier)
        z = [interface_state(omega, m['vector']) for m in modes]
        lams = [m['lambda'] for m in modes]
        h = [zz[1] - yH * zz[0] for zz in z]
        o = [(zz[1] - yO * zz[0]) * np.exp(lam * L_OCEAN) for zz, lam in zip(z, lams)]
        Fraw = h[0] * o[1] - h[1] * o[0]
        scale = max(np.linalg.norm(h) * np.linalg.norm(o), 1e-300)
        F = Fraw / scale
        self.last = {'omega': omega, 'inner_robin': yH, 'outer_robin': yO, 'inner_collocation_residual': rin, 'outer_collocation_residual': rout, 'lambda1': lams[0], 'lambda2': lams[1], 'lambda_sum_error': abs(lams[0] + lams[1] + 2.0 / 9.0), 'max_bulk_relative_sv': max((m['relative_sv'] for m in modes)), 'chiJ': ch['chiJ'], 'chiX': ch['chiX'], 'matching_abs': abs(F)}
        return F
    def real2(self, x):
        omega = complex(float(x[0]), float(x[1]))
        try:
            F = self.value(omega)
            if not np.isfinite(F.real) or not np.isfinite(F.imag):
                return np.array([1000.0, 1000.0])
            return np.array([F.real, F.imag])
        except Exception:
            return np.array([1000.0, 1000.0])

def solve_root(matcher, seed):
    ans = least_squares(matcher.real2, [seed.real, seed.imag], bounds=([0.34, -0.06], [0.42, -0.002]), xtol=5e-11, ftol=5e-11, gtol=5e-11, max_nfev=100)
    w = complex(ans.x[0], ans.x[1])
    F = matcher.value(w)
    return (w, F, ans, dict(matcher.last))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--family', default='phase_prx_final_q_family_output.zip')
    ap.add_argument('--outdir', default='phase_independent_pseudospectral_qnm_output')
    ap.add_argument('--quick', action='store_true')
    args = ap.parse_args()
    out = Path(args.outdir)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    radial, carrier, stored = load_publication_spectral_cache(args.family)
    stored_w = complex(float(stored['root_re']), float(stored['root_im']))
    if abs(stored_w - REFERENCE_W) > 5e-06:
        raise RuntimeError(f'The supplied family ZIP does not contain the publication-resolution q=0.1 root. Its stored root is {stored_w}. Use the FULL phase_prx_final_q_family_output.zip.')
    Mfun, row_names = build_bulk_matrix()
    resolutions = [32, 48, 64] if args.quick else [48, 64, 80]
    rows = []
    seed_rows = []
    errors = []
    print('=' * 78, flush=True)
    print('INDEPENDENT PSEUDOSPECTRAL q=0.1 QNM REPRODUCTION', flush=True)
    print(f'resolutions = {resolutions}', flush=True)
    print('=' * 78, flush=True)
    previous = ROOT_SEEDS[0]
    for ncoll in resolutions:
        print(f'\n[Ncoll={ncoll}]', flush=True)
        matcher = PseudospectralMatcher(Mfun, row_names, radial, carrier, ncoll)
        try:
            w, F, ans, info = solve_root(matcher, previous)
            previous = w
            row = {'Ncoll': ncoll, 'root_re': w.real, 'root_im': w.imag, 'matching_abs': abs(F), 'difference_from_reference': abs(w - REFERENCE_W), 'inner_collocation_residual': info['inner_collocation_residual'], 'outer_collocation_residual': info['outer_collocation_residual'], 'lambda1_re': info['lambda1'].real, 'lambda1_im': info['lambda1'].imag, 'lambda2_re': info['lambda2'].real, 'lambda2_im': info['lambda2'].imag, 'lambda_sum_error': info['lambda_sum_error'], 'max_bulk_relative_sv': info['max_bulk_relative_sv'], 'chiJ_abs': abs(info['chiJ']), 'chiX_abs': abs(info['chiX']), 'solver_success': bool(ans.success)}
            rows.append(row)
            print(f'  Mw={w.real:.12f}{w.imag:+.12f}i  |F|={abs(F):.3e}  |dw_ref|={abs(w - REFERENCE_W):.3e}', flush=True)
        except Exception as exc:
            errors.append({'Ncoll': ncoll, 'stage': 'resolution', 'error': repr(exc), 'traceback': traceback.format_exc()})
            print(f'  FAIL: {exc}', flush=True)
            break
    df = pd.DataFrame(rows)
    df.to_csv(out / 'pseudospectral_resolution_convergence.csv', index=False)
    if len(df) and not args.quick:
        nbest = int(df.iloc[-1].Ncoll)
        print(f'\n[independent seed audit at Ncoll={nbest}]', flush=True)
        for seed in ROOT_SEEDS:
            matcher = PseudospectralMatcher(Mfun, row_names, radial, carrier, nbest)
            try:
                w, F, ans, info = solve_root(matcher, seed)
                seed_rows.append({'seed_re': seed.real, 'seed_im': seed.imag, 'root_re': w.real, 'root_im': w.imag, 'matching_abs': abs(F), 'difference_from_reference': abs(w - REFERENCE_W), 'max_bulk_relative_sv': info['max_bulk_relative_sv'], 'inner_collocation_residual': info['inner_collocation_residual'], 'outer_collocation_residual': info['outer_collocation_residual']})
                print(f'  seed={seed.real:.3f}{seed.imag:+.3f}i -> {w.real:.12f}{w.imag:+.12f}i', flush=True)
            except Exception as exc:
                errors.append({'Ncoll': nbest, 'stage': 'seed_audit', 'seed': str(seed), 'error': repr(exc), 'traceback': traceback.format_exc()})
    sdf = pd.DataFrame(seed_rows)
    sdf.to_csv(out / 'pseudospectral_seed_audit.csv', index=False)
    if len(df) >= 3:
        final = df.iloc[-1]
        prev = df.iloc[-2]
        resolution_shift = abs(complex(final.root_re, final.root_im) - complex(prev.root_re, prev.root_im))
        match_reference = final.difference_from_reference < 0.0001
        converged_resolution = resolution_shift < 7.5e-05
        matching_good = final.matching_abs < 2e-06
        bulk_good = final.max_bulk_relative_sv < 3e-06 and final.lambda_sum_error < 2e-05
        collocation_good = final.inner_collocation_residual < 5e-07 and final.outer_collocation_residual < 5e-07
        if args.quick:
            seed_good = True
        else:
            seed_good = len(sdf) == len(ROOT_SEEDS) and max(abs(complex(r.root_re, r.root_im) - complex(final.root_re, final.root_im)) for _, r in sdf.iterrows()) < 0.0001
        passed = all([match_reference, converged_resolution, matching_good, bulk_good, collocation_good, seed_good])
    else:
        resolution_shift = np.inf
        match_reference = False
        converged_resolution = False
        matching_good = False
        bulk_good = False
        collocation_good = False
        seed_good = False
        passed = False
    if passed and args.quick:
        verdict = 'QUICK_PSEUDOSPECTRAL_PASS__RUN_FULL'
    elif passed:
        verdict = 'INDEPENDENT_PSEUDOSPECTRAL_REPRODUCTION_PASS'
    else:
        verdict = 'INDEPENDENT_PSEUDOSPECTRAL_REPRODUCTION_NOT_ESTABLISHED'
    summary = {'verdict': verdict, 'reference_root': [REFERENCE_W.real, REFERENCE_W.imag], 'resolutions': resolutions, 'final_resolution_shift': float(resolution_shift), 'match_reference_within_1e-4': bool(match_reference), 'resolution_converged': bool(converged_resolution), 'matching_residual_gate': bool(matching_good), 'bulk_7x4_gate': bool(bulk_good), 'pseudospectral_BVP_gate': bool(collocation_good), 'three_seed_same_root_gate': bool(seed_good), 'scope_note': 'This is an independent numerical-method validation of the global boundary-value problem, using the same physical Einstein--Vlasov constitutive response. It is not an external reproduction by another group.'}
    (out / 'INDEPENDENT_PSEUDOSPECTRAL_SUMMARY.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    if errors:
        pd.DataFrame(errors).to_csv(out / 'failed_cases.csv', index=False)
    lines = ['INDEPENDENT PSEUDOSPECTRAL q=0.1 QNM REPRODUCTION', '=' * 78, '', f'VERDICT: {verdict}', '']
    if len(df):
        f = df.iloc[-1]
        lines += [f'highest-resolution root = {f.root_re:.12f}{f.root_im:+.12f} i', f'difference from C9 publication root = {f.difference_from_reference:.6e}', f'last resolution shift = {resolution_shift:.6e}', f'collocation matching |F| = {f.matching_abs:.6e}', f'inner/outer BVP residuals = {f.inner_collocation_residual:.3e}, {f.outer_collocation_residual:.3e}', f'max full 7x4 relative SV = {f.max_bulk_relative_sv:.3e}', '']
    lines += ['INTERPRETATION', '-' * 78, 'The vacuum boundary conditions are computed by independent Chebyshev pseudospectral boundary-value solvers after explicit ingoing/outgoing phase factorization.', 'No C9 transfer matrix, C9 global determinant, C9 homotopy tracker, or ECS root solver is used.', 'A PASS is therefore an independent numerical-method reproduction of the global q=0.1 pole, while retaining the same physical kinetic constitutive response.']
    report = '\n'.join(lines)
    (out / 'INDEPENDENT_PSEUDOSPECTRAL_VERDICT.txt').write_text(report, encoding='utf-8')
    zp = zip_folder(out)
    print('\n' + report, flush=True)
    print('\nDONE', flush=True)
    print('Send back:', zp, flush=True)
if __name__ == '__main__':
    main()
