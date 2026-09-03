"""
Independent derivation of the l=2 even-parity linearized Einstein tensor in
the F(X)=exp(2X/9)/3 background. No supplied Einstein-tensor matrix is
imported.
"""
import sympy as sp


def build_geometry(w, lam):
    T, X, th, ph = sp.symbols("T X theta phi", real=True)
    coords = [T, X, th, ph]
    H0, H2, J, K = sp.symbols("H0 H2 J K")

    F = sp.exp(sp.Rational(2, 9) * X) / 3
    Y = (3 * sp.cos(th) ** 2 - 1) / 2
    E = sp.exp(lam * X - sp.I * w * T)

    g0 = sp.diag(-F, F, 27 * F, 27 * F * sp.sin(th) ** 2)
    h = sp.zeros(4)
    h[0, 0] = F * H0 * E * Y
    h[0, 1] = h[1, 0] = F * J * E * Y
    h[1, 1] = F * H2 * E * Y
    h[2, 2] = 27 * F * K * E * Y
    h[3, 3] = 27 * F * sp.sin(th) ** 2 * K * E * Y

    ginv0 = g0.inv()
    ginv1 = -ginv0 * h * ginv0

    def d_(expr, i):
        return sp.diff(expr, coords[i])

    Gam0 = [[[0] * 4 for _ in range(4)] for _ in range(4)]
    Gam1 = [[[0] * 4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(4):
                acc0 = acc1 = 0
                for e in range(4):
                    dg0 = d_(g0[e, c], b) + d_(g0[e, b], c) - d_(g0[b, c], e)
                    dh = d_(h[e, c], b) + d_(h[e, b], c) - d_(h[b, c], e)
                    acc0 += sp.Rational(1, 2) * ginv0[a, e] * dg0
                    acc1 += sp.Rational(1, 2) * (ginv1[a, e] * dg0 + ginv0[a, e] * dh)
                Gam0[a][b][c] = sp.cancel(acc0)
                Gam1[a][b][c] = sp.expand(acc1)

    R0 = sp.zeros(4)
    R1 = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            v0 = v1 = 0
            for c in range(4):
                v0 += d_(Gam0[c][a][b], c) - d_(Gam0[c][a][c], b)
                v1 += d_(Gam1[c][a][b], c) - d_(Gam1[c][a][c], b)
                for e in range(4):
                    v0 += Gam0[c][a][b] * Gam0[e][c][e] - Gam0[e][a][c] * Gam0[c][b][e]
                    v1 += (
                        Gam1[c][a][b] * Gam0[e][c][e]
                        + Gam0[c][a][b] * Gam1[e][c][e]
                        - Gam1[e][a][c] * Gam0[c][b][e]
                        - Gam0[e][a][c] * Gam1[c][b][e]
                    )
            R0[a, b] = sp.cancel(v0)
            R1[a, b] = sp.expand(v1)

    subs0 = {T: 0, X: 0, th: sp.pi / 3}
    g0e, he = g0.subs(subs0), h.subs(subs0)
    ginv0e, ginv1e = ginv0.subs(subs0), ginv1.subs(subs0)
    R0e, R1e = R0.subs(subs0), R1.subs(subs0)

    Rs0 = sp.cancel(sum(ginv0e[a, b] * R0e[a, b] for a in range(4) for b in range(4)))
    Rs1 = sp.expand(sum(ginv1e[a, b] * R0e[a, b] + ginv0e[a, b] * R1e[a, b]
                        for a in range(4) for b in range(4)))

    G1 = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            G1[a, b] = sp.expand(R1e[a, b] - sp.Rational(1, 2) * he[a, b] * Rs0
                                 - sp.Rational(1, 2) * g0e[a, b] * Rs1)

    Y0 = sp.simplify(Y.subs(th, sp.pi / 3))
    Y1 = sp.simplify(sp.diff(Y, th).subs(th, sp.pi / 3))
    Y2 = sp.simplify(sp.diff(Y, th, 2).subs(th, sp.pi / 3))

    return dict(symbols=(H0, H2, J, K), g0=g0e, h=he, G1=G1,
                Y=Y0, Yth=Y1, Ythth=Y2, theta=sp.pi / 3)


def kinetic_coefficients(w):
    O = sp.sqrt(sp.Rational(1, 27))
    dm, dp = 2 * O - w, 2 * O + w
    return dict(
        O=O,
        d_S=9 * O**2 / (2 * dm * dp),
        d_A=(4 * O**4 + O**2 * w**2 + w**4) / (w**2 * dm**2 * dp**2),
        d_C=sp.I * (O - w) * (O + w) / (w * dm * dp),
        f_SX=sp.I * (O - w) * (O + w) / (2 * w * dm * dp),
        e_S=-3 * sp.I * O * w / (4 * dm * dp),
        e_A=-sp.I * O * w / (dm**2 * dp**2),
        e_C=-O / (2 * dm * dp),
        p_SX=-O / (4 * dm * dp),
        s_S=-3 * O**2 / (8 * dm * dp),
        s_A=O**2 * (4 * O**2 - 3 * w**2) / (4 * w**2 * dm**2 * dp**2),
        s_C=sp.I * O**2 / (4 * w * dm * dp),
    )


def build_stress(geom, w, lam, chiJ, chiX):
    H0, H2, J, K = geom["symbols"]
    theta = geom["theta"]
    Y, Yth, Ythth = geom["Y"], geom["Yth"], geom["Ythth"]
    kin = kinetic_coefficients(w)
    O = kin["O"]

    S = H0 + K
    SX = lam * S
    SXX = lam**2 * S
    A = SXX / 2 + SX / 9
    C = lam * J + 2 * J / 9

    d = sp.expand(kin["d_S"] * S + kin["d_A"] * A + kin["d_C"] * C)
    f_base = sp.expand(kin["f_SX"] * SX)
    e_base = sp.expand(kin["e_S"] * S + kin["e_A"] * A + kin["e_C"] * C)
    p_base = sp.expand(kin["p_SX"] * SX)
    s_base = sp.expand(kin["s_S"] * S + kin["s_A"] * A + kin["s_C"] * C)

    delta_f = (chiX - kin["f_SX"]) * SX + chiJ * J
    delta_p = -sp.I * w * delta_f / (6 * O)
    delta_e = (lam + sp.Rational(2, 9)) * delta_f / (6 * O)
    delta_s = -sp.Rational(9, 2) * sp.I * w * (lam + sp.Rational(2, 9)) * delta_f

    f_tot = f_base + delta_f
    e_tot = e_base + delta_e
    p_tot = p_base + delta_p
    s_tot = s_base + delta_s

    aamp = d + H0 - H2 / 2 - 3 * K / 2
    bamp = f_tot - J
    Damp = d - H2 / 2 - 5 * K / 2

    pref = sp.Rational(1, 36) / sp.pi
    s2 = sp.sin(theta) ** 2
    ct = sp.cos(theta)
    TFtt = Ythth + 3 * Y
    TFpp = (sp.sin(theta) * ct * Yth + 3 * s2 * Y) / s2**2

    T1up = sp.zeros(4)
    T1up[0, 0] = pref * aamp * Y
    T1up[0, 1] = T1up[1, 0] = pref * bamp * Y
    T1up[0, 2] = T1up[2, 0] = pref * O * e_tot * Yth
    T1up[1, 2] = T1up[2, 1] = pref * O * p_tot * Yth
    T1up[2, 2] = pref / sp.Integer(27) * (Damp * Y / 2 + s_tot * TFtt)
    T1up[3, 3] = pref / sp.Integer(27) * (Damp * Y / (2 * s2) + s_tot * TFpp)

    T0up = sp.zeros(4)
    T0up[0, 0] = pref
    T0up[2, 2] = pref / sp.Integer(54)
    T0up[3, 3] = pref / (sp.Integer(54) * s2)

    g0, h = geom["g0"], geom["h"]
    T1cov = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            val = 0
            for c in range(4):
                for dd in range(4):
                    val += g0[a, c] * g0[b, dd] * T1up[c, dd]
                    val += h[a, c] * g0[b, dd] * T0up[c, dd]
                    val += g0[a, c] * h[b, dd] * T0up[c, dd]
            T1cov[a, b] = sp.expand(val)
    return dict(T1cov=T1cov)


def projected_equations(geom, stress):
    G, T = geom["G1"], stress["T1cov"]
    Y, Yth, th = geom["Y"], geom["Yth"], geom["theta"]
    s2 = sp.sin(th) ** 2

    R = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            R[a, b] = sp.expand(G[a, b] - 8 * sp.pi * T[a, b])

    eq = [
        ("TT", sp.cancel(R[0, 0] / Y)),
        ("TX", sp.cancel(R[0, 1] / Y)),
        ("XX", sp.cancel(R[1, 1] / Y)),
        ("Ttheta", sp.cancel(R[0, 2] / Yth)),
        ("Xtheta", sp.cancel(R[1, 2] / Yth)),
    ]
    eq.append(("AngTrace", sp.cancel((R[2, 2] + R[3, 3] / s2) / Y)))
    TFtt = geom["Ythth"] + 3 * Y
    TFpp_over_s2 = (sp.sin(th) * sp.cos(th) * Yth + 3 * s2 * Y) / s2
    TFshape = sp.simplify(TFtt - TFpp_over_s2)
    eq.append(("AngTF", sp.cancel((R[2, 2] - R[3, 3] / s2) / TFshape)))
    return eq


def coefficient_matrix(eq, symbols):
    rows = []
    for name, expr in eq:
        expr = sp.expand(expr)
        row = [sp.expand(sp.diff(expr, z)) for z in symbols]
        rem = sp.cancel(sp.expand(expr - sum(row[j] * symbols[j] for j in range(4))))
        if rem != 0:
            print(f"WARNING nonhomogeneous remainder in {name}: {rem}")
        rows.append(row)
    return sp.Matrix(rows)
