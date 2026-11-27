#!/usr/bin/env python3
"""The growing mode of the c = -1 black universe in closed form, and the slope of the rate through c = -1 (#262).

#212's master equation for Y = dR in the gauge dphi = 0, e^{W t} time dependence, member c of the family (b = 1, K = -2c/pi, 3m = K):

    Y'' = (2/rho - A'/A) Y' + (1/R^4 - A'/(A rho R^2) - 2/(A R^2) + W^2/A^2) Y,     A = 1 + c R^2 + K (rho + R^2 arctan rho),  R^2 = 1 + rho^2.

Result 1. Y = R e^{W x}, x = int d rho / A (the tortoise coordinate), solves it exactly when W = K = -2c/pi, for every c (S1). At c = -1, K = kappa = 2/pi: the
growing mode Y_ref = R e^{kappa x}, Psi = Y_ref e^{-kappa x} = R.  The second solution is a quadrature, Y_inf = Y_ref int_rho^inf (2K rho^2 e^{-2 kappa x}/(A R^2)) d rho, and
has no closed form of the same kind (the formal series of e^{-kappa x} Phi at rho = infinity diverges, N3).
Result 2. The first-order monodromy of the equation about (c, W) = (-1, kappa), computed from residues at rho = 0 (S3), gives d(W/kappa)/dc = -(pi/2) beta*,
beta* = pi/4 - 1/I0 the admixture of the rho^{3/2} solution in the solution decaying at infinity, I0 = int_0^inf 2K rho^2 e^{-2 kappa x}/(A R^2) d rho (N4); so
dW/dc = 1/I0 - 2/pi (N5, N6).

    cminus1_mode.py derive [--mutate=NAME]       sympy checks S1 - S4
    cminus1_mode.py check  [--only=..] [--mutate=NAME]   numerical checks N0 - N7; exit 1 when one fails
    cminus1_mode.py table                        the constants of the note
    cminus1_mode.py mutations DIR PREFIX         fired sets of stored logs DIR/PREFIX<name>.log against the declared ones

Status words as in AGENTS.md. mpmath for every extended-precision number (np.longdouble is float64 on macOS).
"""
import ast
import math
import os
import re
import sys
import time

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
ONLY = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--only=")), None)
RESULTS = []


def want(cid):
    return ONLY is None or cid in ONLY or any(cid.startswith(o) for o in ONLY)


def check(cid, text, ok, detail=""):
    RESULTS.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


# ----------------------------------------------------------------------------------------------------------------------
# truncated power series (lists of mpf)
# ----------------------------------------------------------------------------------------------------------------------
def s_mul(a, b, n):
    out = [mp.mpf(0)] * n
    for i, ai in enumerate(a[:n]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[:n - i]):
            out[i + j] += ai * bj
    return out


def s_inv(a, n):
    out = [mp.mpf(0)] * n
    out[0] = 1 / a[0]
    for m in range(1, n):
        acc = mp.mpf(0)
        for j in range(1, min(m, len(a) - 1) + 1):
            acc += a[j] * out[m - j]
        out[m] = -acc / a[0]
    return out


def s_add(a, b, n):
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def s_scale(a, c, n):
    return [c * (a[i] if i < len(a) else 0) for i in range(n)]


def s_der(a):
    return [i * a[i] for i in range(1, len(a))]


def bg_series0(K, n):
    """A = -rho^2 + K (rho + (1 + rho^2) arctan rho) at rho = 0, n + 2 coefficients"""
    T = [mp.mpf(0)] * (n + 2)
    for m in range(0, n + 2):
        if 2 * m + 1 < n + 2:
            T[2 * m + 1] = mp.mpf((-1) ** m) / (2 * m + 1)
    A = [K * x for x in s_mul([mp.mpf(1), mp.mpf(0), mp.mpf(1)], T, n + 2)]
    A[1] += K
    A[2] -= 1
    return A[:n + 2]


def Aof(K, r):
    return -r * r + K * (r + (1 + r * r) * mp.atan(r))


def Apof(K, r):
    return -2 * r + 2 * K * (1 + r * mp.atan(r))


# ----------------------------------------------------------------------------------------------------------------------
# the equation
# ----------------------------------------------------------------------------------------------------------------------
def PQ(K, r, W):
    """P, Q of Y'' = P Y' + Q Y (mpmath scalars). eq_P: the sign of A'/A in P flipped (mutation)"""
    A = Aof(K, r)
    Ap = Apof(K, r)
    R2 = 1 + r * r
    sg = -1 if MUT == "eq_P" else 1
    P = 2 / r - sg * Ap / A
    Q = 1 / R2 ** 2 - Ap / (A * r * R2) - 2 / (A * R2) + W * W / A ** 2
    return P, Q


def frobenius_s(K, w, s, g1, N):
    """Y = rho^s sum g_k rho^k at W = w K, from the polynomial form of the equation,
    sum_{k+j=n} g_k [C2_j (s+k)(s+k-1) + C1_j (s+k) + C0_j] = 0, with A = rho a, C2 = a^2 R^4, C1 = -a R^4 (a - rho a'),
    C0 = -[rho^2 a^2 - a (a + rho a') R^2 - 2 rho a R^2 + W^2 R^4]. s = 1/2: g1 is the free resonance coefficient; s = 3/2: unique."""
    A = bg_series0(K, N + 4)
    a = A[1:]
    n = N + 3
    ap = s_der(a)
    R2 = [mp.mpf(1), mp.mpf(0), mp.mpf(1)]
    R4 = s_mul(R2, R2, n)
    C2 = s_mul(s_mul(a, a, n), R4, n)
    rap = [mp.mpf(0)] + ap[:n - 1]
    t1 = s_add(a, [-x for x in rap], n)
    C1 = [-x for x in s_mul(s_mul(a, R4, n), t1, n)]
    t2 = s_mul(s_mul(a, s_add(a, rap, n), n), R2, n)
    rho2a2 = [mp.mpf(0), mp.mpf(0)] + s_mul(a, a, n)[:n - 2]
    t3 = [mp.mpf(0)] + s_scale(s_mul(a, R2, n), 2, n)[:n - 1]
    W2R4 = [K * K * w * w * x for x in R4]
    C0 = [-(rho2a2[i] - t2[i] - t3[i] + W2R4[i]) for i in range(n)]
    I = lambda sig: C2[0] * sig * (sig - 1) + C1[0] * sig + C0[0]
    s = mp.mpf(s)
    g = [mp.mpf(1)]
    start = 1
    if s == mp.mpf(1) / 2:
        g.append(mp.mpf(g1))
        start = 2
    for m in range(start, N + 1):
        tot = mp.mpf(0)
        for k in range(m):
            j = m - k
            tot += g[k] * (C2[j] * (s + k) * (s + k - 1) + C1[j] * (s + k) + C0[j])
        g.append(-tot / I(s + m))
    return g


def Y_frob(g, s, rho):
    s = mp.mpf(s)
    ar = abs(rho)
    S = sum(gk * rho ** k for k, gk in enumerate(g))
    dS = sum(k * gk * rho ** (k - 1) for k, gk in enumerate(g) if k > 0)
    return ar ** s * S, ar ** s * (s * S / rho + dS)


def local_series(K, rho0, n):
    r = [rho0, mp.mpf(1)] + [mp.mpf(0)] * (n - 2)
    R2 = [1 + rho0 ** 2, 2 * rho0, mp.mpf(1)] + [mp.mpf(0)] * (n - 3)
    iR2 = s_inv(R2, n)
    T = [mp.atan(rho0)] + [iR2[i - 1] / i for i in range(1, n)]
    A = s_add(s_scale(s_mul(R2, T, n), K, n), s_add(s_scale(r, K, n), [-x for x in s_mul(r, r, n)], n), n)
    return r, R2, iR2, T, A


def step(K, w, rho0, Y0, Yp0, h, n=50):
    """one Taylor step of the Y equation (W = w K) from rho0 to rho0 + h"""
    r, R2, iR2, T, A = local_series(K, rho0, n + 2)
    Ap = s_der(A)
    iA = s_inv(A, n + 2)
    ir = s_inv(r, n + 2)
    iR4 = s_mul(iR2, iR2, n + 2)
    AA = s_mul(Ap, iA, n + 2)
    sg = -1 if MUT == "eq_P" else 1
    P = s_add([2 * x for x in ir], [-sg * x for x in AA], n + 2)
    t1 = s_mul(s_mul(AA, ir, n + 2), iR2, n + 2)
    t2 = s_scale(s_mul(iA, iR2, n + 2), 2, n + 2)
    t3 = s_scale(s_mul(iA, iA, n + 2), K * K * w * w, n + 2)
    Q = [iR4[i] - t1[i] - t2[i] + t3[i] for i in range(n + 2)]
    y = [Y0, Yp0] + [mp.mpf(0)] * n
    for m in range(0, n - 1):
        acc = mp.mpf(0)
        for j in range(m + 1):
            acc += P[j] * (m - j + 1) * y[m - j + 1] + Q[j] * y[m - j]
        y[m + 2] = acc / ((m + 2) * (m + 1))
    Y = sum(y[i] * h ** i for i in range(n))
    Yp = sum(i * y[i] * h ** (i - 1) for i in range(1, n))
    return Y, Yp


def integrate(K, w, rho0, Y0, Yp0, rho1, n=60, frac=mp.mpf("0.25"), hmax=mp.mpf(3), record=()):
    """Taylor-step integration along the real axis from rho0 to rho1 (same sign; the equation is singular at 0 and +-i, so h = frac * distance to the nearest of them)"""
    rho = mp.mpf(rho0)
    Y, Yp = Y0, Yp0
    d = 1 if rho1 > rho0 else -1
    out = {}
    recs = sorted(record, reverse=(d < 0))
    ri = 0
    while (rho1 - rho) * d > 0:
        dist = min(abs(rho), mp.sqrt(1 + rho ** 2))
        h = min(frac * dist, hmax, abs(rho1 - rho))
        while ri < len(recs) and (recs[ri] - rho) * d > 0 and (recs[ri] - rho) * d < h:
            h = abs(recs[ri] - rho)
        Y, Yp = step(K, w, rho, Y, Yp, d * h, n=n)
        rho = rho + d * h
        if ri < len(recs) and abs(rho - recs[ri]) < mp.mpf(10) ** (-mp.mp.dps + 5):
            out[recs[ri]] = (Y, Yp)
            ri += 1
    return rho, Y, Yp, out


# ----------------------------------------------------------------------------------------------------------------------
# the closed form  Y_ref = R e^{kappa x},  e^{kappa x} = |rho|^{1/2} exp(H),  H = int_0^rho (K/A - 1/(2 rho')) d rho'
# ----------------------------------------------------------------------------------------------------------------------
_HS = {}


def h_series(K, N):
    """Taylor coefficients of h = K/A - 1/(2 rho) = -(sum_{j >= 1} a_j rho^{j-1}) / (2 a), a = A/rho"""
    A = bg_series0(K, N + 4)
    a = A[1:N + 3]
    num = [-a[j + 1] for j in range(N + 1)]
    hs = s_mul(num, s_inv(a, N + 1), N + 1)
    return [x / 2 for x in hs]


def H_of(K, rho, radius=mp.mpf("0.25")):
    key = (mp.mp.dps, K)
    if key not in _HS:
        _HS[key] = h_series(K, 100)
    hs = _HS[key]
    Hs = lambda r: sum(hk * r ** (k + 1) / (k + 1) for k, hk in enumerate(hs))
    rho = mp.mpf(rho)
    if abs(rho) <= radius:
        return Hs(rho)
    sgn = 1 if rho > 0 else -1
    r0 = sgn * radius
    hf = lambda r: K / Aof(K, r) - 1 / (2 * r)
    pts = [r0]
    r = abs(r0)
    while r < abs(rho):
        r = min(r * 2, abs(rho))
        pts.append(sgn * r)
    return Hs(r0) + mp.quad(hf, pts)


def Yref_closed(K, rho):
    """Y_ref = R e^{kappa x}. Mutations: cf_R2 (R^2 -> 1 + 1.01 rho^2), cf_k (e^{1.01 kappa x}), cf_xsign (e^{-kappa x})"""
    rho = mp.mpf(rho)
    R2 = 1 + rho * rho
    kk = 1
    if MUT == "cf_R2":
        R2 = 1 + mp.mpf("1.01") * rho * rho
    if MUT == "cf_k":
        kk = mp.mpf("1.01")
    if MUT == "cf_xsign":
        kk = -1
    return mp.sqrt(R2) * abs(rho) ** (mp.mpf(kk) / 2) * mp.exp(kk * H_of(K, rho))


# ----------------------------------------------------------------------------------------------------------------------
# connection constant beta* (Y_inf = Y_nat + beta* Y_b decays at +infinity), two ways
# ----------------------------------------------------------------------------------------------------------------------
def beta_shoot(K, rho_end=100, N=60, n=70):
    """beta* = -lim Y_nat / Y_b at large rho; Y_nat: g1 = -1/(4K) (the continuation in W of the Frobenius series), Y_b = rho^{3/2}(1 + ...); Taylor integration from rho = 0.1.
    Mutation beta_g1: Y_a started with g1 = 0."""
    g1 = mp.mpf(0) if MUT == "beta_g1" else -1 / (4 * K)
    ga = frobenius_s(K, 1, mp.mpf(1) / 2, g1, N)
    gb = frobenius_s(K, 1, mp.mpf(3) / 2, None, N)
    r0 = mp.mpf("0.1")
    Ya = Y_frob(ga, 0.5, r0)
    Yb = Y_frob(gb, 1.5, r0)
    re_ = mp.mpf(rho_end)
    _, A_, Ap_, _ = integrate(K, 1, r0, Ya[0], Ya[1], re_, n=n)
    _, B_, Bp_, _ = integrate(K, 1, r0, Yb[0], Yb[1], re_, n=n)
    return -A_ / B_, -Ap_ / Bp_


def I0_quad(K, log=None):
    """I0 = int_0^inf 2K rho e^{-2H} / (A R^2) d rho  (= int 2K rho^2 e^{-2 kappa x}/(A R^2) d rho, e^{kappa x} = rho^{1/2} e^H); composite panels, H by quadrature"""
    hs = h_series(K, 100)
    Hs = lambda r: sum(hk * r ** (k + 1) / (k + 1) for k, hk in enumerate(hs))
    a_s = bg_series0(K, 100)[1:99]
    a_of = lambda r: sum(c * r ** j for j, c in enumerate(a_s))
    hf = lambda r: K / Aof(K, r) - 1 / (2 * r)
    q = mp.mpf("0.25")
    total = mp.quad(lambda r: 2 * K * mp.exp(-2 * Hs(r)) / (a_of(r) * (1 + r * r)), [0, q])
    Hleft = Hs(q)
    L = q
    while L < 110:
        Rr = 2 * L

        def f(r, L=L, Hleft=Hleft):
            H = Hleft + mp.quad(hf, [L, r])
            return 2 * K * r * mp.exp(-2 * H) / (Aof(K, r) * (1 + r * r))
        total += mp.quad(f, [L, Rr])
        Hleft = Hleft + mp.quad(hf, [L, Rr])
        L = Rr
    return total


# ----------------------------------------------------------------------------------------------------------------------
# symbolic checks S1 - S4
# ----------------------------------------------------------------------------------------------------------------------
_JC = {}


def J_entries():
    """First-order monodromy J = oint Phi0^{-1} B1 Phi0 d rho about (c, W) = (-1, kappa), exactly (K = 2/pi), basis Y_a = rho^{1/2}(1 - rho/(4K) + ...) (g1 = -1/(4K)), Y_b = rho^{3/2}(1 + ...):
    Y_j -> sum_i Y_i J_ij under rho -> rho e^{2 pi i}; B1 = d_eps B, c = -1 + eps, W fixed. Residues at rho = 0 by series. Mutations: j_sign (B1 -> -B1), j_g1 (Y_a with g1 = 0)."""
    key = MUT if MUT in ("j_sign", "j_g1") else None
    if key in _JC:
        return _JC[key]
    import sympy as sp
    NORD = 8
    rho = sp.Symbol("rho")
    K = 2 / sp.pi
    NT = NORD + 8
    T = sum((-1) ** n * rho ** (2 * n + 1) / sp.Integer(2 * n + 1) for n in range(NT))
    R2 = 1 + rho ** 2
    tr = lambda e, n=NT: sp.expand(sum(sp.expand(e).coeff(rho, k) * rho ** k for k in range(n)))
    A = tr(-rho ** 2 + K * (rho + R2 * T))
    A1 = tr(R2 - K * (rho + R2 * T))                      # d_c A at c = -1 (K(c) = -2c/pi, K' = -2/pi = -K)
    aA = sp.expand(sp.cancel(A / rho))
    W = K
    ser = lambda e, n: sp.expand(sp.series(e, rho, 0, n).removeO())
    ell = ser(A1 / (rho * aA), NORD)                      # d_eps ln A = A1/A
    ellp = sp.diff(ell, rho)
    iR2 = ser(1 / R2, NORD + 4)
    P1 = -ellp                                            # d_eps P = -(A1/A)'
    Q1 = ser(-(1 / rho) * iR2 * ellp + 2 * iR2 * A1 / (rho ** 2 * aA ** 2) - 2 * W ** 2 * A1 / (rho ** 3 * aA ** 3), NORD)   # d_eps Q at fixed W
    if MUT == "j_sign":
        P1, Q1 = -P1, -Q1
    C2 = tr(aA * aA * R2 ** 2)
    C1 = tr(-aA * R2 ** 2 * (aA - rho * sp.diff(aA, rho)))
    C0 = tr(-(rho ** 2 * aA * aA - aA * (aA + rho * sp.diff(aA, rho)) * R2 - 2 * rho * aA * R2 + W ** 2 * R2 ** 2))
    cf = lambda P, j: P.coeff(rho, j)

    def frob(s, g1, N):
        I = lambda sig: cf(C2, 0) * sig * (sig - 1) + cf(C1, 0) * sig + cf(C0, 0)
        g = [sp.Integer(1)]
        start = 1
        if s == sp.Rational(1, 2):
            g.append(g1)
            start = 2
        for n in range(start, N + 1):
            tot = sum(g[k] * (cf(C2, n - k) * (s + k) * (s + k - 1) + cf(C1, n - k) * (s + k) + cf(C0, n - k)) for k in range(n))
            g.append(sp.simplify(-tot / I(s + n)))
        return g
    sa, sb = sp.Rational(1, 2), sp.Rational(3, 2)
    g1n = sp.Integer(0) if MUT == "j_g1" else -1 / (4 * K)
    ga = frob(sa, g1n, NORD)
    gb = frob(sb, None, NORD)
    ya = sum(ga[k] * rho ** k for k in range(NORD + 1))
    yb = sum(gb[k] * rho ** k for k in range(NORD + 1))
    dY = lambda s, y: s * y + rho * sp.diff(y, rho)
    Wab = sp.expand(ya * dY(sb, yb) - dY(sa, ya) * yb)     # Wronskian = rho * Wab
    invW = ser(1 / Wab, NORD)
    out = {}
    for i in ("a", "b"):
        for j in ("a", "b"):
            yj, sj = (ya, sa) if j == "a" else (yb, sb)
            yo, so = (yb, sb) if i == "a" else (ya, sa)
            sign = -1 if i == "a" else 1
            inner = sp.expand(yo * (rho * Q1 * yj + P1 * dY(sj, yj)) * invW)
            full = sp.expand(inner * rho ** int(so + sj - 2))
            out[i + j] = sp.simplify(sign * 2 * sp.pi * sp.I * full.coeff(rho, -1))
    _JC[key] = out
    return out


def kdot_symbolic():
    """d kappa/dc at c = -1 by implicit differentiation: A(rho_h(c); c) = 0, kappa = A'(rho_h; c)/2"""
    import sympy as sp
    rho, c = sp.symbols("rho c")
    K = -2 * c / sp.pi
    A = 1 + c * (rho ** 2 + 1) + K * (rho + (rho ** 2 + 1) * sp.atan(rho))
    Ar = sp.diff(A, rho)
    at = {rho: 0, c: -1}
    rhoh_dot = -sp.diff(A, c).subs(at) / Ar.subs(at)
    kd = (sp.diff(Ar, c).subs(at) + sp.diff(Ar, rho).subs(at) * rhoh_dot) / 2
    return sp.simplify(kd), sp.simplify(rhoh_dot)


def derive_S1():
    import sympy as sp
    rho, c, W = sp.symbols("rho c W", real=True)
    K = -2 * c / sp.pi
    R = sp.sqrt(1 + rho ** 2)
    A = 1 + c * R ** 2 + K * (rho + R ** 2 * sp.atan(rho))
    Ap = sp.diff(A, rho)
    P = 2 / rho - Ap / A
    Q = 1 / R ** 4 - Ap / (A * rho * R ** 2) - 2 / (A * R ** 2) + W ** 2 / A ** 2
    E = sp.Function("E")(rho)                                   # e^{W x}, E' = W E / A
    pw = sp.Rational(101, 100) if MUT == "sym_pw" else 1
    Y = R ** pw * E
    d1 = lambda f: sp.diff(f, rho).subs(sp.Derivative(E, rho), W * E / A)
    Yp = d1(Y)
    Ypp = d1(Yp)
    res = sp.simplify((Ypp - P * Yp - Q * Y) / E)
    at_K = sp.simplify(res.subs(W, K))
    off = sp.simplify(res.subs({W: K * 2, c: -1}).subs(rho, sp.Rational(1, 2)).evalf(30))
    check("S1", "Y = R^{pw} e^{W x} (E' = W E/A) inserted in #212's Y equation for general c: the residual divided by e^{Wx} is 2 (pi W + 2c) / (rho R (2c rho^2 atan rho - pi c rho^2 + 2c rho + 2c atan rho - pi c - pi)), zero exactly at W = -2c/pi = K",
          at_K == 0 and off != 0, f"residual at W = K: {at_K}; at W = 2K, c = -1, rho = 1/2: {off}")


def derive_S2():
    import sympy as sp
    rho, b = sp.symbols("rho b", positive=True)
    R = sp.sqrt(rho ** 2 + b ** 2)
    phi = sp.sqrt(2) * sp.atan(rho / b)
    gfac = rho * R / (2 if MUT == "sym_gauge" else sp.sqrt(2))
    Y = sp.simplify((sp.diff(R, b) - gfac * sp.diff(phi, b)).subs(b, 1))
    check("S2", "c = 0 (Ellis wormhole, K = 0, W = 0): Y_K = R is the variation of the throat radius, Y = d_b R - (rho R/sqrt2) d_b phi at b = 1 with phi = sqrt2 arctan(rho/b)",
          sp.simplify(Y - sp.sqrt(rho ** 2 + 1)) == 0, f"Y = {Y}")


def derive_S3_S4():
    import sympy as sp
    J = J_entries()
    kd, rd = kdot_symbolic()
    kap = 2 / sp.pi
    Kdot = -2 / sp.pi
    exp_aa = sp.simplify(sp.I * sp.pi * kd / kap)
    ok3 = (sp.simplify(J["aa"] - exp_aa) == 0 and sp.simplify(J["bb"] + J["aa"]) == 0 and sp.simplify(J["ba"]) == 0 and sp.simplify(J["ab"] + sp.I * sp.pi ** 2) == 0)
    check("S3", "first-order monodromy about (c, W) = (-1, kappa) from the residues at rho = 0 (series to order 8, exact in pi): J_aa = i pi (d kappa/dc)/kappa = -J_bb, J_ba = 0, J_ab = -i pi^2; the prediction J_aa = i pi kdot/kappa, J_ba = 0 follows from det M = 1 and the exponents 1 -+ w/2",
          ok3, "J = [[%s, %s], [%s, %s]]; kdot = %s, d rho_h/dc = %s" % (J["aa"], J["ab"], J["ba"], J["bb"], kd, rd))
    # ingoing line at W = K(c): (v_b / v_a) = 2 pi i (Kdot - kdot) / (J_ab kappa); Y_+ = Y_nat + (g1(Y_+) - g1(nat)) Y_b, g1(Y_+) = 1/(4K)
    if sp.simplify(J["ab"]) != 0:
        ratio = sp.simplify(2 * sp.pi * sp.I * (Kdot - kd) / (J["ab"] * kap))
    else:
        ratio = sp.nan
    target = sp.simplify(1 / (4 * kap) + 1 / (4 * kap))
    check("S4", "the ingoing line at W = K(c), predicted by J, is Y_ref = R e^{kappa x} itself: v_b/v_a = 2 pi i (Kdot - kappadot)/(J_ab kappa) equals g1(Y_ref) - g1(Y_nat) = 1/(2K) = pi/4 (an exact solution that J did not use)",
          sp.simplify(ratio - target) == 0, f"predicted {ratio}, Y_ref {target}")


# ----------------------------------------------------------------------------------------------------------------------
# numerical checks N0 - N7
# ----------------------------------------------------------------------------------------------------------------------
def N0():
    import numpy as np
    import stability as st
    fam = st.Fam(-1.0)
    mp.mp.dps = 30
    K = 2 / mp.pi
    worst = 0.0
    for r in (-3.0, -0.7, -0.2, 0.3, 1.1, 4.0):
        A = float(fam.A(r))
        Ap = float(fam.Ap(r))
        Cx_s, Dx_s = st.coef_x(fam, r, A, A / r, Ap)
        P, Q = PQ(K, mp.mpf(r), 0)                              # Q at W = 0 is Q_0
        Am = Aof(K, mp.mpf(r))
        Apm = Apof(K, mp.mpf(r))
        Cx = Am * P + Apm
        Dx = Am ** 2 * Q
        worst = max(worst, abs(float(Cx) - Cx_s) / (1e-12 + abs(Cx_s)), abs(float(Dx) - Dx_s) / (1e-12 + abs(Dx_s)))
    check("N0", "the Y equation as coded here (P, Q; Y_xx = (A P + A') Y_x + A^2 (Q - W^2/A^2) Y + W^2 Y) equals #212's coef_x (stability.py) at 6 points of rho in [-3, 4]", worst < 1e-12, f"max relative difference {worst:.1e}")


def N1():
    mp.mp.dps = 60
    K = 2 / mp.pi
    Nser = 90
    g = frobenius_s(K, 1, mp.mpf(1) / 2, 1 / (4 * K), Nser)
    worst = mp.mpf(0)
    for s in (1, -1):
        for x in ("0.01", "0.05", "0.1"):
            r = s * mp.mpf(x)
            ys = Y_frob(g, 0.5, r)[0]
            worst = max(worst, abs(Yref_closed(K, r) - ys) / abs(ys))
    check("N1", "Y_ref = R e^{kappa x} against the Frobenius series rho^{1/2} sum g_k rho^k (g1 = 1/(4K) = pi/8, g_k, k >= 2, from the recursion, 90 terms, 60 digits) at rho = +-0.01, +-0.05, +-0.1",
          worst < mp.mpf(10) ** -30, f"max relative difference {mp.nstr(worst, 3)} (bar 1e-30); g2 = {mp.nstr(g[2], 20)} = 5/12 + 3 pi^2/128: {mp.nstr(g[2] - (mp.mpf(5) / 12 + 3 * mp.pi ** 2 / 128), 3)}")


def N2():
    mp.mp.dps = 45
    K = 2 / mp.pi
    g = frobenius_s(K, 1, mp.mpf(1) / 2, 1 / (4 * K), 60)
    pts = [mp.mpf(x) for x in ("0.3", "1", "2", "5", "10", "20", "40")]
    worst = mp.mpf(0)
    for sgn in (1, -1):
        r0 = sgn * mp.mpf("0.1")
        Y0, Yp0 = Y_frob(g, 0.5, r0)
        _, _, _, out = integrate(K, 1, r0, Y0, Yp0, sgn * mp.mpf(40), n=70, record=[sgn * p for p in pts])
        for p in pts:
            r = sgn * p
            worst = max(worst, abs(Yref_closed(K, r) - out[r][0]) / abs(out[r][0]))
    check("N2", "Y_ref against the numerical integration of the master equation at W = kappa (Taylor steps of 70 terms, 45 digits, from the series at rho = +-0.1) at rho = +-0.3, 1, 2, 5, 10, 20, 40: the exterior and the interior",
          worst < mp.mpf(10) ** -20, f"max relative difference {mp.nstr(worst, 3)} (bar 1e-20)")


def N3():
    """formal series at rho = +infinity: Y = e^{eps kappa x} Phi, Phi = rho sum c_n rho^{-n}. The recursion amplifies rounding like the divergent series it also
    produces (a factor n!/(2K)^n), so it runs at 300 digits."""
    mp.mp.dps = 300
    K = 2 / mp.pi
    n = 120
    M = n + 8
    S = [mp.mpf(0)] * (M + 2)
    for m in range(0, (M + 2) // 2):
        S[2 * m] = mp.mpf((-1) ** m) / (2 * m + 1)
    P1 = [S[i] + (S[i - 2] if i >= 2 else 0) for i in range(M + 2)]
    A = [mp.mpf(0)] * M
    A[0] = mp.mpf(1)
    for i in range(M):
        A[i] += -K * P1[i + 1]                                   # A = 1 - K (P1(s) - 1)/s, s = 1/rho
    dA = s_der(A)
    Ap = [mp.mpf(0), mp.mpf(0)] + [-x for x in dA][:M - 2]
    iA = s_inv(A, M)
    AA = s_mul(Ap, iA, M)
    iR2 = s_mul([mp.mpf(0), mp.mpf(0), mp.mpf(1)], s_inv([mp.mpf(1), mp.mpf(0), mp.mpf(1)], M), M)
    iR4 = s_mul(iR2, iR2, M)
    sser = [mp.mpf(0), mp.mpf(1)] + [mp.mpf(0)] * (M - 2)
    res = {}
    sg = -1 if MUT == "inf_P" else 1
    for eps in (1, -1):
        Pe = [sg * eps * 2 * K * iA[i] + AA[i] - (2 if i == 1 else 0) for i in range(M)]
        P = [(2 if i == 1 else 0) - AA[i] for i in range(M)]
        dinvA = [-x for x in s_mul(Ap, s_mul(iA, iA, M), M)]
        PA = s_mul(P, iA, M)
        t1 = s_mul(s_mul(AA, sser, M), iR2, M)
        t2 = s_scale(s_mul(iA, iR2, M), 2, M)
        Q0 = [iR4[i] - t1[i] - t2[i] for i in range(M)]
        Qe = [eps * K * dinvA[i] - eps * K * PA[i] - Q0[i] for i in range(M)]
        phi = {-1: mp.mpf(1)}
        for m in range(1, n + 1):
            tot = mp.mpf(0)
            if (m - 2) in phi:
                tot += phi[m - 2] * (m - 2) * (m - 1)
            for j in range(1, m + 2):
                k = m - 1 - j
                if k in phi:
                    tot -= Pe[j] * k * phi[k]
            for j in range(2, m + 2):
                k = m - j
                if k in phi:
                    tot += Qe[j] * phi[k]
            phi[m - 1] = tot / (2 * eps * K * m)
        res[eps] = phi
    # eps = +1: coefficients of rho (1 + rho^{-2})^{1/2} = R
    worst = mp.mpf(0)
    for k in range(-1, n - 1):
        j = k + 1                                                 # phi_k multiplies rho^{-k}; R = sum binom(1/2, m) rho^{1 - 2m}
        ref = mp.binomial(mp.mpf(1) / 2, j // 2) if j % 2 == 0 else mp.mpf(0)
        worst = max(worst, abs(res[1][k] - ref))
    k0 = 100
    rat = abs(res[-1][k0 + 1] / res[-1][k0]) / k0
    pred = 1 / (2 * K)
    ok = worst < mp.mpf(10) ** -40 and abs(rat / pred - 1) < 0.05
    check("N3", "formal series at rho = +infinity, Y = e^{+-kappa x} rho sum c_n rho^{-n}: the e^{+kappa x} series is that of R = rho sqrt(1 + rho^{-2}) (120 terms, difference below 1e-40); the e^{-kappa x} series diverges, |c_{n+1}/c_n|/n at n = 100 within 5% of 1/(2K) = pi/4",
          ok, f"max difference from R's series {mp.nstr(worst, 3)}; e^(-kappa x): |c_101/c_100|/100 = {mp.nstr(rat, 6)} against pi/4 = {mp.nstr(pred, 6)}")


def N4():
    mp.mp.dps = 50
    K = 2 / mp.pi
    t0 = time.time()
    b1, b2 = beta_shoot(K)
    I0 = get_I0()
    bq = mp.pi / 4 - 1 / I0
    d = max(abs(b1 - bq), abs(b2 - bq)) / abs(bq)
    check("N4", "beta* (Y_inf = Y_nat + beta* Y_b decays at +infinity), two ways: (i) -Y_nat/Y_b at rho = 100 from Taylor integration of the master equation (Y_nat, Y_b started from their series at rho = 0.1; from Y and from Y'); (ii) beta* = pi/4 - 1/I0, I0 = int_0^inf 2K rho^2 e^{-2 kappa x}/(A R^2) d rho by quadrature of the closed form",
          d < mp.mpf(10) ** -30, f"beta* = {mp.nstr(bq, 30)}; I0 = {mp.nstr(I0, 30)}; 1/I0 = {mp.nstr(1 / I0, 30)}; relative difference {mp.nstr(d, 3)} (bar 1e-30); {time.time() - t0:.0f} s")
    return I0, bq


def rates(res, cs):
    import numpy as np
    from scipy.optimize import brentq
    import stability as st
    out = []
    for c in cs:
        fam = st.Fam(c, res)
        kap = fam.A1 / 2
        f = lambda W: st.mismatch(fam, W)
        lo, hi = 0.97 * kap, 1.03 * kap
        if np.sign(f(lo)) == np.sign(f(hi)):
            lo, hi = 0.9 * kap, 1.1 * kap
        W = brentq(f, lo, hi, xtol=1e-15, rtol=1e-14)
        out.append((c, float(kap), float(W)))
    return out


FIT_C = (-1.01, -1.001, -1.0001, -0.9999, -0.999, -0.99)


def fit_eps(eps, y):
    """least squares of y/eps on 1, eps (second order in eps) and on 1, eps, eps^2 (third order); returns the coefficient vectors. Mutation fit_eps: eps -> -eps"""
    import numpy as np
    if MUT == "fit_eps":
        eps = -eps
    M2 = np.vstack([np.ones_like(eps), eps]).T
    c2 = np.linalg.lstsq(M2, y / eps, rcond=None)[0]
    M3 = np.vstack([np.ones_like(eps), eps, eps ** 2]).T
    c3 = np.linalg.lstsq(M3, y / eps, rcond=None)[0]
    return c2, c3


def fit_slope(res):
    import numpy as np
    r = rates(res, FIT_C)
    eps = np.array([c + 1.0 for c, _, _ in r])
    y = np.array([W / k - 1.0 for _, k, W in r])
    c2, c3 = fit_eps(eps, y)
    W = {round(c, 4): w for c, _, w in r}
    cd = (W[-0.9999] - W[-1.0001]) / 2e-4
    return c2, c3, cd, r


_I0 = {}


def get_I0():
    if "v" not in _I0:
        mp.mp.dps = 50
        _I0["v"] = I0_quad(2 / mp.pi)
    return _I0["v"]


def N5_N6():
    mp.mp.dps = 50
    K = 2 / mp.pi
    I0 = get_I0()
    beta = mp.pi / 4 - 1 / I0
    J = J_entries()
    Jab = complex(J["ab"].evalf(30))
    Jab_mp = mp.mpc(Jab.real, Jab.imag)
    if MUT == "s_g1":
        beta = beta - mp.pi / 8                                   # beta* in the g1 = 0 basis
    s = (Jab_mp * beta / (2j * mp.pi)).real
    kdot = mp.pi / 4 - 2 / mp.pi
    if MUT == "s_kappa":
        s = s + kdot / K                                          # d W/dc / kappa instead of d (W/kappa)/dc
    fits = {}
    for res in ("R2", "R1"):
        fits[res] = fit_slope(res)
    c2, c3, cd, r = fits["R2"]
    rel2 = abs(float(s) - c2[0]) / abs(c2[0])
    rel3 = abs(float(s) - c3[0]) / abs(c3[0])
    check("N5", "d(W/kappa)/dc at c = -1 derived by first-order monodromy, s = J_ab beta*/(2 pi i) = -(pi/2) beta*, against the fit of W/kappa(c) - 1 = s1 (c+1) + s2 (c+1)^2 to #212's rates (stability.py, R2) at c = -1 +- 1e-2, 1e-3, 1e-4, to 1e-3",
          rel2 < 1e-3,
          f"s = {mp.nstr(s, 15)}; fit s1 = {c2[0]:.8f} (relative difference {rel2:.1e}), s2 = {c2[1]:.5f}; third-order fit s1 = {c3[0]:.8f} ({rel3:.1e}); R1 fit s1 = {fits['R1'][0][0]:.8f}")
    dwdc = 1 / I0 - 2 / mp.pi
    if MUT == "dwdc_K":
        dwdc = 1 / I0
    rel = abs(float(dwdc) - cd) / abs(cd)
    check("N6", "dW/dc at c = -1, 1/I0 - 2/pi, against the central difference of #212's W (R2) at c = -0.9999 and -1.0001, to 1e-4",
          rel < 1e-4, f"1/I0 - 2/pi = {mp.nstr(dwdc, 15)}; central difference {cd:.10f}; relative difference {rel:.1e}")
    return s, dwdc


PSLQ_B1 = ("1", "pi", "pi^2", "pi^3", "1/pi", "1/pi^2")
PSLQ_B2 = ("1", "pi", "pi^2", "1/pi", "ln2", "ln pi", "G", "gamma")
PSLQ_TOL = mp.mpf(10) ** -33
PSLQ_MAXC = 1000


def pslq_basis(names):
    val = {"1": mp.mpf(1), "pi": mp.pi, "pi^2": mp.pi ** 2, "pi^3": mp.pi ** 3, "1/pi": 1 / mp.pi, "1/pi^2": 1 / mp.pi ** 2, "ln2": mp.log(2), "ln pi": mp.log(mp.pi),
           "G": mp.catalan, "gamma": mp.euler}
    return [val[n] for n in names]


def pslq_try(t, names):
    return mp.pslq([t] + pslq_basis(names), tol=PSLQ_TOL, maxcoeff=PSLQ_MAXC, maxsteps=200000)


def N7():
    mp.mp.dps = 50
    planted = mp.mpf(3) / 2 - 5 * mp.pi / 3 + 2 * mp.pi ** 2
    if MUT == "pslq_plant":
        planted += mp.mpf(10) ** -20
    rel = pslq_try(planted, PSLQ_B1)
    ok_pl = False
    if rel is not None:
        gcd = math.gcd(*[abs(int(x)) for x in rel if x != 0])
        rn = [int(x) // gcd for x in rel]
        rn = [-x for x in rn] if rn[0] < 0 else rn
        ok_pl = rn == [6, -9, 10, -12, 0, 0, 0]
    control = mp.e ** (mp.mpf(1) / 3) + mp.sqrt(2)
    rc = pslq_try(control, PSLQ_B1)
    rc2 = pslq_try(control, PSLQ_B2)
    check("N7", "PSLQ at 50 digits (tolerance 1e-33, max coefficient 1000) recovers the planted number 3/2 - 5 pi/3 + 2 pi^2 (relation 6 t - 9 + 10 pi - 12 pi^2 = 0) in the basis {1, pi, pi^2, pi^3, 1/pi, 1/pi^2}, and returns none for exp(1/3) + sqrt2 in both declared bases",
          ok_pl and rc is None and rc2 is None, f"planted: {rel}; control: {rc}, {rc2}")


def N8():
    """planted: synthetic W/kappa - 1 = s1 e + s2 e^2 + s3 e^3 with known coefficients and noise 1e-10 (the size of #212's root-finding noise) at the six c, through the same fit"""
    import numpy as np
    rs = np.random.RandomState(262)
    eps = np.array([c + 1.0 for c in FIT_C])
    s1, s2, s3 = 0.314029018887, 0.304, 0.27
    y = s1 * eps + s2 * eps ** 2 + s3 * eps ** 3 + 1e-10 * rs.uniform(-1, 1, size=eps.shape) / 0.6366
    c2, c3 = fit_eps(eps, y)
    d2, d3 = abs(c2[0] - s1) / s1, abs(c3[0] - s1) / s1
    check("N8", "planted: synthetic W/kappa - 1 = s1 e + s2 e^2 + s3 e^3 (s1 = 0.314029, s2 = 0.304, s3 = 0.27) with noise 1e-10 at the six c of N5, through the same fit: s1 recovered to 1e-3 (second and third order)",
          d2 < 1e-3 and d3 < 1e-3, f"second order {c2[0]:.8f} ({d2:.1e}), third order {c3[0]:.8f} ({d3:.1e})")


def pslq_real(I0):
    out = []
    mp.mp.dps = 50
    tau = 1 / I0
    s = -(mp.pi / 2) * (mp.pi / 4 - 1 / I0)
    for tname, t in (("1/I0 = dW/dc + 2/pi", tau), ("d(W/kappa)/dc", s)):
        for bname, names in (("B1", PSLQ_B1), ("B2", PSLQ_B2)):
            out.append((tname, bname, pslq_try(t, names)))
    return out


# ----------------------------------------------------------------------------------------------------------------------
MUTATIONS = {
    # name: checks that must fire (development runs of negative-mass-raw I262-mac/dev, before they were declared)
    "sym_pw": ["S1"], "sym_gauge": ["S2"], "j_sign": ["S3", "S4", "N5"], "j_g1": ["S3"],
    "eq_P": ["N0", "N2", "N4"], "cf_R2": ["N1", "N2"], "cf_k": ["N1", "N2"], "cf_xsign": ["N1", "N2"],
    "inf_P": ["N3"], "beta_g1": ["N4"], "s_g1": ["N5"], "s_kappa": ["N5"], "dwdc_K": ["N6"], "pslq_plant": ["N7"], "fit_eps": ["N5", "N8"],
}
DERIVE = [("S1", derive_S1), ("S2", derive_S2), ("S3", derive_S3_S4)]
NUMERIC = [("N0", N0), ("N1", N1), ("N2", N2), ("N3", N3), ("N4", N4), ("N5", N5_N6), ("N7", N7), ("N8", N8)]


def run_all(items):
    t0 = time.time()
    for cid, fn in items:
        if want(cid) or want(cid + "a") or (cid == "S3" and want("S4")) or (cid == "N5" and want("N6")):
            t1 = time.time()
            fn()
            print(f"    ({cid}: {time.time() - t1:.0f} s)", flush=True)
    bad = [cid for cid, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS)} checks, {len(RESULTS) - len(bad)} pass; failed: {bad}; mutation {MUT}; {time.time() - t0:.0f} s")
    if MUT is not None and MUT in MUTATIONS:
        want_ = MUTATIONS[MUT]
        print(f"mutation {MUT}: fired {bad}, declared {want_}: {'equal' if sorted(bad) == sorted(want_) else 'DIFFERENT'}")
    return bad


def table():
    mp.mp.dps = 50
    K = 2 / mp.pi
    g = frobenius_s(K, 1, mp.mpf(1) / 2, -1 / (4 * K), 8)
    g0 = frobenius_s(K, 1, mp.mpf(1) / 2, 0, 4)
    gp = frobenius_s(K, 1, mp.mpf(1) / 2, 1 / (4 * K), 4)
    print("Frobenius coefficients of Y at W = kappa, exponent 1/2 (g1 free; K = 2/pi):")
    print("  g1 = -1/(4K) = -pi/8 (Y_nat):  g2 =", mp.nstr(g[2], 25), "= 5/12 - pi^2/128")
    print("  g1 = 0 (#248 D7):             g2 =", mp.nstr(g0[2], 25), "= 5/12 + pi^2/128")
    print("  g1 = 1/(4K) = pi/8 (Y_ref):   g2 =", mp.nstr(gp[2], 25), "= 5/12 + 3 pi^2/128")
    I0 = get_I0()
    beta = mp.pi / 4 - 1 / I0
    b1, b2 = beta_shoot(K)
    s = -(mp.pi / 2) * beta
    kdot = mp.pi / 4 - 2 / mp.pi
    print("\nconstants (50 digits):")
    print("  I0                    =", mp.nstr(I0, 40))
    print("  1/I0 = d(W - K)/dc    =", mp.nstr(1 / I0, 40))
    print("  beta* (quadrature)    =", mp.nstr(beta, 40))
    print("  beta* (shooting)      =", mp.nstr(b2, 40))
    print("  d(W/kappa)/dc = -(pi/2) beta*  =", mp.nstr(s, 40))
    print("  dW/dc = 1/I0 - 2/pi            =", mp.nstr(1 / I0 - 2 / mp.pi, 40))
    print("  dK/dc = -2/pi, dkappa/dc = pi/4 - 2/pi =", mp.nstr(kdot, 20), "; (dW/dc - dkappa/dc)/kappa =", mp.nstr((1 / I0 - 2 / mp.pi - kdot) / K, 30))
    print("  Y_inf = Y_ref I(rho)/I0,  I(rho) = int_rho^inf 2K r^2 e^{-2 kappa x}/(A R^2) dr;  Y_inf = Y_ref - Y_b/I0 = Y_nat + beta* Y_b")
    print("  e^{2 kappa x} = rho exp(2H) at rho = 1 (the normalisation here; #248's chart has Pi(rho = 1) = 1):", mp.nstr(mp.exp(2 * H_of(K, mp.mpf(1))), 25))
    print("\nrates and fits (#212's code):")
    for res in ("R1", "R2"):
        c2, c3, cd, r = fit_slope(res)
        print(f"  {res}: s1 = {c2[0]:.9f}, s2 = {c2[1]:.6f} (second order); s1 = {c3[0]:.9f}, s2 = {c3[1]:.6f}, s3 = {c3[2]:.4f} (third order); central difference dW/dc = {cd:.10f}")
        for c, k, W in r:
            print(f"      c = {c:8.4f}  kappa = {k:.12f}  W = {W:.13f}  W/kappa - 1 = {W / k - 1:.6e}")
    print("\nW against K = -2c/pi (3m) at #212's members (README.md Sec. 4 table; W is the rate of the one growing mode):")
    below, above = [], []
    sec, seen = False, False
    for line in open(os.path.join(HERE, "README.md")):
        if line.startswith("## "):
            if sec:
                seen = True
            sec = line.startswith("## 4.") and not seen
        cells = [x.strip().replace("\u2212", "-") for x in line.split("|")]
        if sec and len(cells) > 6 and re.fullmatch(r"-?\d+(\.\d+)?", cells[1]) and re.fullmatch(r"\d+\.\d+", cells[5] or ""):
            c, W = float(cells[1]), float(cells[5])
            if abs(c + 1) < 1e-12:
                continue
            (below if W < -2 * c / math.pi else above).append((c, W, -2 * c / math.pi))
    print("  W < K at c =", [c for c, _, _ in below])
    print("  W > K at c =", [c for c, _, _ in above])
    print("\nPSLQ (50 digits, tolerance 1e-33, max coefficient 1000), basis declared before the run:")
    for tname, bname, rel in pslq_real(I0):
        print(f"  target {tname:22s} basis {bname}: {rel}")


def report_mutations(raw, prefix):
    print("== Mutations: fired sets in %s<name>.log against the declared ones ==" % prefix)
    n_eq = 0
    for name, declared in MUTATIONS.items():
        fn = os.path.join(raw, f"{prefix}{name}.log")
        if not os.path.exists(fn):
            print(f"  {name:11s} NO LOG")
            continue
        m = re.findall(r"failed: (\[.*?\]); mutation", open(fn).read())
        fired = sorted(set(x for grp in m for x in ast.literal_eval(grp))) if m else None
        eq = fired is not None and sorted(fired) == sorted(declared)
        n_eq += eq
        print(f"  {name:11s} fired {fired}  declared {declared}  {'equal' if eq else 'DIFFERENT'}")
    print(f"{n_eq} of {len(MUTATIONS)} equal")
    return n_eq == len(MUTATIONS)


def main():
    sub = next((a for a in sys.argv[1:] if not a.startswith("--")), "check")
    if sub == "derive":
        bad = run_all(DERIVE)
    elif sub == "check":
        bad = run_all(NUMERIC)
    elif sub == "table":
        table()
        return
    elif sub == "mutations":
        args = [a for a in sys.argv[2:] if not a.startswith("--")]
        sys.exit(0 if report_mutations(args[0], args[1] if len(args) > 1 else "") else 1)
    else:
        sys.exit(f"unknown command {sub}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
