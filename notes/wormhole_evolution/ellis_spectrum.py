"""Linear spherically symmetric spectrum of the symmetric Ellis wormhole (D = 5 and D = 4) beyond its unstable mode (#272).

The equations are those of notes/wormhole_largeD (derive.py, master.py): gauge delta phi = 0, h = diag(-e^{2nu}, e^{2lam}),
r = r0 + eps Y, with the background r0'^2 = 1 - r0^{-2(n-1)}, n = D - 2, throat radius 1. This file re-derives the time-dependent
linear system, tests the zero modes as exact linear solutions, finds the poles of the resolvent of the master equation E2 by two
eigen-solvers, and evolves E2 in time.

    python3 ellis_spectrum.py run [--only=planted,zero,spectrum,time] [--mutate=NAME] [--crate=PATH]

Every check prints PASS or FAIL with its numbers; the exit status is 1 when any check fails. --mutate=NAME runs a control (see MUTATIONS).
"""
import json
import os
import subprocess
import sys
import time as _time
import warnings

import numpy as np
import scipy.linalg as sl
import scipy.sparse as sp_
import scipy.sparse.linalg as sla
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
MUTATIONS = ("wrong_p", "pt_V", "ray_flip", "no_projection", "wrong_secular", "ecs_noscale", "pole_off")
MUTATE = None
ONLY = None
CRATE = os.path.join(HERE, "..", "wormhole_largeD", "w2_crate", "target", "release", "w2_crate")
RESULTS = []


def check(name, ok, text):
    RESULTS.append((name, ok))
    print(("PASS    " if ok else "FAIL    ") + name + ": " + text, flush=True)


# ------------------------------------------------------------------ power series for the throat (Frobenius at l = 0)
def smul(a, b, K):
    return np.convolve(a, b)[:K]


def spow(a, p, K):
    b = np.zeros(K, dtype=complex)
    b[0] = a[0] ** p
    for k in range(1, K):
        b[k] = sum((p * j - (k - j)) * a[j] * b[k - j] for j in range(1, k + 1)) / (k * a[0])
    return b


_CACHE = {}


def throat_series(n, K=44):
    if n in _CACHE:
        return _CACHE[n]
    r = np.zeros(K, dtype=complex)
    r[0] = 1.0
    for k in range(0, K - 2):
        r[k + 2] = (n - 1) * spow(r, 1 - 2 * n, K)[k] / ((k + 2) * (k + 1))
    rho = np.array([(j + 2) * r[j + 2] for j in range(K - 2)] + [0, 0], dtype=complex)
    rpp = np.array([(k + 2) * (k + 1) * r[k + 2] for k in range(K - 2)] + [0, 0], dtype=complex)
    rho_inv = np.concatenate([spow(rho[:K - 2], -1, K - 2), np.zeros(2)])
    r_inv = spow(r, -1, K)
    l2 = np.zeros(K, dtype=complex)
    l2[2] = 1
    alpha = 2 * smul(rpp, rho_inv, K) + n * smul(smul(l2, rho, K), r_inv, K)
    beta = n * smul(smul(l2, smul(rho, rho, K), K), smul(r_inv, r_inv, K), K) + smul(rpp, r_inv, K)
    _CACHE[n] = (r, rho, alpha, beta)
    return _CACHE[n]


def frob(n, W2, parity, l0, K=44):
    """r, r', Y, Y' at l0: Y(0) = 1 (even) or Y ~ l^3 (odd). k(k-3) y_k = sum_j alpha_j (k-j) y_{k-j} - sum_j gamma_j y_{k-2-j}"""
    r, rho, alpha, beta = throat_series(n, K)
    gam = beta.copy()
    gam[0] -= W2
    y = np.zeros(K, dtype=complex)
    y[0 if parity == "even" else 3] = 1
    for k in range(1, K):
        if k == 3 or (parity == "odd" and k < 3):
            continue
        s = sum(alpha[j] * (k - j) * y[k - j] for j in range(2, k + 1)) - sum(gam[j] * y[k - 2 - j] for j in range(0, k - 1))
        y[k] = s / (k * (k - 3))
    pw = l0 ** np.arange(K)
    d = np.arange(1, K)
    return (np.dot(r, pw), np.dot(r[1:] * d, pw[:-1]), np.dot(y, pw), np.dot(y[1:] * d, pw[:-1]))


def rhs_e2(n, W2):
    shift = 1 if MUTATE == "wrong_p" else 0

    def f(l, s):
        r, rp, y, yp = s
        rpp = (n - 1) * r ** (1 - 2 * n)
        return [rp, rpp, yp, (2 * rpp / rp + (n + shift) * rp / r) * yp - (n * rp ** 2 / r ** 2 + rpp / r) * y + W2 * y]
    return f


def path(f, s0, start, d, length, rtol=1e-12, atol=1e-15):
    sol = solve_ivp(lambda s, y: [d * v for v in f(start + s * d, y)], (0, length), np.array(s0, dtype=complex),
                    method="DOP853", rtol=rtol, atol=atol)
    return sol.y[:, -1]


def _outgoing_e2(n, W, parity, theta, sfar, lm, l0=0.25, rtol=1e-12):
    f = rhs_e2(n, W * W)
    R, Rp, Y, Yp = frob(n, W * W, parity, l0)
    s = path(f, [R, Rp, Y, Yp], 0.0, 1.0, lm - l0, rtol)
    d = np.exp(1j * theta)
    bg = path(lambda l, y: [y[1], (n - 1) * y[0] ** (1 - 2 * n)], [s[0], s[1]], lm, d, sfar, rtol)
    lfar = lm + d * sfar
    back = path(f, [bg[0], bg[1], 1.0 + 0j, -W + n / (2 * lfar)], lm, d, -sfar, rtol)
    return s, back, lfar


def Ma_e2(n, W, parity, theta, sfar, lm=2.0, rtol=1e-12):
    """analytic in W: Wronskian of the throat solution and the outgoing solution at lm, times e^{-W(l_far - lm)} (l_far/lm)^{-n/2}"""
    s, back, lfar = _outgoing_e2(n, W, parity, theta, sfar, lm, rtol=rtol)
    return (s[3] * back[2] - s[2] * back[3]) * np.exp(-W * (lfar - lm)) * (lfar / lm) ** (-n / 2)


def logder_e2(n, W, parity, theta, sfar, lm=2.0, rtol=1e-13):
    s, back, _ = _outgoing_e2(n, W, parity, theta, sfar, lm, rtol=rtol)
    return s[3] / s[2], back[3] / back[2]


def sech2(x):
    e = np.exp(-x if x.real >= 0 else x)
    return 4 * e * e / (1 + e * e) ** 2


PT_V0, PT_AL = 2.5, 1.0


def Ma_pt(W, parity, theta, sfar, rtol=1e-12):
    """Poeschl-Teller u'' = (W^2 + V0 sech^2(al l)) u, parity at l = 0, outgoing along the ray from 0"""
    V0 = PT_V0 + (0.1 if MUTATE == "pt_V" else 0.0)
    d = np.exp(1j * theta)
    sol = solve_ivp(lambda s, y: [d * y[1], d * (W * W + V0 * sech2(PT_AL * s * d)) * y[0]], (sfar, 0.0), np.array([1.0 + 0j, -W]),
                    method="DOP853", rtol=rtol, atol=1e-30)
    u, up = sol.y[0, -1], sol.y[1, -1]
    return (up if parity == "even" else u) * np.exp(-W * d * sfar)


def ray_angle(W, margin=0.3, tmin=-2.2):
    th = float(np.clip(np.pi / 2 - margin - np.angle(W), tmin, 0.0))
    return -th if MUTATE == "ray_flip" else th


def secant(F, x0, x1, tol=1e-14, maxit=60):
    f0, f1 = F(x0), F(x1)
    for it in range(maxit):
        if f1 == f0:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1 = x1, f1, x2
        f1 = F(x1)
        if abs(x1 - x0) < tol * max(1, abs(x1)):
            return x1, it + 1
    return x1, maxit


def box_params(x1, x2, y1, y2, margin=0.30, Tsup=20.0, thmin=-2.2):
    corners = [complex(x, y) for x in (x1, x2) for y in (y1, y2)]
    th = float(np.clip(np.pi / 2 - margin - max(np.angle(c) for c in corners), thmin, 0.0))
    if MUTATE == "ray_flip":
        th = -th
    m = min((c * np.exp(1j * th)).real for c in corners)
    return th, Tsup / max(m, 1e-3), m


def winding(F, x1, x2, y1, y2, h=0.1):
    pts = []
    for a, b in [((x1, y1), (x2, y1)), ((x2, y1), (x2, y2)), ((x2, y2), (x1, y2)), ((x1, y2), (x1, y1))]:
        k = max(int(np.ceil(abs(complex(b[0] - a[0], b[1] - a[1])) / h)), 4)
        pts += [complex(a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])) for t in np.linspace(0, 1, k, endpoint=False)]
    v = np.array([F(w) for w in pts])
    ph = np.angle(np.roll(v, -1) / v)
    return ph.sum() / (2 * np.pi), np.abs(ph).max()


def scan(n, parity, X=8.0, Y=8.0, tile=2.0, h=0.12, log=None):
    """argument-principle scan of Re W in [-X-0.2, -0.2], Im W in [0, Y]; returns the zeros of M and the number of skipped cells"""
    roots, skipped = [], []

    def box(x1, x2, y1, y2):
        th, sf, m = box_params(x1, x2, y1, y2, Tsup=14.0 if x2 > -2 else 20.0)
        if m < 0.15:
            if max(x2 - x1, y2 - y1) > 0.3:
                xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
                for q in [(x1, xm, y1, ym), (xm, x2, y1, ym), (x1, xm, ym, y2), (xm, x2, ym, y2)]:
                    box(*q)
            else:
                skipped.append((x1, x2, y1, y2))
            return
        F = lambda W: Ma_e2(n, W, parity, th, sf)
        w, mx = winding(F, x1, x2, y1, y2, h)
        if abs(w) < 0.5:
            return
        if max(x2 - x1, y2 - y1) <= 0.2:
            c = complex((x1 + x2) / 2, (y1 + y2) / 2)
            r, _ = secant(F, c, c * 1.002 + 0.001j, maxit=50)
            roots.append((r, w, mx))
            return
        xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
        for q in [(x1, xm, y1, ym), (xm, x2, y1, ym), (x1, xm, ym, y2), (xm, x2, ym, y2)]:
            box(*q)
    xs = np.concatenate([np.arange(-X - 0.2, -0.2, tile), [-0.2]])
    ys = np.arange(0.0, Y + 1e-9, tile)
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            box(xs[i], xs[i + 1], ys[j], ys[j + 1])
    return roots, skipped


def polish_stable(n, parity, w0):
    """the root over ray angles and far points at lm = 2 (a root that moves is not one)"""
    out = []
    for th in (-1.8, -1.87, -2.0, -2.1):
        for T in (12.0, 20.0):
            if (w0 * np.exp(1j * th)).real < 0.1:
                continue
            F = lambda w: Ma_e2(n, w, parity, th, T / ((w * np.exp(1j * th)).real))
            try:
                r, _ = secant(F, w0, w0 * 1.001 + 0.001j, maxit=40)
                out.append(r)
            except Exception:
                pass
    return np.array(out)


# ------------------------------------------------------------------ second eigen-solver: Chebyshev collocation on a complex-scaled contour
def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    dX = x[:, None] - x[None, :]
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    return D - np.diag(D.sum(axis=1)), x


def _contour(s, s0, w, theta):
    c = np.exp(1j * theta) - 1
    x, x0 = (np.abs(s) - s0) / w, -s0 / w
    Fz = lambda x: x + np.log(np.cosh(x))
    return (np.sign(s) * (np.abs(s) + c * w / 2 * (Fz(x) - Fz(x0))), 1 + c * 0.5 * (1 + np.tanh(x)),
            np.sign(s) * c * 0.5 / np.cosh(x) ** 2 / w)


def _grid(N, S, amap=2.5):
    D, x = cheb(N)
    s = S * np.sinh(amap * x) / np.sinh(amap)
    sx = S * amap * np.cosh(amap * x) / np.sinh(amap)
    D1 = D / sx[:, None]
    return s, D1, D1 @ D1


def ecs_operator(kind, n, N=481, S=26.0, s0=1.5, w=0.3, theta=-0.5):
    """A Y = W^2 Y on the contour l = z(s): kind 'e2' (E2 of the wormhole, n = D - 2) or 'pt' (u'' - V u, V = V0 sech^2)"""
    if MUTATE == "ecs_noscale":
        theta = 0.0
    s, D1, D2 = _grid(N, S)
    z, dz, ddz = _contour(s, s0, w, theta)
    Dl = D1 / dz[:, None]
    Dll = (D2 - (ddz / dz)[:, None] * D1) / (dz ** 2)[:, None]
    if kind == "pt":
        V = np.array([PT_V0 * sech2(PT_AL * zz) for zz in z])
        return Dll - np.diag(V), s, N
    c = np.exp(1j * theta) - 1
    uniq = np.unique(np.abs(s))
    sol = solve_ivp(lambda ss, y: [(1 + c * 0.5 * (1 + np.tanh((ss - s0) / w))) * v for v in (y[1], (n - 1) * y[0] ** (1 - 2 * n))],
                    (0, uniq[-1]), np.array([1.0 + 0j, 0j]), t_eval=uniq, method="DOP853", rtol=1e-13, atol=1e-15)
    idx = np.searchsorted(uniq, np.abs(s))
    r, rp = sol.y[0][idx], sol.y[1][idx] * np.sign(s)
    rpp = (n - 1) * r ** (1 - 2 * n)
    a = 2 * rpp / rp + (n + (1 if MUTATE == "wrong_p" else 0)) * rp / r
    return Dll - a[:, None] * Dl + np.diag(n * rp ** 2 / r ** 2 + rpp / r), s, N


def ecs_eigs(kind, n, parity, **kw):
    A, s, N = ecs_operator(kind, n, **kw)
    half = np.arange(1, N)[: (N - 1) // 2]
    mirror = N - half
    B = A[np.ix_(half, half)] + (1 if parity == "even" else -1) * A[np.ix_(half, mirror)]
    return sl.eigvals(B)


def ecs_isolated(n, parity, N=481, ths=(-0.7, -0.85, -1.0), wmax=80.0):
    """eigenvalues of the rotated operator present at all rotation angles (the continuum moves with the angle, a mode does not)"""
    sets = [ecs_eigs("e2", n, parity, N=N, theta=th) for th in ths]
    sets = [e[np.abs(e) < wmax] for e in sets]
    return [e for e in sets[0] if all(np.min(np.abs(S - e)) < 2e-3 * max(1, abs(e)) for S in sets[1:])], [len(S) for S in sets]


# ------------------------------------------------------------------ time domain: Y_tt = Y'' - p Y' + q Y on l > 0 with parity at l = 0
class Evo:
    def __init__(self, n, h, lmax, parity="even"):
        self.n, self.h = n, h
        J = int(round(lmax / h))
        self.l = (np.arange(J) + 0.5) * h
        sol = solve_ivp(lambda l, y: [y[1], (n - 1) * y[0] ** (1 - 2 * n)], (0, self.l[-1]), [1.0, 0.0], t_eval=self.l,
                        method="DOP853", rtol=1e-13, atol=1e-15)
        self.r, self.rp = sol.y
        rpp = (n - 1) * self.r ** (1 - 2 * n)
        self.p = 2 * rpp / self.rp + (n + (1 if MUTATE == "wrong_p" else 0)) * self.rp / self.r
        self.q = n * self.rp ** 2 / self.r ** 2 + rpp / self.r
        sgn = 1.0 if parity == "even" else -1.0
        lo, di, up = 1 / h ** 2 + self.p / (2 * h), -2 / h ** 2 + self.q, 1 / h ** 2 - self.p / (2 * h)
        diag = di.copy()
        diag[0] += sgn * lo[0]
        diag[-1] += 2 * up[-1]
        lower, upper = lo[1:].copy(), up[:-1].copy()
        lower[-1] -= up[-1]
        self.A = sp_.diags([lower, diag, upper], [-1, 0, 1], format="csr")

    def unstable(self, guess):
        vals, vecs = sla.eigs(self.A, k=1, sigma=guess, which="LM")
        _, vecsT = sla.eigs(self.A.T.tocsr(), k=1, sigma=guess, which="LM")
        return vals[0].real, vecs[:, 0].real, vecsT[:, 0].real

    def run(self, Y, S, tmax, dt, proj=None, dtout=0.1):
        A = self.A
        ts, th = [0.0], [(9 * Y[0] - Y[1]) / 8]
        every = int(round(dtout / dt))
        assert abs(every * dt - dtout) < 1e-12
        for k in range(int(round(tmax / dt))):
            k1y, k1s = S, A @ Y
            k2y, k2s = S + 0.5 * dt * k1s, A @ (Y + 0.5 * dt * k1y)
            k3y, k3s = S + 0.5 * dt * k2s, A @ (Y + 0.5 * dt * k2y)
            k4y, k4s = S + dt * k3s, A @ (Y + dt * k3y)
            Y = Y + dt / 6 * (k1y + 2 * k2y + 2 * k3y + k4y)
            S = S + dt / 6 * (k1s + 2 * k2s + 2 * k3s + k4s)
            if proj is not None:
                W, v, lv = proj
                c = lv @ (S + W * Y) / (2 * W * (lv @ v))
                Y, S = Y - c * v, S - c * W * v
            if (k + 1) % every == 0:
                ts.append((k + 1) * dt)
                th.append((9 * Y[0] - Y[1]) / 8)
        return np.array(ts), np.array(th)


def evolve(n, h, tmax, S0f, W2guess):
    ev = Evo(n, h, tmax + 40.0)
    w2, v, lv = ev.unstable(W2guess)
    W = np.sqrt(w2)
    Y0, S0 = np.zeros_like(ev.l), S0f(ev.l)
    proj = None if MUTATE == "no_projection" else (W, v, lv)
    if proj is not None:
        c = lv @ (S0 + W * Y0) / (2 * W * (lv @ v))
        Y0, S0 = Y0 - c * v, S0 - c * W * v
    wt = 1 / (ev.rp ** 2 * ev.r ** n)
    P = float(np.sum(wt * ev.r * S0f(ev.l)) * h)
    nst = int(np.ceil(0.1 / (0.4 * h)))
    t, y = ev.run(Y0, S0, tmax, 0.1 / nst, proj)
    return t, y, P, w2


# ------------------------------------------------------------------ the time-dependent linear system (sympy)
def sym_system():
    import sympy as sp
    t, l, eps = sp.symbols("t l epsilon")
    n = sp.Symbol("n", positive=True)
    r0 = sp.Function("r0")(l)
    N, L, Y = (sp.Function(s)(t, l) for s in ("N", "L", "Y"))
    nu, lam, r = eps * N, eps * L, r0 + eps * Y
    x = [t, l]
    h = sp.diag(-sp.exp(2 * nu), sp.exp(2 * lam))
    hi = h.inv()
    G = [[[sum(hi[a, d] * (sp.diff(h[d, b], x[c]) + sp.diff(h[d, c], x[b]) - sp.diff(h[b, c], x[d])) for d in range(2)) / 2
           for c in range(2)] for b in range(2)] for a in range(2)]
    ricci = lambda b, c: sum(sp.diff(G[a][b][c], x[a]) - sp.diff(G[a][b][a], x[c])
                             + sum(G[a][a][d] * G[d][b][c] - G[a][c][d] * G[d][b][a] for d in range(2)) for a in range(2))
    Rh = sum(hi[a, b] * ricci(a, b) for a in range(2) for b in range(2))
    hess = lambda A, B, f: sp.diff(f, x[A], x[B]) - sum(G[C][A][B] * sp.diff(f, x[C]) for C in range(2))
    grad2 = sum(hi[a, b] * sp.diff(r, x[a]) * sp.diff(r, x[b]) for a in range(2) for b in range(2))
    box = sum(hi[a, b] * hess(a, b, r) for a in range(2) for b in range(2))
    kphi2 = n * sp.diff(r0, l, 2) / r0
    E1 = {k: Rh * h[k] / 2 - (n / r) * hess(*k, r) + (kphi2 if k == (1, 1) else 0) for k in [(0, 0), (0, 1), (1, 1)]}
    E2 = (n - 1) * (1 - grad2) - r * box
    S = sp.diff(sp.exp(nu + lam) * r ** n * hi[1, 1] * sp.sqrt(kphi2), l)
    first = lambda e: sp.simplify(sp.diff(e, eps).subs(eps, 0))
    lin = {"E1tt": first(E1[(0, 0)]), "E1tl": first(E1[(0, 1)]), "E1ll": first(E1[(1, 1)]), "E2": first(E2), "S": first(S)}
    return sp, t, l, n, r0, (N, L, Y), lin


def sym_checks(nval):
    sp, t, l, n_sym, r0, (N, L, Y), lin = sym_system()
    n = sp.Integer(nval)
    rr = sp.Function("rr")(l)
    q = sp.diff(rr, l)
    rpp = (n - 1) * rr ** (1 - 2 * n)
    I_ = sp.Function("I_")(l)       # Y2 = rr * I_, I_' = r'^2 r^(n-2)
    Z_ = sp.Function("Z_")(l)       # generalized zero mode: D Z = Y1 = rr
    R, Q, Is, Zs, Zp = sp.symbols("R Q Is Zs Zp")
    p = 2 * rpp / q + n * q / rr
    qq = n * q ** 2 / rr ** 2 + rpp / rr
    Ider = {1: q ** 2 * rr ** (n - 2)}
    Zder = {2: p * sp.diff(Z_, l) - qq * Z_ + rr}

    def bgrules(e):
        for k in (5, 4, 3):
            e = e.subs(sp.Derivative(rr, (l, k)), sp.diff((n - 1) * rr ** (1 - 2 * n), l, k - 2))
        return e.subs(sp.Derivative(rr, (l, 2)), rpp)

    def reduce(e):
        e = bgrules(sp.expand(e.doit()))
        for k in (2, 3, 4):
            if k not in Ider:
                Ider[k] = bgrules(sp.diff(Ider[k - 1], l))
        for k in (4, 3, 2, 1):
            e = bgrules(e.subs(sp.Derivative(I_, (l, k)), Ider[k]))
        for k in (3, 4):
            if k not in Zder:
                d = sp.diff(Zder[k - 1], l)
                for j in range(k - 1, 1, -1):
                    d = d.subs(sp.Derivative(Z_, (l, j)), Zder[j])
                Zder[k] = bgrules(d)
        for k in (4, 3, 2):
            e = bgrules(e.subs(sp.Derivative(Z_, (l, k)), Zder[k]))
        e = e.subs(sp.Derivative(Z_, l), Zp).subs(Z_, Zs).subs(I_, Is).subs(sp.Derivative(rr, l), Q).subs(rr, R)
        return sp.simplify(sp.together(e).subs(Q, sp.sqrt(1 - R ** (2 - 2 * n))))

    def residuals(Yc, g=0):
        Lc = sp.diff(Yc, l) / q + g
        Nc = Lc - n * Yc / rr
        out = {}
        for k, e in lin.items():
            ee = e.subs(n_sym, n).subs(r0, rr).doit()
            ee = ee.subs({N: Nc, L: Lc, Y: Yc}).doit()
            out[k] = reduce(ee)
        return out
    Y1, Y2 = rr, rr * I_
    Yphi, gphi = rr ** n * q, -(n - 1) * rr ** (1 - n) / q
    cands = {"Y1 (scale)": (Y1, 0), "t Y1": (t * Y1, 0), "Y2 (asymmetry)": (Y2, 0), "t Y2": (t * Y2, 0),
             "t^2/2 Y1 + Y^(1)": (t ** 2 / 2 * Y1 + Z_, 0), "Yphi (translation, g != 0)": (Yphi, gphi)}
    if MUTATE == "wrong_secular":
        cands["t Y1"] = (t * Y1, 0)
        Yc = t * Y1
        Lc = sp.Integer(0) * t                        # L dropped: not the solution of the (t,l) equation
        Nc = Lc - n * Yc / rr
        bad = {k: reduce(e.subs(n_sym, n).subs(r0, rr).doit().subs({N: Nc, L: Lc, Y: Yc}).doit()) for k, e in lin.items()}
        return {"t Y1": bad}, {}
    res = {k: residuals(*v) for k, v in cands.items()}
    fail = {"t Yphi (no secular partner)": residuals(t * Yphi, gphi)}
    return res, fail


def sym_reduction(nval):
    """with L = Y'/r0' + g and N = L - nY/r0 + c(t): E2 reduces to Y_tt = D Y - 2(n-1) g / r0, and E1tt, E1ll contain g only"""
    sp, t, l, n_sym, r0, (N, L, Y), lin = sym_system()
    n = sp.Integer(nval)
    g, c = sp.Function("g")(l), sp.Function("c")(t)
    rr = sp.Function("rr")(l)
    q = sp.diff(rr, l)
    rpp = (n - 1) * rr ** (1 - 2 * n)
    Yf = sp.Function("Yf")(t, l)
    Lc = sp.diff(Yf, l) / q + g
    Nc = Lc - n * Yf / rr + c
    sub = lambda e: e.subs(n_sym, n).subs(r0, rr).doit().subs({N: Nc, L: Lc, Y: Yf}).doit()
    Ytt = sp.solve(sub(lin["E2"]), sp.diff(Yf, t, 2))[0]
    R, Q = sp.symbols("R Q")
    Y0, Y1, Y2, Y3, g0, g1, g2 = sp.symbols("Y0 Y1 Y2 Y3 g0 g1 g2")

    def bgrules(e):
        for k in (4, 3):
            e = e.subs(sp.Derivative(rr, (l, k)), sp.diff((n - 1) * rr ** (1 - 2 * n), l, k - 2))
        return e.subs(sp.Derivative(rr, (l, 2)), rpp)

    def fin(e):
        e = sp.expand(e)
        e = e.subs(sp.Derivative(Yf, (t, 2), l), sp.diff(Ytt, l)).subs(sp.Derivative(Yf, l, (t, 2)), sp.diff(Ytt, l))
        e = bgrules(sp.expand(e.subs(sp.Derivative(Yf, (t, 2)), Ytt)))
        for k, s in [(sp.Derivative(Yf, (l, 3)), Y3), (sp.Derivative(Yf, (l, 2)), Y2), (sp.Derivative(Yf, l), Y1),
                     (sp.Derivative(g, (l, 2)), g2), (sp.Derivative(g, l), g1)]:
            e = e.subs(k, s)
        e = e.subs(Yf, Y0).subs(g, g0).subs(sp.Derivative(rr, l), Q).subs(rr, R)
        return sp.simplify(sp.together(e).subs(Q, sp.sqrt(1 - R ** (2 - 2 * n))))
    ytt = fin(Ytt)
    pp = 2 * (n - 1) * R ** (1 - 2 * n) / sp.sqrt(1 - R ** (2 - 2 * n)) + n * sp.sqrt(1 - R ** (2 - 2 * n)) / R
    qs = n * (1 - R ** (2 - 2 * n)) / R ** 2 + (n - 1) * R ** (-2 * n)
    pred = Y2 - pp * Y1 + qs * Y0 - 2 * (n - 1) * g0 / R
    e1, e2 = fin(sub(lin["E1tt"])), fin(sub(lin["E1ll"]))
    z = {g0: 0, g1: 0, g2: 0}
    return sp.simplify(ytt - pred), (e1.subs(z), e2.subs(z)), (e1, e2)


# ------------------------------------------------------------------ the crate's W^2 (planted)
def crate_w2(D, guess):
    out = subprocess.run([CRATE, "w2", "--lane=qd", f"--D={D}", f"--W2={guess}"], capture_output=True, text=True, timeout=600).stdout
    rec = json.loads([ln for ln in out.splitlines() if ln.startswith("RESULT")][0][len("RESULT "):])
    return float(rec["W2"]), rec["converged"]


# ------------------------------------------------------------------ the checks
def stage_planted():
    print("== planted cases (same code paths, known answers)")
    exact = {k: complex(-PT_AL * (k + 0.5), np.sqrt(PT_V0 - PT_AL ** 2 / 4)) for k in range(3)}
    worst = 0
    parts = []
    for k, W in exact.items():
        par = "even" if k % 2 == 0 else "odd"
        th = ray_angle(W)
        F = lambda w: Ma_pt(w, par, ray_angle(w), 20.0 / ((w * np.exp(1j * ray_angle(w))).real))
        r, it = secant(F, W * 1.01 + 0.01, W * 0.99 - 0.01j)
        worst = max(worst, abs(r - W))
        parts.append("n=%d %s %.2e" % (k, par, abs(r - W)))
    check("PT1", worst <= 1e-8, "Poeschl-Teller V0 = 2.5, alpha = 1, W = -i omega, first three overtones, worst |W - exact| = %.2e (bar 1e-8); %s" % (worst, "; ".join(parts)))
    ok = True
    parts = []
    for (box_, par, want) in [((-1.0, 0.0, 1.0, 2.0), "even", 1), ((-2.0, -1.0, 1.0, 2.0), "odd", 1), ((-3.0, -2.0, 1.0, 2.0), "even", 1),
                              ((-2.0, -1.0, 1.0, 2.0), "even", 0), ((-1.0, 0.0, 1.0, 2.0), "odd", 0), ((-1.0, 0.0, 2.0, 3.0), "even", 0)]:
        th, sf, m = box_params(*box_)
        w, mx = winding(lambda W: Ma_pt(W, par, th, sf), *box_, h=0.1)
        ok &= abs(w - want) < 0.05
        parts.append("%s %s: %.3f (want %d)" % (box_, par, w, want))
    check("PT2", ok, "argument-principle windings, same scan routine; " + "; ".join(parts))
    # ECS on PT
    worst = 0
    parts = []
    for k, W in exact.items():
        par = "even" if k % 2 == 0 else "odd"
        ev = ecs_eigs("pt", 0, par, N=481, S=24.0, theta=-1.2)
        d = np.min(np.abs(ev - W * W))
        worst = max(worst, d / abs(W * W))
        parts.append("n=%d %.1e" % (k, d / abs(W * W)))
    check("PT3", worst <= 1e-6, "second eigen-solver (Chebyshev on a complex-scaled contour) on Poeschl-Teller: worst relative error of W^2 = %.2e (bar 1e-6); %s" % (worst, "; ".join(parts)))
    for D, guess, n in [(5, 2.985, 3), (4, 1.397, 2)]:
        w2c, conv = crate_w2(D, guess)
        W = np.sqrt(w2c)
        F = lambda w: Ma_e2(n, w, "even", 0.0, 20.0 / w.real)
        r, it = secant(F, W * 1.01, W * 0.99)
        e1 = abs(r * r - w2c) / w2c
        check("W1-D%d" % D, conv and e1 <= 1e-8, "unstable W^2 at D = %d from the shooting solver %.14f against #195's crate (qd) %.14f: relative %.2e (bar 1e-8)" % (D, (r * r).real, w2c, e1))
        ev = ecs_eigs("e2", n, "even", N=481, theta=-0.7)
        k = np.argmin(np.abs(ev - w2c))
        e2 = abs(ev[k] - w2c) / w2c
        check("W2-D%d" % D, e2 <= 1e-4, "second eigen-solver (ECS, theta = -0.7, N = 481) %.8f%+.1ei: relative %.2e (bar 1e-4)" % (ev[k].real, ev[k].imag, e2))
        errs = []
        for h in (0.02, 0.01):
            e = Evo(n, h, 60.0)
            errs.append(abs(e.unstable(guess)[0] - w2c) / w2c)
        check("W3-D%d" % D, errs[1] <= 3e-4 and 3.0 <= errs[0] / errs[1] <= 5.0, "finite-difference operator of the time-domain code: relative error %.2e (h = 0.02), %.2e (0.01), ratio %.2f (second order: 4)" % (errs[0], errs[1], errs[0] / errs[1]))


def stage_zero():
    print("== zero modes, secular solutions, generalized kernel")
    for nval in (2, 3):
        d, (e1, e2), _ = sym_reduction(nval)
        check("Z1-n%d" % nval, d == 0 and e1 == 0 and e2 == 0,
              "n = %d (D = %d): with L = Y'/r0' + g(l), N = L - nY/r0 + c(t), E2 gives Y_tt = D Y - 2(n-1) g/r0 (difference %s); E1tt, E1ll at g = 0: %s, %s" % (nval, nval + 2, d, e1, e2))
        t0 = _time.time()
        res, fail = sym_checks(nval)
        for k, v in res.items():
            ok = all(x == 0 for x in v.values())
            check("Z2-n%d %s" % (nval, k), ok, "residuals of (E1tt, E1tl, E1ll, E2, S) = %s" % ", ".join(str(x) for x in v.values()))
        for k, v in fail.items():
            ok = any(x != 0 for x in v.values())
            check("Z2-n%d %s" % (nval, k), ok, "must NOT be a solution; residuals: %s" % ", ".join(str(x)[:60] for x in v.values()))
    # generalized zero mode: D Z = Y1, asymptotics
    for n, pred in [(3, 0.5), (2, 0.5)]:
        l0 = 0.05
        zs = []
        for W2 in (1e-4, -1e-4):
            R, Rp, Y, Yp = frob(n, W2, "even", l0)
            zs.append((Y, Yp))
        z0, zp0 = [(a - b) / 2e-4 for a, b in zip(*zs)]
        g = lambda l, s: [s[1], (n - 1) * s[0] ** (1 - 2 * n), s[3],
                          (2 * (n - 1) * s[0] ** (1 - 2 * n) / s[1] + n * s[1] / s[0]) * s[3] - (n * s[1] ** 2 / s[0] ** 2 + (n - 1) * s[0] ** (-2 * n)) * s[2] + s[0]]
        R, Rp, _, _ = frob(n, 0.0, "even", l0)
        ls = [1e2, 1e3, 1e4]
        sol = solve_ivp(g, (l0, ls[-1]), [R.real, Rp.real, z0.real, zp0.real], t_eval=ls, method="DOP853", rtol=1e-12, atol=1e-12)
        zl = [sol.y[2][i] / ls[i] ** 3 for i in range(3)]
        slope = (zl[2] - zl[1]) / np.log(ls[2] / ls[1])
        if n == 3:
            ok = abs(slope - pred) < 0.01 * pred
            txt = "Z/l^3 at l = 1e2, 1e3, 1e4: %s; per e-fold over the last decade %.4f (predicted 1/2: Z ~ (1/2) l^3 ln l)" % (", ".join("%.4f" % x for x in zl), slope)
        else:
            ok = abs(zl[2] - pred) < 0.02 * pred and abs(slope) < 0.01
            txt = "Z/l^3 at l = 1e2, 1e3, 1e4: %s (predicted 1/2: Z ~ l^3/(6 - 2n))" % ", ".join("%.4f" % x for x in zl)
        check("Z3-n%d" % n, ok, "generalized zero mode D Z = Y1 (even, Z(0) = 0), unbounded as l^3 (ln l at D = 5): " + txt)
    # threshold behaviour of the mismatch Delta(W)
    for n, lo, hi, what in [(2, 0.9, 1.2, "simple pole of the resolvent at W = 0"), (3, 1.55, 2.1, "double pole with a logarithm")]:
        ws = [0.01, 0.02, 0.04, 0.08]
        D = []
        for w in ws:
            a, b = logder_e2(n, complex(w, 0), "odd" if MUTATE == "pole_off" else "even", 0.0, 30.0 / w)
            D.append(abs(a - b))
        sl_ = np.diff(np.log(D)) / np.diff(np.log(ws))
        ok = bool(np.all((sl_ > lo) & (sl_ < hi)))
        check("Z4-n%d" % n, ok, "d ln|Delta| / d ln W at W = 0.01 ... 0.08: %s (band [%.2f, %.2f]): %s at D = %d" % (np.round(sl_, 3), lo, hi, what, n + 2))


def stage_spectrum():
    print("== poles of the resolvent: argument-principle scan (shooting) and complex-scaled collocation")
    for D, n in [(5, 3), (4, 2)]:
        F = lambda w: Ma_e2(n, w, "even", 0.0, 20.0 / w.real)
        w2c = crate_w2(D, {5: 2.985, 4: 1.397}[D])[0]
        W = np.sqrt(w2c)
        a = winding(F, W - 0.5, W + 0.5, -0.5, 0.5, 0.1)[0]
        b = winding(F, W + 0.7, W + 1.7, -0.5, 0.5, 0.1)[0]
        c = winding(F, 0.2, W - 0.7, -0.5, 0.5, 0.1)[0] if W > 0.9 else 0.0
        check("S1-D%d" % D, abs(a - 1) < 0.05 and abs(b) < 0.05 and abs(c) < 0.05, "Re W > 0, boxes around and beside the unstable root W = %.6f: windings %.3f, %.3f, %.3f (want 1, 0, 0)" % (W, a, b, c))
        for parity in ("even", "odd"):
            t0 = _time.time()
            roots, skipped = scan(n, parity)
            lst = []
            for r, w, mx in roots:
                pol = polish_stable(n, parity, r)
                lst.append((r, pol))
            print("   D = %d %s: Re W in [-8.2, -0.2], Im W in [0, 8], %d zero(s), %d cell(s) skipped near the origin, %.0f s" % (D, parity, len(roots), len(skipped), _time.time() - t0))
            for r, pol in lst:
                spread = max(np.abs(pol - r)) if len(pol) else float("nan")
                print("      zero at W = %.5f %+.5fi; over ray angle and far point (n = %d): spread %.1e" % (r.real, r.imag, len(pol), spread))
            incore = [r for r, _ in lst if -1.4 <= r.real <= -0.2 and 0 <= r.imag <= 1.0]
            check("S2-D%d-%s" % (D, parity), len(incore) == 0, "no zero with Re W in [-1.4, -0.2], Im W in [0, 1] (the region of kappa = 0.495 and 0.244): %d found" % len(incore))
            if lst:
                sp_ = max(max(np.abs(p - r)) if len(p) else 1.0 for r, p in lst)
                check("S3-D%d-%s" % (D, parity), sp_ <= 1e-3, "the %d zero(s) near the cut: worst spread over four ray angles and two far points %.1e (bar 1e-3)" % (len(lst), sp_))
        for parity in ("even", "odd"):
            iso, counts = ecs_isolated(n, parity)
            real_pos = [e for e in iso if abs(e.imag) < 1e-3 and e.real > 0.1]
            check("S4-D%d-%s" % (D, parity), (len(iso) - len([e for e in iso if abs(e) < 0.05])) == (1 if parity == "even" else 0),
                  "complex scaling at theta = -0.7, -0.85, -1.0 (N = 481): eigenvalues |W^2| < 80 present at all three: %s; counts %s" % (
                      ["%.5f%+.5fi" % (e.real, e.imag) for e in iso], counts))


def stage_time():
    print("== time-domain linear response of the throat radius (even sector, unstable mode projected out)")
    gauss = lambda c, w: (lambda l: np.exp(-((l - c) / w) ** 2 / 2))
    cases = {"c=3 w=0.5": gauss(3, 0.5), "c=6 w=0.5": gauss(6, 0.5), "c=3 w=1.0": gauss(3, 1.0), "c=1.5 w=0.5": gauss(1.5, 0.5)}
    tmax = 40.0
    for D, n, guess, kappa in [(4, 2, 1.397, 0.244), (5, 3, 2.985, 0.495)]:
        ext = {}
        for name, f in cases.items():
            ra, rb = evolve(n, 0.01, tmax, f, guess), evolve(n, 0.005, tmax, f, guess)
            ext[name] = (ra[0], (4 * rb[1] - ra[1]) / 3, rb[2], rb[1])
        r2 = evolve(n, 0.02, tmax, cases["c=3 w=0.5"], guess)
        ya, yb = evolve(n, 0.01, tmax, cases["c=3 w=0.5"], guess)[1], ext["c=3 w=0.5"][3]
        orders = [np.log2((r2[1][i] - ya[i]) / (ya[i] - yb[i])) for i in (100, 200)]
        corr = max(abs(ext["c=3 w=0.5"][1][i] - yb[i]) / abs(ext["c=3 w=0.5"][1][i]) for i in (100, 200))
        check("T1-D%d" % D, all(1.8 <= o <= 2.3 for o in orders) and corr <= 1e-2 and np.all(np.isfinite(yb)),
              "h = 0.02, 0.01, 0.005: observed order of Y(t, 0) at t = 10, 20: %.2f, %.2f (band [1.8, 2.3]); the Richardson correction is %.1e of Y (bar 1e-2); Y_ext(10), Y_ext(20) = %.5f, %.5f"
              % (orders[0], orders[1], corr, ext["c=3 w=0.5"][1][100], ext["c=3 w=0.5"][1][200]))
        t, yy, P, _ = ext["c=3 w=0.5"]
        v = np.gradient(yy, t)
        loc = lambda tt, vv: (np.log(abs(vv[int(round(tt / 0.1)) + 5])) - np.log(abs(vv[int(round(tt / 0.1)) - 5])))
        rows = ["t=%d: Y %.4f, v %.4e, rate %.3f" % (tt, yy[int(round(tt / 0.1))], v[int(round(tt / 0.1))], -loc(tt, v)) for tt in (8, 10, 12, 15, 20, 25)]
        print("   D = %d (c=3 w=0.5, Richardson of h = 0.01, 0.005): %s" % (D, " | ".join(rows)))
        m = (t >= 6.5) & (t <= 12)
        fit = lambda tt, a, c, k: a + c * np.exp(-k * (tt - 6.5))
        pf, cov = curve_fit(fit, t[m], v[m], p0=[v[m][-1], v[m][0] - v[m][-1], 0.4], maxfev=40000)
        late = []
        for name in cases:
            tn, yn, _, _ = ext[name]
            late.append(-loc(20, np.gradient(yn, tn)))
        spread = max(late) - min(late)
        print("   D = %d: the #255 estimator v = v_inf + c exp(-k (t - 6.5)) on t in [6.5, 12]: v_inf = %.4e, k = %.4f; late rates at t = 20 over four data: %s" % (D, pf[0], pf[2], np.round(late, 3)))
        tol = 0.01 if D == 5 else 0.024
        check("T2-D%d" % D, abs(pf[2] - kappa) <= tol, "k of the #255 estimator on the linear response = %.4f against kappa = %.3f (bar +-%.3f)" % (pf[2], kappa, tol))
        if D == 4:
            lo, hi = 0.205, 0.222      # #267: -ln(0.649 ... 0.669) / (0.975 * 2.00)
            rate = float(np.mean(late))
            check("T2m-D4", lo * 0.9 <= rate <= hi * 1.1, "late rate of the throat speed at t = 20 (mean over four data) %.3f against #267's mass-excess rate %.3f to %.3f per proper time (-ln of 0.649 ... 0.669 per decade over 1.95; bar +-10%%)" % (rate, lo, hi))
        check("T3-D%d" % D, spread <= 0.02, "late decay rate of the throat speed at t = 20 over four data (centers 3, 6, 1.5, widths 0.5, 1): %s, spread %.3f (bar 0.02)" % (np.round(late, 3), spread))
        a, b = ext["c=3 w=0.5"], ext["c=6 w=0.5"]
        col = abs(a[1][200] / a[2] - b[1][200] / b[2]) / abs(a[1][200] / a[2])
        check("T4-D%d" % D, col <= 0.01, "Y(t, 0)/P for two data far from the throat (P = integral w Y1 sigma0, w = 1/(r'^2 r^n)) at t = 20 differ by %.2e (bar 1e-2): the response is that of one conjugate momentum, Y/P = %.4f" % (col, a[1][200] / a[2]))
        pt = [-loc(tt, v) * tt for tt in (15, 20, 25)]
        sp_pt = (max(pt) - min(pt)) / abs(np.mean(pt))
        check("T7-D%d" % D, sp_pt <= 0.12 and np.mean(pt) > 0, "the late throat speed is a power law, not an exponential: -d ln v / d ln t = rate x t at t = 15, 20, 25: %s, spread %.1f %% of the mean (bar 12 %%; an exponential's grows as t, 67 %% over this range), exponent %.2f" % (np.round(pt, 3), 100 * sp_pt, np.mean(pt)))
        if D == 5:
            r15, r25 = -loc(15, v), -loc(25, v)
            check("T5-D5", abs(r25) <= 0.005 and abs(r15) <= 0.005, "the secular speed is constant: d ln v/dt at t = 15, 25 = %.4f, %.4f (a constant speed has 0; bar 0.005); v/P at t = 10, 20 = %.4f, %.4f" % (-r15, -r25, v[100] / P, v[200] / P))
        else:
            sat = (yy[400] - yy[300]) / yy[300]
            check("T5-D4", abs(sat) <= 1e-2, "the shift saturates: Y(40)/Y(30) - 1 = %.2e (bar 1e-2); Y(30)/P = %.4f" % (sat, yy[300] / P))
    B = np.array([0.02, 0.04, 0.07, 0.10])
    rdot = np.array([0.00021, 0.00143, 0.00898, 0.02601])
    Enet = np.array([0.01548, 0.06186, 0.18881, 0.38277])
    ratio = rdot / Enet
    check("T6", ratio.max() / ratio.min() <= 1.2, "rdot_inf / E_net at B = 0.02, 0.04, 0.07, 0.1 (throat_drift.md's table): %s, max/min %.2f (a speed linear in the injected energy needs <= 1.2)" % (np.round(ratio, 4), ratio.max() / ratio.min()))


def main():
    global MUTATE, ONLY, CRATE
    args = sys.argv[1:]
    for a in args:
        if a.startswith("--mutate="):
            MUTATE = a.split("=", 1)[1]
            assert MUTATE in MUTATIONS, MUTATE
        elif a.startswith("--only="):
            ONLY = a.split("=", 1)[1].split(",")
        elif a.startswith("--crate="):
            CRATE = a.split("=", 1)[1]
    rev = subprocess.run(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    print("ellis_spectrum.py  code %s  mutate=%s  crate=%s\n" % (rev, MUTATE, os.path.relpath(CRATE, HERE)), flush=True)
    stages = {"planted": stage_planted, "zero": stage_zero, "spectrum": stage_spectrum, "time": stage_time}
    for k, f in stages.items():
        if ONLY is None or k in ONLY:
            try:
                f()
            except Exception as e:
                check("stage-" + k, False, "stopped by %s: %s" % (type(e).__name__, str(e)[:200]))
    bad = [n for n, ok in RESULTS if not ok]
    print("\n%d checks, %d pass, %d fail%s" % (len(RESULTS), len(RESULTS) - len(bad), len(bad), (": " + ", ".join(bad)) if bad else ""))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
