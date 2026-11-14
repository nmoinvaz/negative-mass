#!/usr/bin/env python3
"""Issue #66: #58's analysis repeated for Natario's zero-expansion warp drive [Nat02].

The t = 0 slice of ds^2 = -dt^2 + (dx^i - X^i dt)^2 is flat, with lapse 1, shift beta = -X and
K_ij = -(d_i X_j + d_j X_i)/2 (the sign of #58). X = v curl(n(r) e_x x r) is divergence free, so tr K = 0,
with n = 1/2 inside the bubble and 0 outside. Here n = f/2 with f the tanh profile of #58 (R = 1), so that
(sigma R, v) match #58's table. G = c = 1.

    python notes/warp/natario.py [--part=all|sym|adm|split|energy|identity|mutations] [--mutate=NAME] [--quick]

Parts:
  sym        mu and J from the 3+1 constraints and from G_ab of the 4-metric, against closed forms and [Nat02]'s rho.
  adm        E and P as surface integrals, with #58's controls.
  split      the pointwise minimal split against a numerical minimisation, on random and on Natario (mu, J).
  energy     quadrature controls, then E_- and the Eulerian energy against v and wall width, beside #58's.
  identity   the Witten identity with #58's FFT solver (ported, K as an argument), controls first.
  mutations  every mutation in a subprocess, each of which must make its named check fail.

--quick shrinks the energy table and the identity grids (used by the mutation runs).
Exits 1 when any check fails.
"""
import argparse
import concurrent.futures
import contextlib
import importlib.util
import io
import math
import os
import re
import subprocess
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--part", default="all")
ap.add_argument("--mutate", default="")
ap.add_argument("--quick", action="store_true")
args = ap.parse_args()
MUT = args.mutate
QUICK = args.quick
FAILED = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILED.append(name)


def load_58():
    """#58's script, unmodified and not run (--part=none), for its profile, nodes, ADM and solver code."""
    path = os.path.join(HERE, "alcubierre_bounds.py")
    spec = importlib.util.spec_from_file_location("alcubierre_bounds", path)
    mod = importlib.util.module_from_spec(spec)
    saved, sys.argv = sys.argv, [path, "--part=none"]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            spec.loader.exec_module(mod)
    except SystemExit:
        pass
    finally:
        sys.argv = saved
    return mod


AB = load_58()

# check names, shared by the checks and by the mutation table
C_DIV = "X is divergence free"
C_NATX = "X = [Nat02] Sec. 2 field + v e_x"
C_MU3 = "3+1 mu = closed form"
C_J3 = "3+1 J = closed form"
C_MU4 = "4D mu = closed form"
C_J4 = "4D J = closed form"
C_NATRHO = "closed-form mu = [Nat02] rho"
C_NUM = "numeric spherical fields = symbolic Cartesian fields"
C_PROF = "tanh profile f', f'', f''' against finite differences and #58's"
C_GAM = "gamma algebra: sum_i |K_ij g^j g0 e|^2 = |K|^2 |e|^2"
C_ADM = "Natario E = 0 and |P| < 1e-12 at R = 8"
C_SCHW = "control: Schwarzschild E = m psi(R)^3"
C_BY = "control: Bowen-York P returned"
C_SPL1 = "minimal mu_- = max(0, -mu, (|J| - mu)/2) at random (mu, J)"
C_SPL2 = "minimal mu_- at Natario (mu, J)"
C_SPLA = "the named split is admissible"
C_GAUSS = "control: Gaussian Eulerian energy"
C_REPRO = "control: this quadrature reproduces #58's table"
C_INTJ = "int J dV = 0 (quadrature)"
C_CONV = "E_- converged (relative change < 1e-4 on doubling the nodes)"
C_THIN = "thin-wall hand estimates (sigma R = 128)"
C_PORT = "control: ported solver reproduces #58's"
C_EPS0 = "solver returns eps0, residual at rounding level"
C_ALG = "pointwise |D~eps0|^2/4pi + mu |eps0|^2 = 0 at every node"
C_WITTEN = "Witten identity G + int T = E - P.xi0 = 0"
C_INTJG = "int J dV = 0 on the grid"
C_TPLUS = "split sectors T+ >= 0 pointwise"
C_GEUL = "grid E_Eulerian = quadrature E_Eulerian"
C_GEM = "grid E_- = quadrature E_-"


# ------------------------------------------------------------------------------------------------ profiles
def tanh_profile(sig, R=1.0):
    """f and f', f'', f''' for #58's profile (f = 1 inside, 0 outside). The warp function is n = f / 2."""
    T = math.tanh(sig * R)

    def f(r):
        ta, tb = np.tanh(sig * (r + R)), np.tanh(sig * (r - R))
        sa, sb = 1 - ta**2, 1 - tb**2
        F3 = sig**3 * (sa * (6 * ta**2 - 2) - sb * (6 * tb**2 - 2)) / (2 * T)
        if MUT == "prof3":
            F3 = -F3
        return (ta - tb) / (2 * T), sig * (sa - sb) / (2 * T), sig**2 * (-2 * sa * ta + 2 * sb * tb) / (2 * T), F3
    return f


def gauss_profile(r):
    e = np.exp(-r**2)
    return e, -2 * r * e, (4 * r**2 - 2) * e, (12 * r - 8 * r**3) * e


def fields_nat(f, v, r, u):
    """mu, J_x and J_perp (cylindrical) at radius r and u = cos(theta), about the x axis."""
    F0, F1, F2, F3 = f(r)
    n1, n2, n3 = F1 / 2, F2 / 2, F3 / 2
    s2 = 1 - u**2
    mu = -(v**2 / (16 * np.pi)) * ((5 if MUT == "mu6" else 6) * n1**2 * u**2 + 0.5 * (2 * n1 + r * n2) ** 2 * s2)
    w = 4 * n1 + r * n2
    rw1 = 4 * n1 + 6 * r * n2 + (0 if MUT == "rw1" else 1) * r**2 * n3
    p = v / (16 * np.pi * r)
    return mu, -p * (2 * w * u**2 + rw1 * s2), p * u * np.sqrt(s2) * (rw1 - 2 * w)


def minus_density(mu, Jn):
    """Pointwise minimal mu_- over splits with mu_+- >= |J_+-|, and the |J_-| that attains it."""
    if MUT == "split_old":
        m = np.maximum(0.0, (Jn - mu) / 2)
        return m, m
    m = np.maximum.reduce([np.zeros_like(mu), -mu, (Jn - mu) / 2])
    s = np.where(Jn >= -mu, m, Jn)
    return m, s


def gammas():
    g0, g = AB.gammas_np()
    if MUT == "gamdup":
        g = [g[0], g[1], g[0]]
    return g0, g


# ------------------------------------------------------------------------------------------------ sym
SYM = {}


def stack(o, shape):
    if isinstance(o, (list, tuple)):
        return np.stack([stack(c, shape) for c in o])
    return np.broadcast_to(np.asarray(o, dtype=float), shape)


def ev(fn, X, r, n1, n2, n3, v):
    return stack(fn(X[0], X[1], X[2], r, n1, n2, n3, v), r.shape)


def get_sym():
    """X, the closed forms and the 3+1 fields as symbols, built once. N[k] stands for n^(k)(r_s)."""
    if SYM:
        return SYM
    t, x, y, z = sp.symbols("t x y z", real=True)
    v, r = sp.symbols("v r", positive=True)
    N = sp.symbols("N0:7", real=True)
    xs = x - v * t
    dr = {t: -v * xs / r, x: xs / r, y: y / r, z: z / r}

    def d(e, c):
        # r and N[k] stand for r_s and n^(k)(r_s), so the chain rule is applied by hand
        out = sp.diff(e, c) + sp.diff(e, r) * dr[c]
        for k in range(6):
            out += sp.diff(e, N[k]) * N[k + 1] * dr[c]
        return out

    def is_zero(e):
        e = sp.together(sp.expand(e.subs(t, 0)))
        num = sp.expand(sp.numer(e))
        rem = sp.Poly(num, z).rem(sp.Poly(z**2 - (r**2 - x**2 - y**2), z))
        return sp.expand(rem.as_expr()) == 0

    rho2 = y**2 + z**2
    cy = sp.Rational(1, 2) if MUT == "div" else 1
    X = [v * (2 * N[0] + N[1] * rho2 / r), -cy * v * N[1] * xs * y / r, -cy * v * N[1] * xs * z / r]
    X3 = (x, y, z)
    # closed forms, hand-derived from the rate-of-strain tensor of X in spherical components
    w = 4 * N[1] + r * N[2]
    rw1 = 4 * N[1] + 6 * r * N[2] + r**2 * N[3]
    mu_c = -(v**2 / (16 * sp.pi)) * (6 * N[1] ** 2 * xs**2 / r**2
                                     + sp.Rational(1, 2) * (2 * N[1] + r * N[2]) ** 2 * rho2 / r**2)
    J_c = [(v / (16 * sp.pi * r**3)) * (-2 * w * xs * (xs, y, z)[i] + rw1 * (-rho2, xs * y, xs * z)[i])
           for i in range(3)]
    # route 1 at t = 0: 3+1 constraints on the flat slice
    beta = [-c for c in X]
    K = [[(d(beta[j], X3[i]) + d(beta[i], X3[j])) / 2 for j in range(3)] for i in range(3)]
    args_ = (x, y, z, r, N[1], N[2], N[3], v)
    t0 = {t: 0}
    SYM.update(t=t, x=x, y=y, z=z, v=v, r=r, N=N, xs=xs, d=d, is_zero=is_zero, X=X, X3=X3, mu_c=mu_c, J_c=J_c,
               K=K, w=w, rw1=rw1,
               K_f=sp.lambdify(args_, [[K[i][j].subs(t0) for j in range(3)] for i in range(3)], "numpy"),
               mu_f=sp.lambdify(args_, mu_c.subs(t0), "numpy"),
               J_f=sp.lambdify(args_, [c.subs(t0) for c in J_c], "numpy"))
    return SYM


def part_sym():
    print("== sym: mu and J by two routes, and against [Nat02]")
    S = get_sym()
    t, x, y, z, v, r, N, xs = (S[k] for k in ("t", "x", "y", "z", "v", "r", "N", "xs"))
    d, is_zero, X, X3, mu_c, J_c = S["d"], S["is_zero"], S["X"], S["X3"], S["mu_c"], S["J_c"]
    rho2 = y**2 + z**2
    check(C_DIV, is_zero(sum(d(X[i], X3[i]) for i in range(3))))

    # [Nat02] Sec. 2, in the frame comoving with the bubble: X_N = -2 v f cos(th) e_r + v (2f + r f') sin(th) e_th,
    # with f = 0 inside and 1/2 outside. X_N = X - v e_x with n = 1/2 - f.
    FN = sp.symbols("FN0:4", real=True)
    cth_er = [xs * c / r**2 for c in (xs, y, z)]
    sth_eth = [-rho2 / r**2, xs * y / r**2, xs * z / r**2]
    sg = -1 if MUT == "natX" else 1
    XN = [-2 * v * FN[0] * cth_er[i] + sg * v * (2 * FN[0] + r * FN[1]) * sth_eth[i] for i in range(3)]
    to_f = {N[0]: sp.Rational(1, 2) - FN[0], N[1]: -FN[1], N[2]: -FN[2], N[3]: -FN[3]}
    check(C_NATX, all(is_zero(X[i].subs(to_f) - (v if i == 0 else 0) - XN[i]) for i in range(3)))
    c3 = 2 if MUT == "natrho" else 3
    rho_N = -(v**2 / (8 * sp.pi)) * (c3 * FN[1] ** 2 * xs**2 / r**2 + (FN[1] + r * FN[2] / 2) ** 2 * rho2 / r**2)
    check(C_NATRHO, is_zero(mu_c.subs(to_f) - rho_N))

    bscale = 2 if MUT == "beta2" else 1
    beta = [-bscale * c for c in X]
    K = [[(d(beta[j], X3[i]) + d(beta[i], X3[j])) / 2 for j in range(3)] for i in range(3)]
    trK = sum(K[i][i] for i in range(3))
    mu3 = (trK**2 - sum(K[i][j] ** 2 for i in range(3) for j in range(3))) / (16 * sp.pi)
    sJ = -1 if MUT == "Jsign" else 1
    J3 = [sJ * (sum(d(K[j][i], X3[j]) for j in range(3)) - d(trK, X3[i])) / (8 * sp.pi) for i in range(3)]
    print(f"  tr K = 0 identically: {is_zero(trK)}")
    check(C_MU3, is_zero(mu3 - mu_c))
    check(C_J3, all(is_zero(J3[i] - J_c[i]) for i in range(3)))

    # route 2: Einstein tensor of the 4-metric (X depends on t through x - v t), Eulerian projections
    X4 = (t, x, y, z)
    X2 = sum(c**2 for c in X)
    g = sp.zeros(4, 4)
    gi = sp.zeros(4, 4)
    g[0, 0], gi[0, 0] = -1 + X2, -1
    for i in range(3):
        g[0, i + 1] = g[i + 1, 0] = -X[i]
        gi[0, i + 1] = gi[i + 1, 0] = -X[i]
        g[i + 1, i + 1] = 1
        for j in range(3):
            gi[i + 1, j + 1] = (1 if i == j else 0) - X[i] * X[j]
    assert (g * gi).applyfunc(sp.expand) == sp.eye(4)
    dg = [[[d(g[a, b], X4[c]) for c in range(4)] for b in range(4)] for a in range(4)]
    Gam = [[[sp.expand(sum(gi[a, e] * (dg[e][b][c] + dg[e][c][b] - dg[b][c][e]) for e in range(4)) / 2)
             for c in range(4)] for b in range(4)] for a in range(4)]
    Ric = sp.zeros(4, 4)
    for b in range(4):
        for c in range(b, 4):
            e = sum(d(Gam[a][b][c], X4[a]) - d(Gam[a][b][a], X4[c]) for a in range(4))
            e += sum(Gam[a][a][q] * Gam[q][b][c] - Gam[a][c][q] * Gam[q][b][a] for a in range(4) for q in range(4))
            Ric[b, c] = Ric[c, b] = sp.expand(e)
    Rs = sp.expand(sum(gi[a, b] * Ric[a, b] for a in range(4) for b in range(4)))
    G = (Ric - g * Rs / 2).applyfunc(sp.expand)
    sn = -1 if MUT == "nusign" else 1
    nu = [1, sn * X[0], sn * X[1], sn * X[2]]
    mu4 = sum(G[a, b] * nu[a] * nu[b] for a in range(4) for b in range(4)) / (8 * sp.pi)
    J4 = [-sum(nu[a] * G[a, i + 1] for a in range(4)) / (8 * sp.pi) for i in range(3)]
    check(C_MU4, is_zero(mu4 - mu_c))
    check(C_J4, all(is_zero(J4[i] - J_c[i]) for i in range(3)))

    # the numpy forms of the energy table against the symbolic Cartesian forms at random points
    rng = np.random.default_rng(66)
    P = rng.uniform(-1.7, 1.7, size=(3, 600))
    rr = np.sqrt(np.sum(P**2, axis=0))
    keep = (rr > 0.3) & (rr < 1.8)
    P, rr = P[:, keep], rr[keep]
    f, vv = tanh_profile(4.0), 0.7
    F = f(rr)
    mu_s = ev(S["mu_f"], P, rr, F[1] / 2, F[2] / 2, F[3] / 2, vv)
    J_s = ev(S["J_f"], P, rr, F[1] / 2, F[2] / 2, F[3] / 2, vv)
    rho = np.hypot(P[1], P[2])
    mu_n, Jx_n, Jp_n = fields_nat(f, vv, rr, P[0] / rr)
    Jp_s = (J_s[1] * P[1] + J_s[2] * P[2]) / rho
    err = max(np.max(np.abs(mu_n - mu_s)) / np.max(np.abs(mu_s)), np.max(np.abs(Jx_n - J_s[0])) / np.max(np.abs(J_s)),
              np.max(np.abs(Jp_n - Jp_s)) / np.max(np.abs(J_s)))
    check(C_NUM, err < 1e-12, f"max relative difference {err:.1e} at {rr.size} points")

    # on the x axis inside the wall
    on_axis = {y: 0, z: 0, r: x, xs: x}
    print(f"  on the x axis (theta = 0): mu = {sp.simplify(mu_c.subs(t, 0).subs(on_axis))}")

    # profile derivatives against finite differences, and f, f', f'' against #58's profile
    r0, h = np.array([0.3, 0.9, 1.1, 1.7]), 1e-5
    ok = True
    for sig in (2.0, 4.0, 16.0):
        fp, fa = tanh_profile(sig), AB.tanh_profile(sig)
        fd = [(fp(r0 + h)[k] - fp(r0 - h)[k]) / (2 * h) for k in range(3)]
        ok = ok and all(np.max(np.abs(fd[k] - fp(r0)[k + 1])) < 1e-5 * sig ** (k + 1) for k in range(3))
        ok = ok and all(np.max(np.abs(np.array(fp(r0)[k]) - np.array(fa(r0)[k]))) < 1e-13 * sig**k for k in range(3))
    check(C_PROF, ok)

    # K_ij g^j g0 and the Witten gradient: sum_i |K_ij g^j g0 e|^2 = tr(K^2) |e|^2 for symmetric K
    g0s, gs = AB.gammas_sym()
    if MUT == "gamdup":
        gs = [gs[0], gs[1], gs[0]]
    k = sp.symbols("k11 k12 k13 k22 k23 k33", real=True)
    er, ei = sp.symbols("er0:4", real=True), sp.symbols("ei0:4", real=True)
    Ks = sp.Matrix([[k[0], k[1], k[2]], [k[1], k[3], k[4]], [k[2], k[4], k[5]]])
    eps = sp.Matrix([er[a] + sp.I * ei[a] for a in range(4)])
    tot = 0
    for i in range(3):
        A = sum((Ks[i, j] * gs[j] for j in range(3)), sp.zeros(4, 4)) * g0s * eps
        tot += (A.H * A)[0, 0]
    want = sum(Ks[i, j] ** 2 for i in range(3) for j in range(3)) * (eps.H * eps)[0, 0]
    check(C_GAM, sp.expand(tot - want) == 0)


# ------------------------------------------------------------------------------------------------ adm
def by_K(Pv, scale=1.0):
    def K(X):
        r = np.sqrt(np.sum(X**2, axis=0))
        n = X / r
        Pn = sum(Pv[i] * n[i] for i in range(3))
        return np.array([[scale * 1.5 / r**2 * (Pv[i] * n[j] + Pv[j] * n[i] - ((1 if i == j else 0) - n[i] * n[j]) * Pn)
                          for j in range(3)] for i in range(3)])
    return K


def part_adm():
    print("== adm: E and P as surface integrals")
    S = get_sym()
    f, Pv = tanh_profile(4.0), np.array([0.3, -0.2, 0.5])
    tail = by_K(Pv)

    def K_nat(X):
        r = np.sqrt(np.sum(X**2, axis=0))
        F = f(r)
        K = ev(S["K_f"], X, r, F[1] / 2, F[2] / 2, F[3] / 2, 1.0)
        return K + tail(X) if MUT == "adm_tail" else K
    flat = lambda X: np.broadcast_to(np.eye(3)[:, :, None, None], (3, 3) + X.shape[1:])  # noqa: E731
    for R in (2.0, 4.0, 8.0):
        E, P = AB.adm_E(flat, R), AB.adm_P(K_nat, R)
        print(f"  Natario (sigma R = 4, v = 1), R = {R}: E = {E:.3e}, P = {np.array2string(P, precision=3)}")
    check(C_ADM, abs(E) < 1e-12 and np.max(np.abs(P)) < 1e-12)
    m = 0.7

    def schw(X):
        r = np.sqrt(np.sum(X**2, axis=0))
        return np.eye(3)[:, :, None, None] * (1 + m / (2 * r) * (1.1 if MUT == "adm_m" else 1.0)) ** 4
    E = AB.adm_E(schw, 5.0)
    want = m * (1 + m / 10) ** 3
    print(f"  control, Schwarzschild m = {m}, R = 5: E = {E:.8f}, exact m psi(R)^3 = {want:.8f}")
    check(C_SCHW, abs(E - want) < 1e-6)
    P = AB.adm_P(by_K(Pv, 1.1 if MUT == "adm_p" else 1.0), 3.0)
    print(f"  control, Bowen-York P = {Pv}: surface integral {np.array2string(P, precision=10)}")
    check(C_BY, np.max(np.abs(P - Pv)) < 1e-10)


# ------------------------------------------------------------------------------------------------ split
def part_split():
    print("== split: minimal mu_- against a numerical minimisation")
    from scipy.optimize import minimize
    rng = np.random.default_rng(66)

    def worst_of(points):
        worst, bad = 0.0, 0
        for mu, J in points:
            Jn = np.linalg.norm(J)

            def obj(Jm):
                return max(np.linalg.norm(Jm), np.linalg.norm(J + Jm) - mu)
            opts = {"xatol": 1e-12, "fatol": 1e-12, "maxiter": 20000}
            best = min((minimize(obj, x0, method="Nelder-Mead", options=opts)
                        for x0 in [np.zeros(3), -J / 2, -J, rng.normal(size=3)]), key=lambda o: o.fun)
            for _ in range(3):  # restarts, since Nelder-Mead stalls on the kink of the max
                best = min(best, minimize(obj, best.x, method="Nelder-Mead", options=opts), key=lambda o: o.fun)
            m, s = minus_density(np.array(mu), np.array(Jn))
            worst = max(worst, abs(best.fun - float(m)) / (1 + abs(best.fun)))
            Jm = -float(s) * J / Jn
            adm = float(m) >= np.linalg.norm(Jm) - 1e-12 and mu + float(m) >= np.linalg.norm(J + Jm) - 1e-12
            bad += (not adm)
        return worst, bad
    pts = [(rng.normal() * rng.choice([0.1, 1, 3]), rng.normal(size=3) * rng.choice([0.1, 1, 3])) for _ in range(150)]
    nneg = sum(-mu > np.linalg.norm(J) for mu, J in pts)
    w1, b1 = worst_of(pts)
    print(f"  150 random (mu, J), {nneg} with -mu > |J|: max relative difference {w1:.2e}")
    check(C_SPL1, w1 < 1e-6)
    pts, nneg = [], 0
    for _ in range(150):
        sig = rng.choice([4.0, 8.0, 16.0])
        r = 1 + rng.uniform(-5, 5) / sig
        u = rng.uniform(-1, 1)
        v = 10 ** rng.uniform(-2, 1)
        mu, Jx, Jp = fields_nat(tanh_profile(sig), v, np.array(r), np.array(u))
        pts.append((float(mu), np.array([float(Jx), float(Jp), 0.0])))
        nneg += (-float(mu) > np.hypot(float(Jx), float(Jp)))
    w2, b2 = worst_of(pts)
    print(f"  150 Natario (mu, J), sigma R in 4, 8, 16, v in 0.01-10, {nneg} with -mu > |J|: max relative "
          f"difference {w2:.2e}")
    check(C_SPL2, w2 < 1e-6)
    check(C_SPLA, b1 + b2 == 0)


# ------------------------------------------------------------------------------------------------ energy
def volume_integrals(fields, f, v, sig, npan=40, nper=16, nu=64):
    r, wr = AB.radial_nodes(sig, npan=npan, nper=nper)
    u, wu = np.polynomial.legendre.leggauss(nu)
    if MUT == "quad":
        wu = wu * (1 + 0.01 * u**2)
    Rr, U = np.meshgrid(r, u, indexing="ij")
    W = 2 * np.pi * np.outer(wr * r**2, wu)
    mu, Jx, Jp = fields(f, v, Rr, U)
    Jn = np.hypot(Jx, Jp)
    m, _ = minus_density(mu, Jn)
    return {"E_eul": np.sum(mu * W), "E_minus": np.sum(m * W), "deficit": np.sum(np.maximum(0, Jn - mu) * W),
            "int_absJ": np.sum(Jn * W), "int_Jx": np.sum(Jx * W)}


def parse_58_table():
    num = r"\|\s+([-\d.e+]+)\s+"
    pat = re.compile(r"^\s+(\d+\.\d) \|\s+(\d+\.\d) " + num * 4 + r"\|\s+([-\d.e+]+)$")
    rows = {}
    with open(os.path.join(HERE, "alcubierre_bounds_output.txt")) as fh:
        for line in fh:
            m = pat.match(line.rstrip("\n"))
            if m:
                rows[(float(m[1]), float(m[2]))] = (float(m[3]), float(m[4]))
    return rows


def part_energy():
    print("== energy: quadrature controls, then the tanh profile")
    # planted truth 1, exact Eulerian energy for n = exp(-r^2) / 2
    rr = sp.symbols("r", positive=True)
    n = sp.exp(-rr**2) / 2
    n1, n2 = sp.diff(n, rr), sp.diff(n, rr, 2)
    vv = 0.8
    shell = rr**2 * (4 * n1**2 + sp.Rational(2, 3) * (2 * n1 + rr * n2) ** 2)
    exact = float(-sp.Rational(1, 8) * vv**2 * sp.integrate(shell, (rr, 0, sp.oo)))
    got = volume_integrals(fields_nat, gauss_profile, vv, 1.0)["E_eul"]
    print(f"  control, f = exp(-r^2), v = {vv}: quadrature {got:.12f}, exact -(v^2/8) int r^2 [4 n'^2 + "
          f"(2/3)(2n' + r n'')^2] dr = {exact:.12f}")
    check(C_GAUSS, abs(got - exact) < 1e-10 * abs(exact))
    # planted truth 2, the same quadrature on Alcubierre's fields returns #58's stored table
    ref = parse_58_table()
    worst = 0.0
    for (sig, v), (eul, em) in sorted(ref.items()):
        o = volume_integrals(AB.fields_sph, AB.tanh_profile(sig), v, sig, npan=80, nper=24, nu=96)
        worst = max(worst, abs(o["E_eul"] - eul) / abs(eul), abs(o["E_minus"] - em) / abs(em))
    print(f"  control, Alcubierre fields through this quadrature, {len(ref)} rows of alcubierre_bounds_output.txt: "
          f"max relative difference {worst:.1e}")
    check(C_REPRO, len(ref) == 20 and worst < 1e-5)

    sigs = (4.0, 16.0) if QUICK else (2.0, 4.0, 8.0, 16.0, 32.0)
    vs = (0.1, 1.0) if QUICK else (0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0)
    print("  sigma R |   v   | E_Eulerian   | E_- (= W-)   | int(|J|-mu)+ | E_-/(v sigma^2) | E_-/|E_Eul| | "
          "E_- of #58   | ratio Nat/Alc | rel. change, nodes doubled")
    chs, old, rows, intj = [], [], {}, 0.0
    for sig in sigs:
        f = tanh_profile(sig)
        for v in vs:
            nodes = dict(npan=4, nper=4, nu=8) if MUT == "coarse" else dict(npan=80, nper=24, nu=96)
            b = volume_integrals(fields_nat, f, v, sig, **nodes)
            d = volume_integrals(fields_nat, f, v, sig, npan=160, nper=32, nu=192)
            c = volume_integrals(AB.fields_sph, AB.tanh_profile(sig), v, sig, npan=80, nper=24, nu=96)
            ch = max(abs(b[k] - d[k]) / abs(d[k]) for k in ("E_eul", "E_minus", "deficit"))
            chs.append(ch)
            a0 = volume_integrals(fields_nat, f, v, sig)
            old.append(max(abs(a0[k] - b[k]) / abs(b[k]) for k in ("E_eul", "E_minus", "deficit")))
            intj = max(intj, abs(b["int_Jx"]) / b["int_absJ"])
            rows[(sig, v)] = b["E_minus"]
            print(f"  {sig:7.1f} | {v:5.3f} | {b['E_eul']:12.5e} | {b['E_minus']:12.5e} | {b['deficit']:12.5e} | "
                  f"{b['E_minus'] / (v * sig**2):15.5e} | {b['E_minus'] / abs(b['E_eul']):11.4e} | "
                  f"{c['E_minus']:12.5e} | {b['E_minus'] / c['E_minus']:13.4e} | {ch:.1e}")
    check(C_CONV, max(chs) < 1e-4, f"worst {max(chs):.1e}")
    print(f"  for reference, the nodes (40, 16, 64) of #58's protocol differ from the reported (80, 24, 96) by up "
          f"to {max(old):.1e}")
    check(C_INTJ, intj < 1e-10, f"max |int J_x dV| / int |J| dV = {intj:.1e}")
    print("  local exponent p = d ln E_- / d ln v between successive v")
    print("  sigma R | " + " | ".join(f"{a:g}-{b:g}" for a, b in zip(vs[:-1], vs[1:])))
    for sig in sigs:
        print(f"  {sig:7.1f} | " + " | ".join(f"{math.log(rows[(sig, b)] / rows[(sig, a)]) / math.log(b / a):.3f}"
                                           for a, b in zip(vs[:-1], vs[1:])))
    # thin-wall hand estimates: E_eul -> -v^2 sigma^3 / 180 and, for v sigma << 1, E_- -> pi v sigma^2 / (24 sqrt 3)
    sig = 128.0
    f = tanh_profile(sig)
    e1 = volume_integrals(fields_nat, f, 1.0, sig, npan=80, nper=24, nu=96)["E_eul"] / (-sig**3 / 180)
    e2 = volume_integrals(fields_nat, f, 1e-4, sig, npan=80, nper=24, nu=96)["E_minus"] / (
        (math.pi / (24 if MUT == "thin" else 24 * math.sqrt(3))) * 1e-4 * sig**2)
    print(f"  thin wall, sigma R = {sig:g}: E_eul / (-v^2 sigma^3 / 180) = {e1:.4f}, "
          f"E_- / (pi v sigma^2 / (24 sqrt 3)) at v = 1e-4 = {e2:.4f}")
    check(C_THIN, abs(e1 - 1) < 0.02 and abs(e2 - 1) < 0.02)


# ------------------------------------------------------------------------------------------------ identity
def alc_fields(f, v):
    """#58's K, mu and J for Alcubierre's drive on a Cartesian grid (a node at r = 0 is allowed)."""
    def fields(X):
        r = np.sqrt(np.sum(X**2, axis=0))
        rs = np.where(r == 0, 1.0, r)
        F0, F1, F2 = f(rs)
        n = X / rs
        F1r = np.where(r == 0, f(np.array([1e-8]))[2][0], F1 / rs)
        Kij = np.zeros((3, 3) + r.shape)
        for i in range(3):
            Kij[0, i] += -(v / 2) * F1 * n[i]
            Kij[i, 0] += -(v / 2) * F1 * n[i]
        trK = Kij[0, 0] + Kij[1, 1] + Kij[2, 2]
        mu = (trK**2 - np.sum(Kij**2, axis=(0, 1))) / (16 * np.pi)
        hess = [[F2 * n[i] * n[j] + F1r * ((1 if i == j else 0) - n[i] * n[j]) for j in range(3)] for i in range(3)]
        lap = hess[0][0] + hess[1][1] + hess[2][2]
        J = np.stack([(v / (16 * np.pi)) * (hess[0][i] - (lap if i == 0 else 0)) for i in range(3)])
        return Kij, mu, J
    return fields


def nat_fields(f, v):
    """K from the symbolic 3+1 route, mu and J from the closed forms (the two share no derivation)."""
    S = get_sym()

    def fields(X):
        r = np.sqrt(np.sum(X**2, axis=0))
        F = f(r)
        a = (r, F[1] / 2, F[2] / 2, F[3] / 2, v)
        J = ev(S["J_f"], X, *a)
        if MUT == "Jgrid":
            J = J * (1 + 0.1 * r)
        return ev(S["K_f"], X, *a), ev(S["mu_f"], X, *a), J
    return fields


def solve_spinor(fields, L, N, eps0s, shift=0.0, tol=1e-11, maxit=400):
    """#58's Witten solve (Picard iteration of the FFT inverse of D), with K, mu and J as an argument."""
    g0, g = gammas()
    xg = -L + 2 * L / N * (np.arange(N) + shift)
    X = np.stack(np.meshgrid(xg, xg, xg, indexing="ij"))
    Kij, mu, J = fields(X)
    trK = Kij[0, 0] + Kij[1, 1] + Kij[2, 2]
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=2 * L / N)
    k = np.stack(np.meshgrid(k1, k1, k1, indexing="ij"))
    k2 = np.sum(k**2, axis=0)
    k2[0, 0, 0] = 1.0

    def apply(Mat, sp_):
        return np.einsum("ab,b...->a...", Mat, sp_)

    def deriv(sp_):
        h = np.fft.fftn(sp_, axes=(1, 2, 3))
        return [np.fft.ifftn(1j * k[j] * h, axes=(1, 2, 3)) for j in range(3)]

    def Dinv(src):
        h = np.fft.fftn(src, axes=(1, 2, 3))
        out = sum(apply(-1j * g[j], h * (k[j] / k2)) for j in range(3))
        out[:, 0, 0, 0] = 0
        return np.fft.ifftn(out, axes=(1, 2, 3))

    sB = -1.0 if MUT == "port" else 1.0
    sM = -1.0 if MUT == "musign" else 1.0
    B = lambda s_: sB * 0.5 * trK * apply(g0, s_)  # noqa: E731  D~ = D + (1/2) tr K g0 on the flat slice
    dV = (2 * L / N) ** 3
    Jn = np.sqrt(np.sum(J**2, axis=0))
    m_minus, s = minus_density(mu, Jn)
    Jm = -s * J / np.where(Jn == 0, 1, Jn)
    Krms = np.sqrt(np.mean(np.sum(Kij**2, axis=(0, 1))))
    common = {"E_minus": np.sum(m_minus) * dV, "E_eul": np.sum(mu) * dV, "intJ": np.sum(J, axis=(1, 2, 3)) * dV,
              "absJ": np.sum(Jn) * dV, "Krms": Krms, "trK_max": float(np.max(np.abs(trK)))}
    outs = []
    for eps0 in eps0s:
        e0 = np.asarray(eps0, dtype=complex)[:, None, None, None] * np.ones((1,) + mu.shape)
        psi = np.zeros_like(e0)
        for it in range(maxit):
            new = Dinv(-B(e0 + psi))
            dpsi = np.sqrt(np.mean(np.abs(new - psi) ** 2))
            psi = new
            if dpsi < tol:
                break
        eps = e0 + psi
        de = deriv(psi)
        resid = sum(apply(g[j], de[j]) for j in range(3)) + B(eps)
        b0 = np.sqrt(np.mean(np.abs(B(e0)) ** 2))
        A = [de[i] + 0.5 * sum(Kij[i, j] * apply(g[j] @ g0, eps) for j in range(3)) for i in range(3)]
        gpt = sum(np.sum(np.abs(a) ** 2, axis=0) for a in A) / (2 * np.pi if MUT == "norm" else 4 * np.pi)
        Gterm = np.sum(gpt) * dV
        e2 = np.sum(np.abs(eps) ** 2, axis=0)
        xi = np.stack([np.real(np.einsum("a...,ab,b...->...", eps.conj(), g0 @ g[j], eps)) for j in range(3)])
        T = sM * mu * e2 - np.sum(J * xi, axis=0)  # T(xi, n) = mu |eps|^2 - J.xi
        Tm = m_minus * e2 - np.sum(Jm * xi, axis=0)
        outs.append({"it": it + 1, "res": np.sqrt(np.mean(np.abs(resid) ** 2)) / b0 if b0 > 0 else float("nan"),
                     "res0": b0 / Krms, "dEps": np.sqrt(np.mean(np.abs(psi) ** 2)),
                     "alg": float(np.max(np.abs(gpt + sM * mu * e2)) / np.max(np.abs(mu))),
                     "G": Gterm, "intT": np.sum(T) * dV, "intTm": np.sum(Tm) * dV, "intTp": np.sum(T + Tm) * dV,
                     "minTp": float((T + Tm).min()), "xi": np.mean(xi, axis=(1, 2, 3)),
                     "eps2_min": float(e2.min()), "eps2_max": float(e2.max()), **common})
    return outs


E1 = [1, 0, 0, 0]
_a = np.array([1, 1]) / np.sqrt(2)
EXPLUS = np.concatenate([_a, 1j * _a]) / np.sqrt(2)
EXMINUS = np.concatenate([_a, -1j * _a]) / np.sqrt(2)


def part_identity():
    print("== identity: the Witten solve with #58's solver, ported with K as an argument")
    sig, v, L, N = 2.0, 0.5, 8.0, 64
    f = AB.tanh_profile(sig)
    ref = AB.solve_spinor(f, v, L, N, E1)
    got = solve_spinor(alc_fields(f, v), L, N, [E1])[0]
    rel = {k: abs(got[k] - ref[k]) / abs(ref[k]) for k in ("G", "intT", "intTm", "intTp", "E_minus", "E_eul")}
    print(f"  control, Alcubierre (sigma R = {sig:g}, v = {v}, L = {L:g}, N = {N}): #58's G = {ref['G']:.6e}, "
          f"ported {got['G']:.6e}, iterations {ref['it']} and {got['it']}, max relative difference "
          f"{max(rel.values()):.1e}")
    check(C_PORT, max(rel.values()) < 1e-9 and ref["it"] == got["it"])

    cfgs = [(2.0, 0.5, 8.0, (128,))] if QUICK else [
        (2.0, 0.5, 8.0, (64, 96, 128)), (2.0, 1.0, 8.0, (96, 128)), (4.0, 0.5, 5.0, (96, 128)),
        (4.0, 2.0, 5.0, (128,)), (8.0, 0.1, 3.0, (128,))]
    print("  Natario data, eps0 = e1 (xi = 0) and the null spinors with xi = -+ e_x; the grid is offset by half a cell")
    print("  sigma R |  v   |  L | N   | eps0 | iters | |D~eps0|/|K| | |eps-eps0| | pointwise |      G      | "
          "int T(xi,n) | |G+int T|/G | |G+int T|/terms | |int J|/int|J| | int T-(xi,n) | int T+(xi,n) | min T+")
    ok_tp = True
    worst = {"eps": 0.0, "res": 0.0, "alg": 0.0, "wit": 0.0, "intj": 0.0, "eul": 0.0, "em": 0.0}
    for sig, v, L, Ns in cfgs:
        f = tanh_profile(sig)
        q = volume_integrals(fields_nat, f, v, sig, npan=80, nper=24, nu=96)
        for N in Ns:
            outs = solve_spinor(nat_fields(f, v), L, N, [E1, EXPLUS, EXMINUS], shift=0.5)
            for name, o in zip(("e1", "e+", "e-"), outs):
                terms = o["G"] + np.linalg.norm(o["xi"]) * o["absJ"]
                rel, relt = abs(o["G"] + o["intT"]) / o["G"], abs(o["G"] + o["intT"]) / terms
                intj = np.max(np.abs(o["intJ"])) / o["absJ"]
                print(f"  {sig:7.1f} | {v:4.2f} | {L:2.0f} | {N:3d} | {name:4s} | {o['it']:5d} | {o['res0']:11.2e} | "
                      f"{o['dEps']:10.2e} | {o['alg']:9.2e} | {o['G']:11.5e} | {o['intT']:12.5e} | {rel:11.2e} | "
                      f"{relt:14.2e} | {intj:14.2e} | {o['intTm']:12.5e} | {o['intTp']:12.5e} | {o['minTp']:.1e}")
                if N == Ns[-1]:
                    worst["eps"] = max(worst["eps"], o["dEps"] / o["Krms"])
                    worst["res"] = max(worst["res"], o["res0"])
                    worst["alg"] = max(worst["alg"], o["alg"])
                    worst["wit"] = max(worst["wit"], relt)
                    worst["intj"] = max(worst["intj"], intj)
                    ok_tp = ok_tp and o["minTp"] > -1e-12 * o["E_minus"]
            o = outs[0]
            if N == Ns[-1]:
                worst["eul"] = max(worst["eul"], abs(o["E_eul"] - q["E_eul"]) / abs(q["E_eul"]))
                worst["em"] = max(worst["em"], abs(o["E_minus"] - q["E_minus"]) / q["E_minus"])
            print(f"          xi(e+) = {np.array2string(outs[1]['xi'], precision=4)}, xi(e-) = "
                  f"{np.array2string(outs[2]['xi'], precision=4)}, max|tr K| = {o['trK_max']:.1e}, E_- grid "
                  f"{o['E_minus']:.6e} vs quadrature {q['E_minus']:.6e}, E_Eul grid {o['E_eul']:.8e} vs "
                  f"{q['E_eul']:.8e}, G/|E_Eul| = {outs[0]['G'] / abs(q['E_eul']):.10f}")
    check(C_EPS0, worst["eps"] < 1e-12 and worst["res"] < 1e-12,
          f"worst |eps - eps0|/|K| {worst['eps']:.1e}, |D~ eps0|/|K| {worst['res']:.1e}")
    check(C_ALG, worst["alg"] < 1e-12, f"worst {worst['alg']:.1e} of max |mu|")
    check(C_WITTEN, worst["wit"] < 1e-8,
          f"worst |G + int T| / (G + |xi0| int|J|) {worst['wit']:.1e} at the finest grid")
    check(C_INTJG, worst["intj"] < 1e-8, f"worst {worst['intj']:.1e}")
    check(C_TPLUS, ok_tp)
    check(C_GEUL, worst["eul"] < 1e-6, f"worst {worst['eul']:.1e}")
    check(C_GEM, worst["em"] < 1e-3, f"worst {worst['em']:.1e}")


# ------------------------------------------------------------------------------------------------ mutations
# mutation -> [(part, checks that must fail)]
MUT_EXPECT = {
    "div": [("sym", [C_DIV, C_NATX, C_MU3, C_J3, C_MU4, C_J4]), ("identity", [C_EPS0, C_ALG, C_WITTEN])],
    "beta2": [("sym", [C_MU3, C_J3])],
    "Jsign": [("sym", [C_J3])],
    "nusign": [("sym", [C_MU4, C_J4])],
    "natX": [("sym", [C_NATX])],
    "natrho": [("sym", [C_NATRHO])],
    "rw1": [("sym", [C_NUM])],
    "mu6": [("sym", [C_NUM]), ("identity", [C_GEUL])],
    "prof3": [("sym", [C_PROF])],
    "gamdup": [("sym", [C_GAM]), ("identity", [C_ALG])],
    "adm_tail": [("adm", [C_ADM])],
    "adm_m": [("adm", [C_SCHW])],
    "adm_p": [("adm", [C_BY])],
    "split_old": [("split", [C_SPL1, C_SPL2, C_SPLA]), ("identity", [C_TPLUS])],
    "quad": [("energy", [C_GAUSS, C_REPRO, C_INTJ])],
    "coarse": [("energy", [C_CONV])],
    "thin": [("energy", [C_THIN])],
    "port": [("identity", [C_PORT])],
    "musign": [("identity", [C_ALG, C_WITTEN])],
    "norm": [("identity", [C_PORT, C_ALG, C_WITTEN])],
    "Jgrid": [("identity", [C_INTJG, C_GEM])],
}
FLAGS = {"identity": ["--quick"], "energy": ["--quick"]}


def run_sub(part, mut):
    cmd = [sys.executable, os.path.abspath(__file__), f"--part={part}"] + FLAGS.get(part, []) + (
        [f"--mutate={mut}"] if mut else [])
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, [ln for ln in p.stdout.splitlines() if ln.startswith("  [FAIL] ")]


def part_mutations():
    print("== mutations: each must make its named checks fail, and the unmutated quick runs must pass")
    allc = {val for key, val in globals().items() if key.startswith("C_")}
    covered = {n for lst in MUT_EXPECT.values() for _, names in lst for n in names}
    check("every check is named by a mutation", allc <= covered,
          f"{len(allc)} checks, unnamed {sorted(allc - covered)}")
    jobs = [("", part, []) for part in ("identity", "energy")]
    jobs += [(m, part, names) for m, lst in MUT_EXPECT.items() for part, names in lst]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        res = list(ex.map(lambda j: run_sub(j[1], j[0]), jobs))
    for (m, part, names), (code, fails) in zip(jobs, res):
        if not m:
            check(f"unmutated {part} --quick passes", code == 0 and not fails, f"exit {code}, {len(fails)} failed")
            continue
        missing = [n for n in names if not any(ln.startswith(f"  [FAIL] {n}") for ln in fails)]
        check(f"--mutate={m} ({part}) fails its {len(names)} named check(s)", not missing and code == 1,
              f"exit {code}, {len(fails)} failed" + (f", NOT failed: {missing}" if missing else ""))


PARTS = {"sym": part_sym, "adm": part_adm, "split": part_split, "energy": part_energy, "identity": part_identity,
         "mutations": part_mutations}
if MUT and MUT not in MUT_EXPECT:
    sys.exit(f"unknown mutation {MUT}")
for name, fn in PARTS.items():
    if args.part in ("all", name):
        fn()
print("FAILED: " + ", ".join(FAILED) if FAILED else "all checks passed")
sys.exit(1 if FAILED else 0)
