#!/usr/bin/env python3
"""Issue #58: the hypotheses of phase 4's signed-source mass bounds on Alcubierre initial data.

The t = 0 slice of ds^2 = -dt^2 + (dx - v f(r_s) dt)^2 + dy^2 + dz^2 is flat (g_ij = delta_ij), with lapse 1, shift
beta_x = -v f and K_ij = (d_i beta_j + d_j beta_i)/2. G = c = 1, bubble radius R = 1.

    python notes/warp/alcubierre_bounds.py [--part=all|sym|adm|split|energy|spinor] [--mutate=NAME]

Parts:
  sym     mu and J from the 3+1 constraints and from G_ab of the 4-metric, against the closed forms; the spherical
          forms used below; K_ij g^i g^j = tr K for the Witten operator.
  adm     E and P as surface integrals (zero here); controls with known E (Schwarzschild) and P (Bowen-York).
  split   the pointwise minimal split T = T+ - T- (mu+- >= |J+-|) against a numerical minimisation.
  energy  E- (= W-, psi = 1), the Eulerian energy and the DEC deficit for the tanh profile, by quadrature.
  spinor  the existence constant kappa_K = ||tr K||_{L^3} / (2 sqrt(S)), and the Witten spinor by FFT with its identity.

Mutations (each must make its check fail): beta2, adm_m, adm_p, split_old, quad, Jsign, norm.
Exits 1 when any check fails.
"""
import argparse
import math
import sys

import numpy as np
import sympy as sp

ap = argparse.ArgumentParser()
ap.add_argument("--part", default="all")
ap.add_argument("--mutate", default="")
args = ap.parse_args()
MUT = args.mutate
FAILED = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILED.append(name)


# ------------------------------------------------------------------------------------------------ sym
def part_sym():
    print("== sym: mu and J by two routes")
    t, x, y, z = sp.symbols("t x y z", real=True)
    v, r = sp.symbols("v r", positive=True)
    F = sp.symbols("F0:7", real=True)  # F[k] = f^(k)(r_s)
    xs = x - v * t
    dr = {t: -v * xs / r, x: xs / r, y: y / r, z: z / r}

    def d(e, c):
        # r and F[k] are symbols standing for r_s and f^(k)(r_s); apply the chain rule by hand
        out = sp.diff(e, c) + sp.diff(e, r) * dr[c]
        for k in range(6):
            out += sp.diff(e, F[k]) * F[k + 1] * dr[c]
        return out

    def is_zero(e):
        e = sp.together(sp.expand(e.subs(t, 0)))
        num = sp.expand(sp.numer(e))
        rem = sp.Poly(num, z).rem(sp.Poly(z**2 - (r**2 - x**2 - y**2), z))
        return sp.expand(rem.as_expr()) == 0

    # closed forms at t = 0 (the issue's)
    n = (x / r, y / r, z / r)
    X3 = (x, y, z)
    hess = [[F[2] * n[i] * n[j] + F[1] / r * ((1 if i == j else 0) - n[i] * n[j]) for j in range(3)] for i in range(3)]
    lap = F[2] + 2 * F[1] / r
    mu_c = -(v**2 / (32 * sp.pi)) * (y**2 + z**2) / r**2 * F[1] ** 2
    J_c = [(v / (16 * sp.pi)) * (hess[0][i] - (lap if i == 0 else 0)) for i in range(3)]

    # route 1: 3+1 constraints on the flat slice
    bscale = 2 if MUT == "beta2" else 1
    beta = [-bscale * v * F[0], 0, 0]
    K = [[(d(beta[j], X3[i]) + d(beta[i], X3[j])) / 2 for j in range(3)] for i in range(3)]
    trK = sum(K[i][i] for i in range(3))
    mu3 = (trK**2 - sum(K[i][j] ** 2 for i in range(3) for j in range(3))) / (16 * sp.pi)
    J3 = [(sum(d(K[j][i], X3[j]) for j in range(3)) - d(trK, X3[i])) / (8 * sp.pi) for i in range(3)]
    check("3+1 mu = closed form", is_zero(mu3 - mu_c))
    check("3+1 J = closed form", all(is_zero(J3[i] - J_c[i]) for i in range(3)))

    # route 2: Einstein tensor of the 4-metric, Eulerian projections
    X4 = (t, x, y, z)
    g = sp.Matrix([[-1 + v**2 * F[0] ** 2, -v * F[0], 0, 0], [-v * F[0], 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    gi = sp.simplify(g.inv())
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
    nu = [1, v * F[0], 0, 0]
    mu4 = sum(G[a, b] * nu[a] * nu[b] for a in range(4) for b in range(4)) / (8 * sp.pi)
    J4 = [-sum(nu[a] * G[a, i + 1] for a in range(4)) / (8 * sp.pi) for i in range(3)]
    check("4D mu = closed form", is_zero(mu4 - mu_c))
    check("4D J = closed form", all(is_zero(J4[i] - J_c[i]) for i in range(3)))

    # spherical forms used by the numerics (u = cos theta, s = sin theta, about the x axis)
    c_, s_ = x / r, sp.sqrt(y**2 + z**2) / r
    Jx_sph = -(v / (16 * sp.pi)) * (F[2] * s_**2 + F[1] / r * (1 + c_**2))
    Jy_sph = (v / (16 * sp.pi)) * (F[2] - F[1] / r) * x * y / r**2
    check("spherical J_x, J_y", is_zero(Jx_sph - J_c[0]) and is_zero(Jy_sph - J_c[1]))
    check("tr K = -v f' cos(theta)", is_zero(trK + bscale * v * F[1] * x / r) if MUT != "beta2" else True)

    # on the x axis inside the wall: mu = 0, J_x = -v f'/(8 pi r)
    on_axis = {y: 0, z: 0, r: x}
    print(f"  on the x axis: mu = {sp.simplify(mu_c.subs(on_axis))}, J_x = {sp.simplify(J_c[0].subs(on_axis))}")

    # Witten operator on the flat slice: K_ij g^i g^j = (tr K) 1 for symmetric K
    gam0, gam = gammas_sym()
    k = sp.symbols("k11 k12 k13 k22 k23 k33", real=True)
    Ks = sp.Matrix([[k[0], k[1], k[2]], [k[1], k[3], k[4]], [k[2], k[4], k[5]]])
    S = sum((Ks[i, j] * gam[i] * gam[j] for i in range(3) for j in range(3)), sp.zeros(4, 4))
    check("K_ij g^i g^j = tr K", (S - Ks.trace() * sp.eye(4)).applyfunc(sp.expand) == sp.zeros(4, 4))
    cl = all((gam[i] * gam[j] + gam[j] * gam[i] - 2 * (1 if i == j else 0) * sp.eye(4)) == sp.zeros(4, 4)
             for i in range(3) for j in range(3))
    cl = cl and gam0 * gam0 == -sp.eye(4) and all(gam0 * gm + gm * gam0 == sp.zeros(4, 4) for gm in gam)
    cl = cl and gam0.H == -gam0 and all(gm.H == gm for gm in gam)
    check("gamma algebra: g0^2 = -1 anti-Hermitian, g^i Hermitian, anticommuting", cl)


def gammas_sym():
    I2, Z2 = sp.eye(2), sp.zeros(2, 2)
    sig = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
    gam0 = sp.I * sp.diag(I2, -I2)
    gam = [sp.Matrix(sp.BlockMatrix([[Z2, s], [s, Z2]])) for s in sig]
    return gam0, gam


def gammas_np():
    g0, g = gammas_sym()
    return np.array(g0.tolist(), dtype=complex), [np.array(m.tolist(), dtype=complex) for m in g]


# ------------------------------------------------------------------------------------------------ profiles
def tanh_profile(sig, R=1.0):
    T = math.tanh(sig * R)

    def f(r):
        ta, tb = np.tanh(sig * (r + R)), np.tanh(sig * (r - R))
        sa, sb = 1 - ta**2, 1 - tb**2
        return (ta - tb) / (2 * T), sig * (sa - sb) / (2 * T), sig**2 * (-2 * sa * ta + 2 * sb * tb) / (2 * T)
    return f


def gauss_profile(r):
    e = np.exp(-r**2)
    return e, -2 * r * e, (4 * r**2 - 2) * e


def fields_sph(f, v, r, u):
    """mu, J_x, J_perp (cylindrical) at radius r and u = cos(theta)."""
    F0, F1, F2 = f(r)
    s2 = 1 - u**2
    mu = -(v**2 / (32 * np.pi)) * s2 * F1**2
    Jx = -(v / (16 * np.pi)) * (F2 * s2 + F1 / r * (1 + u**2))
    Jp = (v / (16 * np.pi)) * (F2 - F1 / r) * u * np.sqrt(s2)
    return mu, Jx, Jp


def minus_density(mu, Jn):
    """Pointwise minimal mu_- over splits with mu_+- >= |J_+-|, and the |J_-| that attains it."""
    if MUT == "split_old":
        m = np.maximum(0.0, (Jn - mu) / 2)
        return m, m
    m = np.maximum.reduce([np.zeros_like(mu), -mu, (Jn - mu) / 2])
    s = np.where(Jn >= -mu, m, Jn)  # J_- = -s J/|J|
    return m, s


# ------------------------------------------------------------------------------------------------ adm
def sphere_quad(nu=48, nphi=96):
    u, wu = np.polynomial.legendre.leggauss(nu)
    ph = np.arange(nphi) * 2 * np.pi / nphi
    U, PH = np.meshgrid(u, ph, indexing="ij")
    W = np.outer(wu, np.full(nphi, 2 * np.pi / nphi))
    S = np.sqrt(1 - U**2)
    return np.stack([U, S * np.cos(PH), S * np.sin(PH)]), W


def adm_E(gfun, R, h=1e-4):
    n, W = sphere_quad()
    X = R * n
    # E = (1/16 pi) oint (d_j g_ij - d_i g_jj) n_i dS
    tot = np.zeros_like(W)
    for i in range(3):
        for j in range(3):
            ej = np.eye(3)[j][:, None, None] * h
            ei = np.eye(3)[i][:, None, None] * h
            dj_gij = (gfun(X + ej)[i, j] - gfun(X - ej)[i, j]) / (2 * h)
            di_gjj = (gfun(X + ei)[j, j] - gfun(X - ei)[j, j]) / (2 * h)
            tot += (dj_gij - di_gjj) * n[i]
    return float(np.sum(tot * W) * R**2 / (16 * np.pi))


def adm_P(Kfun, R):
    n, W = sphere_quad()
    Kij = Kfun(R * n)
    tr = Kij[0, 0] + Kij[1, 1] + Kij[2, 2]
    return np.array([np.sum(sum((Kij[i, j] - (tr if i == j else 0)) * n[j] for j in range(3)) * W) * R**2 / (8 * np.pi)
                     for i in range(3)])


def K_alcubierre(f, v):
    def K(X):
        r = np.sqrt(np.sum(X**2, axis=0))
        F0, F1, F2 = f(r)
        n = X / r
        out = np.zeros((3, 3) + r.shape)
        for i in range(3):
            out[0, i] += -(v / 2) * F1 * n[i]
            out[i, 0] += -(v / 2) * F1 * n[i]
        return out
    return K


def part_adm():
    print("== adm: E and P as surface integrals")
    f = tanh_profile(4.0)
    flat = lambda X: np.broadcast_to(np.eye(3)[:, :, None, None], (3, 3) + X.shape[1:])  # noqa: E731
    for R in (2.0, 4.0, 8.0):
        E = adm_E(flat, R)
        P = adm_P(K_alcubierre(f, 1.0), R)
        print(f"  Alcubierre (sigma R = 4, v = 1), R = {R}: E = {E:.3e}, P = {np.array2string(P, precision=3)}")
    check("Alcubierre E = 0 and |P| < 1e-12 at R = 8", abs(E) < 1e-12 and np.max(np.abs(P)) < 1e-12)
    m = 0.7

    def schw(X):
        r = np.sqrt(np.sum(X**2, axis=0))
        psi = 1 + m / (2 * r) * (1.1 if MUT == "adm_m" else 1.0)
        return np.eye(3)[:, :, None, None] * psi**4
    for R in (5.0, 20.0):
        E = adm_E(schw, R)
        want = m * (1 + m / (2 * R)) ** 3
        print(f"  control, Schwarzschild m = {m}, R = {R}: E = {E:.8f}, exact m psi(R)^3 = {want:.8f}")
    check("control: Schwarzschild E = m psi(R)^3", abs(E - want) < 1e-6)
    Pv = np.array([0.3, -0.2, 0.5])

    def by(X):
        r = np.sqrt(np.sum(X**2, axis=0))
        n = X / r
        Pn = sum(Pv[i] * n[i] for i in range(3))
        sc = 1.1 if MUT == "adm_p" else 1.0
        return np.array([[sc * 1.5 / r**2 * (Pv[i] * n[j] + Pv[j] * n[i] - ((1 if i == j else 0) - n[i] * n[j]) * Pn)
                          for j in range(3)] for i in range(3)])
    P = adm_P(by, 3.0)
    print(f"  control, Bowen-York P = {Pv}: surface integral {np.array2string(P, precision=10)}")
    check("control: Bowen-York P returned", np.max(np.abs(P - Pv)) < 1e-10)


# ------------------------------------------------------------------------------------------------ split
def part_split():
    print("== split: minimal mu_- against a numerical minimisation")
    from scipy.optimize import minimize
    rng = np.random.default_rng(58)
    worst, worst_sum, n_negmu = 0.0, 0.0, 0
    for _ in range(300):
        mu = rng.normal() * rng.choice([0.1, 1, 3])
        J = rng.normal(size=3) * rng.choice([0.1, 1, 3])
        Jn = np.linalg.norm(J)
        n_negmu += (-mu > Jn)
        obj = lambda Jm: max(np.linalg.norm(Jm), np.linalg.norm(J + Jm) - mu)  # noqa: E731
        best = min((minimize(obj, x0, method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-12, "maxiter": 20000})
                    for x0 in [np.zeros(3), -J / 2, -J, rng.normal(size=3)]), key=lambda o: o.fun)
        m, s = minus_density(np.array(mu), np.array(Jn))
        worst = max(worst, abs(best.fun - float(m)) / (1 + abs(best.fun)))
        # the split it names must be admissible: mu_- >= |J_-|, mu_+ = mu + mu_- >= |J + J_-|
        Jm = -float(s) * J / Jn
        adm = float(m) >= np.linalg.norm(Jm) - 1e-12 and mu + float(m) >= np.linalg.norm(J + Jm) - 1e-12
        worst_sum = max(worst_sum, 0 if adm else 1)
    print(f"  300 random (mu, J), {n_negmu} with -mu > |J|: max relative difference {worst:.2e}")
    check("minimal mu_- = max(0, -mu, (|J| - mu)/2)", worst < 1e-6)
    check("the named split is admissible", worst_sum == 0)


# ------------------------------------------------------------------------------------------------ energy
def radial_nodes(sig, R=1.0, npan=40, nper=16, outer=40.0):
    a, b = max(0.0, R - 14 / sig), R + 14 / sig
    edges = list(np.linspace(0, a, 6)[:-1]) if a > 0 else []
    edges += list(np.linspace(a, b, npan + 1)) + list(np.linspace(b, R + outer / sig, 9)[1:])
    x, w = np.polynomial.legendre.leggauss(nper)
    rs, ws = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        rs.append((hi - lo) / 2 * x + (hi + lo) / 2)
        ws.append((hi - lo) / 2 * w)
    return np.concatenate(rs), np.concatenate(ws)


def volume_integrals(f, v, sig, npan=40, nper=16, nu=64):
    r, wr = radial_nodes(sig, npan=npan, nper=nper)
    u, wu = np.polynomial.legendre.leggauss(nu)
    if MUT == "quad":
        wu = wu * (1 + 0.01 * u**2)
    Rr, U = np.meshgrid(r, u, indexing="ij")
    W = 2 * np.pi * np.outer(wr * r**2, wu)
    mu, Jx, Jp = fields_sph(f, v, Rr, U)
    Jn = np.hypot(Jx, Jp)
    m, _ = minus_density(mu, Jn)
    return {"E_eul": np.sum(mu * W), "E_minus": np.sum(m * W), "deficit": np.sum(np.maximum(0, Jn - mu) * W),
            "int_absJ": np.sum(Jn * W)}


def part_energy():
    print("== energy: quadrature control, then the tanh profile")
    rr = sp.symbols("r", positive=True)
    vv = 0.8
    exact = float(-(vv**2 / 12) * sp.integrate((sp.diff(sp.exp(-rr**2), rr)) ** 2 * rr**2, (rr, 0, sp.oo)))
    got = volume_integrals(gauss_profile, vv, 1.0)["E_eul"]
    print(f"  control, f = exp(-r^2), v = {vv}: quadrature {got:.12f}, exact -(v^2/12) int f'^2 r^2 dr = {exact:.12f}")
    check("control: Gaussian Eulerian energy", abs(got - exact) < 1e-10 * abs(exact))
    # profile derivatives against finite differences
    f = tanh_profile(4.0)
    r0, h = np.array([0.3, 0.9, 1.1, 1.7]), 1e-5
    fd1 = (f(r0 + h)[0] - f(r0 - h)[0]) / (2 * h)
    fd2 = (f(r0 + h)[1] - f(r0 - h)[1]) / (2 * h)
    check("tanh profile f', f'' against finite differences",
          np.max(np.abs(fd1 - f(r0)[1])) < 1e-6 and np.max(np.abs(fd2 - f(r0)[2])) < 1e-5)
    print("  sigma R |   v  | E_Eulerian   | E_- (= W-)   | int(|J|-mu)+ | E_-/(v R^2 sigma) | rel. change, 2x nodes")
    rows = []
    for sig in (2.0, 4.0, 8.0, 16.0):
        f = tanh_profile(sig)
        for v in (0.1, 0.5, 1.0, 2.0, 5.0):
            a = volume_integrals(f, v, sig)
            b = volume_integrals(f, v, sig, npan=80, nper=24, nu=96)
            ch = max(abs(a[k] - b[k]) / abs(b[k]) for k in ("E_eul", "E_minus", "deficit"))
            rows.append(ch)
            print(f"  {sig:7.1f} | {v:4.1f} | {b['E_eul']:12.5e} | {b['E_minus']:12.5e} | {b['deficit']:12.5e} | "
                  f"{b['E_minus'] / (v * sig):17.5e} | {ch:.1e}")
    check("E_- converged (relative change < 1e-4 on doubling)", max(rows) < 1e-4)


# ------------------------------------------------------------------------------------------------ spinor
S_SOB = 3 * (np.pi / 2) ** (4 / 3)  # sharp Sobolev constant on R^3: int |grad h|^2 >= S ||h||_6^2


def kappa_K(f, v, sig):
    r, wr = radial_nodes(sig, npan=80, nper=24)
    F0, F1, F2 = f(r)
    L3 = (v**3 * np.pi * np.sum(wr * r**2 * np.abs(F1) ** 3)) ** (1 / 3)  # int |cos|^3 sin dtheta = 1/2
    return L3 / (2 * math.sqrt(S_SOB))


def solve_spinor(f, v, L, N, eps0, tol=1e-11, maxit=400):
    g0, g = gammas_np()
    xg = -L + 2 * L / N * np.arange(N)
    X = np.stack(np.meshgrid(xg, xg, xg, indexing="ij"))
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
    k1 = 2 * np.pi * np.fft.fftfreq(N, d=2 * L / N)
    k = np.stack(np.meshgrid(k1, k1, k1, indexing="ij"))
    k2 = np.sum(k**2, axis=0)
    k2[0, 0, 0] = 1.0

    def apply(Mat, sp_):  # (4,4) @ (4,...)
        return np.einsum("ab,b...->a...", Mat, sp_)

    def deriv(sp_):
        h = np.fft.fftn(sp_, axes=(1, 2, 3))
        return [np.fft.ifftn(1j * k[j] * h, axes=(1, 2, 3)) for j in range(3)]

    def Dinv(src):
        h = np.fft.fftn(src, axes=(1, 2, 3))
        out = sum(apply(-1j * g[j], h * (k[j] / k2)) for j in range(3))
        out[:, 0, 0, 0] = 0
        return np.fft.ifftn(out, axes=(1, 2, 3))

    e0 = np.asarray(eps0, dtype=complex)[:, None, None, None] * np.ones((1,) + r.shape)
    psi = np.zeros_like(e0)
    B = lambda s_: 0.5 * trK * apply(g0, s_)  # noqa: E731  D~ = D + (1/2) tr K g0 on the flat slice
    for it in range(maxit):
        new = Dinv(-B(e0 + psi))
        dpsi = np.sqrt(np.mean(np.abs(new - psi) ** 2))
        psi = new
        if dpsi < tol:
            break
    eps = e0 + psi
    de = deriv(psi)
    resid = sum(apply(g[j], de[j]) for j in range(3)) + B(eps)
    res_rel = np.sqrt(np.mean(np.abs(resid) ** 2)) / np.sqrt(np.mean(np.abs(B(e0)) ** 2))
    dV = (2 * L / N) ** 3
    A = [de[i] + 0.5 * sum(Kij[i, j] * apply(g[j] @ g0, eps) for j in range(3)) for i in range(3)]
    grad2 = sum(np.sum(np.abs(a) ** 2, axis=0) for a in A)
    Gterm = np.sum(grad2) * dV / (2 * np.pi if MUT == "norm" else 4 * np.pi)
    e2 = np.sum(np.abs(eps) ** 2, axis=0)
    xi = np.stack([np.real(np.einsum("a...,ab,b...->...", eps.conj(), g0 @ g[j], eps)) for j in range(3)])
    sJ = 1.0 if MUT == "Jsign" else -1.0  # T(xi, n) = mu |eps|^2 - J.xi in these conventions
    T = mu * e2 + sJ * np.sum(J * xi, axis=0)
    Jn = np.sqrt(np.sum(J**2, axis=0))
    m_minus, s = minus_density(mu, Jn)
    Jm = -s * J / np.where(Jn == 0, 1, Jn)
    Tm = m_minus * e2 + sJ * np.sum(Jm * xi, axis=0)
    out = {"it": it + 1, "res": res_rel, "G": Gterm, "intT": np.sum(T) * dV, "intTm": np.sum(Tm) * dV,
           "intTp": np.sum(T + Tm) * dV, "E_minus": np.sum(m_minus) * dV, "E_eul": np.sum(mu) * dV,
           "eps2_min": float(e2.min()), "eps2_max": float(e2.max())}
    return out


def part_spinor():
    print(f"== spinor: existence constant kappa_K = ||tr K||_L3 / (2 sqrt(S)), S = 3 (pi/2)^(4/3) = {S_SOB:.6f}")
    print("  sigma R | v_c (kappa_K = 1) | kappa_K at v = 0.5, 1, 2")
    for sig in (2.0, 4.0, 8.0, 16.0):
        f = tanh_profile(sig)
        k1 = kappa_K(f, 1.0, sig)
        print(f"  {sig:7.1f} | {1 / k1:17.5f} | " + ", ".join(f"{kappa_K(f, v, sig):.4f}" for v in (0.5, 1.0, 2.0)))
    print("  Witten spinor by FFT (eps0 = (1,0,0,0)); identity G + int T(xi, n) = E - P.xi0 = 0")
    print("  sigma R |  v  |  L  |  N  | kappa_K | iters | |D~eps|/|B eps0| |      G      |  int T(xi,n) | "
          "rel. resid | int T-(xi,n) | int T+(xi,n) |   E_-   | |eps|^2 range")
    worst = 0.0
    for sig, v, L, Ns in ((2.0, 0.5, 8.0, (64, 96, 128)), (2.0, 1.0, 8.0, (96, 128)), (4.0, 0.5, 5.0, (96, 128))):
        f = tanh_profile(sig)
        kap = kappa_K(f, v, sig)
        for N in Ns:
            o = solve_spinor(f, v, L, N, [1, 0, 0, 0])
            rel = abs(o["G"] + o["intT"]) / o["G"]
            if N == Ns[-1]:
                worst = max(worst, rel)
            print(f"  {sig:7.1f} | {v:3.1f} | {L:3.0f} | {N:3d} | {kap:7.4f} | {o['it']:5d} | {o['res']:16.2e} | "
                  f"{o['G']:11.5e} | {o['intT']:12.5e} | {rel:10.2e} | {o['intTm']:12.5e} | {o['intTp']:12.5e} | "
                  f"{o['E_minus']:.5f} | {o['eps2_min']:.4f}-{o['eps2_max']:.4f}")
    check("Witten identity G + int T(xi, n) = 0 at the finest grid (relative 1e-4)", worst < 1e-4,
          f"worst {worst:.2e}")


PARTS = {"sym": part_sym, "adm": part_adm, "split": part_split, "energy": part_energy, "spinor": part_spinor}
for name, fn in PARTS.items():
    if args.part in ("all", name):
        fn()
print("FAILED: " + ", ".join(FAILED) if FAILED else "all checks passed")
sys.exit(1 if FAILED else 0)
