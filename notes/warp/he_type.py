#!/usr/bin/env python3
"""Hawking-Ellis type of the Alcubierre Einstein tensor, in closed form (issue #64).

G^a_b and R^a_b differ by a multiple of the identity, so they have the same type. In the orthonormal frame
(E_0 = n, E_1 = d_x, e_rho, e_varphi) of phantom_source.py, R(e_varphi, .) = 0 (I1). So e_varphi is an eigenvector of
R^a_b with eigenvalue 0, and the type is that of the 3x3 mixed block M = eta R on (E_0, E_1, e_rho). Its characteristic
polynomial is the cubic lam^3 + c2 lam^2 + c1 lam + c0, with discriminant D. D > 0 gives Type I, D < 0 Type IV, and on
D = 0 the rank of M - lam I at the repeated root decides between Types I, II and III (`he_type_exact`).

    python3 notes/warp/he_type.py                 all checks, then every control (exit 1 if any check fails)
    python3 notes/warp/he_type.py --mutate=NAME   one control (its target checks must fail, exit 1)
    python3 notes/warp/he_type.py --no-controls   the checks only

Routes. (C) The closed form. phantom_source.sph_tetrad_ricci gives the frame Ricci tensor in comoving spherical
coordinates for a general w(r) = v (1 - f(r)), and `mixed_block`, `cubic` and `disc` turn its block into D, as a sympy
expression in (v, f, f', f'', r, cos theta) and, lambdified, at grid points. (N) The numerical eigenvalues. `einstein_N`
builds the Einstein tensor of the lab metric at t = 0 from g and its first and second derivatives, written out by the
chain rule, then numpy's general eigensolver labels G^a_b (mpmath's at 50 digits for a recheck). C and N share the
metric's definition, written separately in each, and the profile jets `tanh_jets`, which are checked against mpmath's
numerical derivative. The planted tensors go through the same `mixed_block`, `cubic`, `disc` and `he_type_exact`.
"""
import os
import subprocess
import sys
import time

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phantom_source as ps  # noqa: E402

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "alc-drop-w2": (["decouple", "structure", "closed_vs_numeric"],
                    "drops the shift's square from g_tt in route C (phantom_source's control)"),
    "num-drop-w2": (["numN", "closed_vs_numeric"], "drops the shift's square from g_tt in route N"),
    "disc-sign": (["planted"], "flips the sign of the cubic discriminant"),
    "no-rank": (["planted"], "labels every D = 0 tensor Type II, without the rank test"),
    "jets-f2": (["jets"], "scales f'' of the tanh profile by 1 + 1/1000"),
    "le-speed": (["le_fraction"], "evaluates Le's setup at v = 1 instead of 1/2"),
    "rodal-speed": (["rodal_cuts"], "evaluates Rodal's setup at v = 1/4 instead of 1"),
}

RESULTS = []
MUT = None


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# The cubic and its discriminant, for entries of any arithmetic type (sympy, rationals, numpy arrays, mpf)
# ---------------------------------------------------------------------------------------------------------------
def mixed_block(R):
    """M^a_b = eta^ac R_cb on (E_0, E_1, e_rho), from the lower-index 3x3 block R (nested lists)."""
    return [[-R[0][0], -R[0][1], -R[0][2]], [R[1][0], R[1][1], R[1][2]], [R[2][0], R[2][1], R[2][2]]]


def cubic(M):
    """(c2, c1, c0) of det(lam I - M) = lam^3 + c2 lam^2 + c1 lam + c0."""
    c2 = -(M[0][0] + M[1][1] + M[2][2])
    c1 = ((M[0][0] * M[1][1] - M[0][1] * M[1][0]) + (M[0][0] * M[2][2] - M[0][2] * M[2][0])
          + (M[1][1] * M[2][2] - M[1][2] * M[2][1]))
    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
           + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    return c2, c1, -det


def disc(c2, c1, c0):
    """Discriminant of the monic cubic: > 0 three distinct real roots, < 0 a complex pair."""
    D = 18 * c2 * c1 * c0 - 4 * c2**3 * c0 + c2**2 * c1**2 - 4 * c1**3 - 27 * c0**2
    return -D if MUT == "disc-sign" else D


def he_type_exact(R4):
    """Hawking-Ellis type of a 4x4 exact symmetric frame tensor R_ab whose axis 3 decouples."""
    assert all(R4[a, 3] == 0 for a in range(3)), "e_3 does not decouple"
    M = mixed_block([[R4[a, b] for b in range(3)] for a in range(3)])
    c2, c1, c0 = cubic(M)
    D = sp.nsimplify(disc(c2, c1, c0))
    if D > 0:
        return "I", D
    if D < 0:
        return "IV", D
    if MUT == "no-rank":
        return "II", D
    lam = sp.Symbol("lam")
    roots = sp.roots(sp.Poly(lam**3 + c2 * lam**2 + c1 * lam + c0, lam))
    lr, mult = max(roots.items(), key=lambda kv: kv[1])
    rk = (sp.Matrix(M) - lr * sp.eye(3)).rank()
    # Repeated real root: diagonalizable (rank 3 - mult) is Type I; one 2x2 null Jordan block is II; 3x3 is III.
    if mult == 3:
        return {0: "I", 1: "II", 2: "III"}[rk], D
    return {1: "I", 2: "II"}[rk], D


def numeric_type(Gmix, tol=1e-8):
    """'IV' where an eigenvalue of the (n, 4, 4) mixed tensor has |Im| > tol * scale, else 'I' (numpy)."""
    ev = np.linalg.eigvals(Gmix)
    scale = np.maximum(np.abs(ev).max(axis=1), np.abs(Gmix).max(axis=(1, 2)))
    return np.abs(ev.imag).max(axis=1) > tol * scale


# ---------------------------------------------------------------------------------------------------------------
# Route C: the frame Ricci tensor of phantom_source.py, with w = v (1 - f)
# ---------------------------------------------------------------------------------------------------------------
V_, C_, S_ = sp.symbols("v c s", real=True)
F0, F1, F2 = sp.symbols("f f1 f2", real=True)  # f, f', f'' at r_s
r_ = ps.r_
FRAME = ["0", "1", "rho", "varphi"]


def reduce_s(E):
    """Remainder of E modulo s^2 - (1 - c^2), as a polynomial in s."""
    num, den = sp.fraction(sp.together(sp.expand(E)))
    return sp.expand(sp.rem(sp.expand(num), S_**2 - 1 + C_**2, S_)) / den


def route_C():
    Rsph = ps.sph_tetrad_ricci("alc-drop-w2" if MUT == "alc-drop-w2" else None)
    trig = {sp.sin(ps.th): S_, sp.cos(ps.th): C_}
    jets = {ps.W[2]: -V_ * F2, ps.W[1]: -V_ * F1, ps.W[0]: V_ * (1 - F0)}
    return Rsph.applyfunc(lambda e: sp.expand(sp.expand_trig(e).subs(trig).subs(jets)))


# ---------------------------------------------------------------------------------------------------------------
# Profile jets, shared by both routes
# ---------------------------------------------------------------------------------------------------------------
def tanh_jets(r, sig, Rb, lib=np):
    """(f, f', f'') of Alcubierre's profile f = [tanh(sig (r + Rb)) - tanh(sig (r - Rb))] / (2 tanh(sig Rb))."""
    tp, tm = lib.tanh(sig * (r + Rb)), lib.tanh(sig * (r - Rb))
    N = 2 * lib.tanh(sig * Rb)
    f2 = 2 * sig**2 * (tm * (1 - tm**2) - tp * (1 - tp**2)) / N
    return (tp - tm) / N, sig * (tm**2 - tp**2) / N, f2 * (1 + 1e-3 if MUT == "jets-f2" else 1)


def wall_radii(sig, Rb, lo=0.1, hi=0.9):
    """Radii where f = hi and f = lo (Le's wall mask f in [lo, hi])."""
    fn = lambda r, t: (mp.tanh(sig * (r + Rb)) - mp.tanh(sig * (r - Rb))) / (2 * mp.tanh(sig * Rb)) - t
    return float(mp.findroot(lambda r: fn(r, hi), Rb - 0.1 / sig)), float(mp.findroot(lambda r: fn(r, lo), Rb + 0.1 / sig))


# ---------------------------------------------------------------------------------------------------------------
# Route N: the Einstein tensor of the lab metric, by hand, for float arrays or object arrays of mpf
# ---------------------------------------------------------------------------------------------------------------
def einstein_N(X, y, z, r, v, f, f1, f2):
    """Mixed G^a_b and R^a_b, shape (n, 4, 4), of ds^2 = -dt^2 + (dx + b dt)^2 + dy^2 + dz^2, b = -v f(r_s), at t = 0.
    X = x - v t. b depends on t only through X, so d_t = -v d_x on it."""
    n = X.shape[0]
    zero = X * 0
    one = zero + 1
    xs = [X, y, z]
    b = -v * f
    bi = [-v * f1 * xi / r for xi in xs]
    bij = [[-v * (f2 * xs[i] * xs[j] / r**2 + f1 * ((one if i == j else zero) / r - xs[i] * xs[j] / r**3))
            for j in range(3)] for i in range(3)]
    db = [-v * bi[0]] + bi
    ddb = [[None] * 4 for _ in range(4)]
    ddb[0][0] = v**2 * bij[0][0]
    for j in range(3):
        ddb[0][j + 1] = ddb[j + 1][0] = -v * bij[0][j]
        for i in range(3):
            ddb[i + 1][j + 1] = bij[i][j]
    keep = 0 if MUT == "num-drop-w2" else 1  # g_tt = -1 + keep * b^2
    dt = object if X.dtype == object else float
    g = np.zeros((n, 4, 4), dtype=dt)
    dg = np.zeros((n, 4, 4, 4), dtype=dt)  # dg[:, a, i, j] = d_a g_ij
    ddg = np.zeros((n, 4, 4, 4, 4), dtype=dt)  # ddg[:, a, c, i, j] = d_a d_c g_ij
    if dt is object:
        g[...] = zero[:, None, None]
        dg[...] = zero[:, None, None, None]
        ddg[...] = zero[:, None, None, None, None]
    g[:, 0, 0] = -1 + keep * b**2
    g[:, 0, 1] = g[:, 1, 0] = b
    for i in (1, 2, 3):
        g[:, i, i] = one
    for a in range(4):
        dg[:, a, 0, 0] = keep * 2 * b * db[a]
        dg[:, a, 0, 1] = dg[:, a, 1, 0] = db[a]
        for c in range(4):
            ddg[:, a, c, 0, 0] = keep * 2 * (db[a] * db[c] + b * ddb[a][c])
            ddg[:, a, c, 0, 1] = ddg[:, a, c, 1, 0] = ddb[a][c]
    det2 = g[:, 0, 0] * g[:, 1, 1] - g[:, 0, 1]**2  # the (t, x) block; y and z are Euclidean
    gi = np.zeros_like(g)
    if dt is object:
        gi[...] = zero[:, None, None]
    gi[:, 0, 0], gi[:, 1, 1] = g[:, 1, 1] / det2, g[:, 0, 0] / det2
    gi[:, 0, 1] = gi[:, 1, 0] = -g[:, 0, 1] / det2
    gi[:, 2, 2] = gi[:, 3, 3] = one
    dgi = -np.einsum("nik,nakl,nlj->naij", gi, dg, gi)
    T = np.einsum("nilj->nlij", dg) + np.einsum("njli->nlij", dg) - dg  # d_i g_lj + d_j g_li - d_l g_ij
    dT = np.einsum("nmilj->nmlij", ddg) + np.einsum("nmjli->nmlij", ddg) - ddg
    Gam = np.einsum("nkl,nlij->nkij", gi, T) / 2
    dGam = (np.einsum("nmkl,nlij->nmkij", dgi, T) + np.einsum("nkl,nmlij->nmkij", gi, dT)) / 2
    Ric = (np.einsum("nkkij->nij", dGam) - np.einsum("njkik->nij", dGam)
           + np.einsum("nkkl,nlij->nij", Gam, Gam) - np.einsum("nkjl,nlik->nij", Gam, Gam))
    Rs = np.einsum("nij,nij->n", gi, Ric)
    G = Ric - Rs[:, None, None] * g / 2
    return np.einsum("nik,nkj->nij", gi, G), np.einsum("nik,nkj->nij", gi, Ric), G


def numN_points(sig, Rb, v, X, y, z):
    r = np.sqrt(X**2 + y**2 + z**2)
    f, f1, f2 = tanh_jets(r, sig, Rb)
    Gm, Rm, G = einstein_N(X, y, z, r, v, f, f1, f2)
    return Gm, Rm, G, f


# ---------------------------------------------------------------------------------------------------------------
# Route C at grid points
# ---------------------------------------------------------------------------------------------------------------
class ClosedForm:
    def __init__(self, Rc):
        args = (V_, F0, F1, F2, r_, C_, S_)
        self.fn = {(a, b): sp.lambdify(args, Rc[a, b], "numpy") for a in range(3) for b in range(a, 3)}
        self.fm = {(a, b): sp.lambdify(args, Rc[a, b], "mpmath") for a in range(3) for b in range(a, 3)}

    def D(self, v, sig, Rb, r, c, mpm=False):
        """Discriminant D and the block's largest entry, at radii r and cos(theta) = c, through `cubic` and `disc`."""
        if mpm:
            f, f1, f2 = tanh_jets(r, sig, Rb, lib=mp)
            s = mp.sqrt(1 - c**2)
            fn = self.fm
        else:
            f, f1, f2 = tanh_jets(r, sig, Rb)
            s = np.sqrt(np.clip(1 - c**2, 0, None))
            fn = self.fn
        comp = {k: fn[k](v, f, f1, f2, r, c, s) for k in fn}
        R = [[comp[(min(a, b), max(a, b))] for b in range(3)] for a in range(3)]
        M = mixed_block(R)
        scale = max(abs(x) for x in comp.values()) if mpm else np.max([np.abs(x) for x in comp.values()], axis=0)
        return disc(*cubic(M)), scale, comp

    def le_split_IV(self, v, sig, Rb, r, c):
        """Le's momentum discriminant (rho_n + S_par)^2 - 4 |j|^2 < 0, from the Einstein tensor in the frame."""
        _, _, k = self.D(v, sig, Rb, r, c)
        Rs = -k[(0, 0)] + k[(1, 1)] + k[(2, 2)]  # R_varphi,varphi = 0
        rho = k[(0, 0)] + Rs / 2
        j1, j2 = k[(0, 1)], k[(0, 2)]
        S11, S12, S22 = k[(1, 1)] - Rs / 2, k[(1, 2)], k[(2, 2)] - Rs / 2
        jj = j1**2 + j2**2
        Spar = np.where(jj > 0, (S11 * j1**2 + 2 * S12 * j1 * j2 + S22 * j2**2) / np.where(jj > 0, jj, 1), 0)
        return (rho + Spar)**2 - 4 * jj < 0


def polar_grid(r1, r2, n):
    """Midpoints of an n x n grid on [r1, r2] x [0, pi] and the proper-volume weights r^2 sin(theta)."""
    rr = r1 + (r2 - r1) * (np.arange(n) + 0.5) / n
    th = np.pi * (np.arange(n) + 0.5) / n
    R, T = np.meshgrid(rr, th, indexing="ij")
    return R.ravel(), np.cos(T).ravel(), (R**2 * np.sin(T)).ravel()


def fraction_I(cf, v, sig, Rb, r1, r2, n):
    r, c, w = polar_grid(r1, r2, n)
    D, _, _ = cf.D(v, sig, Rb, r, c)
    return float((w * (D > 0)).sum() / w.sum())


def compare(cf, sig, Rb, v, r, c, label):
    """Closed-form sign of D against route N's eigenvalues at the points (r, cos theta), random azimuths."""
    rng = np.random.default_rng(64)
    phi = rng.uniform(0, 2 * np.pi, r.size)
    s = np.sqrt(np.clip(1 - c**2, 0, None))
    X, y, z = r * c, r * s * np.cos(phi), r * s * np.sin(phi)
    D, scale, _ = cf.D(v, sig, Rb, r, c)
    ivC = D < 0
    ivN = np.zeros(r.size, dtype=bool)
    for k in range(0, r.size, 20000):
        Gm, _, _, _ = numN_points(sig, Rb, v, X[k:k + 20000], y[k:k + 20000], z[k:k + 20000])
        ivN[k:k + 20000] = numeric_type(Gm)
    near = np.abs(D) < 1e-10 * scale**6
    closest = np.argsort(np.abs(D) / scale**6)[:3]  # the three points nearest the boundary are always rechecked
    redo = np.union1d(closest, np.where((ivC != ivN) | near)[0])
    unresolved = [None] * max(0, redo.size - 100)  # beyond 100 rechecks, count the rest as unresolved
    mp.mp.dps = 50
    for i in redo[:100]:
        rm, cm, ph = mp.mpf(float(r[i])), mp.mpf(float(c[i])), mp.mpf(float(phi[i]))
        sm = mp.sqrt(1 - cm**2)
        Dm, scm, _ = cf.D(mp.mpf(v), sig, Rb, rm, cm, mpm=True)
        Xa, ya, za = [np.array([e], dtype=object) for e in (rm * cm, rm * sm * mp.cos(ph), rm * sm * mp.sin(ph))]
        ra = np.array([rm], dtype=object)
        jm = [np.array([e], dtype=object) for e in tanh_jets(rm, sig, Rb, lib=mp)]
        Gm, _, _ = einstein_N(Xa, ya, za, ra, mp.mpf(v), *jm)
        ev = mp.eig(mp.matrix(Gm[0].tolist()))[0]
        sc = max(abs(e) for e in ev)
        ivNm = max(abs(mp.im(e)) for e in ev) > mp.mpf(10)**-30 * sc
        if (Dm < 0) != ivNm:
            unresolved.append((float(r[i]), float(c[i]), float(Dm), float(max(abs(mp.im(e)) for e in ev))))
    mp.mp.dps = 15
    nd = int((ivC != ivN).sum())
    print(f"  {label}: {r.size} points, Type IV by D < 0 at {int(ivC.sum())}, by eigenvalues at {int(ivN.sum())}; "
          f"float64 disagreements {nd}, near D = 0 {int(near.sum())}, rechecked at 50 digits {min(redo.size, 100)}, "
          f"unresolved {len(unresolved)}")
    return r.size, nd, unresolved


# ---------------------------------------------------------------------------------------------------------------
# Planted truth: tensors of known type, in a generic Lorentz frame, through he_type_exact
# ---------------------------------------------------------------------------------------------------------------
def planted_cases():
    Q = sp.Rational
    eta3 = sp.diag(-1, 1, 1)
    B1 = sp.Matrix([[Q(5, 3), Q(4, 3), 0], [Q(4, 3), Q(5, 3), 0], [0, 0, 1]])  # boost along e_1, gamma = 5/3
    B2 = sp.Matrix([[Q(13, 12), 0, Q(5, 12)], [0, 1, 0], [Q(5, 12), 0, Q(13, 12)]])  # boost along e_2
    Rot = sp.Matrix([[1, 0, 0], [0, Q(3, 5), -Q(4, 5)], [0, Q(4, 5), Q(3, 5)]])
    L = B1 * Rot * B2
    assert L.T * eta3 * L == eta3
    u = sp.Matrix([-1, 0, 0])  # covector of the rest frame
    k = sp.Matrix([1, 1, 0])  # null covector
    m = sp.Matrix([0, 0, 1])
    cases = [
        ("perfect fluid mu = 3/2, p = 1/5", (Q(3, 2) + Q(1, 5)) * u * u.T + Q(1, 5) * eta3, "I"),
        ("dust mu = 2", 2 * u * u.T, "I"),
        ("anisotropic fluid (2, 1/3, -1/7)", sp.diag(2, Q(1, 3), -Q(1, 7)), "I"),
        ("Type IV canonical (rho, nu, p) = (1/2, 2/3, 1/7)",
         sp.Matrix([[Q(1, 2), Q(2, 3), 0], [Q(2, 3), -Q(1, 2), 0], [0, 0, Q(1, 7)]]), "IV"),
        ("Type IV canonical (rho, nu, p) = (-3, 1/9, 5)",
         sp.Matrix([[-3, Q(1, 9), 0], [Q(1, 9), 3, 0], [0, 0, 5]]), "IV"),
        ("null dust + lam eta + q e2 e2", Q(3, 4) * k * k.T + Q(1, 3) * eta3 + Q(2, 5) * m * m.T, "II"),
        ("null dust + lam eta", Q(3, 4) * k * k.T + Q(1, 3) * eta3, "II"),
        ("Le's Type III example l m + m l + lam eta", k * m.T + m * k.T - Q(1, 6) * eta3, "III"),
        ("lam eta", Q(2, 7) * eta3, "I"),
    ]
    out = []
    for name, T3, want in cases:
        T3 = L.T * T3 * L  # T'_ab = L^c_a L^d_b T_cd
        T4 = sp.zeros(4, 4)
        T4[:3, :3] = T3
        T4[3, 3] = Q(-5, 11)
        out.append((name, T4, want))
    return out


# ---------------------------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------------------------
def main():
    global MUT
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    MUT = opts.get("mutate")
    if MUT is not None and MUT not in MUTATIONS:
        sys.exit(f"unknown mutation {MUT!r}. Known: {', '.join(MUTATIONS)}")
    t0 = time.time()
    print(f"he_type.py, issue #64. sympy {sp.__version__}, numpy {np.__version__}, mpmath {mp.__version__}, python "
          f"{sys.version.split()[0]}. mutation: {MUT or 'none'}" + (f" ({MUTATIONS[MUT][1]})" if MUT else ""))

    # 1. Planted truth through the classifier.
    print("\n== 1. Planted truth: tensors of known type in a boosted and rotated frame, e_3 decoupled ==")
    bad = []
    for name, T4, want in planted_cases():
        got, D = he_type_exact(T4)
        Gm = np.array((sp.diag(-1, 1, 1, 1) * T4).tolist(), dtype=float)[None]
        num = ("IV" if numeric_type(Gm)[0] else "I") if D != 0 else "-"
        ok = got == want and (num == "-" or num == want)
        bad += [] if ok else [name]
        print(f"  {name:>46}: D = {sp.nsimplify(D)!s:>14}, closed form {got:>3}, eigenvalues {num:>2}, "
              f"expected {want}")
    check("planted", not bad, f"every planted tensor gets its type from D and the rank test; wrong: {bad}")

    # 2. Route C: decoupling and the closed form.
    print("\n== 2. Route C: the frame Ricci tensor with w = v (1 - f), and the cubic of its (E_0, E_1, e_rho) block ==")
    Rc = route_C()
    names = {0: "0", 1: "1", 2: "rho", 3: "varphi"}
    for a in range(4):
        for b in range(a, 4):
            print(f"  R_({names[a]},{names[b]}) = {sp.factor(Rc[a, b])}")
    dec = [reduce_s(Rc[a, 3]) for a in range(4)] + [reduce_s(Rc[0, 0] + Rc[1, 1] + 2 * Rc[2, 2])]
    check("decouple", all(e == 0 for e in dec), f"R(e_varphi, .) and R_00 + R_11 + 2 R_rho,rho reduce to {dec}")

    a_, b_, F_ = F1, F2, 1 - F0
    s2 = 1 - C_**2
    J2 = s2 * (a_ + r_ * b_)**2 + 4 * C_**2 * a_**2
    p6 = -J2**3 / (16 * r_**6)
    Q1 = -C_**2 * F_ * (a_ - r_ * b_)**2 - 4 * r_ * a_**3
    Q2 = a_**2 * (F_ + r_ * a_)**2 + C_**2 * F_ * (a_ + r_ * b_) * (F_ * (r_ * b_ - a_) - 2 * r_ * a_**2)
    p12 = s2**2 * F_**2 * Q1**2 * Q2 / (16 * r_**6)
    Qeq = 4 * V_**2 * a_**2 * F_ * (F_ + r_ * a_) + (a_ + r_ * b_)**2
    Qax = a_**2 + r_**2 * V_**2 * (F_ * b_ - a_**2)**2
    eqD = V_**6 * (r_**2 * V_**2 * a_**4 - (a_ + r_ * b_)**2) * Qeq**2 / (16 * r_**6)
    axD = -4 * a_**2 * V_**6 * Qax**2 / r_**6

    # Screen at exact rational points, with cos and sin from tan(theta / 2) = t. D / v^6 at v = 1..5 must be one cubic in
    # v^2 with the expected p6 and p12, and D must equal the equator and axis forms. The symbolic proof runs only after.
    Qr = sp.Rational
    bad = []
    for (fv, f1v, f2v, rv, tv) in [(Qr(1, 3), Qr(-5, 2), Qr(7, 4), Qr(9, 10), Qr(2, 7)),
                                   (Qr(4, 5), Qr(-1, 3), Qr(-2, 9), Qr(6, 5), Qr(5, 3)),
                                   (Qr(1, 9), Qr(-3), Qr(11, 2), Qr(1), Qr(1, 4))]:
        for (cv, sv) in [((1 - tv**2) / (1 + tv**2), 2 * tv / (1 + tv**2)), (0, 1), (1, 0)]:
            pt = {F0: fv, F1: f1v, F2: f2v, r_: rv, C_: cv, S_: sv}
            Dv = []
            for vv in range(1, 6):
                Rn = Rc.subs(pt).subs(V_, vv)
                Dv.append(disc(*cubic(mixed_block([[Rn[a, b] for b in range(3)] for a in range(3)]))) / vv**6)
            co = sp.Matrix([[vv**(2 * k) for k in range(4)] for vv in range(1, 5)]).solve(sp.Matrix(Dv[:4]))
            ok = sum(co[k] * 25**k for k in range(4)) == Dv[4] and co[0] == p6.subs(pt) and co[3] == p12.subs(pt)
            if cv == 0:
                ok = ok and Dv[1] * 2**6 == eqD.subs(pt).subs(V_, 2)
            if cv == 1:
                ok = ok and Dv[1] * 2**6 == axD.subs(pt).subs(V_, 2)
            bad += [] if ok else [(fv, f1v, f2v, rv, cv)]
    print(f"  screen at 9 exact points (f, f', f'', r, cos theta): failures at {bad}")
    p = None
    res = {"exact-point screen": not bad}
    if not bad:
        M = mixed_block([[Rc[a, b] for b in range(3)] for a in range(3)])
        c2, c1, c0 = [reduce_s(e) for e in cubic(M)]
        print(f"  c2 = {sp.factor(c2)}")
        D = reduce_s(disc(c2, c1, c0))
        num, den = sp.fraction(sp.together(D))
        Pv = sp.Poly(sp.expand(num), V_)
        powers = sorted(mm[0] for mm in Pv.monoms())
        p = {k: sp.expand(Pv.coeff_monomial(V_**k)) / den for k in powers}
        print(f"  D has the powers v^{powers} (numerator {len(sp.expand(num).args)} terms, denominator {den})")
        jj = reduce_s(Rc[0, 1]**2 + Rc[0, 2]**2)
        # Static part: the v^2 terms. At v = 1 it is block diagonal, so its D is chi(lam0)^2 (mu1 - mu2)^2 >= 0.
        st = Rc.applyfunc(lambda e: sp.expand(e).coeff(V_, 2))
        lam0, S11, S12, S22 = -st[0, 0], st[1, 1], st[1, 2], st[2, 2]
        chi = lam0**2 - (S11 + S22) * lam0 + S11 * S22 - S12**2
        res.update({
            "powers": powers == [6, 8, 10, 12],
            "p6 = -(J^2)^3/(16 r^6)": sp.expand(p.get(6, 0) - p6) == 0,
            "p6 = -4 |j|^6 / v^6": sp.expand(reduce_s(p.get(6, 0) * V_**6 + 4 * jj**3)) == 0,
            "p12 = s^4 F^2 Q1^2 Q2/(16 r^6)": sp.expand(p.get(12, 0) - p12) == 0,
            "p12 = chi(lam0)^2 [(S11 - S22)^2 + 4 S12^2]":
                sp.expand(reduce_s(p.get(12, 0) - chi**2 * ((S11 - S22)**2 + 4 * S12**2))) == 0,
            "equator": sp.expand(D.subs(C_, 0) - eqD) == 0,
            "axis": sp.expand(D.subs(C_, 1) - axD) == 0,
        })
    for k, ok in res.items():
        print(f"  {k}: {ok}")
    print(f"  Q1 = {Q1}\n  Q2 = {Q2}\n  Q_eq = {Qeq}\n  Q_ax = {Qax}")
    for k in (8, 10):
        if p is not None and k in p:
            print(f"  p{k} (coefficient of v^{k} in D) = {sp.factor(p[k])}")
    check("structure", all(res.values()), "D = v^6 (p6 + p8 v^2 + p10 v^4 + p12 v^6) with p6 = -4|j|^6/v^6 < 0, "
          "p12 >= 0, and the equator and axis factorizations: " + ", ".join(k for k, ok in res.items() if not ok))

    # 3. Route N, and the profile jets.
    print("\n== 3. Route N (the lab Einstein tensor by hand) and the tanh jets ==")
    mp.mp.dps = 30
    worst = 0
    for sig, Rb, rs in ((8, 1, (0.3, 0.9, 1.0, 1.1, 1.7)), (4, 5, (4.5, 5.0, 5.2))):
        fm = lambda r: (mp.tanh(sig * (r + Rb)) - mp.tanh(sig * (r - Rb))) / (2 * mp.tanh(sig * Rb))
        for r0 in rs:
            jn = tanh_jets(np.array([r0]), sig, Rb)
            for kk in range(3):
                ex = mp.diff(fm, mp.mpf(r0), kk)
                worst = max(worst, float(abs(jn[kk][0] - ex)) / sig**kk)  # relative to the jet's scale sig^k
    mp.mp.dps = 15
    check("jets", worst < 1e-9, f"tanh_jets against mpmath's numerical derivatives at 8 radii, worst relative error "
          f"{worst:.1e}")
    rng = np.random.default_rng(57)
    Xs, ys, zs = rng.uniform(-1.3, 1.3, 40), rng.uniform(-1.3, 1.3, 40), rng.uniform(-1.3, 1.3, 40)
    Gm, Rm, G, f = numN_points(8, 1, 0.7, Xs, ys, zs)
    r = np.sqrt(Xs**2 + ys**2 + zs**2)
    _, f1, _ = tanh_jets(r, 8, 1)
    b = -0.7 * f
    Gnn = G[:, 0, 0] - 2 * b * G[:, 0, 1] + b**2 * G[:, 1, 1]  # n = d_t - b d_x
    rhoE = -(0.7**2 / (32 * np.pi)) * (ys**2 + zs**2) / r**2 * f1**2
    err_rho = float(np.max(np.abs(Gnn / (8 * np.pi) - rhoE)) / np.max(np.abs(rhoE)))
    _, Rax, _, _ = numN_points(8, 1, 1.0, np.array([1.0, -1.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]))
    evs = [np.sort_complex(np.linalg.eigvals(Rax[i])) for i in range(2)]
    fa, fa1, fa2 = tanh_jets(np.array([1.0]), 8, 1)
    mR00, R01 = float(fa1[0]**2 - (1 - fa[0]) * fa2[0]), float(fa1[0])  # -R_00 and R_01 on the axis, v = r_s = 1
    want = np.sort_complex(np.array([0, 0, mR00 - 1j * R01, mR00 + 1j * R01]))
    err_ax = float(max(np.max(np.abs(e - want)) for e in evs)) / 16
    print(f"  R^a_b on the axis at x - v t = +-1 (sig = 8, R = 1, v = 1): {evs[0]}, {evs[1]}; the closed form's "
          f"-R_00 +- i R_01 = {mR00:.9f} +- {abs(R01):.9f} i")
    check("numN", err_rho < 1e-8 and err_ax < 1e-10,
          f"route N's G(n, n)/8pi equals -(v^2/32pi)(rho^2/r_s^2) f'^2 at 40 points (relative error {err_rho:.1e}), "
          f"and its R^a_b on the axis has eigenvalues 0, 0, -R_00 +- i R_01 (v^2 (f'^2 - (1 - f) f'') +- i v f'/r_s, "
          f"phantom_source's 16 +- 4i to 6 digits; relative error {err_ax:.1e})")

    # 4. The closed form against numerical eigenvalues at every grid point.
    print("\n== 4. Closed form (sign of D) against route N's eigenvalues at every grid point ==")
    cf = ClosedForm(Rc)
    total, nd_tot, unres = 0, 0, []
    le = (8, 1)
    r1, r2 = wall_radii(*le)
    print(f"  Le's wall, sig = 8, R_b = 1: f in [0.1, 0.9] is r_s in [{r1:.6f}, {r2:.6f}], width {r2 - r1:.4f}")
    for v in (0.5, 1.0, 2.0, 2.5):
        rr, cc, _ = polar_grid(r1, r2, 240)
        n, nd, un = compare(cf, *le, v, rr, cc, f"Le's wall, polar 240 x 240, v = {v}")
        total, nd_tot, unres = total + n, nd_tot + nd, unres + un
    g1 = (np.arange(100) + 0.5) * 0.06 - 3  # Le's box +-3 R_b, N = 100, uniform here
    Xg, yg, zg = np.meshgrid(g1, g1, g1, indexing="ij")
    rg = np.sqrt(Xg**2 + yg**2 + zg**2).ravel()
    inw = (rg >= r1) & (rg <= r2)
    n, nd, un = compare(cf, *le, 0.5, rg[inw], Xg.ravel()[inw] / rg[inw], "Le's wall, uniform 100^3 on +-3, v = 0.5")
    total, nd_tot, unres = total + n, nd_tot + nd, unres + un
    ro = (4, 5)
    vr = 0.25 if MUT == "rodal-speed" else 1.0
    rr, cc, _ = polar_grid(5 - 3 / 4, 5 + 3 / 4, 240)
    n, nd, un = compare(cf, *ro, vr, rr, cc, f"Rodal's band r in [rho - 3/sig, rho + 3/sig], polar 240 x 240, v = {vr}")
    total, nd_tot, unres = total + n, nd_tot + nd, unres + un
    rax = np.linspace(5 - 3 / 4, 5 + 3 / 4, 301)
    n, nd, un = compare(cf, *ro, vr, np.concatenate([rax, rax]), np.concatenate([np.ones(301), -np.ones(301)]),
                        "Rodal's axis, theta = 0 and pi")
    total, nd_tot, unres = total + n, nd_tot + nd, unres + un
    if unres:
        print(f"  unresolved at 50 digits (r, cos theta, D, max|Im lam|), first 5 rechecked: "
              f"{[u for u in unres if u is not None][:5]}")
    check("closed_vs_numeric", not unres, f"the sign of D gives the eigenvalues' type at all {total} grid points "
          f"({nd_tot} float64 disagreements; {len(unres)} unresolved after up to 100 rechecks per grid at 50 digits)")

    # 5. Type I fractions for Le's setup.
    print("\n== 5. Le's setup: R_b = 1, sig = 8, t = 0, wall f in [0.1, 0.9], proper-volume weighted ==")
    print("  | v | Type I %, polar 1000^2 | Type I %, polar 2000^2 | uniform 100^3 | uniform 200^3 | Le's Table 3 | "
          "Le's split-discriminant agreement % |")
    le_pub = {0.5: 1.3, 1.0: 5.7, 2.0: 16.7, 2.5: 24.9}
    fr = {}
    for v in le_pub:
        fa = fraction_I(cf, v, *le, r1, r2, 1000)
        fb = fraction_I(cf, v, *le, r1, r2, 2000)
        cart = []
        for N in (100, 200):
            g1 = (np.arange(N) + 0.5) * (6 / N) - 3
            Xg, yg, zg = np.meshgrid(g1, g1, g1, indexing="ij")
            rg = np.sqrt(Xg**2 + yg**2 + zg**2).ravel()
            inw = (rg >= r1) & (rg <= r2)
            Dg, _, _ = cf.D(v, *le, rg[inw], Xg.ravel()[inw] / rg[inw])
            cart.append(float((Dg > 0).mean()))
        rr, cc, ww = polar_grid(r1, r2, 1000)
        Dp, _, _ = cf.D(v, *le, rr, cc)
        agree = float((ww * ((Dp < 0) == cf.le_split_IV(v, *le, rr, cc))).sum() / ww.sum())
        fr[v] = fb
        print(f"  | {v} | {100 * fa:.3f} | {100 * fb:.3f} | {100 * cart[0]:.2f} | {100 * cart[1]:.2f} | "
              f"{le_pub[v]} | {100 * agree:.2f} |")
    vle = 1.0 if MUT == "le-speed" else 0.5
    fle = fr[vle]
    check("le_fraction", abs(100 * fle - 1.3) <= 2.0, f"Type I fraction at v = {vle}: {100 * fle:.2f} % against Le's "
          f"1.3 % (tolerance 2 percentage points)")

    # 6. Rodal's setup.
    print("\n== 6. Rodal's setup: rho = 5, sig = 4, v = 1, cuts at theta = 0 and pi/2 ==")
    rr = np.linspace(5 - 3 / 4, 5 + 3 / 4, 150001)
    Dax, _, _ = cf.D(vr, *ro, rr, np.ones_like(rr))
    Deq, _, _ = cf.D(vr, *ro, rr, np.zeros_like(rr))
    f, f1, f2 = tanh_jets(rr, *ro)
    eq_rule = rr * vr * f1**2 > np.abs(f1 + rr * f2)  # Type I on the equator by the factorized form
    typeI = Deq > 0
    edges = np.flatnonzero(np.diff(typeI.astype(int)))
    fD = lambda x: float(cf.D(vr, *ro, np.array([x]), np.array([0.0]))[0][0])
    segs = [brentq(fD, rr[e], rr[e + 1], xtol=1e-12) for e in edges]
    print(f"  axis: Type IV at {int((Dax < 0).sum())} of {rr.size} radii in [4.25, 5.75]")
    print(f"  equator: sign changes of D at r = {[round(x, 5) for x in segs]}; the rule r v f'^2 > |f' + r f''| agrees "
          f"with D > 0 at {int((eq_rule == typeI).sum())} of {rr.size} radii")
    Delta = 0.014 * 5
    one_ann = len(segs) == 2 and segs[0] < 5 < segs[1]
    width = segs[1] - segs[0] if one_ann else float("nan")
    print(f"  equatorial Type I annulus: r in [{segs[0] if segs else 'none'}, {segs[-1] if segs else 'none'}], width "
          f"{width:.4f} = {width / (2 * Delta):.2f} x Rodal's 2 Delta = {2 * Delta:.3f}")
    for lab, (a1, a2) in (("Le's mask f in [0.1, 0.9]", wall_radii(*ro)), ("Rodal's band [rho - 3/sig, rho + 3/sig]",
                                                                          (5 - 3 / 4, 5 + 3 / 4))):
        print(f"  Type I fraction, {lab}, r in [{a1:.4f}, {a2:.4f}]: polar 1000^2 "
              f"{100 * fraction_I(cf, vr, *ro, a1, a2, 1000):.3f} %, 2000^2 {100 * fraction_I(cf, vr, *ro, a1, a2, 2000):.3f} %")
    check("rodal_cuts", bool((Dax < 0).all()) and one_ann and 0.5 <= width / (2 * Delta) <= 2
          and bool((eq_rule == typeI).all()),
          f"Type IV on the whole axis, and one equatorial Type I annulus about r = rho of width {width:.4f}, within a "
          f"factor 2 of Rodal's 2 Delta = {2 * Delta:.2f}")

    # 7. The speed structure: D / v^6 is a cubic in V = v^2 with p6 < 0 and p12 >= 0, so it has one or three positive
    # roots V1 <= V2 <= V3. Type IV for v^2 < V1, Type I above V1 with one root; with three, Type I on (V1, V2), Type IV
    # again on (V2, V3) and Type I above V3.
    print("\n== 7. Speed structure: the positive roots of p12 V^3 + p10 V^2 + p8 V + p6, V = v^2 (numerical evidence) ==")
    pf = {k: sp.lambdify((F0, F1, F2, r_, C_), p[k] * 16 * r_**6, "numpy", cse=True) for k in p} if p else {}
    for lab, (sig, Rb), (a1, a2) in (("Le", le, (r1, r2)), ("Rodal", ro, wall_radii(*ro))) if p else ():
        rr, cc, ww = polar_grid(a1, a2, 400)
        f, f1, f2 = tanh_jets(rr, sig, Rb)
        co = np.array([pf[k](f, f1, f2, rr, cc) * np.ones_like(rr) for k in sorted(p)])
        cm = np.zeros((rr.size, 3, 3))  # companion matrix of V^3 + (p10 V^2 + p8 V + p6) / p12
        cm[:, 0, :] = -(co[[2, 1, 0]] / co[3]).T
        cm[:, 1, 0] = cm[:, 2, 1] = 1
        roots = np.linalg.eigvals(cm)
        isreal = np.abs(roots.imag) <= 1e-9 * np.abs(roots).max(axis=1, keepdims=True)
        npos = ((roots.real > 0) & isreal).sum(axis=1)
        posr = np.sort(np.where((roots.real > 0) & isreal, roots.real, np.inf), axis=1)
        v1 = np.sqrt(posr[npos >= 1, 0])
        i3 = np.flatnonzero(npos == 3)
        print(f"  {lab}'s wall, polar 400^2: one positive root on {100 * ww[npos == 1].sum() / ww.sum():.3f} % of the "
              f"volume, three on {100 * ww[npos == 3].sum() / ww.sum():.3f} %, none resolved by the companion matrix at "
              f"{int((npos == 0).sum())} points next to the axis, where p12 is of order sin(theta)^4. The lowest root "
              f"v1 = sqrt(V1) ranges over [{v1.min():.4f}, {v1.max():.3g}], and the re-entrant Type IV band, where "
              f"there is one, starts at v2 = sqrt(V2) >= {np.sqrt(posr[i3, 1]).min() if i3.size else np.nan:.3f}")
        if i3.size:
            k = i3[np.argmax(posr[i3, 1] - posr[i3, 0])]  # the widest first Type I window
            vv = np.sqrt(posr[k])
            probes = [vv[0] / 2, (vv[0] + vv[1]) / 2, (vv[1] + vv[2]) / 2, 2 * vv[2]]
            s_ = np.sqrt(1 - cc[k]**2)
            labs = []
            for x in probes:
                Gm, _, _, _ = numN_points(sig, Rb, x, np.array([rr[k] * cc[k]]), np.array([rr[k] * s_]), np.array([0.0]))
                labs.append("IV" if numeric_type(Gm)[0] else "I")
            print(f"    at r = {rr[k]:.5f}, cos theta = {cc[k]:.5f}: roots v = {np.round(vv, 4)}; route N's type at "
                  f"v = {np.round(probes, 4)} is {labs}")
        vs = np.array([0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 5, 10])
        fI = [fraction_I(cf, x, sig, Rb, a1, a2, 600) for x in vs]
        print("    v:        " + " ".join(f"{x:6g}" for x in vs))
        print("    Type I %: " + " ".join(f"{100 * x:6.2f}" for x in fI))

    failed = [nm for nm, ok in RESULTS if not ok]
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - t0:.0f} s")
    if MUT is not None:
        want = MUTATIONS[MUT][0]
        print(f"control {MUT}: target checks {want}, fired {[n for n in want if n in failed]}")
    return failed


if __name__ == "__main__":
    failed = main()
    controls_ok = True
    if "--mutate" not in " ".join(sys.argv) and "--no-controls" not in sys.argv:
        print("\n== Controls: each mutation must make its target checks fail ==", flush=True)
        for name, (want, what) in MUTATIONS.items():
            pr = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
            fails = [ln.split(":")[0][7:] for ln in pr.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = all(w in fails for w in want) and pr.returncode == 1
            controls_ok &= fired
            print(f"[{'FIRED' if fired else 'SILENT'}] --mutate={name} ({what}): exit {pr.returncode}, failed {fails}",
                  flush=True)
            shown = [ln.split(": ", 1)[1] for ln in pr.stdout.splitlines() if "] le_fraction: " in ln]
            print(f"    its map: {shown[0] if shown else 'no le_fraction line'}")
            if pr.returncode not in (0, 1):
                print(pr.stderr[-2000:])
        print(f"controls: {'all fired' if controls_ok else 'SOME DID NOT FIRE'}")
    sys.exit(0 if not failed and controls_ok else 1)
