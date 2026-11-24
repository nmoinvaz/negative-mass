#!/usr/bin/env python3
"""Hawking-Ellis type of the accelerating Alcubierre bubble, in closed form (issue #73).

    ds^2 = -dt^2 + (dx + b dt)^2 + dy^2 + dz^2,   b = -v_s(t) f(r_s),   r_s^2 = (x - x_s(t))^2 + y^2 + z^2,

with v_s(t) and x_s(t) independent free functions. Notation: v = v_s, al = v_s', u = x_s' (Alcubierre's journey is u = v).
As in he_type.py, G^a_b and R^a_b have the same type, (I1) R(e_varphi, .) = 0 makes e_varphi decouple, and the type is
that of the mixed block M = eta R on (E_0 = n, E_1 = d_x, e_rho), decided by the sign of the discriminant D of its cubic
(D > 0 Type I, D < 0 Type IV) and, on D = 0, by the rank test (`ht.he_type_exact`).

    python3 notes/warp/he_accel.py                 all checks, then every control (exit 1 if any check fails)
    python3 notes/warp/he_accel.py --mutate=NAME   one control (its target checks must fail, exit 1)
    python3 notes/warp/he_accel.py --no-controls   the checks only
    python3 notes/warp/he_accel.py --stop=7        sections 1 to 7 only (the symbolic derivation, no grids)
    python3 notes/warp/he_accel.py --quick         the grids of the README's `--quick` (a subset of the full ones)
    python3 notes/warp/he_accel.py --controls-quick   as the default run, with `--quick` for each control

Imported unchanged. From phantom_source.py (`ps`): `ricci`, `shift_metric`, `shift_tetrad_ricci`, `bubble_rule`, `JetD`,
`ricci_at_point`, the jet symbols VS, XS, F and X, y, z, R (route D, #61's frame Ricci tensor of the bubble, and route E,
its exact-point planted truth). From he_type.py (`ht`): `mixed_block`, `cubic`, `disc`, `he_type_exact`, `planted_cases`,
`route_C`, `ClosedForm`, `numeric_type`, `einstein_N`, `tanh_jets`, `wall_radii`, `polar_grid`, `reduce_s` and its symbols.

Routes. (D) `ps.shift_tetrad_ricci` with the jets of v_s, x_s and f as free symbols, at z = 0, y > 0, where
(E_2, E_3) = (e_rho, e_varphi); (E) `ps.ricci_at_point`, exact rationals, sympy's own derivatives of explicit functions;
(N') `einstein_b`, a numerical Einstein tensor for an explicit b(t, x, y, z) whose derivatives are sympy's `diff` of the
explicit function, with the Christoffel and Ricci einsums of he_type's `einstein_N` (copied, tested by `routeN_planted`).
D and E share the metric's definition `ps.shift_metric`. N' shares with D only the metric's form and the profile f.
"""
import os
import re
import subprocess
import sys
import time

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import he_type as ht  # noqa: E402
import phantom_source as ps  # noqa: E402

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "disc-sign": (["planted_cls", "closed_vs_numeric"], "flips the sign of the cubic discriminant (he_type's switch)"),
    "no-rank": (["planted_cls"], "labels every D = 0 tensor Type II, without the rank test (he_type's switch)"),
    "num-drop-w2": (["routeN_planted", "routeN_vs_E", "closed_vs_numeric"], "drops b^2 from g_tt in route N'"),
    "num-static": (["routeN_planted", "routeN_vs_E", "closed_vs_numeric"], "zeroes every t-derivative of b in route N'"),
    "accel-static": (["entries_E", "const_limit", "closed_vs_numeric"],
                     "phantom_source's control: drops every t-derivative of b, v_s and x_s in route D's jets"),
    "accel-lapse": (["entries_E", "decouple"], "phantom_source's control: lapse 1 + b in route D"),
    "alc-drop-w2": (["entries_E", "decouple", "const_limit", "axis", "closed_vs_numeric"],
                    "drops b^2 from g_tt in route D"),
    "drop-accel": (["closed_vs_numeric", "drop_detected"], "evaluates D at al = 0, u = v (the constant-speed D)"),
    "time-lemma-sign": (["time_lemma"], "flips the sign of the DeltaT_rho term in the time lemma's claim"),
    "weight-break": (["structure"], "counts al as weight 1, not 2, in the scaling of the entries"),
    "parity-break": (["parity"], "flips al but not c in the parity identity"),
    "axis-flip": (["axis"], "claims D = +4 R_01^2 (R_00^2 + R_01^2)^2 on the axis"),
    "rest-v": (["rest"], "takes the rest slice at t = 0.05 (v = tanh 0.1), and v = 1/10 at the exact points"),
}

RESULTS = []
MUT = None
QUICK = False
T0 = time.time()


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# Symbols. he_type's V_, C_, S_, F0, F1, F2 and phantom_source's r_ are used so the entries compare with route C.
# ---------------------------------------------------------------------------------------------------------------
V_, C_, S_, F0, F1, F2, r_ = ht.V_, ht.C_, ht.S_, ht.F0, ht.F1, ht.F2, ps.r_
AL, U = sp.symbols("alpha u", real=True)
JETS = (V_, AL, U, F0, F1, F2, r_, C_, S_)
reduce_s = ht.reduce_s


def route_D(mut):
    """Frame Ricci tensor (E_0 = n, E_1, E_2, E_3) = (n, d_x, d_y, d_z) of the accelerating bubble, jets as symbols."""
    b = -ps.VS[0] * ps.F[0]
    if mut == "accel-static":
        return ps.shift_tetrad_ricci(b, ps.bubble_rule(True), False)
    if mut == "accel-lapse":
        return ps.shift_tetrad_ricci(b, ps.bubble_rule(False), True)
    if mut == "alc-drop-w2":
        g, E = ps.shift_metric(b, False)
        g[0, 0] = -1
        Ric = ps.ricci(g, [0, 1, 2, 3], ps.JetD(ps.bubble_rule(False)), sp.expand)
        return (E.T * Ric * E).applyfunc(sp.expand)
    return ps.shift_tetrad_ricci(b, ps.bubble_rule(False), False)


def block_z0(Ra):
    """Ra at z = 0, y > 0 in (c, s) = (cos theta, sin theta), s^2 reduced to 1 - c^2. At z = 0, e_rho = E_2, e_varphi = E_3."""
    out = Ra.subs({ps.z: 0}).subs({ps.X: ps.R * C_, ps.y: ps.R * S_})
    out = out.subs({ps.R: r_, ps.VS[0]: V_, ps.VS[1]: AL, ps.XS[1]: U, ps.F[0]: F0, ps.F[1]: F1, ps.F[2]: F2})
    return out.applyfunc(reduce_s)


def closed_pieces():
    """The closed forms of the block's entries, in the jets (v, al, u, f, f', f'')."""
    f, a, b, v, al, u, r, c, s = F0, F1, F2, V_, AL, U, r_, C_, S_
    Tx = -al * c * a + v * u * (s**2 * a / r + c**2 * b)  # d_t d_x b
    Tr = -al * s * a + v * u * c * s * (b - a / r)  # d_t d_rho b
    A = v**2 / (2 * r) * (2 * f * (a * s**2 + r * b * c**2) + r * a**2 * (1 + c**2))
    B = s * c * v**2 / (2 * r) * (2 * a**2 * r + f * b * r - f * a)
    q = v**2 * a**2 * s**2
    j1 = v * (a * (1 + c**2) + r * b * s**2) / (2 * r)
    jr = v * c * s * (a - r * b) / (2 * r)
    return dict(Tx=Tx, Tr=Tr, A=A, B=B, q=q, j1=j1, jr=jr)


def closed_block(P):
    """R_00 = Tx - A, R_01 = j1, R_0rho = jr, R_11 = A - q - Tx, R_1rho = B - Tr/2, R_rhorho = q/2, e_varphi decoupled."""
    h = P["B"] - P["Tr"] / 2
    return sp.Matrix([[P["Tx"] - P["A"], P["j1"], P["jr"], 0], [P["j1"], P["A"] - P["q"] - P["Tx"], h, 0],
                      [P["jr"], h, P["q"] / 2, 0], [0, 0, 0, 0]])


def D_shift(j1, jr, q, h, p):
    """Discriminant of the cubic of M - M_00 I, M = eta R, with q = 2 R_rhorho, h = R_1rho, p = R_rhorho + R_00 and (I2):
    M - M_00 I = [[0, -j1, -jr], [j1, -q, h], [jr, h, p]], cubic l^3 + (q - p) l^2 + (j1^2 + jr^2 - p q - h^2) l
    + (2 j1 jr h - p j1^2 + q jr^2)."""
    c2 = q - p
    c1 = j1**2 + jr**2 - p * q - h**2
    c0 = 2 * j1 * jr * h - p * j1**2 + q * jr**2
    return 18 * c2 * c1 * c0 - 4 * c2**3 * c0 + c2**2 * c1**2 - 4 * c1**3 - 27 * c0**2


def shifted_args(blk):
    """(j1, jr, q, h, p) read off a block with (I2)."""
    return blk[0, 1], blk[0, 2], 2 * blk[2, 2], blk[1, 2], blk[2, 2] + blk[0, 0]


def D_direct(blk):
    """D of the unshifted block through he_type's mixed_block, cubic and disc, as an expression."""
    return ht.disc(*ht.cubic(ht.mixed_block([[blk[a, b] for b in range(3)] for a in range(3)])))


def exact_R4(blk, vals):
    return blk.subs(vals).applyfunc(sp.cancel)


# ---------------------------------------------------------------------------------------------------------------
# Profiles of the motion and the tanh profile of the wall (as he_type)
# ---------------------------------------------------------------------------------------------------------------
T_ = sp.Symbol("t", real=True)
XX, YY, ZZ = sp.symbols("x y z", real=True)
PROFILES = {
    "P1": (1 / sp.cosh(2 * T_)**2, sp.tanh(2 * T_) / 2),
    "P2": (sp.tanh(2 * T_), sp.Rational(3, 10) * sp.log(sp.cosh(2 * T_)) - sp.Rational(2, 15) * sp.cos(3 * T_) + T_ / 4),
}
LE = (8, 1)
RODAL = (4, 5)
TIMES = (-0.6, -0.3, 0.0, 0.3, 0.6)


def motion(name, t0):
    """(v, al, u, x_s) at time t0, floats."""
    vs, xs = PROFILES[name]
    ev = lambda e: float(e.subs(T_, sp.Float(t0, 30)))  # noqa: E731
    return ev(vs), ev(sp.diff(vs, T_)), ev(sp.diff(xs, T_)), ev(xs)


def motion_mp(name, t0):
    vs, xs = PROFILES[name]
    fs = [sp.lambdify(T_, e, "mpmath") for e in (vs, sp.diff(vs, T_), sp.diff(xs, T_), xs)]
    return [fn(mp.mpf(t0)) for fn in fs]


def tanh_f(rr, sig, Rb):
    return (sp.tanh(sig * (rr + Rb)) - sp.tanh(sig * (rr - Rb))) / (2 * sp.tanh(sig * Rb))


def slices():
    ts = (-0.3, 0.3) if QUICK else TIMES
    return [(p, t) for p in ("P1", "P2") for t in ts]


def pgrid(r1, r2, n=240):
    """ht.polar_grid, or every 2nd index of both axes of it in quick mode (a subset of the full grid's points)."""
    rr, cc, ww = ht.polar_grid(r1, r2, n)
    if QUICK:
        idx = np.arange(n * n).reshape(n, n)[::2, ::2].ravel()
        rr, cc, ww = rr[idx], cc[idx], ww[idx]
    return rr, cc, ww


# ---------------------------------------------------------------------------------------------------------------
# Closed form at grid points, through he_type's mixed_block, cubic and disc
# ---------------------------------------------------------------------------------------------------------------
class ClosedFormAccel:
    def __init__(self, blk):
        self.fn = {(a, b): sp.lambdify(JETS, blk[a, b], "numpy") for a in range(3) for b in range(a, 3)}
        self.fm = {(a, b): sp.lambdify(JETS, blk[a, b], "mpmath") for a in range(3) for b in range(a, 3)}

    def _D(self, v, al, u, sig, Rb, r, c, mpm):
        if mpm:
            f, f1, f2 = ht.tanh_jets(r, sig, Rb, lib=mp)
            s = mp.sqrt(1 - c**2)
            fn = self.fm
            comp = {k: fn[k](v, al, u, f, f1, f2, r, c, s) for k in fn}
            scale = max(abs(x) for x in comp.values())
        else:
            f, f1, f2 = ht.tanh_jets(r, sig, Rb)
            s = np.sqrt(np.clip(1 - c**2, 0, None))
            fn = self.fn
            one = np.ones_like(r)
            comp = {k: fn[k](v, al, u, f, f1, f2, r, c, s) * one for k in fn}
            scale = np.max([np.abs(x) for x in comp.values()], axis=0)
        R = [[comp[(min(a, b), max(a, b))] for b in range(3)] for a in range(3)]
        return ht.disc(*ht.cubic(ht.mixed_block(R))), scale, comp

    def D(self, v, al, u, sig, Rb, r, c, mpm=False):
        """Discriminant D, the block's largest entry and its entries, at radii r and cos(theta) = c."""
        if MUT == "drop-accel":
            al, u = 0 * al, v
        return self._D(v, al, u, sig, Rb, r, c, mpm)

    def D_dropped(self, v, al, u, sig, Rb, r, c):
        """D with the acceleration terms dropped: al = 0, u = v at the same v. Not changed by any mutation."""
        return self._D(v, 0 * al, v, sig, Rb, r, c, False)


def fraction_I(cf, v, al, u, sig, Rb, r1, r2, n, chunk=100):
    """Proper-volume weighted Type I fraction (D > 0) of the wall r in [r1, r2], midpoint rule on n x n (r, theta)."""
    rr = r1 + (r2 - r1) * (np.arange(n) + 0.5) / n
    th = np.pi * (np.arange(n) + 0.5) / n
    num = den = 0.0
    for k in range(0, n, chunk):
        Rg, Tg = np.meshgrid(rr[k:k + chunk], th, indexing="ij")
        D, _, _ = cf.D(v, al, u, sig, Rb, Rg.ravel(), np.cos(Tg).ravel())
        w = (Rg**2 * np.sin(Tg)).ravel()
        num += float((w * (D > 0)).sum())
        den += float(w.sum())
    return num / den


# ---------------------------------------------------------------------------------------------------------------
# Route N': the Einstein tensor of an explicit shift b(t, x, y, z), numerically
# ---------------------------------------------------------------------------------------------------------------
def einstein_b(b, db, ddb):
    """Mixed G^a_b and R^a_b, shape (n, 4, 4), of ds^2 = -dt^2 + (dx + b dt)^2 + dy^2 + dz^2 at n points, from b, its first
    derivatives db[a] and its second derivatives ddb[a][c] with respect to (t, x, y, z). The tensor algebra is that of
    ht.einstein_N, which this generalizes from b = -v f(x - v t, y, z) to any b."""
    n = b.shape[0]
    zero = b * 0
    one = zero + 1
    if MUT == "num-static":
        db = [zero] + list(db[1:])
        ddb = [[zero if (i == 0 or j == 0) else ddb[i][j] for j in range(4)] for i in range(4)]
    keep = 0 if MUT == "num-drop-w2" else 1  # g_tt = -1 + keep * b^2
    dt = object if b.dtype == object else float
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
    det2 = g[:, 0, 0] * g[:, 1, 1] - g[:, 0, 1]**2
    gi = np.zeros_like(g)
    if dt is object:
        gi[...] = zero[:, None, None]
    gi[:, 0, 0], gi[:, 1, 1] = g[:, 1, 1] / det2, g[:, 0, 0] / det2
    gi[:, 0, 1] = gi[:, 1, 0] = -g[:, 0, 1] / det2
    gi[:, 2, 2] = gi[:, 3, 3] = one
    dgi = -np.einsum("nik,nakl,nlj->naij", gi, dg, gi)
    T = np.einsum("nilj->nlij", dg) + np.einsum("njli->nlij", dg) - dg
    dT = np.einsum("nmilj->nmlij", ddg) + np.einsum("nmjli->nmlij", ddg) - ddg
    Gam = np.einsum("nkl,nlij->nkij", gi, T) / 2
    dGam = (np.einsum("nmkl,nlij->nmkij", dgi, T) + np.einsum("nkl,nmlij->nmkij", gi, dT)) / 2
    Ric = (np.einsum("nkkij->nij", dGam) - np.einsum("njkik->nij", dGam)
           + np.einsum("nkkl,nlij->nij", Gam, Gam) - np.einsum("nkjl,nlik->nij", Gam, Gam))
    Rs = np.einsum("nij,nij->n", gi, Ric)
    G = Ric - Rs[:, None, None] * g / 2
    return np.einsum("nik,nkj->nij", gi, G), np.einsum("nik,nkj->nij", gi, Ric)


class RouteN:
    """Route N' for the explicit shift b_expr(t, x, y, z): sympy differentiates b_expr, lambdify evaluates."""

    def __init__(self, b_expr):
        co = (T_, XX, YY, ZZ)
        self.exprs = ([b_expr] + [sp.diff(b_expr, c) for c in co]
                      + [sp.diff(b_expr, co[i], co[j]) for i in range(4) for j in range(i, 4)])
        self.fn = sp.lambdify(co, self.exprs, modules="numpy", cse=True)
        self._fm = None

    @staticmethod
    def _unpack(vals, ones):
        vals = [x * ones for x in vals]
        ddb = [[None] * 4 for _ in range(4)]
        k = 5
        for i in range(4):
            for j in range(i, 4):
                ddb[i][j] = ddb[j][i] = vals[k]
                k += 1
        return vals[0], vals[1:5], ddb

    def mixed(self, t0, x, y, z):
        b, db, ddb = self._unpack(self.fn(t0, x, y, z), np.ones_like(x))
        return einstein_b(b, db, ddb)

    def mixed_mp(self, t0, x, y, z):
        """At one point, mpmath scalars t0, x, y, z, in whatever precision mp.mp.dps is set to."""
        if self._fm is None:
            self._fm = sp.lambdify((T_, XX, YY, ZZ), self.exprs, modules="mpmath", cse=True)
        ones = np.array([mp.mpf(1)], dtype=object)
        b, db, ddb = self._unpack(self._fm(t0, x, y, z), ones)
        return einstein_b(b, db, ddb)


# ---------------------------------------------------------------------------------------------------------------
# Closed form against route N' at grid points
# ---------------------------------------------------------------------------------------------------------------
def compare(cf, rn, name, t0, sig, Rb, r, c, label, seed):
    """Sign of D (cf) against the eigenvalue type of route N''s G^a_b at the points (r, cos theta), random azimuths. Points
    where float64 disagrees, or |D| < 1e-10 scale^6, and the three nearest D = 0 are rechecked at 50 digits (both routes in
    mpmath), up to 1000 per grid. A point is D = 0 only by an exact factor: at v = 0 exactly, D = al^6 f'^6 s^4 / 16 (the
    `structure` check), and s^4 = (1 - c^2)^2 is evaluated in exact rationals from the float c. Then D = 0, which is not Type
    IV (a real repeated root, the rank test decides among I, II and III), and route N' must show no complex eigenvalue. No
    tolerance decides D = 0. Returns (points, float64 disagreements, unresolved, D, scale)."""
    v, al, u, xs0 = motion(name, t0)
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0, 2 * np.pi, r.size)
    s = np.sqrt(np.clip(1 - c**2, 0, None))
    X, y, z = r * c, r * s * np.cos(phi), r * s * np.sin(phi)
    D, scale, _ = cf.D(v, al, u, sig, Rb, r, c)
    ivC = D < 0
    ivN = np.zeros(r.size, dtype=bool)
    for k in range(0, r.size, 20000):
        Gm, _ = rn.mixed(t0, xs0 + X[k:k + 20000], y[k:k + 20000], z[k:k + 20000])
        ivN[k:k + 20000] = ht.numeric_type(Gm)
    near = np.abs(D) < 1e-10 * scale**6
    closest = np.argsort(np.abs(D) / scale**6)[:3]
    redo = np.union1d(closest, np.where((ivC != ivN) | near)[0])
    unresolved = [None] * max(0, redo.size - 1000)
    mp.mp.dps = 50
    vm, alm, um, xsm = motion_mp(name, t0)
    t0m = mp.mpf(t0)
    for i in redo[:1000]:
        rm, cm, ph = mp.mpf(float(r[i])), mp.mpf(float(c[i])), mp.mpf(float(phi[i]))
        sm = mp.sqrt(1 - cm**2)
        Dm, _, _ = cf.D(vm, alm, um, sig, Rb, rm, cm, mpm=True)
        exact_zero = v == 0 and (1 - sp.Rational(float(c[i]))**2)**2 == 0
        ivCm = bool(Dm < 0) and not exact_zero
        Gm, _ = rn.mixed_mp(t0m, xsm + rm * cm, rm * sm * mp.cos(ph), rm * sm * mp.sin(ph))
        ev = mp.eig(mp.matrix(Gm[0].tolist()))[0]
        sc = max(abs(e) for e in ev)
        ivNm = max(abs(mp.im(e)) for e in ev) > mp.mpf(10)**-30 * sc
        if ivCm != ivNm:
            unresolved.append((float(r[i]), float(c[i]), float(Dm), float(max(abs(mp.im(e)) for e in ev))))
    mp.mp.dps = 15
    nd = int((ivC != ivN).sum())
    print(f"  {label}: {r.size} points, Type IV by D < 0 at {int(ivC.sum())}, by eigenvalues at {int(ivN.sum())}; "
          f"float64 disagreements {nd}, near D = 0 {int(near.sum())}, rechecked at 50 digits {min(redo.size, 1000)}, "
          f"unresolved {len(unresolved)}", flush=True)
    return r.size, nd, unresolved, D, scale


# ---------------------------------------------------------------------------------------------------------------
# Exact points
# ---------------------------------------------------------------------------------------------------------------
QR = sp.Rational
JET_PTS = [(QR(1, 3), QR(-5, 2), QR(7, 4), QR(9, 10)), (QR(4, 5), QR(-1, 3), QR(-2, 9), QR(6, 5)),
           (QR(1, 9), QR(-3), QR(11, 2), QR(1))]  # (f, f', f'', r), as he_type's screen


def cs_from_tan(t):
    """(cos, sin) of an angle with tan(angle / 2) = t."""
    return (1 - t**2) / (1 + t**2), 2 * t / (1 + t**2)


def entries_vs_E(blk):
    """Route D's block at z = 0 against route E, at exact points with explicit rational v_s(t), x_s(t), f(r)."""
    tt, xx = ps.TX
    rr = sp.Symbol("rr", positive=True)
    f_ex, vs_ex, xs_ex = (2 + rr**2) / (3 + rr**4), (3 + 2 * tt) / (2 + tt**2), tt**3 / 5 - tt**2 / 3 + tt / 2
    pyth = [(QR(3, 5), QR(4, 5)), (QR(-4, 5), QR(3, 5)), (QR(5, 13), QR(12, 13)), (QR(-12, 13), QR(5, 13)),
            (QR(8, 17), QR(15, 17)), (QR(-15, 17), QR(8, 17)), (QR(7, 25), QR(24, 25)), (QR(-20, 29), QR(21, 29))]
    ks = [1, QR(3, 2), 2, QR(1, 2), QR(4, 3), QR(5, 3), QR(3, 4), 1]
    tvs = [QR(1, 3), QR(-2, 5), QR(3, 2), QR(1, 7), QR(-5, 4), 2, QR(-1, 2), QR(3, 4)]
    bad = 0
    stray = blk.free_symbols - set(JETS)
    for (cc, ss), k, tv in zip(pyth, ks, tvs):
        Xv, yv, rv = k * cc, k * ss, sp.Integer(1) * k
        sub = {V_: vs_ex.subs(tt, tv), AL: sp.diff(vs_ex, tt).subs(tt, tv), U: sp.diff(xs_ex, tt).subs(tt, tv),
               F0: f_ex.subs(rr, rv), F1: sp.diff(f_ex, rr).subs(rr, rv), F2: sp.diff(f_ex, rr, 2).subs(rr, rv),
               r_: rv, C_: cc, S_: ss}
        mine = blk.subs(sub).applyfunc(sp.cancel)
        b_ex = -vs_ex * f_ex.subs(rr, sp.sqrt((xx - xs_ex)**2 + ps.y**2 + ps.z**2))
        ex = ps.ricci_at_point(b_ex, (tv, xs_ex.subs(tt, tv) + Xv, yv, 0), False)
        bad += any(sp.simplify(e) != 0 for e in (mine - ex))
    return bad, len(pyth), stray


# ---------------------------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------------------------
def main():
    global MUT, QUICK
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    MUT = opts.get("mutate")
    QUICK = "quick" in opts
    stop = int(opts.get("stop", 99))
    if MUT is not None and MUT not in MUTATIONS:
        sys.exit(f"unknown mutation {MUT!r}. Known: {', '.join(MUTATIONS)}")
    if MUT in ("disc-sign", "no-rank"):
        ht.MUT = MUT
    t0 = time.time()
    print(f"he_accel.py, issue #73. sympy {sp.__version__}, numpy {np.__version__}, mpmath {mp.__version__}, python "
          f"{sys.version.split()[0]}. mutation: {MUT or 'none'}" + (f" ({MUTATIONS[MUT][1]})" if MUT else "")
          + f". grids: {'quick (every 2nd index of the 240 x 240 grids, t = +-0.3)' if QUICK else 'full'}")
    le_r = ht.wall_radii(*LE)
    ro_r = (5 - 3 / 4, 5 + 3 / 4)

    # 1. Planted classifications through he_type's classifier.
    print("\n== 1. Planted truth: #64's nine tensors of known type through he_type_exact ==")
    bad = []
    for name, T4, want in ht.planted_cases()[:9]:
        got, D = ht.he_type_exact(T4)
        bad += [] if got == want else [name]
        print(f"  {name:>46}: D = {D!s:>14}, type {got:>3}, expected {want}")
    check("planted_cls", not bad, f"every planted tensor gets its type from D and the rank test; wrong: {bad}")

    # 2. Route N' planted: the constant-speed case against he_type's route N, and the accelerating case against route E.
    print("\n== 2. Route N': planted against he_type's route N (constant speed) and against route E (accelerating) ==")
    V0 = sp.Rational(3, 4)
    rs_c = sp.sqrt((XX - V0 * T_)**2 + YY**2 + ZZ**2)
    rn_c = RouteN(-V0 * tanh_f(rs_c, *LE))
    rng = np.random.default_rng(73)
    Xp, yp, zp = rng.uniform(-1.3, 1.3, 40), rng.uniform(-1.3, 1.3, 40), rng.uniform(-1.3, 1.3, 40)
    tp = 0.35
    Gn, _ = rn_c.mixed(tp, Xp + float(V0) * tp, yp, zp)
    rp = np.sqrt(Xp**2 + yp**2 + zp**2)
    f_, f1_, f2_ = ht.tanh_jets(rp, *LE)
    Gh, _, _ = ht.einstein_N(Xp, yp, zp, rp, float(V0), f_, f1_, f2_)
    err = float(np.max(np.abs(Gn - Gh)) / np.max(np.abs(Gh)))
    check("routeN_planted", err < 1e-9, f"route N' for b = -v f(x - v t, y, z), v = 3/4, equals ht.einstein_N's mixed G at 40 "
          f"random points, t = {tp} (relative error {err:.1e})")

    tt, xx = ps.TX
    rr = sp.Symbol("rr", positive=True)
    f_ex, vs_ex, xs_ex = (2 + rr**2) / (3 + rr**4), (3 + 2 * tt) / (2 + tt**2), tt**3 / 5 - tt**2 / 3 + tt / 2
    b_ex = -vs_ex * f_ex.subs(rr, sp.sqrt((xx - xs_ex)**2 + ps.y**2 + ps.z**2))
    rn_e = RouteN(b_ex.subs({tt: T_, xx: XX, ps.y: YY, ps.z: ZZ}))
    worst = 0.0
    pts = [(12, 3, 4, 13), (-12, 4, 3, 13), (4, 3, 0, 5), (8, 9, 12, 17), (0, 7, 24, 25), (-15, 12, 16, 25)]
    times = [QR(1, 3), QR(-2, 5), QR(3, 2), QR(1, 7), QR(-5, 4), 2]
    eta = sp.diag(-1, 1, 1, 1)
    for (x0, y0, z0, r0), tv in zip(pts, times):
        Xv, yv, zv = [QR(c_, 7) for c_ in (x0, y0, z0)]
        exr = ps.ricci_at_point(b_ex, (tv, xs_ex.subs(tt, tv) + Xv, yv, zv), False)
        Me = np.array((eta * exr).tolist(), dtype=float)
        xl = float(xs_ex.subs(tt, tv) + Xv)
        _, Rm = rn_e.mixed(float(tv), np.array([xl]), np.array([float(yv)]), np.array([float(zv)]))
        pe, pn = np.poly(Me), np.poly(Rm[0])
        worst = max(worst, float(np.max(np.abs(pe - pn) / np.maximum(np.abs(pe), 1e-3 * np.abs(pe).max()))))
    check("routeN_vs_E", worst < 1e-8, f"route N''s mixed Ricci tensor has the characteristic polynomial of route E's "
          f"(exact rationals, explicit v_s, x_s, f) at 6 points with t != 0 (worst relative coefficient error {worst:.1e})")

    # 3. Route D: the block, its entries against route E, and the identities.
    print("\n== 3. Route D: the frame Ricci tensor of the bubble at z = 0, and the block (E_0, E_1, e_rho) ==")
    Ra = route_D(MUT)
    jets_in = sorted(str(x) for x in Ra.free_symbols if x in ps.VS or x in ps.XS or x in ps.F)
    print(f"  jets of v_s, x_s and f in the frame Ricci tensor: {jets_in}")
    blk = block_z0(Ra)
    names = {0: "0", 1: "1", 2: "rho", 3: "varphi"}
    for a in range(4):
        for b in range(a, 4):
            print(f"  R_({names[a]},{names[b]}) = {sp.factor(blk[a, b])}")
    nbad, npts, stray = entries_vs_E(blk)
    check("entries_E", nbad == 0 and not stray, f"route D's block at z = 0 equals route E at {npts - nbad} of {npts} exact "
          f"points (explicit v_s = (3 + 2t)/(2 + t^2), x_s = t^3/5 - t^2/3 + t/2, f = (2 + r^2)/(3 + r^4)); symbols other "
          f"than (v, al, u, f, f', f'', r, c, s) in the block: {sorted(map(str, stray))}")
    dec = [reduce_s(blk[a, 3]) for a in range(4)] + [reduce_s(blk[0, 0] + blk[1, 1] + 2 * blk[2, 2])]
    check("decouple", all(e == 0 for e in dec), f"(I1) R(e_varphi, .) and (I2) R_00 + R_11 + 2 R_rho,rho reduce to {dec}")

    # 4. The time lemma: the entries in closed form, and acceleration entering only through two derivatives of b.
    print("\n== 4. The entries in closed form, and the time lemma ==")
    P = closed_pieces()
    for k, e in P.items():
        print(f"  {k} = {e}")
    cl = closed_block(P)
    ok_cl = (blk - cl).applyfunc(reduce_s) == sp.zeros(4, 4)
    dTx = P["Tx"] - P["Tx"].subs({AL: 0, U: V_})
    dTr = P["Tr"] - P["Tr"].subs({AL: 0, U: V_})
    sg = -1 if MUT == "time-lemma-sign" else 1
    E = lambda i, j: sp.Matrix(4, 4, lambda a, b: 1 if (a, b) == (i, j) else 0)  # noqa: E731
    claim = dTx * (E(0, 0) - E(1, 1)) - sg * dTr / 2 * (E(1, 2) + E(2, 1))
    ok_tl = (blk - blk.subs({AL: 0, U: V_}) - claim).applyfunc(reduce_s) == sp.zeros(4, 4)
    print(f"  dTx = {sp.factor(dTx)}\n  dTr = {sp.factor(dTr)}")
    check("time_lemma", ok_cl and ok_tl, f"the block equals R_00 = Tx - A, R_01 = j1, R_0rho = jr, R_11 = A - q - Tx, "
          f"R_1rho = B - Tr/2, R_rhorho = q/2 [{ok_cl}], and R(v, u, al) - R(v, v, 0) = DeltaTx (E_00 - E_11) - DeltaTr/2 "
          f"(E_1rho + E_rho1) [{ok_tl}]")

    # 5. Structure of D.
    print("\n== 5. D: the shifted cubic, the weights and the powers of v, al and u ==")
    j1, jr, q, h, p, lam = sp.symbols("j1 jr q h p lam")
    Dc = sp.expand(D_shift(j1, jr, q, h, p))
    print(f"  D = 18 c2 c1 c0 - 4 c2^3 c0 + c2^2 c1^2 - 4 c1^3 - 27 c0^2, c2 = q - p, c1 = j1^2 + jr^2 - p q - h^2, "
          f"c0 = 2 j1 jr h - p j1^2 + q jr^2: {len(Dc.args)} terms, degree {sp.degree(Dc, p)} in p")
    R00, R01, R0r, Rrr, R1r = sp.symbols("R00 R01 R0r Rrr R1r")
    Rg = sp.Matrix([[R00, R01, R0r, 0], [R01, -R00 - 2 * Rrr, R1r, 0], [R0r, R1r, Rrr, 0], [0, 0, 0, 0]])
    ok_shift = sp.expand(D_direct(Rg) - Dc.subs(dict(zip((j1, jr, q, h, p), shifted_args(Rg))), simultaneous=True)) == 0
    wl = sp.Symbol("wl", positive=True)
    Dw = sp.expand(Dc.subs({j1: wl * j1, jr: wl * jr, q: wl**2 * q, h: wl**2 * h, p: wl**2 * p}, simultaneous=True))
    wpow = sorted(m[0] for m in sp.Poly(Dw, wl).monoms())
    ok_w6 = sp.expand(Dw.coeff(wl, 6) + 4 * (j1**2 + jr**2)**3) == 0
    ok_w12 = sp.expand(Dw.coeff(wl, 12) - (p * q + h**2)**2 * ((p + q)**2 + 4 * h**2)) == 0
    for k in (8, 10):
        print(f"  weight {k} part: {sp.factor(Dw.coeff(wl, k))}")
    alw = 1 if MUT == "weight-break" else 2
    ok_bw = True
    for a in range(3):
        for b in range(a, 3):
            w_ = 1 if (a == 0 and b > 0) else 2
            sc = reduce_s(blk[a, b].subs({V_: wl * V_, AL: wl**alw * AL, U: wl * U}, simultaneous=True))
            ok_bw &= sp.expand(sc - wl**w_ * blk[a, b]) == 0
    jb, jrb, qb, hb, pb = shifted_args(blk)
    Db = D_shift(jb, jrb, qb, hb, pb)
    num, den = sp.fraction(sp.together(Db))
    Pn = sp.Poly(sp.expand(num), V_, AL, U)
    mons = Pn.monoms()
    wts = sorted({i + k + 2 * j for (i, j, k) in mons})
    print(f"  D of the bubble: 64 r^6 D is a polynomial in (v, al, u) with {len(mons)} monomials, denominator {den}; "
          f"weights i + k + 2j of v^i al^j u^k: {wts}; highest power of v {max(m[0] for m in mons)}, of al "
          f"{max(m[1] for m in mons)}, of u {max(m[2] for m in mons)}")
    ok_pw = (wts == [6, 8, 10, 12] and all(k <= i for (i, j, k) in mons) and max(m[1] for m in mons) == 6)
    lead = reduce_s(Pn.coeff_monomial(AL**6) / den)
    D0 = reduce_s(sp.together(Db.subs(V_, 0)))
    ok_lead = (sp.cancel(sp.together(lead - F1**6 * (1 - C_**2)**2 / 16)) == 0
               and sp.cancel(sp.together(D0 - AL**6 * F1**6 * (1 - C_**2)**2 / 16)) == 0)
    print(f"  al^6 coefficient of D = f'^6 s^4/16 [{ok_lead}]; D at v = 0 equals al^6 f'^6 s^4/16 [{ok_lead}]")
    print(f"  D at (al, u) = (0, v) equals he_type's D, entry by entry, in `const_limit`")
    check("structure", ok_shift and wpow == [6, 8, 10, 12] and ok_w6 and ok_w12 and ok_bw and ok_pw and ok_lead,
          f"D of the block = D(j1, jr, q, h, p) of the shifted cubic [{ok_shift}]; the scaling (v, u, al) -> (w v, w u, w^2 al) "
          f"gives R_0i w, the other entries w^2 [{ok_bw}] and D the weights {wpow} with w^6 part -4|j|^6 [{ok_w6}] and w^12 "
          f"part (pq + h^2)^2 ((p + q)^2 + 4h^2) [{ok_w12}]; the bubble's D has weights {wts}, u only with v, al^6 at most "
          f"[{ok_pw}]; al^6 coefficient f'^6 s^4/16 and D(v = 0) = al^6 f'^6 s^4/16 [{ok_lead}]")

    # 5b. Parity: (al, cos theta) -> (-al, -cos theta) leaves D unchanged, for every v and u.
    flip = {AL: -AL} if MUT == "parity-break" else {AL: -AL, C_: -C_}
    ok_par = reduce_s(sp.expand(num) - sp.expand(num.subs(flip, simultaneous=True))) == 0
    rr_, cc_, _ = pgrid(*le_r)
    v, al, u, _ = motion("P1", 0.3)
    cf = ClosedFormAccel(blk)
    D1, scl1, _ = cf._D(v, al, u, *LE, rr_, cc_, False)
    D2, _, _ = cf._D(v, -al, u, *LE, rr_, cc_ if MUT == "parity-break" else -cc_, False)
    par_num = float(np.max(np.abs(D1 - D2) / scl1**6))
    check("parity", ok_par and par_num < 1e-9, f"64 r^6 D(v, -al, u; -c) - 64 r^6 D(v, al, u; c) reduces to 0 as a polynomial in "
          f"(v, al, u, f, f', f'', r, c) [{ok_par}], so the type at (theta, al) is the type at (pi - theta, -al) and every "
          f"wall fraction is even in al; on Le's grid at P1 t = 0.3 the largest |D(-al, -c) - D(al, c)| / scale^6 is "
          f"{par_num:.1e} (< 1e-9)")

    # 6. The constant-speed limit against #64.
    print("\n== 6. Constant-speed limit (al = 0, u = v) against #64: entries, maps and fractions ==")
    Rc = ht.route_C()
    ok_ent = (blk.subs({AL: 0, U: V_}).applyfunc(reduce_s) - Rc).applyfunc(reduce_s) == sp.zeros(4, 4)
    print(f"  the block at (al, u) = (0, v) equals ht.route_C()'s frame Ricci tensor entry by entry: {ok_ent}")
    cfh = ht.ClosedForm(Rc)
    mp.mp.dps = 15
    ndiff = ntot = 0
    unres_c = []
    cases = [(LE, v, "Le", le_r) for v in (0.5, 1.0, 2.0, 2.5)] + [(RODAL, 1.0, "Rodal", ro_r)]
    for (sb, v, nm, (a1, a2)) in cases:
        rr_, cc_, _ = pgrid(a1, a2)
        Da, _, _ = cf.D(v, 0.0, v, *sb, rr_, cc_)
        Dh, _, _ = cfh.D(v, *sb, rr_, cc_)
        dif = np.where((Da < 0) != (Dh < 0))[0]
        for i in dif[:100]:
            mp.mp.dps = 50
            da = cf.D(mp.mpf(v), mp.mpf(0), mp.mpf(v), *sb, mp.mpf(float(rr_[i])), mp.mpf(float(cc_[i])), mpm=True)[0]
            dh = cfh.D(mp.mpf(v), *sb, mp.mpf(float(rr_[i])), mp.mpf(float(cc_[i])), mpm=True)[0]
            mp.mp.dps = 15
            if (da < 0) != (dh < 0):
                unres_c.append((nm, v, float(rr_[i]), float(cc_[i])))
        ndiff += dif.size
        ntot += rr_.size
    print(f"  type maps (D < 0) of the two closed forms at {ntot} points of Le's wall at v = 0.5, 1, 2, 2.5 and Rodal's "
          f"band at v = 1: float64 sign differences {ndiff}, unresolved at 50 digits {len(unres_c)}")
    # P1 at t = 0 is the constant-speed case v = 1: al = 0 and u = v = 1, exactly.
    v1, al1, u1, _ = motion("P1", 0.0)
    rr_, cc_, _ = pgrid(*le_r)
    Da, _, _ = cf.D(v1, al1, u1, *LE, rr_, cc_)
    Dh, _, _ = cfh.D(1.0, *LE, rr_, cc_)
    p1_same = (v1, al1, u1) == (1.0, 0.0, 1.0) and bool(((Da < 0) == (Dh < 0)).all())
    print(f"  P1 at t = 0: (v, al, u) = ({v1}, {al1}, {u1}); its map on Le's wall equals the constant-speed map at v = 1: "
          f"{p1_same}")
    txt = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "he_type_output.txt")).read()
    pub = {float(m[0]): float(m[1]) for m in re.findall(r"^\s*\| ([0-9.]+) \| ([0-9.]+) \| [0-9.]+ \| [0-9.]+ \| [0-9.]+ \| "
                                                         r"[0-9.]+ \| [0-9.]+ \|$", txt, re.M)}
    fr_ok = True
    fr_msg = []
    if not QUICK:
        for v in (0.5, 1.0, 2.0, 2.5):
            mine = 100 * fraction_I(cf, v, 0.0, v, *LE, *le_r, 1000)
            fr_ok &= v in pub and abs(mine - pub[v]) < 1.01e-3
            fr_msg.append(f"v = {v}: {mine:.3f} % against he_type_output.txt's {pub.get(v)} %")
        print("  Type I fractions of Le's wall, polar 1000^2, (al, u) = (0, v): " + "; ".join(fr_msg))
    check("const_limit", ok_ent and not unres_c and p1_same and fr_ok,
          f"at (al, u) = (0, v) the block equals #64's entry by entry [{ok_ent}]; its type map equals #64's at {ntot} grid "
          f"points ({ndiff} float64 sign differences, {len(unres_c)} unresolved at 50 digits); P1 at t = 0 gives the v = 1 "
          f"map [{p1_same}]; Le's Type I fractions equal he_type_output.txt's to 3 digits "
          f"[{'not run in quick mode' if QUICK else fr_ok}]")

    # 7. The axis and the rest slice.
    print("\n== 7. The axis (Type IV for every al and u) and the instant of rest (Type I) ==")
    ax = [blk.subs({C_: cs_, S_: 0}) for cs_ in (1, -1)]
    sgn_ax = -1 if MUT == "axis-flip" else 1
    ok_axs = True
    for B_ in ax:
        Dax = sp.expand(D_direct(B_))
        claim_ax = sp.expand(-4 * sgn_ax * B_[0, 1]**2 * (B_[0, 0]**2 + B_[0, 1]**2)**2)
        ok_axs &= sp.expand(Dax - claim_ax) == 0
    ex_ok = True
    for vv, aa, uu in ((QR(1, 2), QR(3, 2), QR(1, 3)), (QR(-1, 3), QR(-2), QR(1, 5))):
        for (f, f1, f2, rv) in JET_PTS:
            for cs_ in (1, -1):
                R4 = exact_R4(blk, {V_: vv, AL: aa, U: uu, F0: f, F1: f1, F2: f2, r_: rv, C_: cs_, S_: 0})
                ex_ok &= ht.he_type_exact(R4)[0] == "IV"
    n_ax = 0
    ax_bad = []
    for (name, t_) in slices():
        v, al, u, _ = motion(name, t_)
        if v == 0:
            continue
        for sb, (a1, a2) in ((LE, le_r), (RODAL, ro_r)):
            rax = np.linspace(a1, a2, 301)
            for cs_ in (1.0, -1.0):
                D, scl, _ = cf.D(v, al, u, *sb, rax, cs_ * np.ones_like(rax))
                flag = np.where((D >= 0) | (np.abs(D) < 1e-10 * scl**6))[0]
                n_ax += rax.size
                mp.mp.dps = 50
                vm, alm, um, _ = motion_mp(name, t_)
                for i in flag[:200]:
                    Dm = cf.D(vm, alm, um, *sb, mp.mpf(float(rax[i])), mp.mpf(cs_), mpm=True)[0]
                    if Dm >= 0:
                        ax_bad.append((name, t_, float(rax[i])))
                ax_bad += [(name, t_, "unrechecked")] * max(0, flag.size - 200)
                mp.mp.dps = 15
    check("axis", ok_axs and ex_ok and not ax_bad, f"on c = +-1, D = -4 R_01^2 (R_00^2 + R_01^2)^2 as a polynomial identity in "
          f"(v, al, u, f, f', f'', r) [{ok_axs}], so Type IV wherever v f' != 0 for every al and u; he_type_exact gives IV at "
          f"{2 * 2 * len(JET_PTS)} exact axis points [{ex_ok}]; D < 0 at all {n_ax} axis points of the slices with v != 0, "
          f"both walls, after 50-digit rechecks of the flagged ones (failures {ax_bad[:3]})")
    Br = blk.subs(V_, 0)
    ok_rs = reduce_s(Br[0, 1]) == 0 and reduce_s(Br[0, 2]) == 0
    t_rest = 0.05 if MUT == "rest-v" else 0.0
    v_rest = QR(1, 10) if MUT == "rest-v" else 0
    ex_ok = True
    for aa in (QR(2), QR(-3, 2)):
        for (f, f1, f2, rv) in JET_PTS + [(QR(1, 2), 0, QR(1, 3), QR(1))]:
            for cs_, ss_ in (cs_from_tan(QR(2, 7)), (1, 0), (0, 1), (-1, 0)):
                R4 = exact_R4(blk, {V_: v_rest, AL: aa, U: QR(1, 4), F0: f, F1: f1, F2: f2, r_: rv, C_: cs_, S_: ss_})
                ex_ok &= ht.he_type_exact(R4)[0] == "I"
    v, al, u, _ = motion("P2", t_rest)
    nr = 0
    rest_min = 0.0
    for sb, (a1, a2) in ((LE, le_r), (RODAL, ro_r)):
        rr_, cc_, _ = pgrid(a1, a2)
        D, scl, _ = cf.D(v, al, u, *sb, rr_, cc_)
        rest_min = min(rest_min, float((D / scl**6).min()))
        nr += rr_.size
    check("rest", ok_rs and ex_ok and rest_min >= -1e-12, f"at v = 0 the momentum entries R_01, R_0rho vanish identically "
          f"[{ok_rs}]; he_type_exact gives I at {2 * (len(JET_PTS) + 1) * 4} exact points including the axis and D = 0 [{ex_ok}]; "
          f"on {nr} grid points of both walls at P2, t = {t_rest} (v = {v:.4f}, al = {al:.4f}, u = {u:.4f}) the smallest "
          f"D / scale^6 is {rest_min:.2e} (>= -1e-12)")

    if stop < 8:
        return finish(t0)

    # 8. The closed form against route N' at every grid point, for both profiles at every time.
    print("\n== 8. Closed form (sign of D) against route N''s eigenvalues at every grid point ==")
    total = nd_tot = 0
    unres = []
    drop_rows = []
    rns = {}
    seed = 1
    for (name, t_) in slices():
        v, al, u, xs0 = motion(name, t_)
        print(f"  {name}, t = {t_}: v = {v:.6f}, al = {al:.6f}, u = {u:.6f}, x_s = {xs0:.6f}")
        vs, xs = PROFILES[name]
        for wl_, sb, (a1, a2) in (("Le", LE, le_r), ("Rodal", RODAL, ro_r)):
            if (name, wl_) not in rns:
                rs = sp.sqrt((XX - xs.subs(T_, T_))**2 + YY**2 + ZZ**2)
                rns[(name, wl_)] = RouteN(-vs * tanh_f(rs, *sb))
            rn = rns[(name, wl_)]
            rr_, cc_, ww_ = pgrid(a1, a2)
            n, nd, un, Dfull, _ = compare(cf, rn, name, t_, *sb, rr_, cc_, f"{wl_}'s wall, polar grid", seed)
            seed += 1
            total, nd_tot, unres = total + n, nd_tot + nd, unres + un
            Ddrop, _, _ = cf.D_dropped(v, al, u, *sb, rr_, cc_)
            drop_rows.append((name, t_, wl_, v, al, u, int(((Dfull < 0) != (Ddrop < 0)).sum()), rr_.size,
                              float((ww_ * (Dfull > 0)).sum() / ww_.sum()), float((ww_ * (Ddrop > 0)).sum() / ww_.sum())))
        rax = np.linspace(5 - 3 / 4, 5 + 3 / 4, 301)
        n, nd, un, _, _ = compare(cf, rns[(name, "Rodal")], name, t_, *RODAL, np.concatenate([rax, rax]),
                                  np.concatenate([np.ones(301), -np.ones(301)]), "Rodal's axis, theta = 0 and pi", seed)
        seed += 1
        total, nd_tot, unres = total + n, nd_tot + nd, unres + un
    if unres:
        print(f"  unresolved at 50 digits (r, cos theta, D, max|Im lam|), first 5: {[x for x in unres if x is not None][:5]}")
    check("closed_vs_numeric", not unres, f"the sign of D gives the eigenvalues' type at all {total} grid points "
          f"({nd_tot} float64 disagreements; {len(unres)} unresolved after up to 1000 rechecks per grid at 50 digits)")

    # 9. Dropping the acceleration terms changes the map.
    print("\n== 9. The map with the acceleration terms dropped (al = 0, u = v at the same v), against the full one ==")
    print("  | profile | t | wall | v | al | u | points whose type differs | of | Type I % full | Type I % dropped |")
    for row in drop_rows:
        print("  | {} | {} | {} | {:.4f} | {:.4f} | {:.4f} | {} | {} | {:.3f} | {:.3f} |".format(
            row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], 100 * row[8], 100 * row[9]))
    need = [row for row in drop_rows if row[3] != 0 and (row[4] != 0 or abs(row[5] - row[3]) > 1e-12)]
    same0 = [row for row in drop_rows if row[3] == 1.0 and row[4] == 0.0 and row[5] == 1.0]
    ok_diff = all(row[6] >= 1 for row in need) and all(row[6] == 0 for row in same0)
    ok_pp = all(any(abs(row[8] - row[9]) >= 0.01 for row in need if row[0] == pf and row[2] == "Le")
                for pf in ("P1", "P2"))
    check("drop_detected", ok_diff and ok_pp, f"at all {len(need)} (profile, time, wall) slices with v != 0 and (al != 0 or "
          f"u != v) the type differs at >= 1 grid point [{ok_diff}], the Le Type I fraction differs by >= 1 percentage point "
          f"at some slice of each profile [{ok_pp}], and at P1 t = 0 (al = 0, u = v) no point differs "
          f"[{all(row[6] == 0 for row in same0)}]")

    # 10. Tables of the map (numerical evidence, not checks).
    if MUT is None:
        tables(cf, le_r, ro_r)

    return finish(t0)


def finish(t0):
    failed = [nm for nm, ok in RESULTS if not ok]
    rss = [ln for ln in open("/proc/self/status") if ln.startswith("VmHWM")]
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - t0:.0f} s, "
          f"{rss[0].split(':')[1].strip() if rss else '?'} peak")
    if MUT is not None:
        want = MUTATIONS[MUT][0]
        print(f"control {MUT}: target checks {want}, fired {[n for n in want if n in failed]}")
    return failed


def tables(cf, le_r, ro_r):
    """Type I fractions and the equatorial annulus as functions of the acceleration (numerical evidence)."""
    print("\n== 10. The map as the acceleration changes (numerical evidence, polar 1000^2, u = v) ==")
    als = [-32, -16, -8, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 8, 16, 32]
    print("  Type I % of Le's wall (sigma = 8, R_b = 1, f in [0.1, 0.9]) at speed v and acceleration al:")
    print("  | v \\ al | " + " | ".join(f"{a:g}" for a in als) + " |")
    print("  |---|" + "---|" * len(als))
    for v in (0.5, 1.0, 2.0):
        print(f"  | {v} | " + " | ".join(f"{100 * fraction_I(cf, v, a, v, *LE, *le_r, 1000):.2f}" for a in als) + " |", flush=True)
    print("  Type I % of Rodal's band (sigma = 4, rho = 5, r in [4.25, 5.75]) at v = 1:")
    print("  | al | " + " | ".join(f"{a:g}" for a in als) + " |")
    print("  |---|" + "---|" * len(als))
    print("  | Type I % | " + " | ".join(f"{100 * fraction_I(cf, 1.0, a, 1.0, *RODAL, *ro_r, 1000):.2f}" for a in als) + " |",
          flush=True)
    print("  Type I % of Le's wall at v = 0.5 for u != v (al = 0): u = 0, 0.25, 0.5, 1, 2:  " + ", ".join(
        f"{100 * fraction_I(cf, 0.5, 0.0, u, *LE, *le_r, 1000):.2f}" for u in (0.0, 0.25, 0.5, 1.0, 2.0)))
    print("  Type I % of Le's wall along Alcubierre's journey P1 (u = v), full against dropped:")
    print("  | t | v | al | Type I % full | Type I % dropped |")
    print("  |---|---|---|---|---|")
    for t_ in (-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9):
        v, al, u, _ = motion("P1", t_)
        fa = fraction_I(cf, v, al, u, *LE, *le_r, 1000)
        fd = fraction_I(cf, v, 0.0, v, *LE, *le_r, 1000)
        print(f"  | {t_} | {v:.4f} | {al:.4f} | {100 * fa:.3f} | {100 * fd:.3f} |", flush=True)
    print("  Equatorial cut of Rodal's band (theta = pi/2), v = 1, u = v: the radii where D changes sign, so the Type I "
          "intervals, for each al:")
    rr = np.linspace(5 - 3 / 4, 5 + 3 / 4, 150001)
    for a in (-4.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 4.0):
        Deq, _, _ = cf.D(1.0, a, 1.0, *RODAL, rr, np.zeros_like(rr))
        typeI = Deq > 0
        edges = np.flatnonzero(np.diff(typeI.astype(int)))
        fD = lambda x: float(cf.D(1.0, a, 1.0, *RODAL, np.array([x]), np.array([0.0]))[0][0])  # noqa: E731
        segs = [brentq(fD, rr[e], rr[e + 1], xtol=1e-12) for e in edges]
        print(f"    al = {a:+.1f}: Type I at {int(typeI.sum())} of {rr.size} radii, D changes sign at "
              f"{[round(x, 4) for x in segs]}", flush=True)


if __name__ == "__main__":
    failed = main()
    controls_ok = True
    if ("--mutate" not in " ".join(sys.argv) and "--no-controls" not in sys.argv and "--quick" not in sys.argv
            and "--stop" not in " ".join(sys.argv)):
        quick_ctl = "--controls-quick" in sys.argv
        print("\n== Controls: each mutation must make its target checks fail ==", flush=True)
        for name, (want, what) in MUTATIONS.items():
            t1 = time.time()
            pr = subprocess.run([sys.executable, __file__, f"--mutate={name}"] + (["--quick"] if quick_ctl else []),
                                capture_output=True, text=True)
            fails = [ln.split(":")[0][7:] for ln in pr.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = all(w in fails for w in want) and pr.returncode == 1
            controls_ok &= fired
            print(f"[{'FIRED' if fired else 'SILENT'}] --mutate={name} ({what}): exit {pr.returncode}, failed {fails}, "
                  f"{time.time() - t1:.0f} s", flush=True)
            if pr.returncode not in (0, 1):
                print(pr.stderr[-2000:])
        pr = subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "phantom_source.py"),
                             "--mutate=accel-static"], capture_output=True, text=True)
        fails = [ln.split(":")[0][7:] for ln in pr.stdout.splitlines() if ln.startswith("[FAIL]")]
        fired = pr.returncode == 1 and all(w in fails for w in ("gen_xcheck", "accel_xcheck"))
        controls_ok &= fired
        print(f"[{'FIRED' if fired else 'SILENT'}] phantom_source.py --mutate=accel-static: exit {pr.returncode}, failed "
              f"{fails}", flush=True)
        print(f"controls: {'all fired' if controls_ok else 'SOME DID NOT FIRE'}")
    sys.exit(0 if not failed and controls_ok else 1)
