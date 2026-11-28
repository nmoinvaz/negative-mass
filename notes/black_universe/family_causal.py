#!/usr/bin/env python3
"""The causal structure of the black universes' exact time-dependent family (#278).

Member c (b = 1, K = -2c/pi), ingoing null chart (v, rho), phi = sqrt2 arctan rho (#269, exact_mode.py):

    ds^2 = -F dv^2 + 2 dv drho + r^2 dOmega^2,   r = lam(v) R(rho),   F = A - 2 rho lam'/lam + 1/lam^2 - 1,   lam^2 = 1 + s,   s = beta e^{K v},

R^2 = 1 + rho^2, A = 1 + K rho - K R^2 arccot(rho).  v-translation maps beta to beta e^{K d}: only the sign of beta matters; beta = -1 puts v* = 0.
v = const are the ingoing null rays (rho is an affine parameter on them); the outgoing ones solve d rho/dv = F/2.  Every claim about the diagram is
a claim about these two families of curves in the (v, rho) plane.

    family_causal.py derive [--mutate=NAME]               sympy checks S0 - S6
    family_causal.py check  [--only=..] [--mutate=NAME]   numerical checks N0 - N13 (mpmath 40 digits and float64 ray integrations); exit 1 when one fails
    family_causal.py figures [OUTDIR]                     the four conformal diagrams, family_causal.svg
    family_causal.py mutations DIR PREFIX                 fired sets of stored logs DIR/PREFIX<name>.log against the declared ones

Status words as in AGENTS.md.  Ray integrations are numerical evidence.  mpmath for every extended-precision number.
"""
import ast
import os
import re
import sys
import time

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exact_mode as em  # noqa: E402  the closed forms of #269 (A, the exact family, the Kretschmann routine)

sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2", "nonradial"))
try:
    from data_graphs import GRID, INK, INK_SOFT, SERIES, SURFACE, Svg, footnote  # noqa: E402  the repository's vendored SVG helpers (standard library only)
except ImportError:                                                          # the checks do not need them
    Svg = None

MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
ONLY = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--only=")), None)
RESULTS = []
mp.mp.dps = 40

# one perturbation per mutation: the claim or the integrator it changes
AH_COEF = 1.01 if MUT == "ah_coef" else 1.0            # coefficient of the s term of (grad r)^2 and of s_AH
AH_SIGN = -1.0 if MUT == "ah_sign" else 1.0
ELL_SIGN = -1.0 if MUT == "ell_sign" else 1.0          # sign of the outgoing expansion in the region labels
THROAT = -1.0 if MUT == "throat" else 1.0              # F(v, 0) = c + THROAT / lam^2
EH_SHIFT = 5e-3 if MUT == "eh_shift" else 0.0          # the claimed event horizon, displaced
PSTAR = 1.01 if MUT == "pstar" else 1.0                # the claimed end point P* = (v*, -1/(PSTAR K))
BETA_SIGN = -1.0 if MUT == "beta_sign" else 1.0        # the expansion case run with the collapse sign
RATE_RAY = 1.01 if MUT == "rate_ray" else 1.0          # the integrator's rate K (the closed forms and the metric jets keep K)
HZ_SHIFT = 0.01 if MUT == "hz_shift" else 0.0          # the claimed static horizon, displaced
INV_K = 0.99 if MUT == "inv_k" else 1.0                # the v-shift of the |beta| invariance check
ELLIS = -1.0 if MUT == "ellis" else 1.0               # the planted Ellis family with the sign of beta in F reversed
AFF = -1.0 if MUT == "aff_sign" else 1.0              # the affine-parameter integrand e^{+a} with the sign of a reversed
KA = 1.01 if MUT == "kA" else 1.0                     # the coefficient of rho + R^2 arctan(rho) in A, against the rate K
R_SCALE = 1.01 if MUT == "r_scale" else 1.0           # the areal radius of the planted static member
F2 = 2.02 if MUT == "f2" else 2.0                     # the coefficient of rho lam'/lam in the large-rho metric
GPP = 0.01 if MUT == "gpp" else 0.0                   # a g_{rho rho} term in the metric of the affine-parameter check
DS_H = 1.01 if MUT == "ds_h" else 1.0                 # the claimed de Sitter rate H^2 = -2c of the left side of the cone


def want(cid):
    return ONLY is None or cid in ONLY or any(cid.startswith(o) for o in ONLY)


def check(cid, text, ok, detail=""):
    RESULTS.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


# ----------------------------------------------------------------------------------------------------------------------
# derive: S0 - S6 (sympy)
# ----------------------------------------------------------------------------------------------------------------------
v, rho = em.v, em.rho
c, beta = sp.symbols("c beta")
Ksym = -2 * c / sp.pi
R2 = 1 + rho ** 2
Asym = 1 + c * R2 + KA * Ksym * (rho + R2 * sp.atan(rho))     # = 1 + K rho - K R^2 (pi/2 - arctan rho)


def zero(e):
    return sp.simplify(e) == 0


def family_sym():
    s = beta * sp.exp(Ksym * v)
    lam = sp.sqrt(1 + s)
    r = lam * sp.sqrt(R2)
    F = Asym - 2 * rho * sp.diff(lam, v) / lam + 1 / lam ** 2 - 1
    return s, lam, r, F


def derive_S0():
    """Planted: the static member (beta = 0).  (grad r)^2 = A R'^2 = rho^2 A / R^2, zero at the static horizon A = 0 and at the throat."""
    r, F = R_SCALE * sp.sqrt(R2), Asym
    g2 = 2 * sp.diff(r, v) * sp.diff(r, rho) + F * sp.diff(r, rho) ** 2
    ok = zero(g2 - rho ** 2 * Asym / R2)
    check("S0", "static member: (grad r)^2 = g^{ab} d_a r d_b r = rho^2 A / R^2 for symbolic c (planted truth: the apparent horizon of the static member is A = 0 and the throat)", ok)


def derive_S1():
    s, lam, r, F = family_sym()
    g2 = 2 * sp.diff(r, v) * sp.diff(r, rho) + F * sp.diff(r, rho) ** 2          # g^{vv} = 0, g^{v rho} = 1, g^{rho rho} = F
    closed = rho / R2 * (rho * Asym + AH_SIGN * AH_COEF * s * R2 * sp.diff(Asym, rho) / 2)
    check("S1", "family: (grad r)^2 = (rho/R^2) [rho A + s R^2 A'/2], symbolic c and beta, from the metric components", zero(g2 - closed))


def derive_S2():
    ok1 = zero(R2 * sp.diff(Asym, rho) - 2 * Ksym - 2 * rho * (Asym - 1))
    arccot = sp.pi / 2 - sp.atan(rho)
    ok2 = zero(sp.diff(Asym, rho) - 2 * Ksym * (1 - rho * arccot))
    check("S2", "R^2 A' = 2K + 2 rho (A - 1) and A' = 2K (1 - rho arccot rho) (K > 0 for c < 0 and rho arccot rho < 1, so A' > 0 for every rho: proved by arctan x < x)", ok1 and ok2,
          f"identity: {'0' if ok1 else 'nonzero'}, derivative: {'0' if ok2 else 'nonzero'}")


def derive_S3():
    s, lam, r, F = family_sym()
    ell = sp.diff(r, v) + F * sp.diff(r, rho) / 2              # outgoing null normal l = d_v + (F/2) d_rho
    nn = -sp.diff(r, rho)                                      # ingoing null normal n = -d_rho (g(l, n) = -1)
    P = rho * Asym + s * R2 * sp.diff(Asym, rho) / 2
    ok_l = zero(ell - ELL_SIGN * P / (2 * sp.sqrt(R2) * lam))
    ok_n = zero(nn + lam * rho / sp.sqrt(R2))
    g2 = 2 * sp.diff(r, v) * sp.diff(r, rho) + F * sp.diff(r, rho) ** 2
    ok_g = zero(g2 + 2 * ell * nn)
    check("S3", "l(r) = P/(2 R lam), n(r) = -lam rho/R with P = rho A + s R^2 A'/2, and (grad r)^2 = -2 l(r) n(r)", ok_l and ok_n and ok_g,
          f"l: {'0' if ok_l else 'nonzero'}, n: {'0' if ok_n else 'nonzero'}, product: {'0' if ok_g else 'nonzero'}")


def derive_S4():
    s, lam, r, F = family_sym()
    sAH = -2 * rho * Asym / (R2 * sp.diff(Asym, rho))
    P = rho * Asym + s * R2 * sp.diff(Asym, rho) / 2
    ok_a = zero(P - R2 * sp.diff(Asym, rho) / 2 * (s - AH_SIGN * AH_COEF * sAH))
    F0 = sp.simplify(F.subs(rho, 0))
    ok_f = zero(F0 - (c + THROAT / lam ** 2))
    ok_k = zero(F.subs(rho, -1 / Ksym) - Asym.subs(rho, -1 / Ksym))
    check("S4", "P = (R^2 A'/2)(s - s_AH(rho)), s_AH = -2 rho A/(R^2 A'): the apparent horizon is s = s_AH(rho) (and rho = 0); F(v, 0) = c + 1/lam^2; F(v, -1/K) = A(-1/K) for every v", ok_a and ok_f and ok_k,
          f"P: {'0' if ok_a else 'nonzero'}, F(v,0): {'0' if ok_f else 'nonzero'}, F(v,-1/K): {'0' if ok_k else 'nonzero'}")


def derive_S5():
    """Large rho: with A = 1 and R = rho the 2D part is -F dv^2 + 2 dv drho, F = 1/lam^2 - 2 rho lam'/lam; r = lam rho, dT = dv/lam gives -dT^2 + 2 dT dr."""
    lamf = sp.Function("lam")(v)
    rr, T = sp.symbols("r T")
    F = 1 / lamf ** 2 - F2 * rho * sp.diff(lamf, v) / lamf
    dv, dr = sp.symbols("dv dr")
    drho = (dr - rho * sp.diff(lamf, v) * dv) / lamf                  # d rho at r = lam rho
    ds2 = sp.expand(-F * dv ** 2 + 2 * dv * drho)
    mink = sp.expand(-(dv / lamf) ** 2 + 2 * (dv / lamf) * dr)         # -dT^2 + 2 dT dr
    check("S5", "A -> 1 at rho -> +infinity: the (v, rho) part of the family is exactly -dT^2 + 2 dT dr, T = int dv/lam, r = lam rho (Minkowski)", zero(ds2 - mink))


def derive_S6():
    """rho is an affine parameter on v = const: the Christoffels Gamma^a_{rho rho} of -F dv^2 + 2 dv drho vanish for any F(v, rho)."""
    Ff = sp.Function("F")(v, rho)
    g = sp.Matrix([[-Ff, 1], [1, GPP * Ff]])
    gi = g.inv()
    X = [v, rho]
    ok = all(zero(sum(gi[a, d] * (sp.diff(g[d, 1], X[1]) + sp.diff(g[d, 1], X[1]) - sp.diff(g[1, 1], X[d])) for d in range(2)) / 2) for a in range(2))
    check("S6", "Gamma^v_{rho rho} = Gamma^rho_{rho rho} = 0 for any F(v, rho): rho is affine on the ingoing rays v = const", ok)


DERIVE = [("S0", derive_S0), ("S1", derive_S1), ("S2", derive_S2), ("S3", derive_S3), ("S4", derive_S4), ("S5", derive_S5), ("S6", derive_S6)]

# ----------------------------------------------------------------------------------------------------------------------
# float64 machinery for the ray integrations
# ----------------------------------------------------------------------------------------------------------------------


def arccot(x):
    return np.arctan2(1.0, x)


class Fam:
    """Member c, s = beta e^{K v} (beta = 0 static).  rate_scale multiplies K inside the integrator only (the mutation rate_ray)."""

    def __init__(self, c, beta, rate_scale=1.0):
        self.c, self.beta = c, beta
        self.K = -2 * c / np.pi
        self.Kr = self.K * rate_scale
        self.vs = -np.log(-beta) / self.Kr if beta < 0 else np.inf

    def A(self, r):
        return 1 + self.K * r - self.K * (1 + r * r) * arccot(r)

    def Ap(self, r):
        return 2 * self.K * (1 - r * arccot(r))

    def sv(self, v):
        if self.beta == 0:
            return 0.0
        return self.beta * np.exp(self.Kr * v) if self.beta > 0 else -np.exp(self.Kr * (v - self.vs))

    def lam2(self, v):
        if self.beta == 0:
            return 1.0
        return 1 + self.beta * np.exp(self.Kr * v) if self.beta > 0 else -np.expm1(self.Kr * (v - self.vs))

    def sfrac(self, v):                   # s/(1+s)
        if self.beta == 0:
            return 0.0
        if self.beta > 0:
            return 1 / (1 + 1 / self.sv(v))
        x = self.Kr * (v - self.vs)
        return np.exp(x) / np.expm1(x)

    def F(self, v, r):
        return self.A(r) - self.sfrac(v) * (1 + self.K * r)

    def Fr(self, v, r):
        return self.Ap(r) - self.sfrac(v) * self.K

    def rhoh(self):
        return brentq(self.A, 0.0, 10) if self.c < -1 else (0.0 if self.c == -1 else brentq(self.A, -10, 0))

    def ydisc(self, r):
        """Near v* the rays are y = C tau^{-1/2} + O(tau) with y = 1 + K rho; the ray into P* has C = 0."""
        return 1 + self.K * r


class Ellis(Fam):
    """c = 0 (K = 0, A = 1): lam^2 = 1 + beta v, F = (1 - beta rho)/lam^2 (#269 Sec. 3).  Exact rays: 1 - beta rho = C (1 + beta v)^{-1/2}.  Planted truth for the ray code."""

    def __init__(self, beta):
        self.c, self.beta, self.K, self.Kr = 0.0, beta, 0.0, 0.0
        self.vs = -1 / beta if beta < 0 else np.inf

    def lam2(self, v):
        return 1 + self.beta * v

    def sv(self, v):
        return self.beta * v

    def F(self, v, r):
        return (1 - ELLIS * self.beta * r) / self.lam2(v)

    def ydisc(self, r):
        return 1 - self.beta * r


def s_AH(f, r):
    """The closed-form apparent horizon: the value of s at which the outgoing expansion of the sphere at rho vanishes."""
    return -AH_SIGN * AH_COEF * 2 * r * f.A(r) / ((1 + r * r) * f.Ap(r))


def P_claim(f, v, r):
    """Sign-bearing form of the claim: l(r) is proportional to (s - s_AH(rho)) (times R^2 A'/2 > 0 and 1/(2 R lam))."""
    return ELL_SIGN * (f.sv(v) - s_AH(f, r))


WMAX = np.pi / 2 - 1e-9


def ray(f, v0, r0, vend=None, rtol=1e-11, hi=None, lo=None, dense=False):
    """Outgoing null ray d rho/dv = F/2 forward from (v0, r0), integrated as w = arctan(rho).  Stops at |w| >= WMAX, at rho = hi or lo if given, or at vend.

    Returns (kind, v_end, rho_end, sol): kind 'fall' (rho -> -infinity), 'top' (rho -> +infinity before vend), 'hi' / 'lo' (the level reached) or 'end' (vend reached)."""
    if vend is None:
        vend = f.vs - 1e-12 if f.beta < 0 else v0 + 150.0

    def rhs(v, y):
        return [f.F(v, np.tan(y[0])) * np.cos(y[0]) ** 2 / 2]

    evs = [lambda v, y: y[0] - WMAX, lambda v, y: y[0] + WMAX]
    names = ["top", "fall"]
    if hi is not None:
        evs.append(lambda v, y: y[0] - np.arctan(hi))
        names.append("hi")
    if lo is not None:
        evs.append(lambda v, y: y[0] - np.arctan(lo))
        names.append("lo")
    for e in evs:
        e.terminal = True
    sol = solve_ivp(rhs, (v0, vend), [np.arctan(r0)], events=evs, rtol=rtol, atol=1e-14, method="DOP853", dense_output=dense)
    kind = "end"
    for nm, te in zip(names, sol.t_events):
        if len(te):
            kind = nm
    return kind, sol.t[-1], float(np.tan(sol.y[0, -1])), sol


def side(f, v0, r0):
    """+1 when the outgoing ray from (v0, r0) ends on the AF side (rho -> +infinity at v*), -1 when it ends at rho = -infinity (at v < v* or on the cone v = v*).
    Near v* the rays are y = 1 + K rho = C tau^{-1/2} + O(tau), tau = v* - v: the sign of y at the last step is the sign of C (the ray into P* has C = 0)."""
    kind, ve, re_, _ = ray(f, v0, r0)
    return -1 if (kind == "fall" or f.ydisc(re_) < 0) else +1


def eh_backward(f, v0, tau0=1e-7, dense=True):
    """beta < 0: the unique outgoing ray into P* = (v*, -1/(PSTAR K)), from y = 1 + K rho = -(K A0/3) tau at tau = v* - v = tau0, integrated back to v0.  Returns rho_EH(v)."""
    K = f.K
    r_p = -1 / (PSTAR * K)
    A0 = f.A(r_p)
    y0 = -(K * A0 / 3) * tau0

    def rhs(tau, y):
        return [-f.F(f.vs - tau, y[0]) / 2]

    sol = solve_ivp(rhs, (tau0, f.vs - v0), [r_p + y0 / K], rtol=1e-12, atol=1e-14, method="DOP853", dense_output=dense)
    return sol


def rho_eh(sol, f, v):
    return float(sol.sol(f.vs - v)[0]) + EH_SHIFT


def grid(a, b, n=30000):
    """n points from a to b (same sign) spaced geometrically in |x|."""
    xs = np.geomspace(abs(a), abs(b), n)
    return np.sign(a) * xs


def ah_roots(f, v):
    """Roots of s_AH(rho) = s(v): dict with 'plus' (rho > 0 outer), 'minus' (rho < 0), and for s > 0 the two roots below rho_h ('in', 'out')."""
    s = f.sv(v)
    out = {}
    g = lambda x: s_AH(f, x) - s  # noqa: E731

    def first_root(a, b, n=30000):
        xs = grid(a, b, n)
        gs = g(xs)
        idx = np.where(gs[:-1] * gs[1:] < 0)[0]
        return [brentq(g, xs[i], xs[i + 1], xtol=1e-14, rtol=1e-14) for i in idx]

    if s < 0:
        rm = first_root(-1e-9, -1e5)
        rp = first_root(1e-9, 50.0)
        out["minus"] = rm[0] if rm else None
        out["plus"] = rp[-1] if rp else None
    elif f.c < -1:
        rh = f.rhoh()
        rr = first_root(1e-9, rh * (1 - 1e-9))
        out["in"], out["out"] = (rr[0], rr[-1]) if len(rr) == 2 else (None, None)
    return out


# ----------------------------------------------------------------------------------------------------------------------
# check: N0 - N10
# ----------------------------------------------------------------------------------------------------------------------
MEMBERS = [(sp.Integer(-1), -1.0), (sp.Rational(-3, 2), -1.5)]


def N0():
    """Planted truth, beta = 0: the event horizon of the static member is rho_h; the ray code's separatrix finds it."""
    worst, rows = 0.0, []
    for cq in (-1.0, -1.5, -0.3):
        f = Fam(cq, 0.0)
        rh = f.rhoh() + HZ_SHIFT                      # the claim (mutation hz_shift displaces it)
        cls = lambda x: +1 if ray(f, 0.0, x, vend=600.0, hi=f.rhoh() + 1, lo=f.rhoh() - 1)[0] == "hi" else -1  # noqa: E731
        lo, hi = f.rhoh() - 0.2, f.rhoh() + 0.2
        for _ in range(55):
            mid = 0.5 * (lo + hi)
            if cls(mid) > 0:
                hi = mid
            else:
                lo = mid
        sep = 0.5 * (lo + hi)
        err = abs(sep - rh)
        ok_rays = cls(rh + 1e-4) > 0 and cls(rh - 1e-4) < 0
        worst = max(worst, err)
        rows.append(f"c = {cq}: separatrix {sep:.12f}, claimed rho_h {rh:.12f}, diff {err:.1e}, rays at rho_h +-1e-4: {'outside escapes, inside falls' if ok_rays else 'WRONG'}")
        if not ok_rays:
            worst = max(worst, 1.0)
    # the time-dependent planted truth: c = 0, exact rays 1 - beta rho = C (1 + beta v)^{-1/2}; for beta = -1 the ray into P* = (1, -1) is rho = -1 (C = 0)
    worst_e = 0.0
    for bq, v0, r0, v1 in ((1.0, 0.0, 0.3, 5.0), (-1.0, 0.0, 0.3, 0.9), (-1.0, 0.0, -2.0, 0.9)):
        fe = Ellis(bq)
        kind, ve, re_, _ = ray(fe, v0, r0, vend=v1, rtol=1e-12)
        C = (1 - bq * r0) * (1 + bq * v0) ** 0.5
        exact = (1 - C * (1 + bq * v1) ** -0.5) / bq
        worst_e = max(worst_e, abs(re_ - exact) / (1 + abs(exact)))
    fe = Ellis(-1.0)
    lo, hi = -1.2, -0.8
    for _ in range(55):
        mid = 0.5 * (lo + hi)
        if side(fe, 0.0, mid) > 0:
            hi = mid
        else:
            lo = mid
    sep_e = 0.5 * (lo + hi)
    worst_e = max(worst_e, abs(sep_e + 1.0))
    rows.append(f"Ellis c = 0: exact rays to {worst_e:.1e}, beta = -1 separatrix {sep_e:.12f} (exact -1)")
    check("N0", "planted: beta = 0, c = -1, -3/2, -3/10: the ray separatrix (55 bisections) equals the static horizon rho_h to 1e-9 and rays at rho_h +- 1e-4 escape / fall; c = 0 (Ellis, time dependent): rays equal 1 - beta rho = C (1 + beta v)^{-1/2} to 1e-9 and the beta = -1 separatrix is rho = -1",
          worst < 1e-9 and worst_e < 1e-9, "; ".join(rows))


_JETS = {}


def jets_num(cq, bq):
    """numpy-lambdified r_v, r_rho, F of the exact family from exact_mode's explicit expressions (L = 1)."""
    key = (cq, bq)
    if key not in _JETS:
        rexpr, Fexpr = em.family_exprs(cq, "1", bq)[:2]
        _JETS[key] = tuple(sp.lambdify((v, rho), e, "numpy") for e in (sp.diff(rexpr, v), sp.diff(rexpr, rho), Fexpr, rexpr)) + (rexpr, Fexpr)
    return _JETS[key]


def N1():
    """(grad r)^2 = 0: the closed form s = s_AH(rho) against a root of (grad r)^2 computed from the metric components, to 1e-12, and the residual of the closed form."""
    worst_root, worst_res, rows = mp.mpf(0), mp.mpf(0), []
    for cq, cf in MEMBERS:
        for bq in ("-1", "1"):
            rexpr, Fexpr = em.family_exprs(cq, "1", bq)[:2]
            rv = sp.lambdify((v, rho), sp.diff(rexpr, v), "mpmath")
            rr = sp.lambdify((v, rho), sp.diff(rexpr, rho), "mpmath")
            Ff = sp.lambdify((v, rho), Fexpr, "mpmath")
            g2 = lambda vv, x: 2 * rv(vv, x) * rr(vv, x) + Ff(vv, x) * rr(vv, x) ** 2  # noqa: E731
            f = Fam(cf, float(bq))
            Kq = -2 * mp.mpf(cf) / mp.pi
            Am = lambda x: 1 + Kq * x - Kq * (1 + x * x) * mp.acot(x) if x > 0 else 1 + Kq * x - Kq * (1 + x * x) * (mp.pi + mp.atan(1 / x))  # noqa: E731
            Apm = lambda x: 2 * Kq * (1 - x * (mp.acot(x) if x > 0 else mp.pi + mp.atan(1 / x)))  # noqa: E731
            sAHm = lambda x: -AH_SIGN * AH_COEF * 2 * x * Am(x) / ((1 + x * x) * Apm(x))  # noqa: E731
            # the line rho = 0: (grad r)^2 = 0 at every v
            z0 = max(abs(g2(mp.mpf(vv), mp.mpf(0))) for vv in (-3, -0.5, 0.3))
            worst_res = max(worst_res, z0)
            targets = []
            if bq == "-1":
                targets = [("plus", s_t) for s_t in (-0.05, -0.3, -0.7, -0.95)] + [("minus", s_t) for s_t in (-0.05, -0.3, -0.7, -0.95)]
            elif cq == sp.Rational(-3, 2):
                targets = [("in", 0.01), ("out", 0.01), ("in", 0.03), ("out", 0.03)]
            for br, s_t in targets:
                vv = mp.log(abs(mp.mpf(s_t)) / abs(mp.mpf(bq))) / Kq
                guess = {"plus": None, "minus": None, "in": None, "out": None}
                # float root of s_AH(rho) = s_t on the branch, then refined in mpmath
                sgn = {"plus": (1e-9, 50.0), "minus": (-1e-9, -1e5), "in": (1e-9, f.rhoh() * (1 - 1e-9)), "out": (1e-9, f.rhoh() * (1 - 1e-9))}[br]
                xs = grid(sgn[0], sgn[1])
                gs = s_AH(f, xs) - s_t
                idx = np.where(gs[:-1] * gs[1:] < 0)[0]
                if len(idx) == 0:
                    worst_root = mp.inf
                    rows.append(f"c = {cf}, beta = {bq}, {br}, s = {s_t}: the closed form has no root")
                    continue
                x0 = xs[idx[-1] if br in ("plus", "out") else idx[0]]
                rs = mp.findroot(lambda x: sAHm(x) - mp.mpf(s_t), mp.mpf(x0))        # closed form: root of s_AH(rho) = s
                a_, b_ = rs - mp.mpf("0.01") * (abs(rs) + mp.mpf("0.1")), rs + mp.mpf("0.01") * (abs(rs) + mp.mpf("0.1"))
                if g2(vv, a_) * g2(vv, b_) >= 0:
                    worst_root = mp.inf
                    rows.append(f"c = {cf}, beta = {bq}, {br}, s = {s_t}: (grad r)^2 has no sign change within 1% of the closed-form root {mp.nstr(rs, 8)}")
                    continue
                rn = mp.findroot(lambda x: g2(vv, x), (a_, b_), solver="anderson", tol=mp.mpf("1e-34"), maxsteps=200)  # numerical root of (grad r)^2 from the metric
                dd = abs(rn - rs)
                worst_root = max(worst_root, dd)
                rows.append(f"c = {cf}, beta = {bq}, {br}, s = {s_t}: rho = {mp.nstr(rs, 14)}, diff {mp.nstr(dd, 2)}")
    check("N1", "(grad r)^2 = 2 r_v r_rho + F r_rho^2 from the metric: its root at fixed v equals the closed-form root of s_AH(rho) = s to 1e-12 (c = -1, -3/2; beta = -1: both branches at four s; beta = 1, c = -3/2: the two roots), and it vanishes on rho = 0 to 1e-30",
          worst_root < mp.mpf("1e-12") and worst_res < mp.mpf("1e-30"), f"worst root difference {mp.nstr(worst_root, 3)}, rho = 0 residual {mp.nstr(worst_res, 3)}; {len(rows)} roots; " + "; ".join(rows[:2]))


def N2():
    """The four regions from the sign of (n(r), l(r)): metric jets and an integrated ray against the closed-form labels."""
    bad, n, maxrel, maxray = [], 0, 0.0, 0.0
    for cq, cf in MEMBERS:
        for bq in (-1, 1):
            rv, rr, Ff, rf, _, _ = jets_num(cq, str(bq))
            f = Fam(cf, float(bq))
            vgrid = (-6.0, -3.0, -1.0, -0.3, -0.05) if bq < 0 else (-6.0, -3.0, -1.0, 0.0, 1.0, 3.0)
            for vv in vgrid:
                for x in (-8.0, -3.0, -1.0, -0.3, -0.05, 0.05, 0.2, 0.4, 0.7, 1.5, 4.0, 10.0):
                    nval, lval = -rr(vv, x), rv(vv, x) + Ff(vv, x) * rr(vv, x) / 2
                    lam = np.sqrt(f.lam2(vv))
                    Pv = x * f.A(x) + f.sv(vv) * (1 + x * x) * f.Ap(x) / 2
                    lclaim = ELL_SIGN * Pv / (2 * np.sqrt(1 + x * x) * lam)
                    sgn_claim = np.sign(P_claim(f, vv, x))
                    sgn_n = -np.sign(x)
                    n += 1
                    if abs(lval) > 1e-9 and abs(Pv) > 1e-9:
                        maxrel = max(maxrel, abs(lval - lclaim) / max(abs(lval), abs(lclaim)))
                        if np.sign(lval) != sgn_claim or np.sign(nval) != sgn_n:
                            bad.append((cf, bq, vv, x))
                        if abs(lval) > 1e-3:                       # an integrated ray: dr/dv along it equals l(r)
                            fr = Fam(cf, float(bq), RATE_RAY)            # the integrator's family; the labels and the jets are the explicit family
                            kind, ve, re_, sol = ray(fr, vv, x, vend=vv + 1e-6, rtol=1e-13)
                            r1 = np.sqrt(fr.lam2(ve)) * np.sqrt(1 + re_ ** 2)
                            r0 = np.sqrt(fr.lam2(vv)) * np.sqrt(1 + x * x)
                            dr = (r1 - r0) / (ve - vv)
                            maxray = max(maxray, abs(dr - lval) / abs(lval))
                            if np.sign(dr) != sgn_claim:
                                bad.append((cf, bq, vv, x, "ray"))
    check("N2", "region labels: at 4 x (5-6 v) x 12 points the signs of n(r) = -r_rho and l(r) = r_v + F r_rho/2 from the metric jets equal the closed-form labels (-sign rho; sign(s - s_AH(rho))), l(r) equals P/(2R lam) (1e-9), and dr/dv along an integrated ray equals l(r) (1e-3)",
          not bad and maxrel < 1e-9 and maxray < 1e-3, f"{n} points; worst relative l(r) difference {maxrel:.1e}, worst ray {maxray:.1e}; disagreements {bad[:3]}")


def N3():
    """beta > 0: every outgoing ray falls to rho = -infinity at a finite v, with r -> infinity, from starts up to rho_0 = 1e8; no ray escapes to the AF side."""
    bad, rows, vfall = [], [], []
    for cq, cf in MEMBERS:
        f = Fam(cf, BETA_SIGN * 1.0)
        for v0 in (-10.0, -3.0, 0.0, 2.0):
            if f.beta < 0 and v0 >= f.vs - 0.01:
                continue                                             # only the mutation beta_sign: no chart beyond v*
            for r0 in (-30.0, -1.0, 0.0, 0.3, 1.0, 10.0, 1e3, 1e5, 1e8):
                kind, ve, re_, _ = ray(f, v0, r0, vend=v0 + 150.0)
                r_end = np.sqrt(f.lam2(ve)) * np.sqrt(1 + re_ ** 2)
                if kind != "fall" or r_end < 1e6:
                    bad.append((cf, v0, r0, kind))
                if v0 == -10.0:
                    vfall.append((cf, r0, ve))
    # the late-time 2D metric: F -> A - 1 - K rho = -K R^2 arccot(rho) < 0 for every rho
    worst_late, neg = 0.0, True
    for cq, cf in MEMBERS:
        f = Fam(cf, 1.0)
        for x in np.linspace(-20, 20, 81):
            Finf = -f.K * (1 + x * x) * arccot(x)
            neg &= Finf < 0
            worst_late = max(worst_late, abs(f.F(60.0, x) - Finf) / (1 + x * x))
    mono = all(vfall[i][2] < vfall[i + 1][2] for i in range(len(vfall) - 1) if vfall[i][0] == vfall[i + 1][0])
    info = "; ".join(f"c = {cf}: rho_0 = {r0:g} falls at v = {ve:.2f}" for cf, r0, ve in vfall if r0 in (0.0, 1e3, 1e8))
    check("N3", "beta = +1, c = -1, -3/2: 72 outgoing rays from v_0 = -10, -3, 0, 2 and rho_0 from -30 to 1e8 all reach rho = -infinity at finite v with r > 1e6, none escapes, and the fall time v_b rises with rho_0 (the curve rho = -infinity is spacelike); F -> -K R^2 arccot(rho) < 0 at v = 60 (1e-9)",
          not bad and neg and worst_late < 1e-9 and mono, f"{len(bad)} rays escape {bad[:2]}; late-time F difference {worst_late:.1e}; fall time rises with rho_0 (I+ spacelike): {mono}; " + info)


def claim_horizons(f, sol, v0):
    """Candidate surfaces at v_0 (beta < 0): name -> rho."""
    d = {"EH": rho_eh(sol, f, v0), "throat rho = 0": 0.0, "static horizon rho_h": f.rhoh()}
    ah = ah_roots(f, v0)
    d["AH+ (rho > 0)"] = ah["plus"]
    d["AH- (rho < 0)"] = ah["minus"]
    return d


def N4():
    """beta = -1: the event horizon is the separatrix of the outgoing rays; backward from P* and forward by bisection agree; of the candidate surfaces only it separates the two ends."""
    worst_sep, rows, bad = 0.0, [], []
    for cq, cf in MEMBERS:
        f = Fam(cf, -1.0 * 1.0, 1.0 * RATE_RAY)
        sol = eh_backward(f, -10.0)
        for v0 in (-10.0, -3.0):
            r_b = float(sol.sol(f.vs - v0)[0])                                       # the unshifted backward ray
            lo, hi = r_b - 0.2, r_b + 0.2
            if not (side(f, v0, lo) < 0 < side(f, v0, hi)):
                bad.append((cf, v0, "bracket"))
                continue
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if side(f, v0, mid) > 0:
                    hi = mid
                else:
                    lo = mid
            sep = 0.5 * (lo + hi)
            worst_sep = max(worst_sep, abs(sep - r_b))
            rows.append(f"c = {cf}, v_0 = {v0}: backward {r_b:.12f}, bisection {sep:.12f}")
        for v0 in (-10.0, -5.0, -2.0):
            hz = claim_horizons(f, sol, v0)
            res = {}
            for nm, rh in hz.items():
                if rh is None:
                    bad.append((cf, v0, nm, "no root"))
                    continue
                d = 1e-3
                inside = side(f, v0, rh - d)
                outside = side(f, v0, rh + d)
                res[nm] = (inside < 0 < outside)
            expect = {"EH": True, "throat rho = 0": False, "static horizon rho_h": False, "AH+ (rho > 0)": False, "AH- (rho < 0)": False}
            if cf == -1.0:
                expect["static horizon rho_h"] = False
            for nm in res:
                if res[nm] != expect[nm]:
                    bad.append((cf, v0, nm, res[nm]))
            rows.append(f"c = {cf}, v_0 = {v0}: " + ", ".join(f"{nm}: {'separates' if r else 'does not'}" for nm, r in res.items()))
    check("N4", "beta = -1, c = -1, -3/2: the ray into P* (backward from v* - 1e-7) equals the forward-bisection separatrix to 1e-8 at v_0 = -10, -3; rays at EH +- 1e-3 escape / fall at v_0 = -10, -5, -2, while for the throat, the static horizon and both apparent-horizon branches the ray on one side reaches the wrong end (the controls)",
          worst_sep < 1e-8 and not bad, f"worst backward-forward {worst_sep:.1e}; unexpected {bad[:3]}; " + rows[-1])
    return rows


def N5():
    """beta = -1: rho_AH-(v) < rho_EH(v) < rho_AH+(v) on the whole of [-10, v*), the strict inequalities with their smallest gaps."""
    bad, rows = [], []
    for cq, cf in MEMBERS:
        f = Fam(cf, -1.0)
        sol = eh_backward(f, -10.0)
        gap_m, gap_p = np.inf, np.inf
        for tau in np.concatenate([np.geomspace(1e-3, 1, 25), np.linspace(1.2, 10, 30)]):
            vv = f.vs - tau
            ah = ah_roots(f, vv)
            e = rho_eh(sol, f, vv)
            if ah["minus"] is None or ah["plus"] is None:
                bad.append((cf, vv, "no root"))
                continue
            gap_m, gap_p = min(gap_m, e - ah["minus"]), min(gap_p, ah["plus"] - e)
            if not (ah["minus"] < e < ah["plus"]):
                bad.append((cf, vv))
        rows.append(f"c = {cf}: smallest gaps rho_EH - rho_AH- = {gap_m:.4f}, rho_AH+ - rho_EH = {gap_p:.4f}")
    check("N5", "beta = -1, c = -1, -3/2: at 55 values of v from v* - 1e-3 to -10 the apparent horizons bracket the event horizon, rho_AH- < rho_EH < rho_AH+", not bad, "; ".join(rows) + (f"; violations {bad[:3]}" if bad else ""))


def N6():
    """beta = -1: the end point P*: slope, r/sqrt(tau), and the radius r_* at which near rays meet the cone v = v*."""
    rows, ok = [], True
    for cq, cf in MEMBERS:
        f = Fam(cf, -1.0, RATE_RAY)
        K, tau = f.K, 1e-5
        sol = eh_backward(f, -10.0)
        e = float(sol.sol(tau)[0]) + EH_SHIFT
        A0 = f.A(-1 / (PSTAR * K))
        slope = (e + 1 / K) / tau
        r_over = np.sqrt(f.lam2(f.vs - tau)) * np.sqrt(1 + e * e) / np.sqrt(tau)
        r_claim = np.sqrt(K) * np.sqrt(1 + 1 / K ** 2)
        okp = abs(slope - (-A0 / 3)) < 1e-3 * abs(A0 / 3) and abs(r_over - r_claim) < 2e-3 * r_claim
        r_b = float(sol.sol(f.vs - -10.0)[0])
        rs = []
        for d in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
            kind, ve, re_, _ = ray(f, -10.0, r_b + d)
            rs.append(np.sqrt(f.lam2(ve)) * np.sqrt(1 + re_ ** 2) if kind in ("end", "top") else np.nan)
        mono = all(rs[i] > rs[i + 1] > 0 for i in range(len(rs) - 1)) and rs[-1] < 0.1 * rs[1]
        ok &= okp and mono
        rows.append(f"c = {cf}: (rho_EH + 1/K)/tau = {slope:.6f} (-A0/3 = {-A0 / 3:.6f}), r/sqrt(tau) = {r_over:.6f} (sqrt(K) R(-1/K) = {r_claim:.6f}), r_*(d) = " + ", ".join(f"{x:.5f}" for x in rs))
    check("N6", "beta = -1, c = -1, -3/2: along the ray into P* (rho_EH + 1/K)/tau -> -A(-1/K)/3 and r/sqrt(tau) -> sqrt(K) R(-1/K) (1e-3, 2e-3 at tau = 1e-5); rays at rho_EH + d (d = 1e-3 ... 1e-7) reach the cone v = v* at radii r_* that fall strictly with d, r_*(1e-7) < 0.1 r_*(1e-4)",
          ok, "; ".join(rows))


def N7():
    """The throat line rho = 0 and the causal character of the apparent horizons, by rays."""
    bad, rows = [], []
    for cq, cf in MEMBERS:
        for bq in (-1.0, 1.0):
            f = Fam(cf, bq)
            vlist = (-6.0, -2.0, -0.5) if bq < 0 else (-6.0, 0.0, 3.0)
            for vv in vlist:
                kind, ve, re_, _ = ray(f, vv, 0.0, vend=vv + 1e-3, rtol=1e-12)
                lam2 = f.lam2(vv)
                claim = cf + THROAT / lam2                             # F(v, 0)
                moved = re_ / (ve - vv)
                if abs(moved - claim / 2) > 1e-2 * max(abs(claim), 1e-3):
                    bad.append((cf, bq, vv, moved, claim / 2))
            # the throat is null where lam^2 = -1/c (needs lam^2 < 1, beta < 0): the ray from rho = 0 stands still there
            if bq < 0 and cf < -1:
                vn = f.vs + np.log(1 + 1 / cf) / f.K                      # lam^2 = -1/c: s = -1/c - 1 = -1/3 at c = -3/2
                kind, ve, re_, _ = ray(f, vn, 0.0, vend=vn + 1e-6, rtol=1e-13)
                null_dev = abs(re_ / (ve - vn))
                rows.append(f"c = {cf}: throat null at v = {vn:.6f}, d rho/dv there {null_dev:.1e}")
                if null_dev > 1e-6:
                    bad.append((cf, "null", null_dev))
    # the apparent horizons: the ray started on the curve ends inside or outside it as the closed-form slopes say
    for cf, bq, brs, vl in [(-1.0, -1.0, ("plus", "minus"), (-6.0, -2.0, -0.5)), (-1.5, -1.0, ("plus", "minus"), (-6.0, -2.0, -0.5)), (-1.5, 1.0, ("in", "out"), (-6.0, -4.0))]:
        f = Fam(cf, bq)
        for br in brs:
            for vv in vl:
                x = ah_roots(f, vv)[br]
                if x is None:
                    bad.append((cf, br, vv, "no root"))
                    continue
                s0 = f.sv(vv)
                h = 1e-6
                sa = lambda y: s_AH(f, y)  # noqa: E731
                dsd = (sa(x + h) - sa(x - h)) / (2 * h)
                rate_ah = f.Kr * s0 / dsd                              # d rho_AH / dv
                slope_ray = f.F(vv, x) / 2
                dv = 0.02
                kind, ve, re_, _ = ray(f, vv, x, vend=vv + dv, rtol=1e-12)
                x1 = ah_roots(f, vv + dv)[br]
                if x1 is None:
                    bad.append((cf, br, vv, "no root"))
                    continue
                pred = (slope_ray - rate_ah) * dv                       # ray minus curve after dv, first order
                obs = re_ - x1
                char = "spacelike" if rate_ah > slope_ray else "timelike"
                if np.sign(pred) != np.sign(obs) or abs(obs - pred) > 0.1 * abs(pred) + 1e-9:
                    bad.append((cf, br, vv, pred, obs))
                rows.append(f"c = {cf}, AH {br}, v = {vv}: {char} (d rho_AH/dv = {rate_ah:.4f}, ray slope {slope_ray:.4f})")
    check("N7", "rays on the throat rho = 0 move at d rho/dv = F(v, 0)/2 = (c + 1/lam^2)/2 (1%), the throat is null at lam^2 = -1/c (1e-6); the ray from a point of an apparent horizon branch (beta = -1: both branches at v = -6, -2, -0.5; beta = 1, c = -3/2: both roots of the lens at v = -6, -4) ends on the side of the curve that the closed-form slope d rho_AH/dv = K s/s_AH'(rho) against F/2 predicts, to 10%",
          not bad, f"{len(bad)} disagreements {bad[:2]}; " + "; ".join(rows))
    return rows


def N8():
    """The family depends on beta only through a v-translation: the event horizon as a function of s is the same for |beta| = 1 and 0.3."""
    worst = 0.0
    for cq, cf in MEMBERS:
        f1 = Fam(cf, -1.0)
        f2 = Fam(cf, -0.3)
        s1 = eh_backward(f1, -10.0)
        s2 = eh_backward(f2, -10.0)
        for s_t in (-0.01, -0.1, -0.5, -0.9):
            v1 = f1.vs + np.log(-s_t) / f1.K
            v2 = f2.vs + np.log(-s_t) / (f2.K * INV_K)
            worst = max(worst, abs(float(s1.sol(f1.vs - v1)[0]) - float(s2.sol(f2.vs - v2)[0])))
    check("N8", "beta = -1 and -0.3: the ray into P* gives the same rho_EH at the same s (v-translation, 1e-8), c = -1, -3/2, s = -0.01 ... -0.9", worst < 1e-8, f"worst {worst:.1e}")


def N9():
    """Info: the Kretschmann invariant on the event horizon as it approaches P*, against ρ = -1/K and a ray that misses P*."""
    for cq, cf in MEMBERS:
        f = Fam(cf, -1.0)
        sol = eh_backward(f, -10.0)
        out = []
        for u in ("1/100", "1/1000", "1/10000", "1/100000"):
            um = mp.mpf(1) / mp.mpf(u.split("/")[1])
            tau = -mp.log(1 - um) / mp.mpf(f.K)
            e = float(sol.sol(float(tau))[0])
            vals = []
            for r_s in (mp.mpf(e), -1 / mp.mpf(f.K), mp.mpf(e) - mp.mpf("0.3")):
                J, V, Vp, v0, _ = em.eval_family(cq, "1", "-1", u, mp.nstr(r_s, 30))
                vals.append(em.kretschmann(J))
            out.append(f"lam^2 = {u}: K on EH {mp.nstr(vals[0], 8)}, on rho = -1/K {mp.nstr(vals[1], 8)}, at rho_EH - 0.3 {mp.nstr(vals[2], 6)}")
        print(f"  info: c = {cf}: " + "; ".join(out))
        A0 = f.A(-1 / f.K)
        print(f"  info: c = {cf}: the line rho = -1/K = {-1 / f.K:.6f} has F = A = {A0:.6f} at every v, so it is {'spacelike' if A0 < 0 else 'timelike'}; the event horizon ends on it at v*")


def N10():
    """rho = -infinity is infinitely far along the outgoing rays (affine parameter grows like |rho|), and rho is affine on ingoing rays (S6)."""
    ok, rows = True, []
    for cq, cf in MEMBERS:
        for bq in (-1.0, 1.0):
            f = Fam(cf, bq)
            v0 = -4.0
            r0 = -2.0
            # y = (w, a, lam_aff): a = (1/2) int F_rho dv, d lam_aff/dv = e^{a}
            def rhs(vv, y):
                x = np.tan(y[0])
                return [f.F(vv, x) * np.cos(y[0]) ** 2 / 2, f.Fr(vv, x) / 2, np.exp(AFF * y[1])]
            e1 = lambda vv, y: y[0] - np.arctan(-1e3)  # noqa: E731
            e2 = lambda vv, y: y[0] - np.arctan(-1e6)  # noqa: E731
            e2.terminal = True
            e1.terminal = False
            sol = solve_ivp(rhs, (v0, v0 + 20), [np.arctan(r0), 0.0, 0.0], events=[e1, e2], rtol=1e-12, atol=1e-14, method="DOP853")
            if len(sol.y_events[0]) == 0 or len(sol.y_events[1]) == 0:
                ok = False
                rows.append(f"c = {cf}, beta = {bq}: no crossing")
                continue
            la1, la2 = sol.y_events[0][0][2], sol.y_events[1][0][2]
            ratio = la2 / la1
            ok &= abs(ratio / 1e3 - 1) < 0.02
            rows.append(f"c = {cf}, beta = {bq:+.0f}: lam_aff(rho = -1e6)/lam_aff(rho = -1e3) = {ratio:.2f}")
    check("N10", "outgoing rays into rho = -infinity (c = -1, -3/2, beta = +-1, from (v, rho) = (-4, -2)): the affine parameter lam_aff = int exp(int F_rho dv/2) dv grows linearly in |rho|, lam_aff(-1e6)/lam_aff(-1e3) = 1e3 within 2%", ok, "; ".join(rows))


def N11():
    """beta = -1: the rays left of the event horizon.  Close to it they reach the cone v = v* at rho = -infinity with a finite r_*; beyond a critical offset they fall to
    rho = -infinity at v_b < v*; v_b -> v* and r_* -> infinity at the critical offset, so I+ (rho = -infinity) meets the cone at a point Q with r = infinity."""
    ok, rows = True, []
    for cq, cf in MEMBERS:
        f = Fam(cf, -1.0, RATE_RAY)
        sol = eh_backward(f, -10.0)
        v0 = -10.0
        rb = float(sol.sol(f.vs - v0)[0]) + EH_SHIFT
        isfall = lambda d: ray(f, v0, rb - d)[0] == "fall"  # noqa: E731
        lo, hi = 1e-7, 1e-1
        if isfall(lo) or not isfall(hi):
            ok = False
            rows.append(f"c = {cf}: bracket failed")
            continue
        for _ in range(80):
            mid = np.sqrt(lo * hi) if hi / lo > 1.5 else 0.5 * (lo + hi)
            if isfall(mid):
                hi = mid
            else:
                lo = mid
        dc = 0.5 * (lo + hi)
        rs, tf = [], []
        for eps in (1e-1, 1e-2, 1e-3):
            k1, v1, r1, _ = ray(f, v0, rb - dc * (1 - eps))
            rs.append(np.sqrt(f.lam2(v1)) * np.sqrt(1 + r1 ** 2) if (k1 == "end" and 1 + f.K * r1 < 0) else np.nan)
            k2, v2, r2, _ = ray(f, v0, rb - dc * (1 + eps))
            tf.append(f.vs - v2 if k2 == "fall" else np.nan)
        inc = rs[0] < rs[1] < rs[2] and rs[2] > 100
        dec = tf[0] > tf[1] > tf[2] > 0 and tf[2] < 1e-6
        small = [ray(f, v0, rb - dc * q)[0] == "end" for q in (1e-1, 1e-3)]
        tau_far = []
        for dd in (3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 2 * dc):
            k3, v3, r3, _ = ray(f, v0, rb - dd)
            tau_far.append(f.vs - v3 if k3 == "fall" else np.nan)
        mono = all(tau_far[i] > tau_far[i + 1] > 0 for i in range(len(tau_far) - 1))
        ok &= inc and dec and all(small) and mono
        rows.append(f"c = {cf}: d_crit = {dc:.6e}; r_*(eps = 0.1, 0.01, 0.001) = {rs[0]:.4g}, {rs[1]:.4g}, {rs[2]:.4g}; tau_fall = {tf[0]:.3e}, {tf[1]:.3e}, {tf[2]:.3e}; tau_fall at d = 0.3 ... 2 d_crit falls monotonically (I+ spacelike): {mono}")
    check("N11", "beta = -1, c = -1, -3/2: rays at rho_EH(-10) - d reach the cone v = v* at rho = -infinity for d < d_crit and fall at v_b < v* for d > d_crit; at d_crit (1 -+ eps), eps = 0.1, 0.01, 0.001, r_* rises (> 100 at 0.001) and tau_fall = v* - v_b falls (< 1e-6): I+ meets the cone at a point with r = infinity; tau_fall falls monotonically along the 8 offsets from 0.3 to 2 d_crit (the curve I+ is spacelike)",
          ok, "; ".join(rows))


def N12():
    """beta = +1: the trapped region.  c = -1: none at any v.  c = -3/2: a lens between the two roots of s_AH(rho) = s that closes at s_max = max s_AH on (0, rho_h)."""
    ok, rows = True, []
    for cq, cf in MEMBERS:
        rv, rr, Ff, rf, _, _ = jets_num(cq, "1")
        f = Fam(cf, 1.0)
        xs = np.concatenate([grid(1e-7, 1e-1, 2000), np.linspace(0.1, 3.0, 20001), grid(3.0, 1e3, 2000)])

        def regions(vv, xx):
            nval, lval = -rr(vv, xx), rv(vv, xx) + Ff(vv, xx) * rr(vv, xx) / 2
            return nval, lval

        if cf == -1.0:
            worst_l = np.inf
            for vv in np.arange(-12.0, 12.01, 0.5):
                nval, lval = regions(vv, xs)
                P = xs * f.A(xs) + f.sv(vv) * (1 + xs ** 2) * f.Ap(xs) / 2
                worst_l = min(worst_l, np.min(lval[P > 1e-12] / (P[P > 1e-12] / (2 * np.sqrt(1 + xs[P > 1e-12] ** 2) * np.sqrt(f.lam2(vv))))))
                if np.any(lval <= 0):
                    ok = False
            rows.append(f"c = {cf}: l(r) > 0 at all 49 v from -12 to 12 and {len(xs)} rho > 0 (no trapped region)")
        else:
            res = minimize_scalar(lambda x: -s_AH(f, x), bounds=(1e-6, f.rhoh() * (1 - 1e-6)), method="bounded", options={"xatol": 1e-12})
            smax, rm = -res.fun, res.x
            for fac, expect in ((0.995, True), (1.005, False)):
                vv = np.log(fac * smax) / f.K
                nval, lval = regions(vv, xs)
                found = bool(np.any((nval < 0) & (lval < 0)))
                ok &= (found == expect)
                rows.append(f"c = {cf}: at s = {fac} s_max (s_max = {smax:.8f} at rho = {rm:.6f}, v_max = {np.log(smax) / f.K:.5f}) trapped points {'found' if found else 'none'}")
        xn = -grid(1e-7, 1e4, 4000)
        for vv in np.arange(-12.0, 12.01, 1.0):
            nval, lval = regions(vv, xn)
            if not (np.all(nval > 0) and np.all(lval > 0)):
                ok = False
    check("N12", "beta = +1: from the metric jets, c = -1 has l(r) > 0 for every rho > 0 at 49 values of v (no trapped region ever); c = -3/2 has trapped points at s = 0.995 s_max and none at 1.005 s_max (s_max = max of s_AH on (0, rho_h)); rho < 0 is anti-trapped (n, l > 0) at all v",
          ok, "; ".join(rows))



def N13():
    """beta = -1: the curvature along the rays that reach the cone v = v*.  On the AF side (rho -> +infinity) Ricci scalar ~ sqrt(tau) and Kretschmann ~ tau tend to 0
    (flat); on the other side (rho -> -infinity) they tend to the de Sitter values 4V(-infinity) = -24c = 12 H^2 and 24 H^4, H^2 = -2c.  Only the vertex P* is singular."""
    ok, rows = True, []
    for cq, cf in MEMBERS:
        f = Fam(cf, -1.0)
        sol = eh_backward(f, -10.0)
        rb = float(sol.sol(f.vs + 10.0)[0])                         # the unshifted ray into P*
        H2 = -2 * cf * DS_H
        res = {}
        for name, dd in (("right", 1e-4), ("left", -3e-5)):
            kind, ve, re_, sl = ray(f, -10.0, rb + dd, dense=True)
            vals = {}
            for tau in (1e-6, 1e-8):
                rho_t = float(np.tan(sl.sol(f.vs - tau)[0]))
                u = mp.nstr(-mp.expm1(-mp.mpf(f.K) * tau), 30)
                J, V, Vp, v0, _ = em.eval_family(cq, "1", "-1", u, mp.nstr(mp.mpf(rho_t), 25))
                ricci = -J["F"] * 2 / (1 + mp.mpf(rho_t) ** 2) ** 2 + 4 * V      # trace of the field equations: -F phi'^2 + 4 V
                vals[tau] = (float(ricci), float(em.kretschmann(J)))
            res[name] = vals
        rr = res["right"]
        ratio_R = rr[1e-6][0] / rr[1e-8][0]
        ratio_K = rr[1e-6][1] / rr[1e-8][1]
        ok_r = 8 < ratio_R < 12 and 80 < ratio_K < 120 and abs(rr[1e-8][0]) < 0.02 and rr[1e-8][1] < 1e-3
        ll = res["left"][1e-8]
        ok_l = abs(ll[0] / (12 * H2) - 1) < 5e-3 and abs(ll[1] / (24 * H2 ** 2) - 1) < 5e-3
        ok &= ok_r and ok_l
        rows.append(f"c = {cf}: right ray (tau = 1e-6 -> 1e-8): Ricci {rr[1e-6][0]:.4g} -> {rr[1e-8][0]:.4g} (ratio {ratio_R:.2f}), Kretschmann {rr[1e-6][1]:.4g} -> {rr[1e-8][1]:.4g} (ratio {ratio_K:.1f}); "
                    f"left ray at tau = 1e-8: Ricci {ll[0]:.5g} (12 H^2 = {12 * H2:.5g}), Kretschmann {ll[1]:.5g} (24 H^4 = {24 * H2 ** 2:.5g})")
    check("N13", "beta = -1, c = -1, -3/2, rays at rho_EH(-10) + 1e-4 (AF side) and - 3e-5 (other side) to the cone: AF side, Ricci(1e-6)/Ricci(1e-8) in (8, 12), Kretschmann ratio in (80, 120), |Ricci| < 0.02 and Kretschmann < 1e-3 at tau = 1e-8; "
          "other side, Ricci and Kretschmann at tau = 1e-8 within 0.5% of 12 H^2 and 24 H^4, H^2 = -2c", ok, "; ".join(rows))


# ----------------------------------------------------------------------------------------------------------------------
# figures: the four conformal diagrams
# ----------------------------------------------------------------------------------------------------------------------
V0 = -7.0                    # the ingoing null cone v = V0 is the initial surface of every diagram
XMIN = 1e-9                  # u = sign(x) log10(1 + |x|/XMIN) / LNORM, x = rho_0 - rho_ref: any monotone relabeling of the outgoing rays is a conformal diagram
XMAX = 10.0                  # |x| >= XMAX is drawn on the border (rho_0 beyond rho_ref +- 10)
LNORM = np.log10(1 + XMAX / XMIN)


def vmap(f, v):
    """The second null coordinate in (0, 1): linear in v up to v* for beta < 0, (2/pi) arctan((v - V0)/4) for beta > 0."""
    if f.beta < 0:
        return (v - V0) / (f.vs - V0)
    return 2 / np.pi * np.arctan((v - V0) / 4.0)


def ulabel(x):
    return np.sign(x) * np.minimum(1.0, np.log10(1 + np.abs(x) / XMIN) / LNORM)


def xlabel(u):
    return np.sign(u) * XMIN * (10 ** (np.abs(u) * LNORM) - 1)


def rho_ref(f, sol):
    return float(sol.sol(f.vs - V0)[0]) if f.beta < 0 else f.rhoh()


def u_of(f, v, r, ref):
    """The outgoing-ray label of the point (v, rho): integrate the ray back to V0."""
    if abs(v - V0) < 1e-12:
        r0 = r
    else:
        def rhs(vv, y):
            return [f.F(vv, np.tan(y[0])) * np.cos(y[0]) ** 2 / 2]
        sol = solve_ivp(rhs, (v, V0), [np.arctan(r)], rtol=1e-11, atol=1e-14, method="DOP853")
        r0 = float(np.tan(sol.y[0, -1]))
    return float(ulabel(r0 - ref))


def curve(f, ref, vs_, rhos):
    """(u, V) points of a curve rho(v) sampled at vs_."""
    return [(u_of(f, vv, rr, ref), float(vmap(f, vv))) for vv, rr in zip(vs_, rhos) if rr is not None and np.isfinite(rr)]


def diagram_data(cq, cf, bq):
    f = Fam(cf, float(bq))
    sol = eh_backward(f, V0 - 1.0) if bq < 0 else None
    ref = rho_ref(f, sol)
    d = {"f": f, "ref": ref, "sol": sol}
    top = f.vs - 1e-3 if bq < 0 else 40.0
    if bq < 0:
        vs_ = np.concatenate([np.linspace(V0, f.vs - 1.0, 60), f.vs - np.geomspace(1.0, 2e-3, 70)[1:]])
    else:
        vs_ = np.concatenate([np.linspace(V0, 0.0, 40), np.linspace(0.0, 12.0, 60)[1:], np.linspace(12.0, 40.0, 30)[1:]])
    d["vs"] = vs_
    d["throat"] = curve(f, ref, vs_, [0.0] * len(vs_))
    if bq < 0:
        roots = [ah_roots(f, vv) for vv in vs_]
        d["AHp"] = curve(f, ref, vs_, [rt["plus"] for rt in roots])
        d["AHm"] = curve(f, ref, vs_, [rt["minus"] for rt in roots])
        d["EH"] = [(0.0, float(vmap(f, vv))) for vv in vs_]
    elif cf < -1:
        res = minimize_scalar(lambda x: -s_AH(f, x), bounds=(1e-6, f.rhoh() * (1 - 1e-6)), method="bounded", options={"xatol": 1e-12})
        smax = -res.fun
        vmax = np.log(smax) / f.K
        vv_ = np.concatenate([np.linspace(V0, vmax - 0.5, 40), vmax - np.geomspace(0.5, 1e-4, 40)[1:]])
        roots = [ah_roots(f, vv) for vv in vv_]
        d["vs_ah"] = vv_
        d["AHin"] = curve(f, ref, vv_, [rt.get("in") for rt in roots])
        d["AHout"] = curve(f, ref, vv_, [rt.get("out") for rt in roots])
        d["vmax"] = vmax
    # I+: rho = -infinity, from the rays that fall
    us = np.linspace(-1, 1, 500) if bq > 0 else np.linspace(-1, 0, 500)
    iplus = []
    for u in us:
        x = float(xlabel(u))
        r0 = ref + x
        kind, ve, re_, _ = ray(f, V0, r0, vend=(None if bq < 0 else V0 + 80.0))
        if kind == "fall":
            iplus.append((float(u), float(vmap(f, ve))))
    d["Iplus"] = iplus
    return d


def draw_panel(svg, x0, y0, W, d, title, sub, bq, cf):
    f = d["f"]
    P = lambda u, V: (x0 + (V + (u + 1) / 2) / 2 * W, y0 + W / 2 - (V - (u + 1) / 2) / 2 * W)  # noqa: E731
    svg.text(x0, y0 - 24, title, size=12, fill=INK, weight="bold")
    svg.text(x0, y0 - 9, sub, size=10)
    # the diamond
    corners = [P(-1, 0), P(1, 0), P(1, 1), P(-1, 1)]
    svg.add('<polygon points="%s" fill="%s" stroke="%s" stroke-width="1"/>' % (" ".join(f"{a:.1f},{b:.1f}" for a, b in corners), SURFACE, INK_SOFT))
    # grid of null lines
    for q in (-0.5, 0.0, 0.5):
        (a1, b1), (a2, b2) = P(q, 0), P(q, 1)
        svg.add(f'<line x1="{a1:.1f}" y1="{b1:.1f}" x2="{a2:.1f}" y2="{b2:.1f}" stroke="{GRID}" stroke-width="1"/>')
    for q in (0.25, 0.5, 0.75):
        (a1, b1), (a2, b2) = P(-1, q), P(1, q)
        svg.add(f'<line x1="{a1:.1f}" y1="{b1:.1f}" x2="{a2:.1f}" y2="{b2:.1f}" stroke="{GRID}" stroke-width="1"/>')

    def ribbon(left, right, fill, op):
        """polygon between two curves u_left(V), u_right(V) (lists of (u, V)); None edge = the diagram border."""
        if not left or not right:
            return
        pts = [P(u, V) for u, V in left] + [P(u, V) for u, V in reversed(right)]
        svg.add('<polygon points="%s" fill="%s" fill-opacity="%s"/>' % (" ".join(f"{a:.1f},{b:.1f}" for a, b in pts), fill, op))

    Vs = [V for u, V in d["throat"]]
    border_r = [(1.0, V) for V in Vs]
    # the I+ curve as u_I(V) on the V grid (monotone: interpolate)
    ip = d["Iplus"]
    if ip:
        iu = np.array([p[0] for p in ip])
        iV = np.array([p[1] for p in ip])
        order = np.argsort(iV)
        left_edge = [(float(np.interp(V, iV[order], iu[order], left=-1.0)), V) for V in Vs if V <= iV.max()]
    else:
        left_edge = [(-1.0, V) for V in Vs]
    Vcut = [V for u, V in left_edge]
    C_ANTI, C_TRAP, C_INNER = SERIES[1], SERIES[0], SERIES[9]
    tr = d["throat"][:len(left_edge)]
    if bq < 0:
        ahm = d["AHm"][:len(left_edge)]
        ahp = d["AHp"][:len(left_edge)]
        ribbon(left_edge, ahm, C_ANTI, "0.30")
        ribbon(ahm, tr, C_INNER, "0.22")
        ribbon(tr, ahp, C_TRAP, "0.30")
    else:
        ribbon(left_edge, tr, C_ANTI, "0.30")
        if cf < -1 and d.get("AHin"):
            nin = min(len(d["AHin"]), len(d["AHout"]))
            lens_top = d["vmax"]
            ribbon(d["AHin"][:nin], d["AHout"][:nin], C_TRAP, "0.30")
    # curves
    def poly(pts, color, width=1.6, dash=None):
        if len(pts) < 2:
            return
        dd = f' stroke-dasharray="{dash}"' if dash else ""
        svg.add('<polyline points="%s" fill="none" stroke="%s" stroke-width="%s"%s/>' % (" ".join("%.1f,%.1f" % P(u, V) for u, V in pts), color, width, dd))
    poly(d["throat"], INK, 1.2, "5 3")
    if bq < 0:
        poly(d["AHp"], SERIES[1], 1.8)
        poly(d["AHm"], SERIES[3], 1.8)
        poly(d["EH"], SERIES[0], 2.4)
    elif cf < -1 and d.get("AHin"):
        poly(d["AHin"], SERIES[1], 1.8)
        poly(d["AHout"], SERIES[1], 1.8)
    poly(ip, SERIES[7], 2.4)

    def tag(pts, Vt, text, dx, dy, anchor="start", color=INK):
        if not pts:
            return
        u, V = min(pts, key=lambda q: abs(q[1] - Vt))
        a, b = P(u, V)
        svg.text(a + dx, b + dy, text, size=10, fill=color, anchor=anchor)
    if bq < 0:
        tag(d["EH"], 0.80, "event horizon", -4, 14, "start", SERIES[0])
        tag(d["AHp"], 0.28, "AH⁺", 8, 14, "start", SERIES[1])
        tag(d["AHm"], 0.30, "AH⁻", 8, 12, "start", SERIES[3])
        tag(ip, 0.62, "I⁺", -8, -4, "end", SERIES[7])
        tag(d["throat"], 0.12, "throat ρ = 0", 10, 4, "start", INK)
        # curves of finite rho reach P* only as v -> v* (the label tends to the event horizon like sqrt(v* - v)): dotted to P*
        for pts, color in ((d["throat"], INK), (d["AHp"], SERIES[1])):
            if pts:
                poly([pts[-1], (0.0, 1.0)], color, 1.2, "1 3")
    else:
        tag(ip, 0.5, "I⁺", -8, -4, "end", SERIES[7])
        tag(d["throat"], 0.15, "throat ρ = 0", 10, 4, "start", INK)
        if cf < -1 and d.get("AHin"):
            tag(d["AHin"], 0.04, "trapped lens", 4, -10, "end", SERIES[0])
    # the cone and the vertex
    if bq < 0:
        uq = ip[-1][0] if ip else -1.0
        a1, b1 = P(uq, 1.0)
        a2, b2 = P(1.0, 1.0)
        svg.add(f'<line x1="{a1:.1f}" y1="{b1:.1f}" x2="{a2:.1f}" y2="{b2:.1f}" stroke="{SERIES[7]}" stroke-width="2.4" stroke-dasharray="1 0"/>')
        a, b = P(0.5, 1.0)
        svg.text(a + 8, b - 4, "flat-like", size=10, fill=INK)
        a, b = P(uq / 2, 1.0)
        svg.text(a + 8, b - 4, "dS-like", size=10, fill=INK)
        px, py = P(0.0, 1.0)
        svg.add(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="4" fill="{INK}"/>')
        svg.text(px + 6, py - 6, "P* (r = 0)", size=10, fill=INK)
        qx, qy = P(uq, 1.0)
        svg.add(f'<circle cx="{qx:.1f}" cy="{qy:.1f}" r="3.5" fill="none" stroke="{SERIES[7]}" stroke-width="1.5"/>')
        svg.text(qx - 4, qy - 7, "Q (r = ∞)", size=10, fill=INK, anchor="end")
    (a, b) = P(1.0, 0.5)
    svg.text(a + 6, b + 2, "I⁻", size=11, fill=INK)
    (a, b) = P(-1.0, 0.02)
    svg.text(a - 3, b + 14, "v = v₀", size=10, fill=INK_SOFT, anchor="end")
    return P


def figures(outdir):
    """The four conformal diagrams of the family in one SVG, family_causal.svg, and the numbers behind their curves."""
    W, width = 330, 820
    svg = Svg(width)
    svg.text(16, 24, "Conformal diagrams of the exact family: the future of the ingoing null cone v = v₀ (v₀ = −7)", size=14, fill=INK, weight="bold")
    svg.text(16, 42, "Null coordinates: the second is v; the first labels the outgoing rays by where they cross v = v₀ (log scale about the reference ray: the event horizon for β < 0, the static horizon for β > 0).", size=10)
    svg.text(16, 56, "Left: collapse, β < 0 (curves of finite ρ are drawn to v* − 2·10⁻³; they reach P* only in the limit). Right: expansion, β > 0. β = ±1; any other |β| is a v-translation.", size=10)
    order = [(sp.Integer(-1), -1.0, -1), (sp.Integer(-1), -1.0, +1), (sp.Rational(-3, 2), -1.5, -1), (sp.Rational(-3, 2), -1.5, +1)]
    for k, (cq, cf, bq) in enumerate(order):
        col, row = k % 2, k // 2
        d = diagram_data(cq, cf, bq)
        x0 = 56 + col * 400
        y0 = 120 + row * (W + 70)
        if bq < 0:
            sub = "β = −1, v* = 0: the cone v = v* ends the diagram at the top right"
        else:
            sub = "β = +1: no event horizon" + (f"; the trapped lens closes at v = {d['vmax']:.2f}" if cf < -1 else "; no trapped region")
        draw_panel(svg, x0, y0, W, d, f"c = {cf:g}, {'β < 0 (collapse)' if bq < 0 else 'β > 0 (expansion)'}", sub, bq, cf)
    ybot = 120 + 2 * (W + 70) - 20
    items = [("event horizon", SERIES[0]), ("AH⁺ (θ_ℓ = 0, ρ > 0)", SERIES[1]), ("AH⁻ (θ_ℓ = 0, ρ < 0)", SERIES[3]), ("I⁺ (ρ → −∞)", SERIES[7]), ("throat ρ = 0 (θ_n = 0)", INK)]
    lx = 16
    for label, color in items:
        svg.add(f'<line x1="{lx}" y1="{ybot - 4}" x2="{lx + 18}" y2="{ybot - 4}" stroke="{color}" stroke-width="2.4"/>')
        svg.text(lx + 22, ybot, label, size=10, fill=INK)
        lx += 30 + 6.2 * len(label)
    ybot += 18
    lx = 16
    for label, color in [("trapped (θ_ℓ, θ_n < 0)", SERIES[0]), ("anti-trapped (θ_ℓ, θ_n > 0)", SERIES[1]), ("θ_ℓ < 0 < θ_n", SERIES[9])]:
        svg.add(f'<rect x="{lx}" y="{ybot - 9}" width="14" height="10" fill="{color}" fill-opacity="0.35"/>')
        svg.text(lx + 20, ybot, label, size=10, fill=INK)
        lx += 36 + 6.2 * len(label)
    h = footnote(svg, ["family_causal.py figures (#278)", "outgoing rays: d rho/dv = F/2, DOP853, rtol 1e-11", "null-ray evidence, see family_causal.md"], "", [], ybot + 6)
    os.makedirs(outdir, exist_ok=True)
    fn = os.path.join(outdir, "family_causal.svg")
    open(fn, "w").write(svg.finish(h))
    print("wrote", fn)

NUMERIC = [("N0", N0), ("N1", N1), ("N2", N2), ("N3", N3), ("N4", N4), ("N5", N5), ("N6", N6), ("N7", N7), ("N8", N8), ("N9", N9), ("N10", N10), ("N11", N11), ("N12", N12), ("N13", N13)]


# name: checks that must fire (development runs of negative-mass-raw I278-mac/dev, before they were declared)
MUTATIONS = {
    "ah_coef": ["S1", "S4", "N1", "N12"],
    "ah_sign": ["S1", "S4", "N1", "N2", "N4", "N5", "N7", "N12"],
    "ell_sign": ["S3", "N2"],
    "throat": ["S4", "N7"],
    "eh_shift": ["N4", "N6", "N11"],
    "pstar": ["N4", "N6"],
    "beta_sign": ["N3"],
    "rate_ray": ["N2", "N6"],
    "hz_shift": ["N0"],
    "inv_k": ["N8"],
    "ellis": ["N0"],
    "aff_sign": ["N10"],
    "kA": ["S1", "S2", "S3"],
    "r_scale": ["S0"],
    "f2": ["S5"],
    "gpp": ["S6"],
    "ds_h": ["N13"],
}


def run_all(items):
    t0 = time.time()
    for cid, fn in items:
        if want(cid):
            t1 = time.time()
            fn()
            print(f"    ({cid}: {time.time() - t1:.0f} s)", flush=True)
    bad = [cid for cid, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS)} checks, {len(RESULTS) - len(bad)} pass; failed: {bad}; mutation {MUT}; {time.time() - t0:.0f} s")
    if MUT is not None and MUT in MUTATIONS:
        want_ = MUTATIONS[MUT]
        print(f"mutation {MUT}: fired {bad}, declared {want_}: {'equal' if sorted(bad) == sorted(want_) else 'DIFFERENT'}")
    return bad


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
    elif sub == "figures":
        args = [a for a in sys.argv[2:] if not a.startswith("--")]
        figures(args[0] if args else HERE)
        sys.exit(0)
    elif sub == "mutations":
        args = [a for a in sys.argv[2:] if not a.startswith("--")]
        sys.exit(0 if report_mutations(args[0], args[1] if len(args) > 1 else "") else 1)
    else:
        sys.exit(f"unknown command {sub}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
