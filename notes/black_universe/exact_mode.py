#!/usr/bin/env python3
"""Why Y = R e^{K(t+x)} solves the black universes' linear radial problem for every member: an exact nonlinear family (#269).

Member c of the family (b = 1, K = -2c/pi = 3m): static metric  ds^2 = -A dv^2 + 2 dv drho + R^2 dOmega^2 (ingoing null chart, v = t + x, dx = drho/A),
R^2 = 1 + rho^2, phi = sqrt2 arctan rho, A = 1 + c R^2 + K (rho + R^2 arctan rho), the BF06 potential, R_ab = eps phi_a phi_b + V g_ab, eps = -1, eps box phi = V'.

The candidate (class (b) of the issue, "b -> b(v)"): in the chart with phi = sqrt2 arctan rho as the radial coordinate and g_{v rho} = 1, g_{rho rho} = 0,

    r = lam(v) R(rho),     F = A - 2 rho lam'/lam + (1/lam^2 - 1),     ds^2 = -F dv^2 + 2 dv drho + r^2 dOmega^2,     phi unchanged.

Result. Four of the five field equations hold for every lam(v) (S1); the fifth is  (lam lam')' = K lam lam'  (S2), so  lam^2 = L + beta e^{K v}  (S3):
an exact two-parameter family for every member, with the rate K the first integral K = (R^2 A' - 2 rho A)/2 + rho of the static equations (S4).
Linearized about beta = 0 it is delta r = (beta/2) R e^{K(t+x)}: #212's Y_K (S5, N2).  beta < 0 collapses at finite v (N4); beta > 0 expands (N5).

    exact_mode.py derive [--mutate=NAME]               sympy checks S0 - S5
    exact_mode.py check  [--only=..] [--mutate=NAME]   numerical checks N0 - N6 at 50 digits; exit 1 when one fails
    exact_mode.py mutations DIR PREFIX                 fired sets of stored logs DIR/PREFIX<name>.log against the declared ones

Status words as in AGENTS.md. mpmath for every extended-precision number.
"""
import ast
import os
import re
import sys
import time

import mpmath as mp
import sympy as sp

MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
ONLY = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--only=")), None)
RESULTS = []
mp.mp.dps = 50

# the candidate, one coefficient or one ingredient per mutation
CAND = {"kf": sp.Integer(2), "kh": sp.Integer(1), "rate": "K", "shape": "R", "pot": 1, "eps": -1, "charge": 1}
if MUT == "kf":
    CAND["kf"] = sp.Rational(101, 50)            # coefficient 2 of rho lam'/lam
elif MUT == "kh":
    CAND["kh"] = sp.Rational(101, 100)           # coefficient 1 of (1/lam^2 - 1)
elif MUT == "rate":
    CAND["rate"] = "K(1+1e-3)"
elif MUT == "kappa":
    CAND["rate"] = "kappa"                       # numeric checks only: the surface gravity in place of K
elif MUT == "shape":
    CAND["shape"] = "rho"                        # r = lam rho in place of lam R
elif MUT == "naive":
    CAND["kf"], CAND["kh"] = sp.Integer(0), sp.Integer(0)       # b -> b(v) with F = A unchanged
elif MUT == "pot":
    CAND["pot"] = -1                             # V -> -V
elif MUT == "eps":
    CAND["eps"] = 1                              # canonical scalar in place of phantom
elif MUT == "charge":
    CAND["charge"] = 0                           # the charge without its + rho
HAND_DROP = MUT == "hand_drop"                   # drops a Christoffel term of the hand-reduced equations
KRET_DROP = MUT == "kret_drop"                   # the same term in the Kretschmann routine


def want(cid):
    return ONLY is None or cid in ONLY or any(cid.startswith(o) for o in ONLY)


def check(cid, text, ok, detail=""):
    RESULTS.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


# ----------------------------------------------------------------------------------------------------------------------
# the background and the brute-force field equations (sympy)
# ----------------------------------------------------------------------------------------------------------------------
v, rho, th, ph = sp.symbols("v rho theta varphi", real=True)
c, Ks = sp.symbols("c Ks")
K = -2 * c / sp.pi
X4 = [v, rho, th, ph]
Rs = sp.sqrt(1 + rho ** 2)                       # R
Phi = sp.sqrt(2) * sp.atan(rho)                  # phi(rho)


def background(Kc):
    """A, V(phi(rho)), dV/dphi(phi(rho)) for the coefficients (c, Kc); the BF06 potential, linear in (c, Kc)."""
    A = 1 + c * (rho ** 2 + 1) + Kc * (rho + (rho ** 2 + 1) * sp.atan(rho))
    psi = sp.Symbol("psi")
    Vf = -c * (3 - 2 * sp.cos(psi) ** 2) - Kc * (3 * sp.sin(psi) * sp.cos(psi) + psi * (3 - 2 * sp.cos(psi) ** 2))
    V0 = CAND["pot"] * Vf.subs(psi, sp.atan(rho))
    V1 = CAND["pot"] * (sp.diff(Vf, psi) / sp.sqrt(2)).subs(psi, sp.atan(rho))
    return A, V0, V1


def null_metric(F, r):
    g = sp.zeros(4, 4)
    g[0, 0], g[0, 1], g[1, 0] = -F, 1, 1
    g[2, 2], g[3, 3] = r ** 2, r ** 2 * sp.sin(th) ** 2
    return g


def einstein(g, phi, V0, V1, eps):
    """The five residuals R_ab - eps phi_a phi_b - V g_ab (vv, v rho, rho rho, theta theta) and eps box phi - V' from a brute-force Ricci tensor."""
    gi = g.inv()
    n = 4
    Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X4[cc]) + sp.diff(g[d, cc], X4[b]) - sp.diff(g[b, cc], X4[d])) for d in range(n)) / 2
             for cc in range(n)] for b in range(n)] for a in range(n)]

    def ric(b, d):
        s = sum(sp.diff(Gam[a][b][d], X4[a]) - sp.diff(Gam[a][b][a], X4[d]) for a in range(n))
        s += sum(Gam[a][a][k] * Gam[k][b][d] - Gam[a][d][k] * Gam[k][b][a] for a in range(n) for k in range(n))
        return s

    dphi = [sp.diff(phi, x) for x in X4]
    E = {}
    for nm, (i, j) in {"vv": (0, 0), "vr": (0, 1), "rr": (1, 1), "thth": (2, 2)}.items():
        E[nm] = ric(i, j) - eps * dphi[i] * dphi[j] - V0 * g[i, j]
    sq = sp.sqrt(-g.det())
    box = sum(sp.diff(sq * gi[a, b] * dphi[b], X4[a]) for a in range(n) for b in range(n)) / sq
    E["scal"] = eps * box - V1
    return E


def cand_F(A, lam):
    return A - CAND["kf"] * rho * sp.diff(lam, v) / lam + CAND["kh"] * (1 / lam ** 2 - 1)


def cand_r(lam):
    return lam * (Rs if CAND["shape"] == "R" else rho)


def rate_sym():
    return K * sp.Rational(1001, 1000) if CAND["rate"] == "K(1+1e-3)" else K


def is_zero(e):
    return sp.simplify(e) == 0


# ----------------------------------------------------------------------------------------------------------------------
# derive: S0 - S5
# ----------------------------------------------------------------------------------------------------------------------
def derive_S0():
    A, V0, V1 = background(K)
    E = einstein(null_metric(A, Rs), Phi, V0, V1, CAND["eps"])
    res = {k: is_zero(e.doit()) for k, e in E.items()}
    check("S0", "the static member solves the five equations in the null chart, symbolic c (convention and potential: planted truth)", all(res.values()),
          ", ".join(f"{k}: {'0' if ok else 'nonzero'}" for k, ok in res.items()))


_E_LAM = {}


def lam_equations(Kc, tag):
    """Residuals of the candidate for an arbitrary lam(v), background coefficients (c, Kc)."""
    if tag not in _E_LAM:
        A, V0, V1 = background(Kc)
        lam = sp.Function("lam")(v)
        E = einstein(null_metric(cand_F(A, lam), cand_r(lam)), Phi, V0, V1, CAND["eps"])
        _E_LAM[tag] = ({k: sp.simplify(e.doit()) for k, e in E.items()}, lam)
    return _E_LAM[tag]


def derive_S1_S2():
    Es, lam = lam_equations(K, "af")
    zero4 = {k: Es[k] == 0 for k in ("rr", "vr", "thth", "scal")}
    check("S1", "for an arbitrary lam(v), r = lam R, F = A - 2 rho lam'/lam + 1/lam^2 - 1: E_rho rho, E_v rho, E_theta theta and the scalar equation vanish identically",
          all(zero4.values()), ", ".join(f"{k}: {'0' if ok else 'nonzero'}" for k, ok in zero4.items()))
    ld, ldd = sp.diff(lam, v), sp.diff(lam, v, 2)
    expect = 2 * (K * lam * ld - lam * ldd - ld ** 2) / ((1 + rho ** 2) * lam ** 2)
    ok = is_zero(Es["vv"] - expect)
    check("S2", "the remaining equation: E_vv = 2 (K lam lam' - lam lam'' - lam'^2) / ((1 + rho^2) lam^2), K = -2c/pi; it is zero iff (lam lam')' = K lam lam'",
          ok, "the rho dependence is the overall factor 1/(1 + rho^2)" if ok else f"E_vv - expected = {sp.simplify(Es['vv'] - expect)}")


def derive_S3():
    A, V0, V1 = background(K)
    Lp, bt = sp.symbols("L beta", positive=True)
    lam = sp.sqrt(Lp + bt * sp.exp(rate_sym() * v))
    E = einstein(null_metric(cand_F(A, lam), cand_r(lam)), Phi, V0, V1, CAND["eps"])
    res = {k: is_zero(e.doit()) for k, e in E.items()}
    check("S3", "lam^2 = L + beta e^{K v}: all five residuals vanish identically for symbolic c, L, beta", all(res.values()),
          ", ".join(f"{k}: {'0' if ok else 'nonzero'}" for k, ok in res.items()))


def derive_S4():
    A, _, _ = background(K)
    Rsq = 1 + rho ** 2
    charge = (Rsq * sp.diff(A, rho) - 2 * rho * A) / 2 + CAND["charge"] * rho
    ok_a = is_zero(charge - K)
    check("S4a", "the rate is the first integral K = (R^2 A' - 2 rho A)/2 + rho of the static equations (rho independent, equal to -2c/pi)", ok_a,
          "" if ok_a else f"charge - K = {sp.simplify(charge - K)}")
    Es, lam = lam_equations(Ks, "ck")                      # c and Ks independent: A = 1 + c R^2 + Ks (rho + R^2 arctan rho), V linear in (c, Ks)
    zero4 = all(Es[k] == 0 for k in ("rr", "vr", "thth", "scal"))
    ld, ldd = sp.diff(lam, v), sp.diff(lam, v, 2)
    expect = 2 * (Ks * lam * ld - lam * ldd - ld ** 2) / ((1 + rho ** 2) * lam ** 2)
    ok_b = zero4 and is_zero(Es["vv"] - expect)
    check("S4b", "with c and the coefficient Ks of rho + R^2 arctan rho independent (no asymptotic flatness): the same four vanish and E_vv has Ks in place of K", ok_b,
          "the rate is the coefficient Ks, whatever c" if ok_b else "failed")


def derive_S5():
    Wr = rate_sym()
    A, _, _ = background(K)
    Rb = Rs if CAND["shape"] == "R" else rho
    Qf = sp.Function("Q")(rho)
    q1 = Wr * Qf / A
    q2 = sp.diff(q1, rho).subs(sp.Derivative(Qf, rho), q1)
    Y = Rb * Qf
    Ap = sp.diff(A, rho)
    Pc = 2 / rho - Ap / A
    Qc = 1 / (1 + rho ** 2) ** 2 - Ap / (A * rho * (1 + rho ** 2)) - 2 / (A * (1 + rho ** 2)) + Wr ** 2 / A ** 2
    LY = (sp.diff(Y, rho, 2) - Pc * sp.diff(Y, rho) - Qc * Y).subs(sp.Derivative(Qf, (rho, 2)), q2).subs(sp.Derivative(Qf, rho), q1)
    res = sp.simplify(LY / Qf)
    check("S5", "the first-order part delta r = (beta/2) R e^{K(t+x)} is Y_K: #212's operator Y'' - P Y' - Q Y of the master equation annihilates R e^{Wx} at W = the candidate's rate",
          res == 0, "" if res == 0 else f"residual {sp.simplify(res)}")


def derive_S6():
    """The hand-reduced equations (the ones of the by-hand derivation in exact_mode.md) with sympy jets: same residuals as the brute-force tensor of S1, S2."""
    A, V0, V1 = background(K)
    lam = sp.Function("lam")(v)
    r_e, F_e = cand_r(lam), cand_F(A, lam)
    part = lambda ex, i, j: sp.diff(ex, *([v] * i + [rho] * j)) if (i or j) else ex  # noqa: E731
    J = {nm: part(r_e, *key) for nm, key in JORD.items()}
    J.update({nm: part(F_e, *key) for nm, key in {"F": (0, 0), "Fv": (1, 0), "Fr": (0, 1), "Frr": (0, 2)}.items()})
    J.update({nm: part(Phi, *key) for nm, key in {"pv": (1, 0), "pr": (0, 1), "pvr": (1, 1), "prr": (0, 2)}.items()})
    res = hand_residuals(J, V0, V1, CAND["eps"])
    ld, ldd = sp.diff(lam, v), sp.diff(lam, v, 2)
    expect = 2 * (K * lam * ld - lam * ldd - ld ** 2) / ((1 + rho ** 2) * lam ** 2)
    ok = {k: is_zero(res[k]) for k in ("rr", "vr", "thth", "scal")}
    ok["vv"] = is_zero(res["vv"] - expect)
    check("S6", "the hand-reduced equations (R_ab = (R2/2) h_ab - 2 H_ab/r, R_theta theta = 1 - r box r - (grad r)^2, box phi) at arbitrary lam(v) give the same residuals as the brute-force tensor: four zero, E_vv as in S2",
          all(ok.values()), ", ".join(f"{k}: {'ok' if o else 'differs'}" for k, o in ok.items()))


# ----------------------------------------------------------------------------------------------------------------------
# numerical machinery (mpmath, 50 digits)
# ----------------------------------------------------------------------------------------------------------------------
def mpf_of(q):
    q = sp.Rational(q)
    return mp.mpf(q.p) / mp.mpf(q.q)


def kappa_of(cq):
    """Surface gravity A'(rho_h)/2 of the horizon of member cq (rho_h from A = 0)."""
    cm = mpf_of(cq)
    Km = -2 * cm / mp.pi
    Af = lambda y: 1 + cm * (y * y + 1) + Km * (y + (y * y + 1) * mp.atan(y))  # noqa: E731
    try:
        rh = mp.findroot(Af, mp.mpf(0) if cq == -1 else mp.mpf("0.5"))
    except ValueError:
        return Km                                # no horizon (c > 0): the mutation leaves the rate at K
    return mp.diff(Af, rh) / 2


def rate_num(cq):
    Km = -2 * mpf_of(cq) / mp.pi
    if CAND["rate"] == "K(1+1e-3)":
        return Km * (1 + mp.mpf("1e-3"))
    if CAND["rate"] == "kappa":
        return kappa_of(cq)
    return Km


_JET = {}


def jets(expr, pt, orders):
    """Values of the partial derivatives (i, j) = (d/dv)^i (d/drho)^j of expr at pt = (v0, rho0)."""
    out = {}
    for (i, j) in orders:
        key = (id(expr), i, j)
        if key not in _JET:
            d = expr
            if i:
                d = sp.diff(d, v, i)
            if j:
                d = sp.diff(d, rho, j)
            _JET[key] = (expr, sp.lambdify((v, rho), d, "mpmath"))      # keeps expr alive so id() stays valid
        out[(i, j)] = _JET[key][1](*pt)
    return out


def hess(J, drop=False):
    """Covariant Hessian of r in the 2D metric -F dv^2 + 2 dv drho (hand-derived Christoffels: G^v_vv = F'/2, G^rho_vv = (F F' - F_v)/2, G^rho_v rho = -F'/2)."""
    F, Fv, Fr = J["F"], J["Fv"], J["Fr"]
    rv, rr = J["rv"], J["rr"]
    Hvv = J["rvv"] - (Fr / 2) * rv - (F * Fr / 2 - Fv / 2) * rr
    Hvr = J["rvr"] + (0 if drop else (Fr / 2) * rr)
    Hrr = J["rrr"]
    return Hvv, Hvr, Hrr


def hand_residuals(J, V, Vp, eps):
    """The reduced equations (hand-derived): R_ab = (R2/2) h_ab - 2 H_ab/r, R_thth = 1 - (r box r + (grad r)^2), R2 = -F_rr, box phi."""
    F, Fr, Frr = J["F"], J["Fr"], J["Frr"]
    r, rv, rr = J["r"], J["rv"], J["rr"]
    pv, pr, pvr, prr = J["pv"], J["pr"], J["pvr"], J["prr"]
    R2 = -Frr
    Hvv, Hvr, Hrr = hess(J, HAND_DROP)
    grad2 = 2 * rv * rr + F * rr ** 2
    boxr = 2 * Hvr + F * Hrr
    boxphi = (2 * r * rv * pr + r ** 2 * pvr + 2 * r * rr * pv + r ** 2 * pvr + (2 * r * rr * F + r ** 2 * Fr) * pr + r ** 2 * F * prr) / r ** 2
    return {
        "vv": (R2 / 2) * (-F) - 2 * Hvv / r - eps * pv ** 2 + V * F,
        "vr": R2 / 2 - 2 * Hvr / r - eps * pv * pr - V,
        "rr": -2 * Hrr / r - eps * pr ** 2,
        "thth": 1 - (r * boxr + grad2) - V * r ** 2,
        "scal": eps * boxphi - Vp,
    }


JORD = {"r": (0, 0), "rv": (1, 0), "rr": (0, 1), "rvv": (2, 0), "rvr": (1, 1), "rrr": (0, 2)}


def collect_jets(rexpr, Fexpr, phiexpr, pt):
    J = {}
    for nm, key in JORD.items():
        J[nm] = jets(rexpr, pt, [key])[key]
    for nm, key in {"F": (0, 0), "Fv": (1, 0), "Fr": (0, 1), "Frr": (0, 2)}.items():
        J[nm] = jets(Fexpr, pt, [key])[key]
    for nm, key in {"pv": (1, 0), "pr": (0, 1), "pvr": (1, 1), "prr": (0, 2)}.items():
        J[nm] = jets(phiexpr, pt, [key])[key]
    return J


def kretschmann(J):
    """K = R2^2 + 8 H_ab H^ab / r^2 + 4 (1 - (grad r)^2)^2 / r^4 for ds^2 = h_ab + r^2 dOmega^2."""
    F, r, rv, rr = J["F"], J["r"], J["rv"], J["rr"]
    R2 = -J["Frr"]
    Hvv, Hvr, Hrr = hess(J, KRET_DROP)
    H = mp.matrix([[Hvv, Hvr], [Hvr, Hrr]])
    hinv = mp.matrix([[0, 1], [1, F]])
    M = hinv * H * hinv
    HH = sum(M[i, j] * H[i, j] for i in range(2) for j in range(2))
    grad2 = 2 * rv * rr + F * rr ** 2
    return R2 ** 2 + 8 * HH / r ** 2 + 4 * (1 - grad2) ** 2 / r ** 4


_FAM = {}


def family_exprs(cq, Lq, bq):
    """Explicit sympy r(v, rho), F(v, rho) of the candidate at member cq, lam^2 = Lq + bq e^{rate v}; and V(phi(rho)), V'(phi) as functions of rho."""
    if (cq, Lq, bq) not in _FAM:
        _FAM[(cq, Lq, bq)] = _family_exprs(cq, Lq, bq)
    return _FAM[(cq, Lq, bq)]


def _family_exprs(cq, Lq, bq):
    cs = sp.Rational(cq)
    Kq = -2 * cs / sp.pi
    A, V0, V1 = background(Kq)
    A, V0, V1 = (e.subs(c, cs) for e in (A, V0, V1))
    Wq = rate_num(cq)
    Wsym = sp.Float(str(Wq), 50)
    lam = sp.sqrt(sp.Rational(Lq) + sp.Rational(bq) * sp.exp(Wsym * v))
    Fexpr = cand_F(A, lam)
    rexpr = cand_r(lam)
    return rexpr, Fexpr, Kq, V0, V1, Wq


def v_at_s(cq, Lq, bq, uq, Wq):
    """v with lam^2 = uq: L + b e^{W v} = u."""
    return mp.log((mpf_of(uq) - mpf_of(Lq)) / mpf_of(bq)) / Wq


def eval_family(cq, Lq, bq, uq, rho0):
    rexpr, Fexpr, Kq, V0, V1, Wq = family_exprs(cq, Lq, bq)
    v0 = v_at_s(cq, Lq, bq, uq, Wq)
    pt = (v0, mp.mpf(rho0))
    J = collect_jets(rexpr, Fexpr, Phi, pt)
    fV = sp.lambdify(rho, V0, "mpmath")
    fVp = sp.lambdify(rho, V1, "mpmath")
    return J, fV(pt[1]), fVp(pt[1]), v0, Wq


MEMBERS = [-1, sp.Rational(-3, 10), sp.Rational(-3, 2), -10, sp.Rational(1, 2)]
RHOS = ["-3", "-0.4", "0.5", "2.2"]
UVALS = {"b+": [("1", "7/10", "6/5"), ("1", "7/10", "3"), ("3/2", "1/3", "2")], "b-": [("1", "-7/10", "7/10"), ("1", "-7/10", "3/10"), ("3/2", "-1/3", "6/5")]}


# ----------------------------------------------------------------------------------------------------------------------
# check: N0 - N6
# ----------------------------------------------------------------------------------------------------------------------
def N0():
    """The hand-reduced equations equal the brute-force Ricci residuals on a generic non-solution."""
    r_g = sp.exp(v / 3) * (2 + rho ** 2 / 5 + v * rho / 7)
    F_g = 1 + v * rho / 3 + rho ** 2 / 11 + v ** 2 / 13
    p_g = rho / 3 + v / 5 + v * rho / 17
    Vf = lambda p: p ** 2 / 3 + p ** 3 / 5 + sp.Rational(1, 2)  # noqa: E731
    pp = sp.Symbol("pp")
    V0 = Vf(p_g)
    V1 = sp.diff(Vf(pp), pp).subs(pp, p_g)
    E = einstein(null_metric(F_g, r_g), p_g, V0, V1, -1)
    pt = (mp.mpf(1) / 3, mp.mpf(2) / 5)
    sub = {v: sp.Rational(1, 3), rho: sp.Rational(2, 5), th: sp.Rational(3, 7)}
    brute = {k: mp.mpf(str(sp.N(e.doit().subs(sub), 50))) for k, e in E.items()}
    J = collect_jets(r_g, F_g, p_g, pt)
    Vn = mp.mpf(str(sp.N(V0.subs(sub), 50)))
    Vpn = mp.mpf(str(sp.N(V1.subs(sub), 50)))
    hand = hand_residuals(J, Vn, Vpn, -1)
    scale = max(abs(x) for x in brute.values())
    diff = max(abs(hand[k] - brute[k]) for k in brute) / scale
    check("N0", "hand-reduced equations against the brute-force Ricci residuals on a generic metric, five components (bar 1e-30; the metric is not a solution: residuals are O(1))",
          diff < mp.mpf("1e-30") and scale > 1e-3, f"{mp.nstr(diff, 3)}; largest residual {mp.nstr(scale, 3)}")


def N1():
    """The five residuals of the exact family, hand-reduced equations, 50 digits, members c = -1, -3/10, -3/2, -10, 1/2."""
    worst, n = mp.mpf(0), 0
    per = {}
    for cq in MEMBERS:
        Kq = -2 * mpf_of(cq) / mp.pi
        mw = mp.mpf(0)
        for (Lq, bq, uq) in UVALS["b+"] + UVALS["b-"]:
            for r0 in RHOS:
                J, V, Vp, v0, Wq = eval_family(cq, Lq, bq, uq, r0)
                if abs(Kq) < 1e-30:
                    continue
                res = hand_residuals(J, V, Vp, CAND["eps"])
                scale = max(max(abs(x) for x in (J["F"], J["Fr"], J["r"], V, Vp)), 1)
                m = max(abs(x) for x in res.values()) / scale
                mw = max(mw, m)
                n += 1
        per[cq] = mw
        worst = max(worst, mw)
    check("N1", "exact family (lam^2 = L + beta e^{K v}, L = 1, 3/2, beta of both signs, lam^2 from 0.3 to 3.0): the five hand-reduced residuals at 50 digits, relative to the largest term (bar 1e-30)",
          worst < mp.mpf("1e-30"), f"worst {mp.nstr(worst, 3)} over {n} points; " + ", ".join(f"c = {cq}: {mp.nstr(m, 2)}" for cq, m in per.items()))


def build_static():
    """The linearized field equations in the static gauge from a brute-force Ricci tensor, with Y, H0, H2 the gauge image of the candidate's first order."""
    t = sp.Symbol("t", real=True)
    e = sp.Symbol("e")
    Ww, kf, kh = sp.symbols("Ww kf kh")
    A, Rb, Pb, V0b, V1b, Qf, Xi = (sp.Function(nm)(rho) for nm in ("A", "Rb", "Pb", "V0b", "V1b", "Q", "Xi"))
    T = sp.exp(Ww * t)
    D = kf * Ww * rho / 2 + kh
    Yf = (Rb if CAND["shape"] == "R" else rho) * Qf * T
    H2 = 2 * D * Qf * T / A
    H0 = T * (-2 * D * Qf / A + 4 * Ww * Xi)
    X = [t, rho, th, ph]
    g = sp.diag(-A * (1 + e * H0), (1 + e * H2) / A, (Rb + e * Yf) ** 2, (Rb + e * Yf) ** 2 * sp.sin(th) ** 2)
    gi = g.inv()
    n = 4
    Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X[cc]) + sp.diff(g[d, cc], X[b]) - sp.diff(g[b, cc], X[d])) for d in range(n)) / 2
             for cc in range(n)] for b in range(n)] for a in range(n)]

    def ric(b, d):
        s = sum(sp.diff(Gam[a][b][d], X[a]) - sp.diff(Gam[a][b][a], X[d]) for a in range(n))
        s += sum(Gam[a][a][k] * Gam[k][b][d] - Gam[a][d][k] * Gam[k][b][a] for a in range(n) for k in range(n))
        return s

    lin = lambda ex: sp.diff(ex, e).subs(e, 0)  # noqa: E731
    eps = CAND["eps"]
    dphi = [sp.diff(Pb, x) for x in X]
    Eq = {}
    for nm, (i, j) in {"tt": (0, 0), "tr": (0, 1), "rr": (1, 1), "thth": (2, 2)}.items():
        Eq[nm] = lin(ric(i, j) - eps * dphi[i] * dphi[j] - V0b * g[i, j])
    sq = sp.sqrt(-g.det())
    box = sum(sp.diff(sq * gi[a, b] * dphi[b], X[a]) for a in range(n) for b in range(n)) / sq
    Eq["scal"] = lin(eps * box - V1b)
    q1 = Ww * Qf / A
    q2 = sp.diff(q1, rho).subs(sp.Derivative(Qf, rho), q1)
    x1 = D * Qf / A ** 2
    x2 = sp.diff(x1, rho).subs(sp.Derivative(Qf, rho), q1)
    out = {}
    for k, ex in Eq.items():
        ex = ex.doit().subs({sp.Derivative(Qf, (rho, 2)): q2, sp.Derivative(Xi, (rho, 2)): x2}).subs({sp.Derivative(Qf, rho): q1, sp.Derivative(Xi, rho): x1})
        out[k] = ex.subs(t, 0)
    return out, (A, Rb, Pb, V0b, V1b, Qf, Xi, Ww, kf, kh)


def N2():
    """Static gauge delta g_{t rho} = 0, delta phi = 0: the linearized equations at the gauge image of the candidate's first-order perturbation."""
    out, (A, Rb, Pb, V0b, V1b, Qf, Xi, Ww, kf, kh) = build_static()
    abstract = ("A", "Rb", "Pb", "V0b", "V1b")
    Aex, V0ex, V1ex = background(K)
    explicit = {"A": Aex, "Rb": Rs, "Pb": Phi, "V0b": V0ex, "V1b": V1ex}
    XiV = sp.Symbol("XiV")
    comp = {}
    for k, ex in out.items():                       # compile each equation once: terms as a list of mpmath functions of the jets
        reps, keys = {}, []
        for d in sorted(ex.atoms(sp.Derivative), key=lambda z: -z.derivative_count):
            f = d.expr
            if f.func.__name__ in abstract:
                sym = sp.Symbol(f"{f.func.__name__}_{d.derivative_count}")
                reps[d] = sym
                keys.append((f.func.__name__, d.derivative_count, sym))
        for f in ex.atoms(sp.Function):
            if f.func.__name__ in abstract and f not in reps:
                sym = sp.Symbol(f"{f.func.__name__}_0")
                reps[f] = sym
                keys.append((f.func.__name__, 0, sym))
        e2 = ex.xreplace(reps).xreplace({Qf: sp.Integer(1), Xi: XiV}).subs(th, sp.Rational(1, 3))
        terms = sp.Add.make_args(sp.expand(e2))
        assert not any(t.has(Qf) or t.has(Xi) for t in terms), k
        comp[k] = (sp.lambdify([sym for _, _, sym in keys] + [rho, Ww, kf, kh, XiV], list(terms), "mpmath"), keys)
    worst, n = mp.mpf(0), 0
    per = {}
    for cq in [-1, sp.Rational(-3, 10), sp.Rational(-3, 2), -10]:
        cs = sp.Rational(cq)
        Wq = rate_num(cq)
        mw = mp.mpf(0)
        fcache = {}
        for r0 in ("-3", "-1.5", "-0.4", "0.3", "0.8", "1.7", "4", "9"):
            x0 = mp.mpf(r0)
            Aval = sp.lambdify(rho, Aex.subs(c, cs), "mpmath")(x0)
            if abs(Aval) < mp.mpf("1e-3"):
                continue
            for xi_val in ("0", "3/7"):
                tot = mp.mpf(0)
                for k, (ft, keys) in comp.items():
                    vals = []
                    for nm, kk, sym in keys:
                        if (nm, kk) not in fcache:
                            fx = explicit[nm].subs(c, cs)
                            fcache[(nm, kk)] = sp.lambdify(rho, sp.diff(fx, rho, kk) if kk else fx, "mpmath")
                        vals.append(fcache[(nm, kk)](x0))
                    tv = ft(*vals, x0, Wq, mpf_of(CAND["kf"]), mpf_of(CAND["kh"]), mpf_of(xi_val))
                    s_abs = sum(abs(z) for z in tv)
                    tot = max(tot, abs(sum(tv)) / max(s_abs, mp.mpf("1e-40")))
                mw = max(mw, tot)
                n += 1
        per[cq] = mw
        worst = max(worst, mw)
    check("N2", "static gauge, brute-force linearized E_tt, E_t rho, E_rho rho, E_theta theta and the scalar equation at Y = R Q T, H2 = 2D Q T/A, H0 = T (-2D Q/A + 4 W Xi) (the gauge image of the candidate, Xi free): residuals relative to the sum of the terms, 50 digits, members c = -1, -3/10, -3/2, -10, 8 points each, Xi = 0 and 3/7 (bar 1e-25)",
          worst < mp.mpf("1e-25"), f"worst {mp.nstr(worst, 3)} over {n} point-gauge pairs; " + ", ".join(f"c = {cq}: {mp.nstr(m, 2)}" for cq, m in per.items()))


def N3():
    """Planted truth for the curvature routine: Schwarzschild, F = 1 - 2M/rho, r = rho, K = 48 M^2/r^6."""
    Mv = sp.Rational(3, 2)
    rexpr, Fexpr = rho, 1 - 2 * Mv / rho
    pt = (mp.mpf(0), mp.mpf(4))
    J = collect_jets(rexpr, Fexpr, sp.Integer(0) * rho, pt)
    kr = kretschmann(J)
    exact = 48 * mpf_of(Mv) ** 2 / pt[1] ** 6
    check("N3", "the Kretschmann routine on Schwarzschild (M = 3/2, r = 4): 48 M^2/r^6 (bar 1e-40)", abs(kr - exact) / exact < mp.mpf("1e-40"), f"{mp.nstr(kr, 15)} against {mp.nstr(exact, 15)}")


def N4():
    """beta < 0: collapse at the finite advanced time v*, a curvature singularity on the null hypersurface v = v* at every rho (K lam^4 -> C(rho))."""
    ok_all, rows = True, []
    for cq in (-1, sp.Rational(-3, 2), -10):
        Wq = rate_num(cq)
        for r0 in ("-2", "0.5", "3"):
            C = []
            for uq in ("1/10000", "1/1000000", "1/100000000"):
                J, V, Vp, v0, _ = eval_family(cq, "1", "-7/10", uq, r0)
                C.append(kretschmann(J) * mpf_of(uq) ** 2)
            vs = mp.log(1 / mp.mpf("0.7")) / Wq
            kr = kretschmann(eval_family(cq, "1", "-7/10", "1/100000000", r0)[0])
            ok = (C[2] > 0) and abs(C[2] - C[1]) / C[2] < mp.mpf("1e-3") and kr > mp.mpf("1e12")
            ok_all &= ok
            rows.append(f"c = {cq}, rho = {r0}: K lam^4 = {mp.nstr(C[1], 5)}, {mp.nstr(C[2], 5)}; K = {mp.nstr(kr, 3)} at lam^2 = 1e-8; v* = {mp.nstr(vs, 6)}")
    check("N4", "beta = -7/10, L = 1: K_Kretschmann * lam^4 converges (1e-3) to a positive limit as lam^2 -> 0 (1e-6, 1e-8) at c = -1, -3/2, -10 and rho = -2, 1/2, 3, and K > 1e12 at lam^2 = 1e-8",
          ok_all, "; ".join(rows[:3]) + f"; ... ({len(rows)} cases)")
    for cq in (-1, sp.Rational(-3, 2)):                         # the line rho = -1/K, where F = A for all v (info)
        Kq = -2 * mpf_of(cq) / mp.pi
        for uq in ("1/1000000", "1/100000000"):
            J, V, Vp, v0, _ = eval_family(cq, "1", "-7/10", uq, mp.nstr(-1 / Kq, 40))
            print(f"  info: c = {cq}, rho = -1/K = {mp.nstr(-1 / Kq, 8)}, lam^2 = {uq}: F = {mp.nstr(J['F'], 10)}, K_Kretschmann = {mp.nstr(kretschmann(J), 8)}")
    return rows


def N5():
    """beta > 0: no singularity, r -> infinity as e^{K v/2}; late slices are everywhere trapped-like (F < 0) for c < 0; curvature has a finite limit."""
    ok_all, rows = True, []
    for cq in (-1, sp.Rational(-3, 2), -10, sp.Rational(-3, 10)):
        Wq = rate_num(cq)
        for r0 in ("-2", "0.5", "3"):
            Ks_ = []
            for uq in ("1000001", "1000000000001"):
                J, V, Vp, v0, _ = eval_family(cq, "1", "7/10", uq, r0)
                Ks_.append(kretschmann(J))
                last = J
            grow = last["rv"] / last["r"]
            ok = abs(Ks_[1] - Ks_[0]) / abs(Ks_[1]) < mp.mpf("1e-5") and last["F"] < 0 and abs(grow - Wq / 2) / (Wq / 2) < mp.mpf("1e-10")
            ok_all &= ok
            rows.append(f"c = {cq}, rho = {r0}: K = {mp.nstr(Ks_[1], 6)}, F = {mp.nstr(last['F'], 5)}, d ln r/dv = {mp.nstr(grow, 8)} (K/2 = {mp.nstr(Wq / 2, 8)})")
    check("N5", "beta = 7/10, L = 1: from lam^2 = 1e6 to 1e12 the Kretschmann invariant changes by < 1e-5 (finite limit), F < 0, and d ln r/dv = K/2 to 1e-10, at c = -1, -3/2, -10, -3/10 and rho = -2, 1/2, 3",
          ok_all, "; ".join(rows[:2]) + f"; ... ({len(rows)} cases)")
    return rows


def N6():
    """beta = 0, L != 1: the static family r = sqrt(L) R, F = A + 1/L - 1 solves the equations; Y = R at W = 0 is not a solution of #212's master equation (info)."""
    worst = mp.mpf(0)
    for cq in MEMBERS:
        for (Lq, r0) in (("3/2", "0.5"), ("1/3", "-2"), ("5", "3")):
            rexpr, Fexpr, Kq, V0, V1, Wq = family_exprs(cq, Lq, "0")
            pt = (mp.mpf(0), mp.mpf(r0))
            J = collect_jets(rexpr, Fexpr, Phi, pt)
            res = hand_residuals(J, sp.lambdify(rho, V0, "mpmath")(pt[1]), sp.lambdify(rho, V1, "mpmath")(pt[1]), CAND["eps"])
            worst = max(worst, max(abs(x) for x in res.values()) / max(abs(J["F"]), 1))
    check("N6", "static family L != 1 (r = sqrt(L) R, F = A + 1/L - 1, the same potential): five residuals at 50 digits, 5 members x 3 values of L (bar 1e-30)", worst < mp.mpf("1e-30"),
          f"worst {mp.nstr(worst, 3)}")
    # info: Y = R, W = 0 in #212's operator
    for cq in (-1, sp.Rational(-3, 10)):
        cs = sp.Rational(cq)
        Af = background(-2 * cs / sp.pi)[0].subs(c, cs)
        Ap = sp.diff(Af, rho)
        Pc = 2 / rho - Ap / Af
        Qc = 1 / (1 + rho ** 2) ** 2 - Ap / (Af * rho * (1 + rho ** 2)) - 2 / (Af * (1 + rho ** 2))
        Y = sp.sqrt(1 + rho ** 2)
        val = (sp.diff(Y, rho, 2) - Pc * sp.diff(Y, rho) - Qc * Y).subs(rho, 1)
        print(f"  info: #212's operator on Y = R at W = 0, c = {cq}, rho = 1: {sp.N(val, 15)}  (nonzero: its derivation divides by W through E_t rho)")


NUMERIC = [("N0", N0), ("N1", N1), ("N2", N2), ("N3", N3), ("N4", N4), ("N5", N5), ("N6", N6)]
DERIVE = [("S0", derive_S0), ("S1", derive_S1_S2), ("S3", derive_S3), ("S4", derive_S4), ("S5", derive_S5), ("S6", derive_S6)]

# name: checks that must fire (development runs of negative-mass-raw I269-mac/dev, before they were declared)
MUTATIONS = {
    "kf": ["S1", "S2", "S3", "S4b", "S6", "N1", "N2"],
    "kh": ["S1", "S2", "S3", "S4b", "S6", "N1", "N2", "N6"],
    "rate": ["S3", "S5", "N1", "N2"],
    "kappa": ["N1", "N2"],
    "shape": ["S1", "S2", "S3", "S4b", "S5", "S6", "N1", "N2", "N4", "N6"],
    "naive": ["S1", "S2", "S3", "S4b", "S6", "N1", "N2", "N4", "N5", "N6"],
    "pot": ["S0", "S1", "S2", "S3", "S4b", "S6", "N1", "N2", "N6"],
    "eps": ["S0", "S1", "S3", "S4b", "S6", "N1", "N6"],
    "charge": ["S4a"],
    "hand_drop": ["S6", "N0", "N1", "N6"],
    "kret_drop": ["N3"],
}


def run_all(items):
    t0 = time.time()
    for cid, fn in items:
        if want(cid) or (cid == "S1" and want("S2")) or (cid == "S4" and (want("S4a") or want("S4b"))):
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
    elif sub == "mutations":
        args = [a for a in sys.argv[2:] if not a.startswith("--")]
        sys.exit(0 if report_mutations(args[0], args[1] if len(args) > 1 else "") else 1)
    else:
        sys.exit(f"unknown command {sub}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
