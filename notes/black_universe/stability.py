#!/usr/bin/env python3
"""Radial (l = 0) stability of the phantom "black universes" of Bronnikov and Fabris [BF06, gr-qc/0511109],
the family of Bronnikov, Konoplya and Zhidenko [BKZ12, arXiv:1205.2224].  See README.md in this directory.

    stability.py check    [--mutate=M]     background field equations, derived master equation, comparison with [BKZ12]
    stability.py planted  [--mutate=M]     Ellis wormhole (c = 0) through the same code: one growing mode, W^2 a^2
    stability.py count c [--res=R1|R2] [--out=FILE.json] [--mutate=M]    growing modes of one member
    stability.py report DIR                table, comparison with [BKZ12], from DIR/member_*.json
    stability.py spectral DIR              independent Chebyshev eigenvalue check of the counts in DIR/member_*R2.json
    stability.py controls                  runs every --mutate=M and asserts that its target check fires

Units: b = min R = 1.  Family parameter c; 3m = -2c/pi (asymptotic flatness at rho = +inf, [BKZ12] (5)).
Exit status 1 when a check fails (a mutated run is meant to exit 1).
"""
import json
import math
import re
import subprocess
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = __import__("os").path.dirname(__import__("os").path.abspath(__file__))

# ------------------------------------------------------------------------------------------------ declared settings
PLANTED_FILE = HERE + "/../wormhole_largeD/w2_output.txt"      # crate (qd) value of W^2 a^2 at D = 4, row "TS13 Table I"
PLANTED_TOL = {"R1": 1e-7, "R2": 1e-9}                          # absolute, on W^2 a^2, declared before the runs of record
W_LO, W_HI = 0.01, 8.0                                          # range searched, units 1/b
RES = {
    "R1": dict(rtol=1e-9, z0=1e-9, rmax0=150.0, rmaxW=30.0, delta=0.3, nW=60),
    "R2": dict(rtol=1e-12, z0=1e-13, rmax0=300.0, rmaxW=45.0, delta=0.2, nW=150),
}
ROOT_AGREE = 1e-6                                               # relative, R1 against R2
MEMBERS = [-10.0, -3.0, -2.0, -1.5, -1.1, -1.01, -1.001, -1.0001, -1.0, -0.9999, -0.999, -0.99, -0.95, -0.9, -0.7, -0.5,
           -0.3, -0.2, -0.1, -0.01]
# what [BKZ12] shows for a member (Sec. V.C, Figs 5-8): "time domain" = a time-domain profile that grows (or decays at c = -1),
# "potential" = only the effective potential with a negative region and the statement of Conclusion 2 (stable only at c = -1)
BKZ_EVIDENCE = {-0.1: "time domain, Fig 7", -0.2: "time domain, Fig 7", -0.3: "time domain, Fig 7", -0.9: "time domain, Fig 8",
                -0.95: "time domain, Fig 8", -0.99: "time domain, Fig 8", -1.0: "time domain, Fig 5 (decays)",
                -1.001: "time domain, Fig 5", -1.01: "potential, Fig 5", -1.5: "potential, Fig 6", -2.0: "potential, Fig 6",
                -3.0: "potential, Fig 6"}
MUT = set()


def mutating(name):
    return name in MUT


def planted_value():
    for line in open(PLANTED_FILE):
        m = re.match(r"\s+4\s+TS13 Table I\s+\S+\s+(\S+)", line)
        if m:
            return m.group(1)
    raise RuntimeError("planted value not found in " + PLANTED_FILE)


# ------------------------------------------------------------------------------------------------ symbolic checks
def sym_family(mut=()):
    """BKZ12 (3), (4) in closed form (b symbolic). Returns sympy objects."""
    import sympy as sp
    rho, b = sp.symbols("rho b", positive=True)
    c, m = sp.symbols("c m", real=True)
    R = sp.sqrt(rho**2 + b**2)
    T = sp.atan(rho / b)
    m_A = 2 * m if "bg_A" in mut else 3 * m                       # control: wrong coefficient in A
    A = R**2 * (c / b**2 + 1 / R**2 + m_A / b**3 * (b * rho / R**2 + T))
    phi = sp.sqrt(2) * T
    psi = sp.Symbol("psi")
    s3 = 2 if "bg_V" in mut else 3                                # control: wrong coefficient in V
    Vpsi = -(c / b**2) * (3 - 2 * sp.cos(psi)**2) - 3 * m / b**3 * (s3 * sp.sin(psi) * sp.cos(psi) + psi * (3 - 2 * sp.cos(psi)**2))
    return dict(sp=sp, rho=rho, b=b, c=c, m=m, R=R, T=T, A=A, phi=phi, psi=psi, Vpsi=Vpsi)


def check_background_symbolic(mut=()):
    """[BF06] (3)-(6), eps = -1, on the closed forms: residuals must be 0."""
    F = sym_family(mut)
    sp, rho, R, A, phi, T, psi, Vpsi = (F[k] for k in ("sp", "rho", "R", "A", "phi", "T", "psi", "Vpsi"))
    eps = -1
    V = Vpsi.subs(psi, T)
    Vphi = sp.diff(Vpsi, psi).subs(psi, T) / sp.sqrt(2)
    Ts = sp.Symbol("Ts")
    eqs = {
        "(4) (A' R^2)' + 2 R^2 V": sp.diff(sp.diff(A, rho) * R**2, rho) + 2 * R**2 * V,
        "(5) 2 R''/R + eps phi'^2": 2 * sp.diff(R, rho, 2) / R + eps * sp.diff(phi, rho)**2,
        "(6) A (R^2)'' - R^2 A'' - 2": A * sp.diff(R**2, rho, 2) - R**2 * sp.diff(A, rho, 2) - 2,
        "(3) eps (A R^2 phi')' - R^2 V_phi": eps * sp.diff(A * R**2 * sp.diff(phi, rho), rho) - R**2 * Vphi,
        "theta-theta: 1 - A (R R')' - R^2 V - R A' R'": 1 - A * sp.diff(R * sp.diff(R, rho), rho) - R**2 * V - R * sp.diff(A, rho) * sp.diff(R, rho),
    }
    out = {}
    for k, e in eqs.items():
        e = sp.simplify(e.doit().subs(sp.atan(rho / F["b"]), Ts))
        out[k] = e
    return out


def family_mp(c, mut=()):
    """closed-form background in mpmath, b = 1, 3m = -2c/pi; returns dict of callables of rho."""
    import mpmath as mp
    K = -2 * mp.mpf(c) / mp.pi
    cc = mp.mpf(c)
    Km = (2 if "bg_A" in mut else 3) * K / 3                      # control: same wrong coefficient as in sym_family
    A = lambda r: 1 + cc * (r * r + 1) + Km * (r + (r * r + 1) * mp.atan(r))
    V = lambda r: -cc * (1 + 2 * r * r / (r * r + 1)) - K * ((2 if "bg_V" in mut else 3) * r / (r * r + 1) + (1 + 2 * r * r / (r * r + 1)) * mp.atan(r))
    return dict(A=A, V=V, K=K)


def check_background_numeric(mut=(), cs=(-3.0, -1.0, -0.3, 0.0, 0.5), pts=(-3.0, -1.0, -0.3, 0.2, 0.9, 2.5, 7.0)):
    """the same equations at points, 40 digits, derivatives by mpmath; also A'(h) = 2 (rho_0 - h)/R^2 at the horizon."""
    import mpmath as mp
    mp.mp.dps = 40
    worst = mp.mpf(0)
    for c in cs:
        f = family_mp(c, mut)
        A, V = f["A"], f["V"]
        for r0 in pts:
            r0 = mp.mpf(r0)
            R = lambda r: mp.sqrt(r * r + 1)
            Ap = mp.diff(A, r0); App = mp.diff(A, r0, 2)
            Rr, Rp, Rpp = R(r0), mp.diff(R, r0), mp.diff(R, r0, 2)
            phip = mp.sqrt(2) / (r0 * r0 + 1)
            Vv = V(r0)
            Vphi = mp.diff(V, r0) / phip
            tt = App * Rr**2 + 2 * Ap * Rr * Rp + 2 * Rr**2 * Vv
            r5 = 2 * Rpp / Rr - phip**2
            sc = -mp.diff(lambda r: A(r) * R(r)**2 * mp.sqrt(2) / (r * r + 1), r0) - Rr**2 * Vphi
            for e in (tt, r5, sc):
                worst = max(worst, abs(e))
    # horizon surface gravity: A'(h) = 2 (rho_0 - h)/R(h)^2, rho_0 = 3m = K  ([BF06] (7); (8) as printed lacks the 2)
    horizon_err = mp.mpf(0)
    for c in (-3.0, -1.5, -0.3, -0.1):
        f = family_mp(c, mut)
        A = f["A"]
        try:
            h = mp.findroot(A, (mp.mpf(-1000), f["K"]), solver="anderson", tol=1e-35)
        except ValueError:
            horizon_err = mp.inf                                     # a mutated background has no horizon in the bracket
            break
        lhs = mp.diff(A, h)
        rhs = 2 * (f["K"] - h) / (h * h + 1)
        if "bg_kappa" in mut:
            rhs = rhs / 2                                          # control: the printed (8), without the 2
        horizon_err = max(horizon_err, abs(lhs - rhs))
    return worst, horizon_err


def derive_master(mut=()):
    """Linearize R_mn - eps d_m phi d_n phi - V g_mn = 0 and the scalar equation about the static background (generic A, R, phi0),
    perturbation ~ exp(W t). Returns residual dictionary and the derived master equations."""
    import sympy as sp
    t, rho, th, ph, eta, W = sp.symbols("t rho theta varphi eta W")
    eps = sp.Symbol("epsilon")
    A = sp.Function("A")(rho); R = sp.Function("R")(rho); P0 = sp.Function("P0")(rho)
    Vs = sp.Function("Vs")(rho); Vp = sp.Function("Vp")(rho); Vpp = sp.Function("Vpp")(rho)
    n1, l1, y1, f1 = [sp.Function(s)(t, rho) for s in ("n1", "l1", "y1", "f1")]
    X = [t, rho, th, ph]
    nu = sp.log(A) / 2 + eta * n1; lam = -sp.log(A) / 2 + eta * l1; r = R + eta * y1; phi = P0 + eta * f1
    g = sp.diag(-sp.exp(2 * nu), sp.exp(2 * lam), r**2, r**2 * sp.sin(th)**2)
    gi = g.inv()
    G = [[[sum(gi[a, d] * (sp.diff(g[d, b_], X[c_]) + sp.diff(g[d, c_], X[b_]) - sp.diff(g[b_, c_], X[d])) for d in range(4)) / 2
           for c_ in range(4)] for b_ in range(4)] for a in range(4)]
    def Ric(b_, c_):
        s = 0
        for a in range(4):
            s += sp.diff(G[a][b_][c_], X[a]) - sp.diff(G[a][b_][a], X[c_])
            for d in range(4):
                s += G[a][a][d] * G[d][b_][c_] - G[a][c_][d] * G[d][b_][a]
        return s
    Vexp = Vs + eta * Vp * f1
    dphi = [sp.diff(phi, x) for x in X]
    E = {(a, b_): Ric(a, b_) - eps * dphi[a] * dphi[b_] - Vexp * g[a, b_] for (a, b_) in [(0, 0), (0, 1), (1, 1), (2, 2)]}
    sq = sp.sqrt(-g.det())
    E["S"] = eps * sum(sp.diff(sq * gi[a, a] * dphi[a], X[a]) for a in range(4)) / sq - (Vp + eta * Vpp * f1)
    Nn, Ll, Yy, Ff = [sp.Function(s)(rho) for s in "NLYF"]
    Rpp_bg = -eps * sp.diff(P0, rho)**2 * R / 2
    Ar2 = (A * (2 * sp.diff(R, rho)**2 + 2 * R * Rpp_bg) - 2) / R**2
    def bg(e):
        e = e.subs(Vpp, sp.diff(Vp, rho) / sp.diff(P0, rho)) if e.has(Vpp) else e
        e = e.subs(Vp, eps * (A * sp.diff(P0, rho, 2) + 2 * A * sp.diff(P0, rho) * sp.diff(R, rho) / R + sp.diff(A, rho) * sp.diff(P0, rho))).doit()
        e = e.subs(Vs, -(sp.diff(A, rho, 2) * R**2 + 2 * sp.diff(A, rho) * R * sp.diff(R, rho)) / (2 * R**2))
        for _ in range(4):
            e = e.subs(sp.diff(R, rho, 4), sp.diff(Rpp_bg, rho, 2)).subs(sp.diff(R, rho, 3), sp.diff(Rpp_bg, rho)).subs(sp.diff(R, rho, 2), Rpp_bg)
            e = e.subs(sp.diff(A, rho, 4), sp.diff(Ar2, rho, 2)).subs(sp.diff(A, rho, 3), sp.diff(Ar2, rho)).subs(sp.diff(A, rho, 2), Ar2).doit()
        return sp.simplify(e)
    res = {}
    for k in [(0, 0), (0, 1), (1, 1), (2, 2), "S"]:
        res["background " + str(k)] = bg(E[k].subs(eta, 0).doit())
    def lin(k, **kw):
        e1 = sp.diff(E[k], eta).subs(eta, 0).doit()
        e1 = e1.subs({n1: kw.get("N", Nn) * sp.exp(W * t), l1: kw.get("L", Ll) * sp.exp(W * t), y1: kw.get("Y", Yy) * sp.exp(W * t),
                      f1: kw.get("F", Ff) * sp.exp(W * t)}).doit()
        return bg(sp.simplify(e1 * sp.exp(-W * t)))
    Rp, Ap, Pp = sp.diff(R, rho), sp.diff(A, rho), sp.diff(P0, rho)
    zero = sp.Integer(0)
    # --- gauge dphi = 0 (G1)
    E1 = {k: lin(k, F=zero) for k in [(0, 0), (0, 1), (1, 1), (2, 2), "S"]}
    Ls = (sp.diff(Yy, rho) - Ap * Yy / (2 * A)) / Rp                               # from the t-rho equation
    Np = sp.Symbol("Np")
    Npsol = sp.solve(E1["S"].subs(sp.diff(Nn, rho), Np), Np)[0]
    e22 = bg(E1[(2, 2)].subs(sp.diff(Nn, rho), Np).subs(Np, Npsol).subs(Ll, Ls).doit())
    Ypp = sp.simplify(sp.solve(e22, sp.diff(Yy, rho, 2))[0])
    if "master_p" in mut:                                                          # control: a wrong coefficient of Y'
        Ypp = Ypp + Rp / R * sp.diff(Yy, rho)
    tr = {}
    for key in [(0, 0), (1, 1)]:
        e = E1[key].subs(sp.diff(Nn, rho, 2), sp.diff(Npsol, rho)).subs(sp.diff(Nn, rho), Npsol).subs(Ll, Ls).doit()
        for _ in range(3):
            e = e.subs(sp.diff(Yy, rho, 3), sp.diff(Ypp, rho)).subs(sp.diff(Yy, rho, 2), Ypp).subs(Ll, Ls).doit()
        tr[key] = bg(e)
    res["G1 t-t residual after the master equation"] = tr[(0, 0)]
    res["G1 rho-rho residual after the master equation"] = tr[(1, 1)]
    # --- gauge dR = 0 ([BKZ12] gauge, BK)
    E2 = {k: lin(k, Y=zero) for k in [(0, 0), (0, 1), (1, 1), (2, 2), "S"]}
    Lb = eps * Pp * R * Ff / (2 * Rp)
    e22 = E2[(2, 2)].subs(sp.diff(Nn, rho), Np).subs(Ll, Lb).doit()
    Npb = sp.solve(e22, Np)[0]
    S2 = bg(E2["S"].subs(sp.diff(Nn, rho), Np).subs(Np, Npb).subs(Ll, Lb).doit())
    Fpp = sp.simplify(sp.solve(S2, sp.diff(Ff, rho, 2))[0])
    cF1 = sp.simplify(Fpp.coeff(sp.diff(Ff, rho)))
    cF0 = sp.simplify(sp.expand(Fpp - cF1 * sp.diff(Ff, rho)).coeff(Ff))
    Ut = cF0 - W**2 / A**2
    for key in [(0, 0), (1, 1)]:
        e = E2[key].subs(sp.diff(Nn, rho, 2), sp.diff(Npb, rho)).subs(sp.diff(Nn, rho), Npb).subs(Ll, Lb).doit()
        for _ in range(3):
            e = e.subs(sp.diff(Ff, rho, 3), sp.diff(Fpp, rho)).subs(sp.diff(Ff, rho, 2), Fpp).subs(Ll, Lb).doit()
        res["BK %s residual after the master equation" % (key,)] = bg(e)
    # --- comparison with [BKZ12] (40), (41):  A * Ut = eps phi'^2/R'^2 (R^2 V - 1) + 2 phi' R V_phi/R' + eps V_phiphi
    sgn = -1 if "bkz41_sign" in mut else 1                                         # control: sign of the V_phi term
    rhs41 = eps * Pp**2 / Rp**2 * (R**2 * Vs - 1) + sgn * 2 * Pp * R * Vp / Rp + eps * Vpp
    res["[BKZ12] (41): A*Utilde - rhs"] = bg((A * Ut - rhs41).subs(eps**2, 1)).subs(eps, -1)
    res["  ... and eps = +1"] = bg((A * Ut - rhs41).subs(eps**2, 1)).subs(eps, 1)
    # --- gauge map F = -phi0' Y / R' maps the dphi = 0 master equation onto the dR = 0 one
    Ym = sp.Function("Y")(rho)
    Fm = -Pp * Ym / Rp
    ex = (sp.diff(Fm, rho, 2) - Fpp.subs(Ff, Fm).subs(sp.diff(Ff, rho), sp.diff(Fm, rho))).doit()
    Yp2 = Ypp.subs(Yy, Ym)
    for _ in range(3):
        ex = ex.subs(sp.diff(Ym, rho, 3), sp.diff(Yp2, rho)).subs(sp.diff(Ym, rho, 2), Yp2).doit()
    res["gauge map F = -phi0' Y/R' residual"] = bg(ex)
    # --- the derived coefficients for the numerical code, family version (b = 1), against the generic ones, at 30 digits
    return res, dict(Ypp=Ypp, Fpp=Fpp, Ut=Ut, cF1=cF1, A=A, R=R, P0=P0, rho=rho, W=W, eps=eps, Yy=Yy, Vs=Vs, Vp=Vp, Vpp=Vpp)


def family_vs_generic(D, mut=()):
    """the simplified family coefficients used by the numerics, against the generic derived Y'' equation, mpmath 30 digits."""
    import sympy as sp
    import mpmath as mp
    mp.mp.dps = 30
    rho = D["rho"]; A = D["A"]; R = D["R"]; P0 = D["P0"]; W = D["W"]; Yy = D["Yy"]
    c = sp.Rational(-3, 10); K = -2 * c / sp.pi
    Af = 1 + c * (rho**2 + 1) + K * (rho + (rho**2 + 1) * sp.atan(rho))
    Rf = sp.sqrt(rho**2 + 1); Pf = sp.sqrt(2) * sp.atan(rho)
    Ypp = D["Ypp"].subs(D["eps"], -1)
    Y0, Y1 = sp.symbols("Y0 Y1")
    ex = Ypp.subs(sp.diff(Yy, rho), Y1).subs(Yy, Y0)
    ex = ex.subs({A: Af, R: Rf, P0: Pf}).doit()
    # family version (rho-form): Y'' = (2/rho - A'/A) Y' + (1/R^4 - A'/(A rho R^2) - 2/(A R^2)) Y + W^2 Y/A^2
    q_extra = 0 if "family_q" not in mut else 1 / (Af * Rf**2)                       # control: wrong term in the Y coefficient
    fam = (2 / rho - sp.diff(Af, rho) / Af) * Y1 + (1 / Rf**4 - sp.diff(Af, rho) / (Af * rho * Rf**2) - 2 / (Af * Rf**2) + q_extra + W**2 / Af**2) * Y0
    d = sp.lambdify((rho, Y0, Y1, W), ex - fam, "mpmath")
    worst = max(abs(d(mp.mpf(r), mp.mpf("0.37"), mp.mpf("-0.81"), mp.mpf("0.7"))) for r in ("-2.3", "-0.6", "0.35", "1.4", "5.0"))
    return worst


def check_command():
    ok_all = True
    def line(tag, ok, msg):
        nonlocal ok_all
        ok_all &= ok
        print(("PASS  " if ok else "FAIL  ") + tag + ": " + msg, flush=True)
    t0 = time.time()
    res = check_background_symbolic(MUT)
    for k, e in res.items():
        line("S1 background (sympy, closed form)", e == 0, f"{k}: residual {e}")
    worst, herr = check_background_numeric(MUT)
    line("S2 background (mpmath 40 digits, 5 members, 7 points)", worst < 1e-30, f"max residual {float(worst):.3e} (declared <= 1e-30)")
    line("S3 horizon A'(h) = 2 (rho_0 - h)/R^2 ([BF06] (7); (8) as printed has no 2)", herr < 1e-30, f"max abs difference {float(herr):.3e}")
    print(f"      elapsed {time.time() - t0:.0f} s", flush=True)
    res, D = derive_master(MUT)
    for k, e in res.items():
        line("S4 linearized equations (sympy, generic A, R, phi0)", e == 0, f"{k}: {e if e != 0 else 0}")
    w = family_vs_generic(D, MUT)
    line("S5 family coefficients vs generic derived equation (mpmath 30 digits)", w < 1e-25, f"max abs difference {float(w):.3e}")
    print(f"      elapsed {time.time() - t0:.0f} s", flush=True)
    w = check_vhat_vs_bkz(MUT)
    line("S8 Vhat of u = Y/rho against [BKZ12] (40), (41) (mpmath 30 digits, members c = 0, -0.3, -1.5)", w < 1e-25, f"max abs difference {float(w):.3e} (declared <= 1e-25)")
    lim, resid = check_pt_limit(MUT)
    line("S7 Poeschl-Teller limit c -> -1 (sympy)", resid == 0, f"limit of Vhat = {lim}; residual of u = sech(A1 x/2), W = A1/2: {resid}")
    err = check_pt_numeric(MUT)
    line("S7b Vhat of the member c = -1.0001 against the Poeschl-Teller well", err < 3e-3, f"max |Vhat - V_PT|/(A1^2/2) over the well = {err:.3e} (declared <= 3e-3)")
    lim, resid = check_pt_limit_plus(MUT)
    line("S7c c -> -1 from above (sympy)", resid == 0, f"limit of Vhat = {lim}; residual of u = 1/sinh(A1 x/2), W = A1/2: {resid}")
    err = check_pt_numeric_plus(MUT)
    line("S7d Vhat of the member c = -0.9999 against (A1^2/2)/sinh^2(A1 x/2)", err < 3e-3, f"max relative difference, -6 <= A1 x <= 3, |A1 x| >= 0.5 = {err:.3e} (declared <= 3e-3)")
    printed_vs_derived()
    return ok_all


def check_vhat_vs_bkz(mut=(), cs=(0.0, -0.3, -1.5), pts=("-1.3", "0.7", "2.4")):
    """Vhat of the Schroedinger form (u = Y/rho) against [BKZ12] (40), (41): V_eff = U + R_{r* r*}/R with
    U/A = eps phi'^2/R'^2 (R^2 V - 1) + 2 phi' R V_phi/R' + eps V_phiphi; d/dr* = A d/d rho. mpmath, 30 digits."""
    import mpmath as mp
    mp.mp.dps = 30
    worst = mp.mpf(0)
    sgn = -1 if "veff_sign" in mut else 1                          # control: the sign of the V_phi term, as in the printed (37)
    for c in cs:
        f = family_mp(c)
        A, V = f["A"], f["V"]
        for r0 in pts:
            r0 = mp.mpf(r0)
            R = lambda r: mp.sqrt(r * r + 1)
            phip = lambda r: mp.sqrt(2) / (r * r + 1)
            Vphi = lambda r: mp.diff(V, r) / phip(r)
            Rr, Rp, Rpp = R(r0), mp.diff(R, r0), mp.diff(R, r0, 2)
            Ar, Ap = A(r0), mp.diff(A, r0)
            Vv = V(r0)
            Vpp = mp.diff(Vphi, r0) / phip(r0)
            e = -1
            U = Ar * (e * phip(r0)**2 / Rp**2 * (Rr**2 * Vv - 1) + sgn * 2 * phip(r0) * Rr * Vphi(r0) / Rp + e * Vpp)
            Veff = U + (Ar * Ap * Rp + Ar**2 * Rpp) / Rr
            R2 = r0 * r0 + 1
            Vhat = 2 * Ar**2 / r0**2 - (1 + 1 / R2) * Ar * Ap / r0 + Ar**2 / R2**2 - 2 * Ar / R2
            worst = max(worst, abs(Veff - Vhat))
    return worst


def check_pt_limit(mut=()):
    """c -> -1: in the scaling rho = rho_h t, t fixed, rho_h -> 0, the potential Vhat of the Schroedinger form tends in x to the
    reflectionless Poeschl-Teller well with l(l+1) = 2, whose bound state has W = A1/2 = kappa. Checked in sympy."""
    import sympy as sp
    rh, t, A1, y, x = sp.symbols("rho_h t A1 y x", positive=True)
    rho = rh * t
    A = A1 * rh * (t - 1); Ap = A1                                  # leading orders of A and A' at the horizon
    R2 = 1 + rho**2
    two = 3 if "pt_ell" in mut else 2                                # control: a wrong coefficient of the wall term
    Vh = two * A**2 / rho**2 - (1 + 1 / R2) * A * Ap / rho + A**2 / R2**2 - 2 * A / R2
    lim = sp.simplify(sp.limit(Vh, rh, 0))
    V = lim.subs(t, 1 + sp.exp(A1 * x))                              # dx = d rho / A = dt / (A1 (t - 1))
    u = 1 / sp.cosh(A1 * x / 2)                                      # l = 1 ground state
    resid = sp.simplify(sp.expand((sp.diff(u, x, 2) - V * u - (A1 / 2)**2 * u).rewrite(sp.exp)))   # u_xx - V u = W^2 u, W = kappa
    return lim, resid


def check_pt_limit_plus(mut=()):
    """c -> -1 from above (rho_h = -s < 0): rho = s t, the wall 2A^2/rho^2 of the throat now sits at t = 0, one horizon-throat distance
    from the horizon at t = -1. In x, Vhat tends to (A1^2/2)/sinh^2(A1 x/2) and u = 1/sinh(A1 x/2), with a pole at the throat
    (Y = rho u regular), solves u_xx - Vhat u = kappa^2 u and decays at both ends. Checked in sympy."""
    import sympy as sp
    s_, t, A1, x = sp.symbols("s t A1 x", positive=True)
    rho = s_ * t
    A = A1 * s_ * (t + 1); Ap = A1
    R2 = 1 + rho**2
    two = 3 if "pt_ell" in mut else 2
    Vh = two * A**2 / rho**2 - (1 + 1 / R2) * A * Ap / rho + A**2 / R2**2 - 2 * A / R2
    lim = sp.simplify(sp.limit(Vh, s_, 0))
    V = lim.subs(t, sp.exp(A1 * x) - 1)                              # dx = d rho / A = dt / (A1 (t + 1))
    u = 1 / sp.sinh(A1 * x / 2)
    resid = sp.simplify(sp.expand((sp.diff(u, x, 2) - V * u - (A1 / 2)**2 * u).rewrite(sp.exp)))
    return lim, resid


def check_pt_numeric_plus(mut=()):
    """c = -0.9999: Vhat(rho) against (A1^2/2)/sinh^2(A1 x/2), relative difference over -6 <= y <= 3, |y| >= 0.5, y = A1 x
    (beyond y = 4 the O(rho) corrections to the limit exceed 0.3 %, at y = 6 they are 10 %)"""
    fam = Fam(-0.9999, "R2")
    s_ = -fam.rh
    err = 0.0
    for yv in list(np.linspace(-6, -0.5, 12)) + list(np.linspace(0.5, 3, 12)):
        r = s_ * (math.exp(yv) - 1)
        pt = fam.A1**2 / 2 / math.sinh(yv / 2)**2
        err = max(err, abs(fam.Vhat(r) - pt) / pt)
    return err


def check_pt_numeric(mut=()):
    """c = -1.0001: Vhat(rho) near the horizon against the Poeschl-Teller shape, max relative difference of the well"""
    fam = Fam(-1.0001, "R2")
    ys = np.linspace(-6, 6, 25)
    err = 0.0
    for yv in ys:
        r = fam.rh * (1 + math.exp(yv))
        pt = -2 * fam.A1**2 * math.exp(yv) / (1 + math.exp(yv))**2
        err = max(err, abs(fam.Vhat(r) - pt) / (fam.A1**2 / 2))
    return err


def printed_vs_derived():
    """[BKZ12] (36), (37) as printed against the derived equation (the printed forms are reported, not corrected)."""
    import mpmath as mp
    mp.mp.dps = 30
    c = mp.mpf("-0.3"); K = -2 * c / mp.pi
    A = lambda r: 1 + c * (r * r + 1) + K * (r + (r * r + 1) * mp.atan(r))
    V = lambda r: -c * (1 + 2 * r * r / (r * r + 1)) - K * (3 * r / (r * r + 1) + (1 + 2 * r * r / (r * r + 1)) * mp.atan(r))
    print("      [BKZ12] (37) as printed: U = e^{-2 mu2} ( eps (e^{2 mu3} - V) phi'^2/mu3'^2 - 2 phi' V_phi/mu3' + eps V_phiphi ),")
    print("      against (41) U/A = eps phi'^2/R'^2 (R^2 V - 1) + 2 phi' R V_phi/R' + eps V_phiphi, eps = -1, member c = -0.3:")
    print("        rho      U/A from (41)        U/A from (37) as printed     ratio")
    for r in ("-1.3", "0.7", "2.4"):
        r = mp.mpf(r)
        R = mp.sqrt(r * r + 1); Rp = r / R; phip = mp.sqrt(2) / R**2
        Vv = V(r); Vphi = mp.diff(V, r) / phip; Vpp = mp.diff(lambda s: mp.diff(V, s) / (mp.sqrt(2) / (s * s + 1)), r) / phip
        e = -1
        u41 = e * phip**2 / Rp**2 * (R**2 * Vv - 1) + 2 * phip * R * Vphi / Rp + e * Vpp
        mu3p = Rp / R
        u37 = e * (R**2 - Vv) * phip**2 / mu3p**2 - 2 * phip * Vphi / mu3p + e * Vpp
        print(f"        {mp.nstr(r, 4):>6}   {mp.nstr(u41, 14):>18}   {mp.nstr(u37, 14):>18}   {mp.nstr(u37 / u41, 8):>12}")
    print("      (36) as printed: e^{-2 mu2 - 2 nu} ddot(dphi) - dphi'' - dphi' (nu' - 2 mu3' + mu2') + U dphi = 0 gives, with (7), the")
    print("      coefficients 1 and -2 R'/R. The derived equation has 1/A^2 and A'/A + 2 R'/R (= nu' - mu2' + 2 mu3').")


# ------------------------------------------------------------------------------------------------ numerics
class Fam:
    """member c of the family, b = 1: A, A', horizon; eps = -1 phantom (eps = +1 is the control with the kinetic sign flipped)."""

    def __init__(self, c, res="R2"):
        self.c = float(c); self.K = -2 * self.c / math.pi
        self.eps = +1 if mutating("eps") else -1
        self.res = dict(RES[res])
        if mutating("rmax"):
            self.res.update(rmax0=4.0, rmaxW=0.0)
        if mutating("hz_start"):
            self.res.update(z0=0.05)                              # control: the horizon start too far from the horizon
        self.rh = self.A1 = None
        if self.c < 0:
            self.rh = self._horizon(); self.A1 = float(self.Ap(self.rh))

    def A(self, r):
        return 1 + self.c * (r * r + 1) + self.K * (r + (r * r + 1) * np.arctan(r))

    def Ap(self, r):
        return 2 * self.c * r + 2 * self.K * (1 + r * np.arctan(r))

    def _horizon(self):
        if self.c == -1.0:
            return 0.0
        lo = -1.0
        while self.A(lo) > 0:
            lo *= 2
        return brentq(self.A, lo, self.K + 1.0, xtol=1e-15, rtol=1e-15)

    def a_over_z(self, z):
        """A(rh + z)/z"""
        if abs(z) < 1e-3:
            rh, K, c, A1 = self.rh, self.K, self.c, self.A1
            s = 1 + rh * rh
            A2 = 2 * c + 2 * K * (np.arctan(rh) + rh / s)
            A3 = 4 * K / s**2
            A4 = -16 * K * rh / s**3
            A5 = -16 * K * (1 - 5 * rh * rh) / s**4
            return A1 + A2 * z / 2 + A3 * z**2 / 6 + A4 * z**3 / 24 + A5 * z**4 / 120
        return self.A(self.rh + z) / z

    def Vhat(self, r):
        """potential of the Schroedinger form u_xx = (V + W^2) u, u = Y/rho (the wall 2 A^2/rho^2 at the throat); eps = -1."""
        A = self.A(r); Ap = self.Ap(r); R2 = r * r + 1
        return 2 * A * A / r**2 - (1 + 1 / R2) * A * Ap / r + A * A / R2**2 - 2 * A / R2


def coef_x(fam, r, A, Ar, Ap):
    """Y_xx = Cx Y_x + Dx Y + W^2 Y with dx = d rho / A; Ar = A/rho"""
    R2 = r * r + 1
    e = fam.eps
    return -2 * e * Ar / R2 + 2 * A * r / R2, -e * A * A / R2**2 + e * Ar * Ap / R2 - 2 * A / R2


def rhs_x(fam, W):
    def f(x, s):
        r, Y, Yx = s
        A = fam.A(r); Ap = fam.Ap(r)
        Cx, Dx = coef_x(fam, r, A, A / r, Ap)
        return [A, Yx, Cx * Yx + Dx * Y + W * W * Y]
    return f


def rhs_zeta(fam, W):
    def f(x, s):
        zeta, Y, Yx = s
        z = np.exp(zeta); r = fam.rh + z
        at = fam.a_over_z(z); A = z * at
        Ar = at if fam.rh == 0.0 else A / r
        Cx, Dx = coef_x(fam, r, A, Ar, fam.Ap(r))
        return [at, Yx, Cx * Yx + Dx * Y + W * W * Y]
    return f


def rhs_arc(fam, W, d):
    """rho-form along rho = d exp(i theta): Y'' = P1 Y' + P0 Y + W^2 Y / A^2"""
    e = fam.eps
    def f(th, s):
        Y, Yr = s
        r = d * np.exp(1j * th)
        A = fam.A(r); Ap = fam.Ap(r); R2 = r * r + 1
        P1 = -2 * e / (r * R2) + 2 * r / R2 - Ap / A
        P0 = -e / R2**2 + e * Ap / (A * r * R2) - 2 / (A * R2)
        return [1j * r * Yr, 1j * r * (P1 * Yr + P0 * Y + W * W * Y / A**2)]
    return f


def run_chunks(f, s, direction, W, rtol, atol, stop_ev, node_ev=None):
    """integrate in x with renormalization of (Y, Y_x); returns state, accumulated log norm, nodes"""
    x = 0.0; ln = 0.0; nodes = 0
    evs = [stop_ev] + ([node_ev] if node_ev else [])
    while True:
        L = 6.0 / max(W, 1e-3)
        sol = solve_ivp(f, [x, x + direction * L], s, method="DOP853", rtol=rtol, atol=atol, events=evs)
        s = sol.y[:, -1].copy(); x = sol.t[-1]
        if node_ev:
            nodes += len(sol.t_events[1])
        n = np.hypot(s[1], s[2]); s[1:] /= n; ln += math.log(n)
        if sol.status == 1:
            return s, ln, nodes
        if abs(x) > 1e7:
            raise RuntimeError("integration did not reach the stop event")


def make_ev(fn, direction, terminal):
    fn.terminal = terminal; fn.direction = direction
    return fn


def match_point(fam):
    return 1.0 if (fam.rh is None or fam.rh < 0.5) else fam.rh + 1.0


def left_solution(fam, W, rm, nodes=False):
    """(rho, Y, Y_x, nodes, signflag) at rho = rm of the solution regular at the horizon (ingoing), continued through the throat."""
    p = fam.res; rtol = p["rtol"]
    node_ev = make_ev(lambda x, s: s[1], 0, False) if nodes else None
    n_nodes = 0; flag = None
    if mutating("throat_wall") and (fam.rh is None or fam.rh < 0):
        # control: independent sides, Dirichlet wall at the throat (Y ~ rho^3), right side only
        r0 = 1e-4
        s = np.array([r0, r0**3, fam.A(r0) * 3 * r0**2])
        s[1:] /= np.hypot(s[1], s[2])
        stop = make_ev(lambda x, ss: ss[0] - rm, 1, True)
        s, ln, nn = run_chunks(rhs_x(fam, W), s, +1, W, rtol, [1e-12, 1e-30, 1e-30], stop, node_ev)
        return s[0], s[1], s[2], nn, None
    if fam.rh is None:                                       # c = 0: second spatial infinity, decaying toward rho -> -inf
        r0 = -max(p["rmax0"], p["rmaxW"] / W)
        lg = W + fam.Vhat(r0) / (2 * W)
        s = np.array([r0, 1.0, 1.0 / r0 + lg]); s[1:] /= np.hypot(s[1], s[2])
        stop = make_ev(lambda x, ss: ss[0] + p["delta"], 1, True)
        s, ln, nn = run_chunks(rhs_x(fam, W), s, +1, W, rtol, [1e-12, 1e-30, 1e-30], stop, node_ev)
        n_nodes += nn
        rstart = s[0]
    else:
        if fam.rh == 0.0:
            s = np.array([math.log(p["z0"]), 1.0, fam.A1 + W])         # Y ~ rho^(1 + W/A1), the ingoing solution
        else:
            s = np.array([math.log(p["z0"]), 1.0, W])                  # Y ~ z^(W/A1)
        s[1:] /= np.hypot(s[1], s[2])
        rstop = rm if fam.rh >= 0 else -min(p["delta"], 0.5 * abs(fam.rh))
        stop = make_ev(lambda x, ss: fam.rh + np.exp(ss[0]) - rstop, 1, True)
        s, ln, nn = run_chunks(rhs_zeta(fam, W), s, +1, W, rtol, [1e-12, 1e-30, 1e-30], stop, node_ev)
        n_nodes += nn
        if fam.rh >= 0:
            return fam.rh + np.exp(s[0]), s[1], s[2], n_nodes, None
        s = np.array([fam.rh + np.exp(s[0]), s[1], s[2]])
        rstart = s[0]
    # throat: complex arc, upper half plane
    d = -rstart
    Yl, Yxl = s[1], s[2]
    sol = solve_ivp(lambda th, y: rhs_arc(fam, W, d)(th, y), [np.pi, 0.0], [Yl + 0j, Yxl / fam.A(rstart) + 0j], method="DOP853",
                    rtol=rtol, atol=1e-30)
    Yc, Yrc = sol.y[:, -1]
    flag = (np.sign(Yl) == np.sign(Yc.real), abs(Yc.imag) / abs(Yc))
    s = np.array([d, Yc.real, (Yrc * fam.A(d)).real]); s[1:] /= np.hypot(s[1], s[2])
    stop = make_ev(lambda x, ss: ss[0] - rm, 1, True)
    s, ln, nn = run_chunks(rhs_x(fam, W), s, +1, W, rtol, [1e-12, 1e-30, 1e-30], stop, node_ev)
    return s[0], s[1], s[2], n_nodes + nn, flag


def right_solution(fam, W, rm, nodes=False):
    p = fam.res; rtol = p["rtol"]
    r0 = max(p["rmax0"], p["rmaxW"] / W)
    lg = -W - fam.Vhat(r0) / (2 * W)
    s = np.array([r0, 1.0, fam.A(r0) / r0 + lg]); s[1:] /= np.hypot(s[1], s[2])
    node_ev = make_ev(lambda x, ss: ss[1], 0, False) if nodes else None
    stop = make_ev(lambda x, ss: ss[0] - rm, -1, True)
    s, ln, nn = run_chunks(rhs_x(fam, W), s, -1, W, rtol, [1e-12, 1e-30, 1e-30], stop, node_ev)
    return s[0], s[1], s[2], nn


def mismatch(fam, W, nodes=False):
    rm = match_point(fam)
    _, Yl, Ylx, nl, flag = left_solution(fam, W, rm, nodes)
    _, Yr, Yrx, nr = right_solution(fam, W, rm, nodes)
    D = Yl * Yrx - Ylx * Yr
    return (D, nl + nr, flag) if nodes else D


def count_member(c, res="R2", W_lo=W_LO, W_hi=W_HI, log=print):
    fam = Fam(c, res)
    if mutating("scan_range"):
        W_lo = 3.0                                          # control: a search range that excludes the growing modes
    p = fam.res
    Ws = np.exp(np.linspace(math.log(W_lo), math.log(W_hi), p["nW"]))
    t0 = time.time()
    D = np.array([mismatch(fam, W) for W in Ws])
    roots = []
    for i in range(len(Ws) - 1):
        if D[i] == 0 or np.sign(D[i]) != np.sign(D[i + 1]):
            roots.append(brentq(lambda W: mismatch(fam, W), Ws[i], Ws[i + 1], xtol=1e-14, rtol=1e-13))
    info = []
    for W in roots:
        Dn, nn, flag = mismatch(fam, W, nodes=True)
        info.append(dict(W=W, W2=W * W, W_m13=W * fam.K, nodes=int(nn), arc_same_sign=None if flag is None else bool(flag[0]),
                         arc_imag=None if flag is None else float(flag[1])))
    out = dict(c=c, res=res, rh=fam.rh, A1=fam.A1, kappa=None if fam.A1 is None else fam.A1 / 2, K=fam.K, n_roots=len(roots), roots=info,
               D_min_abs=float(np.min(np.abs(D))), D_first=float(D[0]), D_last=float(D[-1]),
               sign_changes=int(np.sum(np.sign(D[:-1]) != np.sign(D[1:]))), W_lo=W_lo, W_hi=W_hi, nW=p["nW"],
               seconds=time.time() - t0, mutate=sorted(MUT))
    return out


def planted_command(res_list=("R1", "R2")):
    ok_all = True
    ref = planted_value()
    import mpmath as mp
    mp.mp.dps = 30
    for res in res_list:
        o = count_member(0.0, res)
        n = o["n_roots"]
        d = abs(mp.mpf(o["roots"][0]["W2"]) - mp.mpf(ref)) if n else None
        ok = (n == 1) and d is not None and d <= PLANTED_TOL[res]
        ok_all &= ok
        print(("PASS  " if ok else "FAIL  ") + f"P1 planted Ellis (c = 0) {res}: growing modes {n} (declared 1); W^2 a^2 = "
              f"{o['roots'][0]['W2']:.15f}" if n else "FAIL  P1 planted: no growing mode found",
              f"; reference {ref}; |difference| = {float(d):.3e} (declared <= {PLANTED_TOL[res]:g}); nodes {o['roots'][0]['nodes'] if n else None}"
              if n else "", flush=True)
    return ok_all


# ------------------------------------------------------------------------------------------------ independent spectral check
SPEC = {"S-A": dict(N=401, X=45.0, tol_throat=1e-2, tol_nothroat=1e-2),
        "S-B": dict(N=801, X=90.0, tol_throat=2e-2, tol_nothroat=1e-4)}


def cheb(N):
    xc = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1.0) ** np.arange(N + 1)
    dX = xc[:, None] - xc[None, :]
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    return D - np.diag(D.sum(axis=1)), xc


def spectral_member(c, which):
    """eigenvalues lambda = W^2 of  u_xx - Vhat u = lambda u  (u = Y/rho, c <= -1: no throat) or of the Y-form
    Y_xx = Cx Y_x + Dx Y + lambda Y (c > -1), by Chebyshev collocation on x in [-X, X] with Dirichlet ends. Shares with the shooting
    code only the background (A, A') and the coefficients; nothing of the boundary conditions, throat crossing or matching."""
    fam = Fam(c, "R2")
    N, X = SPEC[which]["N"], SPEC[which]["X"]
    D, xc = cheb(N)
    x = X * xc
    if fam.rh is None:                                        # Ellis, A = 1, rho = x
        rho = x.copy(); A = np.ones_like(x); Ar = 1 / rho; Ap = np.zeros_like(x)
    else:
        rref = 0.0 if fam.rh < 0 else fam.rh + 1.0
        zref = rref - fam.rh
        neg = x[x < 0]; pos = x[x >= 0]
        z = np.empty_like(x)
        if len(neg):
            sol = solve_ivp(lambda xx, q: [fam.a_over_z(np.exp(q[0]))], [0, -X * 1.0001], [math.log(zref)], method="DOP853",
                            rtol=1e-13, atol=1e-13, t_eval=np.sort(neg)[::-1])
            z[x < 0] = np.exp(sol.y[0])[np.argsort(np.argsort(-neg))]
        sol = solve_ivp(lambda xx, q: [fam.A(fam.rh + np.exp(q[0])) / np.exp(q[0])], [0, X * 1.0001], [math.log(zref)], method="DOP853",
                        rtol=1e-13, atol=1e-13, t_eval=np.sort(pos))
        z[x >= 0] = np.exp(sol.y[0])[np.argsort(np.argsort(pos))]
        rho = fam.rh + z
        at = np.array([fam.a_over_z(v) for v in z])
        A = z * at
        Ar = at if fam.rh == 0.0 else A / rho
        Ap = fam.Ap(rho)
    D1 = D / X
    D2 = D1 @ D1
    sgn = -1.0 if mutating("spec_sign") else 1.0               # control: wrong sign of the potential term
    if fam.c <= -1.0:
        R2 = rho * rho + 1
        Vh = 2 * Ar**2 - (1 + 1 / R2) * Ar * Ap + A * A / R2**2 - 2 * A / R2
        Lop = D2 - sgn * np.diag(Vh)
        form = "u"
    else:
        Cx, Dx = coef_x(fam, rho, A, Ar, Ap)
        Lop = D2 - Cx[:, None] * D1 - sgn * np.diag(Dx)
        form = "Y"
    ev = np.linalg.eigvals(Lop[1:-1, 1:-1])
    good = ev[(ev.real > 1e-3) & (ev.real < 100)]
    return form, np.sort_complex(good), float(np.min(np.abs(rho)))


def spectral_command(d):
    import glob
    ok_all = True
    rows = {}
    for fn in glob.glob(d + "/member_*.json"):
        o = json.load(open(fn))
        if o["res"] == "R2" or o["c"] not in rows:
            rows[o["c"]] = o
    for c in sorted(rows):
        for which in ("S-A", "S-B"):
            form, ev, minrho = spectral_member(c, which)
            realpos = [e.real for e in ev if abs(e.imag) < 1e-6 * abs(e)]
            cplx = [e for e in ev if abs(e.imag) >= 1e-6 * abs(e)]
            ref = [r["W"] for r in rows[c]["roots"]]
            tol = SPEC[which]["tol_nothroat"] if c <= -1.0 else SPEC[which]["tol_throat"]
            ok = len(realpos) == len(ref) and not cplx and all(abs(math.sqrt(a) - b) / b <= tol for a, b in zip(sorted(realpos), sorted(ref)))
            ok_all &= ok
            rel = [abs(math.sqrt(a) - b) / b for a, b in zip(sorted(realpos), sorted(ref))]
            print(("PASS  " if ok else "FAIL  ") + f"S6 spectral {which} ({form}-form) c = {c:g}: positive real eigenvalues {len(realpos)} "
                  f"(shooting {len(ref)}), complex with Re > 1e-3: {len(cplx)}, W = {[round(math.sqrt(a), 9) for a in sorted(realpos)]}, "
                  f"|dW|/W = {['%.2e' % r for r in rel]} (declared <= {tol:g})", flush=True)
    return ok_all



def agree_resolutions(o1, o2):
    """R1 against R2: same number of growing modes, and every rate within ROOT_AGREE (relative)"""
    same = o1["n_roots"] == o2["n_roots"]
    rel = max([abs(a["W"] - b["W"]) / b["W"] for a, b in zip(o1["roots"], o2["roots"])] or [0.0])
    return same and rel <= ROOT_AGREE, rel


def bkz_verdict(c, n):
    """[BKZ12]: black universes are stable only at c = -1 and unstable at every other c (Sec. V.C, Conclusions 2)"""
    return (n == 0) == (c == -1.0)


def report_command(d):
    import glob
    rows = {}
    for fn in glob.glob(d + "/member_*.json"):
        o = json.load(open(fn))
        rows.setdefault(o["c"], {})[o["res"]] = o
    ok_all = True
    print("Black universes, radial modes, b = 1 (W in units 1/b; the bracket is W in BKZ12's units m = 1/3, i.e. W * 3m)")
    print("%8s %10s %9s | %-4s %-4s | %-30s %-5s %-8s | %-12s | %-9s %s" % (
        "c", "rho_h", "kappa", "nR1", "nR2", "W [W in m=1/3 units]", "nodes", "W/kappa", "R1 v R2", "BKZ12", "evidence in BKZ12"))
    one = []
    for c in sorted(rows):
        r1, r2 = rows[c].get("R1"), rows[c].get("R2")
        agree = "-"
        if r1 and r2:
            ok, rel = agree_resolutions(r1, r2)
            ok_all &= ok
            agree = f"{rel:.1e} {'ok' if ok else 'FAIL'}"
        ref = r2 or r1
        nr = ref["n_roots"]
        agrees = bkz_verdict(c, nr)
        ok_all &= agrees
        Ws = "; ".join(f"{x['W']:.10f} [{x['W_m13']:.6f}]" for x in ref["roots"]) or "none"
        nodes = ",".join(str(x["nodes"]) for x in ref["roots"]) or "-"
        wk = ",".join("%.5f" % (x["W"] / ref["kappa"]) for x in ref["roots"]) if ref["kappa"] else "-"
        if nr == 1:
            one.append((c, ref["roots"][0]["W"], ref["roots"][0]["W_m13"]))
        print("%8g %10.6f %9s | %-4s %-4s | %-30s %-5s %-8s | %-12s | %-9s %s" % (
            c, ref["rh"], "-" if ref["kappa"] is None else "%.6f" % ref["kappa"], r1["n_roots"] if r1 else "-",
            r2["n_roots"] if r2 else "-", Ws, nodes, wk, agree, "agrees" if agrees else "DISAGREES",
            BKZ_EVIDENCE.get(c, "Conclusion 2 only")))
    print("members with exactly one growing mode: %d of %d" % (len(one), len(rows)))
    return ok_all


CONTROLS = {
    # mutation: (kind, target); kind: check/planted run the command and need a FAIL line with the target;
    # agree = R1 with the mutation against an unmutated R2 must fail agree_resolutions; bkz = the count must fail bkz_verdict
    "bg_V": ("check", "S1 background"), "bg_A": ("check", "S1 background"), "bg_kappa": ("check", "S3 horizon"),
    "master_p": ("check", "S4 linearized"), "bkz41_sign": ("check", "S4 linearized"), "family_q": ("check", "S5 family"),
    "eps": ("planted", "P1 planted"), "throat_wall": ("planted", "P1 planted"), "rmax": ("planted", "P1 planted"),
    "hz_start": ("agree", "R1 v R2"), "scan_range": ("bkz", "[BKZ12] verdict"), "spec_sign": ("spectral", "S6 spectral"),
    "pt_ell": ("check", "S7 Poeschl"), "veff_sign": ("check", "S8 Vhat"),
}
CONTROL_MEMBER = -0.3


def run_count(c, res, mut=None, out=None):
    args = [sys.executable, __file__, "count", str(c), "--res=" + res] + (["--mutate=" + mut] if mut else []) + (["--out=" + out] if out else [])
    r = subprocess.run(args, capture_output=True, text=True)
    return json.loads(r.stdout[r.stdout.index("{"):])


def controls_command():
    ok_all = True
    base = None
    for m, (kind, target) in CONTROLS.items():
        if kind == "spectral":
            import tempfile
            tmp = tempfile.mkdtemp(prefix="nm_i212_ctl_")
            run_count(-1.5, "R1", None, out=tmp + "/member_-1.5_R1.json")
            r = subprocess.run([sys.executable, __file__, "spectral", tmp, "--mutate=" + m], capture_output=True, text=True)
            txt = r.stdout + r.stderr
            fired = r.returncode == 1 and any(target in l and l.startswith("FAIL") for l in txt.splitlines())
            detail = "; ".join(l[:140] for l in txt.splitlines() if l.startswith("FAIL"))[:260]
            ret = r.returncode
        elif kind in ("check", "planted"):
            args = [sys.executable, __file__, kind] + (["--res=R1"] if kind == "planted" else []) + ["--mutate=" + m]
            r = subprocess.run(args, capture_output=True, text=True)
            txt = r.stdout + r.stderr
            fired = r.returncode == 1 and any(target in l and l.startswith("FAIL") for l in txt.splitlines())
            detail = "; ".join(l[:140] for l in txt.splitlines() if l.startswith("FAIL"))[:260]
            ret = r.returncode
        else:
            o = run_count(CONTROL_MEMBER, "R1", m)
            if kind == "agree":
                base = base or run_count(CONTROL_MEMBER, "R2")
                good, rel = agree_resolutions(o, base)
                fired = not good
                detail = f"c = {CONTROL_MEMBER}: n_roots {o['n_roots']} against {base['n_roots']}; max |dW|/W = {rel:.3e} (declared <= {ROOT_AGREE:g})"
            else:
                fired = not bkz_verdict(CONTROL_MEMBER, o["n_roots"])
                detail = f"c = {CONTROL_MEMBER}: n_roots = {o['n_roots']} (unmutated: 1; [BKZ12]: unstable)"
            ret = 0
        print(("fires   " if fired else "SILENT  ") + f"--mutate={m} (target {target}): exit {ret}; {detail}", flush=True)
        ok_all &= fired
    return ok_all


def main():
    global MUT
    argv = [a for a in sys.argv[1:] if not a.startswith("--mutate=")]
    for a in sys.argv[1:]:
        if a.startswith("--mutate="):
            MUT |= set(a.split("=", 1)[1].split(","))
    opts = {a.split("=")[0]: a.split("=", 1)[1] for a in argv if a.startswith("--") and "=" in a}
    pos = [a for a in argv if not a.startswith("--")]
    cmd = pos[0]
    if cmd == "check":
        sys.exit(0 if check_command() else 1)
    if cmd == "planted":
        sys.exit(0 if planted_command((opts["--res"],) if "--res" in opts else ("R1", "R2")) else 1)
    if cmd == "count":
        o = count_member(float(pos[1]), opts.get("--res", "R2"))
        s = json.dumps(o, indent=1)
        if "--out" in opts:
            open(opts["--out"], "w").write(s)
        print(s)
        sys.exit(0)
    if cmd == "report":
        sys.exit(0 if report_command(pos[1]) else 1)
    if cmd == "spectral":
        sys.exit(0 if spectral_command(pos[1]) else 1)
    if cmd == "controls":
        sys.exit(0 if controls_command() else 1)
    raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
