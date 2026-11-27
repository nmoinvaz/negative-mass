#!/usr/bin/env python3
"""Polar (even-parity) linearization about the phantom black universes / M-AdS wormholes / Ellis wormhole of
stability.py ([BF06] family, b = 1), derived in sympy.  See polar.md.

    polar_derive.py check [--mutate=M]   every symbolic and numeric check; exit status 1 when one fails
    polar_derive.py controls             runs every --mutate=M and asserts that its target check fires
    polar_derive.py emit                 Python source of the coefficient functions that polar.py carries

Family.  ds^2 = -A dt^2 + drho^2/A + R^2 dOmega^2, R^2 = rho^2 + 1, phi0' = sqrt2/R^2, eps = -1 (phantom).  The background equations
leave A free up to A'' = 2 (A - 1)/R^2 (two constants, c and m); V, V_phi, V_phiphi follow from A (bg_rules).  Equations:
E_ab = R_ab - eps d_a phi d_b phi - V g_ab = 0 and S = eps Box phi - V_phi = 0, perturbation ~ exp(W t) Y_l(theta).
Fields: h_tt = A H0 Y, h_t rho = H1 Y, h_rho rho = H2 Y/A, h_thth = R^2 K Y, h_phph = R^2 sin^2 K Y, delta phi = F Y (Regge-Wheeler gauge).

Gauge.  l >= 2: Regge-Wheeler (no h_tA, h_rhoA, no G): fixed completely and regular where phi0' != 0 and R != 0, so it has no singular point at
the throat (the delta R = 0 gauge of [BKZ12] does, at R' = 0).  l = 1: h_tA = h_rhoA = 0 and h_t rho = 0 (the gauge of [ABCKY23] Sec. 3.3), which leaves
the one pure-gauge solution of the checks below.  At a horizon the RW fields are singular like z^(W/A1 -/+ 1) (z = rho - rho_h): that is the gauge, not
the solution, and polar.py takes the ingoing solutions by Frobenius series.
"""
import math
import pickle
import sys
import time

import sympy as sp

t, rho, th, ph = sp.symbols("t rho theta varphi", real=True)
W, L = sp.symbols("W L")                       # L = l (l + 1)
eta = sp.Symbol("eta")
X = [t, rho, th, ph]
A = sp.Function("A")(rho)
R2 = rho**2 + 1
phi0 = sp.sqrt(2) * sp.atan(rho)
Yf = sp.Function("Y")(th)
Ys, Yps = sp.symbols("Ys Yps")
ex = sp.exp(W * t)
As, A1, lam = sp.symbols("As A1 lam")
MUT = set()
N = 4
HERE = __import__('os').path.dirname(__import__('os').path.abspath(__file__))


def mutating(name):
    return name in MUT


# ------------------------------------------------------------------------------------------------ background, linearization
def bg_rules():
    """V, V_phi, V_phiphi along the family as functions of rho, from A (the tt and scalar equations, checked in D1)"""
    Ap = sp.diff(A, rho)
    App = 2 * (A - 1) / R2
    Vs = -(A - 1 + rho * Ap) / R2
    if mutating("bg_V"):
        Vs = -(A - 1 + 2 * rho * Ap) / R2                    # control: wrong coefficient in V
    Vp = -sp.sqrt(2) * Ap / R2
    Vpp = -App + 2 * rho * Ap / R2
    return Vs, Vp, Vpp


def Asub(z, depth=4):
    App = 2 * (A - 1) / R2
    for _ in range(depth):
        z = z.subs(sp.diff(A, rho, 4), sp.diff(App, rho, 2)).subs(sp.diff(A, rho, 3), sp.diff(App, rho)).subs(sp.diff(A, rho, 2), App)
    return z


def linearize(h, dphi, eps=-1):
    """first-order variation of E_ab and S for the perturbation h_ab (4x4 matrix), delta phi; background residuals E0, S0"""
    g0 = sp.diag(-A, 1 / A, R2, R2 * sp.sin(th)**2)
    g0i = g0.inv()
    g = g0 + eta * h
    gi = g0i - eta * g0i * h * g0i
    Vs, Vp, Vpp = bg_rules()
    phi = phi0 + eta * dphi

    def Gam(gm, gim):
        return [[[sum(gim[a, d] * (sp.diff(gm[d, b], X[c]) + sp.diff(gm[d, c], X[b]) - sp.diff(gm[b, c], X[d])) for d in range(N)) / 2
                  for c in range(N)] for b in range(N)] for a in range(N)]
    Gfull = Gam(g, gi)
    G0 = [[[sp.simplify(Gfull[a][b][c].subs(eta, 0)) for c in range(N)] for b in range(N)] for a in range(N)]
    dG = [[[sp.expand(sp.diff(Gfull[a][b][c], eta).subs(eta, 0)) for c in range(N)] for b in range(N)] for a in range(N)]

    def ric(Gm):
        Rm = sp.zeros(N, N)
        for a in range(N):
            for b in range(N):
                Rm[a, b] = sum(sp.diff(Gm[c][a][b], X[c]) - sp.diff(Gm[c][a][c], X[b])
                               + sum(Gm[c][c][d] * Gm[d][a][b] - Gm[c][b][d] * Gm[d][a][c] for d in range(N)) for c in range(N))
        return Rm
    R0 = ric(G0).applyfunc(sp.simplify)
    dR = sp.zeros(N, N)
    for a in range(N):
        for b in range(N):
            dR[a, b] = sum(sp.diff(dG[c][a][b], X[c]) - sp.diff(dG[c][a][c], X[b])
                           + sum(dG[c][c][d] * G0[d][a][b] + G0[c][c][d] * dG[d][a][b]
                                 - dG[c][b][d] * G0[d][a][c] - G0[c][b][d] * dG[d][a][c] for d in range(N)) for c in range(N))
    dph = [sp.diff(phi, x) for x in X]
    dphi0 = [sp.diff(phi0, x) for x in X]
    sgn = -1 if mutating("lin_sign") else 1                  # control: wrong sign of the V_phi delta phi g_ab term
    dE = sp.zeros(N, N)
    E0 = sp.zeros(N, N)
    for a in range(N):
        for b in range(N):
            dpp = sp.diff(dph[a] * dph[b], eta).subs(eta, 0)
            dE[a, b] = dR[a, b] - eps * dpp - (sgn * Vp * dphi * g0[a, b] + Vs * h[a, b])
            E0[a, b] = R0[a, b] - eps * dphi0[a] * dphi0[b] - Vs * g0[a, b]
    sq = R2 * sp.sin(th) * (1 + eta * sum(g0i[a, b] * h[a, b] for a in range(N) for b in range(N)) / 2)
    box = sum(sp.diff(sq * sum(gi[a, b] * dph[b] for b in range(N)), X[a]) for a in range(N)) / sq
    S = eps * box - (Vp + eta * Vpp * dphi)
    return dE, sp.diff(S, eta).subs(eta, 0), E0, S.subs(eta, 0)


def angsub(z):
    """Y'' -> -cot Y' - L Y (Legendre), then Y, Y' -> symbols"""
    LL = L + 1 if mutating("ang_L") else L                    # control: wrong angular eigenvalue
    Ypp = -sp.cos(th) / sp.sin(th) * sp.diff(Yf, th) - LL * Yf
    for _ in range(3):
        z = z.subs(sp.diff(Yf, th, 4), sp.diff(Ypp, th, 2)).subs(sp.diff(Yf, th, 3), sp.diff(Ypp, th)).subs(sp.diff(Yf, th, 2), Ypp)
    return sp.expand(z.subs(sp.diff(Yf, th), Yps).subs(Yf, Ys))


def clean(z):
    return Asub(angsub(Asub(sp.expand(sp.simplify(z / ex)))))


def separate(dE, dS):
    """angular separation of the linearized equations: dict of ODE-level equations"""
    P = {k: clean(dE[k]) for k in [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2), (3, 3)]}
    P["S"] = clean(dS)
    Ev = {"tt": sp.simplify(P[(0, 0)].coeff(Ys)), "tr": sp.simplify(P[(0, 1)].coeff(Ys)), "rr": sp.simplify(P[(1, 1)].coeff(Ys)),
          "tth": sp.simplify(P[(0, 2)].coeff(Yps)), "rth": sp.simplify(P[(1, 2)].coeff(Yps)), "S": sp.simplify(P["S"].coeff(Ys))}
    Rthth, Rphph = P[(2, 2)], sp.expand(sp.simplify(P[(3, 3)] / sp.sin(th)**2))
    Ev["ang_tr"] = sp.simplify(sp.expand(Rthth + Rphph).coeff(Ys))
    Ev["ang_tl"] = sp.simplify(sp.expand(Rthth - Rphph).coeff(Yps) / (sp.cos(th) / sp.sin(th)))
    Ev["trace_has_Yp"] = sp.simplify(sp.expand(Rthth + Rphph).coeff(Yps))
    return Ev


class Jets:
    """jet symbols f_k for the perturbation functions and total derivative along rho (A'' = 2 (A - 1)/R^2)"""

    def __init__(self, names):
        self.names = names
        self.fn = {n: sp.Function(n)(rho) for n in names}
        self.jet = {n: [sp.Symbol(f"{n}_{k}") for k in range(5)] for n in names}

    def tojet(self, z, extra=None):
        z = z.subs(extra or {}).doit()
        for n in self.names:
            for k in (4, 3, 2, 1):
                z = z.subs(sp.diff(self.fn[n], rho, k), self.jet[n][k])
            z = z.subs(self.fn[n], self.jet[n][0])
        return z.subs(sp.diff(A, rho), A1).subs(A, As)

    def D(self, z):
        out = sp.diff(z, rho)
        for n in self.names:
            for k in range(4):
                out += sp.diff(z, self.jet[n][k]) * self.jet[n][k + 1]
        return out + sp.diff(z, As) * A1 + sp.diff(z, A1) * 2 * (As - 1) / R2


def fields_h(H0, H1, H2, K):
    h = sp.zeros(4, 4)
    h[0, 0] = A * H0 * Yf * ex
    h[0, 1] = h[1, 0] = H1 * Yf * ex
    h[1, 1] = H2 * Yf * ex / A
    h[2, 2] = R2 * K * Yf * ex
    h[3, 3] = R2 * sp.sin(th)**2 * K * Yf * ex
    return h


# ------------------------------------------------------------------------------------------------ l >= 2
def reduce_l2():
    """RW gauge, l >= 2.  Returns the first-order 4-system for (H0, H1, K, F), the algebraic relation, the jets, residual checks"""
    H0, H1, H2, K, F = [sp.Function(n)(rho) for n in ("H0", "H1", "H2", "K", "F")]
    dE, dS, E0, S0 = linearize(fields_h(H0, H1, H2, K), F * Yf * ex)
    Ev = separate(dE, dS)
    jt = Jets(["H0", "H1", "K", "F"])
    H0j = jt.fn["H0"]
    J = {k: sp.simplify(jt.tojet(v.subs(sp.Function("H2")(rho), H0j))) for k, v in Ev.items() if k not in ("ang_tl", "trace_has_Yp")}
    j = jt.jet
    sol = sp.solve([J["tth"], J["tr"], J["rth"]], [j["H1"][1], j["K"][1], j["H0"][1]], dict=True)
    assert len(sol) == 1
    sol = {k: sp.simplify(v) for k, v in sol[0].items()}

    def elim(z):
        z = z.subs({j["H1"][2]: jt.D(sol[j["H1"][1]]), j["K"][2]: jt.D(sol[j["K"][1]]), j["H0"][2]: jt.D(sol[j["H0"][1]])}).subs(sol)
        return sp.simplify(z)
    sec = {k: elim(J[k]) for k in ("tt", "rr", "ang_tr", "S")}
    Al = sp.numer(sp.together(sec["tt"]))
    F1 = sp.simplify(sp.solve(Al, j["F"][1])[0])
    if mutating("alg_coef"):
        F1 = F1 + j["F"][0] * rho                           # control: a wrong term in the algebraic relation
    sysd = {j["H0"][1]: sol[j["H0"][1]], j["H1"][1]: sol[j["H1"][1]], j["K"][1]: sol[j["K"][1]], j["F"][1]: F1}
    return dict(Ev=Ev, jt=jt, J=J, sol=sol, sec=sec, Al=Al, F1=F1, sysd=sysd, E0=E0, S0=S0)


def matrix_l2(r):
    """4x4 matrix for (H0, Pi = H1/W, K, F) as a function of lam = W^2"""
    j = r["jt"].jet
    vars4 = [j["H0"][0], j["H1"][0], j["K"][0], j["F"][0]]
    M = sp.Matrix([[sp.simplify(sp.expand(r["sysd"][j[n][1]]).coeff(v)) for v in vars4] for n in ("H0", "H1", "K", "F")])
    T = sp.diag(1, W, 1, 1)
    Mn = sp.simplify(T.inv() * M * T)
    return Mn.subs(W, sp.sqrt(lam))


# ------------------------------------------------------------------------------------------------ l = 1
def reduce_l1():
    """l = 1, h_tA = h_rhoA = h_t rho = 0 (so K = -H2): 3-system (H0, H2, F), algebraic relation, pure-gauge solution"""
    H0, H2, K, F = [sp.Function(n)(rho) for n in ("H0", "H2", "K", "F")]
    h = sp.zeros(4, 4)
    h[0, 0] = A * H0 * Yf * ex
    h[1, 1] = H2 * Yf * ex / A
    h[2, 2] = R2 * K * Yf * ex
    h[3, 3] = R2 * sp.sin(th)**2 * K * Yf * ex
    dE, dS, E0, S0 = linearize(h, F * Yf * ex)
    Ev = separate(dE, dS)
    jt = Jets(["H0", "H2", "F"])
    Kf = sp.Function("K")(rho)
    J = {k: sp.simplify(jt.tojet(v.subs(L, 2).subs(Kf, -H2))) for k, v in Ev.items() if k not in ("ang_tl", "trace_has_Yp")}
    j = jt.jet
    sol = sp.solve([J["tr"], J["rth"]], [j["H2"][1], j["H0"][1]], dict=True)
    assert len(sol) == 1
    sol = {k: sp.simplify(v) for k, v in sol[0].items()}

    def elim(z):
        z = z.subs({j["H0"][2]: jt.D(sol[j["H0"][1]]), j["H2"][2]: jt.D(sol[j["H2"][1]])}).subs(sol)
        return sp.simplify(z)
    sec = {k: elim(J[k]) for k in ("tt", "rr", "ang_tr", "S")}
    Al = sp.numer(sp.together(sec["tt"]))
    F1 = sp.simplify(sp.solve(Al, j["F"][1])[0])
    sysd = {j["H0"][1]: sol[j["H0"][1]], j["H2"][1]: sol[j["H2"][1]], j["F"][1]: F1}
    return dict(Ev=Ev, jt=jt, J=J, sol=sol, sec=sec, Al=Al, F1=F1, sysd=sysd)


def matrix_l1(r):
    j = r["jt"].jet
    vars3 = [j["H0"][0], j["H2"][0], j["F"][0]]
    M = sp.Matrix([[sp.simplify(sp.expand(r["sysd"][j[n][1]]).coeff(v)) for v in vars3] for n in ("H0", "H2", "F")])
    return M.subs(W, sp.sqrt(lam))


def gauge_l1():
    """the pure-gauge solution of the l = 1 system: xi^theta = c(rho) = sqrt(A)/R, xi^rho = -A R^2 c', xi^t = R^2 W c/A (h_t theta = h_rho theta = h_t rho = 0)"""
    c = sp.sqrt(A) / (R2 if mutating("l1_gauge") else sp.sqrt(R2))   # control: a wrong gauge function
    xr = -A * R2 * sp.diff(c, rho)
    xt = R2 * W * c / A
    htt = -2 * A * W * xt + xr * sp.diff(-A, rho)
    hrr = 2 * sp.diff(xr, rho) / A + xr * sp.diff(1 / A, rho)
    htr = -A * sp.diff(xt, rho) + W * xr / A

    def toS(z):
        z = z.doit()
        z = z.subs(sp.diff(A, rho, 3), sp.diff(2 * (A - 1) / R2, rho)).subs(sp.diff(A, rho, 2), 2 * (A - 1) / R2)
        return sp.simplify(z.subs(sp.diff(A, rho), A1).subs(A, As))
    g = sp.Matrix([toS(htt / A), toS(hrr * A), toS(xr * sp.sqrt(2) / R2)])
    return g, toS(htr)


# ------------------------------------------------------------------------------------------------ code emission
def emit(name, args, mat, cplx=False):
    from sympy.printing.pycode import pycode
    flat = [sp.simplify(e) for e in mat]
    reps, exprs = sp.cse(flat, optimizations="basic")
    lines = ["def %s(%s):" % (name, ", ".join(args))]
    for s, e in reps:
        lines.append("    %s = %s" % (s, pycode(e)))
    lines.append("    return (%s)" % ", ".join(("(%s)+0j" if cplx else "1.0*(%s)") % pycode(e) for e in exprs))
    return "\n".join(lines)


# ------------------------------------------------------------------------------------------------ checks
def run_checks():
    ok_all = True
    T0 = time.time()

    def line(tag, ok, msg):
        nonlocal ok_all
        ok_all &= ok
        print(("PASS  " if ok else "FAIL  ") + tag + ": " + msg, flush=True)

    try:
        # D1 background
        Vs, Vp, Vpp = bg_rules()
        dE0, _, E0, S0 = linearize(sp.zeros(4, 4), 0)
        resid = [sp.simplify(Asub(E0[a, b])) for a in range(4) for b in range(4)] + [sp.simplify(Asub(S0))]
        line("D1 background E_ab = 0, S = 0 on the family, generic A with A'' = 2(A - 1)/R^2 (sympy)", all(e == 0 for e in resid),
             "residuals " + ("all 0" if all(e == 0 for e in resid) else str([e for e in resid if e != 0][:2])))
        # V against stability.py's closed form at 40 digits
        import mpmath as mp
        sys.path.insert(0, HERE)
        import stability as st
        mp.mp.dps = 40
        worst = mp.mpf(0)
        for cv in (-3, -0.3, 0.5):
            fm = st.family_mp(cv)
            Vme = sp.lambdify((rho, A, sp.Symbol("Ap")), Vs.subs(sp.diff(A, rho), sp.Symbol("Ap")), "mpmath")
            for r in (mp.mpf(-2), mp.mpf("0.7"), mp.mpf(3)):
                worst = max(worst, abs(fm["V"](r) - Vme(r, fm["A"](r), mp.diff(fm["A"], r))))
        line("D1b V(rho) of the family against stability.py's family_mp (mpmath 40 digits, 3 members)", worst < 1e-30, f"max difference {float(worst):.2e} (declared <= 1e-30)")

        # D2 pure gauge
        xt, xr, xs = [sp.Function(n)(rho) for n in ("xt", "xr", "xs")]
        xi = [xt * Yf * ex, xr * Yf * ex, xs * sp.diff(Yf, th) * ex, 0]
        g0 = sp.diag(-A, 1 / A, R2, R2 * sp.sin(th)**2)
        h = sp.zeros(4, 4)
        for a in range(4):
            for b in range(4):
                h[a, b] = sum(xi[c] * sp.diff(g0[a, b], X[c]) + g0[c, b] * sp.diff(xi[c], X[a]) + g0[a, c] * sp.diff(xi[c], X[b]) for c in range(4))
        dE, dS, _, _ = linearize(h, sum(xi[c] * sp.diff(phi0, X[c]) for c in range(4)))
        res = [sp.simplify(Asub(angsub(Asub(sp.simplify(dE[k] / ex))))) for k in [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2), (3, 3)]] + [sp.simplify(Asub(angsub(Asub(sp.simplify(dS / ex)))))]
        line("D2 linearized E_ab, S on a pure-gauge perturbation (Lie derivative of g and phi along xi = (T Y, X Y, S Y'), sympy)", all(e == 0 for e in res),
             "8 residuals " + ("all 0" if all(e == 0 for e in res) else str([e for e in res if e != 0][:1])[:120]))
        print("      elapsed %.0f s" % (time.time() - T0), flush=True)

        # D3 l >= 2
        r = reduce_l2()
        Ev = r["Ev"]
        line("D3a traceless angular equation of RW gauge", sp.simplify(Ev["ang_tl"] - (sp.Function("H2")(rho) - sp.Function("H0")(rho))) == 0 and Ev["trace_has_Yp"] == 0,
             "ang_tl = H2 - H0 and the trace part has no Y'")
        j = r["jt"].jet
        prop = [sp.simplify(sp.cancel(r["sec"][k] / r["sec"]["tt"])) for k in ("rr", "ang_tr")]
        line("D3b tt, rr and trace-angular equations, after the three first-order ones, are one algebraic relation", all(not p.has(j["H0"][0]) for p in prop),
             "ratios rr/tt, ang_tr/tt = %s, %s (functions of rho, A only)" % (sp.factor(prop[0]), sp.factor(prop[1])))
        jt = r["jt"]
        sysd = r["sysd"]
        F2S = sp.solve(r["sec"]["S"], j["F"][2])[0]
        sub_sys = lambda z: z.subs({j["H1"][1]: sysd[j["H1"][1]], j["K"][1]: sysd[j["K"][1]], j["H0"][1]: sysd[j["H0"][1]]}).subs(j["F"][1], r["F1"])
        cons = sp.simplify(sub_sys(jt.D(r["F1"])) - sub_sys(F2S))
        line("D3c constraint propagation: d/drho of F' (from the algebraic relation) equals F'' of the scalar equation on the 4-system", cons == 0, f"difference {cons}")
        resids = []
        for k in ("tt", "rr", "ang_tr", "S"):
            z = r["J"][k].subs({j["H1"][2]: jt.D(sysd[j["H1"][1]]), j["K"][2]: jt.D(sysd[j["K"][1]]), j["H0"][2]: jt.D(sysd[j["H0"][1]]), j["F"][2]: jt.D(r["F1"])})
            for _ in range(3):
                z = z.subs({j["H1"][1]: sysd[j["H1"][1]], j["K"][1]: sysd[j["K"][1]], j["H0"][1]: sysd[j["H0"][1]], j["F"][1]: r["F1"]})
            resids.append(sp.simplify(z))
        line("D3d the second-order equations tt, rr, trace-angular, S on the 4-system (second derivatives from the system)", all(e == 0 for e in resids), "4 residuals " + ("all 0" if all(e == 0 for e in resids) else str([e for e in resids if e != 0])[:100]))
        M4 = matrix_l2(r)
        deg = max(sp.degree(sp.numer(sp.together(e)), lam) for e in M4 if e != 0)
        line("D3e the matrix for (H0, H1/W, K, F) is polynomial in W^2 = lam, no odd power of W", all(not e.has(sp.sqrt(lam)) for e in M4) and deg <= 1, f"degree in lam {deg}")
        print("      elapsed %.0f s" % (time.time() - T0), flush=True)

        # D4 Ellis limit against [ABCKY23] (KC0), (H0C0), (H1C0), (upC0): u = F/sqrt2, omega = i W
        Me = M4.subs({As: 1, A1: 0})
        H0s, Ps, Ks, Fs = sp.symbols("H0 Pi K F")
        u = Fs / sp.sqrt(2)
        v = sp.Matrix([H0s, Ps, Ks, Fs])
        mine = Me * v
        Wv = sp.sqrt(lam)
        H1v = Wv * Ps
        om = sp.I * Wv                                         # exp(W t) = exp(-i omega t)
        Lq = L
        sgn = -1 if mutating("ellis_sign") else 1                 # control: wrong sign of the scalar coupling in the comparison
        ab = {
            "K'": (sp.I * Lq * H1v - 2 * om * (rho * Ks - rho * H0s - 2 * sgn * u)) / (2 * om * R2),
            "H0'": ((-2 * sp.I * om**2 * R2 + sp.I * Lq) * H1v - 2 * om * (rho * Ks - rho * H0s + 2 * sgn * u)) / (2 * om * R2),
            "H1'": -sp.I * om * (H0s + Ks),
            "u'": -((2 * om**2 * R2 + (Lq + 0) - 2) * Ks + (Lq - 2) * H0s - 2 * sp.I * om * rho * H1v) / 4 - 2 * rho / R2 * u,
        }
        mp_ = {"K'": mine[2], "H0'": mine[0], "H1'": Wv * mine[1], "u'": mine[3] / sp.sqrt(2)}
        diffs = {k: sp.simplify(mp_[k] - ab[k]) for k in ab}
        line("D4a Ellis (c = 0) limit of the 4-system against [ABCKY23] (KC0), (H0C0), (H1C0) with u = F/sqrt2, omega = iW (sympy)", all(diffs[k] == 0 for k in ("K'", "H0'", "H1'")),
             "differences K': %s, H0': %s, H1': %s" % (diffs["K'"], diffs["H0'"], diffs["H1'"]))
        expect = Ks * (Lq - 2) / 2
        line("D4b (upC0) as printed differs from the derived u' by exactly K (l+2)(l-1)/2, i.e. the sign of its (l+2)(l-1) K term (reported, not corrected; D5 and D4a hold)",
             sp.simplify(diffs["u'"] - expect) == 0, "derived minus printed = %s" % sp.factor(diffs["u'"]))
        # D5: at c = 0 the scalar F decouples; its potential is that of [ABCKY23] (eq_Schr)
        row = Me[3, :]
        Fd = (row * v)[0]                                       # F'
        Fdd = sp.diff(Fd, rho) + sum(sp.diff(Fd, v[i]) * (Me * v)[i] for i in range(4))
        f = [sp.diff(Fd, v[i]) for i in range(4)]
        gcoef = [sp.simplify(sp.diff(Fdd, v[i])) for i in range(4)]
        a_ = sp.simplify(gcoef[0] / f[0])
        decoupled = all(sp.simplify(gcoef[i] - a_ * f[i]) == 0 for i in range(3))
        b_ = sp.simplify(gcoef[3] - a_ * f[3])
        line("D5a at c = 0 the scalar F decouples: F'' = a F' + b F, no H0, H1, K left", decoupled, "a = %s, b = %s" % (sp.factor(a_), sp.factor(b_)))
        Rr = sp.sqrt(R2)
        Z = sp.Function("Z")(rho)
        ode = sp.simplify(sp.diff(Z / Rr, rho, 2) - a_ * sp.diff(Z / Rr, rho) - b_ * Z / Rr)
        Zpp = sp.simplify(sp.solve(ode, sp.diff(Z, rho, 2))[0])
        V_ab = L / R2 - (2 if mutating("pot") else 3) / R2**2
        line("D5b Z = R F obeys Z'' = (V + lam) Z with V = L/R^2 - 3/R^4, the master equation of [ABCKY23] (eq_Schr)", not Zpp.has(sp.diff(Z, rho)) and sp.simplify(sp.diff(Zpp, Z) - (V_ab + lam)) == 0,
             "derived V = %s" % sp.factor(sp.simplify(sp.diff(Zpp, Z) - lam)))
        # D5c with (47) as printed, F does not decouple at c = 0: the (ell+2)(ell-1) K term of the F' equation has the opposite sign, and F'' keeps H1 and K terms
        printed_corr = 0 if mutating("print47") else sp.sqrt(2) * Ks * (Lq - 2) / 2        # control: the derived row in place of the printed one
        prow = sp.expand((Me * v)[3] - printed_corr)
        rws = [(Me * v)[0], (Me * v)[1], (Me * v)[2], prow]
        f2 = [sp.diff(prow, v[i]) for i in range(4)]
        g2 = [sp.simplify(sp.diff(sp.diff(prow, rho) + sum(sp.diff(prow, v[i]) * rws[i] for i in range(4)), v[j])) for j in range(4)]
        a2 = sp.simplify(g2[0] / f2[0])
        res2 = [sp.simplify(g2[i] - a2 * f2[i]) for i in range(3)]
        line("D5c with (47) as printed the scalar does not decouple at c = 0 (F'' - a F' keeps H1 and K terms)", any(e != 0 for e in res2),
             "residuals of the H0, H1, K coefficients: %s" % [sp.factor(e) for e in res2])
        print("      elapsed %.0f s" % (time.time() - T0), flush=True)

        # D6 l = 1
        r1 = reduce_l1()
        j1 = r1["jt"].jet
        sysd1 = r1["sysd"]
        jt1 = r1["jt"]
        F2S1 = sp.solve(r1["sec"]["S"], j1["F"][2])[0]
        sub1 = lambda z: z.subs({j1["H0"][1]: sysd1[j1["H0"][1]], j1["H2"][1]: sysd1[j1["H2"][1]]}).subs(j1["F"][1], r1["F1"])
        cons1 = sp.simplify(sub1(jt1.D(r1["F1"])) - sub1(F2S1))
        line("D6a l = 1 (h_t rho = 0 gauge): D(F') from the algebraic relation equals F'' of the scalar equation", cons1 == 0, f"difference {cons1}")
        resids = []
        for k in ("tt", "rr", "ang_tr", "S"):
            z = r1["J"][k].subs({j1["H0"][2]: jt1.D(sysd1[j1["H0"][1]]), j1["H2"][2]: jt1.D(sysd1[j1["H2"][1]]), j1["F"][2]: jt1.D(r1["F1"])})
            for _ in range(3):
                z = z.subs({j1["H0"][1]: sysd1[j1["H0"][1]], j1["H2"][1]: sysd1[j1["H2"][1]], j1["F"][1]: r1["F1"]})
            resids.append(sp.simplify(z))
        line("D6b l = 1: tt, rr, trace-angular, S on the 3-system", all(e == 0 for e in resids), "4 residuals " + ("all 0" if all(e == 0 for e in resids) else "nonzero"))
        M3 = matrix_l1(r1)
        g, htr = gauge_l1()
        Dv = lambda z: sp.diff(z, rho) + sp.diff(z, As) * A1 + sp.diff(z, A1) * 2 * (As - 1) / R2
        gres = (g.applyfunc(Dv) - M3.subs(lam, W**2) * g).applyfunc(sp.simplify)
        line("D6c l = 1: the pure-gauge solution xi^theta = sqrt(A)/R (with h_t rho = 0) solves the 3-system", htr == 0 and all(e == 0 for e in gres), f"h_t rho of the gauge field {htr}; residual {list(gres)}")
        # D6d at c = 0 the gauge-invariant scalar chi = F - (g_F/g_H2) H2 of l = 1 obeys chi'' + (2 rho/R^2) chi' - lam chi = 0, i.e. Z = R chi: Z'' = (1/R^4 + lam) Z
        Me1 = M3.subs({As: 1, A1: 0})
        ge = g.subs({As: 1, A1: 0})
        H0q, H2q, Fq, q1, q2 = sp.symbols("H0q H2q Fq q1 q2")
        vv = sp.Matrix([H0q, H2q, Fq])
        lamW = W**2
        Mw = Me1.subs(lam, lamW)
        chi = Fq - ge[2] / ge[1] * H2q
        psi = H0q - ge[0] / ge[1] * H2q
        Dd = lambda ex_: sp.diff(ex_, rho) + sum(sp.diff(ex_, vv[i]) * (Mw * vv)[i] for i in range(3))
        rep = {Fq: q1 + ge[2] / ge[1] * H2q, H0q: q2 + ge[0] / ge[1] * H2q}
        dchi = sp.simplify(Dd(chi).subs(rep))
        dpsi = sp.simplify(Dd(psi).subs(rep))
        N2 = sp.Matrix([[sp.diff(dchi, q1), sp.diff(dchi, q2)], [sp.diff(dpsi, q1), sp.diff(dpsi, q2)]]).applyfunc(sp.simplify)
        closed = not dchi.has(H2q) and not dpsi.has(H2q)
        a_, b_, c_, d_ = N2[0, 0], N2[0, 1], N2[1, 0], N2[1, 1]
        Rr = sp.sqrt(R2)
        chif = sp.Function("chi")(rho)
        psi_expr = (sp.diff(chif, rho) - a_ * chif) / b_
        ode = sp.simplify(sp.diff(psi_expr, rho) - c_ * chif - d_ * psi_expr)
        target = R2 * sp.diff(chif, rho, 2) + 2 * rho * sp.diff(chif, rho) - W**2 * R2 * chif
        ratio = sp.simplify(ode / target)
        line("D6d at c = 0 the gauge-invariant scalar of l = 1 obeys chi'' + (2 rho/R^2) chi' - W^2 chi = 0 (2-system closed, H2 drops out)", closed and not ratio.has(chif),
             "ode/(R^2 chi'' + 2 rho chi' - W^2 R^2 chi) = %s" % ratio)
        # D6e the potentials 1/R^4 (l = 1) and L/R^2 - 3/R^4 at L = 2 (l >= 2 master equation) are Darboux partners with w = -rho/R^2
        w_ = -rho / R2 * (-1 if mutating("susy") else 1)       # control: the wrong sign of the superpotential
        Vb = sp.simplify(w_**2 - sp.diff(w_, rho))
        Va = sp.simplify(w_**2 + sp.diff(w_, rho))
        line("D6e V(l = 1) = 1/R^4 = w^2 - w' and V(l >= 2 at L = 2) = 2/R^2 - 3/R^4 = w^2 + w', w = -rho/R^2: isospectral partners", sp.simplify(Vb - 1 / R2**2) == 0 and sp.simplify(Va - (2 / R2 - 3 / R2**2)) == 0,
             "w^2 - w' = %s, w^2 + w' = %s" % (sp.factor(Vb), sp.factor(Va)))
        # D6f the Darboux map Q = d/drho + w keeps the quasinormal-mode conditions: w is regular on the real line and tends to 0 at both ends, so Q maps
        # exp(-+W rho) to (-+W + w) exp(-+W rho) ~ -+W exp(-+W rho), nonzero for W != 0, and the decaying and the outgoing behavior are preserved
        wf = -rho / (sp.sqrt(R2) if mutating("susy_end") else R2)   # control: a superpotential that does not vanish at infinity
        poles = sp.solve(sp.denom(sp.together(sp.simplify(wf**2))), rho)
        regular = all(not z.is_real for z in poles)
        lim_p, lim_m = sp.limit(wf, rho, sp.oo), sp.limit(wf, rho, -sp.oo)
        line("D6f the Darboux superpotential w = -rho/R^2 is regular on the real line (poles at rho = +-i) and w -> 0 at both ends, so Q = d/drho + w keeps exp(-+W rho)", lim_p == 0 and lim_m == 0 and regular,
             "w(+inf) = %s, w(-inf) = %s, poles of w at %s" % (lim_p, lim_m, poles))
        print("      elapsed %.0f s" % (time.time() - T0), flush=True)

        # D7 emitted code in polar.py against the derived matrices
        sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
        import importlib
        pol = importlib.import_module("polar")
        import random
        random.seed(1)
        worst = 0.0
        for _ in range(6):
            rv, lv, Lv, Av, A1v = (random.uniform(-3, 3), random.uniform(0.05, 4), random.choice([2.0, 6.0, 12.0]), random.uniform(0.2, 2), random.uniform(-1.5, 1.5))
            ref4 = sp.Matrix(M4.subs({rho: rv, lam: lv, L: Lv, As: Av, A1: A1v})).evalf(20)
            got4 = pol.m4_entries(rv, lv, Lv, Av, A1v)
            ref3 = sp.Matrix(M3.subs({rho: rv, lam: lv, As: Av, A1: A1v})).evalf(20)
            got3 = pol.m3_entries(rv, lv, Av, A1v)
            refg = sp.Matrix(g.subs(W, sp.sqrt(lam)).subs({rho: rv, lam: lv, As: Av, A1: A1v})).evalf(20)
            gotg = pol.g3_entries(rv, lv, Av, A1v)
            for refm, got in ((ref4, got4), (ref3, got3), (refg, gotg)):
                for e1, e2 in zip(list(refm), got):
                    worst = max(worst, abs(float(e1) - float(e2)) / (1 + abs(float(e1))))
        if mutating("code"):
            worst += 1e-3                                       # control: a code table that differs from the derivation
        line("D7 coefficient code in polar.py (m4_entries, m3_entries, g3_entries) against the derived matrices at 6 random points", worst < 1e-12, f"max relative difference {worst:.2e} (declared <= 1e-12)")
        print("total %.0f s" % (time.time() - T0))
    except Exception as exc:                                 # a crash is a failure of the run, not a silent pass
        line("CRASH", False, repr(exc)[:200])
    return ok_all


CONTROLS = {"bg_V": "D1 background", "lin_sign": "D2 linearized", "ang_L": "D4a Ellis", "alg_coef": "D3c constraint", "ellis_sign": "D4a Ellis",
            "pot": "D5b Z = R F", "print47": "D5c with (47)", "l1_gauge": "D6c l = 1", "susy": "D6e V(l = 1)", "susy_end": "D6f the Darboux", "code": "D7 coefficient code"}


def controls_command():
    import subprocess
    ok_all = True
    for m, target in CONTROLS.items():
        r = subprocess.run([sys.executable, __file__, "check", "--mutate=" + m], capture_output=True, text=True)
        txt = r.stdout + r.stderr
        fails = [l for l in txt.splitlines() if l.startswith("FAIL")]
        fired = r.returncode == 1 and any(target in l for l in fails)
        ok_all &= fired
        print(("fires   " if fired else "SILENT  ") + f"--mutate={m} (target {target}): exit {r.returncode}; " + "; ".join(l[:110] for l in fails)[:300], flush=True)
    return ok_all


def emit_command():
    r = reduce_l2()
    M4 = matrix_l2(r)
    r1 = reduce_l1()
    M3 = matrix_l1(r1)
    g, _ = gauge_l1()
    gs = g.subs(W, sp.sqrt(lam))
    out = [emit("m4_entries", ["rho", "lam", "L", "As", "A1"], M4), emit("m3_entries", ["rho", "lam", "As", "A1"], M3), emit("g3_entries", ["rho", "lam", "As", "A1"], gs),
           emit("m4_entries_c", ["rho", "lam", "L", "As", "A1"], M4, True), emit("m3_entries_c", ["rho", "lam", "As", "A1"], M3, True),
           emit("g3_entries_c", ["rho", "lam", "As", "A1"], gs, True)]
    print("\n\n".join(out))


def main():
    global MUT
    argv = [a for a in sys.argv[1:] if not a.startswith("--mutate=")]
    for a in sys.argv[1:]:
        if a.startswith("--mutate="):
            MUT |= set(a.split("=", 1)[1].split(","))
    if not argv:
        raise SystemExit(__doc__)
    if argv[0] == "check":
        sys.exit(0 if run_checks() else 1)
    if argv[0] == "controls":
        sys.exit(0 if controls_command() else 1)
    if argv[0] == "emit":
        emit_command()
        return
    raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
