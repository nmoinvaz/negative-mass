#!/usr/bin/env python3
"""Linearization of the double-null equations with the potential about a static member, checked symbolically (#232).

    python3 linear_derive.py                  all checks; exit 1 if a check fails
    python3 linear_derive.py --mutate=NAME    break one step on purpose (the checks in MUTATIONS must fire)

Equations of #217 (evolve_derive.py): ds^2 = -2 e^{-f} dp dm + r^2 dOmega^2, eps = -1,

    r_pm   = -(r_p r_m + e^{-f}/2)/r + V r e^{-f}/2
    f_pm   = -2 (r_p r_m + e^{-f}/2)/r^2 + eps phi_p phi_m
    phi_pm = -(r_p phi_m + r_m phi_p)/r - (eps/2) e^{-f} V'(phi)
    C_pp   = r_pp + f_p r_p + (r/2) eps phi_p^2 = 0,   C_mm likewise.

Part A  the linearization U_pm = B0 U + Bp U_p + Bm U_m, U = (dr, df, dphi), against the derivative of the nonlinear right-hand sides (40 digits)
Part B  pure-gauge perturbations L_xi(member) solve the linearized evolution equations and constraints on the exact member (six members, 40 digits);
        Y = dr - (rho R/sqrt2) dphi is invariant under them
Part C  the master equation of #212 for u = Y/rho: u_pm + Vhat u/2 = M, M is a first-order gauge-invariant combination, static on solutions
Part D  the horizon: Raychaudhuri on the horizon, exponents of the Y equation there, the logarithm at c = -1, W = kappa
Part E  V'' of the potential
"""
import math
import os
import sys
import time

import mpmath as mp
import sympy as sp

MUT = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--mutate=")), None)
MUTATIONS = {
    # name: checks that must fire (computed in the development runs dev_70 to dev_74 of negative-mass-raw I232-mac, before they were declared)
    "lin_Vpp": ["B1", "C2"],        # V'' dropped from the potential's expansion (linearization and nonlinear right-hand side alike, so A1 stays consistent)
    "lin_cross": ["A1", "B1"],      # sign of the dphi_p phi_m cross terms of the f equation in the linearization only
    "gauge_sign": ["B1", "B2"],     # sign of xi' + eta' in df of the gauge perturbation
    "frob_log": ["D3"],             # the resonance coefficient of the Frobenius recursion evaluated at the wrong exponent
    "V2_coef": ["E1"],              # V'' = V_psipsi (the factor 1/2 dropped)
    "lin_const": ["A0"],            # a term independent of the perturbation added to the dphi row of the linearization
    "Y_coef": ["B3"],               # Y = dr - (rho R/2) dphi in the invariance check
    "M_ep": ["C1"],                 # a dphi_p term added to the right-hand side of the dphi equation in the master-equation check
    "hz_rp": ["D1"],                # r_p not set to zero in the horizon constraint
    "exp_W": ["D2"],                # the expected exponents of Y at the c = -1.5 horizon scaled by 1.1
    "exp_one": ["D3a"],             # the expected exponents at c = -1 without the 1
}
if MUT is not None and MUT not in MUTATIONS:
    sys.exit(f"unknown mutation {MUT}; known: {', '.join(MUTATIONS)}")

results = []


def check(cid, text, ok, detail=""):
    results.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


mp.mp.dps = 40
EPS = -1
s2 = sp.sqrt(2)

# ---------------------------------------------------------------------------------------------------------------------
# the member as functions of rho
# ---------------------------------------------------------------------------------------------------------------------
rho = sp.Symbol("rho", real=True)


def member(c):
    c = sp.nsimplify(c)
    K = -2 * c / sp.pi
    A = 1 + c * (rho ** 2 + 1) + K * (rho + (rho ** 2 + 1) * sp.atan(rho))
    R = sp.sqrt(rho ** 2 + 1)
    Ph = s2 * sp.atan(rho)
    psi = sp.atan(rho)
    V = -c * (2 - sp.cos(2 * psi)) - K * (sp.Rational(3, 2) * sp.sin(2 * psi) + psi * (2 - sp.cos(2 * psi)))
    V1 = sp.diff(V, rho) / sp.diff(Ph, rho)
    V2 = sp.diff(V1, rho) / sp.diff(Ph, rho)
    return dict(c=c, K=K, A=A, R=R, Ph=Ph, V=V, V1=V1, V2=V2)


# ---------------------------------------------------------------------------------------------------------------------
# the linearization with generic background jets (module level, silent: linear.py imports derived_matrices)
# ---------------------------------------------------------------------------------------------------------------------
import random

r, rp, rm, E, php, phm, V0, V1, V2 = sp.symbols("r rp rm E php phm V0 V1 V2")        # E = e^{-f}
a, ap, am, d, dp, dm, e, ep, em = sp.symbols("a ap am d dp dm e ep em")              # dr, df, dphi and their p and m derivatives
x_, s = sp.symbols("x s")                                                            # dphi at the point (for V(phi), V'(phi)), the linearization parameter
Vq = V0 + V1 * x_ + (0 if MUT == "lin_Vpp" else 1) * V2 * x_ ** 2 / 2
Vpq = V1 + (0 if MUT == "lin_Vpp" else 1) * V2 * x_
Fr = -(rp * rm + E / 2) / r + sp.Rational(1, 2) * Vq * r * E
Ff = -2 * (rp * rm + E / 2) / r ** 2 + EPS * php * phm
Fp = -(rp * phm + rm * php) / r - sp.Rational(EPS, 2) * E * Vpq
pert = {r: r + s * a, rp: rp + s * ap, rm: rm + s * am, E: E * sp.exp(-s * d), php: php + s * ep, phm: phm + s * em, x_: s * e}


def lin(F):
    return sp.expand(sp.diff(F.subs(pert, simultaneous=True), s).subs(s, 0))


Lin = [lin(Fr), lin(Ff), lin(Fp)]
if MUT == "lin_cross":
    Lin[1] = sp.expand(Lin[1].subs({em: -em, ep: -ep}, simultaneous=True))
if MUT == "lin_const":
    Lin[2] = sp.expand(Lin[2] + E * V0)
Uj = [a, d, e]
Upj = [ap, dp, ep]
Umj = [am, dm, em]
B0 = sp.Matrix(3, 3, lambda i, j: Lin[i].coeff(Uj[j]))
Bp = sp.Matrix(3, 3, lambda i, j: Lin[i].coeff(Upj[j]))
Bm = sp.Matrix(3, 3, lambda i, j: Lin[i].coeff(Umj[j]))
fn_nl = sp.lambdify((r, rp, rm, E, php, phm, V0, V1, V2, x_), [Fr, Ff, Fp], "mpmath")
fn_B = sp.lambdify((r, rp, rm, E, php, phm, V0, V1, V2), [B0, Bp, Bm], "mpmath")
fn_coef = sp.lambdify((r, rp, rm, E, php, phm, V0, V1, V2), [B0, Bp, Bm], "numpy")


def derived_matrices(rv, rpv, rmv, Ev, phpv, phmv, v0, v1, v2):
    """(B0, Bp, Bm) as 3 x 3 numpy arrays from the derivation, for linear.py's check of its coefficient code."""
    import numpy as np
    return [np.array(M, float) for M in fn_coef(rv, rpv, rmv, Ev, phpv, phmv, v0, v1, v2)]


def nl_rhs(st, dstate, scale):
    r0, rp0, rm0, E0, php0, phm0, v0, v1, v2 = st
    ar, arp, arm, ad, aphp, aphm, ae = dstate
    return fn_nl(r0 + scale * ar, rp0 + scale * arp, rm0 + scale * arm, E0 * mp.exp(-scale * ad), php0 + scale * aphp, phm0 + scale * aphm, v0, v1, v2, scale * ae)


def part_a():
    print("== Part A: the linearization ==")
    random.seed(232)
    for nm, M in (("B0", B0), ("Bp", Bp), ("Bm", Bm)):
        print(f"  {nm} =")
        for i in range(3):
            print("     [" + ",  ".join(str(sp.simplify(M[i, j])) for j in range(3)) + "]")
    rest = [sp.simplify(Lin[i] - sum(B0[i, j] * Uj[j] + Bp[i, j] * Upj[j] + Bm[i, j] * Umj[j] for j in range(3))) for i in range(3)]
    check("A0", "the linearization is a linear combination of (dr, df, dphi) and their first derivatives (no remainder)", all(t == 0 for t in rest))
    worst = mp.mpf(0)
    for _ in range(25):
        st = [mp.mpf(random.uniform(0.4, 3)), mp.mpf(random.uniform(-1, 1)), mp.mpf(random.uniform(-1, 1)), mp.mpf(random.uniform(0.1, 3)),
              mp.mpf(random.uniform(-1, 1)), mp.mpf(random.uniform(-1, 1)), mp.mpf(random.uniform(-2, 2)), mp.mpf(random.uniform(-2, 2)), mp.mpf(random.uniform(-2, 2))]
        dstate = [mp.mpf(random.uniform(-1, 1)) for _ in range(7)]
        ar, arp, arm, ad, aphp, aphm, ae = dstate
        h_ = mp.mpf("1e-15")
        fp_, fm_ = nl_rhs(st, dstate, h_), nl_rhs(st, dstate, -h_)
        fd = [(fp_[i] - fm_[i]) / (2 * h_) for i in range(3)]
        Bn = fn_B(*st)
        U = [ar, ad, ae]
        Up = [arp, mp.mpf(random.uniform(-1, 1)), aphp]      # df_p does not enter any right-hand side
        Um = [arm, mp.mpf(random.uniform(-1, 1)), aphm]
        lin_v = [sum(Bn[0][i, j] * U[j] + Bn[1][i, j] * Up[j] + Bn[2][i, j] * Um[j] for j in range(3)) for i in range(3)]
        worst = max(worst, max(abs(fd[i] - lin_v[i]) for i in range(3)))
    check("A1", "the linearized right-hand sides equal the derivative of the nonlinear ones (central difference, h = 1e-15, 40 digits, 25 random states)",
          worst < mp.mpf("1e-20"), f"max difference {mp.nstr(worst, 3)}")


# ---------------------------------------------------------------------------------------------------------------------
# Part B: pure gauge perturbations on the exact member
# ---------------------------------------------------------------------------------------------------------------------
NJ = 6
xi = sp.symbols(f"xi0:{NJ}")        # xi(p) and its derivatives
et = sp.symbols(f"eta0:{NJ}")       # eta(m) and its derivatives


def make_ops(A):
    def Dp(F):
        return sp.diff(F, rho) * A / s2 + sum(sp.diff(F, xi[k]) * xi[k + 1] for k in range(NJ - 1))

    def Dm(F):
        return -sp.diff(F, rho) * A / s2 + sum(sp.diff(F, et[k]) * et[k + 1] for k in range(NJ - 1))
    return Dp, Dm


def gauge_residuals(c):
    M = member(c)
    A, R, Ph = M["A"], M["R"], M["Ph"]
    Dp, Dm = make_ops(A)
    fp_ = -Dp(sp.log(A))
    fm_ = -Dm(sp.log(A))
    r_p, r_m, ph_p, ph_m = Dp(R), Dm(R), Dp(Ph), Dm(Ph)
    sg = -1 if MUT == "gauge_sign" else 1
    dr_ = xi[0] * r_p + et[0] * r_m
    dph_ = xi[0] * ph_p + et[0] * ph_m
    df_ = xi[0] * fp_ + et[0] * fm_ - sg * (xi[1] + et[1])
    back = (R, r_p, r_m, A, ph_p, ph_m, M["V"], M["V1"], M["V2"])
    Bn = [Mx.subs(dict(zip((r, rp, rm, E, php, phm, V0, V1, V2), back)), simultaneous=True) for Mx in (B0, Bp, Bm)]
    Uv = sp.Matrix([dr_, df_, dph_])
    Upv = Uv.applyfunc(Dp)
    Umv = Uv.applyfunc(Dm)
    Upm = Upv.applyfunc(Dm)
    res_ev = Upm - (Bn[0] * Uv + Bn[1] * Upv + Bn[2] * Umv)
    Cpp = Dp(Dp(dr_)) + r_p * Dp(df_) + fp_ * Dp(dr_) + sp.Rational(EPS, 2) * ph_p ** 2 * dr_ + EPS * R * ph_p * Dp(dph_)
    Cmm = Dm(Dm(dr_)) + r_m * Dm(df_) + fm_ * Dm(dr_) + sp.Rational(EPS, 2) * ph_m ** 2 * dr_ + EPS * R * ph_m * Dm(dph_)
    Y = dr_ - (rho * R / 2 if MUT == "Y_coef" else rho * R / s2) * dph_
    return list(res_ev) + [Cpp, Cmm], Y


def part_b():
    print("\n== Part B: gauge perturbations of the exact member ==")
    members = [-1.5, -1.01, -1.0, -0.99, -0.3, 0.0]
    rng = random.Random(7)
    worst_ev, worst_c, worst_Y = mp.mpf(0), mp.mpf(0), mp.mpf(0)
    t0 = time.time()
    for cv in members:
        exprs, Y = gauge_residuals(cv)
        syms = (rho,) + tuple(xi) + tuple(et)
        fn = sp.lambdify(syms, exprs + [Y], "mpmath")
        for _ in range(3):
            rv = mp.mpf(rng.choice([-0.6, -0.15, 0.35, 0.9, 1.7, 3.2])) + mp.mpf(rng.uniform(-0.05, 0.05))
            jets = [mp.mpf(rng.uniform(-1, 1)) for _ in range(2 * NJ)]
            vals = fn(rv, *jets)
            worst_ev = max(worst_ev, max(abs(v) for v in vals[:3]))
            worst_c = max(worst_c, max(abs(v) for v in vals[3:5]))
            worst_Y = max(worst_Y, abs(vals[5]))
    check("B1", "L_xi(member) solves the three linearized evolution equations: six members, three points, random jets of xi(p), eta(m) (40 digits)",
          worst_ev < mp.mpf("1e-25"), f"max residual {mp.nstr(worst_ev, 3)}, {time.time() - t0:.0f} s")
    check("B2", "and both linearized constraints", worst_c < mp.mpf("1e-25"), f"max residual {mp.nstr(worst_c, 3)}")
    check("B3", "Y = dr - (rho R/sqrt2) dphi vanishes on them (gauge invariance)", worst_Y < mp.mpf("1e-30"), f"max |Y| {mp.nstr(worst_Y, 3)}")


# ---------------------------------------------------------------------------------------------------------------------
# Part C: the master equation up to a static term
# ---------------------------------------------------------------------------------------------------------------------
def master_coeffs(c):
    """M = u_pm + Vhat u/2 on solutions of the evolution equations, u = Y/rho: a linear form in (a, e, ap, am, ep, em, d) with coefficients in rho."""
    Mm = member(c)
    A, R, Ph = Mm["A"], Mm["R"], Mm["Ph"]
    Dp, Dm = make_ops(A)
    r_p, r_m, ph_p, ph_m = Dp(R), Dm(R), Dp(Ph), Dm(Ph)
    back = (R, r_p, r_m, A, ph_p, ph_m, Mm["V"], Mm["V1"], Mm["V2"])
    Bn = [Mx.subs(dict(zip((r, rp, rm, E, php, phm, V0, V1, V2), back)), simultaneous=True) for Mx in (B0, Bp, Bm)]
    Uv = sp.Matrix([a, d, e])
    Upv = sp.Matrix([ap, dp, ep])
    Umv = sp.Matrix([am, dm, em])
    rhs = Bn[0] * Uv + Bn[1] * Upv + Bn[2] * Umv                      # (a_pm, d_pm, e_pm)
    if MUT == "M_ep":
        rhs = rhs + sp.Matrix([0, 0, ep])
    g = rho * R / s2
    Y = a - g * e
    Yp = ap - g * ep - Dp(g) * e
    Ym = am - g * em - Dm(g) * e
    Ypm = rhs[0] - g * rhs[2] - Dp(g) * em - Dm(g) * ep - Dp(Dm(g)) * e
    iq = 1 / rho
    u_pm = Ypm * iq + Yp * Dm(iq) + Ym * Dp(iq) + Y * Dp(Dm(iq))
    Vhat = 2 * A ** 2 / rho ** 2 - (1 + 1 / R ** 2) * A * sp.diff(A, rho) / rho + A ** 2 / R ** 4 - 2 * A / R ** 2
    Mx = sp.expand(u_pm + Vhat / 2 * Y / rho)
    return Mx, dict(A=A, R=R, Ph=Ph, Dp=Dp, Dm=Dm, back=back)


def part_c():
    print("\n== Part C: the master equation of #212 ==")
    worst_zero, worst_static = mp.mpf(0), mp.mpf(0)
    nonzero = []
    apm = sp.Symbol("apm")
    for cv in (-1.5, -1.0, -0.3, 0.0):
        Mx, info = master_coeffs(cv)
        Ac = info["A"]
        c_ep, c_em = Mx.coeff(ep), Mx.coeff(em)
        c_a, c_e, c_ap, c_am, c_d = (Mx.coeff(q) for q in (a, e, ap, am, d))
        f_all = sp.lambdify(rho, [c_ep, c_em, c_a, c_e, c_ap, c_am, c_d], "mpmath")
        Dp, Dm = info["Dp"], info["Dm"]
        R_, r_p, r_m, E_, ph_p, ph_m, V_, V1_, V2_ = info["back"]
        fp_ = -Dp(sp.log(Ac))
        fm_ = -Dm(sp.log(Ac))
        # the constraints solved for the second derivatives dr_pp, dr_mm
        app_s = -(r_p * dp + fp_ * ap + sp.Rational(EPS, 2) * ph_p ** 2 * a + EPS * R_ * ph_p * ep)
        amm_s = -(r_m * dm + fm_ * am + sp.Rational(EPS, 2) * ph_m ** 2 * a + EPS * R_ * ph_m * em)
        # (Dp + Dm) of M: the coefficients are functions of rho and (Dp + Dm) rho = 0; a -> ap + am, e -> ep + em, ap -> app + apm, am -> apm + amm, d -> dp + dm
        tot = sp.expand(c_a * (ap + am) + c_e * (ep + em) + c_ap * (app_s + apm) + c_am * (apm + amm_s) + c_d * (dp + dm))
        cl = [tot.coeff(q) for q in (ap, am, ep, em, dp, dm, a, apm)]
        fl = sp.lambdify(rho, cl, "mpmath")
        for rv in (mp.mpf("0.7"), mp.mpf("1.9"), mp.mpf("-0.4"), mp.mpf("2.6")):
            vals = f_all(rv)
            worst_zero = max(worst_zero, abs(vals[0]), abs(vals[1]))
            worst_static = max(worst_static, max(abs(v) for v in fl(rv)))
            nonzero.append(max(abs(v) for v in vals[2:]))
    check("C1", "M = u_pm + Vhat u/2 on solutions of the evolution equations contains dr, dphi, dr_p - dr_m and df only, no dphi_p, dphi_m: four members, four points",
          worst_zero < mp.mpf("1e-25"), f"max coefficient of dphi_p, dphi_m {mp.nstr(worst_zero, 3)}")
    check("C2", "(d/dt) M = 0 on solutions of the evolution equations and both constraints, so M is a function of x alone: the master equation holds up to a static term",
          worst_static < mp.mpf("1e-25"), f"max residual coefficient {mp.nstr(worst_static, 3)}; M itself is not zero (largest coefficient {mp.nstr(max(nonzero), 3)})")


# ---------------------------------------------------------------------------------------------------------------------
# Part D: the horizon
# ---------------------------------------------------------------------------------------------------------------------
W, w = sp.symbols("W w", positive=True)


def y_equation(c):
    Mm = member(c)
    A, R = Mm["A"], Mm["R"]
    Ap = sp.diff(A, rho)
    P = 2 / rho - Ap / A
    Q = 1 / R ** 4 - Ap / (A * rho * R ** 2) - 2 / (A * R ** 2) + W ** 2 / A ** 2          # Y'' = P Y' + Q Y (README.md Sec. 2)
    return A, Ap, P, Q


def part_d():
    print("\n== Part D: the horizon ==")
    # D1: the constraint C_pp on a null surface along which r_p = phi_p = 0 (the horizon): the equation for dr along it
    rq, rp_, php_, fp_s, a_, ap_, app_, dp_, ep_ = sp.symbols("rq rp_ php_ fp_s a_ ap_ app_ dp_ ep_")
    dCpp = app_ + dp_ * rp_ + fp_s * ap_ + sp.Rational(EPS, 2) * php_ ** 2 * a_ + EPS * rq * php_ * ep_
    hz = sp.simplify(dCpp.subs({php_: 0} if MUT == "hz_rp" else {rp_: 0, php_: 0}))
    check("D1", "on a null surface with r_p = phi_p = 0 the linearized constraint is dr_pp + f_p dr_p = 0, so dr = c1 + c2 e^{b p} when f_p = -b",
          sp.simplify(hz - (app_ + fp_s * ap_)) == 0)
    # D2: exponents of the Y equation at a horizon, c = -1.5
    A15, Ap15, P15, Q15 = y_equation(sp.Rational(-3, 2))
    rh15 = sp.Float(sp.nsolve(A15, rho, 0.34), 30)
    Pn = sp.lambdify(rho, (rho - rh15) * P15, "mpmath")
    Qn = sp.lambdify((rho, W), (rho - rh15) ** 2 * Q15, "mpmath")
    rhm = mp.mpf(str(rh15))
    eps_ = mp.mpf("1e-14")
    Ap_h = sp.lambdify(rho, Ap15, "mpmath")(rhm)
    W0 = mp.mpf("0.4941711675")
    p0h = Pn(rhm + eps_)
    q0h = Qn(rhm + eps_, W0)
    roots = mp.polyroots([1, -1 - p0h, -q0h])                  # s(s-1) - p0 s - q0 = 0
    want = sorted([-W0 / Ap_h * (mp.mpf('1.1') if MUT == 'exp_W' else 1), W0 / Ap_h * (mp.mpf('1.1') if MUT == 'exp_W' else 1)])
    got = sorted([mp.re(x) for x in roots])
    check("D2", "c = -1.5: the exponents of Y at the horizon are +-W/A'(rho_h), so Y ~ (rho - rho_h)^{W/A'} e^{Wt} = e^{Wv}, regular and nonzero on the horizon",
          max(abs(got[0] - want[0]), abs(got[1] - want[1])) < mp.mpf("1e-8"), f"{mp.nstr(got[0], 8)}, {mp.nstr(got[1], 8)} against {mp.nstr(want[0], 8)}, {mp.nstr(want[1], 8)}")
    # D3: c = -1
    c1 = -1
    K1 = -2 * sp.Integer(c1) / sp.pi
    A1 = 1 + c1 * (rho ** 2 + 1) + K1 * (rho + (rho ** 2 + 1) * sp.atan(rho))
    Ap1 = sp.diff(A1, rho)
    kap1 = sp.simplify(Ap1.subs(rho, 0) / 2)
    Wk = w * kap1
    pf = -2 + rho * Ap1 / A1                                                       # rho^2 Y'' + rho p Y' + q Y = 0
    qf = -rho ** 2 * (1 / (rho ** 2 + 1) ** 2 - Ap1 / (A1 * rho * (rho ** 2 + 1)) - 2 / (A1 * (rho ** 2 + 1)) + Wk ** 2 / A1 ** 2)
    ps_ = sp.series(pf, rho, 0, 3).removeO()
    qs_ = sp.series(qf, rho, 0, 3).removeO()
    p0f, p1f = (sp.simplify(ps_.coeff(rho, k)) for k in (0, 1))
    q0f, q1f = (sp.simplify(qs_.coeff(rho, k)) for k in (0, 1))
    sx = sp.Symbol("s")
    ind = sp.expand(sx * (sx - 1) + p0f * sx + q0f)
    ex1 = sp.solve(ind, sx)
    okex = sorted([sp.simplify(t) for t in ex1], key=str) == sorted([-w / 2, w / 2] if MUT == "exp_one" else [1 - w / 2, 1 + w / 2], key=str)
    check("D3a", "c = -1: the exponents of Y at the horizon (which is the throat) are 1 -+ W/(2 kappa), kappa = 2/pi (sympy, W = w kappa)", okex, f"kappa = {kap1}, exponents {ex1}")
    s_lo = (1 - w / 2) if MUT != "frob_log" else (1 + w / 2)
    res_c = sp.simplify((p1f * s_lo + q1f).subs(w, 1))
    check("D3", "c = -1: the logarithm of the exponent 1 - W/(2 kappa) branch has the coefficient p1 s + q1 = pi w (1 - w)/8, which vanishes at W = kappa: both branches are smooth there",
          res_c == 0 and sp.simplify(p1f * (1 - w / 2) + q1f - sp.pi * w * (1 - w) / 8) == 0, f"p1 s + q1 at w = 1: {res_c}")


# ---------------------------------------------------------------------------------------------------------------------
# Part E: V''
# ---------------------------------------------------------------------------------------------------------------------
def part_e():
    print("\n== Part E: the potential ==")
    worstE = mp.mpf(0)
    psi = sp.Symbol("psi")
    for cv in (-1.5, -1.0, -0.3, 0.0):
        Mm = member(cv)
        c_, K_ = Mm["c"], Mm["K"]
        Vpsipsi = -4 * c_ * sp.cos(2 * psi) + 2 * K_ * sp.sin(2 * psi) - 4 * K_ * psi * sp.cos(2 * psi)         # linear.py pot3
        Vpp_claim = Vpsipsi / (1 if MUT == "V2_coef" else 2)
        fn1 = sp.lambdify(rho, Mm["V2"], "mpmath")
        fn2 = sp.lambdify(psi, Vpp_claim, "mpmath")
        for rv in (mp.mpf("-0.8"), mp.mpf("0.3"), mp.mpf("1.4")):
            worstE = max(worstE, abs(fn1(rv) - fn2(mp.atan(rv))))
    check("E1", "V''(phi) of linear.py (V_psipsi/2, psi = arctan rho) equals the second derivative of V along the member: four members, three points, 40 digits",
          worstE < mp.mpf("1e-30"), f"max difference {mp.nstr(worstE, 3)}")


def main():
    t0 = time.time()
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    bad = [cid for cid, ok in results if not ok]
    print(f"\n{len(results)} checks, {len(results) - len(bad)} pass; failed: {bad}; mutation {MUT}; {time.time() - t0:.0f} s")
    if MUT is not None:
        want = MUTATIONS[MUT]
        print(f"mutation {MUT}: fired {bad}, declared {want}: {'equal' if sorted(bad) == sorted(want) else 'DIFFERENT'}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
