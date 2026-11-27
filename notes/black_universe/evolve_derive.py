#!/usr/bin/env python3
"""Double-null equations of D = 4 gravity with a phantom scalar and a potential V(phi) (#217, design).

    python3 evolve_derive.py                  all checks; exit 1 if a check fails
    python3 evolve_derive.py --mutate=NAME    break one claimed equation on purpose (the checks in MUTATIONS must fire)

Convention of #212: R_mn - eps d_m phi d_n phi - V(phi) g_mn = 0 and eps box phi = V'(phi), eps = -1 phantom, 8 pi G = 1, which is
G_mn = T_mn with L_matter = -eps (d phi)^2 / 2 - V. Metric of #39: ds^2 = -2 e^{-f} dx+ dx- + r^2 dOmega_2^2; below p = x+, m = x-.

Part A  brute force: the Ricci tensor of the full 4-dimensional metric and the Laplacian of phi, the field equations solved for the
        mixed derivatives, against the claimed forms; V = 0 against dn.py's right-hand side (dn.F_np, 40 random states)
Part B  constraint propagation with V (the potential must drop out of d_m C_pp)
Part C  the static member of the family of [BF06] in the tortoise null gauge: generic reduction to #212's equations, then the
        closed forms at 40 digits
"""
import math
import os
import sys
import time

import sympy as sp

MUT = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--mutate=")), None)
MUTATIONS = {
    # name: checks that must fire (computed in the development runs dev_01 to dev_05 of negative-mass-raw I217-mac, before they were declared here)
    "V_sign_r": ["A1", "B1", "C1", "C2"],       # sign of V in r_pm
    "V_coef_r": ["A1", "B1", "C1", "C2"],       # coefficient 1/2 of V r e^{-f} in r_pm replaced by 1
    "V_sign_f": ["A1", "B1", "C1", "C2"],       # sign of V e^{-f} in f_pm
    "Vp_sign": ["A1", "B1", "C1", "C2"],        # sign of the V' term in phi_pm
    "eps_sign": ["A1", "A2", "A3", "C1", "C2"], # kinetic sign eps -> -eps in the claims (B1 holds for either sign: it is a consistency check)
}
if MUT is not None and MUT not in MUTATIONS:
    sys.exit(f"unknown mutation {MUT}; known: {', '.join(MUTATIONS)}")

HERE = os.path.dirname(os.path.abspath(__file__))
results = []


def check(cid, text, ok, detail=""):
    results.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def zero2(e):
    """the expression vanishes for eps = +1 and for eps = -1"""
    return all(sp.simplify(e.subs(EPS, s)) == 0 for s in (1, -1))


p, m = sp.symbols("p m")
th, vp = sp.symbols("theta varphi")
r = sp.Function("r")(p, m)
f = sp.Function("f")(p, m)
ph = sp.Function("phi")(p, m)
EPS = sp.Symbol("epsilon")                       # +1 normal, -1 phantom; kept symbolic with EPS**2 = 1
a0, a1, a2 = sp.symbols("a0 a1 a2")              # V(phi) = a0 + a1 phi + a2 phi^2: generic V and V' at a point (the equations use V, V' only)
V = a0 + a1 * ph + a2 * ph ** 2
Vp = a1 + 2 * a2 * ph

rp, rm_, rpp, rmm, rpm = (sp.diff(r, *a) for a in ((p,), (m,), (p, p), (m, m), (p, m)))
fp, fm, fpm = sp.diff(f, p), sp.diff(f, m), sp.diff(f, p, m)
pp_, pm_, ppm = sp.diff(ph, p), sp.diff(ph, m), sp.diff(ph, p, m)

# ---------------------------------------------------------------------------------------------------------------------
# claims (the equations to be proved); a mutation edits exactly one
# ---------------------------------------------------------------------------------------------------------------------
sV_r = -1 if MUT == "V_sign_r" else 1
cV_r = sp.Integer(1) if MUT == "V_coef_r" else sp.Rational(1, 2)
sV_f = -1 if MUT == "V_sign_f" else 1
sVp = -1 if MUT == "Vp_sign" else 1
eps_c = -EPS if MUT == "eps_sign" else EPS

n = 2
rpm_claim = -(rp * rm_ + sp.exp(-f) / 2) / r + sV_r * cV_r * V * r * sp.exp(-f)
fpm_claim = (n / r) * rpm_claim + eps_c * pp_ * pm_ - sV_f * V * sp.exp(-f)
ppm_claim = -(rp * pm_ + rm_ * pp_) / r - sVp * (eps_c / 2) * sp.exp(-f) * Vp
Cpp = rpp + fp * rp + (r / n) * eps_c * pp_ ** 2
Cmm = rmm + fm * rm_ + (r / n) * eps_c * pm_ ** 2


def ricci_full():
    th_ = [th, vp]
    X = [p, m] + th_
    g = sp.zeros(4, 4)
    g[0, 1] = g[1, 0] = -sp.exp(-f)
    g[2, 2] = r ** 2
    g[3, 3] = r ** 2 * sp.sin(th) ** 2
    gi = sp.zeros(4, 4)
    gi[0, 1] = gi[1, 0] = -sp.exp(f)
    gi[2, 2] = 1 / g[2, 2]
    gi[3, 3] = 1 / g[3, 3]
    G = [[[sp.simplify(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(4)) / 2)
           for c in range(4)] for b in range(4)] for a in range(4)]

    def ric(b, c):
        return sp.simplify(sum(sp.diff(G[a][b][c], X[a]) - sp.diff(G[a][b][a], X[c])
                               + sum(G[a][a][d] * G[d][b][c] - G[a][c][d] * G[d][b][a] for d in range(4)) for a in range(4)))
    R = sp.zeros(4, 4)
    for b in range(4):
        for c in range(b, 4):
            R[b, c] = R[c, b] = ric(b, c)
    return X, g, gi, R


print("== Part A: brute force ==")
t0 = time.time()
X, g, gi, R = ricci_full()
ok_off = all(R[a, b] == 0 for a in range(4) for b in range(4) if a != b and (a, b) not in ((0, 1), (1, 0)))
print(f"  Ricci tensor of the D = 4 metric: {time.time() - t0:.1f} s, off-diagonal components other than (p, m) vanish: {ok_off}")

# the field equations E_ab = R_ab - eps d_a phi d_b phi - V g_ab = 0, written out
dphi = [sp.diff(ph, x) for x in X]
E = {(a, b): R[a, b] - EPS * dphi[a] * dphi[b] - V * g[a, b] for (a, b) in ((0, 0), (1, 1), (0, 1), (2, 2))}
sqrtg = sp.exp(-f) * r ** 2 * sp.sin(th)
box = sum(sp.diff(sqrtg * gi[a, a2_] * dphi[a2_], X[a]) for a in range(4) for a2_ in range(4) if gi[a, a2_] != 0) / sqrtg
ES = EPS * box - Vp

# solve for the mixed derivatives (the unknowns are the derivative objects)
S_rpm, S_fpm, S_ppm = sp.symbols("S_rpm S_fpm S_ppm")
sub_unk = {rpm: S_rpm, fpm: S_fpm, ppm: S_ppm}
rpm_sol = sp.simplify(sp.solve(E[(2, 2)].subs(sub_unk), S_rpm)[0])
fpm_sol = sp.simplify(sp.solve(E[(0, 1)].subs(sub_unk).subs(S_rpm, rpm_sol), S_fpm)[0])
ppm_sol = sp.simplify(sp.solve(sp.simplify(ES).subs(sub_unk), S_ppm)[0])
A1 = zero2(rpm_sol - rpm_claim) and zero2(fpm_sol - fpm_claim) and zero2(ppm_sol - ppm_claim)
print("  r_pm   = -(r_p r_m + e^{-f}/2)/r + V r e^{-f}/2")
print("  f_pm   = (2/r) r_pm + eps phi_p phi_m - V e^{-f}")
print("  phi_pm = -(r_p phi_m + r_m phi_p)/r - (eps/2) e^{-f} V'(phi)")
print("  (with r_pm inserted: f_pm = -2 (r_p r_m + e^{-f}/2)/r^2 + eps phi_p phi_m, the V terms cancel)")
check("A1", "the brute-force R_theta_theta, R_pm and eps box phi = V' solved for r_pm, f_pm, phi_pm equal the claimed forms (V general)", A1)

# constraints: g_pp = 0, so V does not enter
Cpp_bf = sp.simplify(-(r / n) * E[(0, 0)] + 0)
Cmm_bf = sp.simplify(-(r / n) * E[(1, 1)] + 0)
A2 = zero2(Cpp_bf - Cpp) and zero2(Cmm_bf - Cmm)
print("  C_pp = r_pp + f_p r_p + (r/2) eps phi_p^2 = 0,  C_mm likewise (no V: g_pp = g_mm = 0)")
check("A2", "-(r/n) E_pp and -(r/n) E_mm are C_pp and C_mm of #39, with no V term", A2)

# A3: V = 0 is dn.py's right-hand side (the kernel's F_np), at random states; also dn's constraint form
sys.path.insert(0, os.path.join(HERE, "..", "wormhole_evolution"))
import numpy as np   # noqa: E402
import dn            # noqa: E402

rng = np.random.default_rng(217)
sy = sp.symbols("r0 f0 rp0 rm0 pp0 pm0")
rr0, ff0, rpp0, rmm0, ppp0, pmm0 = sy
lam = sp.lambdify((rr0, ff0, rpp0, rmm0, ppp0, pmm0), [
    rpm_claim.subs({rp: rpp0, rm_: rmm0, pp_: ppp0, pm_: pmm0}).subs({r: rr0, f: ff0}).subs({a0: 0, a1: 0, a2: 0}).subs(EPS, -1),
    fpm_claim.subs({rp: rpp0, rm_: rmm0, pp_: ppp0, pm_: pmm0}).subs({r: rr0, f: ff0}).subs({a0: 0, a1: 0, a2: 0}).subs(EPS, -1),
    ppm_claim.subs({rp: rpp0, rm_: rmm0, pp_: ppp0, pm_: pmm0}).subs({r: rr0, f: ff0}).subs({a0: 0, a1: 0, a2: 0}).subs(EPS, -1)])
par = dn.kernel_params(2, [-1.0])
worst = 0.0
for _ in range(40):
    c = np.array([[rng.uniform(0.3, 3.0)], [rng.uniform(-1, 1)], [rng.uniform(-1, 1)]])
    up = rng.normal(size=(3, 1))
    um = rng.normal(size=(3, 1))
    gg = dn.F_np(c, up, um, par, np.array([-1.0]))[:, 0]
    mine = np.array(lam(c[0, 0], c[1, 0], up[0, 0], um[0, 0], up[2, 0], um[2, 0]), float)
    worst = max(worst, float(np.max(np.abs(gg - mine))))
check("A3", "with V = 0 the claimed right-hand sides equal dn.F_np (the #39 kernel, ghost, n = 2) at 40 random states", worst <= 1e-14,
      f"max difference {worst:.1e}")

# ---------------------------------------------------------------------------------------------------------------------
print("\n== Part B: constraint propagation ==")


def on_shell(expr):
    """replace mixed derivatives by the evolution equations (third-order ones by differentiating them)"""
    rules = {rpm: rpm_claim, fpm: fpm_claim, ppm: ppm_claim}
    d3 = {}
    for var, rhs in ((r, rpm_claim), (f, fpm_claim), (ph, ppm_claim)):
        d3[sp.diff(var, p, p, m)] = sp.diff(rhs, p)
        d3[sp.diff(var, p, m, m)] = sp.diff(rhs, m)
    for _ in range(3):
        expr = expr.subs(d3).doit()
        expr = expr.subs(rules).doit()
    return expr


Csym = sp.Symbol("C")
dC = on_shell(sp.diff(Cpp, m)).subs(rpp, Csym - fp * rp - (r / n) * eps_c * pp_ ** 2)
kap = dC.subs(Csym, 0)
lin = sp.diff(dC, Csym)
print(f"  d_m C_pp on solutions of the evolution equations = ({sp.simplify(lin.subs(EPS, -1))}) C_pp + ({sp.simplify(kap.subs(EPS, -1))}) at eps = -1")
check("B1", "d_m C_pp = -(r_m/r) C_pp exactly on solutions of the evolution equations with V general (the potential drops out; at V = 0 it is -(n-1) r_m/r)",
      zero2(kap) and zero2(lin + rm_ / r))

# ---------------------------------------------------------------------------------------------------------------------
print("\n== Part C: the static member, tortoise null gauge ==")
# x = (p - m)/sqrt2, f = -ln A, r = R(rho(x)), phi = Phi(rho(x)), d rho/dx = A, so e^{-f} = A and ds^2 = A(-dt^2 + dx^2) with t = (p + m)/sqrt2
rho = sp.Symbol("rho")
A = sp.Function("A")(rho)
Rf = sp.Function("R")(rho)
Ph = sp.Function("Phi")(rho)
Vr = sp.Function("Vr")(rho)                  # V along the solution; V'(phi) = Vr'/Phi'
SQ = sp.sqrt(2)


def dx(F):
    return A * sp.diff(F, rho)                # d/dx of a function of rho


def static_derivs(F):
    """(F_p, F_m, F_pp, F_mm, F_pm) of a function of rho(x), x = (p - m)/sqrt2"""
    Fx = dx(F)
    Fxx = dx(Fx)
    return Fx / SQ, -Fx / SQ, Fxx / 2, Fxx / 2, -Fxx / 2


r_p, r_m, r_pp, r_mm, r_pm = static_derivs(Rf)
f_all = static_derivs(-sp.log(A))
ph_all = static_derivs(Ph)
emf = A
Vpr = sp.diff(Vr, rho) / sp.diff(Ph, rho)
Eq_r = r_pm - (-(r_p * r_m + emf / 2) / Rf + sV_r * cV_r * Vr * Rf * emf)
Eq_f = f_all[4] - ((n / Rf) * (-(r_p * r_m + emf / 2) / Rf + sV_r * cV_r * Vr * Rf * emf) + eps_c * ph_all[0] * ph_all[1] - sV_f * Vr * emf)
Eq_p = ph_all[4] - (-(r_p * ph_all[1] + r_m * ph_all[0]) / Rf - sVp * (eps_c / 2) * emf * Vpr)
Eq_C = r_pp + f_all[0] * r_p + (Rf / n) * eps_c * ph_all[0] ** 2
# #212: (3) eps (A R^2 Phi')' - R^2 V_phi = 0, (4) (A' R^2)' + 2 R^2 V = 0, (5) 2 R''/R + eps Phi'^2 = 0, (6) A (R^2)'' - R^2 A'' - 2 = 0, and
# the angular equation 1 - A (R R')' - R^2 V - R A' R' = 0
e3 = EPS * sp.diff(A * Rf ** 2 * sp.diff(Ph, rho), rho) - Rf ** 2 * Vpr
e4 = sp.diff(sp.diff(A, rho) * Rf ** 2, rho) + 2 * Rf ** 2 * Vr
e5 = 2 * sp.diff(Rf, rho, 2) / Rf + EPS * sp.diff(Ph, rho) ** 2
e6 = A * sp.diff(Rf ** 2, rho, 2) - Rf ** 2 * sp.diff(A, rho, 2) - 2
eth = 1 - A * sp.diff(Rf * sp.diff(Rf, rho), rho) - Rf ** 2 * Vr - Rf * sp.diff(A, rho) * sp.diff(Rf, rho)
# the residuals as the multiples of #212's equations that the algebra gives (declared here, checked as identities for eps = +1 and -1)
mult_r = A / (2 * Rf)                                        # r_pm = (A/2R) [angular]
mult_C = A ** 2 * Rf / 4                                     # C_pp = (A^2 R/4) (5)
mult_p = -EPS * A / (2 * Rf ** 2)                            # phi_pm = -(eps A/2R^2) (3)
mult_f = (-A / (2 * Rf ** 2), A ** 2 / 2)                    # f_pm = -(A/2R^2) (6) + (A^2/2) (5)
ok_r = zero2(Eq_r - mult_r * eth)
ok_C = zero2(Eq_C - mult_C * e5)
ok_p = zero2(Eq_p - mult_p * e3)
ok_f = zero2(Eq_f - (mult_f[0] * e6 + mult_f[1] * e5))
print("  static residuals: r_pm = (A/2R) * [1 - A (R R')' - R^2 V - R A' R'],  C_pp = (A^2 R/4) * [2R''/R + eps Phi'^2],")
print("                    phi_pm = -(eps A/2R^2) * [eps (A R^2 Phi')' - R^2 V_phi],  f_pm = -(A/2R^2) * [A (R^2)'' - R^2 A'' - 2] + (A^2/2) * [2R''/R + eps Phi'^2]")
check("C1", "on the static ansatz (r = R(rho), e^{-f} = A, drho/dx = A) the r_pm, f_pm, phi_pm equations and C_pp are the multiples above of #212's "
      "angular equation, (6) and (5), (3) and (5)", ok_r and ok_C and ok_p and ok_f)

# C2: the closed forms of the family at 40 digits (b = 1, 3m = K = -2c/pi), residuals of the five equations
import mpmath as mp   # noqa: E402
mp.mp.dps = 40


def family(c):
    c = mp.mpf(c)
    K = -2 * c / mp.pi
    Af = lambda x: 1 + c * (x * x + 1) + K * (x + (x * x + 1) * mp.atan(x))
    Rr = lambda x: mp.sqrt(x * x + 1)
    Phi = lambda x: mp.sqrt(2) * mp.atan(x)
    Vv = lambda x: -c * (1 + 2 * x * x / (x * x + 1)) - K * (3 * x / (x * x + 1) + (1 + 2 * x * x / (x * x + 1)) * mp.atan(x))
    return Af, Rr, Phi, Vv


def static_residuals(c, x):
    Af, Rr, Phi, Vv = family(c)
    x = mp.mpf(x)
    ep = -1 if MUT != "eps_sign" else 1
    # derivatives by the chain rule d/dx = A d/d rho (exact), with rho-derivatives by mpmath at 40 digits
    D1 = lambda F: (lambda y: mp.diff(F, y))
    def Fx(F):
        return lambda y: Af(y) * mp.diff(F, y)
    def derivs(F):
        F1 = Fx(F)
        F2 = Fx(F1)
        return F1(x) / mp.sqrt(2), -F1(x) / mp.sqrt(2), F2(x) / 2, F2(x) / 2, -F2(x) / 2
    rP, rM, rPP, rMM, rPM = derivs(Rr)
    fP, fM, fPP, fMM, fPM = derivs(lambda y: -mp.log(Af(y)))
    pP, pM, pPP, pMM, pPM = derivs(Phi)
    r0, A0, V0 = Rr(x), Af(x), Vv(x)
    Vphi = mp.diff(Vv, x) / mp.diff(Phi, x)
    emf0 = A0
    cr = mp.mpf(1) / 2 if MUT != "V_coef_r" else mp.mpf(1)
    sr = sV_r
    rpm_c = -(rP * rM + emf0 / 2) / r0 + sr * cr * V0 * r0 * emf0
    fpm_c = (2 / r0) * rpm_c + ep * pP * pM - sV_f * V0 * emf0
    ppm_c = -(rP * pM + rM * pP) / r0 - sVp * (ep / 2) * emf0 * Vphi
    Cp = rPP + fP * rP + (r0 / 2) * ep * pP ** 2
    Cm = rMM + fM * rM + (r0 / 2) * ep * pM ** 2
    return [rPM - rpm_c, fPM - fpm_c, pPM - ppm_c, Cp, Cm]


worst = mp.mpf(0)
for c in (-1.5, -1.0, -0.3, 0.0):
    for x in ("-0.8", "-0.2", "0.35", "1.7", "4.0"):
        worst = max(worst, max(abs(e) for e in static_residuals(c, x)))
check("C2", "the family of [BF06] (b = 1, 3m = -2c/pi; c = -1.5, -1, -0.3, 0) satisfies the three evolution equations and both constraints, 40 digits",
      worst <= mp.mpf("1e-30"), f"max residual {mp.nstr(worst, 3)}")

bad = [cid for cid, ok in results if not ok]
print(f"\n{len(results)} checks, {len(results) - len(bad)} pass; failed: {bad}; mutation {MUT}")
if MUT is not None:
    want = MUTATIONS[MUT]
    print(f"mutation {MUT}: fired {bad}, declared {want}: {'equal' if sorted(bad) == sorted(want) else 'DIFFERENT'}")
sys.exit(1 if bad else 0)
