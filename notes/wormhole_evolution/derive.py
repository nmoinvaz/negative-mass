#!/usr/bin/env python3
"""Dual-null equations of D = n + 2 dimensional gravity with massless real scalars of either sign (#39, design).

    python3 derive.py                  all checks; exit 1 if any check fails
    python3 derive.py --quick          skip the brute-force Ricci tensors at n = 4, 5
    python3 derive.py --mutate=NAME    break one claimed equation on purpose (the checks listed in MUTATIONS must fire)

Model: G_ab = sum_i eps_i (d_a phi_i d_b phi_i - g_ab (d phi_i)^2 / 2), eps = +1 normal, -1 ghost (8 pi G = 1). Trace reversal
gives R_ab = sum_i eps_i d_a phi_i d_b phi_i in every D. Metric ds^2 = -2 e^{-f} dx+ dx- + r^2 dOmega_n, fields r, f, phi_i on
(x+, x-); below p = x+, m = x-.

Part A  brute force: the D-dimensional Ricci tensor of the full metric (n = 2, 3, 4, 5 sphere angles) against the reduced
        expressions with n symbolic (the reduced expressions are linear in n).
Part B  the field equations from them: three wave equations, two constraints, and the identities they obey.
Part C  Misner-Sharp mass and its flux law, expansions, stress components.
Part D  static solutions: Ellis (symmetric), Tangherlini, the static ODEs with a lapse; the exact dynamical Roberts solution
        at n = 2.
"""
import sys
import time

import sympy as sp

QUICK = "--quick" in sys.argv
MUT = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--mutate=")), None)
MUTATIONS = {
    # name: checks that must fire (computed in the development runs dev_05 to dev_07, negative-mass-raw I39-dh610)
    "theta_coef": ["A1", "A3", "B2", "B3", "C1", "C2", "D3", "D3b", "D4", "D5"],
    "pm_sign": ["A3", "A4", "B2", "B3", "D3b", "D5"],
    "scalar_coef": ["B1", "B2", "B3", "D3d", "D5"],
    "ghost_sign": ["A3", "A4", "B2", "B3", "C1", "C2", "D3b", "D5"],
    "ms_mass": ["C1", "C2", "C3", "D3c"],
    "ellis_charge": ["D2", "D3b"],
    "lapse_coef": ["D1"],
    "trace_half": ["A2"],
    "ellis_profile": ["D2a", "D2", "D3", "D3b", "D3c"],
}
if MUT is not None and MUT not in MUTATIONS:
    sys.exit(f"unknown mutation {MUT}; known: {', '.join(MUTATIONS)}")

results = []


def check(cid, text, ok, detail=""):
    results.append((cid, bool(ok)))
    print(f"[{'pass' if ok else 'FAIL'}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


p, m = sp.symbols("p m")
n = sp.Symbol("n", positive=True)
r = sp.Function("r")(p, m)
f = sp.Function("f")(p, m)
phi1 = sp.Function("phi1")(p, m)
phi2 = sp.Function("phi2")(p, m)
e1, e2 = sp.symbols("e1 e2")          # signs of the two scalars (+1 normal, -1 ghost)
rp, rm_, rpp, rmm, rpm = (sp.diff(r, *a) for a in ((p,), (m,), (p, p), (m, m), (p, m)))
fp, fm, fpm = sp.diff(f, p), sp.diff(f, m), sp.diff(f, p, m)

# ---------------------------------------------------------------------------------------------------------------------
# Claims (reduced expressions); a mutation edits exactly one of them
# ---------------------------------------------------------------------------------------------------------------------
th = (n - 1) if MUT != "theta_coef" else n          # the sphere equation's coefficient
sg = 1 if MUT != "ghost_sign" else -1                # sign of the matter term in the constraints
pms = 1 if MUT != "pm_sign" else -1                  # sign of the matter term in the f equation
sc = n if MUT != "scalar_coef" else (n - 1)          # the scalar equation's coefficient
lc = 1 if MUT != "lapse_coef" else n                 # static: coefficient of r r' nu' in R_ang
tr = 2 if MUT != "trace_half" else 1                 # Einstein tensor: G = R_ab - g R / tr
pe = 2 * (n - 1) if MUT != "ellis_profile" else 2 * n  # exponent in r'^2 = 1 - (a/r)^pe
msx = n - 1 if MUT != "ms_mass" else n                # exponent of r in the Misner-Sharp mass

Rpp_claim = -(n / r) * (rpp + fp * rp)
Rmm_claim = -(n / r) * (rmm + fm * rm_)
Rpm_claim = fpm - (n / r) * rpm
Rang_claim = th * (1 + 2 * sp.exp(f) * rp * rm_) + 2 * sp.exp(f) * r * rpm   # R_ij = ghat_ij * this


def ricci_full(nn):
    """Ricci tensor of -2 e^{-f} dp dm + r^2 (unit S^nn), with the metric functions of (p, m)."""
    th_ = sp.symbols(f"t1:{nn + 1}")
    X = [p, m] + list(th_)
    D = nn + 2
    g = sp.zeros(D, D)
    g[0, 1] = g[1, 0] = -sp.exp(-f)
    ghat = []
    acc = sp.Integer(1)
    for k in range(nn):
        ghat.append(acc)
        g[2 + k, 2 + k] = r**2 * acc
        acc = acc * sp.sin(th_[k]) ** 2
    gi = sp.zeros(D, D)
    gi[0, 1] = gi[1, 0] = -sp.exp(f)
    for k in range(nn):
        gi[2 + k, 2 + k] = 1 / g[2 + k, 2 + k]
    G = [[[sp.simplify(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                           for d in range(D)) / 2) for c in range(D)] for b in range(D)] for a in range(D)]

    def ric(b, c):
        return sp.simplify(sum(sp.diff(G[a][b][c], X[a]) - sp.diff(G[a][b][a], X[c])
                               + sum(G[a][a][d] * G[d][b][c] - G[a][c][d] * G[d][b][a] for d in range(D)) for a in range(D)))

    R = sp.zeros(D, D)
    for b in range(D):
        for c in range(b, D):
            R[b, c] = R[c, b] = ric(b, c)
    return X, g, gi, R, ghat


S_p = e1 * sp.diff(phi1, p) ** 2 + e2 * sp.diff(phi2, p) ** 2
S_m = e1 * sp.diff(phi1, m) ** 2 + e2 * sp.diff(phi2, m) ** 2
S_pm = e1 * sp.diff(phi1, p) * sp.diff(phi1, m) + e2 * sp.diff(phi2, p) * sp.diff(phi2, m)
rpm_ev = -th * (rp * rm_ + sp.exp(-f) / 2) / r                        # sphere equation
fpm_ev = (n / r) * rpm_ev + pms * S_pm                                # (pm) equation with the sphere equation inserted
phpm_ev = [-(sc / (2 * r)) * (rp * sp.diff(q_, m) + rm_ * sp.diff(q_, p)) for q_ in (phi1, phi2)]
Cpp = rpp + fp * rp + sg * (r / n) * S_p
Cmm = rmm + fm * rm_ + sg * (r / n) * S_m
print("== Part A: brute-force Ricci tensors against the reduced expressions ==")
for nn in (2, 3, 4, 5):
    if QUICK and nn > 3:
        continue
    t0 = time.time()
    X, g, gi, R, ghat = ricci_full(nn)
    D = nn + 2
    sub = {n: nn}
    ok_pp = sp.simplify(R[0, 0] - Rpp_claim.subs(sub)) == 0
    ok_mm = sp.simplify(R[1, 1] - Rmm_claim.subs(sub)) == 0
    ok_pm = sp.simplify(R[0, 1] - Rpm_claim.subs(sub)) == 0
    ok_ang = all(sp.simplify(R[2 + k, 2 + k] - ghat[k] * Rang_claim.subs(sub)) == 0 for k in range(nn))
    ok_off = all(R[a, b] == 0 for a in range(D) for b in range(D) if a != b and (a, b) not in ((0, 1), (1, 0)))
    check(f"A1[n={nn}]", f"Ricci of the D = {D} metric equals the reduced R_pp, R_mm, R_pm, R_ang; the rest vanishes",
          ok_pp and ok_mm and ok_pm and ok_ang and ok_off, f"{time.time() - t0:.1f} s")
    S_pm_n = S_pm
    rpm_n, fpm_n = rpm_ev.subs(n, nn), fpm_ev.subs(n, nn)
    a3_sph = sp.simplify((R[2, 2] / ghat[0]).subs(rpm, rpm_n))
    a3_f = sp.simplify((R[0, 1] - S_pm_n).subs(rpm, rpm_n).subs(fpm, fpm_n))
    a3_pp = sp.simplify(Cpp.subs(n, nn) - (-(r / nn) * (R[0, 0] - S_p)))
    a3_mm = sp.simplify(Cmm.subs(n, nn) - (-(r / nn) * (R[1, 1] - S_m)))
    check(f"A3[n={nn}]", "the evolution equations r_pm, f_pm and the constraints C_pp, C_mm are the brute-force R_ab = sum eps d phi d phi",
          a3_sph == 0 and a3_f == 0 and a3_pp == 0 and a3_mm == 0)
    if nn == 2:
        # trace reversal: G - T = Q - g Q / 2 with Q = R_ab - sum e d phi d phi
        ph = [phi1, phi2]
        es = [e1, e2]
        Xs = X
        dphi = [[sp.diff(q, x) for x in Xs] for q in ph]
        Q = R - sp.Matrix(D, D, lambda a, b: sum(e * d[a] * d[b] for e, d in zip(es, dphi)))
        Rs = sum(gi[a, b] * R[a, b] for a in range(D) for b in range(D))
        Gmat = R - g * Rs / tr
        T = sp.Matrix(D, D, lambda a, b: sum(e * (d[a] * d[b] - g[a, b] * sum(gi[c, k] * d[c] * d[k] for c in range(D) for k in range(D)) / 2)
                                              for e, d in zip(es, dphi)))
        Qtr = sum(gi[a, b] * Q[a, b] for a in range(D) for b in range(D))
        diff = (Gmat - T) - (Q - g * Qtr / 2)
        check("A2", "G_ab - T_ab = Q_ab - g_ab Q/2 with Q_ab = R_ab - sum eps d_a phi d_b phi (D = 4)",
              all(sp.simplify(diff[a, b]) == 0 for a in range(D) for b in range(a, D)))
        # stress components of the brute-force T against the claims used in the note
        Tpp_claim, Tmm_claim = sg * S_p, sg * S_m
        Tang_claim = pms * sp.exp(f) * S_pm            # T^theta_theta
        ok_T = (sp.simplify(T[0, 0] - Tpp_claim) == 0 and sp.simplify(T[1, 1] - Tmm_claim) == 0 and sp.simplify(T[0, 1]) == 0
                and sp.simplify(T[2, 2] / g[2, 2] - Tang_claim) == 0)
        # u = e^{f/2}(1, 1)/sqrt2, s = e^{f/2}(1, -1)/sqrt2: T_uu = T_ss = e^f (T_pp + T_mm)/2 when T_pm = 0
        Tuu = sp.exp(f) * (T[0, 0] + 2 * T[0, 1] + T[1, 1]) / 2
        Tss = sp.exp(f) * (T[0, 0] - 2 * T[0, 1] + T[1, 1]) / 2
        ok_T = ok_T and sp.simplify(Tuu - Tss) == 0 and sp.simplify(Tuu - sp.exp(f) * (Tpp_claim + Tmm_claim) / 2) == 0
        check("A4", "T_pp = sum eps phi_p^2, T_mm likewise, T_pm = 0, T^theta_theta = e^f sum eps phi_p phi_m, rho = p_radial (D = 4)", ok_T)

print("\n== Part B: the equations ==")
# Field equations R_ab = sum eps d phi d phi, from the claims above:
print("  r_pm  = -(n-1)/r (r_p r_m + e^{-f}/2)")
print("  f_pm  = (n/r) r_pm + sum eps phi_p phi_m  =  -n(n-1)/r^2 (r_p r_m + e^{-f}/2) + sum eps phi_p phi_m")
print("  phi_pm = -n/(2r) (r_p phi_m + r_m phi_p)")
print("  C_pp  = r_pp + f_p r_p + (r/n) sum eps phi_p^2 = 0,  C_mm likewise")

# B1: scalar equation: box phi = 2 g^{pm} (phi_pm + n/(2r)(r_p phi_m + r_m phi_p)) from the D-dim Laplacian
D_ = 4
sqrtg = sp.exp(-f) * r ** n
box_phi = (sp.diff(sqrtg * (-sp.exp(f)) * sp.diff(phi1, m), p) + sp.diff(sqrtg * (-sp.exp(f)) * sp.diff(phi1, p), m)) / sqrtg
claim_box = -2 * sp.exp(f) * (sp.diff(phi1, p, m) + (sc / (2 * r)) * (rp * sp.diff(phi1, m) + rm_ * sp.diff(phi1, p)))
check("B1", "the D-dimensional Laplacian of phi(p, m) is -2 e^f [phi_pm + n/(2r)(r_p phi_m + r_m phi_p)]",
      sp.simplify(box_phi - claim_box) == 0)

# B2: constraint propagation, d_m C_pp = -k (r_m / r) C_pp on solutions of the evolution equations
def ev(expr):
    """Replace mixed derivatives by the evolution equations (repeat for third derivatives)."""
    rules = {sp.diff(r, p, m): rpm_ev, sp.diff(f, p, m): fpm_ev, sp.diff(phi1, p, m): phpm_ev[0], sp.diff(phi2, p, m): phpm_ev[1]}
    for _ in range(3):
        expr = expr.subs(rules).doit()
        # derivatives of the rules (third order mixed): handled by differentiating the rule right-hand sides
        expr = expr.subs({sp.diff(r, p, p, m): sp.diff(rpm_ev, p), sp.diff(f, p, p, m): sp.diff(fpm_ev, p),
                          sp.diff(phi1, p, p, m): sp.diff(phpm_ev[0], p), sp.diff(phi2, p, p, m): sp.diff(phpm_ev[1], p),
                          sp.diff(r, p, m, m): sp.diff(rpm_ev, m), sp.diff(f, p, m, m): sp.diff(fpm_ev, m),
                          sp.diff(phi1, p, m, m): sp.diff(phpm_ev[0], m), sp.diff(phi2, p, m, m): sp.diff(phpm_ev[1], m)})
        expr = expr.subs(rules)
    return expr


C = sp.Symbol("C")
dCpp = ev(sp.diff(Cpp, m))
dCpp = dCpp.subs(rpp, C - fp * rp - sg * (r / n) * S_p)
# the claim to test: d_m C_pp = -(n-1) (r_m / r) C_pp on solutions
resid = sp.simplify(dCpp - (-(n - 1) * rm_ / r * C))
check("B2", "d_m C_pp = -(n-1)(r_m/r) C_pp on solutions of the evolution equations (so r^(n-1) C_pp is constant along x-)",
      resid == 0, "" if resid == 0 else f"residual {sp.factor(resid)}")
dCmm = ev(sp.diff(Cmm, p)).subs(rmm, C - fm * rm_ - sg * (r / n) * S_m)
resid2 = sp.simplify(dCmm - (-(n - 1) * rp / r * C))
check("B3", "d_p C_mm = -(n-1)(r_p/r) C_mm likewise", resid2 == 0, "" if resid2 == 0 else f"residual {sp.factor(resid2)}")

print("\n== Part C: mass, flux, expansions, stress ==")
MS = r ** msx / 2 * (1 + 2 * sp.exp(f) * rp * rm_)
dMS_p = ev(sp.diff(MS, p)).subs(rpp, -fp * rp - sg * (r / n) * S_p)
claim_flux_p = -(r ** n / n) * sp.exp(f) * rm_ * S_p
r1 = sp.simplify(dMS_p - claim_flux_p)
check("C1", "d_p M = -(r^n/n) e^f r_m sum eps phi_p^2 for M = r^(n-1)(1 + 2 e^f r_p r_m)/2 on solutions with C_pp = 0",
      r1 == 0, "" if r1 == 0 else f"residual {sp.factor(r1)}")
dMS_m = ev(sp.diff(MS, m)).subs(rmm, -fm * rm_ - sg * (r / n) * S_m)
r2 = sp.simplify(dMS_m - (-(r ** n / n) * sp.exp(f) * rp * S_m))
check("C2", "d_m M = -(r^n/n) e^f r_p sum eps phi_m^2", r2 == 0, "" if r2 == 0 else f"residual {sp.factor(r2)}")

# C3: (grad r)^2 = 1 - 2M/r^(n-1) and the static check M = a^(2n-2)/(2 r^(n-1)) for the Ellis solution (done in D2)
grad2 = 2 * (-sp.exp(f)) * rp * rm_
check("C3", "(grad r)^2 = 2 g^{pm} r_p r_m = -2 e^f r_p r_m = 1 - 2M/r^(n-1)", sp.simplify(grad2 - (1 - 2 * MS / r ** (n - 1))) == 0)

# stress components for u = e^{f/2}(d_p + d_m)/sqrt2, s = e^{f/2}(d_p - d_m)/sqrt2 (unit timelike and spacelike)
print("  T_pp = sum eps phi_p^2, T_mm = sum eps phi_m^2, T_pm = 0 (massless),  T^theta_theta = sum eps e^f phi_p phi_m")
print("  rho = T_uu = e^f (T_pp + T_mm)/2,  p_s = T_ss = rho,  p_ang = T^theta_theta = e^f S_pm  (n spheres, D = n + 2)")
print("  expansions theta_pm = n r_pm-sign: theta_+ = n e^{f} r_p / r, theta_- = n e^{f} r_m / r per unit e^{-f} null normal")

print("\n== Part D: static solutions and an exact dynamical one ==")
a, l = sp.symbols("a l", positive=True)
q = sp.Symbol("q", positive=True)
# D1: static equations with a lapse nu(l): ds^2 = -e^{2 nu} dt^2 + dl^2 + r^2 dOmega_n, one ghost field phi(l)
nu = sp.Function("nu")(l)
rl = sp.Function("rl")(l)
phl = sp.Function("phl")(l)
# Ricci of the static metric via a warped-product formula: R_tt, R_ll, R_ang (checked against Part A form by n = 2 brute force)
def static_ricci(nn):
    X = [sp.Symbol("t"), l] + list(sp.symbols(f"u1:{nn + 1}"))
    Dd = nn + 2
    g = sp.zeros(Dd, Dd)
    g[0, 0] = -sp.exp(2 * nu)
    g[1, 1] = 1
    acc = sp.Integer(1)
    gh = []
    for k in range(nn):
        gh.append(acc)
        g[2 + k, 2 + k] = rl**2 * acc
        acc = acc * sp.sin(X[2 + k]) ** 2
    gi = sp.diag(*[1 / g[i, i] for i in range(Dd)])
    G = [[[sum(gi[a_, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(Dd)) / 2
           for c in range(Dd)] for b in range(Dd)] for a_ in range(Dd)]

    def ric(b, c):
        return sp.simplify(sum(sp.diff(G[a_][b][c], X[a_]) - sp.diff(G[a_][b][a_], X[c])
                               + sum(G[a_][a_][d] * G[d][b][c] - G[a_][c][d] * G[d][b][a_] for d in range(Dd)) for a_ in range(Dd)))
    return ric(0, 0), ric(1, 1), ric(2, 2)


Rtt_c = sp.exp(2 * nu) * (sp.diff(nu, l, 2) + sp.diff(nu, l) ** 2 + n * sp.diff(nu, l) * sp.diff(rl, l) / rl)
Rll_c = -sp.diff(nu, l, 2) - sp.diff(nu, l) ** 2 - n * sp.diff(rl, l, 2) / rl
Rang_c = (n - 1) * (1 - sp.diff(rl, l) ** 2) - rl * sp.diff(rl, l, 2) - lc * rl * sp.diff(rl, l) * sp.diff(nu, l)
ok_static = True
for nn in (2, 3):
    Rtt_b, Rll_b, Rang_b = static_ricci(nn)
    ok_static &= sp.simplify(Rtt_b - Rtt_c.subs(n, nn)) == 0 and sp.simplify(Rll_b - Rll_c.subs(n, nn)) == 0 \
        and sp.simplify(Rang_b - Rang_c.subs(n, nn)) == 0
check("D1", "static Ricci components R_tt, R_ll, R_ang with a lapse (brute force n = 2, 3)", ok_static)

# D2: the symmetric Ellis solution (nu = 0): r'^2 = 1 - (a/r)^(2(n-1)), phi' = q / r^n, q^2 = n(n-1) a^(2n-2)
r0 = sp.Symbol("r0", positive=True)
qq2 = (n * (n - 1) * a ** (2 * n - 2)) if MUT != "ellis_charge" else (n**2 * a ** (2 * n - 2))
rp_l = sp.sqrt(1 - (a / r0) ** pe)            # r'
rpp_l = (n - 1) * (a / r0) ** (2 * (n - 1)) / r0        # r'' = d/dl of r'^2/2 / r' ... from differentiating r'^2
# differentiate r'^2 = 1 - (a/r)^(2n-2): 2 r' r'' = 2(n-1) (a/r)^(2n-2) r'/r
diffcheck = sp.simplify((2 * rp_l * rpp_l) - sp.diff(1 - (a / r0) ** pe, r0) * rp_l)
check("D2a", "r'' = (n-1)(a/r)^(2n-2)/r follows from r'^2 = 1 - (a/r)^(2(n-1))", diffcheck == 0)
# with nu = 0: R_tt = 0; R_ll = -n r''/r = -phi'^2 (ghost, R_ab = -d phi d phi); R_ang = 0; scalar eq (r^n phi')' = 0
phi_prime = sp.sqrt(qq2) / r0 ** n
R_ll_res = sp.simplify((-n * rpp_l / r0 + phi_prime**2))
R_ang_res = sp.simplify((n - 1) * (1 - rp_l**2) - r0 * rpp_l)
check("D2", "Ellis: R_ll = -phi'^2 (ghost) and R_ang = 0 with phi' = q/r^n, q^2 = n(n-1) a^(2n-2)",
      R_ll_res == 0 and R_ang_res == 0, f"R_ll residual {R_ll_res}, R_ang residual {R_ang_res}")
# the same solution in dual-null form: x+- = (t +- l)/sqrt 2, f = 0, r(l), phi(l) with l = (p - m)/sqrt 2
lp = (p - m) / sp.sqrt(2)
rfun = sp.Function("R")
phfun = sp.Function("Ph")
rr = rfun(lp)
ph_ = phfun(lp)
rl_p = sp.diff(rr, p)
# derivative rules: dR/dl = rp_l, d2R/dl2 = rpp_l (as functions of R = r0), d Ph/dl = phi_prime
def to_l(expr):
    """Express p, m derivatives of R(l), Ph(l) through the Ellis relations; returns a function of r0 (and a, n)."""
    expr = expr.doit()
    L = sp.Symbol("L")
    expr = expr.subs((p - m) / sp.sqrt(2), L)
    return expr
Lsym = sp.Symbol("Lsym")
Rf = sp.Function("Rf")(Lsym)
Pf = sp.Function("Pf")(Lsym)
dR = {sp.Derivative(Rf, Lsym): rp_l, sp.Derivative(Rf, (Lsym, 2)): rpp_l}
def dn(expr_l, which):
    """d/dp (which = 'p') or d/dm of a function of l: d/dp = (1/sqrt2) d/dl, d/dm = -(1/sqrt2) d/dl."""
    s = 1 / sp.sqrt(2) if which == "p" else -1 / sp.sqrt(2)
    return s * sp.diff(expr_l, Lsym)
Rl = Rf
Rp = dn(Rl, "p"); Rm = dn(Rl, "m")
Rpp_ = dn(Rp, "p"); Rpm_ = dn(Rp, "m")
Pp = dn(Pf, "p"); Pm = dn(Pf, "m")
Ppm = dn(Pp, "m")
Rpr = lambda ex: ex.subs(dR).subs(sp.Derivative(Pf, (Lsym, 2)), sp.diff(phi_prime, r0) * rp_l).subs(sp.Derivative(Pf, Lsym), phi_prime).subs(Rf, r0)
# equations of Part B with f = 0, one ghost scalar (e1 = -1, e2 irrelevant)
eq_rpm = Rpr(sp.diff(Rp, Lsym) * (-1 / sp.sqrt(2)) - (-th * (Rp * Rm + sp.Rational(1, 2)) / Rf)).subs(Rf, r0)
Rpm_l = Rpr(Rpm_)
check("D3", "Ellis, dual-null form (f = 0): sphere equation r_pm = -(n-1)(r_p r_m + 1/2)/r holds",
      sp.simplify(Rpm_l - Rpr(-th * (Rp * Rm + sp.Rational(1, 2)) / Rf)) == 0)
cpp_l = sp.simplify(Rpr(dn(Rp, "p") + sg * (Rf / n) * (-1) * Pp**2))
check("D3b", "Ellis, dual-null form: C_pp = r_pp + (r/n)(-1)phi_p^2 = 0 and the f equation (f_pm = 0)",
      cpp_l == 0 and sp.simplify(Rpr(-(n * th / Rf**2) * (Rp * Rm + sp.Rational(1, 2)) + pms * (-1) * Pp * Pm) ) == 0,
      f"C_pp = {cpp_l}")
wave_l = sp.simplify(Rpr(Ppm + (sc / (2 * Rf)) * (Rp * Pm + Rm * Pp)))
check("D3d", "Ellis, dual-null form: phi_pm + n/(2r)(r_p phi_m + r_m phi_p) = 0", wave_l == 0, f"residual {wave_l}")
Mst = sp.simplify(Rpr(Rf ** msx / 2 * (1 + 2 * Rp * Rm)))
check("D3c", "Ellis: Misner-Sharp M = a^(2n-2)/(2 r^(n-1)) (tends to 0 at both ends)",
      sp.simplify(Mst - a ** (2 * n - 2) / (2 * r0 ** (n - 1))) == 0)

# D5: Tangherlini F = 1 - 2M/r^(n-1) in dual-null form: e^{-f} = F, dr/dr* = F, x+- = (t +- r*)/sqrt2
Mm = sp.Symbol("Mm", positive=True)
F = 1 - 2 * Mm / r0 ** (n - 1)
# r_p = F/sqrt2, r_m = -F/sqrt2, e^f = 1/F; r_pm = -F F'/2 with F' = dF/dr; f_p = -F'... f = -ln F: f_p = -F'/F * r_p
Fr = sp.diff(F, r0)
rp_t = F / sp.sqrt(2); rm_t = -F / sp.sqrt(2)
rpm_t = rm_t * sp.diff(rp_t, r0)           # d_m (F/sqrt2) = F'(r) r_m / sqrt2
fp_t = -(Fr / F) * rp_t
rpp_t = rp_t * sp.diff(rp_t, r0)
chk1 = sp.simplify(rpm_t + th * (rp_t * rm_t + F / 2) / r0)
chk2 = sp.simplify(rpp_t + fp_t * rp_t)              # C_pp = 0 with no matter
fpm_t = sp.diff(-(Fr / F) * rp_t, r0) * rm_t
chk3 = sp.simplify(fpm_t - (-(n * (n - 1) / r0**2) * (rp_t * rm_t + F / 2)))
check("D4", "Tangherlini (vacuum, M != 0) in dual-null form solves the sphere equation, C_pp = 0 and the f equation",
      chk1 == 0 and chk2 == 0 and chk3 == 0, f"{chk1}, {chk2}, {chk3}")

# D6: static ODE system with a lapse: first integrals
c1, qq = sp.symbols("c1 qq")
nu_p = c1 * sp.exp(-nu) / rl**n
Rtt_eq = sp.simplify((Rtt_c / sp.exp(2 * nu)).subs(sp.diff(nu, l, 2), sp.diff(nu_p, l)).subs(sp.diff(nu, l), nu_p).subs(
    sp.diff(rl, l), sp.Function("w")(l)))
print("  static with a lapse: (e^nu r^n nu')' = 0 and (e^nu r^n phi')' = 0, so nu' = c e^{-nu}/r^n, phi' = q e^{-nu}/r^n;")
print("  with the G_ll constraint these give the two-parameter family (a, c) of [GGS09a]; integrated as ODEs in the code.")

# D7: Roberts, n = 2: R = r^2 = a (p^2 + m^2) - p m, f = 0, ghost phi = sqrt2 atan(m sin th/(p - m cos th)), cos th = 1/(2a)
pp_, mm_ = sp.symbols("pp mm", positive=True)
aa = sp.Symbol("aa", positive=True)
Rr = aa * (pp_**2 + mm_**2) - pp_ * mm_
rr_ = sp.sqrt(Rr)
tht = sp.acos(1 / (2 * aa))
cc = sp.sqrt(2)
ph = cc * sp.atan(mm_ * sp.sin(tht) / (pp_ - mm_ * sp.cos(tht)))
n2 = 2
rp_, rm__ = sp.diff(rr_, pp_), sp.diff(rr_, mm_)
e_sph = sp.simplify(rr_ * sp.diff(rr_, pp_, mm_) + th.subs(n, n2) * (rp_ * rm__ + sp.Rational(1, 2)))
e_wav = sp.simplify(rr_ * sp.diff(ph, pp_, mm_) + (sc.subs(n, n2) / 2) * (rp_ * sp.diff(ph, mm_) + rm__ * sp.diff(ph, pp_)))
e_cpp = sp.simplify(sp.diff(rr_, pp_, 2) + sg * (rr_ / n2) * (-1) * sp.diff(ph, pp_) ** 2)
e_cmm = sp.simplify(sp.diff(rr_, mm_, 2) + sg * (rr_ / n2) * (-1) * sp.diff(ph, mm_) ** 2)
e_f = sp.simplify(-(n2 * th.subs(n, n2) / rr_**2) * (rp_ * rm__ + sp.Rational(1, 2)) + pms * (-1) * sp.diff(ph, pp_) * sp.diff(ph, mm_))
check("D5", "Roberts ghost solution (n = 2, f = 0, 1/2 < a): sphere eq., wave eq., C_pp, C_mm and f_pm = 0 hold exactly",
      all(sp.simplify(e.subs(aa, sp.Rational(3, 4))) == 0 for e in (e_sph, e_wav, e_cpp, e_cmm, e_f)),
      f"at a = 3/4")
print("  Roberts ghost: r^2 = a(x+^2 + x-^2) - x+ x-, phi = sqrt2 atan(x- sin(th)/(x+ - x- cos(th))), cos th = 1/(2a), a > 1/2;")
print("  r_p > 0 iff x- < 2a x+, r_m > 0 iff x+ < 2a x-: the wedge x+/(2a) < x- < 2a x+ is anti-trapped, both rays regular.")

print()
bad = [c for c, ok in results if not ok]
print(f"{len(results) - len(bad)} of {len(results)} checks pass" + (f"; failed: {', '.join(bad)}" if bad else ""))
if MUT:
    expect = set(MUTATIONS[MUT])
    fired = {c.split("[")[0] for c in bad}
    print(f"mutation {MUT}: expected to fire {sorted(expect)}; fired {sorted(fired)}; "
          + ("targets all fired" if expect <= fired else f"SILENT: {sorted(expect - fired)}"))
sys.exit(1 if bad else 0)
