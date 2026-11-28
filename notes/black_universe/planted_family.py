#!/usr/bin/env python3
"""The exact time-dependent family of #269 as planted truth for #217's nonlinear double-null code (#277).

Member c (b = 1, K = -2c/pi), ingoing null chart  ds^2 = -F dv^2 + 2 dv drho + r^2 dOmega^2,  v = t + x,  r = lam(v) R(rho),  phi = sqrt2 arctan(rho),
F = A - 2 rho lam'/lam + 1/lam^2 - 1,  lam^2 = L + beta e^{K v}   (exact_mode.md).  The code (evolve.py) is double null, ds^2 = -2 e^{-f} dx+ dx- + r^2 dOmega^2.

Chart map.  x+ = (v - v0)/sqrt2 on the whole rectangle (x+ is constant on ingoing rays, as v is).  x- labels the outgoing null curves, drho/dv = F/2.  On the ingoing
ray v = v0 the label is the tortoise-like x- with drho/dx- = -F(v0, rho)/sqrt2, rho(0) = rho_c, i.e. e^{-f} = F(v0, rho) (the static member's gauge e^{-f} = A at beta = 0).
The Jacobian w = d rho/d x- (fixed v) obeys w_v = F' w/2 and e^{-f} = -sqrt2 w, so along an outgoing ray
    d rho/dv = F/2,      d ln(e^{-f})/dv = F'/2,      F' = dF/drho = A' - 2 mu,   mu = lam'/lam.
The data are r, f, phi on x- = 0 (the ray through rho_c, parametrized by x+) and on x+ = 0 (v = v0); the kernel is evolve.row_step, unchanged.  The gauge-invariant
Lambda^2 = r^2 cos^2(phi/sqrt2) = r^2/(1 + rho^2) equals lam(v)^2 at every x- on the family; v* = ln(-L/beta)/K is where it vanishes (beta < 0).

    planted_family.py check [--only=D1,D2,...] [--mutate=NAME] [--out=DIR]    all checks and the runs; exit 1 when one fails (a mutated run is meant to exit 1)
    planted_family.py case --c=C --beta=B --h=H [--Xp=X --Xm=X --mutate=NAME]  one run
    planted_family.py mutations DIR PREFIX      fired sets of the stored mutated logs DIR/PREFIX<name>.log against the declared ones
    planted_family.py rays                      reported diagnostic (not a check): the outgoing rays of the family near v* (beta < 0) and at their escape to rho = -infinity (beta > 0)
Status words as in AGENTS.md: every number is numerical evidence; D3's identities are checked at 50 digits.  No sentence here says a result is new.
"""
import json
import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "wormhole_evolution"))
MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
ONLY = next((a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--only=")), None)
OUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--out=")), None)
_argv = sys.argv
sys.argv = [sys.argv[0]]            # exact_mode.py and evolve.py read sys.argv at import
import evolve as ev                 # noqa: E402  (the code under test: kernel, params, run_grid)
import evolve_potential as ep       # noqa: E402
import dn                           # noqa: E402
sys.argv = _argv

SQ2 = math.sqrt(2.0)
KERNEL_MUTS = ("th_coef", "vterm_sign", "vterm_coef", "Vp_sign", "V_sign_f")        # evolve.py's own mutations, applied to the kernel the runs use
if MUT in KERNEL_MUTS:
    ev.MUT = MUT
EPS = np.array([1.0 if MUT == "eps_sign" else -1.0])      # eps = -1 (phantom); the mutation runs the kernel with the canonical sign

# --------------------------------------------------------------------------------------------------------------------
# declared design (README.md of I277-mac)
# --------------------------------------------------------------------------------------------------------------------
CASES = [(-1.5, +0.05), (-1.5, -0.05), (-1.0, +0.05), (-1.0, -0.05)]       # (c, beta), L = 1
XP = {-1.5: 2.16, -1.0: 3.20}       # x+ extent: v_end = 3.055 (v* = 3.137) and 4.525 (v* = 4.706); lam^2 there 0.076 and 0.108 at beta = -0.05
XM = 16.0                           # x- extent
RHO_C, V0 = 1.0, 0.0                # the corner: rho_c = 1 on the throat's outside, v0 = 0
HS = (0.02, 0.01, 0.005, 0.0025, 0.00125)
D4_BAR, D4_FALL = 1e-3, 2.5       # the right-ray constraint residual: bar at the finest h, least factor per halving
RATIO = (4.0, 0.3)                  # the error ratio per halving of h: 4.0 +- 0.3
VSTAR_TOL = 1e-3
BLK = 1024

DATA_MUTS = {"kh": dict(kh=1.01), "kf": dict(kf=1.01), "rate": dict(rate_scale=1.001)}      # F's (1/lam^2 - 1) coefficient, 2 rho lam'/lam coefficient (x1.01), rate K (x1.001)
KRN = ["D2", "O", "V"]
MUTATIONS = {         # name: the checks that must fire in `check`; computed in the development runs dev_11_mut_<name> of negative-mass-raw I277-mac before they were declared
    "beta_sign": ["O", "V"],                                      # data from -beta, compared with +beta (the issue's sign mutation)
    "kh": ["D1", "D1b", "D2", "D3", "D4", "O", "V"],              # coefficient 1 of (1/lam^2 - 1) in F -> 1.01
    "kf": ["D1", "D1b", "D2", "D3", "D4", "O", "V"],              # coefficient 2 of rho lam'/lam in F -> 2.02
    "rate": ["D1", "D1b", "D3", "D4"],                            # rate K -> 1.001 K inside lam^2 only: solves the evolution equations, not the constraints
    "gauge_frozen": ["D4", "O", "V"],                             # right ray's e^{-f} = F(v0, rho(v)) instead of the Jacobian flow
    "th_coef": KRN,                                               # evolve.py's own kernel mutations, applied to the kernel the runs use
    "vterm_sign": KRN,
    "vterm_coef": KRN,
    "Vp_sign": KRN,
    "V_sign_f": KRN,
    "eps_sign": KRN,                                              # the kernel run with eps = +1 (canonical scalar)
    "vfit_K": ["PL2", "V"],                                       # the v* finder with the rate scaled by 1.01
}


class Family:
    """lam^2 = L + beta e^{rate v};  F = A - kf 2 rho mu + kh (1/lam^2 - 1);  mu = lam'/lam.  kf = kh = 1, rate = K is the family of #269; others are mutations."""

    def __init__(self, c, beta, L=1.0, rate_scale=1.0, kf=1.0, kh=1.0):
        self.c = float(c)
        self.K = -2.0 * self.c / math.pi
        self.beta, self.L = float(beta), float(L)
        self.rate = self.K * rate_scale
        self.kf, self.kh = kf, kh

    def A(self, r):
        return 1.0 + self.c * (r * r + 1.0) + self.K * (r + (r * r + 1.0) * np.arctan(r))

    def Ap(self, r):
        return 2.0 * self.c * r + 2.0 * self.K * (1.0 + r * np.arctan(r))

    def App(self, r):
        return 2.0 * self.c + 2.0 * self.K * (np.arctan(r) + r / (1.0 + r * r))

    def s(self, v):
        return self.beta * np.exp(self.rate * v)

    def lam2(self, v):
        return self.L + self.s(v)

    def mu(self, v):
        return 0.5 * self.rate * self.s(v) / self.lam2(v)

    def mudot(self, v):
        return 0.5 * self.rate ** 2 * self.s(v) * self.L / self.lam2(v) ** 2

    def F(self, v, r):
        return self.A(r) - self.kf * 2.0 * r * self.mu(v) + self.kh * (1.0 / self.lam2(v) - 1.0)

    def Fp(self, v, r):
        return self.Ap(r) - self.kf * 2.0 * self.mu(v)

    def Fv(self, v, r):
        return -self.kf * 2.0 * r * self.mudot(v) - self.kh * self.rate * self.s(v) / self.lam2(v) ** 2

    def vstar(self):
        return math.log(-self.L / self.beta) / self.rate if self.beta < 0 else None

    def jets(self, v, rho, w, q):
        """r, f, phi and their x+ and x- derivatives at (v, rho) for e^{-f} = -sqrt2 w, w = d rho/d x-, q = d w/d x- (gauge data); closed forms, by hand.
        r+ = sqrt2 D r with D = d_v + (F/2) d_rho along the outgoing ray; r- = w r_rho; w+ = F' w/sqrt2; q+ = (F'' w^2 + F' q)/sqrt2."""
        lam = np.sqrt(self.lam2(v))
        mu, mud = self.mu(v), self.mudot(v)
        R = np.sqrt(1.0 + rho * rho)
        R1, R2 = rho / R, 1.0 / R ** 3
        F, Fr, Fv, Frr = self.F(v, rho), self.Fp(v, rho), self.Fv(v, rho), self.App(rho)
        P1 = SQ2 / (1.0 + rho * rho)
        P2 = -2.0 * SQ2 * rho / (1.0 + rho * rho) ** 2
        J = {}
        J["r"], J["f"], J["phi"] = lam * R, -np.log(-SQ2 * w), SQ2 * np.arctan(rho)
        Dr = lam * mu * R + lam * R1 * F / 2.0
        J["rp"] = SQ2 * Dr
        J["rm"] = w * lam * R1
        J["rpm"] = SQ2 * w * lam * (mu * R1 + F * R2 / 2.0 + Fr * R1 / 2.0)
        J["rpp"] = 2.0 * lam * ((mu * mu + mud) * R + mu * R1 * F / 2.0 + R1 * Fv / 2.0 + F / 2.0 * (mu * R1 + R2 * F / 2.0 + R1 * Fr / 2.0))
        J["rmm"] = q * lam * R1 + w * w * lam * R2
        J["fp"] = -Fr / SQ2
        J["fm"] = -q / w
        J["fpm"] = Frr * (-SQ2 * w) / 2.0
        J["pp"] = F / (1.0 + rho * rho)
        J["pm"] = w * P1
        J["ppm"] = SQ2 * (F / 2.0 * w * P2 + Fr * w / 2.0 * P1)
        return J


# --------------------------------------------------------------------------------------------------------------------
# data and exact fields
# --------------------------------------------------------------------------------------------------------------------
def left_ray(fam, rho_c, v0, N, h):
    """rho on the ingoing ray v = v0 at x- = j h, drho/dx- = -F(v0, rho)/sqrt2, rho(0) = rho_c."""
    from scipy.integrate import solve_ivp
    xs = np.arange(N + 1) * h
    sol = solve_ivp(lambda x, y: -fam.F(v0, y) / SQ2, [0.0, xs[-1]], [rho_c], method="DOP853", rtol=1e-13, atol=1e-15, t_eval=xs)
    if not sol.success or sol.y.shape[1] != len(xs):
        raise RuntimeError("left ray: " + sol.message)
    return sol.y[0]


def flow(fam, rho0, v0, vnodes, rtol=1e-12, atol=1e-13):
    """Outgoing rays labelled by their rho0 at v0: (rho, ln e^{-f}) at the v nodes, arrays (len(vnodes), len(rho0))."""
    from scipy.integrate import solve_ivp
    n = len(rho0)
    y0 = np.concatenate([rho0, np.log(fam.F(v0, rho0))])

    def rhs(v, y):
        r = y[:n]
        return np.concatenate([0.5 * fam.F(v, r), 0.5 * fam.Fp(v, r)])
    sol = solve_ivp(rhs, [v0, vnodes[-1]], y0, method="DOP853", rtol=rtol, atol=atol, t_eval=vnodes)
    if not sol.success or sol.y.shape[1] != len(vnodes):
        raise RuntimeError("flow: " + sol.message)
    return sol.y[:n].T, sol.y[n:].T


def fields(fam, v, rho, lnE):
    """(r, f, phi) of the family at (v, rho) with e^{-f} = exp(lnE); v broadcast against rho."""
    return np.sqrt(fam.lam2(v)) * np.sqrt(1.0 + rho * rho), -lnE, SQ2 * np.arctan(rho)


def ray_data(fam, rho_c, v0, N, h, Nl, frozen=False):
    """Characteristic data (Sp on x- = 0, Sm on x+ = 0) of the family: arrays (3, N + 1) and (3, Nl + 1) of (r, f, phi), and the left ray's rho0.
    frozen (a mutation): the right ray's e^{-f} taken as F(v0, rho(v)), the static gauge rule, in place of the Jacobian flow."""
    vn = v0 + SQ2 * h * np.arange(N + 1)
    rho_r, lnE_r = flow(fam, np.array([rho_c]), v0, vn)
    if frozen:
        lnE_r = np.log(fam.F(v0, rho_r))
    Sp = np.array(fields(fam, vn[:, None], rho_r, lnE_r))[:, :, 0]
    rho0 = left_ray(fam, rho_c, v0, Nl, h)
    Sm = np.array(fields(fam, v0, rho0, np.log(fam.F(v0, rho0))))
    return Sp, Sm, rho0


def rows(Sp, Sm, h, par, eps, rcut=0.0, fcut=1e300, excise=1):
    """The code's row sweep, yielding (j, row j of (r, f, phi) over x+ = i h, dead mask); row 0 is Sp.  The kernel is evolve.row_step, unchanged;
    no cell is excised unless a value is not finite (rcut = 0, fcut = 1e300)."""
    nf, n1 = Sp.shape
    old = np.ascontiguousarray(Sp.copy())
    new = np.empty_like(old)
    dold = np.zeros(n1, np.uint8)
    dnew = np.zeros(n1, np.uint8)
    yield 0, old, dold
    for j in range(Sm.shape[1] - 1):
        new[:, 0] = Sm[:, j + 1]
        dnew[0] = 0
        ev.row_step(old, new, dold, dnew, h, h, par, eps, rcut, fcut, excise)
        yield j + 1, new, dnew
        old, new = new, old
        dold, dnew = dnew, dold


def lam2_of(G):
    """Lambda^2 = r^2 cos^2(phi/sqrt2) of the evolved fields (first axis r, f, phi)."""
    return G[0] ** 2 * np.cos(G[2] / SQ2) ** 2


def vstar_fit(vn, lam2, K):
    """v* from a least-squares fit of Lambda^2(v) = a + b e^{Kv}: v* = ln(-a/b)/K (nan when -a/b <= 0); returns (v*, a, b)."""
    X = np.stack([np.ones_like(vn), np.exp(K * vn)], axis=1)
    (a, b), *_ = np.linalg.lstsq(X, lam2, rcond=None)
    return (math.log(-a / b) / K if -a / b > 0 else float("nan")), float(a), float(b)


def families(c, beta, mut):
    """(data family, truth family) of a run; the mutations of the data path act here."""
    if mut == "beta_sign":
        return Family(c, -beta), Family(c, beta)
    kw = DATA_MUTS.get(mut, {})
    return Family(c, beta, **kw), Family(c, beta, **kw)


def run_case(c, beta, h, Xp=None, Xm=XM, rho_c=RHO_C, v0=V0, mut=None):
    """One evolution of the family's data by the code and the comparison with the exact family at every node of the rectangle [0, Xp] x [0, Xm] (streamed)."""
    Xp = XP[c] if Xp is None else Xp
    data, truth = families(c, beta, mut)
    mem = ep.Member(c)
    Np, Nm = int(round(Xp / h)), int(round(Xm / h))
    t0 = time.time()
    Sp, Sm, rho0_d = ray_data(data, rho_c, v0, Np, h, Nm, frozen=(mut == "gauge_frozen"))
    rho0_t = left_ray(truth, rho_c, v0, Nm, h) if truth is not data else rho0_d
    vn = v0 + SQ2 * h * np.arange(Np + 1)
    par, eps = ev.params(mem), EPS
    emax, amax = np.zeros(3), [(0, 0)] * 3
    runmax = np.zeros((3, Np + 1))
    relmax = 0.0
    lsum, ndead = np.zeros(Np), 0
    ldev = 0.0
    lam2_ex = truth.lam2(vn[1:])
    blk = None
    for j, row, dead in rows(Sp, Sm, h, par, eps):
        if j % BLK == 0:
            rho_b, lnE_b = flow(truth, rho0_t[j:j + BLK], v0, vn)
            blk = fields(truth, vn[:, None], rho_b, lnE_b)
        jj = j % BLK
        ex = np.stack([blk[0][:, jj], blk[1][:, jj], blk[2][:, jj]])
        d = np.abs(row - ex)
        ndead += int(dead.sum())
        runmax = np.maximum(runmax, d)
        for q in range(3):
            m = float(d[q].max())
            if m > emax[q]:
                emax[q], amax[q] = m, (int(d[q].argmax()), j)
        relmax = max(relmax, float((d[0] / ex[0]).max()))
        if j >= 1:
            L2 = lam2_of(row)[1:]
            lsum += L2
            ldev = max(ldev, float(np.abs(L2 - lam2_ex).max()))
    Lm = lsum / Nm
    K = ep.Member(c).K * (1.01 if mut == "vfit_K" else 1.0)
    vf, a, b = vstar_fit(vn[1:], Lm, K)
    vend = float(vn[-1])
    l2e = float(Lm[-1])
    vlast = vend - math.log(1.0 - l2e) / K if l2e < 1 else float("nan")
    quarter = [float(runmax[0][:int(round(Np * x)) + 1].max()) for x in (0.25, 0.5, 0.75, 1.0)]
    return dict(c=c, beta=beta, h=h, Xp=Xp, Xm=Xm, Np=Np, Nm=Nm, v_end=vend, emax=[float(x) for x in emax], amax=[[int(a_[0]), int(a_[1])] for a_ in amax],
                rel_r=relmax, quarter_r=quarter, lam2_dev=ldev, vfit=vf, fit_a=a, fit_b=b, vlast=vlast, lam2_end=l2e, ndead=ndead, vstar=truth.vstar(),
                time=time.time() - t0)


# --------------------------------------------------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------------------------------------------------
RESULTS = []


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)
    RESULTS.append((cid, ok))
    tag = "SKIP" if ok is None else ("pass" if ok else "FAIL")
    print(f"[{tag}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def want(cid):
    return ONLY is None or cid in ONLY


def load_em():
    """exact_mode.py (#269): its sympy closed forms of the family, imported read-only; it reads sys.argv at import."""
    global _EM
    if _EM is None:
        keep = sys.argv
        sys.argv = [sys.argv[0]]
        import exact_mode as em
        sys.argv = keep
        _EM = em
    return _EM


_EM = None


def rng_points(rng, n, vmax):
    return rng.uniform(0.0, vmax, n), rng.uniform(-1.2, 6.0, n), -rng.uniform(0.2, 3.0, n), rng.normal(size=n)       # v, rho, w < 0, q


def sympy_jets(c, beta):
    """The jets of exact_mode's closed forms by the chain rule in sympy: the operators d+ and d- of the module docstring, applied to r(v, rho), F(v, rho), phi(rho)."""
    import sympy as sp
    em = load_em()
    cq = sp.Rational(str(c))
    bq = sp.Rational(str(beta))
    rexpr, Fexpr, Kq, V0, V1, Wq = em.family_exprs(cq, sp.Integer(1), bq)
    v, rho = em.v, em.rho
    w, q = sp.symbols("w q")
    S2 = sp.sqrt(2)
    Fr, Frr = sp.diff(Fexpr, rho), sp.diff(Fexpr, rho, 2)

    def dp(g):
        return S2 * (sp.diff(g, v) + Fexpr / 2 * sp.diff(g, rho) + Fr * w / 2 * sp.diff(g, w) + (Frr * w ** 2 + Fr * q) / 2 * sp.diff(g, q))
    Phi = S2 * sp.atan(rho)
    f = -sp.log(-S2 * w)
    rp = dp(rexpr)
    rm = w * sp.diff(rexpr, rho)
    pm = w * sp.diff(Phi, rho)
    J = dict(r=rexpr, f=f, phi=Phi, rp=rp, rm=rm, rpm=dp(rm), rpp=dp(rp), rmm=q * sp.diff(rexpr, rho) + w ** 2 * sp.diff(rexpr, rho, 2), fp=dp(f), fm=-q / w, fpm=dp(-q / w),
             pp=dp(Phi), pm=pm, ppm=dp(pm))
    return {k: sp.lambdify((v, rho, w, q), e, "mpmath") for k, e in J.items()}


def check_D1():
    """closed forms r, F, V, V' of the module against exact_mode's sympy forms at 50 digits"""
    import mpmath as mp
    import sympy as sp
    em = load_em()
    rng = np.random.default_rng(277)
    worst = 0.0
    for c, beta in CASES:
        fam = Family(c, beta, **DATA_MUTS.get(MUT, {}))
        cq, bq = sp.Rational(str(c)), sp.Rational(str(beta))
        rexpr, Fexpr, Kq, V0, V1, Wq = em.family_exprs(cq, sp.Integer(1), bq)
        fr = sp.lambdify((em.v, em.rho), [rexpr, Fexpr], "mpmath")
        fV = sp.lambdify(em.rho, [V0, V1], "mpmath")
        vend = XP[c] * SQ2
        for _ in range(20):
            v, rho = float(rng.uniform(0.0, vend)), float(rng.uniform(-1.2, 6.0))
            ra, Fa = fr(mp.mpf(v), mp.mpf(rho))
            Va, Vpa = fV(mp.mpf(rho))
            mem = ep.Member(c)
            Vb, Vpb = ep.pot(SQ2 * math.atan(rho), c, mem.K)
            mine = [float(np.sqrt(fam.lam2(v)) * math.sqrt(1 + rho * rho)), float(fam.F(v, rho)), Vb, Vpb]
            ref = [ra, Fa, Va, Vpa]
            for x, y in zip(mine, ref):
                worst = max(worst, float(abs(x - y) / (1 + abs(y))))
    check("D1", "the module's closed forms r = lam R, F, and the potential V(phi(rho)), V'(phi) (evolve_potential.pot) equal exact_mode's sympy forms at 50 digits, 4 cases x 20 points",
          worst <= 1e-13, f"max relative difference {worst:.1e}, bar 1e-13")


def check_D1b():
    """hand-written jets against the chain-rule jets of exact_mode's closed forms"""
    import mpmath as mp
    rng = np.random.default_rng(278)
    worst, which = 0.0, ""
    for c, beta in CASES:
        fam = Family(c, beta, **DATA_MUTS.get(MUT, {}))
        sj = sympy_jets(c, beta)
        v, rho, w, q = rng_points(rng, 12, XP[c] * SQ2)
        J = fam.jets(v, rho, w, q)
        for i in range(len(v)):
            for k, fn in sj.items():
                ref = float(fn(mp.mpf(float(v[i])), mp.mpf(float(rho[i])), mp.mpf(float(w[i])), mp.mpf(float(q[i]))))
                e = abs(J[k][i] - ref) / (1 + abs(ref))
                if e > worst:
                    worst, which = float(e), f"{k} at c = {c}, beta = {beta}"
    check("D1b", "the module's hand-written jets (r, f, phi and their x+, x- derivatives, mixed and second) equal the chain-rule jets of exact_mode's closed forms (sympy, 50 digits), 4 cases x 12 points x 13 jets",
          worst <= 1e-12, f"max relative difference {worst:.1e} ({which}), bar 1e-12")


def kernel_residuals(c, J):
    """F_np and the numba _F of the code on the jets: residuals against r+-, f+-, phi+-, for arrays of points."""
    mem = ep.Member(c)
    par, eps = ev.params(mem), EPS
    n = len(J["r"])
    cv = np.stack([J["r"], J["f"], J["phi"]])
    up = np.stack([J["rp"], np.zeros(n), J["pp"]])
    um = np.stack([J["rm"], np.zeros(n), J["pm"]])
    g = ev.F_np(cv, up, um, par, eps)
    gn = np.empty_like(g)
    tmp = np.empty(3)
    for i in range(n):
        ev._F(cv[:, i].copy(), up[:, i].copy(), um[:, i].copy(), par, eps, tmp)
        gn[:, i] = tmp
    ref = np.stack([J["rpm"], J["fpm"], J["ppm"]])
    return np.abs(g - ref) / (1 + np.abs(ref)), np.abs(gn - ref) / (1 + np.abs(ref)), np.abs(g - gn)


def check_D2():
    rng = np.random.default_rng(279)
    worst = 0.0
    lines = []
    for c, beta in CASES:
        fam = Family(c, beta, **DATA_MUTS.get(MUT, {}))
        v, rho, w, q = rng_points(rng, 400, XP[c] * SQ2)
        J = fam.jets(v, rho, w, q)
        rn, rb, dd = kernel_residuals(c, J)
        worst = max(worst, float(rn.max()), float(rb.max()))
        lines.append(f"c = {c}, beta = {beta:+}: F_np residual r+- {rn[0].max():.1e} f+- {rn[1].max():.1e} phi+- {rn[2].max():.1e}; _F (numba) {rb.max():.1e}")
    check("D2", "the code's equations (F_np, and the numba _F of the kernel) on the family's jets: the residuals of r+-, f+-, phi+- vanish, 4 cases x 400 random points (v, rho, e^{-f}, gauge q)",
          worst <= 1e-12, f"max relative residual {worst:.1e}, bar 1e-12; " + "; ".join(lines))


def check_D3():
    rng = np.random.default_rng(280)
    worst = 0.0
    lines = []
    for c, beta in CASES:
        fam = Family(c, beta, **DATA_MUTS.get(MUT, {}))
        v, rho, w, q = rng_points(rng, 400, XP[c] * SQ2)
        J = fam.jets(v, rho, w, q)
        # C++ = r++ + f+ r+ + (r/2) eps phi+^2, C-- likewise, eps = -1; relative to the largest term
        tp = np.stack([J["rpp"], J["fp"] * J["rp"], -0.5 * J["r"] * J["pp"] ** 2])
        tm = np.stack([J["rmm"], J["fm"] * J["rm"], -0.5 * J["r"] * J["pm"] ** 2])
        cp = np.abs(tp.sum(axis=0)) / np.abs(tp).max(axis=0)
        cm = np.abs(tm.sum(axis=0)) / np.abs(tm).max(axis=0)
        worst = max(worst, float(cp.max()), float(cm.max()))
        lines.append(f"c = {c}, beta = {beta:+}: C++ {cp.max():.1e}, C-- {cm.max():.1e}")
    check("D3", "the family's jets satisfy the null constraints C++ = r++ + f+ r+ + (r/2) eps phi+^2 and C-- (eps = -1) identically, 4 cases x 400 random points",
          worst <= 1e-12, f"max residual relative to the largest term {worst:.1e}, bar 1e-12; " + "; ".join(lines))


def ratios_ok(vals):
    rs = [vals[i] / vals[i + 1] for i in range(len(vals) - 1)]
    return rs, all(RATIO[0] - RATIO[1] <= r <= RATIO[0] + RATIO[1] for r in rs)


def check_D4():
    """the null constraint C++ on the stored right-ray data (central differences), relative to its largest term: small, and falling with h.  The left ray is not used:
    its derivatives are exponentially small at large x-, so the second differences there sit on the roundoff floor (dev_09)."""
    allok = True
    lines = []
    for c, beta in CASES:
        data, _ = families(c, beta, MUT)
        sup = []
        for h in HS:
            Np = int(round(XP[c] / h))
            Sp, _, _ = ray_data(data, RHO_C, V0, Np, h, 1, frozen=(MUT == "gauge_frozen"))
            r, f, p = Sp
            rpp = (r[2:] - 2 * r[1:-1] + r[:-2]) / h ** 2
            fp = (f[2:] - f[:-2]) / (2 * h)
            rp = (r[2:] - r[:-2]) / (2 * h)
            pp = (p[2:] - p[:-2]) / (2 * h)
            t = np.stack([rpp, fp * rp, -0.5 * r[1:-1] * pp ** 2])
            sup.append(float((np.abs(t.sum(axis=0)) / np.abs(t).max(axis=0)).max()))
        fall = [sup[i] / sup[i + 1] for i in range(len(sup) - 1)]
        ok = sup[-1] <= D4_BAR and all(x >= D4_FALL for x in fall)
        allok &= bool(ok)
        lines.append(f"c = {c}, beta = {beta:+}: sup|C++|/largest term " + ", ".join(f"{x:.2e}" for x in sup) + " (falls by " + " ".join(f"{x:.2f}" for x in fall) + ")")
    check("D4", f"the null constraint C++ on the stored right-ray data (second-order central differences), relative to its largest term, is at most {D4_BAR:g} at h = {HS[-1]} and falls by at least {D4_FALL} "
          f"per halving of h = {HS}", allok, "; ".join(lines))


def check_D5():
    """the exact map's own accuracy: the flow of the outgoing rays with the tolerances tightened"""
    worst = 0.0
    for c, beta in CASES:
        fam = Family(c, beta)
        h = 0.01
        Np, Nm = int(round(XP[c] / h)), 1600
        rho0 = left_ray(fam, RHO_C, V0, Nm, h)
        vn = V0 + SQ2 * h * np.arange(Np + 1)
        a = flow(fam, rho0, V0, vn)
        b = flow(fam, rho0, V0, vn, rtol=1e-13, atol=1e-14)
        fa, fb = fields(fam, vn[:, None], *a), fields(fam, vn[:, None], *b)
        worst = max(worst, max(float(np.abs(x - y).max()) for x, y in zip(fa, fb)))
    check("D5", "the exact map (outgoing-ray flow, rtol 1e-12) changes by less than 1e-10 in r, f, phi when the tolerances are tightened to 1e-13 (the code's errors are 1e-7 and above)",
          worst <= 1e-10, f"max change {worst:.1e}, bar 1e-10")


def check_D6():
    """the module's row sweep equals evolve.run_grid on the family's data"""
    worst, same = 0.0, True
    for c, beta in CASES:
        N, h = 100, 0.02
        data, _ = families(c, beta, MUT)
        Sp, Sm, _ = ray_data(data, RHO_C, V0, N, h, N, frozen=(MUT == "gauge_frozen"))
        d = dn.Data(2, list(EPS), N, h, Sp, Sm, name="family")
        mem = ep.Member(c)
        par = ev.params(mem)
        for kernel in ("numba", "numpy"):
            G, dead = ev.run_grid(d, par, kernel, rcut=0.0, fcut=1e300)
            mine = np.empty_like(G)
            for j, row, dd in rows(Sp, Sm, h, par, EPS):
                mine[:, :, j] = row
            diff = float(np.max(np.abs(G - mine)))
            if kernel == "numba":
                same &= diff == 0.0
            worst = max(worst, diff) if kernel == "numpy" else worst
            same &= (diff <= 1e-12 if kernel == "numpy" else True)
    check("D6", "the module's row sweep equals evolve.run_grid on the family's data: bit-identical with the numba kernel, 1e-12 with the numpy reference kernel (4 cases, N = 100, h = 0.02)",
          same, f"numpy max |dG| {worst:.1e}")


def check_PL1():
    """the ratio estimator on planted errors"""
    hs = np.array(HS)
    ok = []
    notes = []
    for p, expect in ((2.0, True), (1.8, False), (2.2, False), (1.0, False)):
        rs, good = ratios_ok(list(0.37 * hs ** p))
        ok.append(good == expect)
        notes.append(f"h^{p}: ratio {rs[0]:.3f} -> {'accepted' if good else 'rejected'}")
    check("PL1", "planted errors C h^p through the same ratio band: p = 2 accepted; p = 1.8, 2.2, 1 rejected", all(ok), "; ".join(notes))


def check_PL2():
    """the v* finder on planted Lambda^2(v) = 1 + beta e^{Kv}"""
    rng = np.random.default_rng(281)
    ok = True
    notes = []
    for c in (-1.5, -1.0):
        K = ep.Member(c).K * (1.01 if MUT == "vfit_K" else 1.0)
        K0 = ep.Member(c).K
        vn = V0 + SQ2 * 0.01 * np.arange(1, int(XP[c] / 0.01) + 1)
        for beta in (-0.05, -0.2, -0.01):
            vs = math.log(-1.0 / beta) / K0
            vn_ = vn[vn < 0.97 * vs]
            vf, a, b = vstar_fit(vn_, 1.0 + beta * np.exp(K0 * vn_), K)
            ok &= abs(vf - vs) <= 1e-10
            notes.append(f"c = {c}, beta = {beta}: {vf - vs:+.1e}")
        # noise: the finder should still find v* inside the bar for noise 1e-5 (a typical error of the code at h = 0.02 is 1e-4 in Lambda^2)
        beta = -0.05
        vs = math.log(-1.0 / beta) / K0
        errs = []
        for _ in range(40):
            vf, a, b = vstar_fit(vn, 1.0 + beta * np.exp(K0 * vn) + 1e-5 * rng.normal(size=len(vn)), K)
            errs.append(abs(vf - vs))
        ok &= max(errs) <= VSTAR_TOL
        notes.append(f"noise 1e-5, 40 seeds: max |error| {max(errs):.1e}")
        # a control that must fire: the finder with the wrong rate
        vfw, _, _ = vstar_fit(vn, 1.0 + beta * np.exp(K0 * vn), 1.01 * K0)
        ok &= abs(vfw - vs) > VSTAR_TOL
        notes.append(f"control (rate x 1.01): {vfw - vs:+.1e}")
    check("PL2", "the v* finder (fit of Lambda^2 to a + b e^{Kv}) recovers planted v* to 1e-10, stays inside 1e-3 at noise 1e-5 in 40 seeds, and a finder with the rate scaled by 1.01 misses by more than 1e-3",
          ok, "; ".join(notes))


RUN_RESULTS = {}


def check_runs():
    """the runs: O (error ratios), V (v* for beta < 0), and the report lines"""
    allO, allV = True, True
    vnotes = []
    for c, beta in CASES:
        res = []
        for h in HS:
            r = run_case(c, beta, h, mut=MUT)
            res.append(r)
            RUN_RESULTS[(c, beta, h)] = r
            print(f"  run c = {c}, beta = {beta:+}, h = {h}: rectangle {r['Np']} x {r['Nm']} (v_end {r['v_end']:.4f}); sup|dr| {r['emax'][0]:.4e} at (i, j) = {tuple(r['amax'][0])} "
                  f"(relative {r['rel_r']:.2e}); sup|df| {r['emax'][1]:.4e}; sup|dphi| {r['emax'][2]:.4e}; running sup|dr| at 1/4, 1/2, 3/4, 1 of v: "
                  + ", ".join(f"{x:.2e}" for x in r["quarter_r"]) + f"; sup|Lambda^2 - lam^2| {r['lam2_dev']:.2e}; dead cells {r['ndead']}; {r['time']:.1f} s", flush=True)
        rs, ok = ratios_ok([x["emax"][0] for x in res])
        rf, okf = ratios_ok([x["emax"][1] for x in res])
        rp, okp = ratios_ok([x["emax"][2] for x in res])
        allO &= ok
        print(f"  c = {c}, beta = {beta:+}: error ratios per halving of h, r: " + " ".join(f"{x:.3f}" for x in rs) + "; f: " + " ".join(f"{x:.3f}" for x in rf)
              + "; phi: " + " ".join(f"{x:.3f}" for x in rp) + f"   (r {'within' if ok else 'OUTSIDE'} 4.0 +- 0.3; f, phi reported, f {'within' if okf else 'outside'}, phi {'within' if okp else 'outside'})", flush=True)
        if beta < 0:
            vs = res[0]["vstar"]
            errs = [x["vfit"] - vs for x in res]
            okv = all(abs(e) <= VSTAR_TOL for e in errs)       # nan compares False
            allV &= okv
            rv = [abs(errs[i] / errs[i + 1]) if errs[i + 1] else float("nan") for i in range(len(errs) - 1)]
            print(f"  c = {c}, beta = {beta:+}: v* = {vs:.6f}; found (fit of the evolved interior Lambda^2) " + ", ".join(f"{x['vfit']:.6f}" for x in res) + "; errors " + ", ".join(f"{e:+.2e}" for e in errs)
                  + "; error ratios " + " ".join(f"{x:.2f}" for x in rv) + "; fit a, b at the finest h " + f"{res[-1]['fit_a']:.6f}, {res[-1]['fit_b']:.6f}; direct estimate v - ln(1 - Lambda^2)/K from the last row "
                  + ", ".join(f"{x['vlast'] - vs:+.1e}" for x in res), flush=True)
            vnotes.append(f"c = {c}: max |v* error| {max(abs(e) for e in errs):.1e}")
    check("O", f"the evolved r equals lam(v) R(rho) over the whole rectangle: sup error falls by 4.0 +- 0.3 per halving of h = {HS}, four cases (c = -1.5, -1; beta = +-0.05)", allO,
          "per-case lines above")
    check("V", f"for beta < 0, v* = ln(-L/beta)/K is found from the evolved interior Lambda^2 to {VSTAR_TOL:g}, at every h and both members", allV, "; ".join(vnotes))


def rays():
    """Where the exact map ends.  beta < 0: along outgoing rays (labelled by rho0 at v0 = 0) rho -> +infinity as v -> v*, 1 + K rho ~ (v* - v)^{-1/2}, r = lam R tends to a constant
    (printed against eps = v* - v).  beta > 0: the v at which an outgoing ray reaches rho = -1000 (the integration of the exact map stops near there), against the rectangle's v_end."""
    from scipy.integrate import solve_ivp
    for c, beta in CASES:
        fam = Family(c, beta)
        K = fam.K
        h = 0.005
        Nm = int(round(XM / h))
        rho0 = left_ray(fam, RHO_C, V0, Nm, h)
        vend = SQ2 * XP[c]
        picks = [("x- = 0 (rho0 = 1)", rho0[0]), ("x- = 8", rho0[int(8 / h)]), ("x- = 16", rho0[-1])]
        if beta < 0:
            vs = fam.vstar()
            eps = [1.0, 1e-2, 1e-4, 1e-6, 1e-8]
            vn = np.array([vs - e for e in eps])
            print(f"c = {c}, beta = {beta:+}: v* = {vs:.6f}; outgoing rays at v = v* - eps")
            for name, r0 in picks:
                rho, lnE = flow(fam, np.array([r0]), V0, vn, rtol=1e-11, atol=1e-12)
                rr = np.sqrt(fam.lam2(vn)) * np.sqrt(1 + rho[:, 0] ** 2)
                print(f"  {name}: rho0 = {r0:.6f}; " + "; ".join(f"eps {e:g}: rho {x:.5g}, (1 + K rho) sqrt(K eps) {(1 + K * x) * math.sqrt(K * e):.5f}, r {rv:.5f}"
                                                              for e, x, rv in zip(eps, rho[:, 0], rr)))
        else:
            print(f"c = {c}, beta = {beta:+}: the rectangle's v_end = {vend:.4f}; escape of the outgoing rays to rho = -1000")
            for name, r0 in picks:
                ev_ = lambda v, y: y[0] + 1000.0
                ev_.terminal = True
                sol = solve_ivp(lambda v, y: [0.5 * fam.F(v, y[0]), 0.5 * fam.Fp(v, y[0])], [V0, 30.0], [r0, math.log(fam.F(V0, r0))], method="DOP853", rtol=1e-11, atol=1e-12, events=ev_)
                ve = float(sol.t_events[0][0]) if len(sol.t_events[0]) else float("nan")
                print(f"  {name}: rho0 = {r0:.6f}; reaches rho = -1000 at v = {ve:.4f}")


def main():
    cmd = next((a for a in sys.argv[1:] if not a.startswith("--")), "check")
    if cmd == "rays":
        return rays() or 0
    if cmd == "mutations":
        return report_mutations(sys.argv[2], sys.argv[3])
    if cmd == "case":
        kv = dict(a[2:].split("=", 1) for a in sys.argv[2:] if a.startswith("--") and "=" in a)
        c, beta, h = float(kv["c"]), float(kv["beta"]), float(kv["h"])
        r = run_case(c, beta, h, Xp=float(kv["Xp"]) if "Xp" in kv else None, Xm=float(kv.get("Xm", XM)), mut=MUT)
        print(json.dumps(r, indent=1))
        return 0
    print(f"planted_family.py check; mutation {MUT}; python {sys.version.split()[0]}; numpy {np.__version__}; {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    for cid, fn in (("D1", check_D1), ("D1b", check_D1b), ("D2", check_D2), ("D3", check_D3), ("D4", check_D4), ("D5", check_D5), ("D6", check_D6), ("PL1", check_PL1), ("PL2", check_PL2)):
        if want(cid):
            fn()
    if want("O") or want("V"):
        check_runs()
    fails = [cid for cid, ok in RESULTS if ok is False]
    print(f"\n{sum(1 for _, ok in RESULTS if ok)} of {len(RESULTS)} pass; fired: {','.join(fails) if fails else 'none'}", flush=True)
    if OUT and RUN_RESULTS:
        os.makedirs(OUT, exist_ok=True)
        with open(os.path.join(OUT, f"results_{MUT or 'none'}.json"), "w") as fh:
            json.dump([v for v in RUN_RESULTS.values()], fh, indent=1)
    return 1 if fails else 0


def report_mutations(d, prefix):
    """fired sets of the stored mutated logs d/prefix<name>.log (the checks printed [FAIL]) against the declared ones"""
    nok = 0
    for name, declared in MUTATIONS.items():
        path = os.path.join(d, f"{prefix}{name}.log")
        if not os.path.exists(path):
            print(f"[MISSING] {name}: {path}")
            continue
        fired = sorted({m.group(1) for m in re.finditer(r"^\[FAIL\] (\w+):", open(path).read(), re.M)})
        ok = fired == sorted(declared)
        nok += ok
        print(f"[{'pass' if ok else 'FAIL'}] mutation {name}: fired {','.join(fired) or 'none'}; declared {','.join(sorted(declared))}")
    print(f"\n{nok} of {len(MUTATIONS)} mutations fire exactly their declared checks")
    return 0 if nok == len(MUTATIONS) else 1


if __name__ == "__main__":
    sys.exit(main())
