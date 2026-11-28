#!/usr/bin/env python3
"""The asymptotic solution of the mirror-symmetric expanding throat and its exponent k(n) (#230). Uses expansion.py, dn.py and checks.py unchanged.
The declarations (criteria, members, h, X, planted cases, mutations) are CRIT and MUTATIONS below, as committed in negative-mass-raw
I230-dh610/README.md before the runs of record.

    python3 throat_scaling.py derive                                       the closed form against the flat-limit field equations (sympy), the Ellis jump, k(n), c(n)
    python3 throat_scaling.py job --D=6 --A=0.05 --h=0.003 --X=200 --out=DIR    one streamed run of family S (A = from:PATH takes the member of a stored run;
                                                                           --b=0.2 adds a step of 0.2 to the field's end value, the jump test)
    python3 throat_scaling.py analyze --i213=DIR --i230=DIR                every check and table, from the stored runs (exit 1 if a check fails)
    python3 throat_scaling.py mutate --mutate=NAME                         the planted checks and a reduced suite with one input broken
    python3 throat_scaling.py dev                                          the planted checks and the reduced suite, nothing broken

Units a = 1, 8 pi G = 1, tau the proper time along the diagonal x+ = x-, the throat of the mirror-symmetric run (#213). Conventions of dn.py:
ds^2 = -2 e^{-f} dx+ dx- + r^2 dOmega_n^2, R_ab = eps d_a phi d_b phi with eps = -1 (ghost).
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dn            # noqa: E402
import expansion as ex   # noqa: E402

SQ2 = math.sqrt(2.0)
MUT = None          # this file's mutation (dn.MUT holds dn.py's)


def mut(name):
    return MUT == name


CRIT = dict(
    # the closed form (symbolic) and the Ellis jump
    n_sym=(2, 3, 4, 5, 6), n_table=tuple(range(2, 9)), ellis_L=4000.0, ellis_tol=1e-9, charge_tol=1e-12,
    # estimators: the last e-fold of r of #213's fit; the extrapolation of the local exponent over the same window; the window for the stability check
    nfold=1.0, nfold_ws=1.5, ws_rel=1e-3,
    # planted controls
    PF=dict(k=2.0, ts=17.8902, tau_end_gap=0.31, npts=40000, sigma=1e-4, seeds=(1, 2, 3), k_tol=1e-3, ts_tol=1e-3),
    PX=dict(cases=((3, 2.0, 4e-3, 84.0), (2, 2.5, 2e-4, 190.0), (4, 2.4, 5e-3, 50.0)), npts=800, sigma=1e-4, seeds=(1, 2, 3), rel=1e-3),
    PS=dict(n=3, xmax=60.0, xmin=20.0, A_pert=1.02),
    # the acceptance on the stored and the new runs
    kc_rel=1e-3,          # #213's stored runs, D = 4 and 5: k and c within this relative distance of the closed form
    c_radius=50.0,        # |rho|/H^2 read at r = 50a (the radius of #213's table) or at 0.8 r_end if the run stops short
    k6_abs=1e-2,          # D = 6: k within 1e-2 (absolute) of 2(n - 1)
    res6=dict(H_rel=1e-2, k_abs=1e-2, r_eval=(10.0, 25.0, 50.0)),
    x800=dict(ts_ref=17.8902, ts_tol=1e-3, k_tol=1e-3, KT_rel=1e-3),
    sp=dict(xmin=20.0, r_rel=1e-3, phi_abs=1e-2),   # spatial: columns x+ = X/4, X/2, X at x- in [xmin, x+]
    # the runs
    runs=dict(D4_tuned_X800=dict(D=4, h=0.01, X=800.0), D6=dict(D=6, X=200.0, h=(0.003, 0.006), members=(0.05, 0.2, 0.01))),
    # the jump test: the S pulse (A = 0.05) plus a step b in the field's end value on each ray (a 'dphi' pulse, centre 3, width 0.5, odd mirror);
    # the closed form predicts A = (2/pi)(phi_inf + b), k = 2/(n A^2 - n - 1), c = n^2 A^2/2 for b with n A^2 > n + 1
    jump=dict(A=0.05, xc=3.0, w=0.5, tol=1e-3, members={4: dict(h=0.02, X=100.0, b=(0.05, 0.1, 0.2, 0.4, -0.05)),
                                                          5: dict(h=0.00707, X=60.0, b=(0.2, 0.5, -0.03)),
                                                          6: dict(h=0.012, X=60.0, b=(0.1, 0.3, -0.03))}),
    mini=dict(D=4, A=0.05, h=0.02, X=100.0, b=0.2),
)

# the checks each mutation must fire, exactly, on the reduced suite (D = 4, h = 0.02, X = 100: A = 0.05 and the jump test at b = 0.2) and the planted cases;
# computed in the development run before the declaration (negative-mass-raw I230-dh610/dev/mut/dev_05_*.log, BUILD ones rerun as dev_06_*.log)
MUTATIONS = {
    "exp_n": ["CF1", "CF2", "PS", "SP1"],                   # r = (x+ + x-)^{1/(n-1)} in the closed form and the spatial test
    "f_coef": ["CF1", "CF2"],                               # (n-1)/n replaced by 1/n in f
    "c_formula": ["C4", "CF2", "CF4", "JTC"],               # c = n^2 A^2/4 instead of n^2 A^2/2
    "k_formula": ["CF4", "K4F", "K4X"],                     # k = 2n instead of 2(n-1)
    "A_scale": ["C4", "CF3", "CF4", "JTC", "JTK", "PS", "SP2"],   # A^2 1 % high
    "jump_off": ["JTC", "JTK", "SP2"],                      # the step b ignored in the prediction
    "fit_unit": ["K4F", "PF"],                              # tau scaled by sqrt 2 in the fit
    "extrap_mean": ["PX"],                                  # the mean of k_loc over the window instead of the extrapolation
    "scalar_coef": ["BUILD"],                               # dn.py's kernel mutations: the reduced run collapses (trapped cells, r_end < 0.5) and does not expand
    "theta_coef": ["BUILD"],
    "ghost_sign": ["BUILD"],
    "data_constraint_sign": ["BUILD"],
    "ellis_charge": [],                                     # silent by design: q enters the constraint solve for r on the rays, not the end value of the field that sets A
    "no_corrector": [],                                     # silent at h = 0.02: the checks' bars (1e-3) are above the predictor's error here
    "constraint_sign": [],                                  # silent by design: diagnostics this script does not read
    "flux_sign": [],
    "mass_exponent": [],
    "finder_sign": [],
}


results = []      # (id, ok)


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok))
    tag = "SKIP" if ok is None else ("pass" if ok else "FAIL")
    print(f"[{tag}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


# --------------------------------------------------------------------------------------------------------------------
# the claimed values (mutable)
# --------------------------------------------------------------------------------------------------------------------
def A2_pred(n):
    """A^2 of the scaling solution: phi = A arcsin((x+ - x-)/(x+ + x-)), A pi = the jump of the static wormhole's field."""
    return (n / (n - 1.0)) * (1.0 + (0.02 if mut("A_scale") else 0.0))


def A_jump(n, b):
    """A of the scaling solution for data whose field has the end value phi_inf + b: A = (2/pi)(phi_inf + b), phi_inf = (pi/2) sqrt(n/(n-1))."""
    return math.sqrt(A2_pred(n)) + (0.0 if mut("jump_off") else 2.0 * b / math.pi)


def k_jump(n, b):
    den = n * A_jump(n, b) ** 2 - n - 1
    return 2.0 / den if den > 0 else float("inf")


def c_jump(n, b):
    return n * n * A_jump(n, b) ** 2 / (4.0 if mut("c_formula") else 2.0)


def k_pred(n):
    return 2.0 * n if mut("k_formula") else 2.0 * (n - 1)


def c_pred(n):
    return n * n * A2_pred(n) / (4.0 if mut("c_formula") else 2.0)


# --------------------------------------------------------------------------------------------------------------------
# the closed form against the flat-limit field equations (sympy; shares nothing with dn.py)
# --------------------------------------------------------------------------------------------------------------------
def closed_form(nn):
    """r^n = u + v, phi = A asin((u - v)/(u + v)), e^f = ... on u = x+, v = x- in the gauge where r^n is additive.
    Returns sympy objects (u, v, A, r, phi, f)."""
    import sympy as sp
    u, v = sp.symbols("u v", positive=True)
    A = sp.symbols("A", positive=True)
    n = sp.Integer(nn)
    S = u + v
    r = S ** (1 / (n - 1 if mut("exp_n") else n))
    phi = A * sp.asin((u - v) / S)
    f = ((1 / n if mut("f_coef") else (n - 1) / n) - A ** 2) * sp.log(S / 2) + A ** 2 * sp.log(u * v)
    return u, v, A, r, phi, f


def chk_CF1_CF2():
    """CF1: the closed form solves the three evolution equations (without the sphere term e^{-f}/2) and the two constraints, exactly.
    CF2: on the diagonal H = k/(tau_s - tau) with k = 2/(n A^2 - n - 1), |rho|/H^2 = n^2 A^2/2, and the throat relation holds identically."""
    import sympy as sp
    eps = -1
    bad1, bad2, lines = [], [], []
    for nn in CRIT["n_sym"]:
        u, v, A, r, phi, f = closed_form(nn)
        n = sp.Integer(nn)
        d = sp.diff
        res = [
            d(r, u, v) + (n - 1) * d(r, u) * d(r, v) / r,                                    # r_{+-} = -(n-1) r_+ r_-/r (sphere term dropped)
            d(f, u, v) - (n / r) * d(r, u, v) - eps * d(phi, u) * d(phi, v),                  # f_{+-} = (n/r) r_{+-} + eps phi_+ phi_-
            d(phi, u, v) + n * (d(r, u) * d(phi, v) + d(r, v) * d(phi, u)) / (2 * r),         # phi_{+-} = -n (r_+ phi_- + r_- phi_+)/(2 r)
            d(r, u, u) + d(f, u) * d(r, u) + (r / n) * eps * d(phi, u) ** 2,                  # C_{++}
            d(r, v, v) + d(f, v) * d(r, v) + (r / n) * eps * d(phi, v) ** 2,                  # C_{--}
        ]
        res = [sp.simplify(e) for e in res]
        if any(e != 0 for e in res):
            bad1.append(nn)
        # the diagonal: u = v = x
        x = sp.symbols("x", positive=True)
        sub = {u: x, v: x}
        H = sp.exp(f / 2) * (d(r, u) + d(r, v)) / (sp.sqrt(2) * r)
        rho = eps * sp.exp(f) * d(phi, u) ** 2
        dtau_dx = sp.sqrt(2) * sp.exp(-f / 2)
        Hx, rhox, dtx = (e.subs(sub) for e in (H, rho, dtau_dx))
        k_sym = sp.simplify(-1 / (sp.diff(1 / Hx, x) / dtx))                  # d(1/H)/dtau = -1/k
        c_sym = sp.simplify(-rhox / Hx ** 2)
        c_claim = sp.simplify(n ** 2 * A ** 2 / (4 if mut("c_formula") else 2))
        k_claim = 2 / (n * A ** 2 - n - 1)
        ok = (sp.simplify(k_sym - k_claim) == 0 and sp.simplify(c_sym - c_claim) == 0
              and sp.simplify(c_sym - n * (1 / k_sym + (n + 1) / sp.Integer(2))) == 0)
        if not ok:
            bad2.append(nn)
        lines.append(f"n = {nn}: residuals {[int(e == 0) for e in res]} (1 = zero), k = {k_sym}, c = {c_sym}")
    for ln in lines:
        print("    CF " + ln)
    check("CF1", "the closed form (r^n = x+ + x-, phi = A arcsin((x+ - x-)/(x+ + x-)), f in closed form) solves the flat-limit evolution equations and both constraints exactly, A free",
          not bad1, f"n = {list(CRIT['n_sym'])}; failing n: {bad1}")
    check("CF2", "on its diagonal H = k/(tau_s - tau) with k = 2/(nA^2 - n - 1), |rho|/H^2 = n^2 A^2/2, and the throat relation c = n(1/k + (n+1)/2) holds identically",
          not bad2, f"n = {list(CRIT['n_sym'])}; failing n: {bad2}")


def ellis_jump(n):
    """The static wormhole's field at l = L from dn.ellis_dense, the tail -q/((n-1) L^{n-1}) added back, and the charge."""
    L = CRIT["ellis_L"]
    sol, q = dn.ellis_dense(n, L)
    ph = sol(L)[2] + q / ((n - 1) * L ** (n - 1))
    return ph, q


def chk_CF3_CF4():
    """CF3: the static wormhole's end value of phi is (pi/2) sqrt(n/(n-1)) and its charge sqrt(n(n-1)), so A = (2/pi) phi_inf = sqrt(n/(n-1)).
    CF4: with A^2 = n/(n-1) the closed form gives k = 2(n-1) and c = n^3/(2(n-1)), the values the checks below use."""
    import sympy as sp
    worst, worst_q = 0.0, 0.0
    for n in CRIT["n_sym"]:
        ph, q = ellis_jump(n)
        pred = math.pi / 2 * math.sqrt(A2_pred(n))
        worst = max(worst, abs(ph - pred))
        worst_q = max(worst_q, abs(q - math.sqrt(n * (n - 1))))
        print(f"    CF3 n = {n}: phi_inf = {ph:.12f}, (pi/2) A = {pred:.12f}, q = {q:.12f}")
    check("CF3", "the static wormhole's end value phi_inf = (pi/2) sqrt(n/(n-1)) (dn.ellis_dense, tail restored) and q = sqrt(n(n-1))",
          worst <= CRIT["ellis_tol"] and worst_q <= CRIT["charge_tol"], f"max |phi_inf - (pi/2)A| {worst:.2e} (tol {CRIT['ellis_tol']:g}), max |q - sqrt(n(n-1))| {worst_q:.2e}")
    bad = []
    rows = []
    for n in CRIT["n_table"]:
        nn = sp.Integer(n)
        A2 = nn / (nn - 1)
        k = sp.simplify(2 / (nn * A2 - nn - 1))
        c = sp.simplify(nn ** 2 * A2 / 2)
        rows.append((n + 2, n, k, c, sp.simplify(nn * (1 / k + (nn + 1) / 2))))
        if sp.simplify(k - sp.nsimplify(k_pred(n))) != 0 or sp.simplify(c - sp.nsimplify(c_pred(n))) != 0 or sp.simplify(c - nn * (1 / k + (nn + 1) / 2)) != 0:
            bad.append(n)
    print("    CF4  D   n   k = 2(n-1)   c = n^3/(2(n-1))   c from the throat relation")
    for D, n, k, c, c2 in rows:
        print(f"    CF4  {D:<3d}{n:<4d}{str(k):<13s}{str(c):<19s}{c2}")
    check("CF4", "with A^2 = n/(n-1) the closed form gives k = 2(n-1) and c = n^3/(2(n-1)) (exact rationals), equal to the values the checks use",
          not bad, f"n = {list(CRIT['n_table'])}; failing n: {bad}")


# --------------------------------------------------------------------------------------------------------------------
# estimators on a throat series
# --------------------------------------------------------------------------------------------------------------------
def throat(path):
    R = ex.load_run(path)
    R.update(ex.finish(R))
    M = ex.diag_arrays(R)
    return R, M


def fit_k(tau, r, H, nfold=None):
    """#213's fit: 1/H linear in tau over the last e-fold of r. Returns k, tau_s."""
    ix = ex.final_window(tau, r, CRIT["nfold"] if nfold is None else nfold)
    t = tau[ix] * (SQ2 if mut("fit_unit") else 1.0)
    b, a = np.polyfit(t, 1 / H[ix], 1)
    return -1 / b, -a / b


def local_exponent(tau, H):
    """k_loc = H^2/(dH/dtau): the exponent without an extrapolated tau_s."""
    return H ** 2 / np.gradient(H, tau)


def k_inf(r, kloc, nfold=None):
    """Least-squares k_loc(r) = k_inf - a (r/r_end)^-g over the last e-fold of r; returns k_inf, g, rms of the residual.
    extrap_mean (mutation): the mean of k_loc over the window instead of the extrapolation."""
    w = CRIT["nfold"] if nfold is None else nfold
    sel = r >= r[-1] * math.exp(-w)
    rr, kk = r[sel], kloc[sel]
    step = max(1, len(rr) // 400)
    rr, kk = rr[::step], kk[::step]
    if mut("extrap_mean"):
        return float(np.mean(kk)), float("nan"), float(np.std(kk))
    best = None
    for g0 in (0.5, 1.0, 2.0, 3.0, 4.0):
        s = least_squares(lambda p: (p[0] - p[1] * (rr / rr[-1]) ** (-p[2])) - kk, [kk[-1], max(kk[-1] - kk[0], 1e-4), g0],
                          bounds=([0, -1e3, 0.05], [1e3, 1e3, 12]))
        if best is None or s.cost < best.cost:
            best = s
    return float(best.x[0]), float(best.x[2]), math.sqrt(2 * best.cost / len(rr))


def c_at(r, rho, H, rc):
    j = min(int(np.searchsorted(r, rc)), len(r) - 1)
    return float(-rho[j] / H[j] ** 2), float(r[j])


# --------------------------------------------------------------------------------------------------------------------
# planted controls
# --------------------------------------------------------------------------------------------------------------------
def chk_PF():
    """PF: a synthetic H = k/(tau_s - tau) with noise through the fit recovers k and tau_s."""
    c = CRIT["PF"]
    worst_k = worst_t = 0.0
    for seed in c["seeds"]:
        rng = np.random.default_rng(seed)
        tau = np.linspace(c["ts"] - 6.0, c["ts"] - c["tau_end_gap"], c["npts"])
        s = c["ts"] - tau
        r = 7.0 * s ** (-c["k"])
        H = c["k"] / s * (1 + c["sigma"] * rng.standard_normal(len(s)))
        k, ts = fit_k(tau, r, H)
        worst_k, worst_t = max(worst_k, abs(k - c["k"])), max(worst_t, abs(ts - c["ts"]))
    check("PF", "planted: a synthetic H = k/(tau_s - tau) with relative noise through the fit recovers k and tau_s",
          worst_k <= c["k_tol"] and worst_t <= c["ts_tol"], f"worst |k - {c['k']}| {worst_k:.2e} (tol {c['k_tol']:g}), worst |tau_s - {c['ts']}| {worst_t:.2e} (tol {c['ts_tol']:g}); "
          f"noise {c['sigma']:g}, {len(c['seeds'])} seeds, window = last e-fold of r")


def planted_kloc(n, g, dev_end, r_end, npts, sigma, seed):
    """k_loc(r) = k - a r^-g with the deviation dev_end at r_end, plus noise."""
    rng = np.random.default_rng(seed)
    k = 2.0 * (n - 1)
    a = dev_end * r_end ** g
    r = np.exp(np.linspace(math.log(r_end * math.exp(-2.5)), math.log(r_end), npts))
    return r, k - a * r ** (-g) + sigma * rng.standard_normal(npts), k


def chk_PX():
    """PX: a synthetic local exponent with a known limit and a power-law approach through the extrapolation recovers the limit."""
    c = CRIT["PX"]
    worst = 0.0
    lines = []
    for n, g, dev, r_end in c["cases"]:
        for seed in c["seeds"]:
            r, kl, k = planted_kloc(n, g, dev, r_end, c["npts"], c["sigma"], seed)
            ki, gg, rms = k_inf(r, kl)
            worst = max(worst, abs(ki / k - 1))
            lines.append((n, g, dev, seed, ki, gg))
    for n, g, dev, seed, ki, gg in lines[::len(c["seeds"])]:
        print(f"    PX n = {n}: planted g = {g}, deviation at r_end {dev:g}; recovered k_inf = {ki:.6f} (k = {2.0 * (n - 1)}), g = {gg:.2f}")
    check("PX", "planted: a synthetic k_loc(r) = k - a r^-g with noise through the extrapolation recovers k", worst <= c["rel"],
          f"worst relative error {worst:.2e} (tol {c['rel']:g}) over {len(c['cases'])} cases x {len(c['seeds'])} seeds, window = last e-fold of r")


# --------------------------------------------------------------------------------------------------------------------
# the spatial test: r^n = P(x+) + P(x-) and phi = A arcsin((P+ - P-)/(P+ + P-)), P(x) = r_diag(x)^n / 2 from the run's own diagonal
# --------------------------------------------------------------------------------------------------------------------
def spatial_dev(xd, rd, xcol, xm, r_col, phi_col, n, A, xmin):
    """Largest relative error of r^n and largest absolute error of phi against the closed form on x- in [xmin, x+ - 2 node]; P from the diagonal."""
    def interp(x):                              # linear, extended linearly beyond the last cell centre (the column x+ = X is a node past it)
        x = np.asarray(x, float)
        y = np.interp(x, xd, rd)
        return y + np.where(x > xd[-1], (rd[-1] - rd[-2]) / (xd[-1] - xd[-2]) * (x - xd[-1]), 0.0)
    P = lambda x: 0.5 * interp(x) ** n
    sel = (xm >= xmin) & (xm <= xcol - 3 * (xm[1] - xm[0]))
    Pp, Pm = P(xcol), P(xm[sel])
    rn_pred = (Pp + Pm)
    xi = (Pp - Pm) / (Pp + Pm)
    A_ = A * (1.01 if mut("A_scale") else 1.0)
    if mut("exp_n"):
        rn_pred = (Pp + Pm) ** (n / (n - 1.0))
    e_r = float(np.max(np.abs(r_col[sel] ** n / rn_pred - 1)))
    e_p = float(np.max(np.abs(phi_col[sel] - A_ * np.arcsin(xi))))
    return e_r, e_p


def run_spatial(R, n, b=0.0):
    """Spatial deviations on the stored columns of a run: list of (x+, e_r, e_phi)."""
    N, h = R["N"], R["h"]
    xd = (np.arange(N) + 0.5) * h
    rd = np.where(np.isfinite(R["r"]), R["r"], np.nan)
    good = np.isfinite(rd)
    out = []
    for i0 in sorted(R["cols"]):
        lo, arr = R["cols"][i0]
        ii = i0 - lo
        xm = np.arange(N + 1) * h
        if i0 * h < CRIT["sp"]["xmin"] + 10.0:        # the column must hold a range of x- inside the anti-trapped region
            continue
        e_r, e_p = spatial_dev(xd[good], rd[good], i0 * h, xm, arr[0, ii], arr[2, ii], n, A_jump(n, b), CRIT["sp"]["xmin"])
        out.append((i0 * h, e_r, e_p))
    return out


def chk_PS():
    """PS: the spatial test on the closed form itself (exact, it must pass) and with the amplitude off by 2 % (it must fail)."""
    c = CRIT["PS"]
    n = c["n"]
    h = 0.01
    N = int(c["xmax"] / h)
    xd = (np.arange(N) + 0.5) * h
    Pf = lambda x: 0.5 * (1.0 + 3.0 * x ** 2)
    rd = (2 * Pf(xd)) ** (1.0 / n)
    xm = np.arange(N + 1) * h
    A = math.sqrt(A2_pred(n))
    out = {}
    for lab, Ause in (("exact", A), ("pert", A * c["A_pert"])):
        xcol = c["xmax"]
        Pp, Pm = Pf(xcol), Pf(xm)
        r_col = (Pp + Pm) ** (1.0 / n)
        phi_col = Ause * np.arcsin((Pp - Pm) / (Pp + Pm))
        out[lab] = spatial_dev(xd, rd, xcol, xm, r_col, phi_col, n, A, CRIT["sp"]["xmin"])
    ok = (out["exact"][0] <= CRIT["sp"]["r_rel"] and out["exact"][1] <= CRIT["sp"]["phi_abs"]
          and out["pert"][1] > CRIT["sp"]["phi_abs"])
    check("PS", "planted: the spatial test passes on the closed form (n = 3, P(x) = (1 + 3x^2)/2) and fails when the amplitude is 2 % off",
          ok, f"exact: r^n {out['exact'][0]:.1e}, phi {out['exact'][1]:.1e}; amplitude +2 %: phi {out['pert'][1]:.1e}; bars {CRIT['sp']['r_rel']:g}, {CRIT['sp']['phi_abs']:g}")


# --------------------------------------------------------------------------------------------------------------------
# the runs
# --------------------------------------------------------------------------------------------------------------------
def jump_data(n, h, X, b):
    """Family S (ghost phi-Gaussian, amplitude A = CRIT jump A, centre 3, width 0.5) plus a step b in the end value of the field on each ray (odd mirror)."""
    c = CRIT["jump"]
    N = int(round(X / h))
    pulses = [dict(field=0, ray="+", A=c["A"], xc=ex.PULSE["xc"], w=ex.PULSE["w"]),
              dict(field=0, ray="+", A=b / (c["w"] * math.sqrt(math.pi)), xc=c["xc"], w=c["w"], kind="dphi")]
    return dn.ellis_rays(n, N, h, pulses=pulses, symmetric=True)


def run_name(D, member, h, X, b=None):
    if b is not None:
        return f"D{D}_S_b{b:+g}_h{h:g}_X{X:g}"
    return f"D{D}_{'tuned' if member == 'tuned' else 'A' + str(member)}_h{h:g}_X{X:g}"


def job(argv):
    """One streamed run of family S with expansion.py's code. A = from:PATH takes the member (ghost amplitude) of a stored run. The log gets a line per
    tenth of the rows, and the diagonal series so far is saved beside the result (partial.npz) at each."""
    a = dict(x.lstrip("-").split("=", 1) for x in argv)
    D, h, X, out = int(a["D"]), float(a["h"]), float(a["X"]), a["out"]
    n = D - 2
    os.makedirs(out, exist_ok=True)
    b = float(a["b"]) if "b" in a else None           # b: the step of the jump test (A is then the S pulse's amplitude)
    if a["A"].startswith("from:"):
        A = float(ex.load_run(a["A"][5:])["meta"]["A"])
        member = "tuned"
    else:
        A = float(a["A"])
        member = a["A"]
    name = run_name(D, member, h, X, b)
    t0 = time.time()
    partial = os.path.join(out, name + ".partial.npz")

    def log(msg):
        print(f"{time.strftime('%H:%M:%S')} {name}: {msg} ({time.time() - t0:.0f} s)", flush=True)
        try:
            loc = sys._getframe(1).f_locals          # the row sweep's own arrays, read only
            R = loc["R"]
            np.savez_compressed(partial, **{k: R[k] for k in ("r", "f", "rp", "rm", "alive_diag")}, phip=R["phip"], phim=R["phim"], h=h, A=A)
        except Exception as e:                       # a snapshot must never stop the run
            print(f"snapshot failed: {e!r}", flush=True)

    d = ex.s_data(n, h, X, A) if b is None else jump_data(n, h, X, b)
    log(f"start, N = {d.N}, A = {A!r}, b = {b}")
    R = ex.stream(d, tol=ex.CRIT["trap_tol"], log=log)
    secs = time.time() - t0
    meta = dict(kind="S" if b is None else "S+step", n=n, h=h, X=X, member=member, A=A, b=b, secs=secs, mut=None, dn_mut=dn.MUT, script="throat_scaling.py job")
    ex.save_run(os.path.join(out, name + ".npz"), R, meta)
    print(f"{name}: {secs:.1f} s, N = {R['N']}, A = {A!r}, first dead cell {R.get('first_dead')}", flush=True)


# --------------------------------------------------------------------------------------------------------------------
# analysis of the stored runs
# --------------------------------------------------------------------------------------------------------------------
def summarize(path, tag=None, rc=None):
    """One run: the estimators and the spatial test. A dict."""
    R, M = throat(path)
    n = int(R["n"])
    tau, r, H, rho = M["tau"], M["r"], M["H"], M["rho"]
    if len(r) < 100 or r[-1] < 20.0:
        raise ValueError(f"the diagonal series does not reach r = 20 (r_end = {r[-1] if len(r) else float('nan'):.3g}, {len(r)} cells): the run does not expand")
    kl = local_exponent(tau, H)
    kf, ts = fit_k(tau, r, H)
    ki, g, rms = k_inf(r, kl)
    ki2, _, _ = k_inf(r, kl, CRIT["nfold_ws"])
    rc = min(CRIT["c_radius"], 0.8 * r[-1]) if rc is None else rc
    c, rr = c_at(r, rho, H, rc)
    K = H ** 2 + 1 / r ** 2
    s = dict(path=path, name=os.path.basename(path).replace(".npz", ""), n=n, D=n + 2, h=float(R["h"]), X=float(R["meta"].get("X", R["N"] * R["h"])), N=int(R["N"]),
             A=float(R["meta"].get("A", float("nan"))), b=R["meta"].get("b"), r_end=float(r[-1]), tau_end=float(tau[-1]), kfit=kf, ts=ts, kinf=ki, g=g, rms=rms, kinf_ws=ki2,
             c=c, r_c=rr, kl_end=float(kl[-3]), KT=float(K[-1] * (ts - tau[-1]) ** 2), R=R, M=M, kloc=kl, spatial=run_spatial(R, n, R["meta"].get("b") or 0.0))
    return s


def table_line(s):
    return (f"{s['name']:30s} D={s['D']} h={s['h']:<9.5g} N={s['N']:<6d} r_end={s['r_end']:7.1f} tau_s-tau_end={s['ts'] - s['tau_end']:.3f}  k_fit={s['kfit']:.5f}  tau_s={s['ts']:.5f}  "
            f"k_loc(end)={s['kl_end']:.5f}  k_inf={s['kinf']:.5f} (g={s['g']:.2f}, rms {s['rms']:.1e})  c(r={s['r_c']:.0f})={s['c']:.5f}  K(tau_s-tau)^2={s['KT']:.5f}")


def approach_table(s):
    r, M = s["M"]["r"], s["M"]
    n = s["n"]
    cp = c_pred(n)
    out = []
    for r0 in (10.0, 20.0, 50.0, 100.0):
        if r0 <= r[-3]:
            j = int(np.searchsorted(r, r0))
            out.append(f"r = {r[j]:.0f}: k_loc {s['kloc'][j]:.5f}, c/c_pred - 1 {(-M['rho'][j] / M['H'][j] ** 2) / cp - 1:+.2e}")
    return "; ".join(out)


def analyze(runs):
    """runs: dict group -> list of summaries (gated). Groups: 'D4' (#213 h = 0.01), 'D5' (#213 finer h), 'D6' (h = 0.003), 'D6c' (h = 0.006), 'X800'."""
    eps_k = CRIT["kc_rel"]
    for grp, n in (("D4", 2), ("D5", 3)):
        S = runs.get(grp, [])
        if not S:
            for cid in ("K%dF" % (n + 2), "K%dX" % (n + 2), "C%d" % (n + 2)):
                check(cid, "no runs", None)
            continue
        D = n + 2
        kp, cp = k_pred(n), c_pred(n)
        dk = [abs(s["kfit"] / kp - 1) for s in S]
        dx = [abs(s["kinf"] / kp - 1) for s in S]
        dc = [abs(s["c"] / cp - 1) for s in S]
        check(f"K{D}F", f"D = {D}: #213's fitted k (last e-fold of r) within {eps_k:g} (relative) of k = 2(n-1) = {kp:g} for every stored member",
              max(dk) <= eps_k, f"worst {max(dk):.2e}; " + ", ".join(f"{s['name']} {s['kfit']:.5f}" for s in S))
        check(f"K{D}X", f"D = {D}: the local exponent extrapolated to r -> infinity (last e-fold, power-law approach) within {eps_k:g} (relative) of {kp:g}",
              max(dx) <= eps_k, f"worst {max(dx):.2e}; " + ", ".join(f"{s['name']} {s['kinf']:.5f}" for s in S))
        check(f"C{D}", f"D = {D}: |rho|/H^2 at r = 50a within {eps_k:g} (relative) of c = n^3/(2(n-1)) = {cp:.6f}",
              max(dc) <= eps_k, f"worst {max(dc):.2e}; " + ", ".join(f"{s['name']} {s['c']:.5f}" for s in S))
    S = runs.get("D6", [])
    if S:
        kp, cp = k_pred(4), c_pred(4)
        dk = [abs(s["kfit"] - kp) for s in S]
        dx = [abs(s["kinf"] - kp) for s in S]
        dc = [abs(s["c"] / cp - 1) for s in S]
        check("K6F", f"D = 6 at h <= 0.003: the fitted k (last e-fold of r) within {CRIT['k6_abs']:g} of the prediction {kp:g}", max(dk) <= CRIT["k6_abs"],
              f"worst {max(dk):.2e}; " + ", ".join(f"{s['name']} {s['kfit']:.5f}" for s in S))
        check("K6X", f"D = 6: the extrapolated local exponent within {CRIT['k6_abs']:g} of {kp:g}", max(dx) <= CRIT["k6_abs"],
              f"worst {max(dx):.2e}; " + ", ".join(f"{s['name']} {s['kinf']:.5f}" for s in S))
        check("C6", f"D = 6: |rho|/H^2 at r = 50a within {eps_k:g} (relative) of c = {cp:.6f}", max(dc) <= eps_k,
              f"worst {max(dc):.2e}; " + ", ".join(f"{s['name']} {s['c']:.5f}" for s in S))
    else:
        for cid in ("K6F", "K6X", "C6"):
            check(cid, "no runs", None)
    # resolution: D = 6, h = 0.006 against h = 0.003 for the same member
    fine = {s["A"]: s for s in runs.get("D6", [])}
    coarse = runs.get("D6c", [])
    if coarse and fine:
        c = CRIT["res6"]
        worst_H, worst_k = 0.0, 0.0
        for sc in coarse:
            sf = fine.get(sc["A"])
            if sf is None:
                continue
            for r0 in c["r_eval"]:
                Hc = float(np.interp(r0, sc["M"]["r"], sc["M"]["H"]))
                Hf = float(np.interp(r0, sf["M"]["r"], sf["M"]["H"]))
                worst_H = max(worst_H, abs(Hc / Hf - 1))
            worst_k = max(worst_k, abs(sc["kfit"] - sf["kfit"]))
        check("RES6", f"D = 6: H at fixed r (r = 10, 25, 50) within {c['H_rel']:g} (relative) and the fitted k within {c['k_abs']:g} between h = 0.006 and h = 0.003",
              worst_H <= c["H_rel"] and worst_k <= c["k_abs"], f"worst H {worst_H:.2e}, worst k {worst_k:.2e}")
    else:
        check("RES6", "no runs", None)
    # window stability of the extrapolation, on every gated run
    allS = [s for g in ("D4", "D5", "D6", "X800", "JT") for s in runs.get(g, [])]
    if allS:
        ws = max(abs(s["kinf"] - s["kinf_ws"]) / s["kinf"] for s in allS)
        check("WS", f"the extrapolated exponent does not depend on the window (last 1.0 against 1.5 e-folds of r) beyond {CRIT['ws_rel']:g} (relative)", ws <= CRIT["ws_rel"],
              f"worst {ws:.2e}")
    # the closer approach to tau_s at D = 4
    S = runs.get("X800", [])
    if S:
        s = S[0]
        x = CRIT["x800"]
        check("T800", f"D = 4 tuned at X = 800: tau_s from the fit within {x['ts_tol']:g} of {x['ts_ref']}", abs(s["ts"] - x["ts_ref"]) <= x["ts_tol"],
              f"tau_s = {s['ts']:.5f}, last cell at tau_s - tau = {s['ts'] - s['tau_end']:.4f}")
        check("K800", f"D = 4 tuned at X = 800: k = 2.000 +- {x['k_tol']:g}", abs(s["kfit"] - 2.0) <= x["k_tol"], f"k = {s['kfit']:.5f}, k_loc(end) {s['kl_end']:.5f}")
        check("KT800", f"D = 4 tuned at X = 800: K (tau_s - tau)^2 stays at 4.000 within {x['KT_rel']:g} (relative) at the last cell",
              abs(s["KT"] / 4.0 - 1) <= x["KT_rel"], f"K (tau_s - tau)^2 = {s['KT']:.5f}")
    else:
        for cid in ("T800", "K800", "KT800"):
            check(cid, "no runs", None)
    # the jump test: the exponent and c follow the end value of the field through the closed form
    S = runs.get("JT", [])
    if S:
        tol = CRIT["jump"]["tol"]
        bad_k = [(s["name"], s["kinf"], k_jump(s["n"], s["b"])) for s in S if abs(s["kinf"] / k_jump(s["n"], s["b"]) - 1) > tol]
        bad_c = [(s["name"], s["c"], c_jump(s["n"], s["b"])) for s in S if abs(s["c"] / c_jump(s["n"], s["b"]) - 1) > tol]
        wk = max(abs(s["kinf"] / k_jump(s["n"], s["b"]) - 1) for s in S)
        wc = max(abs(s["c"] / c_jump(s["n"], s["b"]) - 1) for s in S)
        check("JTK", f"jump test: the extrapolated exponent of data whose field end value is raised by b is k = 2/(nA^2 - n - 1), A = (2/pi)(phi_inf + b), within {tol:g} (relative)",
              not bad_k, f"worst {wk:.2e} over {len(S)} runs (D = 4, 5, 6)")
        check("JTC", f"jump test: |rho|/H^2 at r = 50a (or 0.8 r_end) is n^2 A^2/2 with the same A, within {tol:g} (relative)", not bad_c, f"worst {wc:.2e} over {len(S)} runs")
    else:
        for cid in ("JTK", "JTC"):
            check(cid, "no runs", None)
    # the spatial test, every gated run
    allS = [s for g in ("D4", "D5", "D6", "D6c", "X800", "JT") for s in runs.get(g, [])]
    if allS:
        sp = CRIT["sp"]
        er = max(e for s in allS for (_, e, _) in s["spatial"])
        ep = max(e for s in allS for (_, _, e) in s["spatial"])
        check("SP1", f"r^n = P(x+) + P(x-) with P = r_diag^n/2 on the stored columns x+ = X/4, X/2, X at x- >= {sp['xmin']:g}, to {sp['r_rel']:g} (relative)", er <= sp["r_rel"],
              f"worst {er:.2e} over {len(allS)} runs")
        check("SP2", f"phi = A arcsin((P+ - P-)/(P+ + P-)) with A = sqrt(n/(n-1)) on the same columns, to {sp['phi_abs']:g} (absolute)", ep <= sp["phi_abs"],
              f"worst {ep:.2e} over {len(allS)} runs")
    else:
        for cid in ("SP1", "SP2"):
            check(cid, "no runs", None)


def reports(runs):
    print("\n--- runs: the throat series and its limits ---")
    for g in ("D4", "D5", "D6", "D6c", "X800"):
        for s in runs.get(g, []):
            print("  " + table_line(s))
            print(f"      approach: {approach_table(s)}")
            print("      spatial (x+, max |r^n/(P+ + P-) - 1|, max |phi - A asin xi|): " + "; ".join(f"{x:.1f}, {er:.1e}, {ep:.1e}" for x, er, ep in s["spatial"]))
    if runs.get("JT"):
        print("\n--- the jump test: A = (2/pi)(phi_inf + b); k = 2/(nA^2 - n - 1), c = n^2 A^2/2 against the runs (k_inf, c at r = 50a or 0.8 r_end) ---")
        print(f"  {'D':<3s}{'b':>7s}{'A':>9s}{'k_pred':>9s}{'k_inf':>9s}{'k_fit':>9s}{'c_pred':>9s}{'c':>9s}{'r_end':>7s}")
        for s in runs["JT"]:
            print(f"  {s['D']:<3d}{s['b']:>+7.3f}{A_jump(s['n'], s['b']):>9.5f}{k_jump(s['n'], s['b']):>9.4f}{s['kinf']:>9.4f}{s['kfit']:>9.4f}{c_jump(s['n'], s['b']):>9.4f}{s['c']:>9.4f}{s['r_end']:>7.1f}")
    print("\n--- the closed-form values ---")
    print(f"  {'D':<3s}{'n':<3s}{'k = 2(n-1)':>12s}{'c':>12s}")
    for n in (2, 3, 4):
        print(f"  {n + 2:<3d}{n:<3d}{2.0 * (n - 1):>12.4f}{c_pred(n):>12.5f}")


def load_runs(i213, i230):
    runs = {}
    if i213:
        d = os.path.join(i213, "runs")
        runs["D4"] = [summarize(os.path.join(d, f)) for f in ("D4_tuned_h0.01.npz", "D4_A0.01_h0.01.npz", "D4_A0.05_h0.01.npz", "D4_A0.2_h0.01.npz") if os.path.exists(os.path.join(d, f))]
        runs["D5"] = [summarize(os.path.join(d, f)) for f in ("D5_tuned_h0.00353553.npz", "D5_A0.05_h0.00353553.npz") if os.path.exists(os.path.join(d, f))]
    if i230:
        d = os.path.join(i230, "runs")
        runs["D6"] = [summarize(os.path.join(d, f"D6_A{a}_h0.003_X200.npz")) for a in CRIT["runs"]["D6"]["members"] if os.path.exists(os.path.join(d, f"D6_A{a}_h0.003_X200.npz"))]
        runs["D6c"] = [summarize(os.path.join(d, f"D6_A{a}_h0.006_X200.npz")) for a in (0.05,) if os.path.exists(os.path.join(d, f"D6_A{a}_h0.006_X200.npz"))]
        runs["JT"] = []
        for D, m in sorted(CRIT["jump"]["members"].items()):
            for b in m["b"]:
                p = os.path.join(d, run_name(D, None, m["h"], m["X"], b) + ".npz")
                if os.path.exists(p):
                    runs["JT"].append(summarize(p))
        p = os.path.join(d, "D4_tuned_h0.01_X800.npz")
        runs["X800"] = [summarize(p)] if os.path.exists(p) else []
    return runs


def mini_runs():
    """The reduced suite for the mutation runs: D = 4, A = 0.05, h = 0.02, X = 100 (about 15 s)."""
    m = CRIT["mini"]
    d = ex.s_data(m["D"] - 2, m["h"], m["X"], m["A"])
    R = ex.stream(d, tol=ex.CRIT["trap_tol"])
    import tempfile
    p = os.path.join(tempfile.mkdtemp(prefix="ts_mini_"), "mini.npz")
    ex.save_run(p, R, dict(kind="S", n=m["D"] - 2, h=m["h"], X=m["X"], member=str(m["A"]), A=m["A"], mut=MUT, dn_mut=dn.MUT))
    d = jump_data(m["D"] - 2, m["h"], m["X"], m["b"])
    R = ex.stream(d, tol=ex.CRIT["trap_tol"])
    p2 = os.path.join(os.path.dirname(p), "mini_jump.npz")
    ex.save_run(p2, R, dict(kind="S+step", n=m["D"] - 2, h=m["h"], X=m["X"], member="S+step", A=CRIT["jump"]["A"], b=m["b"], mut=MUT, dn_mut=dn.MUT))
    return {"D4": [summarize(p)], "JT": [summarize(p2)]}


def main():
    global MUT
    argv = sys.argv[1:]
    if not argv:
        sys.exit(__doc__)
    cmd, rest = argv[0], argv[1:]
    mutname = i213 = i230 = None
    for x in rest:
        if x.startswith("--mutate="):
            mutname = x.split("=", 1)[1]
        elif x.startswith("--i213="):
            i213 = x.split("=", 1)[1]
        elif x.startswith("--i230="):
            i230 = x.split("=", 1)[1]
    if mutname:
        known = set(MUTATIONS) | {"exp_n", "f_coef", "c_formula", "k_formula", "A_scale", "fit_unit", "extrap_mean", "jump_off"} | set(ex.DN_MUT)
        if mutname not in known:
            sys.exit(f"unknown mutation {mutname}; known: {', '.join(sorted(known))}")
        if mutname in ex.DN_MUT:
            dn.MUT = mutname
        else:
            MUT = mutname
    if cmd == "job":
        job([x for x in rest if not x.startswith("--mutate=")])
        return
    print(f"throat_scaling.py {cmd}: mutation {MUT or dn.MUT}, numba {dn.nb.__version__ if dn.HAVE_NUMBA else 'absent'}, numpy {np.__version__}, {time.strftime('%Y-%m-%d %H:%M:%S')}")
    if cmd == "derive":
        chk_CF1_CF2()
        chk_CF3_CF4()
        bad = [c for c, ok in results if ok is False]
        sys.exit(1 if bad else 0)
    if cmd in ("mutate", "dev"):
        try:
            suite = mini_runs()
        except Exception as e:                       # a mutation can break the build of the reduced suite: that is the check BUILD
            suite = {}
            check("BUILD", "the reduced suite builds and its diagnostics are finite", False, repr(e))
    elif cmd == "analyze":
        suite = load_runs(i213, i230)
    else:
        sys.exit(f"unknown command {cmd}")
    for fn in (chk_CF1_CF2, chk_CF3_CF4, chk_PF, chk_PX, chk_PS):
        try:
            fn()
        except Exception as e:                      # an exception inside a check is that check's failure
            import traceback
            traceback.print_exc()
            check(fn.__name__[4:], "check raised an exception", False, repr(e))
    analyze(suite)
    if cmd == "analyze":
        reports(suite)
    bad = [c for c, ok in results if ok is False]
    n_ok = sum(1 for _, ok in results if ok)
    print(f"\n{n_ok} of {sum(1 for _, ok in results if ok is not None)} checks pass; failing: {bad}")
    if cmd == "mutate":
        want = set(MUTATIONS.get(mutname, []))
        fired = set(bad)
        print(f"mutation {mutname}: declared target {sorted(want)}; fired {sorted(fired)}; " + ("targets all fired" if want <= fired else f"SILENT: {sorted(want - fired)}")
              + ("" if fired <= want else f"; extra: {sorted(fired - want)}"))
    sys.exit(1 if (bad and cmd != "mutate") else 0)


if __name__ == "__main__":
    main()
