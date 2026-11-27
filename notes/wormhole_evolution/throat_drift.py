#!/usr/bin/env python3
"""What drifts the throat radius of family G (#255), and the expanding member A = 0.05 at X = 160.
Uses dn.py, checks.py, threshold.py, expansion.py of #39, #40, #213 and throat_locus.py of #225 unchanged (imported, never modified).

    python3 throat_drift.py bisect  --D=5 --B=0.04 --hk=1 --out=DIR          p*_h(B, h) by bisection of #40's events
    python3 throat_drift.py stall   --D=5 --B=0.04 --hk=1 --delta=1e-10 --out=DIR [--pstar=...]   the pair p*_h -/+ delta kept whole
    python3 throat_drift.py bisect_stall ...                                  both
    python3 throat_drift.py expand  --hk=0 --X=160 --out=DIR                 G, D = 5, A = 0.05 (throat_locus.job_expand)
    python3 throat_drift.py analyze --raw=DIR/runs [--mutate=NAME]            every check and table from the stored runs
    python3 throat_drift.py dev | mutate --mutate=NAME                        the development suite (stored #225 pairs, planted cases)

The criteria are CRIT and MUTATIONS below, as declared in negative-mass-raw I255-dh610/README.md before the runs of record.
Conventions as throat_locus.py: a = 1, 8 pi G = 1, t = (x+ + x-)/sqrt2, the throat locus is the minimum of r on each slice t = const, tau its proper time.
Family G of #40 with the normal pulse amplitude B free: ghost pulse A on '+' plus a normal pulse B on the same ray (centre 3, width 0.5)."""
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import checks       # noqa: E402
import dn           # noqa: E402
import threshold    # noqa: E402
import throat_locus as tl   # noqa: E402
import expansion as xp      # noqa: E402

SQ2 = math.sqrt(2.0)
W_EXACT = checks.W_EXACT
LMUT = None                     # this file's mutation
RAW225 = os.path.expanduser("~/Source/negative-mass-raw/I225-dh610/runs")


def lmut(name):
    return LMUT == name


# --------------------------------------------------------------------------------------------------------------------
# declared criteria
# --------------------------------------------------------------------------------------------------------------------
CRIT = dict(
    Bs=(0.04, 0.1, 0.2), B_ref=0.1, hks=(0, 1, 2),               # D = 5: the three B, #40's three h (index of threshold.CRIT['H'])
    D4=dict(Bs=(0.04, 0.1, 0.2), hks=(0, 1)),                    # D = 4: reported, not gated
    delta=1e-10, delta_ctl=1e-9,                                  # the pair p*_h -/+ delta (the control is #225's delta)
    bis_tol=1e-12, bis_tcap=14.0, bracket=(0.3, 3.0),             # bisection to hi - lo <= tol; the triangle covers t <= bis_tcap + ln(1/(hi - lo))/W (x 1.4 on a retry, 3 times); bracket = (0.3, 3) x (B/0.1)^2 p*_h(0.1, h)
    win=(6.5, 12.0), odd_max=3e-4, vwin=0.5,                      # the fit window in t of rdot(tau), the bar of the pair's odd part, the half window (in tau) of the local slope
    ts=(8.0, 10.0, 12.0),                                         # the times at which r_th - 1 is compared
    pw=dict(lo=1.8, hi=2.4),                                      # the band of the exponent between B = 0.04 and 0.1 (first order would give 1)
    kap=dict(lo=0.75, hi=1.25),                                   # kappa(0.04) / kappa(0.1)
    fin=2,                                                        # the gates PS, PR, PV are evaluated at the finest `fin` h
    cv=dict(order=(1.7, 2.3), two_h=0.01),                        # observed order of r_th(10) - 1 over three h; the two finest within 1 %
    qs=dict(bar=5e-3, win=(8.0, 12.0)),                           # |rho_g / rho_eq - 1| on the locus
    ctl=dict(pts=8, bar=0.0, arr=1e-12),                          # CT1: number of #40's points; bit for bit; CT2, CT3: arrays equal to 1e-12
    dl=dict(t=(7.0, 11.0), bar=1e-5),                             # DL: the tuned mean at delta = 1e-9 and 1e-10 (B = 0.1, h = 0.00707) differ by less than this
    pl1=dict(D=5, c=1.2, hs=(0.02, 0.01, 0.005), t=(0.8, 2.7), bar=2e-4, rdot=5e-4, order=(1.75, 2.25)),     # planted static throat of radius c
    pl2=dict(D=5, vs=(0.13, -0.3), hs=(0.02, 0.01, 0.005), t=(0.8, 2.7), bar=5e-4, rdot=1e-3),               # planted moving static throat (r_th = 1, rdot = 0)
    pl3=dict(inf=(0.026, 0.0, -0.01), kappa=0.48, c=-0.013, noise=1e-6, bar=1e-5, kbar=5e-3),                           # synthetic rdot series: recovered rdot_inf, kappa
    px=dict(hk=(0, 1, 2), alphas=(2.0, 1.0, 3.0)),                                                           # planted B families: the gates pass at 2, fail at 1 and 3
    X=160.0, exp_hk=(0, 1), member="0.05",                        # the expanding member: D = 5, G, A = 0.05, two h
    ft_tsmax=1.2,                                                 # #213's bound on tau_s / tau_end (xp.CRIT['CL']['ft_tsmax'])
    res=dict(H_rel=1e-2, k_rel=1e-2, ts_rel=1e-2),
)

MUTATIONS = {   # name: the checks whose verdict changes against the unmutated suite (computed in the development suite, dev/mut in the raw folder, before they were declared)
    "b_label": ["PR", "PS", "PV", "PX1", "PX2"],                  # B labelled sqrt(B x 0.1): the exponents come out doubled
    "rth_diag": ["CV", "PL2", "QS", "RV"],                        # the diagonal's radius instead of the locus's
    "tau_unit": ["PL3", "RV"],                                    # tau in units of sqrt2 t
    "odd_off": ["DL"],                                            # one side of the pair instead of the mean
    "pstar_off": ["CT2"],                                         # run level: p*_h shifted by 3e-9 in the control pair
    "ts_bound": ["KX2"],                                          # the bound on tau_s / tau_end relaxed to 1.32
    "rho_nolapse": ["QS"],                                        # e^f dropped from the stress on the locus
    "cv_swap": ["CV"],                                            # the h in the wrong order in the order estimate
    "pstar_swap": ["PS", "PX2"],                                  # p*_h of B = 0.04 and 0.2 exchanged
    "r_not_r1": ["PR", "PX2"],                                    # r instead of r - 1 in the exponent
    "rdot_swap": ["PV", "PX2"],                                   # the rdot fits of B = 0.04 and 0.2 exchanged
    "ts_scale": ["RES2"],                                         # the coarse run's tau_s scaled by 1.02
    "bis_loose": ["CT4"],                                         # run level: the control bisection stopped at 1e-9
}

results = []


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok))
    tag = "SKIP" if ok is None else ("pass" if ok else "FAIL")
    print(f"[{tag}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def fmt(xs, f="{:.2e}"):
    return "[" + ", ".join(f.format(x) for x in xs) + "]"


def order_of(errs):
    return [math.log2(errs[i] / errs[i + 1]) if errs[i + 1] > 0 and errs[i] > 0 else float("inf") for i in range(len(errs) - 1)]


# --------------------------------------------------------------------------------------------------------------------
# the family with the normal amplitude B free, and the planted static throat of radius c
# --------------------------------------------------------------------------------------------------------------------
def family_B(B, A):
    """(pulses, symmetric, eps) of threshold.family('G', A) with the normal pulse amplitude B instead of threshold.CRIT['G_B']."""
    pu = threshold.CRIT["pulse"]
    return ([dict(field=0, ray="+", A=A, xc=pu["xc"], w=pu["w"]), dict(field=1, ray="+", A=B, xc=pu["xc"], w=pu["w"])], False, (-1.0, 1.0))


def make_data_B(n, h, N, A, B):
    pulses, sym, eps = family_B(B, A)
    return dn.ellis_rays(n, N, h, pulses=pulses, symmetric=sym, eps=eps)


def static_rays(n, N, h, c):
    """The static Ellis wormhole of throat radius c as characteristic data: r(l) = c r_1(l/c), phi(l) = phi_1(l/c), f = 0 (the equations are invariant under
    l -> c l, r -> c r with phi fixed). The ray constraint r'' = -(S/n) r is solved by RK4 as dn.ellis_rays does; c = 1 reproduces dn.ellis_rays."""
    sol, q = dn.ellis_dense(n, N * h / SQ2 + 1.0)
    qc = q * c ** (n - 1)
    rays = []
    nodes = np.arange(N + 1) * h
    for ray in ("+", "-"):
        sg = 1.0 if ray == "+" else -1.0

        # d phi / ds = sg * q c^(n-1) / r^n / sqrt2 with r = c r_1(l/c): written through r_1 so that c = 1 is dn.ellis_rays's
        def phi_s(s):
            r1 = sol(s / SQ2 / c)[0]
            return sg * qc / (c * r1) ** n / SQ2

        S0, Sh, S1 = (-phi_s(nodes[:-1]) ** 2, -phi_s(nodes[:-1] + 0.5 * h) ** 2, -phi_s(nodes[1:]) ** 2)
        r = np.empty(N + 1)
        v = np.empty(N + 1)
        r[0], v[0] = c, 0.0
        for i in range(N):
            y = (r[i], v[i])

            def d_(y, S):
                return (y[1], -(S / n) * y[0])
            k1 = d_(y, S0[i])
            k2 = d_((y[0] + 0.5 * h * k1[0], y[1] + 0.5 * h * k1[1]), Sh[i])
            k3 = d_((y[0] + 0.5 * h * k2[0], y[1] + 0.5 * h * k2[1]), Sh[i])
            k4 = d_((y[0] + h * k3[0], y[1] + h * k3[1]), S1[i])
            r[i + 1] = y[0] + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            v[i + 1] = y[1] + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        U = np.zeros((3, N + 1))
        U[0] = r
        U[2] = sg * sol(nodes / SQ2 / c)[2]
        rays.append(U)
    return dn.Data(n, [-1.0], N, h, rays[0], rays[1], name="ellis_static_scale")


# --------------------------------------------------------------------------------------------------------------------
# jobs: bisection of #40's events at (B, h), the pair p*_h -/+ delta kept whole
# --------------------------------------------------------------------------------------------------------------------
def fate_run(n, B, h, A, smax):
    """One run to #40's event (the sweep reproduces threshold.point's event and time bit for bit): dict(A, fate, T)."""
    d = make_data_B(n, h, smax, A, B)
    L = tl.sweep(d, smax=smax, stop=True, need_loc=False)
    key = next(iter(L["ev40"]), None)
    return dict(A=A, fate=(-1 if key == "collapse" else 1 if key else 0), T=tl.ev40_time(L, key) if key else None, cell=list(L["ev40"][key]) if key else None)


def bracket_of(D, B, hk):
    """(lo, hi) = (0.3, 3) x (B/B_ref)^2 p*_h(B_ref, h) at the same D and h (p*_h of #40's raw)."""
    ps = tl.pstar_raw40(D, "G", hk)[0]
    c = (B / CRIT["B_ref"]) ** 2 * ps
    return CRIT["bracket"][0] * c, CRIT["bracket"][1] * c


def job_bisect(D, B, hk, out):
    """Bisection of #40's events between lo (collapse) and hi (expansion) to hi - lo <= bis_tol. The bracket is declared (bracket_of); if an end does not have its fate it is moved by a
    factor 2 (lo / 2, hi x 2), at most 4 times, and the moves are in summary.json."""
    n = D - 2
    h = threshold.CRIT["H"][n][hk]
    lo, hi = bracket_of(D, B, hk)
    d = os.path.join(out, f"D{D}_B{B:g}_h{h:g}")
    os.makedirs(d, exist_ok=True)
    W = W_EXACT[n]
    sink = tl.Sink(os.path.join(d, "points.jsonl"))
    t0 = time.time()
    nrun = [0]
    br = dict(lo=lo, hi=hi)
    tol = 1e-9 if lmut("bis_loose") else CRIT["bis_tol"]            # bis_loose: the bisection stopped at 1e-9

    def run(p):
        cap = CRIT["bis_tcap"] + math.log(1.0 / max(br["hi"] - br["lo"], 1e-15)) / W
        for k in range(4):
            smax = tl.tcap_smax(h, cap)
            r = fate_run(n, B, h, p, smax)
            if r["fate"] != 0:
                break
            cap *= 1.4
        nrun[0] += 1
        if r["fate"] == -1:
            br["lo"] = max(br["lo"], p)
        elif r["fate"] == 1:
            br["hi"] = min(br["hi"], p)
        print(f"  bisect run {nrun[0]}: A = {p:.15e} fate {r['fate']} T {r['T']} (cap {cap:.1f}, {time.time() - t0:.0f} s)", flush=True)
        return r
    moves = []
    for k in range(5):
        flo, fhi = run(lo), run(hi)
        sink(dict(phase="end", lo=lo, hi=hi, fate_lo=flo["fate"], fate_hi=fhi["fate"]))
        if flo["fate"] == -1 and fhi["fate"] == 1:
            break
        if k == 4:
            break
        if flo["fate"] != -1:
            lo /= 2
        if fhi["fate"] != 1:
            hi *= 2
        moves.append([lo, hi])
    try:
        l, u, stopped = threshold.bisect(run, lo, hi, tol, sink=sink)
        err = None
    except RuntimeError as exc:
        l, u, stopped, err = lo, hi, "bracket", str(exc)
    s = dict(D=D, B=B, h=h, hk=hk, lo=l, hi=u, pstar=0.5 * (l + u), stopped=stopped, error=err, bracket=list(bracket_of(D, B, hk)), moves=moves, tol=tol, runs=nrun[0], secs=time.time() - t0)
    json.dump(s, open(os.path.join(d, "summary.json"), "w"), indent=1)
    print(f"bisect D={D} B={B:g} h={h:g}: p* = {s['pstar']:.15e}, bracket {u - l:.2e}, {nrun[0]} runs, {s['secs']:.0f} s, stopped {stopped}, error {err}, bracket moves {moves}", flush=True)
    return err is None and stopped is None


def pstar_of(D, B, hk, bis_dir=None):
    """p*_h at (D, B, h): #40's for B = B_ref, else the bisection's summary.json under bis_dir."""
    if B == CRIT["B_ref"]:
        return tl.pstar_raw40(D, "G", hk)[0]
    h = threshold.CRIT["H"][D - 2][hk]
    s = json.load(open(os.path.join(bis_dir, f"D{D}_B{B:g}_h{h:g}", "summary.json")))
    return s["pstar"]


def stall_name(D, B, h, delta, sd):
    return f"stall_D{D}_B{B:g}_h{h:g}_d{delta:g}_side{sd}.npz"


def job_stall(D, B, hk, out, delta, pstar=None, bis_dir=None):
    """The pair p*_h -/+ delta of family G with the normal amplitude B, kept whole as throat_locus.job_stall does: the locus series and the profile of every slice t = 2, 3, ..."""
    n = D - 2
    h = threshold.CRIT["H"][n][hk]
    ps = pstar if pstar is not None else pstar_of(D, B, hk, bis_dir)
    smax = tl.tcap_smax(h)
    os.makedirs(out, exist_ok=True)
    tk = np.arange(tl.STALL["t_lo"], tl.TCAP - 3.0, tl.STALL["dt"])
    prof_s = sorted({int(round(t * SQ2 / h)) for t in tk})
    for sd in (-1, 1):
        A = ps + sd * delta
        t0 = time.time()
        d = make_data_B(n, h, smax, A, B)
        L = tl.sweep(d, smax=smax, stop=True, col_levels=tl.LEVELS["col"], exp_levels=tl.LEVELS["exp"], prof_s=prof_s, log=lambda m: print(m, flush=True))
        tl.locus_series(L)
        meta = dict(D=D, B=B, family="G", h=h, A=A, side=sd, delta=delta, pstar=ps, smax=smax, secs=time.time() - t0, rec=tl.record(L))
        tl.save_run(os.path.join(out, stall_name(D, B, h, delta, sd)), L, meta)
        print(f"stall D={D} B={B:g} h={h:g} delta={delta:g} side={sd}: {time.time() - t0:.0f} s, slices {len(L['s'])}, cross {L['cross']}", flush=True)
    return True


# --------------------------------------------------------------------------------------------------------------------
# the tuned throat series from a pair, rdot, the fit, the exponents
# --------------------------------------------------------------------------------------------------------------------
def pair_series(La, Lb):
    """The pair's mean on the slices both runs have: the tuned member to second order in the offset; odd = (r_+ - r_-)/2 is the unstable mode's part."""
    common = np.intersect1d(La["s"], Lb["s"])
    ia, ib = np.searchsorted(La["s"], common), np.searchsorted(Lb["s"], common)

    def mean(k):
        return 0.5 * (La[k][ia] + Lb[k][ib])
    if lmut("odd_off"):                                  # odd_off: one side of the pair instead of the mean
        def mean(k):                                     # noqa: F811
            return La[k][ia]
    out = dict(s=common, t=mean("t"), tau=mean("tau"), r=mean("r"), l=mean("l"), f=mean("f"), v=mean("v"), odd=0.5 * (Lb["r"][ib] - La["r"][ia]),
               php=0.5 * (La["php"][ia] + Lb["php"][ib]), phm=0.5 * (La["phm"][ia] + Lb["phm"][ib]), rp=mean("rp"), rm=mean("rm"))
    if lmut("rth_diag"):                                 # rth_diag: the diagonal's radius instead of the locus's
        out["r"] = 0.5 * (La["diag"]["r"][common // 2] + Lb["diag"]["r"][common // 2])
    if lmut("tau_unit"):                                 # tau_unit: tau in units of sqrt2 t
        out["tau"] = out["tau"] * SQ2
    return out


def rdot_local(tau, r, win=None):
    """dr/dtau by a local linear fit over +- win in tau, on every k/10-th slice (k = the half window in slices); nan elsewhere."""
    win = CRIT["vwin"] if win is None else win
    out = np.full(len(r), np.nan)
    if len(r) < 8:
        return out
    dt = np.median(np.diff(tau))
    k = max(2, int(round(win / dt)))
    for i in range(k, len(r) - k, max(1, k // 10)):
        out[i] = np.polyfit(tau[i - k:i + k + 1], r[i - k:i + k + 1], 1)[0]
    return out


def fit_rdot(t, tau, r, odd, win=None, odd_max=None):
    """rdot(tau) = a + c exp(-kappa (tau - tau_0)) over t in `win` where |odd| <= odd_max. Returns dict(inf, c, kappa, se_inf, se_kappa, rms, n, t0, t1)."""
    from scipy.optimize import curve_fit
    win = CRIT["win"] if win is None else win
    odd_max = CRIT["odd_max"] if odd_max is None else odd_max
    rd = rdot_local(tau, r)
    ok = np.isfinite(rd) & (t >= win[0]) & (t <= win[1]) & (np.abs(odd) <= odd_max)
    if ok.sum() < 20:
        return dict(inf=float("nan"), c=float("nan"), kappa=float("nan"), se_inf=float("nan"), se_kappa=float("nan"), rms=float("nan"), n=int(ok.sum()), t0=float("nan"), t1=float("nan"))
    x, y = tau[ok], rd[ok]
    f = lambda z, a, c, kap: a + c * np.exp(-kap * (z - x[0]))      # noqa: E731
    p, cov = curve_fit(f, x, y, p0=[y[-1], y[0] - y[-1], 0.5], maxfev=20000)
    se = np.sqrt(np.diag(cov))
    return dict(inf=float(p[0]), c=float(p[1]), kappa=float(p[2]), se_inf=float(se[0]), se_kappa=float(se[2]), rms=float(np.sqrt(np.mean((f(x, *p) - y) ** 2))), n=int(ok.sum()),
                t0=float(t[ok][0]), t1=float(t[ok][-1]))


def expo(a, b, Ba, Bb):
    """Exponent alpha of a ~ B^alpha between B = Bb and B = Ba (a, b the values there)."""
    return math.log(a / b) / math.log(Ba / Bb) if a > 0 and b > 0 else float("nan")


def bval(B):
    """The B of a record as the analysis labels it (b_label: a wrong label, B^(1/2) scaling)."""
    return math.sqrt(B * CRIT["B_ref"]) if lmut("b_label") else B


def Rt(S, ts=None):
    """r_th(t) - 1 of a series at the times ts (linear interpolation)."""
    ts = CRIT["ts"] if ts is None else ts
    return [float(np.interp(t, S["t"], S["r"]) - (0.0 if lmut("r_not_r1") else 1.0)) if S["t"][0] <= t <= S["t"][-1] else float("nan") for t in ts]       # r_not_r1: r instead of r - 1


# --------------------------------------------------------------------------------------------------------------------
# the pair's diagnostics on the locus (reported): the ghost's energy density against the static throat's, the charge flux, the lapse
# --------------------------------------------------------------------------------------------------------------------
def locus_diag(S, n):
    """On the pair's mean: rho_g in the locus frame (the ghost, eps = -1), the equilibrium rho_eq = -n(n-1)/(2 r^2) of a static throat of that radius, the flux rate of the ghost's charge per
    unit t through the locus, q_t = r^n (phi_- (1 - v) - phi_+ (1 + v))/sqrt2, and N_th = dtau/dt on the locus."""
    v, f, r = S["v"], S["f"], S["r"]
    Tpp, Tmm = -S["php"][:, 0] ** 2, -S["phm"][:, 0] ** 2
    ef = np.ones_like(f) if lmut("rho_nolapse") else np.exp(f)           # rho_nolapse: e^f dropped from the stress components
    Ttt, Ttl = ef * (Tpp + Tmm) / 2, ef * (Tpp - Tmm) / 2
    rho = (1 + v ** 2) * Ttt + 2 * v * Ttl
    rho = rho / (1 - v ** 2)
    rho_eq = -n * (n - 1) / (2 * r ** 2)
    qt = r ** n * (S["phm"][:, 0] * (1 - v) - S["php"][:, 0] * (1 + v)) / SQ2
    return dict(ratio=rho / rho_eq, qt=qt, Nth=np.exp(-f / 2) * np.sqrt(1 - v ** 2))


def pulse_integrals(n):
    """I_2 = integral r^n (dg/ds)^2 ds and I_1 = integral r^n dg/ds ds of the unit Gaussian pulse (centre 3, width 0.5) on the initial ray, r from the static wormhole."""
    pu = threshold.CRIT["pulse"]
    sol, q = dn.ellis_dense(n, 40.0)
    s = np.linspace(0.0, 40.0, 40001)
    r = sol(s / SQ2)[0]
    g = np.exp(-(s - pu["xc"]) ** 2 / (2 * pu["w"] ** 2))
    gs = -(s - pu["xc"]) / pu["w"] ** 2 * g
    return float(np.trapezoid(r ** n * gs ** 2, s)), float(np.trapezoid(r ** n * gs, s))


# --------------------------------------------------------------------------------------------------------------------
# the predictions (PS, PR, PV) as a function of a family of tuned series, and the second-order convergence (CV)
# --------------------------------------------------------------------------------------------------------------------
def eval_predictions(T, P, hks, tag="", out=print):
    """T[(B, hk)] = pair_series, P[(B, hk)] = p*_h. At each h in hks (the finest ones): exponents of p*_h, of r_th - 1 at CRIT['ts'] and of rdot_inf between B = 0.04 and 0.1
    (gated: in CRIT['pw']), kappa(0.04)/kappa(0.1) (gated), and the same between 0.1 and 0.2 (reported). Returns the verdicts and the numbers."""
    Bs = CRIT["Bs"]
    B1, B2, B3 = (bval(b) for b in Bs)
    okS, okR, okV = True, True, True
    rows = []
    fits = {}
    for hk in hks:
        for B in Bs:
            S = T[(B, hk)]
            fits[(B, hk)] = fit_rdot(S["t"], S["tau"], S["r"], S["odd"])
        Pk = {b: P[(Bs[i if not lmut("pstar_swap") else (2 - i if i != 1 else 1)], hk)] for i, b in enumerate(Bs)}       # pstar_swap: p*_h of B = 0.04 and 0.2 exchanged
        a_p1 = expo(Pk[Bs[1]], Pk[Bs[0]], B2, B1)
        a_p2 = expo(Pk[Bs[2]], Pk[Bs[1]], B3, B2)
        R = {B: Rt(T[(B, hk)]) for B in Bs}
        a_R1 = [expo(R[Bs[1]][i], R[Bs[0]][i], B2, B1) for i in range(len(CRIT["ts"]))]
        a_R2 = [expo(R[Bs[2]][i], R[Bs[1]][i], B3, B2) for i in range(len(CRIT["ts"]))]
        f1, f2, f3 = (fits[(B, hk)] for B in Bs)
        if lmut("rdot_swap"):                                    # rdot_swap: the fits of B = 0.04 and 0.2 exchanged
            f1, f3 = f3, f1
        a_v1, a_v2 = expo(f2["inf"], f1["inf"], B2, B1), expo(f3["inf"], f2["inf"], B3, B2)
        kr = f1["kappa"] / f2["kappa"] if f2["kappa"] else float("nan")
        lo, hi = CRIT["pw"]["lo"], CRIT["pw"]["hi"]
        okS &= lo <= a_p1 <= hi
        okR &= all(lo <= a <= hi for a in a_R1)
        okV &= (lo <= a_v1 <= hi) and (CRIT["kap"]["lo"] <= kr <= CRIT["kap"]["hi"])
        rows.append((hk, a_p1, a_p2, a_R1, a_R2, a_v1, a_v2, kr, f1, f2, f3))
        out(f"    {tag} h index {hk}: p*_h {P[(Bs[0], hk)]:.6e}, {P[(Bs[1], hk)]:.6e}, {P[(Bs[2], hk)]:.6e}; exponent of p*_h {a_p1:.3f} (0.04 to 0.1), {a_p2:.3f} (0.1 to 0.2)")
        out(f"       r_th - 1 at t = {CRIT['ts']}: B=0.04 {fmt(R[Bs[0]], '{:.5f}')}, B=0.1 {fmt(R[Bs[1]], '{:.5f}')}, B=0.2 {fmt(R[Bs[2]], '{:.5f}')}; exponents {fmt(a_R1, '{:.3f}')} (0.04 to 0.1), {fmt(a_R2, '{:.3f}')} (0.1 to 0.2)")
        out(f"       rdot_inf {f1['inf']:+.5f} +- {f1['se_inf']:.1e}, {f2['inf']:+.5f} +- {f2['se_inf']:.1e}, {f3['inf']:+.5f} +- {f3['se_inf']:.1e} (rms residual {f1['rms']:.1e}, {f2['rms']:.1e}, {f3['rms']:.1e}; n {f1['n']}, {f2['n']}, {f3['n']}); "
            f"exponents {a_v1:.3f}, {a_v2:.3f}; kappa {f1['kappa']:.4f}, {f2['kappa']:.4f}, {f3['kappa']:.4f} (ratio 0.04 / 0.1 = {kr:.3f})")
    return dict(okS=okS, okR=okR, okV=okV, rows=rows, fits=fits)


def eval_convergence(T, B, hks, out=print, tag=""):
    """r_th(10) - 1 at the three h: observed order and the two finest within CRIT['cv']['two_h']."""
    R = [Rt(T[(B, hk)], (10.0,))[0] for hk in (hks[::-1] if lmut("cv_swap") else hks)]       # cv_swap: the h in the wrong order
    if not all(np.isfinite(R)) or len(R) < 3:
        return False, float("nan"), float("nan")
    e1, e2 = R[0] - R[1], R[1] - R[2]
    q = math.log2(abs(e1 / e2)) if e1 * e2 > 0 else float("nan")
    two = abs(R[2] / R[1] - 1)
    ok = CRIT["cv"]["order"][0] <= q <= CRIT["cv"]["order"][1] and two <= CRIT["cv"]["two_h"]
    out(f"    {tag} B = {B:g}: r_th(10) - 1 at the three h {fmt(R, '{:.6f}')}, observed order {q:.2f}, finest pair differ by {two:.2e}")
    return ok, q, two


def planted_B(Tref, Pref, alpha):
    """A planted family: r_th - 1 and the odd part scaled by (B/B_ref)^alpha, p*_h scaled by the same, from the stored pairs at B_ref (known exponent alpha)."""
    T, P = {}, {}
    for (B0, hk), S in Tref.items():
        for B in CRIT["Bs"]:
            f = (B / CRIT["B_ref"]) ** alpha
            S2 = dict(S)
            S2["r"] = 1.0 + f * (S["r"] - 1.0)
            S2["odd"] = f * S["odd"]
            T[(B, hk)] = S2
            P[(B, hk)] = Pref[(B0, hk)] * f
    return T, P


# --------------------------------------------------------------------------------------------------------------------
# planted evolutions: a static throat of radius c (PL1), a static throat moving in the grid (PL2), a synthetic rdot series (PL3)
# --------------------------------------------------------------------------------------------------------------------
def planted_errors(kind, n, par, hs, t_lo, t_hi):
    """Errors of the locus's r_th against the known value (kind 'scale': c = par; 'boost': 1, the pair-rescaled throat with drift v = par) and of rdot against 0, per h."""
    out = []
    for h in hs:
        smax = int((t_hi + 2) * SQ2 / h)
        d = static_rays(n, smax, h, par) if kind == "scale" else tl.scaled_rays(n, smax, h, tl.lam_of_v(par))
        L = tl.sweep(d, smax=smax, stop=False)
        tl.locus_series(L)
        t, r = L["t"], L["r"]
        if lmut("rth_diag"):                                # rth_diag: the diagonal's radius instead of the locus's
            r = L["diag"]["r"][L["s"] // 2]
        tau = L["tau"] * (SQ2 if lmut("tau_unit") else 1.0)
        rd = rdot_local(tau, r)
        m = (t >= t_lo) & (t <= t_hi)
        ex = par if kind == "scale" else 1.0
        out.append(dict(r=float(np.max(np.abs(r[m] - ex))), rdot=float(np.nanmax(np.abs(rd[m]))) if np.isfinite(rd[m]).any() else float("nan")))
    return out


def chk_PL12():
    c1, c2 = CRIT["pl1"], CRIT["pl2"]
    n = c1["D"] - 2
    for cid, kind, pars, c in (("PL1", "scale", (c1["c"],), c1), ("PL2", "boost", c2["vs"], c2)):
        ok = True
        for par in pars:
            try:
                errs = planted_errors(kind, n, par, c["hs"], *c["t"])
            except Exception as exc:                        # a mutated evolution may leave no locus
                check(cid, "planted case: the locus code raised an exception", False, repr(exc))
                ok = None
                break
            fin = errs[-1]
            ords = order_of([e["r"] for e in errs]) if kind == "scale" else []
            ok &= fin["r"] <= c["bar"] and fin["rdot"] <= c["rdot"]
            if kind == "scale":
                ok &= all(c["order"][0] <= x <= c["order"][1] for x in ords)
            print(f"    {cid} D = {c['D']}, " + (f"static throat of radius {par:g}" if kind == "scale" else f"static throat moving with v = {par:g} (r_th = 1)") + f", h = {c['hs'][-1]:g}, t in {c['t']}: "
                  f"max |r_th - exact| {fin['r']:.1e}, max |rdot| {fin['rdot']:.1e}" + (f"; orders {fmt(ords, '{:.2f}')}" if ords else ""))
        if ok is not None:
            check(cid, ("planted static throat of known radius c = %g" % c1["c"]) if kind == "scale" else "planted static throat moving in the grid (v = %s): the throat radius is 1 and its rate 0 (the drift of the position is not a drift of the radius)" % (c2["vs"],),
                  ok, f"bars: |r_th - exact| <= {c['bar']:g}, |rdot| <= {c['rdot']:g}" + (f", orders in {c['order']}" if kind == "scale" else ""))


def chk_PL3():
    c = CRIT["pl3"]
    rng = np.random.default_rng(255)
    tau = np.arange(0.0, 14.0, 0.0005)
    worst_a, worst_k = 0.0, 0.0
    for inf in c["inf"]:
        r = 1.0 + inf * tau + (c["c"] / c["kappa"]) * (1.0 - np.exp(-c["kappa"] * tau)) + c["noise"] * rng.standard_normal(len(tau))
        f = fit_rdot(tau, tau * (SQ2 if lmut("tau_unit") else 1.0), r, np.zeros_like(tau))
        worst_a = max(worst_a, abs(f["inf"] - inf))
        worst_k = max(worst_k, abs(f["kappa"] / c["kappa"] - 1) if f["kappa"] == f["kappa"] else 1.0)
        print(f"    PL3 planted rdot_inf = {inf:+.4f}, kappa = {c['kappa']}: recovered {f['inf']:+.5f}, kappa {f['kappa']:.4f}")
    check("PL3", "planted rdot(tau) = rdot_inf + c exp(-kappa tau) through the same rdot and fit code: rdot_inf and kappa recovered", worst_a <= c["bar"] and worst_k <= c["kbar"],
          f"worst |rdot_inf error| {worst_a:.1e} (bar {c['bar']:.0e}), worst relative kappa error {worst_k:.1e} (bar {c['kbar']:.0e})")


def chk_PL5():
    """#213's five synthetic throat series with known classes through the classification used on the locus (as throat_locus.chk_PL)."""
    ok, rows = True, []
    for name, cls_exp, sl_exp, (tau, r, H, rho, w) in xp.planted_series():
        out = xp.classify(tau, r, H, rho, w)
        met, vals = xp.slope_test(r, H, W_EXACT[2])
        ok &= out["class"] == cls_exp and met == sl_exp
        rows.append(f"{name}: class {out['class']} (expected {cls_exp}), slope test {'met' if met else 'not met'} (expected {'met' if sl_exp else 'not met'})")
    for x in rows:
        print("    PL5 " + x)
    check("PL5", "planted series with known classes (#213's five) through the classification used for the expanding member", ok)


def chk_PX(Tref, Pref):
    """The gates PS, PR, PV on planted B families of known exponent: pass at alpha = 2, fail at 1 and 3."""
    for alpha in CRIT["px"]["alphas"]:
        T, P = planted_B(Tref, Pref, alpha)
        ev = eval_predictions(T, P, CRIT["px"]["hk"][-CRIT["fin"]:], tag=f"PX alpha = {alpha:g}")
        want = alpha == 2.0
        ok = (ev["okS"] and ev["okR"] and ev["okV"]) if want else not (ev["okS"] or ev["okR"] or ev["okV"])
        check(f"PX{alpha:g}", f"planted B family r_th - 1, p*_h, rdot ~ B^{alpha:g}: the gates PS, PR, PV " + ("pass" if want else "all fail"), ok, f"PS {ev['okS']}, PR {ev['okR']}, PV {ev['okV']}")


# --------------------------------------------------------------------------------------------------------------------
# controls: #40's events (CT1), #225's stored pair (CT2), #225's X = 100 run (CT3), the offset of the pair (DL)
# --------------------------------------------------------------------------------------------------------------------
def job_ctl(out, D=5, hk=0):
    """CT1: CRIT['ctl']['pts'] of #40's recorded points of G (B = 0.1) rerun through make_data_B / the sweep. CT2: the pair p*_h -/+ 1e-9 at B = 0.1 rerun (#225's stall pair)."""
    os.makedirs(out, exist_ok=True)
    n = D - 2
    h = threshold.CRIT["H"][n][hk]
    old = tl.points40(D, "G", hk)
    keys = sorted(old)
    sel = [keys[int(i)] for i in np.linspace(0, len(keys) - 1, CRIT["ctl"]["pts"])]
    smax = tl.tcap_smax(h, tl.TCAP_BR)
    rows = []
    for A in sel:
        r = fate_run(n, 0.1, h, A, smax)
        o = old[A]
        rows.append(dict(A=A, fate=r["fate"], T=r["T"], cell=r["cell"], old_fate=o["fate"], old_T=o["T"], old_cell=[o["j"], o["i"]]))
        print(f"CT1 A = {A:.15e}: fate {r['fate']} (#40 {o['fate']}), T {r['T']} (#40 {o['T']})", flush=True)
    json.dump(rows, open(os.path.join(out, f"ct1_D{D}_hk{hk}.json"), "w"), indent=1)
    ps = tl.pstar_raw40(D, "G", hk)[0] + (3e-9 if lmut("pstar_off") else 0.0)       # pstar_off: p*_h shifted by 3e-9 (the pair at 1e-9 no longer straddles)
    job_stall(D, 0.1, hk, out, CRIT["delta_ctl"], pstar=ps)
    job_stall(D, 0.1, hk, out, CRIT["delta"], pstar=ps)
    job_bisect(4, CRIT["B_ref"], 0, os.path.join(out, "bisect"))
    return True


def chk_CT4(ctl_dir, D=4, hk=0):
    """CT4: the bisection code at B = B_ref reproduces #40's p*_h (D = 4, h = 0.02) to the tolerance."""
    h = threshold.CRIT["H"][D - 2][hk]
    try:
        s = json.load(open(os.path.join(ctl_dir, "bisect", f"D{D}_B{CRIT['B_ref']:g}_h{h:g}", "summary.json")))
    except OSError:
        check("CT4", "control: the bisection at B = 0.1 reproduces #40's p*_h", False, "the control run is missing")
        return
    p40 = tl.pstar_raw40(D, "G", hk)[0]
    d = abs(s["pstar"] - p40)
    check("CT4", f"control: the bisection of this work at B = 0.1 (D = {D}, h = {h:g}, tolerance {CRIT['bis_tol']:g}) gives #40's p*_h within 2 x tolerance", d <= 2 * CRIT["bis_tol"] and s["stopped"] is None,
          f"p* {s['pstar']:.15e}, #40's {p40:.15e}, difference {d:.1e}, {s['runs']} runs")


def chk_CT12(ctl_dir, D=5, hk=0):
    h = threshold.CRIT["H"][D - 2][hk]
    try:
        rows = json.load(open(os.path.join(ctl_dir, f"ct1_D{D}_hk{hk}.json")))
    except OSError:
        check("CT1", "control: #40's points reproduced", False, "the control run is missing")
        return
    bad = [r for r in rows if r["T"] != r["old_T"] or r["fate"] != r["old_fate"] or r["cell"] != r["old_cell"]]
    check("CT1", f"control: {len(rows)} of #40's recorded points of G (B = 0.1, D = {D}, h = {h:g}) rerun through make_data_B give the same fate, event cell and time bit for bit", not bad,
          f"{len(rows) - len(bad)} of {len(rows)} equal")
    a = [os.path.join(ctl_dir, stall_name(D, 0.1, h, CRIT["delta_ctl"], sd)) for sd in (-1, 1)]
    b = [os.path.join(RAW225, "stall", f"stall_D{D}_h{h:g}_side{sd}.npz") for sd in (-1, 1)]
    try:
        worst = 0.0
        for fa, fb in zip(a, b):
            La, Lb = tl.load_run(fa), tl.load_run(fb)
            if len(La["s"]) != len(Lb["s"]):
                worst = float("inf")
                continue
            for k in ("r", "l", "f", "tau", "v"):
                worst = max(worst, float(np.max(np.abs(La[k] - Lb[k]))))
        fates = [tl.load_run(f)["meta"]["rec"]["fate40"] for f in a]
        check("CT2", f"control: the pair p*_h -/+ 1e-9 at B = 0.1 (D = {D}, h = {h:g}) rerun through job_stall gives #225's stored locus series (r, l, f, tau, v)", worst <= CRIT["ctl"]["arr"] and fates == [-1, 1],
              f"max |difference| {worst:.1e}, bar {CRIT['ctl']['arr']:g}; fates {fates}")
    except OSError as exc:
        check("CT2", "control: the pair rerun against #225's", False, repr(exc))


# --------------------------------------------------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------------------------------------------------
def load_pair(path_fn):
    a, b = path_fn(-1), path_fn(1)
    if not (os.path.exists(a) and os.path.exists(b)):
        return None
    return tl.load_run(a), tl.load_run(b)


def pair_entry(D, B, hk, La, Lb, ps):
    return dict(D=D, B=B, hk=hk, h=La["h"], S=pair_series(La, Lb), La=La, Lb=Lb, pstar=ps, delta=La["meta"]["delta"] if "delta" in La["meta"] else CRIT["delta_ctl"],
                fate=(La["meta"]["rec"]["fate40"], Lb["meta"]["rec"]["fate40"]), fate_loc=(La["meta"]["rec"]["fate_loc"], Lb["meta"]["rec"]["fate_loc"]))


def load_record(raw):
    """raw = .../I255-dh610/runs. Pairs at delta = 1e-10 (D = 5: the three B and h; D = 4: CRIT['D4']), the bisections, #225's stored pairs (B = 0.1, delta = 1e-9), the expansions."""
    rec = dict(mode="record", T5={}, bis={}, T4={}, ref225={}, ctl=os.path.join(raw, "ctl"))
    for D, Bs, hks, store in ((5, CRIT["Bs"], CRIT["hks"], "T5"), (4, CRIT["D4"]["Bs"], CRIT["D4"]["hks"], "T4")):
        for B in Bs:
            for hk in hks:
                h = threshold.CRIT["H"][D - 2][hk]
                pr = load_pair(lambda sd: os.path.join(raw, "stall", stall_name(D, B, h, CRIT["delta"], sd)))
                if pr is None:
                    continue
                if B != CRIT["B_ref"]:
                    sp = os.path.join(raw, "bisect", f"D{D}_B{B:g}_h{h:g}", "summary.json")
                    if os.path.exists(sp):
                        rec["bis"][(D, B, hk)] = json.load(open(sp))
                ps = pstar_of(D, B, hk, os.path.join(raw, "bisect")) if (B == CRIT["B_ref"] or (D, B, hk) in rec["bis"]) else float("nan")
                rec[store][(B, hk)] = pair_entry(D, B, hk, pr[0], pr[1], ps)
    for D in (4, 5):
        for hk in (0, 1, 2):
            h = threshold.CRIT["H"][D - 2][hk]
            pr = load_pair(lambda sd: os.path.join(RAW225, "stall", f"stall_D{D}_h{h:g}_side{sd}.npz"))
            if pr:
                rec["ref225"][(D, hk)] = pair_entry(D, 0.1, hk, pr[0], pr[1], tl.pstar_raw40(D, "G", hk)[0])
                rec["ref225"][(D, hk)]["delta"] = 1e-9
    rec["expand"] = []
    for hk in CRIT["exp_hk"]:
        h = threshold.CRIT["H"][3][hk]
        f = os.path.join(raw, "expand", f"G_D5_{CRIT['member']}_h{h:g}.npz")
        if os.path.exists(f):
            L = tl.load_run(f)
            rec["expand"].append(dict(D=5, fam="G", member=CRIT["member"], h=L["h"], L=L, meta=L["meta"]))
    rec["expand225"] = {}
    for hk in CRIT["exp_hk"]:
        h = threshold.CRIT["H"][3][hk]
        f = os.path.join(RAW225, "expand", f"G_D5_{CRIT['member']}_h{h:g}.npz")
        if os.path.exists(f):
            rec["expand225"][hk] = tl.load_run(f)
    return rec


def dev_suite():
    """Development suite: #225's stored pairs stand for the B family (planted exponent 2 on top of them, PX), the controls and planted evolutions are run, #225's X = 100 runs are the expansions."""
    rec = load_record(os.path.join(os.path.expanduser("~/Source/negative-mass-raw/I255-dh610"), "runs"))
    rec["mode"] = "dev"
    rec["expand"] = [dict(D=5, fam="G", member=CRIT["member"], h=L["h"], L=L, meta=L["meta"]) for L in rec["expand225"].values()]
    rec["ctl"] = os.path.join(os.path.expanduser("~/Source/negative-mass-raw/I255-dh610"), "dev", "ctl")
    h = threshold.CRIT["H"][3][0]
    pr = load_pair(lambda sd: os.path.join(rec["ctl"], stall_name(5, 0.1, h, CRIT["delta"], sd)))
    rec["T5"] = {(0.1, 0): pair_entry(5, 0.1, 0, pr[0], pr[1], tl.pstar_raw40(5, "G", 0)[0])} if pr else {}
    return rec


# --------------------------------------------------------------------------------------------------------------------
# the expanding member
# --------------------------------------------------------------------------------------------------------------------
def exp_summaries(E):
    for e in E:
        e["S"] = tl.expansion_summary(e["L"], W_EXACT[e["D"] - 2])
    return E


def chk_exp(E, E225, full):
    exp_summaries(E)
    ok_kx, ok_res = True, True
    lines = []
    tsmax = CRIT["ft_tsmax"] * (1.1 if lmut("ts_bound") else 1.0)           # ts_bound: the bound on tau_s / tau_end relaxed to 1.32
    for e in sorted(E, key=lambda e: -e["h"]):
        cl, L = e["S"]["cl"], e["L"]
        ratio = cl["tau_s"] / cl["tau_end"]
        isFT = cl["class"] == "FT" or (lmut("ts_bound") and cl["class"] == "GROW" and ratio <= tsmax)
        ok_kx &= isFT and ratio <= tsmax
        print(f"  --- G D={e['D']} A = {e['member']} h = {e['h']:g} (X = {e['meta']['X']:g}, {e['meta'].get('secs', float('nan')):.0f} s) ---")
        for x in tl.expand_lines(e, W_EXACT[e["D"] - 2], e["D"] - 2):
            print(x)
        print(f"    final window (last e-fold of r, {cl['npts']} points, r/a = {cl['r0']:.1f} to {cl['r1']:.1f}): H max/min {cl['Hratio']:.3f}, H end/start {cl['Hgrow']:.3f}, |rho| max/min {cl['rhoratio']:.3f}, "
              f"w_bar in [{cl['wmin']:.4f}, {cl['wmax']:.4f}]; 1/H linear in tau: k = {cl['k']:.4f}, tau_s = {cl['tau_s']:.4f} (last point {cl['tau_end']:.4f}, ratio {ratio:.4f}), max residual {cl['ft_resid']:.1e} -> class {cl['class']}")
        print(f"    slope test: H/W at r/a = " + ", ".join(f"{k:g}: {v:.3f}" for k, v in e["S"]["slope_vals"].items()) + f"; locus drift l at the end {L['l'][e['S']['ix']][-1]:+.3f}")
    check("KX2", f"G, D = 5, A = {CRIT['member']}, X = {CRIT['X']:g}: class FT with tau_s / tau_end <= {CRIT['ft_tsmax']:g} (#213's bound) at both h", ok_kx and len(E) >= 2,
          "; ".join(f"h = {e['h']:g}: {e['S']['cl']['class']}, k {e['S']['cl']['k']:.4f}, ratio {e['S']['cl']['tau_s'] / e['S']['cl']['tau_end']:.4f}" for e in sorted(E, key=lambda e: -e["h"])))
    E = sorted(E, key=lambda e: -e["h"])
    if len(E) >= 2:
        a, b = E[0], E[1]
        ca, cb = a["S"]["cl"], b["S"]["cl"]
        dH = []
        for r0 in (10.0, 25.0, 50.0):
            ra, rb = a["L"]["r"][a["S"]["ix"]], b["L"]["r"][b["S"]["ix"]]
            if (ra >= r0).any() and (rb >= r0).any():
                dH.append(abs(np.interp(r0, ra, a["S"]["obs"]["H"][a["S"]["ix"]]) / np.interp(r0, rb, b["S"]["obs"]["H"][b["S"]["ix"]]) - 1))
        dk, dts = abs(ca["k"] / cb["k"] - 1), abs(ca["tau_s"] * (1.02 if lmut("ts_scale") else 1.0) / cb["tau_s"] - 1)       # ts_scale: the coarse run's tau_s scaled by 1.02
        c = CRIT["res"]
        ok_res = ca["class"] == cb["class"] and (max(dH) if dH else 0) <= c["H_rel"] and dk <= c["k_rel"] and dts <= c["ts_rel"]
        print(f"    RES2 h {a['h']:g} vs {b['h']:g}: class {ca['class']} / {cb['class']}; tau_s {ca['tau_s']:.4f} / {cb['tau_s']:.4f}; k {ca['k']:.4f} / {cb['k']:.4f}; H at fixed r/a {fmt(dH)}")
        check("RES2", f"the two h agree: same class, H at fixed r/a within {c['H_rel']:.0%}, k and tau_s within {c['k_rel']:.0%}", ok_res)
    if full and E225:
        worst = 0.0
        for hk, L225 in E225.items():
            h = L225["h"]
            Ln = [e["L"] for e in E if abs(e["h"] - h) < 1e-12]
            if not Ln:
                continue
            Ln = Ln[0]
            N100 = L225["N"]
            m = Ln["s"] <= N100
            m2 = L225["s"] <= N100
            if m.sum() != m2.sum():
                worst = float("inf")
                continue
            for k in ("r", "l", "f", "tau"):
                worst = max(worst, float(np.max(np.abs(Ln[k][m] - L225[k][m2]))))
        check("CT3", "control: the X = 160 runs reproduce #225's X = 100 runs (locus r, l, f, tau on every slice that lies inside both rectangles)", worst <= CRIT["ctl"]["arr"], f"max |difference| {worst:.1e}, bar {CRIT['ctl']['arr']:g}")
    return ok_kx


# --------------------------------------------------------------------------------------------------------------------
# the analysis
# --------------------------------------------------------------------------------------------------------------------
def report_pair(e, n):
    S = e["S"]
    d = locus_diag(S, n)
    print(f"    pair D={e['D']} B={e['B']:g} h={e['h']:g} delta={e['delta']:g}: p*_h {e['pstar']:.12e}; slices {len(S['s'])}; fates {e['fate']} on the locus {e['fate_loc']}")
    for tt in (6.0, 8.0, 10.0, 12.0, 13.0):
        k = int(np.argmin(np.abs(S["t"] - tt)))
        if abs(S["t"][k] - tt) > 0.1:
            continue
        print(f"       t {S['t'][k]:5.2f} tau {S['tau'][k]:7.3f} r_th {S['r'][k]:.5f} odd {S['odd'][k]:+.1e} l {S['l'][k]:+.3f} v {S['v'][k]:+.4f} f {S['f'][k]:+.4f} N_th {d['Nth'][k]:.4f} "
              f"rho_g/rho_eq - 1 {d['ratio'][k] - 1:+.1e} q_t {d['qt'][k]:+.4f}")


def chk_RV(S, label):
    """Reproduction: rdot_inf of the B = 0.1 tuned series at the finest h is #225-based dev_04's value 0.02602 within 5 %."""
    f = fit_rdot(S["t"], S["tau"], S["r"], S["odd"])
    ok = abs(f["inf"] / 0.02602 - 1) <= 0.05
    check("RV", f"reproduction: rdot_inf of the B = 0.1 tuned series ({label}) at the finest h is dev_04's 0.02602 within 5 %", ok, f"rdot_inf {f['inf']:+.5f} +- {f['se_inf']:.1e}, kappa {f['kappa']:.4f}, n {f['n']}")


def chk_DL(rec):
    """The tuned mean at B = B_ref from #225's pair (delta = 1e-9) and from this work's (delta = 1e-10) on the same slices."""
    for hk in CRIT["hks"]:
        if (5, hk) in rec["ref225"] and (CRIT["B_ref"], hk) in rec["T5"]:
            a, b = rec["ref225"][(5, hk)]["S"], rec["T5"][(CRIT["B_ref"], hk)]["S"]
            t0, t1 = CRIT["dl"]["t"]
            m = (a["t"] >= t0) & (a["t"] <= t1) & np.isin(a["s"], b["s"])
            ia = np.where(m)[0]
            ib = np.searchsorted(b["s"], a["s"][ia])
            dd = float(np.max(np.abs(a["r"][ia] - b["r"][ib]))) if len(ia) else float("inf")
            check(f"DL[h index {hk}]", f"the tuned mean at delta = 1e-9 (#225) and 1e-10 differs by at most {CRIT['dl']['bar']:g} for t in {CRIT['dl']['t']}", dd <= CRIT["dl"]["bar"], f"max |difference| {dd:.1e}")


def analyze_suite(rec):
    print(f"throat_drift.py analyze: mode {rec['mode']}, mutation {LMUT or dn.MUT}, numba {dn.nb.__version__ if dn.HAVE_NUMBA else 'absent'}, numpy {np.__version__}, {time.strftime('%Y-%m-%d %H:%M:%S')}")
    full = rec["mode"] == "record"
    n5, n4 = 3, 2
    W5 = W_EXACT[n5]
    I2, I1 = pulse_integrals(n5)
    print(f"  pulse integrals on the initial ray (D = 5): I_2 = {I2:.5f}, I_1 = {I1:.5f}")
    # planted and controls
    chk_PL12()
    chk_PL3()
    chk_PL5()
    ref = {(hk): rec["ref225"][(5, hk)] for hk in CRIT["hks"] if (5, hk) in rec["ref225"]}
    if len(ref) == 3:
        Tref = {(0.1, hk): ref[hk]["S"] for hk in ref}
        Pref = {(0.1, hk): ref[hk]["pstar"] for hk in ref}
        chk_PX(Tref, Pref)
    chk_CT12(rec["ctl"])
    chk_CT4(rec["ctl"])
    # D = 5 pairs
    if full:
        T5 = rec["T5"]
        for (B, hk), e in sorted(T5.items()):
            report_pair(e, n5)
        miss = [(B, hk) for B in CRIT["Bs"] for hk in CRIT["hks"] if (B, hk) not in T5]
        for (B, hk), e in sorted(T5.items()):
            ok = e["fate"] == (-1, 1) and e["fate_loc"] == (-1, 1)
            bis = rec["bis"].get((5, B, hk))
            okb = True if B == CRIT["B_ref"] else (bis is not None and bis["stopped"] is None and bis["error"] is None and bis["hi"] - bis["lo"] <= CRIT["bis_tol"] * (1 + 1e-9))
            check(f"JOB[D=5 B={B:g} h={e['h']:g}]", "the pair p*_h -/+ delta has fates (collapse, expansion) on #40's events and on the locus's levels" + ("" if B == CRIT["B_ref"] else "; the bisection finished to the declared tolerance"), ok and okb,
                  f"fates {e['fate']}, locus {e['fate_loc']}" + ("" if bis is None else f"; bracket {bis['hi'] - bis['lo']:.1e}, {bis['runs']} runs"))
            S = e["S"]
            w = (S["t"] >= CRIT["win"][0]) & (S["t"] <= CRIT["win"][1])
            odd = float(np.max(np.abs(S["odd"][w]))) if w.any() else float("inf")
            qs = locus_diag(S, n5)["ratio"]
            wq = (S["t"] >= CRIT["qs"]["win"][0]) & (S["t"] <= CRIT["qs"]["win"][1])
            dq = float(np.max(np.abs(qs[wq] - 1))) if wq.any() else float("inf")
            check(f"NO[D=5 B={B:g} h={e['h']:g}]", f"the pair's odd part (the unstable mode) is at most {CRIT['odd_max']:g} for t in {CRIT['win']}: the pair's mean is the tuned member", odd <= CRIT["odd_max"], f"max |odd| {odd:.1e}")
            check(f"QS[D=5 B={B:g} h={e['h']:g}]", f"on the locus the ghost's energy density is a static throat's of the same radius, |rho_g / rho_eq - 1| <= {CRIT['qs']['bar']:g} for t in {CRIT['qs']['win']}", dq <= CRIT["qs"]["bar"], f"max {dq:.1e}")
        if not miss:
            T = {(B, hk): T5[(B, hk)]["S"] for (B, hk) in T5}
            P = {(B, hk): T5[(B, hk)]["pstar"] for (B, hk) in T5}
            hk_fin = CRIT["hks"][-CRIT["fin"]:]
            ev = eval_predictions(T, P, CRIT["hks"], tag="D=5")
            # gates on the finest h only
            okS = all(CRIT["pw"]["lo"] <= r[1] <= CRIT["pw"]["hi"] for r in ev["rows"] if r[0] in hk_fin)
            okR = all(all(CRIT["pw"]["lo"] <= a <= CRIT["pw"]["hi"] for a in r[3]) for r in ev["rows"] if r[0] in hk_fin)
            okV = all((CRIT["pw"]["lo"] <= r[5] <= CRIT["pw"]["hi"]) and (CRIT["kap"]["lo"] <= r[7] <= CRIT["kap"]["hi"]) for r in ev["rows"] if r[0] in hk_fin)
            check("PS", f"prediction: p*_h(B) ~ B^alpha with alpha in {tuple(CRIT['pw'].values())} between B = 0.04 and 0.1 at the finest {CRIT['fin']} h (first order would give 1)", okS, "exponents " + fmt([r[1] for r in ev["rows"] if r[0] in hk_fin], "{:.3f}"))
            check("PR", f"prediction: r_th(t) - 1 ~ B^alpha, alpha in {tuple(CRIT['pw'].values())}, between B = 0.04 and 0.1 at t = {CRIT['ts']} and the finest {CRIT['fin']} h", okR, "exponents " + "; ".join(fmt(r[3], "{:.3f}") for r in ev["rows"] if r[0] in hk_fin))
            check("PV", f"prediction: rdot_inf(B) ~ B^alpha, alpha in {tuple(CRIT['pw'].values())}, and kappa(0.04) / kappa(0.1) in {tuple(CRIT['kap'].values())}", okV,
                  "exponents " + fmt([r[5] for r in ev["rows"] if r[0] in hk_fin], "{:.3f}") + ", kappa ratios " + fmt([r[7] for r in ev["rows"] if r[0] in hk_fin], "{:.3f}"))
            for B in CRIT["Bs"]:
                ok, q, two = eval_convergence(T, B, CRIT["hks"], tag="D=5")
                check(f"CV[B={B:g}]", f"r_th(10) - 1 converges at second order over the three h (order in {CRIT['cv']['order']}) and the two finest agree within {CRIT['cv']['two_h']:.0%}", ok, f"order {q:.2f}, two finest {two:.1e}")
            print("  energies the pulses leave, E_net = (B^2 - A*^2) I_2, against rdot_inf (finest h):")
            for B in CRIT["Bs"]:
                hk = CRIT["hks"][-1]
                En = (B ** 2 - P[(B, hk)] ** 2) * I2
                f = ev["fits"][(B, hk)]
                print(f"    B = {B:g}: A* = {P[(B, hk)]:.6e}, E_net = {En:.5f}, rdot_inf = {f['inf']:+.5f}, rdot_inf / E_net = {f['inf'] / En:.5f}, kappa {f['kappa']:.4f}")
        # masses on the finest h
        hk = CRIT["hks"][-1]
        print("  Misner-Sharp mass at fixed areal radius on the + and - side of the collapse-side run, finest h:")
        for B in CRIT["Bs"]:
            e = T5.get((B, hk))
            if not e:
                continue
            La = e["La"]
            for tt in (8.0, 10.0, 12.0):
                sk = min(La["prof"], key=lambda s: abs(La["t"][np.searchsorted(La["s"], s)] - tt)) if La["prof"] else None
                if sk is None:
                    continue
                Mp, Mm, rn_, ic, M = tl.slice_masses(La, sk, n5, np.array([2.0, 3.0, 4.0]))
                print(f"    B = {B:g} t = {tt:g}: M_+(R = 2, 3, 4) {fmt(Mp, '{:+.4f}')}, M_-(R = 2, 3, 4) {fmt(Mm, '{:+.4f}')}")
        # D = 4, reported
        if rec["T4"]:
            print("  D = 4, reported (not gated):")
            for (B, hk), e in sorted(rec["T4"].items()):
                S = e["S"]
                f = fit_rdot(S["t"], S["tau"], S["r"], S["odd"])
                print(f"    D=4 B={B:g} h={e['h']:g}: p*_h {e['pstar']:.6e}; r_th - 1 at t = {CRIT['ts']}: {fmt(Rt(S), '{:+.5f}')}; rdot_inf {f['inf']:+.5f} +- {f['se_inf']:.1e}, kappa {f['kappa']:.3f}, n {f['n']}")
    else:
        # development suite: the stored pairs at B = 0.1 (delta = 1e-9) as the family with planted exponent 2
        Tref = {(0.1, hk): ref[hk]["S"] for hk in ref}
        Pref = {(0.1, hk): ref[hk]["pstar"] for hk in ref}
        T, P = planted_B(Tref, Pref, 2.0)
        for hk in CRIT["hks"]:
            for B in CRIT["Bs"]:
                if (B, hk) in T:
                    pass
        ev = eval_predictions(T, P, CRIT["hks"], tag="dev planted alpha = 2")
        hk_fin = CRIT["hks"][-CRIT["fin"]:]
        check("PS", "prediction gate on the planted family (exponent 2)", ev["okS"])
        check("PR", "prediction gate on the planted family (exponent 2)", ev["okR"])
        check("PV", "prediction gate on the planted family (exponent 2)", ev["okV"])
        for B in CRIT["Bs"]:
            ok, q, two = eval_convergence(T, B, CRIT["hks"], tag="dev")
            check(f"CV[B={B:g}]", "convergence on the planted family", ok, f"order {q:.2f}, two finest {two:.1e}")
        for hk in CRIT["hks"]:
            e = ref[hk]
            S = e["S"]
            w = (S["t"] >= CRIT["win"][0]) & (S["t"] <= CRIT["qs"]["win"][1])
            qs = locus_diag(S, n5)["ratio"]
            wq = (S["t"] >= CRIT["qs"]["win"][0]) & (S["t"] <= CRIT["qs"]["win"][1])
            dq = float(np.max(np.abs(qs[wq] - 1)))
            check(f"QS[h index {hk}]", "quasi-static diagnosis on the stored pair", dq <= CRIT["qs"]["bar"], f"max {dq:.1e}")
    chk_DL(rec)
    if (CRIT["B_ref"], CRIT["hks"][-1]) in rec["T5"]:
        chk_RV(rec["T5"][(CRIT["B_ref"], CRIT["hks"][-1])]["S"], "delta = 1e-10")
    elif ref:
        chk_RV(ref[CRIT["hks"][-1]]["S"], "#225's pair, delta = 1e-9")
    if rec["expand"]:
        chk_exp(rec["expand"], rec.get("expand225") if full else None, full)
    return {c: ok for c, ok in results}


def main():
    global LMUT
    argv = sys.argv[1:]
    if not argv:
        sys.exit(__doc__)
    cmd, rest = argv[0], argv[1:]
    a, mut = {}, None
    for x in rest:
        if x.startswith("--mutate="):
            mut = x.split("=", 1)[1]
        elif x.startswith("--") and "=" in x:
            k, v = x[2:].split("=", 1)
            a[k] = v
    if mut:
        if mut in ("ellis_charge", "theta_coef", "scalar_coef", "ghost_sign", "no_corrector", "data_constraint_sign", "finder_sign", "constraint_sign"):
            dn.MUT = mut
        else:
            LMUT = mut
    import platform
    la = open("/proc/loadavg").read().split()[:3] if platform.system() == "Linux" else ["?"]
    print(f"throat_drift.py {cmd} {a}: numba {dn.nb.__version__ if dn.HAVE_NUMBA else 'absent'}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"load average at start {' '.join(la)}, mutation {mut}, {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    t0 = time.time()
    D = int(a.get("D", 5))
    if cmd == "bisect":
        job_bisect(D, float(a["B"]), int(a["hk"]), a["out"])
    elif cmd == "stall":
        job_stall(D, float(a["B"]), int(a["hk"]), a["out"], float(a.get("delta", CRIT["delta"])), pstar=float(a["pstar"]) if "pstar" in a else None, bis_dir=a.get("bis"))
    elif cmd == "bisect_stall":
        B, hk = float(a["B"]), int(a["hk"])
        ok = job_bisect(D, B, hk, os.path.join(a["out"], "bisect"))
        if ok:
            job_stall(D, B, hk, os.path.join(a["out"], "stall"), float(a.get("delta", CRIT["delta"])), bis_dir=os.path.join(a["out"], "bisect"))
        else:
            sys.exit("bisection did not finish")
    elif cmd == "expand":
        tl.job_expand(5, int(a["hk"]), CRIT["member"], float(a.get("X", CRIT["X"])), a["out"], fam="G")
    elif cmd == "ctl":
        job_ctl(a["out"], D=D, hk=int(a.get("hk", 0)))
    elif cmd in ("analyze", "dev", "mutate"):
        def one(m):
            global LMUT
            results.clear()
            LMUT = m
            suite = load_record(a["raw"]) if cmd == "analyze" else dev_suite()
            if "ctl" in a and m:                       # the mutated run reads the control directory made under the mutation
                suite["ctl"] = a["ctl"]
            return analyze_suite(suite)
        if cmd == "mutate":
            base = one(None)
            print("\n==== baseline done; now the mutation ====\n", flush=True)
            mm = one(mut)
            changed = sorted(c.split("[")[0] for c in set(base) | set(mm) if base.get(c) != mm.get(c))
            print(f"\nmutation {mut}: verdicts changed against the unmutated suite: {changed}")
            print(f"{time.time() - t0:.1f} s")
            return
        ok = one(mut)
        bad = [c for c, v in ok.items() if v is False]
        print(f"\n{sum(1 for v in ok.values() if v)} of {len(ok)} checks pass; {len(bad)} fail" + (f": {', '.join(bad)}" if bad else ""))
        print(f"{time.time() - t0:.1f} s")
        sys.exit(1 if bad else 0)
    else:
        sys.exit(f"unknown command {cmd}")
    print(f"{time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
