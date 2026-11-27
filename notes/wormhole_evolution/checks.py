#!/usr/bin/env python3
"""Checks of the double-null evolution dn.py (#39, step 2). Every check can fail; each has a --mutate= whose targets are in
MUTATIONS. The criteria are the CRIT table, as declared in negative-mass-raw I39-dh610/README.md ("Step 2 criteria").

    python3 checks.py                         all checks, numba kernel; exit 1 if any check fails
    python3 checks.py --kernel=numpy          the reference kernel (no numba needed); K is skipped
    python3 checks.py --only=P1,P3            selected checks
    python3 checks.py --mutate=NAME           break one input on purpose (the targets in MUTATIONS must fire)
    python3 checks.py --quick                 coarser resolutions (smoke test; not a run of record)
    python3 checks.py --json=out.json         write the numbers

Checks: I1 K P1 P2 P3 C1 C2 G1 S F1 U. Units a = 1, 8 pi G = 1, t = (x+ + x-)/sqrt2, tau the proper time on the diagonal.
"""
import json
import math
import sys
import time

import numpy as np

import dn

SQ2 = math.sqrt(2.0)
W_EXACT = {2: math.sqrt(1.39705243433350684), 3: math.sqrt(2.98495893035176025), 4: math.sqrt(4.68662054311946479)}  # #195

# criteria: declared in the raw README before the first run of record (see "Step 2 criteria" there)
CRIT = dict(
    I1_tol=1e-8,                    # max |r solved - r Ellis| over n = 2, 3, 4 at h = 0.02, N = 400
    K_tol=1e-12,                    # numba against numpy, max abs difference of every field on the declared grids
    order_lo=1.75, order_hi=2.25,   # observed order at each halving, P1 P2 P3 C1 C2
    P1_cases=[(2, 4.0, (0.04, 0.02, 0.01)), (3, 2.7, (0.04 / SQ2, 0.02 / SQ2, 0.01 / SQ2))],   # (n, t_c, h list)
    P2_h=(0.005, 0.0025, 0.00125), P2_X=0.6,
    P3_h=(0.04, 0.02, 0.01), P3_slope_tol=1e-3,
    C_h=(0.04, 0.02, 0.01), C_X=3.0, C_pulse=dict(A=0.05, xc=1.5, w=0.5),
    G1_h={2: (0.02, 0.01), 3: (0.01 / SQ2, 0.005 / SQ2)}, G1_X={2: 20.0, 3: 15.0}, G1_bracket=0.01, G1_tolA=1e-10,
    G1_window=(1e-3, 1e-2),         # |H|/W after the global minimum of |H|; slope averaged over the two sides of p*
    G1_bar=0.01, G1_agree=0.005,    # |W_fit/W - 1| <= bar at each h; |W_h/W_h/2 - 1| <= agree
    S_dt=4.0,                       # horizon-formation time agrees within S_dt * h
    S_rcut=0.025, S_fcut=50.0,      # r_cut halved, F_cut doubled (defaults 0.05, 25)
    F1_h=0.02, F1_X=12.0, F1_pulse=dict(A=0.04, xc=3.0, w=1 / math.sqrt(3.0), kind='dphi'),   # [SH02] (3.1): c_a = 0.1 in r d phi, r = 2.35 at the pulse; c_b = 3
    trap_tol=1e-2,
)

MUTATIONS = {   # name: checks that must fire; computed in development runs dev_s2_16 and dev_s2_17 (negative-mass-raw I39-dh610/dev)
    "theta_coef": ["P1", "P2", "P3", "C1", "C2", "G1", "F1"],
    "scalar_coef": ["P1", "P3", "C1", "C2", "G1", "F1"],
    "ghost_sign": ["P1", "P3", "C1", "C2", "G1", "F1"],
    "no_corrector": ["P2", "P3", "C1", "C2"],
    "numpy_Spert": ["K"],
    "data_constraint_sign": ["I1", "P1", "C1", "C2", "G1", "F1"],
    "ellis_charge": ["I1", "P1", "C1", "C2", "G1", "F1"],
    "constraint_sign": ["C1"],
    "flux_sign": ["C2"],
    "mass_exponent": ["P2", "C2"],
    "finder_sign": ["P3", "G1", "F1", "S"],
    "finder_swap": ["P2", "P3"],
    "time_unit": ["G1"],
    "cut_extreme": ["S"],
}

results = []
INFO = {}


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)          # numpy booleans would not compare `is False`
    results.append((cid, ok))
    tag = "SKIP" if ok is None else ("pass" if ok else "FAIL")
    print(f"[{tag}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def orders(errs):
    return [math.log2(errs[i] / errs[i + 1]) if errs[i + 1] > 0 and errs[i] > 0 else float("inf") for i in range(len(errs) - 1)]


def order_ok(errs):
    o = orders(errs)
    return all(CRIT["order_lo"] <= x <= CRIT["order_hi"] for x in o), o


def fmt(xs, f="{:.3e}"):
    return "[" + ", ".join(f.format(x) for x in xs) + "]"


def ellis(n, h, X, pulses=(), symmetric=False, eps=(-1.0,), **kw):
    N = int(round(X / h))
    d = dn.ellis_rays(n, N, h, pulses=list(pulses), symmetric=symmetric, eps=eps)
    return d, dn.evolve(d, kernel=KERNEL, tol=CRIT["trap_tol"], **kw)


# ----------------------------------------------------------------------------------------------------------------------
def check_I1():
    errs = []
    for n in (2, 3, 4):
        d = dn.ellis_rays(n, 400, 0.02)
        sol, q = dn.ellis_dense(n, 400 * 0.02 / SQ2 + 1)
        r_ex = sol(np.arange(401) * 0.02 / SQ2)[0]
        errs.append(max(np.max(np.abs(d.Sp[0] - r_ex)), np.max(np.abs(d.Sm[0] - r_ex))))
    check("I1", "the constraint solve on the rays returns the Ellis r(l), n = 2, 3, 4", max(errs) <= CRIT["I1_tol"],
          f"max error {max(errs):.2e}, per n {fmt(errs)}, tol {CRIT['I1_tol']:g}")


def check_K():
    if KERNEL == "numpy" and not dn.HAVE_NUMBA:
        check("K", "numba and numpy kernels agree", None, "numba absent: skipped")
        return
    if not dn.HAVE_NUMBA:
        check("K", "numba and numpy kernels agree", None, "numba absent: skipped")
        return
    worst = 0.0
    same_dead = True
    cases = []
    d1 = dn.ellis_rays(2, 150, 0.02, pulses=[dict(field=0, ray="+", A=0.05, xc=1.5, w=0.5)], symmetric=True)
    cases.append(("ellis D=4 pulse N=150", d1))
    d2 = dn.tangherlini_rays(3, 60, 0.01)
    cases.append(("kruskal D=5 N=60", d2))
    d3 = dn.ellis_rays(2, 150, 0.02, pulses=[dict(field=0, ray="+", A=-0.1, xc=1.5, w=0.5)], symmetric=True)   # collapses: excision
    cases.append(("ellis D=4 collapse N=150", d3))
    for name, d in cases:
        a = dn.evolve(d, kernel="numba", tol=CRIT["trap_tol"])
        b = dn.evolve(d, kernel="numpy", tol=CRIT["trap_tol"])
        live = ~(a.dead.astype(bool) | b.dead.astype(bool))
        diff = np.max(np.abs(a.G[:, live] - b.G[:, live])) if live.any() else 0.0
        same_dead &= bool(np.array_equal(a.dead, b.dead))
        worst = max(worst, diff)
        print(f"    K {name}: max |numba - numpy| {diff:.2e}, dead cells {int(a.dead.sum())} / {int(b.dead.sum())}")
    check("K", "numba and numpy kernels agree to roundoff on three declared grids, with the same dead cells",
          worst <= CRIT["K_tol"] and same_dead, f"max diff {worst:.2e}, tol {CRIT['K_tol']:g}, same dead masks {same_dead}")


def check_P1():
    allok = True
    for n, tc, hs in CRIT["P1_cases"]:
        X = tc / SQ2
        errs = []
        for h in hs:
            N = int(round(X / h))
            d = dn.ellis_rays(n, N, h)
            res = dn.evolve(d, kernel=KERNEL, tol=CRIT["trap_tol"])
            r_ex, ph_ex = dn.ellis_exact_nodes(d)
            errs.append(max(np.max(np.abs(res.G[0] - r_ex)), np.max(np.abs(res.G[1])), np.max(np.abs(res.G[2] - ph_ex))))
        ok, o = order_ok(errs)
        allok &= ok
        print(f"    P1 D={n + 2} t_c={tc}: h {fmt(hs, '{:.5f}')} errors {fmt(errs)} orders {fmt(o, '{:.3f}')}")
    check("P1", "static Ellis wormhole: error against the exact static solution converges at second order for t <= t_c",
          allok, f"orders within [{CRIT['order_lo']}, {CRIT['order_hi']}]")


def check_P2():
    allok = True
    for n in (2, 3, 4):
        e_rf, e_M, e_hx, e_hr = [], [], [], []
        for h in CRIT["P2_h"]:
            N = int(round(CRIT["P2_X"] / h))
            d = dn.tangherlini_rays(n, N, h)
            res = dn.evolve(d, kernel=KERNEL, excise=False, tol=CRIT["trap_tol"])
            r_ex, f_ex = dn.tangherlini_exact_nodes(d)
            e_rf.append(max(np.max(np.abs(res.G[0] - r_ex)), np.max(np.abs(res.G[1] - f_ex))))
            M = dn.mass_nodes(res)
            e_M.append(np.max(np.abs(M[1:-1, 1:-1] - d.M0)))
            xs, ms, rs, sl = dn.contour(res, "+", axis="j")
            e_hx.append(np.max(np.abs(ms - (-d.U0))) if len(ms) else float("inf"))
            e_hr.append(np.max(np.abs(rs - d.rh)) if len(rs) else float("inf"))
        for name, errs in (("(r, f) vs exact", e_rf), ("M - M0", e_M), ("horizon x-", e_hx), ("r on horizon", e_hr)):
            ok, o = order_ok(errs)
            allok &= ok
            print(f"    P2 D={n + 2} {name}: errors {fmt(errs)} orders {fmt(o, '{:.3f}')}")
    check("P2", "Tangherlini in Kruskal coordinates, D = 4, 5, 6: (r, f), M, and the horizon at x- = -U0 with r = (2M)^(1/(n-1)) "
          "converge at second order", allok, f"orders within [{CRIT['order_lo']}, {CRIT['order_hi']}]")


def check_P3():
    a = 0.75
    allok = True
    e_all, e_c1, e_c2 = [], [], []
    mism = []
    slopes = []
    for h in CRIT["P3_h"]:
        N = int(round(2.0 / h))
        d = dn.roberts_rays(N, h, a)
        res = dn.evolve(d, kernel=KERNEL, tol=CRIT["trap_tol"])
        i = np.arange(N + 1)[:, None]
        j = np.arange(N + 1)[None, :]
        r_ex, ph_ex = dn.roberts_exact(1.0 + i * h, 1.0 + j * h, a)
        e_all.append(max(np.max(np.abs(res.G[0] - r_ex)), np.max(np.abs(res.G[1])), np.max(np.abs(res.G[2] - ph_ex))))
        xs, ms, rs, s1 = dn.contour(res, "+")
        xs2, ms2, rs2, s2 = dn.contour(res, "-")
        e_c1.append(np.max(np.abs(xs - ms / (2 * a))) if len(xs) else float("inf"))
        e_c2.append(np.max(np.abs(xs2 - 2 * a * ms2)) if len(xs2) else float("inf"))
        slopes.append((s1, s2))
        # classification of cell centres against the exact r_+, r_-
        c, up, um, alive = dn.cell_fields(res)
        p = 1.0 + (np.arange(N)[:, None] + 0.5) * h
        m = 1.0 + (np.arange(N)[None, :] + 0.5) * h
        rr = np.sqrt(a * (p ** 2 + m ** 2) - p * m)
        rp_ex, rm_ex = (2 * a * p - m) / (2 * rr), (2 * a * m - p) / (2 * rr)

        mism.append(int(np.sum(dn.cell_class(up[0], um[0], 0.0) != np.where((rp_ex > 0) & (rm_ex > 0), 1, np.where((rp_ex < 0) & (rm_ex < 0), -1, 0)))))
    for name, errs in (("(r, f, phi) vs exact", e_all), ("contour r_+ = 0 position", e_c1), ("contour r_- = 0 position", e_c2)):
        ok, o = order_ok(errs)
        allok &= ok
        print(f"    P3 {name}: errors {fmt(errs)} orders {fmt(o, '{:.3f}')}")
    s_ok = all(abs(s[0] - 1 / (2 * a)) <= CRIT["P3_slope_tol"] and abs(s[1] - 2 * a) <= CRIT["P3_slope_tol"] for s in slopes[-1:])
    N_fin = int(round(2.0 / CRIT["P3_h"][-1]))
    m_ok = all(mm <= 3 * int(round(2.0 / h)) for mm, h in zip(mism, CRIT["P3_h"]))
    print(f"    P3 contour slopes dx+/dx- at the finest h: r_+ = 0: {slopes[-1][0]:.6f} (exact {1 / (2 * a):.6f}); r_- = 0: {slopes[-1][1]:.6f} (exact {2 * a:.6f}); "
          f"cells misclassified (trapped / anti-trapped / untrapped) {mism}")
    check("P3", "Roberts ghost solution (D = 4, exact, nonlinear): second-order error, the two trapping-horizon lines and their slopes, and the "
          "cell classification", allok and s_ok and m_ok,
          f"slope tol {CRIT['P3_slope_tol']}, at most 3 misclassified cells per row")


def pulse_run(h, X=None):
    cp = CRIT["C_pulse"]
    X = X or CRIT["C_X"]
    return ellis(2, h, X, pulses=[dict(field=0, ray="+", A=cp["A"], xc=cp["xc"], w=cp["w"])], symmetric=True)


def check_C():
    cmax, fres = [], []
    for h in CRIT["C_h"]:
        d, res = pulse_run(h)
        Cpp, Cmm = dn.constraints(res)
        cmax.append(max(np.nanmax(np.abs(Cpp)), np.nanmax(np.abs(Cmm))))
        fres.append(dn.edge_flux_residual(res))
    e_flux = [f[0] for f in fres]
    scale = [f[1] for f in fres]
    ok1, o1 = order_ok(cmax)
    ok2, o2 = order_ok(e_flux)
    print(f"    C1 max |C_++|, |C_--| {fmt(cmax)} orders {fmt(o1, '{:.3f}')}")
    print(f"    C2 edge flux residual {fmt(e_flux)} orders {fmt(o2, '{:.3f}')}; largest M change on the edge {fmt(scale)}")
    check("C1", "constraints C_++, C_-- of a pulse run fall at second order", ok1, f"max norm at t <= {SQ2 * CRIT['C_X']:.2f}")
    check("C2", "mass change along the edge x+ = X equals the integral of the flux at second order, with a nonzero flux", ok2 and min(scale) > 1e-5,
          "pulse crossing the edge")


# ----------------------------------------------------------------------------------------------------------------------
def fate_run(n, h, A, X, kernel=None):
    N = int(round(X / h))
    d = dn.ellis_rays(n, N, h, pulses=[dict(field=0, ray="+", A=A, xc=3.0, w=0.5)], symmetric=True)
    return dn.evolve(d, kernel=kernel or KERNEL, stop=True, store=False, tol=CRIT["trap_tol"])


def fate(res):
    return 1 if "expand" in res.events else (-1 if "collapse" in res.events else 0)


def bisect_A(n, h, X, lo, hi, tolA):
    f_lo, f_hi = fate(fate_run(n, h, lo, X)), fate(fate_run(n, h, hi, X))
    if not (f_lo == -1 and f_hi == 1):
        raise RuntimeError(f"bracket [{lo}, {hi}] does not straddle the threshold (fates {f_lo}, {f_hi})")
    k = 2
    while hi - lo > tolA:
        mid = 0.5 * (lo + hi)
        f = fate(fate_run(n, h, mid, X))
        k += 1
        if f == 1:
            hi = mid
        elif f == -1:
            lo = mid
        else:
            raise RuntimeError(f"static through X = {X} at A = {mid}")
    return lo, hi, k


def fit_growth(res, W, window):
    """Slope of ln|H| against tau (tau * sqrt2 if the time_unit mutation is on) over the points after the global minimum of |H|
    with |H|/W in the window."""
    H, tau = res.diag["H"], res.diag["tau"]
    ev = next(iter(res.events.values()), None)
    if ev is not None:                              # the series ends at the event, whichever kernel ran
        H, tau = H[: ev[0] + 1], tau[: ev[0] + 1]
    ok = np.isfinite(H) & (H != 0)
    H, tau = H[ok], tau[ok]
    if dn.mut("time_unit"):
        tau = tau * SQ2
    i0 = int(np.argmin(np.abs(H)))
    lvl = np.abs(H) / W
    m = (np.arange(len(H)) >= i0) & (lvl >= window[0]) & (lvl <= window[1])
    if m.sum() < 10:
        return float("nan"), int(m.sum())
    return float(np.polyfit(tau[m], np.log(np.abs(H[m])), 1)[0]), int(m.sum())


def check_G1(cases=None):
    allok = True
    for n in (cases or (2, 3)):
        W = W_EXACT[n]
        out = []
        for h in CRIT["G1_h"][n]:
            X = CRIT["G1_X"][n]
            t0 = time.time()
            lo, hi, k = bisect_A(n, h, X, -CRIT["G1_bracket"], CRIT["G1_bracket"], CRIT["G1_tolA"])
            sl = []
            for A in (lo, hi):
                res = fate_run(n, h, A, X)
                s, npts = fit_growth(res, W, CRIT["G1_window"])
                sl.append((s, npts))
            Wfit = 0.5 * (sl[0][0] + sl[1][0])
            # untuned run: the departure time and the seed constant
            r0 = fate_run(n, h, 0.0, X)
            Hh, tt = r0.diag["H"], r0.diag["tau"]
            idx = np.where(np.abs(Hh) / W >= 2e-3)[0]
            tdep = float(tt[idx[0]]) if len(idx) else float("nan")
            kappa = 2e-3 * W * math.exp(-W * tdep) / h ** 2 if len(idx) else float("nan")
            out.append(dict(h=h, Astar=0.5 * (lo + hi), runs=k, slopes=[s[0] for s in sl], npts=[s[1] for s in sl], Wfit=Wfit,
                            tdep=tdep, kappa=kappa, secs=time.time() - t0))
            print(f"    G1 D={n + 2} h={h:.5f}: A*_h = {0.5 * (lo + hi):.10e} ({k} runs, {time.time() - t0:.1f} s); slopes {sl[0][0]:.6f}, {sl[1][0]:.6f} "
                  f"({sl[0][1]}, {sl[1][1]} pts) -> W_fit {Wfit:.6f} vs W {W:.6f} (rel {Wfit / W - 1:+.2e}); "
                  f"untuned departure |H| = 2e-3 W at tau = {tdep:.3f}, seed constant kappa = {kappa:.3g}")
        INFO[f"G1_D{n + 2}"] = out
        ok_bar = all(abs(o["Wfit"] / W - 1) <= CRIT["G1_bar"] for o in out)
        ok_agree = abs(out[0]["Wfit"] / out[1]["Wfit"] - 1) <= CRIT["G1_agree"]
        pr = out[0]["Astar"] / out[1]["Astar"]
        dt = out[1]["tdep"] - out[0]["tdep"]
        print(f"    G1 D={n + 2}: A*_h ratio between the two h {pr:.3f} (h^2 seed predicts 4); untuned departure shifts by {dt:.3f} "
              f"(ln4/W = {math.log(4) / W:.3f})")
        check(f"G1[D={n + 2}]", f"tune p to p*_h, then the departure grows at W = {W:.6f} within {CRIT['G1_bar']:.0%}, at two resolutions that agree "
              f"within {CRIT['G1_agree']:.1%}", ok_bar and ok_agree,
              f"W_fit {out[0]['Wfit']:.6f}, {out[1]['Wfit']:.6f}; rel {out[0]['Wfit'] / W - 1:+.2e}, {out[1]['Wfit'] / W - 1:+.2e}")
        allok &= ok_bar and ok_agree


# ----------------------------------------------------------------------------------------------------------------------
def f1_cases():
    h, X, cp = CRIT["F1_h"], CRIT["F1_X"], CRIT["F1_pulse"]
    N = int(round(X / h))
    cases = []
    for name, field, A, eps in (("ghost, +A", 0, cp["A"], (-1.0,)), ("ghost, -A", 0, -cp["A"], (-1.0,)), ("normal, +A", 1, cp["A"], (-1.0, 1.0))):
        d = dn.ellis_rays(2, N, h, pulses=[dict(field=field, ray="+", A=A, xc=cp["xc"], w=cp["w"], kind=cp["kind"])], eps=eps)
        cases.append((name, d))
    return cases


def f1_describe(name, d, res):
    out = dict(name=name)
    ev = res.events
    key = next(iter(ev), None)
    out["fate"] = key or "static"
    if key:
        xp, xm, t = dn.event_time(res, key)
        out["t_event"] = t
    M = dn.mass_nodes(res)
    sol, _q = dn.ellis_dense(d.n, d.N * d.h / SQ2 + 1.0)
    M_static = 1.0 / (2.0 * sol(d.N * d.h / SQ2)[0] ** (d.n - 1))   # the unperturbed Ellis value at the same x+ on the ray
    out["E0"] = float(M[-1, 0] - M_static)         # the change of the Misner-Sharp mass along the ray (Bondi-type, finite X)
    th = dn.throat_min(res)
    dg = res.diag
    out["throat_min_end"] = (float(th[-1, 0]), float(th[-1, 1])) if len(th) else None
    j = np.where(np.isfinite(dg["r"]))[0][-1]
    out["diag_end"] = (float(dg["t"][j]), float(dg["r"][j]), float(dg["H"][j]))
    if key == "collapse":
        j0, i0 = ev["collapse"]
        c, up, um, alive = dn.cell_fields(res)
        # the black hole's mass at the first trapped cell: M = r^(n-1)/2 on r_+ = 0 (D = 4: r/2)
        out["M_first_trapped"] = float(c[0][i0, j0] ** (d.n - 1) / 2)
    return out


def check_F1():
    rows = []
    for name, d in f1_cases():
        res = dn.evolve(d, kernel=KERNEL, tol=CRIT["trap_tol"])
        info = f1_describe(name, d, res)
        rows.append(info)
        print(f"    F1 {name}: fate {info['fate']}" + (f" at t = {info['t_event']:.3f}" if "t_event" in info else "")
              + f"; E0 = {info['E0']:+.4f}; diagonal throat at the last live cell t = {info['diag_end'][0]:.2f} r = {info['diag_end'][1]:.4f} H = {info['diag_end'][2]:+.4f}"
              + (f"; min-r slice at the end t = {info['throat_min_end'][0]:.2f} r = {info['throat_min_end'][1]:.4f}" if info['throat_min_end'] else "")
              + (f"; M at the first trapped cell {info['M_first_trapped']:.4f}" if "M_first_trapped" in info else ""))
    INFO["F1"] = rows
    f = {r["name"]: r["fate"] for r in rows}
    # a qualitative report: gated only on [SH02]'s two statements (a normal pulse collapses; the two ghost signs give different fates)
    ok = f["normal, +A"] == "collapse" and f["ghost, +A"] != f["ghost, -A"] and {f["ghost, +A"], f["ghost, -A"]} == {"collapse", "expand"}
    check("F1", "one pulse of each sign at D = 4: a normal pulse collapses; ghost pulses of the two signs end in a black hole and in an expansion "
          "([SH02]'s two fates). Reported as seen, parameters not tuned", ok, ", ".join(f"{k}: {v}" for k, v in f.items()))


def check_S():
    rows = []
    allok = True
    for name, d in f1_cases():
        if name == "normal, +A":
            continue
        base = dn.evolve(d, kernel=KERNEL, tol=CRIT["trap_tol"])
        rc, fc = CRIT["S_rcut"], CRIT["S_fcut"]
        if dn.mut("cut_extreme"):
            rc, fc = 1.5, 1.5
        var = dn.evolve(d, kernel=KERNEL, tol=CRIT["trap_tol"], rcut=rc, fcut=fc)
        kb, kv = next(iter(base.events), None), next(iter(var.events), None)
        tb = dn.event_time(base, kb)[2] if kb else float("nan")
        tv = dn.event_time(var, kv)[2] if kv else float("nan")
        same = kb == kv and kb is not None
        close = same and abs(tb - tv) <= CRIT["S_dt"] * d.h * SQ2
        allok &= close
        print(f"    S {name}: fate {kb} -> {kv}; event time {tb:.4f} -> {tv:.4f} (|dt| = {abs(tb - tv):.2e}, tolerance {CRIT['S_dt'] * d.h * SQ2:.2e}); "
              f"dead cells {int(base.dead.sum())} -> {int(var.dead.sum())}")
    check("S", "sensitivity: with r_cut halved and F_cut doubled, one classified collapse and one classified expansion keep their fate and their "
          f"event time within {CRIT['S_dt']} h", allok)


CHECKS = dict(I1=check_I1, K=check_K, P1=check_P1, P2=check_P2, P3=check_P3, C=check_C, G1=check_G1, F1=check_F1, S=check_S)
KERNEL = "numba"


def main():
    global KERNEL
    argv = sys.argv[1:]
    only = None
    jpath = None
    for a in argv:
        if a.startswith("--kernel="):
            KERNEL = a.split("=", 1)[1]
        elif a.startswith("--only="):
            only = a.split("=", 1)[1].split(",")
        elif a.startswith("--mutate="):
            dn.MUT = a.split("=", 1)[1]
        elif a.startswith("--json="):
            jpath = a.split("=", 1)[1]
        elif a == "--quick":
            CRIT["P3_h"] = (0.08, 0.04, 0.02)
            CRIT["C_h"] = (0.08, 0.04, 0.02)
            CRIT["G1_h"] = {2: (0.04, 0.02), 3: (0.02 / SQ2, 0.01 / SQ2)}
        elif a.startswith("--crit="):
            k, v = a[7:].split("=", 1)
            CRIT[k] = eval(v)
        else:
            sys.exit(f"unknown argument {a}")
    if dn.MUT is not None and dn.MUT not in MUTATIONS:
        sys.exit(f"unknown mutation {dn.MUT}; known: {', '.join(MUTATIONS)}")
    if KERNEL == "numba" and not dn.HAVE_NUMBA:
        sys.exit("numba is not available; run with --kernel=numpy")
    t0 = time.time()
    import platform
    la = open("/proc/loadavg").read().split()[:3] if platform.system() == "Linux" else ["?"]
    nbv = dn.nb.__version__ if dn.HAVE_NUMBA else "absent"
    print(f"checks.py: kernel {KERNEL}, numba {nbv}, numpy {np.__version__}, python {platform.python_version()}, {platform.node()}, "
          f"load average at start {' '.join(la)}, mutation {dn.MUT}, only {only}, {time.strftime('%Y-%m-%d %H:%M:%S')}")
    for name, fn in CHECKS.items():
        if only and name not in only:
            continue
        try:
            fn()
        except Exception as e:                       # an exception inside a check is that check's failure
            import traceback
            traceback.print_exc()
            check(name, "check raised an exception", False, repr(e))
    bad = [c for c, ok in results if ok is False]
    skipped = [c for c, ok in results if ok is None]
    print(f"\n{sum(1 for _, ok in results if ok)} of {len(results)} checks pass; {len(bad)} fail" + (f": {', '.join(bad)}" if bad else "")
          + (f"; skipped: {', '.join(skipped)}" if skipped else "") + f"; {time.time() - t0:.1f} s")
    if dn.MUT:
        want = set(MUTATIONS[dn.MUT])
        fired = {c.split("[")[0] for c in bad}
        print(f"mutation {dn.MUT}: expected to fire {sorted(want)}; fired {sorted(fired)}; "
              + ("targets all fired" if want <= fired else f"SILENT: {sorted(want - fired)}"))
    if jpath:
        json.dump(dict(info=INFO, results=results, crit={k: str(v) for k, v in CRIT.items()}), open(jpath, "w"), indent=1, default=float)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
