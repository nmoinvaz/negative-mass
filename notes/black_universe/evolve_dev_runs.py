#!/usr/bin/env python3
"""Development runs of #217 (design stage), on the prototype evolve_dev.py. Every number quoted in evolution_design.md comes from one of these.
Not the build and not a run of record.

    evolve_dev_runs.py bitid                 prototype kernel against dn.evolve at c = 0 (V flag off and on): max |G - G_dn|
    evolve_dev_runs.py static                static members against the exact solution at h = 0.04, 0.02, 0.01, with orders
    evolve_dev_runs.py rates                 growth rate of Y_inv (= dr - (R'/phi0') dphi, gauge invariant) from truncation seeds, against #212's W
    evolve_dev_runs.py pulse                 one-sided pulse scans: finder events at c = -0.3 and -1.5, both rays
    evolve_dev_runs.py mode C RHOC H X BR [Y0]     mode family (extracted from the code's own near-threshold run): bisection, series and fits by departure
    evolve_dev_runs.py timing                full-run times at several N
    evolve_dev_runs.py data0                 the c = 0 member's rays against dn.ellis_rays
"""
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "wormhole_evolution"))
import dn                      # noqa: E402
import evolve_dev as ed        # noqa: E402
import evolve_potential as ep  # noqa: E402

SQ2 = math.sqrt(2.0)
W212 = {-0.3: 0.9700800261, -1.5: 0.4941711675, -1.01: 0.6331476235, -0.99: 0.6401216241, -1.0: None}   # notes/black_universe/README.md, Sec. 4
KAPPA = {c: ep.Member(c).kappa for c in (-0.3, -1.5, -1.0, -1.01, -0.99)}


def bitid():
    N, h = 300, 0.02
    d = dn.ellis_rays(2, N, h)
    ref = dn.evolve(d, store=True, stop=False)
    mem = ep.Member(0.0)
    for V in (False, True):
        r = ed.evolve(d, ed.params(mem, V=V), stop=False, store=True)
        print(f"V flag {V}: max |G - G_dn| = {np.abs(r['G'] - ref.G).max():.1e}, dead masks equal {bool((r['dead'] == ref.dead).all())}")


def exact_nodes(mem, rc, N, h):
    ks = np.arange(-N, N + 1)
    rho, A, f = ed.tortoise_rho(mem, rc, ks * h / SQ2)
    i = np.arange(N + 1)[:, None]
    j = np.arange(N + 1)[None, :]
    kk = (i - j) + N
    return rho, f, kk, (i + j) * h / SQ2 + 0 * kk


def static():
    for c, rc in ((0.0, 0.0), (-0.3, 0.0), (-1.5, 1.0), (-1.0, 1.0)):
        mem = ep.Member(c)
        print(f"== c = {c}, corner rho_c = {rc}, X = 10")
        prev = None
        for h in (0.04, 0.02, 0.01):
            N = int(round(10.0 / h))
            d = ed.member_rays(mem, rc, N, h)
            t0 = time.time()
            res = ed.evolve(d, ed.params(mem), stop=False, store=True)
            el = time.time() - t0
            rho, f, kk, t = exact_nodes(mem, rc, N, h)
            G = res["G"]
            e = np.maximum(np.maximum(np.abs(G[0] - mem.R(rho)[kk]), np.abs(G[1] - f[kk])), np.abs(G[2] - mem.phi(rho)[kk]))
            row = []
            for tc in (2.0, 4.0, 8.0):
                row.append(float(e[t <= tc].max()))
            print(f"   h = {h}: max error (r, f, phi) for t <= 2, 4, 8: {row[0]:.3e} {row[1]:.3e} {row[2]:.3e}   {el:.2f} s")
            if prev is not None:
                print(f"           ratio to the previous h: {prev[0] / row[0]:.2f} {prev[1] / row[1]:.2f} {prev[2] / row[2]:.2f}")
            prev = row


def rates():
    h, X = 0.01, 24.0
    print(f"h = {h}, X = {X}; slope of ln|Y_inv| against t over t in [8, 14] at the cell l = l0, from the truncation seed alone (no pulse)")
    for c, rc in ((-0.3, 0.0), (-1.5, 1.0), (-1.0, 1.0), (-1.01, 1.0), (-0.99, 1.0)):
        mem = ep.Member(c)
        N = int(round(X / h))
        d = ed.member_rays(mem, rc, N, h)
        res = ed.evolve(d, ed.params(mem), stop=False, store=True)
        rho, f, kk, t = exact_nodes(mem, rc, N, h)
        G = res["G"]
        dr, dph = G[0] - mem.R(rho)[kk], G[2] - mem.phi(rho)[kk]
        rr = rho[kk]
        Y = dr - (rr * mem.R(rr) / SQ2) * dph
        i = np.arange(N + 1)[:, None]
        j = np.arange(N + 1)[None, :]
        l = (i - j) * h / SQ2 + 0 * kk
        out = []
        for l0 in (-2.0, -0.5, 1.0):
            m = np.abs(l - l0) < 0.6 * h / SQ2 + 1e-12
            tt, yy = t[m], Y[m]
            sel = (tt >= 8) & (tt <= 14)
            out.append(np.polyfit(tt[sel], np.log(np.abs(yy[sel])), 1)[0])
        w = W212[c]
        print(f"   c = {c:+.2f}: Y_inv rate at l = -2, -0.5, +1: {out[0]:.3f} {out[1]:.3f} {out[2]:.3f};  #212 W = {w}  kappa = {KAPPA[c]:.4f}")


def pulse():
    h, X = 0.02, 30.0
    for c, rc in ((-0.3, 0.0), (-1.5, 1.0)):
        mem = ep.Member(c)
        N = int(round(X / h))
        for ray in ("+", "-"):
            print(f"== c = {c}, corner rho_c = {rc}, pulse on the {ray} ray (centre 3, width 0.5), h = {h}, X = {X}: raw finders (tolerance 1e-2 on r_+, r_-)")
            for A in (-0.2, -0.05, -0.01, -0.005, 0.005, 0.01, 0.05, 0.2):
                d = ed.member_rays(mem, rc, N, h, pulse=dict(A=A, xc=3.0, w=0.5, ray=ray))
                res = ed.evolve(d, ed.params(mem), stop=True, tol=1e-2, rexp=2.0)
                ev = res["events"]
                s = "none through t = %.1f" % (SQ2 * X) if not ev else "%s at t = %.2f, l = %+.2f" % ((list(ev)[0],) + ed.event_T(d, ev)[1:])
                print(f"   A = {A:+.3f}: {s}")


def bisect(run, lo, hi, tol=1e-14):
    rl, rh = run(lo), run(hi)
    if rl["fate"] * rh["fate"] != -1:
        raise RuntimeError(f"bracket does not straddle: {rl} {rh}")
    fl, k = rl["fate"], 2
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        r = run(mid)
        k += 1
        if r["fate"] == fl:
            lo = mid
        elif r["fate"] == -fl:
            hi = mid
        else:
            raise RuntimeError(f"no departure at {mid!r}")
    return lo, hi, fl, k


def fit(rows, ps, win, sides=(-1, 1)):
    pts = [(r["side"], abs(r["A"] - ps), r["T"]) for r in rows
           if r["T"] is not None and r["side"] in sides and r["fate"] == -r["side"] * rows[0]["orient"] and win[0] * (1 - 1e-6) <= abs(r["A"] - ps) <= win[1] * (1 + 1e-6)]
    present = sorted({s for s, _, _ in pts})
    M = np.array([[1.0 if s == sd else 0.0 for sd in present] + [math.log(1 / d)] for s, d, _ in pts])
    y = np.array([T for _, _, T in pts])
    coef, *_ = np.linalg.lstsq(M, y, rcond=None)
    res = y - M @ coef
    cov = (res @ res / max(1, len(y) - M.shape[1])) * np.linalg.inv(M.T @ M)
    return float(coef[-1]), float(math.sqrt(cov[-1, -1])), len(y)


def mode(c, rc, h, X, br, Y0=3e-3, Ap=0.0, Yext=2e-3):
    mem = ep.Member(c)
    t0 = time.time()
    dR, dL, info = ed.extract_mode(mem, rc, h, X + 20.0, Ap, Y_ext=Yext)
    print(f"c = {c}, corner rho_c = {rc}, h = {h}, X = {X}: mode extracted at t_s = {info['t_s']:.2f} from the run with A_p = {Ap} (row {info['js']}), Y at the corner node {info['Yc']:.3e}; "
          f"largest |dr|, |df|, |dphi| over the rays per unit Y at the corner: {(info['maxdR'] / abs(info['Yc'])).round(1).tolist()} (right ray), {(info['maxdL'] / abs(info['Yc'])).round(1).tolist()} (left ray)")
    N = int(round(X / h))
    par, st = ed.params(mem), ed.Static(mem, rc, N, h)

    def run(A):
        d = ed.mode_rays(mem, rc, N, h, A, dR, dL)
        r = ed.evolve_dep(d, par, st, Y0=Y0)
        dp = r["dep"]
        return dict(A=A, fate=0 if dp is None else dp["sign"], T=None if dp is None else dp["T"])
    lo, hi, fl, k = bisect(run, -br, br)
    ps = 0.5 * (lo + hi)
    print(f"   A*_h = {ps:.6e} ({k} runs, {time.time() - t0:.1f} s with the extraction); Y_inv departure at |Y| = {Y0} on the diagonal; fate(A < A*) = {fl:+d}")
    rows = []
    for e in np.arange(max(2.5, math.floor(-math.log10(Y0) * 2) / 2), 11.51, 0.5):
        for s in (-1, 1):
            r = run(ps + s * 10.0 ** (-e))
            r.update(e=e, side=s, orient=fl)
            rows.append(r)
    for sd in (-1, 1):
        print(f"   side {sd:+d} (A = A* {'-' if sd < 0 else '+'} 10^-e): " + " ".join(f"{r['e']:.1f}:{r['T']:.2f}" for r in rows if r["side"] == sd and r["T"] is not None))
    w = W212[c] or KAPPA[c]                        # c = -1 has no W in #212: the comparison is with kappa (label below)
    lab = "1/W" if W212[c] else "1/kappa (no growing mode in #212)"
    hiw = Y0 / 2
    for win in [w for w in ((1e-11, hiw), (1e-10, hiw), (1e-9, hiw), (1e-8, hiw), (1e-7, hiw), (1e-6, hiw), (1e-9, 1e-4), (1e-8, 1e-4), (1e-7, 1e-3)) if w[0] < w[1] / 3]:
        try:
            s, se, n = fit(rows, ps, win)
            print(f"   window {win}: slope {s:.5f} +- {se:.5f} (n = {n}); {lab} = {1 / w:.5f}; relative {s * w - 1:+.2e}")
        except Exception as ex:                      # noqa: BLE001
            print("   window", win, "failed:", ex)


def timing():
    c, rc = -1.5, 1.0
    mem = ep.Member(c)
    for h, X in ((0.02, 40.0), (0.01, 40.0), (0.005, 40.0), (0.0025, 40.0)):
        N = int(round(X / h))
        d = ed.member_rays(mem, rc, N, h)
        st = ed.Static(mem, rc, N, h)
        t0 = time.time()
        r = ed.evolve_dep(d, ed.params(mem), st, Y0=1e9, stop_dep=False)
        print(f"c = {c}: h = {h}, X = {X}, N = {N}: full run (no stop) {time.time() - t0:.2f} s, rows {r['rows']}")


def data0():
    """member_rays at c = 0 (no pulse) against dn.ellis_rays: r, f, phi on both rays"""
    N, h = 1000, 0.02
    d1 = dn.ellis_rays(2, N, h)
    d2 = ed.member_rays(ep.Member(0.0), 0.0, N, h)
    for nm, a, b in (("+ ray", d1.Sp, d2.Sp), ("- ray", d1.Sm, d2.Sm)):
        print(f"{nm}: max |r| difference {np.abs(a[0] - b[0]).max():.2e}, |f| {np.abs(a[1] - b[1]).max():.2e}, |phi| {np.abs(a[2] - b[2]).max():.2e}")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "mode":
        mode(float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), *(float(x) for x in sys.argv[7:]))
    else:
        globals()[cmd]()
