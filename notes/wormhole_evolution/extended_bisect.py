#!/usr/bin/env python3
"""Bisection of the Ellis threshold in binary128: is the roundoff floor of #40's series the evolution's? (#227)

#40's double-null threshold search (family S, D = 4, h = 0.02) repeated with every number in IEEE binary128 (np.longdouble on
aarch64 Linux, 113-bit mantissa): the initial data (static Ellis wormhole from a Taylor series, Gaussian pulses, the RK4 constraint
solve of dn.ellis_rays), the evolution (dn.cell_np, the numpy kernel that is bit-identical to the numba row kernel in float64) and the
fates. dn.py, checks.py, threshold.py and constraint_plateau.py are imported and not changed. Criteria, planted cases and mutations:
negative-mass-raw I227-rpi5/README.md, committed before the runs of record.

    python3 extended_bisect.py controls [--mutate=NAME]          K1-K4 and the planted cases PL1-PL4 (no evolution data)
    python3 extended_bisect.py point --D=4 --h=0.02 --A=P        one fate run in binary128 (cost)
    python3 extended_bisect.py bisect --D=4 --h=0.02 --out=DIR   bisection + series (binary128), json
    python3 extended_bisect.py f64 --D=4 --h=0.02 --out=DIR      float64 bisections from K brackets (threshold.bisect, numba)
    python3 extended_bisect.py stall --D=4 --h=0.02 --out=DIR    float64 - binary128 on the stall
    python3 extended_bisect.py analyze --raw=DIR [--mutate=NAME] every check and table from the stored runs

Status words as in AGENTS.md: every result is *numerical evidence*.
"""
import argparse
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dn           # noqa: E402
import threshold    # noqa: E402
import constraint_plateau as cp   # noqa: E402

LD = np.longdouble
SQ2 = math.sqrt(2.0)
CRIT = threshold.CRIT
MUT = None          # --mutate=NAME of this script; dn.MUT carries dn's and threshold's own names
TOL, REXP, RCUT, FCUT = CRIT["trap_tol"], CRIT["rexp"], CRIT["rcut"], CRIT["fcut"]
LOGF = math.log(FCUT)


def mut(name):
    return MUT == name


def ld_sqrt2():
    return np.sqrt(LD(2))


# --------------------------------------------------------------------------------------------------------------------
# static Ellis wormhole in binary128 (a = 1): r'' = (n-1) r^(1-2n), phi' = q r^(-n), r(0) = 1, r'(0) = 0, phi(0) = 0, q = sqrt(n (n-1))
# --------------------------------------------------------------------------------------------------------------------
def _ipow(x, m):
    out = LD(1)
    for _ in range(m):
        out = out * x
    return out


def static_taylor(n, dl, K, order=18):
    """r, r', phi at l_k = k dl, k = 0..K, by Taylor series of the given order about each l_k (power-series recurrences for r^(1-2n)
    and r^(-n)); the truncation error of a step is (dl / radius)^(order + 1) with radius >= 1: 1e-37 at dl = 0.007, order 18."""
    dl = LD(dl)
    order = 8 if mut("taylor_order") else order
    q = np.sqrt(LD(n) * LD(n - 1))
    a1, a2 = -(2 * n - 1), -n
    R, V, P = np.empty(K + 1, LD), np.empty(K + 1, LD), np.empty(K + 1, LD)
    r0, v0, p0 = LD(1), LD(0), LD(0)
    R[0], V[0], P[0] = r0, v0, p0
    for k in range(K):
        r = [LD(0)] * (order + 1)
        r[0], r[1] = r0, v0
        y = [LD(0)] * order
        z = [LD(0)] * order
        y[0] = LD(1) / _ipow(r0, 2 * n - 1)
        z[0] = LD(1) / _ipow(r0, n)
        for m in range(order):                           # y[m], z[m] need r[0..m]; r[m + 2] needs y[m]
            if m >= 1:
                s1 = LD(0)
                s2 = LD(0)
                for j in range(1, m + 1):
                    s1 += (a1 * j - (m - j)) * r[j] * y[m - j]
                    s2 += (a2 * j - (m - j)) * r[j] * z[m - j]
                y[m] = s1 / (m * r0)
                z[m] = s2 / (m * r0)
            if m + 2 <= order:
                r[m + 2] = (n - 1) * y[m] / ((m + 1) * (m + 2))
        rn = LD(0)
        for m in range(order, -1, -1):
            rn = rn * dl + r[m]
        vn = LD(0)
        for m in range(order, 0, -1):
            vn = vn * dl + m * r[m]
        sp = LD(0)
        for m in range(order - 1, -1, -1):
            sp = sp * dl + q * z[m] / (m + 1)
        r0, v0, p0 = rn, vn, p0 + dl * sp
        R[k + 1], V[k + 1], P[k + 1] = r0, v0, p0
    return R, V, P


_STATIC = {}


def static_grid(n, h, Nmax):
    """r(l), phi(l) on l = m dl, dl = h / (2 sqrt 2) (m = 2i the grid nodes s = i h, m odd the RK4 midpoints), m = 0..2 Nmax; cached."""
    key = (n, float(h))
    have = _STATIC.get(key)
    if have is not None and len(have[0]) >= 2 * Nmax + 1:
        return have
    dl = LD(h) / (2 * ld_sqrt2())
    Nmax = int(Nmax * 1.15) + 200                   # headroom: the row cap grows run to run
    R, V, P = static_taylor(n, dl, 2 * Nmax)
    _STATIC[key] = (R, V, P, dl)
    return _STATIC[key]


# --------------------------------------------------------------------------------------------------------------------
# initial data of dn.ellis_rays in binary128 (family S: ghost pulse on '+' mirrored onto '-')
# --------------------------------------------------------------------------------------------------------------------
def pulse_ld(s, A, xc=3.0, w=0.5):
    """(value, d/ds) of the Gaussian ghost pulse A (g - g(0)), g = exp(-((s - xc)/w)^2), as dn._pulse, in binary128."""
    xc, w = LD(xc), LD(w)
    g = np.exp(-((s - xc) / w) ** 2)
    g0 = np.exp(-((LD(0) - xc) / w) ** 2)
    return A * (g - g0), A * (-LD(2) * (s - xc) / w ** 2) * g


def rays_ld(n, Nc, h, A, eta_scale=None, bump=None, dense=False):
    """dn.ellis_rays(n, Nc, h, pulses=[ghost A on '+'], symmetric=True, eps=(-1,)) in binary128: Sp, Sm (3, Nc + 1) = (r, f = 0, phi).
    eta_scale: the data's pulse amplitude is A (1 + eta) (the label A is not scaled). bump: a relative bump eta exp(-((s - 5)/0.5)^2)
    in r on both rays after the constraint solve (planted perturbations). dense: the static r and phi from dn.ellis_dense (float64, rtol 1e-13)
    cast to binary128, as the float64 run reads them."""
    cp.need_binary128(LD)
    hL = LD(h)
    sq2 = ld_sqrt2()
    R, V, P, dl = static_grid(n, h, Nc)
    q = np.sqrt(LD(n) * LD(n - 1))
    m = np.arange(2 * Nc + 1)
    s = m * (hL / 2)
    if dense:
        sol, qf = dn.ellis_dense(n, Nc * float(h) / SQ2 + 1.0)
        lf = (m * (float(h) / 2)) / SQ2
        Rm, Pm = np.asarray(sol(lf)[0], LD), np.asarray(sol(lf)[2], LD)
        qq = LD(qf)
        sq = LD(SQ2)
    else:
        Rm, Pm, qq, sq = R[:2 * Nc + 1], P[:2 * Nc + 1], q, sq2
    stat_s = qq / Rm ** n / sq                            # d phi_static / ds on the '+' ray
    Aeff = A * (LD(1) + LD(eta_scale)) if eta_scale else A
    nodes = np.arange(Nc + 1) * hL
    rays = []
    for ray, sgn in (("+", 1), ("-", -1)):
        pv, pd = pulse_ld(s, sgn * Aeff)                  # the mirrored pulse of field 0 has the opposite sign
        phis = sgn * stat_s + pd
        S = -(phis ** 2)                                  # eps = -1
        S0, Sh, S1 = S[0:2 * Nc:2], S[1:2 * Nc:2], S[2:2 * Nc + 1:2]
        r = np.empty(Nc + 1, LD)
        y0, v0 = LD(1), LD(0)
        r[0] = y0
        half = LD(0.5)
        for i in range(Nc):
            a0, b0 = y0, v0
            k1 = (b0, -(S0[i] / n) * a0)
            a2, b2 = a0 + half * hL * k1[0], b0 + half * hL * k1[1]
            k2 = (b2, -(Sh[i] / n) * a2)
            a3, b3 = a0 + half * hL * k2[0], b0 + half * hL * k2[1]
            k3 = (b3, -(Sh[i] / n) * a3)
            a4, b4 = a0 + hL * k3[0], b0 + hL * k3[1]
            k4 = (b4, -(S1[i] / n) * a4)
            y0 = a0 + hL / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            v0 = b0 + hL / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            r[i + 1] = y0
        if bump:
            r = r * (LD(1) + LD(bump) * np.exp(-((nodes - LD(5)) / LD(0.5)) ** 2))
        U = np.zeros((3, Nc + 1), LD)
        U[0] = r
        U[2] = (sgn * Pm[0:2 * Nc + 1:2]) + pv[0:2 * Nc + 1:2]
        rays.append(U)
    return dn.Data(n, [-1.0], Nc, float(h), rays[0], rays[1], name="ellis_ld")


# --------------------------------------------------------------------------------------------------------------------
# the evolution: anti-diagonal sweep of dn.cell_np with the events of dn.evolve (stop=True)
# --------------------------------------------------------------------------------------------------------------------
def sweep(d, Jcap, dtype=LD, noise=0.0, seed=0, snap_every=0, find_events=True, hmut=1.0, smax=None):
    """Evolve the data d on the rectangle i <= d.N (x+), j <= Jcap (x-) by anti-diagonals i + j = s (three in memory), with dn.cell_np in
    `dtype` (float64: the arithmetic of the numba row kernel; np.longdouble: binary128). Events as dn.evolve(stop=True): per row j the first
    trapped cell (both r_+ and r_- below -tol; column smallest) and the diagonal cell (j, j) anti-trapped with r >= rexp; the first row
    with either wins, a trapped cell wins a tie. Stops when the first event row is complete. noise: relative uniform noise of that amplitude
    on every node as it is computed (the roundoff mutation). snap_every: keep r on every snap_every-th node diagonal; smax: last diagonal.
    Returns dict(fate, j, i, T, rows, snaps)."""
    Nc, h, nf, eps, n = d.N, d.h, d.nf, d.eps, d.n
    par = dn.kernel_params(n, eps, 1)
    hd = dtype(h * hmut) if dtype is np.float64 else LD(h) * LD(hmut)
    Sp, Sm = d.Sp.astype(dtype), d.Sm.astype(dtype)
    rng = np.random.default_rng(seed) if noise else None
    INF = 1 << 60
    ctrap = np.full(Jcap + 1, INF, np.int64)
    cexp = np.zeros(Jcap + 1, bool)
    D2 = None
    D1 = np.zeros((nf, Nc + 1), dtype)
    D1[:, 0] = Sp[:, 0]                                   # s = 0: the corner node
    dd2 = None
    dd1 = np.zeros(Nc + 1, bool)
    snaps = {}
    best = INF
    two_h = dtype(2) * hd
    quarter = dtype(0.25)
    sf = -1.0 if dn.mut("finder_sign") else 1.0
    for s in range(1, (Nc + Jcap + 2) if smax is None else min(Nc + Jcap + 2, smax + 1)):
        D0 = np.zeros((nf, Nc + 1), dtype)
        dd0 = np.zeros(Nc + 1, bool)
        i0, i1 = max(1, s - Jcap), min(Nc, s - 1)
        if i0 <= i1 and D2 is not None:
            Sv = D1[:, i0 - 1:i1] if mut("diag_S") else D2[:, i0 - 1:i1]
            Ev = D1[:, i0:i1 + 1]
            Wv = D1[:, i0 - 1:i1]
            ds = dd2[i0 - 1:i1] | dd1[i0:i1 + 1] | dd1[i0 - 1:i1]
            with np.errstate(all="ignore"):
                out = dn.cell_np(Sv, Ev, Wv, hd, par, eps)
                if noise:
                    out = out * (1.0 + noise * rng.uniform(-1.0, 1.0, out.shape))
                bad = ~np.all(np.isfinite(out), axis=0)
                ds = ds | bad | (out[0] < RCUT) | (out[1] < -LOGF)
            D0[:, i0:i1 + 1] = np.where(ds[None, :], Ev, out)
            dd0[i0:i1 + 1] = ds
        if s <= Nc:
            D0[:, s] = Sp[:, s]
        if s <= Jcap:
            D0[:, 0] = Sm[:, s]
        if snap_every and s % snap_every == 0:
            snaps[s] = D0[0].copy()
        c = s - 2
        if find_events and c >= 0 and D2 is not None:
            ic0, ic1 = max(0, c - (Jcap - 1)), min(Nc - 1, c)
            if ic0 <= ic1:
                Sr, Er, Wr, Nr = D2[0, ic0:ic1 + 1], D1[0, ic0 + 1:ic1 + 2], D1[0, ic0:ic1 + 1], D0[0, ic0 + 1:ic1 + 2]
                up = (Er + Nr - Sr - Wr) / two_h
                um = (Wr + Nr - Sr - Er) / two_h
                alive = ~(dd2[ic0:ic1 + 1] | dd1[ic0 + 1:ic1 + 2] | dd1[ic0:ic1 + 1] | dd0[ic0 + 1:ic1 + 2])
                cc = np.where((sf * up > TOL) & (sf * um > TOL), 1, np.where((sf * up < -TOL) & (sf * um < -TOL), -1, 0))
                ii = np.arange(ic0, ic1 + 1)
                jj = c - ii
                tr = alive & (cc == -1)
                if tr.any():
                    ctrap[jj[tr]] = np.minimum(ctrap[jj[tr]], ii[tr])
                    best = min(best, int(jj[tr].min()))
                if c % 2 == 0 and ic0 <= c // 2 <= ic1:
                    k = c // 2 - ic0
                    cen = quarter * (Sr[k] + Er[k] + Wr[k] + Nr[k])
                    if alive[k] and cc[k] == 1 and cen >= REXP:
                        cexp[c // 2] = True
                        best = min(best, c // 2)
        D2, D1 = D1, D0
        dd2, dd1 = dd1, dd0
        if find_events and best < INF and best <= s - Nc - 1:   # rows <= s - Nc - 1 are complete (cell (Nc - 1, j) is on the cell diagonal Nc - 1 + j = s - 2)
            break
    ev = None
    for j in range(Jcap + 1):
        if ctrap[j] < INF:
            ev = ("collapse", j, int(ctrap[j]))
            break
        if cexp[j]:
            ev = ("expand", j, j)
            break
    out = dict(fate=0, T=None, j=None, i=None, snaps=snaps, rows=best if best < INF else None)
    if find_events and ev is not None:
        key, j, i = ev
        out.update(fate=-1 if key == "collapse" else 1, j=int(j), i=int(i), T=((i + 0.5) * h + (j + 0.5) * h) / SQ2)
    return out


# --------------------------------------------------------------------------------------------------------------------
# fate runs of one problem, bisection and series
# --------------------------------------------------------------------------------------------------------------------
class Runner:
    """Fate runs of family S at (D, h) in binary128 (variant E: the data from the exact static solution; F: the data from dn.ellis_dense, as the
    float64 run reads it). The row cap Jcap follows the longest lifetime seen (+ margin) and doubles-up (x1.4) when a run finds no event.
    noise: relative uniform noise a * 2^-52 (a = noise) on every node (the roundoff mutation); eta_scale and bump: planted data perturbations."""

    def __init__(self, n, h, variant="E", noise=0.0, seed=0, eta_scale=None, bump=None, wcol=40, tmax0=8.0, jmax=5000):
        self.n, self.h, self.variant = n, h * (1.01 if mut("h_wrong") else 1.0), variant       # h_wrong: the whole problem at 1.01 h
        self.noise = (1.0 if mut("noise_f64") else noise) * 2.0 ** -52
        self.seed, self.eta_scale, self.bump, self.wcol, self.jmax = seed, eta_scale, bump, wcol, jmax
        self.tmax = tmax0
        self.count = 0
        self.secs = 0.0
        self.off = 0                                  # max |i - j| of an event cell

    def __call__(self, p):
        t0 = time.time()
        p = LD(p)
        hmut = 1.01 if mut("quad_h") else 1.0
        while True:
            Jc = min(self.jmax, int(1.12 * (self.tmax + 1.5) / (SQ2 * self.h)) + 40)
            Nc = Jc + self.wcol
            if mut("float_cast"):                     # the amplitude and the data's arithmetic fall back to float64
                p = LD(float(p))
            d = rays_ld(self.n, Nc, self.h, p, eta_scale=None if mut("eta_dropped") else self.eta_scale, bump=self.bump, dense=(self.variant == "F"))
            o = sweep(d, Jc, LD, noise=self.noise, seed=self.seed + self.count, hmut=hmut)
            self.count += 1
            if o["fate"] != 0 or Jc >= self.jmax:
                break
            self.tmax *= 1.4
        if o["fate"]:
            self.tmax = max(self.tmax, o["T"])
            self.off = max(self.off, abs(o["i"] - o["j"]))
        dt = time.time() - t0
        self.secs += dt
        return dict(A=p, fate=o["fate"], T=o["T"], j=o["j"], i=o["i"], secs=round(dt, 2), Jcap=Jc)


class LDSink:
    def __init__(self, path):
        self.f = open(path, "a")

    def __call__(self, rec):
        rec = dict(rec)
        if "A" in rec:
            rec["A"] = str(rec["A"])
        self.f.write(json.dumps(rec) + "\n")
        self.f.flush()


def ld_series(run, pstar, es, sides, sink=None):
    """Runs at p = p* + s 10^-e (binary128 arithmetic) for e in es, s in sides; rows as threshold.series."""
    rows = []
    for e in es:
        dl = LD(10) ** (-LD(e))
        for s in sides:
            r = run(pstar + s * dl)
            r.update(phase="series", side=s, e=float(e))
            if sink:
                sink(dict(r))
            rows.append(r)
    return rows


def pstar_of(lo, hi):
    return LD(0.5) * (LD(lo) + LD(hi))


def cmd_bisect(a):
    n = a.D - 2
    os.makedirs(a.out, exist_ok=True)
    sink = LDSink(os.path.join(a.out, "points_bisect.jsonl"))
    run = Runner(n, a.h, a.variant, noise=a.noise, eta_scale=a.eta_scale, bump=a.bump)
    lo, hi = LD(a.lo), LD(a.hi)
    t0 = time.time()
    lo, hi, stopped = threshold.bisect(run, lo, hi, LD(a.tol), sink)
    summ = dict(D=a.D, h=a.h, variant=a.variant, noise=a.noise, eta_scale=a.eta_scale, bump=a.bump, lo=str(lo), hi=str(hi), tol=a.tol,
                stopped=stopped, runs=run.count, secs=round(time.time() - t0, 1), off=run.off, tmax=run.tmax, mutation=MUT or dn.MUT)
    json.dump(summ, open(os.path.join(a.out, "bisect.json"), "w"), indent=1)
    print(json.dumps(summ), flush=True)
    return 0


def cmd_series(a):
    n = a.D - 2
    b = json.load(open(os.path.join(a.out, "bisect.json")))
    ps = pstar_of(b["lo"], b["hi"])
    es = [float(x) for x in a.elist.split(",")] if a.elist else [a.e0 + 0.5 * k for k in range(int(round((a.e1 - a.e0) / 0.5)) + 1)]
    sides = (-1, 1) if a.side == 0 else (a.side,)
    tag = {0: "both", -1: "m", 1: "p"}[a.side]
    sink = LDSink(os.path.join(a.out, f"series_{tag}_{a.e0:g}_{a.e1:g}.jsonl"))
    run = Runner(n, a.h, b["variant"], noise=b["noise"], eta_scale=b["eta_scale"], bump=b["bump"], tmax0=a.tmax0)
    t0 = time.time()
    ld_series(run, ps, es, sides, sink)
    print(json.dumps(dict(series=tag, e0=a.e0, e1=a.e1, runs=run.count, secs=round(time.time() - t0, 1), off=run.off)), flush=True)
    return 0


# --------------------------------------------------------------------------------------------------------------------
# float64: bisections from several brackets (threshold.bisect and threshold.point, numba), and the series at the binary128 p*
# --------------------------------------------------------------------------------------------------------------------
def bracket_k(k):
    """Bracket k of the K float64 bisections: k = 0 is #40's [-0.01, 0.01]; the others are inside it, different, so the bisections visit
    different amplitudes and so different roundoff."""
    return -0.01 * (1.0 - 0.031 * k), 0.01 * (1.0 - 0.023 * k)


def cmd_f64(a):
    n = a.D - 2
    os.makedirs(a.out, exist_ok=True)
    X = a.X or CRIT["X"]["S"][n]
    res = []
    for k in range(a.k0, a.k1):
        lo, hi = bracket_k(k)
        sink = threshold.Sink(os.path.join(a.out, f"f64_points_k{k}.jsonl"))
        run = lambda p: threshold.point(n, "S", a.h, X, p)
        t0 = time.time()
        lo, hi, stopped = threshold.bisect(run, lo, hi, a.f64tol, sink)
        rec = dict(k=k, lo=lo, hi=hi, pstar=0.5 * (lo + hi), stopped=stopped, secs=round(time.time() - t0, 1), X=X, h=a.h, D=a.D, tol=a.f64tol)
        json.dump(rec, open(os.path.join(a.out, f"f64_bisect_k{k}.json"), "w"), indent=1)
        print(json.dumps(rec), flush=True)
    return 0


def cmd_f64series(a):
    """float64 runs (numba) at p = float64(p*_binary128 + s 10^-e), e = e0 .. e1."""
    n = a.D - 2
    b = json.load(open(os.path.join(a.out, "bisect.json")))
    ps = pstar_of(b["lo"], b["hi"])
    X = a.X or 34.0
    sink = threshold.Sink(os.path.join(a.out, f"f64series_{b['variant']}.jsonl"))
    es = [float(x) for x in a.elist.split(",")] if a.elist else [a.e0 + 0.5 * k for k in range(int(round((a.e1 - a.e0) / 0.5)) + 1)]
    for e in es:
        dl = LD(10) ** (-LD(e))
        for s in (-1, 1):
            p = float(ps + s * dl)
            r = threshold.point(n, "S", a.h, X, p)
            r.update(phase="f64series", side=s, e=e, pLD=str(ps + s * dl), A=repr(p))
            sink(r)
    return 0


# --------------------------------------------------------------------------------------------------------------------
# the stall: float64 minus binary128 at the same p, and the pair that gives the unstable mode's profile
# --------------------------------------------------------------------------------------------------------------------
def stall_runs(n, h, pc, delta, smax, snap_every, variant="F"):
    """For p = pc -/+ delta (rounded to float64, the same number in both precisions): the float64 run (dn.ellis_rays, the numpy sweep) and
    the binary128 run (variant F: the data's static solution from dn.ellis_dense, as float64 reads it), both to node diagonal smax with r kept on
    every snap_every-th diagonal. Returns {side: dict(p, r64 = {s: array}, r128 = {s: array})}."""
    out = {}
    Nc = smax
    for side in (-1, 1):
        pf = float(LD(pc) + side * LD(delta))
        pulses, sym, eps = threshold.family("S", pf)
        d64 = dn.ellis_rays(n, Nc, h, pulses=pulses, symmetric=sym, eps=eps)
        o64 = sweep(d64, Nc, np.float64, snap_every=snap_every, find_events=False, smax=smax)
        dld = rays_ld(n, Nc, h, LD(pf), dense=(variant == "F"))
        old = sweep(dld, Nc, LD, snap_every=snap_every, find_events=False, smax=smax)
        out[side] = dict(p=pf, r64=o64["snaps"], r128=old["snaps"])
    return out


def stall_profile(n, h, runs, L0):
    """Per diagonal s (t = s h / sqrt2) over the nodes with |l| = |2 i - s| h / sqrt2 <= L0: rms and max of d = r64 - r128 on each side, the
    pair's profile Q = (r128(+) - r128(-)) / (p(+) - p(-)) and the projection eps = <d, Q> / <Q, Q> (units of p) with the residual fraction."""
    rows = []
    sm, sp = runs[-1], runs[1]
    dp = LD(sp["p"]) - LD(sm["p"])
    for s in sorted(sm["r128"]):
        if s < 8:
            continue
        i = np.arange(0, s + 1)
        l = (2 * i - s) * h / SQ2
        sel = np.abs(l) <= L0
        Q = ((sp["r128"][s][:s + 1] - sm["r128"][s][:s + 1]) / dp)[sel].astype(np.float64)
        rec = dict(s=s, t=s * h / SQ2, nsel=int(sel.sum()), qrms=float(np.sqrt(np.mean(Q ** 2))))
        for side, tag in ((-1, "m"), (1, "p")):
            dlt = (LD(runs[side]["r64"][s][:s + 1].astype(np.float64)) - runs[side]["r128"][s][:s + 1])[sel].astype(np.float64)
            qq = float(Q @ Q)
            eps = float(dlt @ Q) / qq if qq > 0 else float("nan")
            res = dlt - eps * Q
            rec[tag] = dict(rms=float(np.sqrt(np.mean(dlt ** 2))), max=float(np.max(np.abs(dlt))), eps=eps,
                            resid=float(np.sqrt(np.mean(res ** 2)) / max(np.sqrt(np.mean(dlt ** 2)), 1e-300)))
        rows.append(rec)
    return rows


def cmd_stall(a):
    n = a.D - 2
    os.makedirs(a.out, exist_ok=True)
    smax = int(math.floor(a.tmax * SQ2 / a.h))
    snap = max(1, int(round(a.snap * SQ2 / a.h)))
    pc = LD(a.pc)
    allrows = []
    for k in range(a.k0, a.k1):
        dl = LD(a.delta) * (LD(1) + LD("0.0731") * k)
        t0 = time.time()
        runs = stall_runs(n, a.h, pc, dl, smax, snap, a.variant)
        rows = stall_profile(n, a.h, runs, a.L0)
        rec = dict(k=k, pc=str(pc), delta=str(dl), p=[runs[-1]["p"], runs[1]["p"]], variant=a.variant, D=a.D, h=a.h, L0=a.L0, smax=smax, snap=snap,
                   secs=round(time.time() - t0, 1), rows=rows)
        json.dump(rec, open(os.path.join(a.out, f"stall_{a.variant}_k{k}.json"), "w"))
        print(json.dumps(dict(k=k, secs=rec["secs"], p=rec["p"])), flush=True)
    return 0


# --------------------------------------------------------------------------------------------------------------------
# checks, controls (no evolution data of record) and planted cases
# --------------------------------------------------------------------------------------------------------------------
RESULTS = []


ONLY = None         # analyze --only=A2,P1: only the checks whose id starts with one of these (views of the run-level mutations hold part of the data)


def check(cid, text, ok, detail=""):
    if ONLY is not None and cid.split("-")[0] not in ONLY:
        return
    ok = None if ok is None else bool(ok)
    RESULTS.append((cid, ok))
    tag = "SKIP" if ok is None else ("pass" if ok else "FAIL")
    print(f"[{tag}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


W = threshold.W_EXACT
QLAT = SQ2 * 0.02          # the row spacing of the event time at h = 0.02

PLANT = dict(p0="4.7527725350123456789012345678901234e-05", cm=-0.79, cp=4.05, kappa=0.37, floor_rms=3e-14, seeds=8)


def prep(p):
    """The amplitude as the run uses it (mutation float_cast: through float64)."""
    return LD(float(p)) if mut("float_cast") else LD(p)


def planted_run(p0, slope, shift=LD(0), floor=0.0, q=QLAT, seed=0):
    """A fate function with a known threshold p0 + shift and a log-law lifetime T = c_side + slope ln(1/|p - p*|) on the row lattice q
    (c_- = -0.79, c_+ = +4.05 as #40's planted series); floor: pseudo-random offset of the threshold (a hash of p) of that rms."""
    import zlib

    def run(p):
        p = prep(p)
        thr = p0 + shift
        if floor:
            thr = thr + LD(floor) * LD(zlib.crc32((str(p) + '/' + str(seed)).encode()) / 2.0 ** 31 - 1.0)
        side = -1 if p < thr else 1
        dl = abs(p - thr)
        T = (PLANT["cm"] if side < 0 else PLANT["cp"]) + slope * float(np.log(1.0 / dl))
        T = (math.floor(T / q) + 0.5) * q
        return dict(A=p, fate=side, T=T, j=None, i=None)
    return run


def series_checks(rows, pstar, es, tag):
    """FM: every point on its own side; TM: T non-decreasing in e on each side. Returns (n_wrong, first_wrong_e, n_tm_viol)."""
    wrong = [r for r in rows if r["fate"] != r["side"]]
    first = min((r["e"] for r in wrong), default=None)
    viol = 0
    for sd in (-1, 1):
        ts = [r["T"] for r in sorted((r for r in rows if r["side"] == sd and r["T"] is not None), key=lambda r: r["e"])]
        viol += sum(1 for a, b in zip(ts[:-1], ts[1:]) if b < a - 1e-9)
    return len(wrong), first, viol


def fit_window(rows, pstar, lo, hi):
    s, se, rms, n = threshold.fit_joint(rows, pstar, (lo, hi))
    return s, se, rms, n


def ctl_k1():
    """K1: the float64 anti-diagonal sweep reproduces the numba row kernel's events."""
    cnt = bad = 0
    cases = [(2, 0.02, 28.0, 4.7527725350e-05), (3, 0.01 / SQ2, 18.0, 2.8847903691e-04), (4, 0.005, 15.0, 1.1155147990e-03)]
    for n, h, X, ps in cases:
        for p in (ps * 1.5, ps * 0.5, ps + 1e-6 * ps, ps - 1e-6 * ps, ps + 1e-9, ps - 1e-9) + ((0.01, -0.01) if n == 2 else ()):
            ref = threshold.point(n, "S", h, X, p)
            Jc = int(ref["j"] * 1.1) + 30
            pulses, sym, eps = threshold.family("S", p)
            d = dn.ellis_rays(n, Jc + 40, h, pulses=pulses, symmetric=sym, eps=eps)
            o = sweep(d, Jc, np.float64)
            cnt += 1
            same = ref["fate"] == o["fate"] and ref["j"] == o["j"] and ref["i"] == o["i"] and ref["T"] == o["T"]
            bad += 0 if same else 1
    check("K1", "float64 sweep = numba row kernel (fate, event cell, T) at D = 4, 5, 6", cnt > 0 and bad == 0, f"{cnt - bad} of {cnt} equal")


def ctl_k2():
    """K2: the static solution in binary128."""
    dl = LD(0.02) / (2 * ld_sqrt2())
    K = 2000
    R, V, P = static_taylor(2, dl, K)
    l = np.arange(K + 1) * dl
    e_r = float(np.max(np.abs(R - np.sqrt(1 + l ** 2))))
    e_p = float(np.max(np.abs(P - np.sqrt(LD(2)) * np.arctan(l))))
    fi = []
    for n in (3, 4):
        R, V, P = static_taylor(n, dl, K)
        fi.append(float(np.max(np.abs(V ** 2 - (1 - R ** (-2 * (n - 1)))))))
    check("K2", "binary128 static solution: closed form at D = 4, first integral r'^2 = 1 - r^-2(n-1) at D = 5, 6", max(e_r, e_p, *fi) <= 1e-30,
          f"D=4 r {e_r:.1e} phi {e_p:.1e}; first integral D=5 {fi[0]:.1e}, D=6 {fi[1]:.1e}; bar 1e-30")


def ctl_k3():
    """K3: the data against the float64 data."""
    A = 4.7527725350e-05
    Nc = 1200
    pulses, sym, eps = threshold.family("S", A)
    d64 = dn.ellis_rays(2, Nc, 0.02, pulses=pulses, symmetric=sym, eps=eps)
    dE = rays_ld(2, Nc, 0.02, LD(A))
    dF = rays_ld(2, Nc, 0.02, LD(A), dense=True)
    de = (float(np.max(np.abs(dE.Sp[0] - d64.Sp[0]))), float(np.max(np.abs(dE.Sp[2] - d64.Sp[2]))))
    df = (float(np.max(np.abs(dF.Sp[0] - d64.Sp[0]))), float(np.max(np.abs(dF.Sp[2] - d64.Sp[2]))))
    sym_err = float(max(np.max(np.abs(dE.Sp[0] - dE.Sm[0])), np.max(np.abs(dE.Sp[2] + dE.Sm[2]))))
    check("K3", "binary128 data against float64 data: variant E differs by the dense solution's error (rtol 1e-13), variant F by float64 roundoff",
          de[0] <= 3e-13 and de[1] <= 2e-12 and df[0] <= 1e-13 and df[1] <= 1e-14 and sym_err == 0.0,
          f"E: r {de[0]:.1e} phi {de[1]:.1e} (bars 3e-13, 2e-12); F: r {df[0]:.1e} phi {df[1]:.1e} (bars 1e-13, 1e-14); '+' and '-' rays mirror: {sym_err:.0e}")


def ctl_k4():
    """K4: the binary128 engine on float64 data reproduces the float64 events where roundoff cannot move them (e <= 8)."""
    ps = 4.7527725350e-05
    cnt = bad = 0
    for e in (3.0, 5.0, 7.0, 8.0):
        for sd in (-1, 1):
            p = ps + sd * 10.0 ** (-e)
            ref = threshold.point(2, "S", 0.02, 28.0, p)
            Jc = int(ref["j"] * 1.1) + 30
            pulses, sym, eps = threshold.family("S", p)
            d = dn.ellis_rays(2, Jc + 40, 0.02, pulses=pulses, symmetric=sym, eps=eps)
            o = sweep(d, Jc, LD, hmut=1.01 if mut("quad_h") else 1.0)
            cnt += 1
            bad += 0 if (ref["fate"] == o["fate"] and ref["j"] == o["j"] and ref["i"] == o["i"]) else 1
    check("K4", "binary128 sweep on the float64 data = the float64 events at e = 3, 5, 7, 8", cnt > 0 and bad == 0, f"{cnt - bad} of {cnt} equal")


def cmd_controls(a):
    only = set(a.only.split(",")) if a.only else None
    for name, f in (("K1", ctl_k1), ("K2", ctl_k2), ("K3", ctl_k3), ("K4", ctl_k4), ("PL", planted_checks)):
        if only is None or name in only:
            f()
    nfail = sum(1 for _, ok in RESULTS if ok is False)
    print(f"{len(RESULTS)} checks, {nfail} fail", flush=True)
    return 1 if nfail else 0


def planted_checks():
    p0 = LD(PLANT["p0"])
    slope = 1.0 / W[2]
    es = [3.0 + 0.5 * k for k in range(35)]
    # PL1: a planted fate function with a known threshold through the same bisection, series and fit
    run = planted_run(p0, slope)
    lo, hi, stopped = threshold.bisect(run, LD("4.7527e-05"), LD("4.7528e-05"), LD("1e-24"))
    ps = pstar_of(lo, hi)
    err = float(abs(ps - p0))
    rows = ld_series(run, ps, es, (-1, 1))
    nw, first, viol = series_checks(rows, ps, es, "PL1")
    s, se, rms, nn = fit_window(rows, ps, 1e-20, 1e-12)
    check("PL1", "planted threshold recovered to the bisection tolerance; planted series: every point on its side, T monotone, slope on [1e-20, 1e-12]",
          err <= 1e-24 and nw == 0 and viol == 0 and abs(s / slope - 1) <= 3e-3 and stopped is None,
          f"|p* - p0| = {err:.1e} (bar 1e-24), wrong side {nw}, T decreases {viol}, slope {s:.5f} against {slope:.5f} ({s / slope - 1:+.1e}, bar 3e-3), stopped {stopped}")
    # PL2: the perturbation analysis on a planted response p*(eta) = p0 + kappa eta
    kap = LD(PLANT["kappa"])
    ok, det = pert_analysis(lambda eta, lo_, hi_, tol: bisect_planted(p0 + kap * LD(eta), slope, lo_, hi_, tol), p0, kap_hint=float(kap), exact=None)
    check("PL2", "planted linear response p*(eta) = p0 + kappa eta: kappa from eta = 1e-8, 1e-10 predicts the shift at eta = +-1e-16", ok, det)
    # PL3: a planted floor in the fate function: the series check must find the first wrong-side point near the floor
    fl = 0.0 if mut("floor_off") else PLANT["floor_rms"]
    firsts = []
    for sd in range(PLANT["seeds"]):
        run = planted_run(p0, slope, floor=fl, seed=sd)
        rows = ld_series(run, p0, es, (-1, 1))
        nw, first, viol = series_checks(rows, p0, es, "PL3")
        firsts.append(first)
    ok = all(f is not None and 12.5 <= f <= 15.5 for f in firsts)
    check("PL3", f"a planted threshold floor of rms {fl:.0e} puts the first wrong-side point at e in [12.5, 15.5]", ok, f"first wrong e over {PLANT['seeds']} draws: {firsts}")


def bisect_planted(thr, slope, lo, hi, tol):
    """Bisection of a planted fate function whose threshold is thr; returns p* (the bracket's midpoint)."""
    run = planted_run(thr, slope)
    lo, hi, stopped = threshold.bisect(run, LD(lo), LD(hi), LD(tol))
    return pstar_of(lo, hi), stopped


PERT = dict(eta=(1e-8, 1e-10), eta_t=1e-16, lin_bar=1e-3, pred_bar=1e-4, w=(1e-6, 1e-8, 1e-14), tolrel=1e-8, tol_t=1e-24)


def pert_analysis(measure, p0, kap_hint=None, exact=None):
    """Linear response of p* to a planted data perturbation of relative size eta. measure(eta, lo, hi, tol) -> (p*, stopped). kappa from
    eta = 1e-8 and 1e-10 (agreement within lin_bar), the shift predicted at eta = +-1e-16 against the measured one (bar pred_bar relative),
    antisymmetry. exact: the exact prediction function of eta, or None. Returns (ok, detail)."""
    ks = []
    for eta, w in zip(PERT["eta"], PERT["w"]):
        ps, st = measure(eta, p0 - LD(w), p0 + LD(w), eta * PERT["tolrel"])
        ks.append(float((ps - p0) / LD(eta)))
    lin = abs(ks[0] / ks[1] - 1) if ks[1] else float("inf")
    kap = ks[1]
    meas = {}
    for sgn in (1, -1):
        eta = sgn * PERT["eta_t"]
        ps, st = measure(eta, p0 - LD(PERT["w"][2]), p0 + LD(PERT["w"][2]), PERT["tol_t"])
        meas[sgn] = float(ps - p0)
    pred = kap * PERT["eta_t"]
    if mut("pert_sign"):
        pred = -pred
    pr = [meas[1] / pred - 1, meas[-1] / (-pred) - 1] if pred else [float("inf")] * 2
    ok = lin <= PERT["lin_bar"] and all(abs(x) <= PERT["pred_bar"] for x in pr)
    det = (f"kappa {ks[0]:.6g} (eta 1e-8), {ks[1]:.6g} (1e-10), agree to {lin:.1e} (bar {PERT['lin_bar']:.0e}); predicted shift {pred:+.4e}, measured {meas[1]:+.4e} "
           f"and {-meas[-1]:+.4e} for -eta (relative error {pr[0]:+.1e}, {pr[1]:+.1e}; bar {PERT['pred_bar']:.0e})")
    return ok, det


# --------------------------------------------------------------------------------------------------------------------
# analysis of the stored runs
# --------------------------------------------------------------------------------------------------------------------
AN = dict(
    tol=1e-24, e_lo=3.0, e_hi=20.0,
    ext=(1e-20, 1e-12), dbl=(1e-11, 1e-3), full=(1e-20, 1e-3),            # fit windows in delta: the extension, #40's primary, both
    slope_bar=5e-3, agree_bar=5e-3, resid_bar=1.0,                          # slope against 1/W, extension against double window, max residual / row spacing
    dbl_floor=1e-11, evc=5, fs_band=1.5, roundoff_bar=1e-13,
)


def read_jsonl(path):
    out = []
    if os.path.exists(path):
        for line in open(path):
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def load_series(d, variant):
    rows = []
    for f in sorted(os.listdir(d)):
        if f.startswith("series_") and f.endswith(".jsonl"):
            for r in read_jsonl(os.path.join(d, f)):
                r["A"] = LD(r["A"])
                rows.append(r)
    return rows


def bracket_history(d, tol):
    """The bracket after each bisection step, from the stored points (the bisect_loose mutation stops at the first width <= tol)."""
    pts = [r for r in read_jsonl(os.path.join(d, "points_bisect.jsonl")) if r["phase"] == "bisect"]
    lo, hi = LD(pts[0]["A"]), LD(pts[1]["A"])
    for r in pts[2:]:
        A = LD(r["A"])
        if hi - lo <= LD(tol):
            break
        if r["fate"] == 1:
            hi = A
        elif r["fate"] == -1:
            lo = A
    return lo, hi


def analyze_variant(raw, variant, es_all):
    d = os.path.join(raw, "runs", variant)
    b = json.load(open(os.path.join(d, "bisect.json")))
    lo, hi = bracket_history(d, 1e-14 if mut("bisect_loose") else AN["tol"]) if mut("bisect_loose") else (LD(b["lo"]), LD(b["hi"]))
    ps = pstar_of(lo, hi) + (LD("1e-18") if mut("pstar_off") else LD(0))
    rows = load_series(d, variant)
    if mut("drop_point"):
        rows = rows[1:]
    if mut("evc_shift"):
        for r in rows:
            r["i"] = r["i"] + 10 if r["i"] is not None else None
    for r in rows:
        r["side_p"] = -1 if r["A"] < ps else 1          # the side of p* the amplitude is on, as the analysis sees it
    info = dict(b=b, ps=ps, rows=rows, d=d)
    tag = variant
    check(f"BR-{tag}", f"bisection to hi - lo <= {AN['tol']:.0e} with every run an event", b["stopped"] is None and hi - lo <= LD(AN["tol"]) * (1 + 1e-9),
          f"hi - lo = {float(hi - lo):.2e}, stopped {b['stopped']}, runs {b['runs']}, {b['secs']} s")
    if not rows and ONLY is not None:                # a view of a run-level mutation without a series
        return info
    have = sorted({r["e"] for r in rows})
    check(f"SR-{tag}", f"series e = {AN['e_lo']:g} .. {AN['e_hi']:g} step 0.5, both sides, complete", have == es_all and len(rows) == 2 * len(es_all),
          f"{len(rows)} points, e from {min(have) if have else None} to {max(have) if have else None}")
    wrong = [r for r in rows if r["fate"] != r["side_p"]]
    first = min((r["e"] for r in wrong), default=None)
    check(f"FM-{tag}", "fate monotone: every series point p* -/+ 10^-e lands on its own side, e <= 20", not wrong,
          "none wrong" if not wrong else f"{len(wrong)} wrong, first at e = {first}: " + ", ".join(f"e={r['e']:g} side {r['side_p']:+d} fate {r['fate']:+d}" for r in sorted(wrong, key=lambda r: r["e"])[:6]))
    viol = []
    for sd in (-1, 1):
        rs = sorted((r for r in rows if r["side_p"] == sd and r["T"] is not None), key=lambda r: r["e"])
        viol += [(sd, a["e"], b2["e"], a["T"], b2["T"]) for a, b2 in zip(rs[:-1], rs[1:]) if b2["T"] < a["T"] - 1e-9]
    check(f"TM-{tag}", "T non-decreasing as |p - p*| falls, each side", not viol, "none" if not viol else f"{len(viol)} decreases, first {viol[0]}")
    off = max((abs(r["i"] - r["j"]) for r in rows if r["i"] is not None), default=0)
    check(f"EVC-{tag}", f"events within {AN['evc']} cells of the diagonal (the column cap does not bind)", off <= AN["evc"] and b["off"] <= AN["evc"], f"max |i - j| = {max(off, b['off'])}")
    q = SQ2 * b["h"]
    fits = {}
    for name, w in (("ext", AN["ext"]), ("dbl", AN["dbl"]), ("full", AN["full"])):
        fits[name] = threshold.fit_joint([dict(r, side=r["side_p"]) for r in rows], ps, w)
    sx, sex, rmsx, nx = fits["ext"]
    sd_, sed, rmsd, nd = fits["dbl"]
    invW = 1.0 / W[2]
    # the fit spread of #40: every window and side choice, slope +- 2 SE
    allf = []
    for w in (AN["ext"], (1e-20, 1e-14), (1e-18, 1e-12), (1e-16, 1e-12)):
        for sides in ((-1, 1), (-1,), (1,)):
            f = threshold.fit_joint([dict(r, side=r["side_p"]) for r in rows], ps, w, sides)
            if np.isfinite(f[0]) and np.isfinite(f[1]):
                allf.append(f)
    env = (min(f[0] - 2 * f[1] for f in allf), max(f[0] + 2 * f[1] for f in allf)) if allf else (float("nan"), float("nan"))
    check(f"LL1-{tag}", f"slope of T against ln(1/delta) over [1e-20, 1e-12] (the decades the double run could not reach) against 1/W = {invW:.6f}",
          np.isfinite(sx) and abs(sx / invW - 1) <= AN["slope_bar"] and threshold.verdict_inside(env, invW),
          f"{sx:.6f} +- {sex:.6f} (OLS), {sx / invW - 1:+.2e} (bar {AN['slope_bar']:.0e}), fit spread [{env[0]:.6f}, {env[1]:.6f}], n = {nx}")
    check(f"LL2-{tag}", "the extension's slope agrees with the slope over #40's window [1e-11, 1e-3]", np.isfinite(sx) and np.isfinite(sd_) and abs(sx / sd_ - 1) <= AN["agree_bar"],
          f"extension {sx:.6f}, double window {sd_:.6f} +- {sed:.6f}: {sx / sd_ - 1:+.2e} (bar {AN['agree_bar']:.0e}); both windows {fits['full'][0]:.6f}")
    # residuals of the extension fit
    pts = [(r["side_p"], float(abs(r["A"] - ps)), r["T"]) for r in rows if r["T"] is not None and r["fate"] == r["side_p"]
           and AN["ext"][0] * (1 - 1e-6) <= float(abs(r["A"] - ps)) <= AN["ext"][1] * (1 + 1e-6)]
    if pts and np.isfinite(sx):
        c = {sd: np.mean([T - sx * math.log(1.0 / dl) for sg, dl, T in pts if sg == sd]) for sd in (-1, 1) if any(sg == sd for sg, _, _ in pts)}
        mr = max(abs(T - c[sg] - sx * math.log(1.0 / dl)) for sg, dl, T in pts)
    else:
        mr = float("nan")
    check(f"LL3-{tag}", "log law: largest residual of the extension fit <= one row spacing sqrt2 h", np.isfinite(mr) and mr <= AN["resid_bar"] * q,
          f"max residual {mr:.4f} against sqrt2 h = {q:.4f} (rms {rmsx:.4f}, uniform quantization {q / math.sqrt(12):.4f})")
    info.update(fits=fits, env=env, sx=sx)
    return info


def f64_offsets(raw, fdir, ps):
    fs = sorted(f for f in os.listdir(os.path.join(raw, "runs", fdir)) if f.startswith("f64_bisect_k") and f.endswith(".json"))
    p64 = [json.load(open(os.path.join(raw, "runs", fdir, f)))["pstar"] for f in fs]
    if mut("drop_k"):
        p64 = p64[:4]
    return np.array([float(LD(x) - ps) for x in p64]) * (1e-4 if mut("offsets_scaled") else 1.0)


def analyze_f64(raw, infos):
    """A2: the binary128 p* against K float64 bisections from different brackets (threshold.bisect, numba). A3: for variant F (the same inputs in both
    precisions) the offset is roundoff alone. FS: the float64 series about the binary128 p* leaves its side where the offsets say."""
    out = {}
    for tag, info in infos.items():
        dev = f64_offsets(raw, "f64", info["ps"])
        out[tag] = dev
        if tag == "E":
            check("KF", "float64 bisections from K >= 8 different brackets, each to hi - lo <= 1e-16", len(dev) >= 8, f"{len(dev)} bisections")
        check(f"A2-{tag}", f"p* (binary128, variant {tag}) agrees with every float64 p* to the double floor {AN['dbl_floor']:.0e}", len(dev) > 0 and np.max(np.abs(dev)) <= AN["dbl_floor"],
              f"float64 minus binary128: mean {dev.mean():+.2e}, rms {np.sqrt(np.mean(dev ** 2)):.2e}, sd {dev.std(ddof=1):.2e}, max {np.max(np.abs(dev)):.2e}, K = {len(dev)}")
    if "F" in out:
        dev = out["F"]
        check("A3-F", f"variant F (the float64 run's own inputs): the float64 offset from the binary128 p* is roundoff alone, rms <= {AN['roundoff_bar']:.0e}",
              np.sqrt(np.mean(dev ** 2)) <= AN["roundoff_bar"], f"rms {np.sqrt(np.mean(dev ** 2)):.2e}, mean {dev.mean():+.2e}, sd {dev.std(ddof=1):.2e}")
    if "E" in infos and "F" in infos:
        dd = float(infos["E"]["ps"] - infos["F"]["ps"])
        print(f"[info] data effect: p*(E) - p*(F) = {dd:+.3e} (the dense static solution's error, rtol 1e-13, moves p* by this)", flush=True)
    for tag, info in infos.items():
        rows = read_jsonl(os.path.join(info["d"], f"f64series_{tag}.jsonl"))
        if not rows:
            continue
        b128 = {(r["e"], r["side_p"]): r for r in info["rows"]}
        wrong = sorted(r["e"] for r in rows if r["fate"] != r["side"])
        diff = sorted(r["e"] for r in rows if (r["e"], r["side"]) in b128 and (r["T"] != b128[(r["e"], r["side"])]["T"] or r["fate"] != b128[(r["e"], r["side"])]["fate"]))
        mx = float(np.max(np.abs(out[tag]))) if len(out.get(tag, [])) else float("nan")
        epred = -math.log10(mx) if mx > 0 else float("nan")
        first = wrong[0] if wrong else None
        check(f"FS-{tag}", f"the float64 series about the binary128 p* leaves its side near e = -log10(max offset) = {epred:.1f} (+- {AN['fs_band']})",
              first is not None and abs(first - epred) <= AN["fs_band"],
              f"first wrong side e = {first}; {len(wrong)} wrong of {len(rows)}; first e whose event differs from binary128: {diff[0] if diff else None} ({len(diff)} differ)")
    return out


def analyze_pert(raw, info):
    """P1: relative amplitude eta of the pulse: p*(eta) = p*/(1 + eta) exactly. P2: a relative bump eta in r on both rays: kappa from eta = 1e-8, 1e-10
    predicts the shift at +-1e-16."""
    pd = os.path.join(raw, "runs", "pert")
    ps0 = info["ps"]
    got = {}
    for name in sorted(os.listdir(pd)) if os.path.isdir(pd) else []:
        f = os.path.join(pd, name, "bisect.json")
        if os.path.exists(f):
            b = json.load(open(f))
            got[name] = (pstar_of(b["lo"], b["hi"]), b)
    tolsum = LD("1.5e-24")
    for name, eta in (("scaleP16", LD("1e-16")), ("scaleM16", LD("-1e-16"))):
        if name not in got:
            check(f"P1-{name}", "scale perturbation", None, "not run")
            continue
        pred = -eta * ps0 / (1 + eta)
        if mut("pert_sign"):
            pred = -pred
        meas = got[name][0] - ps0
        check(f"P1-{name}", f"amplitude of the data's pulse scaled by 1 + eta, eta = {float(eta):+.0e}: the shift of p* is -eta p*/(1 + eta) (bar 1.5e-24, the bracket widths)",
              abs(meas - pred) <= tolsum, f"predicted {float(pred):+.4e}, measured {float(meas):+.4e}, difference {float(meas - pred):+.2e}")
    need = ("bumpP8", "bumpP10", "bumpP16", "bumpM16")
    if all(k in got for k in need):
        lk = {1e-8: got["bumpP8"][0], 1e-10: got["bumpP10"][0], 1e-16: got["bumpP16"][0], -1e-16: got["bumpM16"][0]}

        def measure(eta, lo_, hi_, tol):
            return lk[round(eta, 20) if abs(eta) < 1 else eta], None
        lk = {float(k): v for k, v in lk.items()}
        ok, det = pert_analysis(lambda eta, lo_, hi_, tol: (lk[float(eta)], None), ps0)
        check("P2", "planted bump of relative size eta in r on both rays: kappa from eta = 1e-8, 1e-10 predicts the shift of p* at eta = +-1e-16", ok, det)
    else:
        check("P2", "planted bump in r", None, "not run")
    return got


THIN = [3.0, 6.0, 9.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0]      # the thinned series of the development set
ST = dict(wt=(10.0, 17.0), rate_bar=0.02, floor_bar=3.0, nmin=4)          # window in W t: the transient e^(-W t) has decayed, the pair is still linear


def analyze_stall(raw, f64dev):
    """The float64 minus binary128 difference of r on the stall (variant F data: the same inputs in both precisions), per D: its growth rate against
    W over the declared window, the pair's own profile Q = dr/dp (growth rate against W), the equivalent offset in p = <d, Q>/<Q, Q> (the floor
    predicted from the seed) against the float64 bisections' offsets (D = 4)."""
    out = {}
    for D in (4, 5, 6):
        d = os.path.join(raw, "runs", "stall", f"D{D}")
        fs = sorted(f for f in os.listdir(d) if f.startswith("stall_F_k") and f.endswith(".json")) if os.path.isdir(d) else []
        if not fs:
            continue
        Wd = W[D - 2]
        ta, tb = ST["wt"][0] / Wd, ST["wt"][1] / Wd
        slopes, qslopes, epss, rmss, maxs, resid = [], [], [], [], [], []
        for f in fs:
            rec = json.load(open(os.path.join(d, f)))
            rows = rec["rows"]
            t = np.array([r["t"] for r in rows]) * (SQ2 if mut("t_unit") else 1.0)
            sel = (t >= ta) & (t <= tb)
            qr = np.array([r["qrms"] for r in rows])
            qslopes.append(np.polyfit(t[sel], np.log(qr[sel]), 1)[0])
            for tag in ("m", "p"):
                rm = np.array([r[tag]["rms"] for r in rows])
                mx = np.array([r[tag]["max"] for r in rows])
                ep = np.array([r[tag]["eps"] for r in rows])
                rs = np.array([r[tag]["resid"] for r in rows])
                slopes.append(np.polyfit(t[sel], np.log(rm[sel]), 1)[0])
                epss.append(float(np.mean(ep[sel])))
                rmss.append(float(np.exp(np.polyval(np.polyfit(t[sel], np.log(rm[sel]), 1), 0.0))))
                maxs.append(float(np.mean(mx[sel] / np.exp(Wd * t[sel]))))
                resid.append(float(np.mean(rs[sel])))
        slopes, qslopes, epss = np.array(slopes), np.array(qslopes), np.array(epss)
        n = len(slopes)
        out[D] = dict(eps_rms=float(np.sqrt(np.mean(epss ** 2))), eps_mean=float(epss.mean()), slope=float(slopes.mean()))
        check(f"S0-D{D}", f"the pair's profile Q grows as W = {Wd:.6f} over t in [{ta:.2f}, {tb:.2f}] (W t in {ST['wt']}; the window is in the linear regime)",
              n >= ST["nmin"] and abs(qslopes.mean() / Wd - 1) <= ST["rate_bar"], f"{qslopes.mean():.5f} ({qslopes.mean() / Wd - 1:+.2e}, bar {ST['rate_bar']:.0e}), {len(qslopes)} pairs")
        check(f"S1-D{D}", f"|r_float64 - r_binary128| (rms over |l| <= 4) on the stall grows as W over t in [{ta:.2f}, {tb:.2f}] (W t in {ST['wt']})",
              n >= ST["nmin"] and abs(slopes.mean() / Wd - 1) <= ST["rate_bar"],
              f"rate {slopes.mean():.5f} +- {slopes.std(ddof=1) / math.sqrt(n):.5f} (sd {slopes.std(ddof=1):.4f}), {slopes.mean() / Wd - 1:+.2e} against W, bar {ST['rate_bar']:.0e}; "
              f"size extrapolated to t = 0: rms {np.exp(np.mean(np.log(rmss))):.2e}; max/e^(Wt) {np.mean(maxs):.2e}; offset in p: rms {out[D]['eps_rms']:.2e}, mean {epss.mean():+.2e}; "
              f"residual after the Q projection {np.mean(resid):.2f}; {n} (pair, side) draws")
        out[D].update(n=n)
    for D, dev in sorted((f64dev or {}).items()):
        if D not in out or dev is None or not len(dev):
            continue
        dev = np.asarray(dev)
        dev_rms = float(np.sqrt(np.mean(dev ** 2)))
        sc = 10.0 if mut("seed_scale") else 1.0
        ratio = out[D]["eps_rms"] * sc / dev_rms
        check(f"S2-D{D}", f"the floor in p* predicted from the seed (offset equivalent to |r_f64 - r_b128|) is within a factor {ST['floor_bar']:g} of the float64 bisections' offsets from the binary128 p*",
              1.0 / ST["floor_bar"] <= ratio <= ST["floor_bar"], f"predicted rms {out[D]['eps_rms'] * sc:.2e} (6+ stall draws), measured rms {dev_rms:.2e} (K = {len(dev)} brackets): ratio {ratio:.2f}")
        pm, mm = -out[D]["eps_mean"] * sc, float(dev.mean())
        check(f"S3-D{D}", "the mean offset predicted from the seed has the measured mean's sign and is within a factor 3 of it (p*_float64 - p*_binary128 = -eps)",
              pm * mm > 0 and 1.0 / ST["floor_bar"] <= pm / mm <= ST["floor_bar"], f"predicted {pm:+.2e}, measured {mm:+.2e} (sd of the measured {dev.std(ddof=1):.1e})")
    return out


def analyze(a):
    global ONLY
    raw = a.raw
    ONLY = set(a.only.split(",")) if a.only else None
    es_all = THIN if a.thin else [AN["e_lo"] + 0.5 * k for k in range(int(round((AN["e_hi"] - AN["e_lo"]) / 0.5)) + 1)]
    infos = {}
    for tag in ("E", "F"):
        if os.path.exists(os.path.join(raw, "runs", tag, "bisect.json")):
            infos[tag] = analyze_variant(raw, tag, es_all)
    f64d = {}
    if infos and os.path.isdir(os.path.join(raw, "runs", "f64")):
        f64d[4] = analyze_f64(raw, infos).get("F")
    for D in (5, 6):                                   # the floors at D = 5, 6 (wave 3), when run
        bd, fd = os.path.join(raw, "runs", f"D{D}F"), os.path.join(raw, "runs", f"f64_D{D}")
        if os.path.exists(os.path.join(bd, "bisect.json")) and os.path.isdir(fd):
            bj = json.load(open(os.path.join(bd, "bisect.json")))
            f64d[D] = f64_offsets(raw, f"f64_D{D}", pstar_of(bj["lo"], bj["hi"]))
    analyze_stall(raw, f64d)
    if "E" in infos:
        analyze_pert(raw, infos["E"])
    nfail = sum(1 for _, ok in RESULTS if ok is False)
    print(f"{len(RESULTS)} checks, {sum(1 for _, ok in RESULTS if ok)} pass, {nfail} fail, {sum(1 for _, ok in RESULTS if ok is None)} skipped", flush=True)
    return 1 if nfail else 0


# --------------------------------------------------------------------------------------------------------------------
# figure (the primitives of phase2/nonradial/data_graphs.py, the arriba graph-data skill's, unchanged)
# --------------------------------------------------------------------------------------------------------------------
def cmd_graph(a):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "phase2", "nonradial"))
    import data_graphs as dg
    raw = a.raw
    es_all = THIN if a.thin else [AN["e_lo"] + 0.5 * k for k in range(35)]
    d = os.path.join(raw, "runs", "E")
    b = json.load(open(os.path.join(d, "bisect.json")))
    ps = pstar_of(b["lo"], b["hi"])
    rows = load_series(d, "E")
    f64 = read_jsonl(os.path.join(d, "f64series_E.jsonl"))
    invW = 1.0 / W[2]
    svg = dg.Svg(760)
    X0, PW, PH = 86, 600, 230
    y = 28
    y = dg.caption(svg, 16, y, "Lifetime T against the distance from the threshold, D = 4, family S, h = 0.02",
                   "binary128 (filled) and float64 at the same amplitudes (rings), p*ₕ ∓ 10^−e; dashed: the slope 1/W of the linear mode through the e = 8 points")
    xm, ym = dg.axes(svg, X0, y, PW, PH, (2.5, 20.5), (0, 48), "e, with |p − p*ₕ| = 10^−e (binary128 p*ₕ)", "lifetime T (a)", xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:g}")
    for sd, col in ((-1, dg.SERIES[0]), (1, dg.SERIES[1])):
        pts = sorted((r["e"], r["T"]) for r in rows if r["T"] is not None and r["fate"] == sd and (r["A"] < ps) == (sd < 0))
        svg.add('<polyline fill="none" stroke="' + col + '" stroke-width="1.4" points="' + " ".join(f"{xm(e):.1f},{ym(T):.1f}" for e, T in pts) + '"/>')
        for e, T in pts:
            dg.marker(svg, "circle", xm(e), ym(T), col, f"binary128, e = {e:g}, T = {T:.3f}")
        t8 = dict(pts).get(8.0)
        if t8 is not None:
            xs = np.array([3.0, 20.0])
            svg.add('<polyline fill="none" stroke="' + dg.SERIES[-1] + '" stroke-width="1" stroke-dasharray="4 3" points="' + " ".join(f"{xm(x):.1f},{ym(t8 + invW * math.log(10) * (x - 8.0)):.1f}" for x in xs) + '"/>')
    for r in f64:
        if r["T"] is not None:
            col = dg.SERIES[0] if r["side"] < 0 else dg.SERIES[1]
            svg.add(f'<circle cx="{xm(r["e"]):.1f}" cy="{ym(r["T"]):.1f}" r="7" fill="none" stroke="{col}" stroke-width="1.3"><title>float64, e = {r["e"]:g}, side {r["side"]:+d}, fate {r["fate"]:+d}, T = {r["T"]:.3f}</title></circle>')
    dg.legend(svg, X0 + PW, y - 10, [("collapse side", dg.SERIES[0]), ("expansion side", dg.SERIES[1]), ("1/W guide", dg.SERIES[-1])])
    y1 = y + PH + 84
    y = dg.caption(svg, 16, y1, "The float64 minus binary128 difference of r on the stall, as an offset in p",
                   "ε(t) = ⟨d, Q⟩/⟨Q, Q⟩ over |l| ≤ 4 for each stall draw (thin) and the mean of |ε| (thick); time in units of 1/W; a plateau is the unstable mode carrying a fixed seed")
    xm2, ym2 = dg.axes(svg, X0, y, PW, PH, (0, 20), (1e-15, 1e-11, "log"), "W t", "|ε| = offset in p equivalent to the roundoff", xfmt=lambda v: f"{v:g}", yfmt=lambda v: f"{v:.0e}")
    for D, col in ((4, dg.SERIES[0]), (5, dg.SERIES[1]), (6, dg.SERIES[2])):
        dd = os.path.join(raw, "runs", "stall", f"D{D}")
        if not os.path.isdir(dd):
            continue
        allc = []
        for f in sorted(os.listdir(dd)):
            if f.startswith("stall_F_k") and f.endswith(".json"):
                rec = json.load(open(os.path.join(dd, f)))
                for tag in ("m", "p"):
                    tt_ = np.array([r["t"] for r in rec["rows"]]) * W[D - 2]
                    ee = np.abs([r[tag]["eps"] for r in rec["rows"]])
                    sel = (tt_ >= 3) & (ee > 1e-15)
                    allc.append((tt_, ee))
                    svg.add('<polyline fill="none" stroke="' + col + '" stroke-opacity="0.35" stroke-width="0.8" points="' + " ".join(f"{xm2(t):.1f},{ym2(min(max(e, 1.01e-15), 0.99e-11)):.1f}" for t, e in zip(tt_[sel], ee[sel])) + '"/>')
        if allc:
            tg = allc[0][0]
            m = np.mean([e for _, e in allc], axis=0)
            sel = (tg >= 3) & (m > 1e-15)
            svg.add('<polyline fill="none" stroke="' + col + '" stroke-width="2" points="' + " ".join(f"{xm2(t):.1f},{ym2(min(max(e, 1.01e-15), 0.99e-11)):.1f}" for t, e in zip(tg[sel], m[sel])) + '"/>')
    svg.line(xm2(10.0), y, xm2(10.0), y + PH, dg.SERIES[-1], 1)
    svg.line(xm2(17.0), y, xm2(17.0), y + PH, dg.SERIES[-1], 1)
    dg.legend(svg, X0 + PW, y - 10, [("D = 4", dg.SERIES[0]), ("D = 5", dg.SERIES[1]), ("D = 6", dg.SERIES[2])])
    try:
        commit = subprocess.check_output(["git", "-C", os.path.dirname(os.path.abspath(__file__)), "rev-parse", "--short", "HEAD"], text=True).strip()
    except (subprocess.SubprocessError, OSError):
        commit = "?"
    height = dg.footnote(svg, [f"extended_bisect.py at {commit}", "evolution #39, threshold #40, plateau #226"], dg.host_specs(),
                         ["window of the growth-rate checks: W t in [10, 17] (grey lines)"], y1 + PH + 100, repo_url="https://github.com/nmoinvaz/negative-mass")
    open(a.out, "w").write(svg.finish(height))
    print("wrote", a.out)
    return 0


def main():
    global MUT
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd")
    ap.add_argument("--D", type=int, default=4)
    ap.add_argument("--h", type=float, default=0.02)
    ap.add_argument("--variant", default="E")
    ap.add_argument("--out")
    ap.add_argument("--raw")
    ap.add_argument("--lo", type=float, default=-0.01)
    ap.add_argument("--hi", type=float, default=0.01)
    ap.add_argument("--tol", type=float, default=1e-24)
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--eta-scale", dest="eta_scale", type=float, default=None)
    ap.add_argument("--bump", type=float, default=None)
    ap.add_argument("--e0", type=float, default=3.0)
    ap.add_argument("--e1", type=float, default=20.0)
    ap.add_argument("--elist", default=None)
    ap.add_argument("--thin", action="store_true")
    ap.add_argument("--side", type=int, default=0)
    ap.add_argument("--tmax0", type=float, default=8.0)
    ap.add_argument("--X", type=float, default=None)
    ap.add_argument("--k0", type=int, default=0)
    ap.add_argument("--k1", type=int, default=1)
    ap.add_argument("--pc", type=float, default=0.0)
    ap.add_argument("--delta", type=float, default=1e-10)
    ap.add_argument("--tmax", type=float, default=16.0)
    ap.add_argument("--snap", type=float, default=0.25)
    ap.add_argument("--L0", type=float, default=4.0)
    ap.add_argument("--only", default=None)
    ap.add_argument("--f64tol", type=float, default=1e-14)
    ap.add_argument("--mutate")
    a = ap.parse_args()
    if a.mutate:
        MUT = a.mutate
        if a.mutate in ("time_unit", "finder_sign", "ellis_charge"):
            dn.MUT = a.mutate
    cmds = dict(bisect=cmd_bisect, series=cmd_series, f64=cmd_f64, f64series=cmd_f64series, stall=cmd_stall, controls=cmd_controls, analyze=analyze, graph=cmd_graph)
    sys.exit(cmds[a.cmd](a))


if __name__ == "__main__":
    main()
