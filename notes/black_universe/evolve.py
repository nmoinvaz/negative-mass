#!/usr/bin/env python3
"""Evolution of the phantom black universes along their growing mode (#217): #39's double-null kernel with the potential V(phi) of [BF06],
the static member on the rays in the tortoise null gauge, the growing mode taken from the member's own linear-regime departure, the departure
of the gauge-invariant Y on the diagonal, bisection to A*_h, the lifetime series, and the fits. Design and rulings: evolution_design.md.

    evolve.py check [--only=E0,E1,...] [--mutate=NAME]          the code checks (E0-E7, planted series PL, PX)
    evolve.py job --c=C --h=H [--X=X] [--rho_c=R] --out=DIR [--family=mode|pulse] [--Y0=Y]    extraction, bisection of the departure sign, series (family pulse: the
                                                                one-sided pulse of #40 and the weighted trapped-cell finder, the collapse side)
    evolve.py analyze --raw=DIR [--mutate=NAME]                 E8, E9, XC, E11 and every table, from the stored jobs (directories c<c>_h<h>, O_c<c>_h<h>)
    evolve.py name --c=C --h=H --raw=DIR --out=DIR              naming runs: finders, diagonal, mass, at |A - A*| = 1e-4, 1e-5, 1e-6
    evolve.py resp --c=C [--Ap=1e-6] [--out=DIR] [--tuned=JOBDIR]    reported diagnostic D1: the linear response to a compact pulse to the end of the grid; with --tuned the
                                                                background is the mode family at that job's A* (D1b, added after D1)
    evolve.py bgdep --c=C [--h=H]                               when the untuned member's |Y| on the diagonal reaches 2e-5, 1e-3, 1e-2, 1e-1 (when D1's background leaves the member)

Units b = 1, 8 pi G = 1, D = 4, eps = -1 (the potential of [BF06], evolve_potential.py). t = (x+ + x-)/sqrt2 is the Killing time of the member (e^{-f} = A), so
T is in #212's t. Exit status 1 when a check fails (a mutated run is meant to exit 1).
"""
import json
import math
import os
import sys
import time

import numba as nb
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "wormhole_evolution"))
import dn                      # noqa: E402
import evolve_potential as ep  # noqa: E402

SQ2 = math.sqrt(2.0)
MUT = None                      # --mutate=NAME; the names are in MUTATIONS

W212 = {-0.3: 0.9700800261, -1.5: 0.4941711675, -1.01: 0.6331476235, -0.99: 0.6401216241, -1.1: 0.6031905199, -2.0: 0.3984523430,
        -1.001: 0.6362712254, -3.0: 0.2827745372, -1.0: None, 0.0: math.sqrt(1.39705243433350684)}     # notes/black_universe/README.md, Sec. 4; c = 0: #195


def mut(name):
    return MUT == name


# --------------------------------------------------------------------------------------------------------------------
# kernel: dn.py's with the potential term (new kernel, dn.py unchanged), separate steps hp (x+) and hm (x-) allowed
# --------------------------------------------------------------------------------------------------------------------
def params(mem, V=True):
    """[n, th, sc, pms, ncorr, hasV, c, K, sign of V in r_pm, coefficient of V in r_pm, sign of the V' term, sign of V in f_pm]; mutable."""
    par = np.zeros(12)
    par[:5] = dn.kernel_params(2, [-1.0])
    if mut("th_coef"):
        par[1] = 2.0
    par[5] = 1.0 if V else 0.0
    par[6], par[7] = mem.c, mem.K
    par[8], par[9], par[10], par[11] = 1.0, 0.5, 1.0, 1.0
    if mut("vterm_sign"):
        par[8] = -1.0
    if mut("vterm_coef"):
        par[9] = 1.0
    if mut("Vp_sign"):
        par[10] = -1.0
    if mut("V_sign_f"):
        par[11] = -1.0
    return par


@nb.njit(cache=True)
def _F(c, up, um, par, eps, g):
    n = par[0]
    th = par[1]
    sc = par[2]
    pms = par[3]
    nf = c.shape[0]
    emf = math.exp(-c[1])
    rpm = -th * (up[0] * um[0] + 0.5 * emf) / c[0]
    hasV = par[5] != 0.0
    Vv = 0.0
    Vp = 0.0
    if hasV:
        Vv, Vp = ep.pot(c[2], par[6], par[7])
        rpm = rpm + par[8] * par[9] * Vv * c[0] * emf
    spm = 0.0
    for k in range(nf - 2):
        spm = spm + eps[k] * up[2 + k] * um[2 + k]
    g[0] = rpm
    g1 = (n / c[0]) * rpm + pms * spm
    if hasV:
        g1 = g1 - par[11] * Vv * emf
    g[1] = g1
    kk = -sc / (2.0 * c[0])
    for k in range(2, nf):
        g[k] = kk * (up[0] * um[k] + um[0] * up[k])
    if hasV:
        g[2] = g[2] - par[10] * 0.5 * eps[0] * emf * Vp


@nb.njit(cache=True)
def row_step(old, new, dead_old, dead_new, hp, hm, par, eps, rcut, fcut, excise):
    """dn.row_step with the steps hp (x+, index i) and hm (x-, row); fills new[:, 1:] from old and new[:, :-1]; new[:, 0], dead_new[0] set by the caller."""
    nf, n1 = old.shape
    hh = hp * hm
    ncorr = int(par[4])
    c = np.empty(nf)
    up = np.empty(nf)
    um = np.empty(nf)
    nh = np.empty(nf)
    g = np.empty(nf)
    for i in range(n1 - 1):
        if dead_old[i] or dead_old[i + 1] or dead_new[i]:
            dead_new[i + 1] = 1
            for q in range(nf):
                new[q, i + 1] = old[q, i + 1]
            continue
        for q in range(nf):
            S = old[q, i]
            E = old[q, i + 1]
            W = new[q, i]
            c[q] = 0.5 * (E + W)
            up[q] = (E - S) / hp
            um[q] = (W - S) / hm
        _F(c, up, um, par, eps, g)
        for q in range(nf):
            nh[q] = old[q, i + 1] + new[q, i] - old[q, i] + hh * g[q]
        for _ in range(ncorr):
            for q in range(nf):
                S = old[q, i]
                E = old[q, i + 1]
                W = new[q, i]
                c[q] = 0.25 * (S + E + W + nh[q])
                up[q] = (E + nh[q] - S - W) / (2 * hp)
                um[q] = (W + nh[q] - S - E) / (2 * hm)
            _F(c, up, um, par, eps, g)
            for q in range(nf):
                nh[q] = old[q, i + 1] + new[q, i] - old[q, i] + hh * g[q]
        bad = False
        for q in range(nf):
            if not (nh[q] == nh[q]) or abs(nh[q]) > 1e300:
                bad = True
        if excise and (bad or nh[0] < rcut or math.exp(-nh[1]) > fcut):
            dead_new[i + 1] = 1
            for q in range(nf):
                new[q, i + 1] = old[q, i + 1]
        else:
            dead_new[i + 1] = 0
            for q in range(nf):
                new[q, i + 1] = nh[q]


def pot_np(phi, c, K):
    psi = phi / SQ2
    s2, c2 = np.sin(2.0 * psi), np.cos(2.0 * psi)
    V = -c * (2.0 - c2) - K * (1.5 * s2 + psi * (2.0 - c2))
    Vp = (-2.0 * c * s2 - 2.0 * K * (1.0 + c2 + psi * s2)) / SQ2
    return V, Vp


def F_np(c, up, um, par, eps):
    """the numpy reference of _F (arrays (nf, L))"""
    n, th, sc, pms = par[0], par[1], par[2], par[3]
    nf = c.shape[0]
    emf = np.exp(-c[1])
    rpm = -th * (up[0] * um[0] + 0.5 * emf) / c[0]
    hasV = par[5] != 0.0
    if hasV:
        Vv, Vp = pot_np(c[2], par[6], par[7])
        rpm = rpm + par[8] * par[9] * Vv * c[0] * emf
    spm = 0.0
    for k in range(nf - 2):
        spm = spm + eps[k] * up[2 + k] * um[2 + k]
    g = np.empty_like(c)
    g[0] = rpm
    g1 = (n / c[0]) * rpm + pms * spm
    if hasV:
        g1 = g1 - par[11] * Vv * emf
    g[1] = g1
    kk = -sc / (2.0 * c[0])
    for k in range(2, nf):
        g[k] = kk * (up[0] * um[k] + um[0] * up[k])
    if hasV:
        g[2] = g[2] - par[10] * 0.5 * eps[0] * emf * Vp
    return g


def run_grid(d, par, kernel="numba", rcut=0.05, fcut=25.0):
    """The whole rectangle, stored: G[q, i, j] and the dead mask [i, j]. kernel 'numba' (rows) or 'numpy' (anti-diagonals, the reference)."""
    N, h, nf, eps = d.N, d.h, d.nf, d.eps
    G = np.empty((nf, N + 1, N + 1))
    dead = np.zeros((N + 1, N + 1), np.uint8)
    G[:, :, 0] = d.Sp
    G[:, 0, :] = d.Sm
    if kernel == "numba":
        old = np.ascontiguousarray(d.Sp.copy())
        new = np.empty_like(old)
        dold = np.zeros(N + 1, np.uint8)
        dnew = np.zeros(N + 1, np.uint8)
        for j in range(N):
            new[:, 0] = d.Sm[:, j + 1]
            dnew[0] = 0
            row_step(old, new, dold, dnew, h, h, par, eps, rcut, fcut, 1)
            G[:, :, j + 1] = new
            dead[:, j + 1] = dnew
            old, new = new, old
            dold, dnew = dnew, dold
        return G, dead
    ncorr = int(par[4])
    h2 = h * h
    spert = 1.0 + 1e-8 if mut("numpy_pert") else 1.0
    for k in range(0, 2 * N - 1):
        lo, hi = max(1, k + 2 - N), min(N, k + 1)
        if hi < lo:
            continue
        I = np.arange(lo, hi + 1)
        J = k + 2 - I
        S = G[:, I - 1, J - 1] * spert
        E = G[:, I, J - 1]
        W = G[:, I - 1, J]
        dd = (dead[I - 1, J - 1] | dead[I, J - 1] | dead[I - 1, J]).astype(bool)
        with np.errstate(all="ignore"):
            c = 0.5 * (E + W)
            up, um = (E - S) / h, (W - S) / h
            out = E + W - S + h2 * F_np(c, up, um, par, eps)
            for _ in range(ncorr):
                c = 0.25 * (S + E + W + out)
                up = (E + out - S - W) / (2 * h)
                um = (W + out - S - E) / (2 * h)
                out = E + W - S + h2 * F_np(c, up, um, par, eps)
            bad = ~np.all(np.isfinite(out), axis=0)
            dd = dd | bad | (out[0] < rcut) | (np.exp(-out[1]) > fcut)
        out = np.where(dd[None, :], E, out)
        G[:, I, J] = out
        dead[I, J] = dd
    return G, dead


def rows_iter(d, par, rcut=0.05, fcut=25.0):
    """Row sweep as a generator: after the update of row j yields (j, old, new, dold, dnew); the caller may break."""
    N, h = d.N, d.h
    old = np.ascontiguousarray(d.Sp.copy())
    new = np.empty_like(old)
    dold = np.zeros(N + 1, np.uint8)
    dnew = np.zeros(N + 1, np.uint8)
    for j in range(N):
        new[:, 0] = d.Sm[:, j + 1]
        dnew[0] = 0
        row_step(old, new, dold, dnew, h, h, par, d.eps, rcut, fcut, 1)
        yield j, old, new, dold, dnew
        old, new = new, old
        dold, dnew = dnew, dold


# --------------------------------------------------------------------------------------------------------------------
# data: the member on the two rays through the corner, tortoise null gauge (e^{-f} = A, l = (x+ - x-)/sqrt2 the tortoise coordinate)
# --------------------------------------------------------------------------------------------------------------------
def tortoise_rho(mem, rho_c, xs):
    """rho, A, f at tortoise coordinates xs (any order) with rho(0) = rho_c. On the horizon side z = rho - rho_h is integrated through zeta = ln z."""
    from scipy.integrate import solve_ivp
    xs = np.asarray(xs, float)
    out_rho = np.empty_like(xs)
    out_A = np.empty_like(xs)
    out_f = np.empty_like(xs)
    zc = rho_c - mem.rh if (mem.rh is not None and rho_c > mem.rh) else None
    for sgn in (+1, -1):
        sel = np.where(sgn * xs > 0)[0]
        if len(sel) == 0:
            continue
        x = xs[sel]
        order = np.argsort(sgn * x)
        x = x[order]
        if zc is not None and sgn < 0:
            sol = solve_ivp(lambda t, y: mem.a_over_z(np.exp(y)), [0.0, x[-1]], [math.log(zc)], method="DOP853", rtol=1e-13, atol=1e-14, t_eval=x)
            zeta = sol.y[0]
            z = np.exp(zeta)
            a = mem.a_over_z(z)
            rho, A, f = mem.rh + z, z * a, -zeta - np.log(a)
        else:
            sol = solve_ivp(lambda t, y: mem.A(y), [0.0, x[-1]], [rho_c], method="DOP853", rtol=1e-13, atol=1e-14, t_eval=x)
            rho = sol.y[0]
            A = mem.A(rho)
            f = -np.log(A)
        if not sol.success or len(rho) != len(x):
            raise RuntimeError("tortoise integration failed: " + sol.message)
        idx = sel[order]
        out_rho[idx], out_A[idx], out_f[idx] = rho, A, f
    z0 = np.where(xs == 0)[0]
    if len(z0):
        out_rho[z0] = rho_c
        out_A[z0] = float(mem.A(rho_c))
        out_f[z0] = -math.log(float(mem.A(rho_c)))
    return out_rho, out_A, out_f


def member_rays(mem, rho_c, N, h, pulse=None):
    """dn.Data of the static member, corner at rho_c (x+ = x- = 0). pulse = dict(A, xc, w, ray) adds A (g - g(0)) to phi on a ray and solves the constraint
    for the change of r (RK4 on the deviation, exactly zero at A = 0), as dn.ellis_rays does."""
    s = np.arange(N + 1) * h
    rays = {}
    for ray, sg in (("+", +1.0), ("-", -1.0)):
        rho, A, f = tortoise_rho(mem, rho_c, sg * s / SQ2)
        U = np.zeros((3, N + 1))
        U[0], U[1], U[2] = mem.R(rho), f, mem.phi(rho)
        if pulse is not None and pulse["ray"] == ray and pulse["A"] != 0.0:
            gc = lambda ss: np.exp(-((ss - pulse["xc"]) / pulse["w"]) ** 2)
            g0 = math.exp(-(pulse["xc"] / pulse["w"]) ** 2)
            U[2] = U[2] + pulse["A"] * (gc(s) - g0)

            def at(ss):
                rh_, A_, _ = tortoise_rho(mem, rho_c, sg * ss / SQ2)
                ph0s = (SQ2 / (rh_ * rh_ + 1.0)) * A_ * sg / SQ2
                ps = ph0s + pulse["A"] * (-2.0 * (ss - pulse["xc"]) / pulse["w"] ** 2) * gc(ss)
                return mem.R(rh_), ph0s, ps, -mem.Ap(rh_) * sg / SQ2
            Rn, ph0n, psn, fsn = at(s)
            Rh, ph0h, psh, fsh = at(s[:-1] + 0.5 * h)
            S0, Sh = 0.5 * Rn * (psn ** 2 - ph0n ** 2), 0.5 * Rh * (psh ** 2 - ph0h ** 2)
            Q0, Qh = 0.5 * psn ** 2, 0.5 * psh ** 2
            dl = np.zeros(N + 1)
            dv = np.zeros(N + 1)

            def dd(y, S, Q, fq):
                return (y[1], S + Q * y[0] - fq * y[1])
            for i in range(N):
                y = (dl[i], dv[i])
                k1 = dd(y, S0[i], Q0[i], fsn[i])
                k2 = dd((y[0] + 0.5 * h * k1[0], y[1] + 0.5 * h * k1[1]), Sh[i], Qh[i], fsh[i])
                k3 = dd((y[0] + 0.5 * h * k2[0], y[1] + 0.5 * h * k2[1]), Sh[i], Qh[i], fsh[i])
                k4 = dd((y[0] + h * k3[0], y[1] + h * k3[1]), S0[i + 1], Q0[i + 1], fsn[i + 1])
                dl[i + 1] = y[0] + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
                dv[i + 1] = y[1] + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            U[0] = U[0] + dl
        rays[ray] = U
    return dn.Data(2, [-1.0], N, h, rays["+"], rays["-"], name="member")


class Static:
    """The exact static member at cell centres, index k = i - j + N, averaged over the four nodes like the code's centre values; Q = R'/phi0' = rho R/sqrt2."""

    def __init__(self, mem, rho_c, N, h):
        ks = np.arange(-N - 1, N + 2)
        rho, A, f = tortoise_rho(mem, rho_c, ks * h / SQ2)
        Rn, Pn = mem.R(rho), mem.phi(rho)
        sl = slice(1, -1)
        self.N = N
        self.R = 0.25 * (2 * Rn[sl] + Rn[2:] + Rn[:-2])
        self.P = 0.25 * (2 * Pn[sl] + Pn[2:] + Pn[:-2])
        self.Q = (rho * Rn / SQ2)[sl]
        self.rho = rho[sl]


def Ycell(c, st, k):
    """the gauge-invariant Y at the cell centre values c (r, f, phi) with static index k"""
    return (c[0] - st.R[k]) - st.Q[k] * (c[2] - st.P[k])


# --------------------------------------------------------------------------------------------------------------------
# streaming evolution: departure of Y on the diagonal, finders, diagonal record
# --------------------------------------------------------------------------------------------------------------------
def cell_arrays(old, new, h):
    c = 0.25 * (old[:, :-1] + old[:, 1:] + new[:, :-1] + new[:, 1:])
    up = (old[:, 1:] + new[:, 1:] - old[:, :-1] - new[:, :-1]) / (2 * h)
    um = (new[:, :-1] + new[:, 1:] - old[:, :-1] - old[:, 1:]) / (2 * h)
    return c, up, um


def evolve_dep(d, par, st, Y0=2e-5, tmin=0.0, stop=True, finders=False, record=False, wtol=0.02, rtol=1e-2, tend=None, stop_event=None):
    """Row sweep. Departure: the first row with t >= tmin whose diagonal cell has |Y| >= Y0 (t = the cell centre's (x+ + x-)/sqrt2; fate = sign of Y).
    With stop the run ends there. finders: the first trapped and anti-trapped cell anywhere (raw: dn.cell_class with tolerance rtol on r+, r-;
    weighted: the signs and 2 e^f r+ r- >= wtol), in row order. record: per row the diagonal cell's Y, r, f, r+, r-, the Misner-Sharp mass M = (r/2)(1 + 2 e^f r+ r-)
    there and at the last column's cell, and t."""
    N, h = d.N, d.h
    dep = None
    ev = dict(raw_trapped=None, raw_anti=None, w_trapped=None, w_anti=None)
    rec = None
    if record:
        rec = {k: np.full(N, np.nan) for k in ("t", "Y", "r", "f", "rp", "rm", "M", "Medge", "ledge")}
    t_of = lambda j, i: (d.x0p + d.x0m + (i + j + 1.0) * h) / SQ2
    rows = 0
    for j, old, new, dold, dnew in rows_iter(d, par):
        c, up, um = cell_arrays(old, new, h)
        alive = ~(dold[:-1].astype(bool) | dold[1:].astype(bool) | dnew[:-1].astype(bool) | dnew[1:].astype(bool))
        rows = j + 1
        tj = t_of(j, j)
        if alive[j]:
            k = N
            Y = (c[0, j] - st.R[k]) - st.Q[k] * (c[2, j] - st.P[k])
            if dep is None and tj >= tmin and abs(Y) >= Y0:
                dep = dict(j=j, T=tj, sign=1 if Y > 0 else -1, Y=float(Y))
            if record:
                rec["t"][j], rec["Y"][j], rec["r"][j], rec["f"][j] = tj, Y, c[0, j], c[1, j]
                rec["rp"][j], rec["rm"][j] = up[0, j], um[0, j]
                rec["M"][j] = 0.5 * c[0, j] * (1 + 2 * math.exp(c[1, j]) * up[0, j] * um[0, j])
        if record and alive[N - 1]:
            i = N - 1
            rec["Medge"][j] = 0.5 * c[0, i] * (1 + 2 * math.exp(c[1, i]) * up[0, i] * um[0, i])
            rec["ledge"][j] = (i - j) * h / SQ2
        if finders:
            w = 2.0 * np.exp(c[1]) * up[0] * um[0]
            for name, mask in (("raw_trapped", dn.cell_class(up[0], um[0], rtol) == -1), ("raw_anti", dn.cell_class(up[0], um[0], rtol) == 1),
                               ("w_trapped", (up[0] < 0) & (um[0] < 0) & (w >= wtol)), ("w_anti", (up[0] > 0) & (um[0] > 0) & (w >= wtol))):
                if ev[name] is None:
                    m = alive & mask
                    if m.any():
                        i = int(np.argmax(m))
                        ev[name] = dict(j=j, i=i, t=t_of(j, i), l=(i - j) * h / SQ2, r=float(c[0, i]))
        if stop and dep is not None and not finders:
            break
        if stop and dep is not None and finders and stop_event is None and all(v is not None for v in ev.values()):
            break
        if stop_event is not None and ev[stop_event] is not None:
            break
        if tend is not None and tj > tend:
            break
    return dict(dep=dep, events=ev, rows=rows, rec=rec)


# --------------------------------------------------------------------------------------------------------------------
# the mode family
# --------------------------------------------------------------------------------------------------------------------
def extract_mode(mem, rho_c, h, X0, Y_ext=2e-4, tmin=6.0):
    """Run the untuned member (the truncation error seeds the unstable mode) until the first row after tmin whose diagonal cell has |Y| >= Y_ext, and take
    the data on the two rays through that diagonal corner, minus the exact static member's, divided by Y at the corner node (so that A = Y at the corner).
    Returns (dR, dL, info), arrays (3, K + 1)."""
    N = int(round(X0 / h))
    d = member_rays(mem, rho_c, N, h)
    st = Static(mem, rho_c, N, h)
    par = params(mem)
    js = None
    col = []
    for j, old, new, dold, dnew in rows_iter(d, par):
        if js is None:
            c, up, um = cell_arrays(old, new, h)
            tj = (2 * j + 1.0) * h / SQ2
            if tj >= tmin and abs(Ycell(c[:, j], st, N)) >= Y_ext:
                js = j
                rowS = old.copy()
        if js is not None:
            col.append(new[:, js].copy())
    if js is None:
        raise RuntimeError("the untuned member did not reach the extraction level")
    K = N - js
    right = rowS[:, js:js + K + 1].copy()
    left = np.array([rowS[:, js]] + col).T[:, :K + 1].copy()
    ds = member_rays(mem, rho_c, K, h)
    dR, dL = right - ds.Sp, left - ds.Sm
    Yc = dR[0, 0] - st.Q[N] * dR[2, 0]
    info = dict(js=js, t_s=2 * js * h / SQ2, K=K, Yc=float(Yc), maxd_R=(np.abs(dR).max(axis=1) / abs(Yc)).tolist(), maxd_L=(np.abs(dL).max(axis=1) / abs(Yc)).tolist())
    return dR / Yc, dL / Yc, info


def mode_rays(mem, rho_c, N, h, A, dR, dL):
    ds = member_rays(mem, rho_c, N, h)
    return dn.Data(2, [-1.0], N, h, ds.Sp + A * dR[:, :N + 1], ds.Sm + A * dL[:, :N + 1], name="mode")


# --------------------------------------------------------------------------------------------------------------------
# bisection, series, fits (the same code for the evolution and for the planted oracle)
# --------------------------------------------------------------------------------------------------------------------
CRIT = dict(
    Y0=2e-5, Y_ext=2e-4, bracket=1e-4, tmin_extract=6.0,                      # departure and extraction levels; ruling 1 of the diagnosis
    e_lo=5.0, e_hi=11.5, e_step=0.25,                                          # the series: A = A* -/+ 10^-e; below 1e-7 (e > 7) run and not fitted
    windows=[(1e-7, 1e-5), (3e-7, 1e-5), (1e-7, 3e-6), (3e-7, 3e-6), (1e-7, 1e-6)],   # the first is primary
    bisect_tol=1e-14,
    bracket_pulse=2e-4, pulse=dict(xc=3.0, w=0.5), pulse_ray={-0.3: "+", -1.5: "-"},   # family O (one-sided pulse, #40's): the ray is the one whose collapse side has finder events
    windows_pulse=[(1e-8, 1e-5), (3e-8, 1e-5), (1e-8, 3e-6), (3e-8, 3e-6), (1e-8, 1e-6)],
    H={-0.3: (0.02, 0.01, 0.005, 0.0025), -1.5: (0.02, 0.01, 0.005, 0.0025), -1.0: (0.02, 0.01, 0.005, 0.0025)},   # coarse member first
    X={-0.3: 40.0, -1.5: 50.0, -1.0: 50.0}, RHO_C={-0.3: 0.0, -1.5: 1.0, -1.0: 1.0},
    order_lo=1.75, order_hi=2.25,
)


def bisect(run, lo, hi, tol, sink=None):
    """run(A) -> dict(A, fate, T). Fate -1 or +1 on the two sides; bisection to hi - lo <= tol. Returns lo, hi, fate(lo), count."""
    r_lo, r_hi = run(lo), run(hi)
    for r in (r_lo, r_hi):
        if sink:
            sink(dict(phase="bisect", **r))
    if r_lo["fate"] * r_hi["fate"] != -1:
        raise RuntimeError(f"bracket [{lo}, {hi}] does not straddle the threshold: fates {r_lo['fate']}, {r_hi['fate']} (T {r_lo['T']}, {r_hi['T']})")
    fl, k = r_lo["fate"], 2
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        r = run(mid)
        k += 1
        if sink:
            sink(dict(phase="bisect", **r))
        if r["fate"] == fl:
            lo = mid
        elif r["fate"] == -fl:
            hi = mid
        else:
            raise RuntimeError(f"no departure at A = {mid!r}")
    return lo, hi, fl, k


def pstar_of(lo, hi):
    return 0.5 * (lo + hi) + (3e-8 if mut("pstar_off") else 0.0)


def series(run, ps, es, sink=None, phase="series"):
    out = []
    for e in es:
        for s in (-1, 1):
            r = run(ps + s * 10.0 ** (-e))
            r.update(phase=phase, side=s, e=float(e))
            if sink:
                sink(dict(r))
            out.append(r)
    return out


def fit_joint(rows, ps, window, orient, sides=(-1, 1)):
    """T = c_side + s ln(1/|A - A*|), common slope, one constant per side, over the points whose fate is the side's own. Returns (s, SE, rms, n)."""
    sc = SQ2 if mut("time_unit") else 1.0
    pts = [(r["side"], abs(r["A"] - ps), r["T"] * sc) for r in rows
           if r.get("T") is not None and r["side"] in sides and r["fate"] == -r["side"] * orient and window[0] * (1 - 1e-6) <= abs(r["A"] - ps) <= window[1] * (1 + 1e-6)]
    present = sorted({s for s, _, _ in pts})
    if len(pts) < 4 or not present:
        return float("nan"), float("nan"), float("nan"), len(pts)
    M = np.array([[1.0 if s == sd else 0.0 for sd in present] + [math.log(1.0 / d)] for s, d, _ in pts])
    y = np.array([T for _, _, T in pts])
    coef, *_ = np.linalg.lstsq(M, y, rcond=None)
    res = y - M @ coef
    dof = max(1, len(y) - M.shape[1])
    cov = (float(res @ res) / dof) * np.linalg.inv(M.T @ M)
    return float(coef[-1]), float(math.sqrt(cov[-1, -1])), float(math.sqrt(float(res @ res) / len(y))), len(y)


def fit_set(rows, ps, orient, windows=None):
    """each declared window, joint and each side alone; the primary is the first window, joint; the fit spread is the envelope of s -/+ 2 SE over the fits."""
    allf = []
    for w in (windows or CRIT["windows"]):
        for sides in ((-1, 1), (-1,), (1,)):
            s, se, rms, n = fit_joint(rows, ps, w, orient, sides)
            allf.append(dict(window=w, sides=sides, slope=s, se=se, rms=rms, n=n))
    prim = allf[0]
    ok = [f for f in allf if np.isfinite(f["slope"]) and np.isfinite(f["se"])]
    env = (min(f["slope"] - 2 * f["se"] for f in ok), max(f["slope"] + 2 * f["se"] for f in ok)) if ok else (float("nan"), float("nan"))
    return allf, prim, env


def es_series():
    return list(np.arange(CRIT["e_lo"], CRIT["e_hi"] + 1e-9, CRIT["e_step"]))


# --------------------------------------------------------------------------------------------------------------------
# job: one (member, h): extraction, bisection of the departure sign, series
# --------------------------------------------------------------------------------------------------------------------
class Sink:
    def __init__(self, path):
        self.f = open(path, "a")

    def __call__(self, rec):
        self.f.write(json.dumps(rec) + "\n")
        self.f.flush()


def family_run(mem, rc, h, X, dR, dL, Y0, tmin=0.0):
    N = int(round(X / h))
    par, st = params(mem), Static(mem, rc, N, h)

    def run(A):
        d = mode_rays(mem, rc, N, h, A, dR, dL)
        r = evolve_dep(d, par, st, Y0=Y0, tmin=tmin)
        dp = r["dep"]
        return dict(A=A, fate=0 if dp is None else dp["sign"], T=None if dp is None else dp["T"])
    return run


def pulse_run(mem, rc, h, X, ray):
    """family O: a one-sided Gaussian pulse of amplitude A in phi on the given ray. fate = -1 if a trapped cell (weighted finder) appears before the end of
    the grid, T = its time (the collapse side's lifetime, as #40's); fate = +1 otherwise (T None)."""
    N = int(round(X / h))
    par, st = params(mem), Static(mem, rc, N, h)

    def run(A):
        d = member_rays(mem, rc, N, h, pulse=dict(A=A, ray=ray, **CRIT["pulse"]))
        r = evolve_dep(d, par, st, Y0=1e9, stop=True, finders=True, stop_event="w_trapped")
        e = r["events"]["w_trapped"]
        return dict(A=A, fate=1 if e is None else -1, T=None if e is None else e["t"])
    return run


def job(c, h, out, X=None, rc=None, Y0=None, family="mode"):
    X = X or CRIT["X"][c]
    rc = CRIT["RHO_C"][c] if rc is None else rc
    Y0 = Y0 or CRIT["Y0"]
    os.makedirs(out, exist_ok=True)
    sink = Sink(os.path.join(out, "points.jsonl"))
    mem = ep.Member(c)
    t0 = time.time()
    summ = dict(c=c, h=h, X=X, rho_c=rc, N=int(round(X / h)), Y0=Y0, Y_ext=CRIT["Y_ext"], mutation=MUT, family=family)
    try:
        if family == "mode":
            dR, dL, info = extract_mode(mem, rc, h, X + 25.0, CRIT["Y_ext"], CRIT["tmin_extract"])
            summ["extraction"] = info
            run = family_run(mem, rc, h, X, dR, dL, Y0)
            br = CRIT["bracket"]
        else:
            run = pulse_run(mem, rc, h, X, CRIT["pulse_ray"][c])
            br = CRIT["bracket_pulse"]
            summ["extraction"] = dict(t_s=float("nan"))
        lo, hi, fl, k = bisect(run, -br, br, CRIT["bisect_tol"], sink)
        summ.update(lo=lo, hi=hi, fate_lo=fl, bisect_runs=k, bisect_secs=round(time.time() - t0, 1))
        json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
        ps = 0.5 * (lo + hi)
        series(run, ps, es_series(), sink)
        summ["status"] = "done"
    except Exception as e:           # recorded; analyze reports the job as failed
        import traceback
        traceback.print_exc()
        summ["status"] = "failed: " + repr(e)
    summ["secs"] = round(time.time() - t0, 1)
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(summ), flush=True)
    return summ["status"] == "done"


# --------------------------------------------------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------------------------------------------------
results = []
ONLY = None


def check(cid, text, ok, detail=""):
    ok = None if ok is None else bool(ok)
    results.append((cid, ok))
    tag = "SKIP" if ok is None else ("pass" if ok else "FAIL")
    print(f"[{tag}] {cid}: {text}" + (f"  ({detail})" if detail else ""), flush=True)


def want(cid):
    return ONLY is None or cid in ONLY


def ellis_grids():
    """three grids of E0: static, a collapsing one-sided ghost pulse (cells excised), an expanding one"""
    out = []
    for name, N, h, A in (("static", 300, 0.02, 0.0), ("collapse", 400, 0.04, -0.05), ("expand", 400, 0.04, 0.05)):
        pulses = [dict(field=0, ray="+", A=A, xc=3.0, w=0.5)] if A != 0.0 else []
        out.append((name, dn.ellis_rays(2, N, h, pulses=pulses)))
    return out


def check_E0():
    worst, notes = 0.0, []
    same = True
    mem = ep.Member(0.0)
    for name, d in ellis_grids():
        for kernel in ("numba", "numpy"):
            ref = dn.evolve(d, kernel=kernel, store=True, stop=False)
            for V in (False, True):
                G, dead = run_grid(d, params(mem, V=V), kernel=kernel)
                diff = float(np.max(np.abs(G - ref.G)))
                eqd = bool((dead == ref.dead).all())
                worst = max(worst, diff)
                same &= eqd
                notes.append(f"{name}/{kernel}/V{'on' if V else 'off'} {diff:.0e}{'' if eqd else ' DEAD'}")
    check("E0", "new kernel with the potential flag off and on at c = 0 (V = 0) equals dn.evolve: three grids (static, collapse with excision, expansion), numba and numpy",
          worst == 0.0 and same, f"max |dG| {worst:.1e}; dead masks equal {same}")


def check_E1():
    import sympy as sp
    rng = np.random.default_rng(217)
    r, f, ph, rp, rm, pp, pm, c, K = sp.symbols("r f phi rp rm pp pm c K")
    psi = ph / sp.sqrt(2)
    V = -c * (3 - 2 * sp.cos(psi) ** 2) - K * (3 * sp.sin(psi) * sp.cos(psi) + psi * (3 - 2 * sp.cos(psi) ** 2))     # [BF06] (14), written out
    Vp = sp.diff(V, ph)
    eps = -1
    emf = sp.exp(-f)
    rpm = -(rp * rm + emf / 2) / r + V * r * emf / 2
    fpm = (2 / r) * rpm + eps * pp * pm - V * emf
    ppm = -(rp * pm + rm * pp) / r - sp.Rational(eps, 2) * emf * Vp
    fn = sp.lambdify((r, f, ph, rp, rm, pp, pm, c, K), [rpm, fpm, ppm], "numpy")
    worst = 0.0
    for cc in (-0.3, -1.5, -1.0, 0.0):
        mem = ep.Member(cc)
        par = params(mem)
        for _ in range(40):
            st = dict(r=rng.uniform(0.3, 3), f=rng.uniform(-1, 3), ph=rng.uniform(-2.2, 2.2), rp=rng.normal(), rm=rng.normal(), pp=rng.normal(), pm=rng.normal())
            cv = np.array([st["r"], st["f"], st["ph"]])
            up = np.array([st["rp"], 0.0, st["pp"]])
            um = np.array([st["rm"], 0.0, st["pm"]])
            g = np.empty(3)
            _F(cv, up, um, par, np.array([-1.0]), g)
            mine = np.array(fn(st["r"], st["f"], st["ph"], st["rp"], st["rm"], st["pp"], st["pm"], mem.c, mem.K), float)
            worst = max(worst, float(np.max(np.abs(g - mine) / (1.0 + np.abs(mine)))))
    check("E1", "the kernel's right-hand side equals the equations with V (V and V' of [BF06] written out and differentiated in sympy) at 40 random states, c = -0.3, -1.5, -1, 0",
          worst <= 1e-13, f"max relative difference {worst:.1e}, tolerance 1e-13")


def check_E2():
    worst = 0.0
    same = True
    for cc, rc in ((-0.3, 0.0), (-1.5, 1.0), (-1.0, 1.0)):
        mem = ep.Member(cc)
        N, h = 250, 0.04
        d = member_rays(mem, rc, N, h, pulse=dict(A=0.02, xc=3.0, w=0.5, ray="+"))
        par = params(mem)
        Ga, da = run_grid(d, par, "numba")
        Gb, db = run_grid(d, par, "numpy")
        worst = max(worst, float(np.max(np.abs(Ga - Gb))))
        same &= bool((da == db).all())
    check("E2", "numba and numpy kernels agree with V on, three members (c = -0.3, -1.5, -1; one-sided pulse A = 0.02, h = 0.04, N = 250), same dead masks",
          worst <= 1e-12 and same, f"max |dG| {worst:.1e}, tolerance 1e-12, dead masks equal {same}")


def check_E3():
    N, h = 1000, 0.02
    d1 = dn.ellis_rays(2, N, h, pulses=[dict(field=0, ray="+", A=0.01, xc=3.0, w=0.5)])
    d2 = member_rays(ep.Member(0.0), 0.0, N, h, pulse=dict(A=0.01, xc=3.0, w=0.5, ray="+"))
    er = max(np.abs(d1.Sp[0] - d2.Sp[0]).max(), np.abs(d1.Sm[0] - d2.Sm[0]).max())
    ep_ = max(np.abs(d1.Sp[2] - d2.Sp[2]).max(), np.abs(d1.Sm[2] - d2.Sm[2]).max())
    ef = max(np.abs(d1.Sp[1] - d2.Sp[1]).max(), np.abs(d1.Sm[1] - d2.Sm[1]).max())
    check("E3", "the c = 0 member with a one-sided pulse, on the rays, equals dn.ellis_rays (N = 1000, h = 0.02)", er <= 2e-9 and ep_ <= 1e-12 and ef == 0.0,
          f"r {er:.1e} (tolerance 2e-9, #39's I1), phi {ep_:.1e} (1e-12), f {ef:.0e} (0)")


def exact_nodes(mem, rc, N, h):
    ks = np.arange(-N, N + 1)
    rho, A, f = tortoise_rho(mem, rc, ks * h / SQ2)
    i = np.arange(N + 1)[:, None]
    j = np.arange(N + 1)[None, :]
    kk = (i - j) + N
    return mem.R(rho)[kk], f[kk], mem.phi(rho)[kk], (i + j) * h / SQ2 + 0 * kk


def check_E4():
    allok, notes = True, []
    for cc, rc in ((0.0, 0.0), (-0.3, 0.0), (-1.5, 1.0), (-1.0, 1.0)):
        mem = ep.Member(cc)
        errs = []
        for h in (0.04, 0.02, 0.01):
            N = int(round(10.0 / h))
            d = member_rays(mem, rc, N, h)
            G, dead = run_grid(d, params(mem))
            R, F, P, t = exact_nodes(mem, rc, N, h)
            e = np.maximum(np.maximum(np.abs(G[0] - R), np.abs(G[1] - F)), np.abs(G[2] - P))
            errs.append(float(e[t <= 8.0].max()))
        o = [math.log2(errs[0] / errs[1]), math.log2(errs[1] / errs[2])]
        ok = all(CRIT["order_lo"] <= x <= CRIT["order_hi"] for x in o)
        allok &= ok
        notes.append(f"c = {cc}: errors {errs[0]:.2e}, {errs[1]:.2e}, {errs[2]:.2e}, orders {o[0]:.2f}, {o[1]:.2f}")
    check("E4", "static members c = 0, -0.3, -1.5, -1 stay static to the truncation error for t <= 8 (X = 10, h = 0.04, 0.02, 0.01): second order", allok,
          "orders within [1.75, 2.25]; " + "; ".join(notes))


E5_WINDOW = (1e-4, 1e-3)             # |Y| levels of the fit: the linear regime after the start-up transient (dev: [1e-5, 1e-4] is still inside it at c = -0.3)
E5_TOL = 0.015


def check_E5():
    allok, notes = True, []
    for cc, rc in ((-0.3, 0.0), (-1.5, 1.0)):
        mem = ep.Member(cc)
        devs = []
        for h in (0.01, 0.005):
            N = int(round(25.0 / h))
            d = member_rays(mem, rc, N, h)
            r = evolve_dep(d, params(mem), Static(mem, rc, N, h), Y0=1e9, stop=False, record=True, tend=24.0)
            t, Y = r["rec"]["t"], np.abs(r["rec"]["Y"])
            m = np.isfinite(Y) & (Y >= E5_WINDOW[0]) & (Y <= E5_WINDOW[1]) & (t >= 4.0)
            s = np.polyfit(t[m], np.log(Y[m]), 1)[0]
            dev = s / W212[cc] - 1
            devs.append(dev)
            allok &= abs(dev) <= E5_TOL
            notes.append(f"c = {cc}, h = {h}: {s:.4f} against W = {W212[cc]:.5f} ({dev:+.1e}, {int(m.sum())} rows)")
        allok &= abs(devs[1]) < abs(devs[0])
    check("E5", f"the growth rate of Y on the diagonal, from the truncation seed of the untuned member over |Y| in [{E5_WINDOW[0]:g}, {E5_WINDOW[1]:g}], equals #212's W at c = -0.3 and -1.5, "
          "and the error falls from h = 0.01 to 0.005", allok, f"tolerance {E5_TOL:.1%}; " + "; ".join(notes))


def check_E6():
    import evolve_master1d as m1
    allok, notes = True, []
    for cc in (-2.0, -1.5, -1.1, -1.01, -1.001):
        ts, ys, mem = m1.run(cc)
        m = (ts >= 45) & (ts <= 60) & (ys > 0)
        s = np.polyfit(ts[m], np.log(ys[m]), 1)[0]
        ok = abs(s / W212[cc] - 1) <= 1e-4
        allok &= ok
        notes.append(f"{cc}: {s:.5f} ({s / W212[cc] - 1:+.0e})")
    ts, ys, mem = m1.run(-1.0)
    m = (ts >= 45) & (ts <= 60) & (ys > 0)
    s = np.polyfit(ts[m], np.log(ys[m]), 1)[0]
    allok &= s < 0
    notes.append(f"-1: {s:+.4f} (must be < 0)")
    check("E6", "the 1D master equation in the time domain (no double-null code; shares only #212's V-hat) reproduces #212's W, and decays at c = -1", allok,
          "rate of ln|u(x = 0)| over t = 45 to 60, tolerance 1e-4: " + "; ".join(notes))


def check_planted():
    import threshold as th
    seeds = range(1, 41)
    es = es_series()
    allpl = allpx = True
    notes = []
    for cc in (-0.3, -1.5):
        W = W212[cc]
        h = 0.005
        q = SQ2 * h
        cov = cp = rej = 0
        errs = []
        for sd in seeds:
            ps0 = -7.0e-7
            orc = th.oracle(1.0 / W, {-1: 5.0, 1: 6.5}, ps0, 0.01, q, sd)
            lo, hi, fl, k = bisect(orc, -CRIT["bracket"], CRIT["bracket"], CRIT["bisect_tol"])
            ps = pstar_of(lo, hi)
            rows = series(orc, 0.5 * (lo + hi), es)
            allf, prim, env = fit_set(rows, ps, fl)
            cov += bool(env[0] <= 1.0 / W <= env[1])
            cp += abs(ps - ps0) <= 2 * CRIT["bisect_tol"]
            errs.append(prim["slope"] * W - 1)
            orx = th.oracle(1.0 / (SQ2 * W), {-1: 5.0, 1: 6.5}, ps0, 0.01, q, sd)
            lo, hi, fl, k = bisect(orx, -CRIT["bracket"], CRIT["bracket"], CRIT["bisect_tol"])
            rows = series(orx, 0.5 * (lo + hi), es)
            allf, prim, env = fit_set(rows, pstar_of(lo, hi), fl)
            rej += not (env[0] <= 1.0 / W <= env[1])
        m = len(list(seeds))
        allpl &= cov >= 0.9 * m and cp == m
        allpx &= rej == m
        notes.append(f"W = {W:.5f}: 1/W inside the spread in {cov} of {m}, A* recovered in {cp} of {m}, primary error mean {np.mean(errs):+.1e}; 1/(sqrt2 W) rejected in {rej} of {m}")
    check("PL", "planted series (slope 1/W, sigma 0.01 plus the row quantization) through the same bisect, series, fit_set: A* recovered, 1/W inside the fit spread in >= 90 % of 40 seeds",
          allpl, "; ".join(notes))
    check("PX", "a planted series with the time unit scaled by sqrt2 (slope 1/(sqrt2 W)) is rejected by the comparison with 1/W in every seed", allpx, "")


def y212_profile(c, W, rho_grid):
    """#212's Y(rho) of the growing mode: the ingoing solution from the horizon and the decaying one from large rho, matched at stability.match_point,
    through its own right-hand sides; for a throat member (rh < 0) the throat is crossed on the upper half circle as stability.left_solution does.
    Returns Y on rho_grid (nan where not covered), normalized to 1 at rho = 1."""
    import stability as sb
    from scipy.integrate import solve_ivp
    from scipy.interpolate import CubicSpline
    fam = sb.Fam(c, "R2")
    p = fam.res
    rm = sb.match_point(fam)
    pieces = []        # (rho array increasing, Y array)

    def dense(f, s0, x1, event, ngrid=6000):
        sol = solve_ivp(f, [0.0, x1], s0, method="DOP853", rtol=1e-12, atol=[1e-12, 1e-300, 1e-300], dense_output=True, events=event)
        return sol
    # left, real part: from the horizon outward
    if fam.rh is None:
        raise ValueError("c = 0 not handled")
    s0 = np.array([math.log(p["z0"]), 1.0, (fam.A1 + W) if fam.rh == 0.0 else W])
    rstop = rm if fam.rh >= 0 else -min(p["delta"], 0.5 * abs(fam.rh))
    ev = lambda x, s: fam.rh + math.exp(s[0]) - rstop
    ev.terminal = True
    ev.direction = 1
    sol = dense(sb.rhs_zeta(fam, W), s0, 1e4, ev)
    xs = np.linspace(0.0, sol.t[-1], 8000)
    ys = sol.sol(xs)
    rho_l = fam.rh + np.exp(ys[0])
    Yl = ys[1]
    pieces.append((rho_l, Yl))
    rho_last, Y_last, Yx_last = rho_l[-1], Yl[-1], ys[2][-1]
    if fam.rh < 0:
        d = -rho_last
        arc = solve_ivp(lambda th_, y: sb.rhs_arc(fam, W, d)(th_, y), [np.pi, 0.0], [Y_last + 0j, Yx_last / fam.A(rho_last) + 0j], method="DOP853", rtol=1e-12, atol=1e-30)
        Yc, Yrc = arc.y[:, -1]
        s1 = np.array([d, Yc.real, (Yrc * fam.A(d)).real])
        evr = lambda x, s: s[0] - rm
        evr.terminal = True
        evr.direction = 1
        sol2 = dense(sb.rhs_x(fam, W), s1, 1e4, evr)
        xs2 = np.linspace(0.0, sol2.t[-1], 6000)
        ys2 = sol2.sol(xs2)
        pieces.append((ys2[0], ys2[1]))
    # right: decaying solution from large rho inward, in x decreasing
    r0 = max(p["rmax0"], p["rmaxW"] / W)
    lg = -W - fam.Vhat(r0) / (2 * W)
    s2 = np.array([r0, 1.0, fam.A(r0) / r0 + lg])
    evl = lambda x, s: s[0] - rm
    evl.terminal = True
    evl.direction = -1
    sol3 = solve_ivp(sb.rhs_x(fam, W), [0.0, -1e4], s2, method="DOP853", rtol=1e-12, atol=[1e-12, 1e-300, 1e-300], dense_output=True, events=evl)
    xs3 = np.linspace(0.0, sol3.t[-1], 12000)
    ys3 = sol3.sol(xs3)
    rho_r, Yr = ys3[0][::-1], ys3[1][::-1]
    # amplitude match at rm
    Yleft_rm = pieces[-1][1][-1]
    Yr_rm = Yr[0]
    sc = Yleft_rm / Yr_rm
    rho_all = np.concatenate([pieces[0][0]] + ([pieces[1][0]] if len(pieces) > 1 else []) + [rho_r[1:]])
    Y_all = np.concatenate([pieces[0][1]] + ([pieces[1][1]] if len(pieces) > 1 else []) + [sc * Yr[1:]])
    o = np.argsort(rho_all)
    rho_all, Y_all = rho_all[o], Y_all[o]
    keep = np.concatenate([[True], np.diff(rho_all) > 0])
    spl = CubicSpline(rho_all[keep], Y_all[keep])
    out = np.full_like(np.asarray(rho_grid, float), np.nan)
    ok = (rho_grid >= rho_all[0]) & (rho_grid <= rho_all[-1])
    if fam.rh < 0:                    # the stretch |rho| < d inside the arc is not covered
        ok &= np.abs(rho_grid) >= d
    out[ok] = spl(np.asarray(rho_grid)[ok])
    return out


E7_BAR = 0.01            # declared: max |profile difference| over l in [-6, 6] (outside the throat's arc), relative to the profile's maximum
E7_LEVELS = (1e-3,)      # the time slice: the first row after t = 4 whose diagonal cell has |Y| >= the level (the linear regime, after the start-up transient)


def check_E7():
    allok, notes = True, []
    for cc, rc in ((-1.5, 1.0), (-1.01, 1.0), (-1.1, 1.0), (-0.3, 0.0)):
        mem = ep.Member(cc)
        h, X = 0.005, 24.0
        N = int(round(X / h))
        d = member_rays(mem, rc, N, h)
        st = Static(mem, rc, N, h)
        G, dead = run_grid(d, params(mem))
        ls = np.linspace(-6.0, 6.0, 49)
        rho_l, A_l, f_l = tortoise_rho(mem, rc, ls)
        Yref = y212_profile(cc, W212[cc], rho_l)
        jj = np.arange(N)
        cd = 0.25 * (G[:, jj, jj] + G[:, jj + 1, jj] + G[:, jj, jj + 1] + G[:, jj + 1, jj + 1])
        Yd = Ycell(cd, st, N)
        worst = 0.0
        for lev in E7_LEVELS:
            jstar = int(np.argmax((np.abs(Yd) >= lev) & ((2 * jj + 1) * h / SQ2 >= 4.0)))
            s_ = 2 * jstar + 1
            prof = []
            for l0 in ls:
                k = int(round(l0 * SQ2 / h))
                i = (s_ + k) // 2
                j = s_ - i
                if not (1 <= j < N and 1 <= i < N):
                    prof.append(np.nan)
                    continue
                c4 = 0.25 * (G[:, i, j] + G[:, i + 1, j] + G[:, i, j + 1] + G[:, i + 1, j + 1])
                prof.append(Ycell(c4, st, (i - j) + N))
            prof = np.array(prof)
            i0 = int(np.argmin(np.abs(ls - 1.0)))
            prof = prof / prof[i0]
            yr = Yref / Yref[i0]
            good = np.isfinite(prof) & np.isfinite(yr)
            worst = max(worst, float(np.max(np.abs(prof[good] - yr[good])) / np.max(np.abs(yr[good]))))
        allok &= worst <= E7_BAR
        notes.append(f"c = {cc}: {worst:.1e}")
    check("E7", f"the shape of the code's growing Y(l) (untuned member, h = 0.005, the slice where the diagonal |Y| first reaches {E7_LEVELS[0]:g}, normalized at l = 1) equals #212's Y(rho(l)) "
          "(ingoing from the horizon, matched to the decaying solution; through the throat's arc at c = -0.3)", allok,
          f"max difference over l in [-6, 6] relative to the maximum, bar {E7_BAR:g}: " + "; ".join(notes))


QUICK = dict(E0=check_E0, E1=check_E1, E2=check_E2, E3=check_E3, E4=check_E4, E5=check_E5, E6=check_E6, E7=check_E7, PL=check_planted)

MUTATIONS = {   # name: the checks that must fire in `check`; computed in the development runs dev_111 to dev_119 of negative-mass-raw I217-mac before they were declared
    "th_coef": ["E0", "E1", "E4", "E5", "E7"],            # dn.py's sphere coefficient n - 1 replaced by n (the kernel's, not dn's)
    "vterm_sign": ["E1", "E4", "E5", "E7"],               # sign of V in r_pm
    "vterm_coef": ["E1", "E4", "E5", "E7"],               # coefficient 1/2 of V r e^{-f} in r_pm replaced by 1
    "Vp_sign": ["E1", "E4", "E5", "E7"],                  # sign of the V' term in phi_pm
    "V_sign_f": ["E1", "E4", "E5", "E7"],                 # sign of the V e^{-f} term in f_pm
    "numpy_pert": ["E0", "E2"],                           # a transcription slip in the numpy kernel only
    "time_unit": ["PL", "PX"],                            # T scaled by sqrt2 in the fits (also fires E9 in `analyze`)
    "pstar_off": ["PL"],                                  # A* shifted by 3e-8 in the fits (also fires E9 in `analyze`)
}


# --------------------------------------------------------------------------------------------------------------------
# analysis of the stored jobs: E8, E9, E11, tables, floors
# --------------------------------------------------------------------------------------------------------------------
TARGET = lambda c: (1.0 / W212[c], "1/W") if W212.get(c) else (1.0 / ep.Member(c).kappa, "1/kappa")


def job_dir(raw, c, h):
    return os.path.join(raw, f"c{c:g}_h{h:g}")


def load_jobs(raw, family="mode"):
    jobs = {}
    for name in sorted(os.listdir(raw)):
        p = os.path.join(raw, name)
        if os.path.isfile(os.path.join(p, "summary.json")) and name.startswith("c" if family == "mode" else "O_c"):
            summ = json.load(open(os.path.join(p, "summary.json")))
            pts = []
            if os.path.exists(os.path.join(p, "points.jsonl")):
                pts = [json.loads(x) for x in open(os.path.join(p, "points.jsonl"))]
            jobs[(summ["c"], summ["h"])] = dict(summ, pts=pts)
    return jobs


def floor_of(rows, ps, orient, prim_slope, window, side):
    """the largest |A - A*| below the window at which the series leaves the log law of the window (offset of the side's points) by more than 0.25 in T; None if it does not"""
    sc = SQ2 if mut("time_unit") else 1.0
    pts = [(abs(r["A"] - ps), r["T"] * sc) for r in rows if r["side"] == side and r["T"] is not None and r["fate"] == -side * orient]
    inw = [(d, T) for d, T in pts if window[0] * (1 - 1e-6) <= d <= window[1] * (1 + 1e-6)]
    if len(inw) < 3 or not np.isfinite(prim_slope):
        return None
    off = float(np.mean([T - prim_slope * math.log(1 / d) for d, T in inw]))
    below = sorted([(d, T) for d, T in pts if d < window[0] * (1 - 1e-6)], reverse=True)
    for d, T in below:
        if abs(T - (off + prim_slope * math.log(1 / d))) > 0.25:
            return d
    return None


def fmtf(x, f="{:.1e}"):
    return "none" if x is None else f.format(x)


def analyze(raw, only=None):
    jobs = load_jobs(raw)
    print(f"{len(jobs)} jobs in {raw}")
    summ_rows = {}
    for (c, h), j in sorted(jobs.items(), key=lambda kv: (-kv[0][0], -kv[0][1])):
        if not str(j.get("status", "")).startswith("done"):
            print(f"  c = {c}, h = {h}: {j.get('status')}")
            continue
        orient = j["fate_lo"]
        ps = pstar_of(j["lo"], j["hi"])
        rows = [p for p in j["pts"] if p.get("phase") == "series"]
        allf, prim, env = fit_set(rows, ps, orient)
        tgt, lab = TARGET(c)
        fl = {s: floor_of(rows, ps, orient, prim["slope"], CRIT["windows"][0], s) for s in (-1, 1)}
        tsat = {s: next((r["T"] for r in reversed(rows) if r["side"] == s and r["T"] is not None), None) for s in (-1, 1)}
        inside = env[0] <= tgt <= env[1]
        summ_rows[(c, h)] = dict(ps=ps, prim=prim, env=env, tgt=tgt, lab=lab, inside=inside, fl=fl, tsat=tsat, allf=allf)
        side_fits = {tuple(f["sides"]): f["slope"] for f in allf if f["window"] == CRIT["windows"][0]}
        print(f"  c = {c:+.2f}, h = {h:g}: A*_h = {ps:.6e} ({j['bisect_runs']} runs, extraction t_s = {j['extraction']['t_s']:.2f}, {j['secs']:.0f} s); primary slope {prim['slope']:.5f} +- {prim['se']:.5f} "
              f"(n = {prim['n']}), {lab} = {tgt:.5f}, relative {prim['slope'] / tgt - 1:+.2e}; side alone: collapse-side {side_fits.get((-1,), float('nan')):.5f}, other {side_fits.get((1,), float('nan')):.5f}; "
              f"fit spread [{env[0]:.5f}, {env[1]:.5f}], {lab} inside: {inside}; floor |A - A*|: {fmtf(fl[-1])}, {fmtf(fl[1])}; T at e = {CRIT['e_hi']}: {fmtf(tsat[-1], '{:.2f}')}, {fmtf(tsat[1], '{:.2f}')}")
    # E8
    for c in sorted({k[0] for k in summ_rows}, reverse=True):
        hs = sorted([h for (cc, h) in summ_rows if cc == c], reverse=True)
        if len(hs) < 3:
            continue
        a = [summ_rows[(c, h)]["ps"] for h in hs]
        o = [math.log2(a[i] / a[i + 1]) for i in range(len(a) - 1)]
        rich = a[-1] + (a[-1] - a[-2]) / 3.0
        ok = all(CRIT["order_lo"] <= x <= CRIT["order_hi"] for x in o)
        tag = "E8" if c != -1.0 else "E8c"
        check(f"{tag}[c = {c:g}]", "A*_h converges at second order in h (successive orders within [1.75, 2.25])" + ("" if c != -1.0 else " (the control; reported)"),
              ok if c != -1.0 else None, f"h {hs}: A*_h " + ", ".join(f"{x:.5e}" for x in a) + f"; orders " + ", ".join(f"{x:.3f}" for x in o) + f"; Richardson (order 2) A* = {rich:+.2e}" + ("" if ok else " (outside the band)"))
    # E9
    for c in (-0.3, -1.5):
        hs = sorted([h for (cc, h) in summ_rows if cc == c], reverse=True)
        if len(hs) < 2:
            continue
        h1, h2 = hs[-2], hs[-1]
        a, b = summ_rows[(c, h1)], summ_rows[(c, h2)]
        ok1 = all(abs(x["prim"]["slope"] - x["tgt"]) <= 2 * x["prim"]["se"] for x in (a, b))
        ok2 = all(x["inside"] for x in (a, b))
        ok3 = abs(a["prim"]["slope"] - b["prim"]["slope"]) <= 2 * math.hypot(a["prim"]["se"], b["prim"]["se"])
        check(f"E9[c = {c:g}]", "the lifetime slope equals 1/W within its fit error at the two finest h, which agree: primary window within 2 OLS SE of 1/W, 1/W inside the fit spread, and "
              "the two primary slopes within 2 sqrt(SE1^2 + SE2^2)", ok1 and ok2 and ok3,
              f"h = {h1:g}, {h2:g}: slopes {a['prim']['slope']:.5f} +- {a['prim']['se']:.5f}, {b['prim']['slope']:.5f} +- {b['prim']['se']:.5f}; 1/W = {a['tgt']:.5f}; "
              f"relative {a['prim']['slope'] / a['tgt'] - 1:+.2e}, {b['prim']['slope'] / b['tgt'] - 1:+.2e}; inside the spread: {a['inside']}, {b['inside']}; agreement {a['prim']['slope'] / b['prim']['slope'] - 1:+.2e}")
    # XC: the one-sided pulse (family O, #40's) as the cross-check that the slope does not depend on the family
    pj = load_jobs(raw, "pulse")
    for (c, h), j in sorted(pj.items(), key=lambda kv: (-kv[0][0], -kv[0][1])):
        if not str(j.get("status", "")).startswith("done"):
            print(f"  family O, c = {c}, h = {h}: {j.get('status')}")
            continue
        orient = j["fate_lo"]
        ps = pstar_of(j["lo"], j["hi"])
        rows = [p for p in j["pts"] if p.get("phase") == "series"]
        allf, prim, env = fit_set(rows, ps, orient, CRIT["windows_pulse"])
        tgt, lab = TARGET(c)
        m = summ_rows.get((c, h))
        print(f"  family O (pulse on the {CRIT['pulse_ray'][c]} ray, collapse side, weighted trapped cell), c = {c:+.2f}, h = {h:g}: A*_h = {ps:.6e}; primary slope {prim['slope']:.5f} +- {prim['se']:.5f} (n = {prim['n']}), "
              f"{lab} = {tgt:.5f}, relative {prim['slope'] / tgt - 1:+.2e}; fit spread [{env[0]:.5f}, {env[1]:.5f}]" + ("" if m is None else f"; mode family at the same h: {m['prim']['slope']:.5f}, relative difference {prim['slope'] / m['prim']['slope'] - 1:+.2e}"))
        if h == 0.005 and m is not None:
            check(f"XC[c = {c:g}]", "the one-sided pulse (family O) gives the same slope as the mode family at h = 0.005: within 5 %", abs(prim["slope"] / m["prim"]["slope"] - 1) <= 0.05,
                  f"{prim['slope']:.5f} against {m['prim']['slope']:.5f}")
    # E11
    c = -1.0
    hs = sorted([h for (cc, h) in summ_rows if cc == c], reverse=True)
    if len(hs) >= 3:
        a, b = summ_rows[(c, hs[-2])], summ_rows[(c, hs[-1])]
        o = [math.log2(summ_rows[(c, hs[i])]["ps"] / summ_rows[(c, hs[i + 1])]["ps"]) for i in range(len(hs) - 1)]
        law = (all(x["prim"]["se"] <= 0.10 * x["prim"]["slope"] for x in (a, b)) and abs(a["prim"]["slope"] - b["prim"]["slope"]) <= 2 * math.hypot(a["prim"]["se"], b["prim"]["se"])
               and all(CRIT["order_lo"] <= x <= CRIT["order_hi"] for x in o))
        check("E11", "control c = -1: #212 has no growing mode, so a family through this member should show no divergent lifetime at A = 0. "
              "MET if no lifetime law is found (slope SE > 10 % of the slope, or the two finest slopes disagree, or A*_h not second order); NOT MET if one is found",
              not law, f"law found: {law}; h = {hs[-2]:g}, {hs[-1]:g}: slopes {a['prim']['slope']:.5f} +- {a['prim']['se']:.5f}, {b['prim']['slope']:.5f} +- {b['prim']['se']:.5f}; "
              f"1/kappa = {a['tgt']:.5f}, relative {a['prim']['slope'] / a['tgt'] - 1:+.2e}, {b['prim']['slope'] / b['tgt'] - 1:+.2e}; A*_h orders " + ", ".join(f"{x:.3f}" for x in o))


# --------------------------------------------------------------------------------------------------------------------
# naming runs
# --------------------------------------------------------------------------------------------------------------------
def name_runs(c, h, raw, out):
    j = load_jobs(raw)[(c, h)]
    mem = ep.Member(c)
    X, rc = j["X"], j["rho_c"]
    dR, dL, info = extract_mode(mem, rc, h, X + 25.0, CRIT["Y_ext"], CRIT["tmin_extract"])
    ps = 0.5 * (j["lo"] + j["hi"])
    N = int(round(X / h))
    par, st = params(mem), Static(mem, rc, N, h)
    m_adm = -2.0 * c / (3.0 * math.pi)
    os.makedirs(out, exist_ok=True)
    lines = []
    print(f"naming runs, c = {c}, h = {h}, X = {X}, A*_h = {ps:.6e}; ADM mass of the member m = {m_adm:.5f}; first trapped (T) and anti-trapped (E) cell: raw = dn.cell_class tolerance 1e-2; weighted = signs and 2 e^f r+ r- >= 0.02")
    for e in (4, 5, 6):
        for s in (-1, 1):
            A = ps + s * 10.0 ** (-e)
            d = mode_rays(mem, rc, N, h, A, dR, dL)
            r = evolve_dep(d, par, st, Y0=CRIT["Y0"], stop=False, finders=True, record=True)
            rec, ev, dep = r["rec"], r["events"], r["dep"]
            fmt = lambda x: "none" if x is None else f"t = {x['t']:.2f} l = {x['l']:+.2f} r = {x['r']:.3f}"
            print(f"  A = A* {'-' if s < 0 else '+'} 1e-{e}: departure {None if dep is None else (round(dep['T'], 2), dep['sign'])}; raw T: {fmt(ev['raw_trapped'])}; raw E: {fmt(ev['raw_anti'])}; "
                  f"weighted T: {fmt(ev['w_trapped'])}; weighted E: {fmt(ev['w_anti'])}")
            row = []
            for tt in (10, 20, 30, 40, 50):
                jj = int(tt * SQ2 / (2 * h))
                if jj < N and np.isfinite(rec["r"][jj]):
                    row.append(f"t = {tt}: r {rec['r'][jj]:.3f}, r+ {rec['rp'][jj]:+.2e}, r- {rec['rm'][jj]:+.2e}, 2M/r {2 * rec['M'][jj] / rec['r'][jj]:.4f}, M {rec['M'][jj]:.4f}")
                else:
                    row.append(f"t = {tt}: excised or outside")
            print("      diagonal: " + "; ".join(row))
            lines.append(dict(e=e, side=s, A=A, dep=dep, events=ev))
    json.dump(lines, open(os.path.join(out, f"naming_c{c:g}_h{h:g}.json"), "w"), indent=1, default=float)


# --------------------------------------------------------------------------------------------------------------------
# reported diagnostic: the linear response to a compact pulse, to the end of the grid
# --------------------------------------------------------------------------------------------------------------------
def resp(c, Ap=1e-6, h=0.01, X=40.0, rc=None, out=None, tuned=None):
    mem = ep.Member(c)
    rc = CRIT["RHO_C"].get(c, 1.0) if rc is None else rc
    N = int(round(X / h))
    st, par = Static(mem, rc, N, h), params(mem)
    m_adm = -2.0 * c / (3.0 * math.pi)
    # the static member's Misner-Sharp mass at the diagonal cell and at the edge cells
    def msm(rho):
        Rr = mem.R(rho)
        return 0.5 * Rr * (1.0 - mem.A(rho) * (rho / Rr) ** 2)

    if tuned is not None:                      # D1b (an addition after D1): the background is the mode family at A = A*, so it stays at the member until the floor
        dR, dL, info = extract_mode(mem, rc, h, X + 25.0, CRIT["Y_ext"], CRIT["tmin_extract"])
        sm = json.load(open(os.path.join(tuned, "summary.json")))
        Abg = 0.5 * (sm["lo"] + sm["hi"])

    def one(A):
        d = member_rays(mem, rc, N, h, pulse=dict(A=A, xc=3.0, w=0.5, ray="+"))
        if tuned is not None:
            d.Sp = d.Sp + Abg * dR[:, :N + 1]
            d.Sm = d.Sm + Abg * dL[:, :N + 1]
        return evolve_dep(d, par, st, Y0=1e9, stop=False, record=True)["rec"]
    res = {}
    for A in (Ap, -Ap, 10 * Ap, -10 * Ap):
        res[A] = one(A)
    t = res[Ap]["t"]
    ok = np.isfinite(t)

    def odd(A, key):
        return (res[A][key] - res[-A][key]) / 2.0 / A
    oY, oY10 = odd(Ap, "Y"), odd(10 * Ap, "Y")
    oM, oME = odd(Ap, "M"), odd(Ap, "Medge")
    # static Misner-Sharp mass at the diagonal (rho_c) and at each row's edge cell (l = ledge)
    Ms_diag = msm(rc)
    le = res[Ap]["ledge"]
    good = np.isfinite(le)
    Ms_edge = np.full_like(le, np.nan)
    rho_e, _, _ = tortoise_rho(mem, rc, le[good])
    Ms_edge[good] = msm(rho_e)
    print(("D1b: background = the mode family at A = A*_h = %.6e (tuned). " % Abg if tuned is not None else "D1: background = the untuned member. ")
          + f"c = {c}, rho_c = {rc}, h = {h}, X = {X}, pulse +-{Ap:g} on the + ray (centre 3, width 0.5); member: kappa = {mem.kappa if mem.kappa else float('nan'):.4f}, #212 W = {W212.get(c)}, ADM mass m = {m_adm:.5f}")
    print("   odd part of Y per unit A_p, and its ratio to the same at 10 A_p (1 = linear):")
    tt = np.arange(6, 55, 4.0)
    for T in tt:
        j = int(np.argmin(np.abs(t - T)))
        if not np.isfinite(oY[j]):
            continue
        print(f"      t = {t[j]:5.1f}: Y/A_p {oY[j]:+.3e}  ratio {oY10[j] / oY[j] if oY[j] != 0 else float('nan'):.4f};  dM_diag/A_p {oM[j]:+.3e}  dM_edge/A_p {oME[j]:+.3e} (edge l = {le[j]:+.1f})")
    print("   local rate of ln|Y| (odd part):")
    rates = []
    for a, b in (((6, 10), (10, 14), (14, 18), (18, 22), (22, 26), (26, 30)) if tuned is not None else ((10, 16), (16, 22), (22, 28), (28, 34), (34, 40), (40, 46), (46, 52))):
        m = ok & (t >= a) & (t <= b) & np.isfinite(oY) & (oY != 0)
        if m.sum() > 5:
            rates.append((a, b, float(np.polyfit(t[m], np.log(np.abs(oY[m])), 1)[0])))
    print("      " + "; ".join(f"[{a}, {b}] {s:+.3f}" for a, b, s in rates))
    last = ok & np.isfinite(oY) & (t >= t[ok].max() - 8.0)
    late = float(np.median(oY[last]))
    rate_late = rates[-1][2] if rates else float("nan")
    absY = np.abs(oY)
    sign_changes = int(np.sum(np.diff(np.sign(oY[ok & np.isfinite(oY)])) != 0))
    # saturation time: the first t after which |Y| stays within 10 % of the late value
    within = np.abs(absY - abs(late)) <= 0.1 * abs(late)
    tsat = None
    idx = np.where(ok & np.isfinite(oY))[0]
    for k in idx:
        if within[idx[idx >= k]].all():
            tsat = float(t[k])
            break
    print(f"   late value of Y/A_p (median over the last 8 time units): {late:+.4e}; peak |Y|/A_p {np.nanmax(absY[ok]):.3e} at t = {t[ok][np.nanargmax(absY[ok])]:.1f}; sign changes {sign_changes}; "
          f"rate of ln|Y| over the last window {rate_late:+.3f}; saturated (rate < 0.05 in magnitude over the last window and within 10 % of the late value from t = {tsat}): {abs(rate_late) < 0.05 and tsat is not None}")
    j = int(np.where(ok)[0][-1]) - 20
    print(f"   at the last rows (t = {t[j]:.1f}): M at the diagonal cell, A = +A_p: {res[Ap]['M'][j]:.5f} (static {Ms_diag:.5f}); at the edge cell (l = {le[j]:+.1f}): {res[Ap]['Medge'][j]:.5f} (static {Ms_edge[j]:.5f}); "
          f"member m = {m_adm:.5f}; odd dM_diag/A_p {oM[j]:+.3e}, dM_edge/A_p {oME[j]:+.3e}; 2M/r at the diagonal {2 * res[Ap]['M'][j] / res[Ap]['r'][j]:.4f}")
    if out:
        os.makedirs(out, exist_ok=True)
        np.savez(os.path.join(out, f"resp_c{c:g}.npz"), t=t, oY=oY, oM=oM, oME=oME, le=le)


def bgdep(c, h=0.01, X=40.0, rc=None):
    """the untuned member's own departure: the times at which |Y| on the diagonal first reaches 2e-5, 1e-3, 1e-2 and 1e-1 (D1b: when the background of D1 leaves the member)"""
    mem = ep.Member(c)
    rc = CRIT["RHO_C"].get(c, 1.0) if rc is None else rc
    N = int(round(X / h))
    r = evolve_dep(member_rays(mem, rc, N, h), params(mem), Static(mem, rc, N, h), Y0=1e9, stop=False, record=True)["rec"]
    t, Y = r["t"], np.abs(r["Y"])
    out = []
    for lev in (2e-5, 1e-3, 1e-2, 1e-1):
        m = np.isfinite(Y) & (Y >= lev) & (t >= 4)
        out.append(f"{lev:g}: t = {t[m][0]:.2f}" if m.any() else f"{lev:g}: not reached")
    print(f"c = {c}, h = {h}: the untuned member's |Y| on the diagonal first reaches " + "; ".join(out))


# --------------------------------------------------------------------------------------------------------------------
def main():
    global MUT, ONLY
    args = sys.argv[1:]
    kw = {a.split("=", 1)[0][2:]: a.split("=", 1)[1] for a in args if a.startswith("--") and "=" in a}
    cmd = args[0]
    MUT = kw.get("mutate")
    if MUT is not None and MUT not in MUTATIONS:
        sys.exit(f"unknown mutation {MUT}; known: {', '.join(MUTATIONS)}")
    ONLY = kw["only"].split(",") if "only" in kw else None
    if cmd == "check":
        for cid, fn in QUICK.items():
            if want(cid):
                t0 = time.time()
                try:
                    fn()
                except Exception as e:           # a check that raises under a mutation has fired
                    check(cid, "raised", False, repr(e)[:200])
                print(f"    ({cid}: {time.time() - t0:.1f} s)", flush=True)
        bad = [cid for cid, ok in results if ok is False]
        print(f"\n{len(results)} checks, {sum(1 for _, ok in results if ok)} pass, {len(bad)} fail: {bad}; mutation {MUT}")
        if MUT is not None:
            fired = sorted({cid.split('[')[0] for cid in bad})
            decl = sorted(MUTATIONS[MUT])
            print(f"mutation {MUT}: fired {fired}, declared {decl}: {'equal' if fired == decl else 'DIFFERENT'}")
        sys.exit(1 if bad else 0)
    if cmd == "job":
        c, h = float(kw["c"]), float(kw["h"])
        ok = job(c, h, kw["out"], X=float(kw["X"]) if "X" in kw else None, rc=float(kw["rho_c"]) if "rho_c" in kw else None, Y0=float(kw["Y0"]) if "Y0" in kw else None,
                 family=kw.get("family", "mode"))
        sys.exit(0 if ok else 1)
    if cmd == "analyze":
        analyze(kw["raw"])
        bad = [cid for cid, ok in results if ok is False]
        print(f"\n{len(results)} criteria, {sum(1 for _, ok in results if ok)} met, {len(bad)} not met: {bad}; mutation {MUT}")
        sys.exit(1 if bad else 0)
    if cmd == "name":
        name_runs(float(kw["c"]), float(kw["h"]), kw["raw"], kw["out"])
        return
    if cmd == "bgdep":
        bgdep(float(kw["c"]), float(kw.get("h", 0.01)))
        return
    if cmd == "resp":
        resp(float(kw["c"]), float(kw.get("Ap", 1e-6)), out=kw.get("out"), tuned=kw.get("tuned"))
        return
    sys.exit(f"unknown command {cmd}")


if __name__ == "__main__":
    main()
