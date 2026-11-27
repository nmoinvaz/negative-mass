#!/usr/bin/env python3
"""Horizon-regular gauge test of #217 (design stage diagnosis, item (a) of the coordinator's ruling 4). Prototype, not the build.

The member (c < 0, corner outside the horizon) is evolved with different steps in the two null directions: x+ uniform with step hp (tortoise,
x+ = v/sqrt2), and the second null coordinate either the tortoise y = u/sqrt2 ('tort': z = y, step hm) or the Kruskal-like
z = -exp(-kappa sqrt2 y)/(kappa sqrt2) ('kru': regular on the future horizon z = 0, which the grid approaches but does not reach).
The row index is z, the column index x+. With e^{-f_old} = A: e^{-f_new} = A dy/dz = A exp(kappa sqrt2 y).
The untuned member's departure is followed by Y = dr - (rho R/sqrt2) dphi (cell centre, static values from splines of the exact member,
node-averaged like the code's centre values) on the diagonal l = 0 (x+ = y). Truncation error seeds the unstable mode.

    evolve_dev_reg.py C RHOC GAUGE HP HM TEND [XP]
"""
import math
import os
import sys
import time

import numba as nb
import numpy as np
from scipy.interpolate import CubicSpline

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "wormhole_evolution"))
import dn                      # noqa: E402
import evolve_dev as ed        # noqa: E402
import evolve_potential as ep  # noqa: E402

SQ2 = math.sqrt(2.0)


@nb.njit(cache=True)
def row_step2(old, new, hp, hm, par, eps):
    """dn.row_step's update with different steps: x+ (index i, step hp) and the second coordinate (rows, step hm); no excision."""
    nf, n1 = old.shape
    hh = hp * hm
    ncorr = int(par[4])
    c = np.empty(nf)
    up = np.empty(nf)
    um = np.empty(nf)
    nh = np.empty(nf)
    g = np.empty(nf)
    for i in range(n1 - 1):
        for q in range(nf):
            S = old[q, i]
            E = old[q, i + 1]
            W = new[q, i]
            c[q] = 0.5 * (E + W)
            up[q] = (E - S) / hp
            um[q] = (W - S) / hm
        ed._F(c, up, um, par, eps, g)
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
            ed._F(c, up, um, par, eps, g)
            for q in range(nf):
                nh[q] = old[q, i + 1] + new[q, i] - old[q, i] + hh * g[q]
        for q in range(nf):
            new[q, i + 1] = nh[q]


def static_splines(mem, rho_c, lmax=0.6, n=6001):
    ls = np.linspace(-lmax, lmax, n)
    rho, A, f = ed.tortoise_rho(mem, rho_c, ls)
    R, P = mem.R(rho), mem.phi(rho)
    Q = rho * R / SQ2
    return CubicSpline(ls, R), CubicSpline(ls, P), CubicSpline(ls, Q)


def left_data(mem, rho_c, ys, kap2):
    """r, f_old, phi at the left-ray nodes with x+ = 0, l = -y/sqrt2 (exterior: zeta = ln(rho - rho_h) is integrated)."""
    lmin = -ys.max() / SQ2 - 0.2
    ls = np.linspace(lmin, 0.0, 40001)
    # zeta(l) by the dn-style integration on the left branch, via tortoise_rho on a sorted grid
    from scipy.integrate import solve_ivp
    z0 = rho_c - mem.rh
    sol = solve_ivp(lambda t, y: mem.a_over_z(np.exp(y)), [0.0, lmin], [math.log(z0)], method="DOP853", rtol=1e-13, atol=1e-14, t_eval=ls[::-1])
    zeta = sol.y[0][::-1]
    spl = CubicSpline(ls, zeta)
    zt = spl(-ys / SQ2)
    z = np.exp(zt)
    a = mem.a_over_z(z)
    rho = mem.rh + z
    return mem.R(rho), -zt - np.log(a), mem.phi(rho)


def run(c, rho_c, gauge, hp, hm, tend, Xp=None):
    mem = ep.Member(c)
    kap2 = mem.kappa * SQ2
    yend = tend / SQ2
    Xp = Xp or (yend + 3.0)
    Np = int(round(Xp / hp))
    zc = -1.0 / kap2 if gauge == "kru" else 0.0
    if gauge == "kru":
        zend = -math.exp(-kap2 * yend) / kap2
        Nm = int(math.ceil((zend - zc) / hm))
    else:
        Nm = int(math.ceil(yend / hm))
    js = np.arange(Nm + 2)
    zj = zc + js * hm
    yj = -np.log(-kap2 * zj) / kap2 if gauge == "kru" else zj
    shift = (kap2 * yj) if gauge == "kru" else np.zeros_like(yj)             # f_new = f_old - kappa sqrt2 y
    # data: right ray (row 0, z = zc, y = 0), left ray (column 0, x+ = 0)
    xs = np.arange(Np + 1) * hp / SQ2
    rho, A, f = ed.tortoise_rho(mem, rho_c, xs)
    Sp = np.zeros((3, Np + 1))
    Sp[0], Sp[1], Sp[2] = mem.R(rho), f, mem.phi(rho)
    rl, fl, pl = left_data(mem, rho_c, yj, kap2)
    Sm = np.zeros((3, Nm + 2))
    Sm[0], Sm[1], Sm[2] = rl, fl - shift, pl
    assert abs(Sm[0, 0] - Sp[0, 0]) < 1e-9 and abs(Sm[1, 0] - Sp[1, 0]) < 1e-9, (Sm[:, 0], Sp[:, 0])
    Rs, Ps, Qs = static_splines(mem, rho_c)
    par = ed.params(mem)
    eps = np.array([-1.0])
    old = np.ascontiguousarray(Sp.copy())
    new = np.empty_like(old)
    ts, Ys = [], []
    t0 = time.time()
    for j in range(Nm):
        new[:, 0] = Sm[:, j + 1]
        row_step2(old, new, hp, hm, par, eps)
        ym = 0.5 * (yj[j] + yj[j + 1])
        i = int(round(ym / hp - 0.5))
        if 0 <= i < Np - 1:
            lc = ((i + 0.5) * hp - ym) / SQ2
            nodes_l = [((i + di) * hp - yj[j + dj]) / SQ2 for di in (0, 1) for dj in (0, 1)]
            Rn = np.mean([Rs(x) for x in nodes_l])
            Pn = np.mean([Ps(x) for x in nodes_l])
            cr = 0.25 * (old[0, i] + old[0, i + 1] + new[0, i] + new[0, i + 1])
            cp = 0.25 * (old[2, i] + old[2, i + 1] + new[2, i] + new[2, i + 1])
            Ys.append((cr - Rn) - float(Qs(lc)) * (cp - Pn))
            ts.append((((i + 0.5) * hp) + ym) / SQ2)
        old, new = new, old
    return np.array(ts), np.array(Ys), dict(Nm=Nm, Np=Np, secs=time.time() - t0, mem=mem)


if __name__ == "__main__":
    c, rc, gauge = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
    hp, hm, tend = float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
    Xp = float(sys.argv[7]) if len(sys.argv) > 7 else None
    ts, Ys, info = run(c, rc, gauge, hp, hm, tend, Xp)
    mem = info["mem"]
    W = {-1.5: 0.4941711675, -1.01: 0.6331476235, -1.1: 0.6031905199, -0.99: 0.6401216241}.get(c)
    print(f"c = {c}, rho_c = {rc}, gauge {gauge}, hp = {hp}, hm = {hm}, rows {info['Nm']}, columns {info['Np']}, {info['secs']:.1f} s, kappa = {mem.kappa:.4f}, #212 W = {W}")
    for ta, tb in ((5, 8), (8, 11), (11, 14), (6, 12), (8, 14)):
        m = (ts >= ta) & (ts <= tb) & (np.abs(Ys) > 0)
        if m.sum() > 5:
            print(f"   t in [{ta}, {tb}]: slope of ln|Y| {np.polyfit(ts[m], np.log(np.abs(Ys[m])), 1)[0]:+.4f}, |Y| from {abs(Ys[m][0]):.2e} to {abs(Ys[m][-1]):.2e}, sign {np.sign(Ys[m][0]):+.0f}")
    for tt in (4, 8, 12):
        k = np.argmin(np.abs(ts - tt))
        print(f"   Y(t = {ts[k]:.2f}) = {Ys[k]:+.3e}")
