#!/usr/bin/env python3
"""Time-domain check of #212's master equation, without the double-null machinery (#217, design stage).

    u_tt = u_xx - Vhat(x) u,   Vhat = 2 A^2/rho^2 - (1 + 1/R^2) A A'/rho + A^2/R^4 - 2 A/R^2,   dx = d rho / A,   u = Y/rho   (README.md, Sec. 2)

Compact initial data, leapfrog, both ends far away. For the members with no throat in the static region (c <= -1) Vhat is smooth. The growth
rate of ln|u| at x = 0 is compared with #212's W. At c = -1 #212 finds no growing mode; the double-null prototype (evolve_dev_runs.py rates)
finds growth at 0.63-0.64 per unit time from truncation seeds. This script says whether the linear master equation grows at c = -1.

    evolve_master1d.py C [--L=90] [--dx=0.02] [--T=60]
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evolve_dev as ed        # noqa: E402
import evolve_potential as ep  # noqa: E402

SQ2 = math.sqrt(2.0)
W212 = {-1.01: 0.6331476235, -1.1: 0.6031905199, -1.5: 0.4941711675, -2.0: 0.3984523430, -3.0: 0.2827745372, -1.001: 0.6362712254}


def vhat_table(mem, rho_c, xs):
    rho, A, f = ed.tortoise_rho(mem, rho_c, xs)
    R2 = rho * rho + 1.0
    Ap = mem.Ap(rho)
    V = 2 * A * A / rho ** 2 - (1 + 1 / R2) * A * Ap / rho + A * A / R2 ** 2 - 2 * A / R2
    return V, rho


def run(c, L=90.0, dx=0.02, T=60.0, rho_c=1.0, x0=2.0, w=0.7):
    mem = ep.Member(c)
    n = int(round(L / dx))
    xs = np.arange(-n, n + 1) * dx
    V, rho = vhat_table(mem, rho_c, xs)
    u = np.exp(-((xs - x0) / w) ** 2)
    dt = 0.5 * dx
    um = u.copy()                                  # u_t = 0 at t = 0
    lap = np.zeros_like(u)
    lap[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / dx ** 2
    up = u + 0.5 * dt * dt * (lap - V * u)
    i0 = n
    ts, ys = [], []
    steps = int(T / dt)
    for s in range(1, steps):
        lap[1:-1] = (up[2:] - 2 * up[1:-1] + up[:-2]) / dx ** 2
        un = 2 * up - u + dt * dt * (lap - V * up)
        u, up = up, un
        if s % 25 == 0:
            ts.append(s * dt)
            ys.append(abs(up[i0]))
    return np.array(ts), np.array(ys), mem


def main():
    c = float(sys.argv[1])
    kw = {a.split("=")[0][2:]: float(a.split("=")[1]) for a in sys.argv[2:] if a.startswith("--")}
    ts, ys, mem = run(c, **kw)
    print(f"c = {c}: kappa = {mem.kappa:.4f}, #212 W = {W212.get(c)}")
    for t0, t1 in ((10, 20), (20, 30), (30, 45), (45, 60)):
        m = (ts >= t0) & (ts <= t1) & (ys > 0)
        if m.sum() > 5:
            print(f"   slope of ln|u(x = 0)| over t in [{t0}, {t1}]: {np.polyfit(ts[m], np.log(ys[m]), 1)[0]:+.4f}   |u| at the ends: {ys[m][0]:.2e} -> {ys[m][-1]:.2e}")


if __name__ == "__main__":
    main()
