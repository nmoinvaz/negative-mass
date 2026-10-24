"""Data for Paper III's figure: the growing mode's profile delta r(r) / delta r(a) at several D, and the table of growth
rates (ours for D >= 20 from modes_r.py; Torii & Shinkai's Table I for D <= 11).

    python3 profiles.py        -> profiles.csv, growth.csv

The regular solution is integrated from the throat to r_m = a e^{X/s} (F = 1 - e^{-X}), and continued beyond with the
exact decaying far solution r^{(n+1)/2} K_nu(W r), nu = (n - 1)/2, scaled to match Y at r_m (mpmath for the Bessel
ratio, which overflows in double precision at large nu).
"""
import csv
import os

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

import modes_r as m

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 30

# Torii & Shinkai, PRD 88, 064027 (2013), Table I: -omega^2 a^2 (their n is the spacetime dimension D)
TS13 = {4: 1.39705243371511, 5: 2.98495893027790, 6: 4.68662054299460, 7: 6.46258414126318, 8: 8.28975936306259,
        9: 10.1535530451867, 10: 12.0442650147438, 11: 13.9552091676647, 20: 31.5751101285105,
        50: 91.3457759137153, 100: 191.283017729717}


def profile(D, W2, rgrid, X=12.0):
    n = D - 2
    s = 2 * (n - 1)
    rho0 = 1e-4 / s
    rho_m = np.expm1(X / s)

    def f(rho, z):
        Y, L = z
        r = 1 + rho
        F, Fp = m.F_of(rho, n), m.Fp_of(rho, n)
        return [L, (W2 * Y * r + L * (r * Fp / 2 + n * F) - n * Y * F / r - Y * Fp / 2) / (F * r)]

    inner = rgrid[rgrid - 1 <= rho_m]
    sol = solve_ivp(f, (rho0, rho_m), [m._YS(rho0, n, W2), m._LS(rho0, n, W2)], method="DOP853", rtol=1e-12,
                    atol=1e-30, dense_output=True)
    Y = np.empty_like(rgrid)
    Y[: len(inner)] = sol.sol(np.maximum(inner - 1, rho0))[0]
    Ym = sol.y[0, -1]
    r_m = 1 + rho_m
    W = mp.sqrt(W2)
    nu = mp.mpf(n - 1) / 2
    ref = mp.power(r_m, (n + 1) / 2) * mp.besselk(nu, W * r_m)
    for i in range(len(inner), len(rgrid)):
        r = rgrid[i]
        Y[i] = float(Ym * mp.power(r, (n + 1) / 2) * mp.besselk(nu, W * r) / ref)
    return Y


if __name__ == "__main__":
    rgrid = np.concatenate([1 + np.geomspace(1e-6, 0.05, 200), np.linspace(1.05, 3.0, 400)[1:]])
    rgrid = np.unique(np.concatenate([[1.0 + 1e-7], rgrid]))
    Ds = (10, 20, 100, 1000)
    cols = {}
    for D in Ds:
        W2 = TS13[D] if D in TS13 and D <= 11 else m.W2_of(D)
        if D == 10:
            W2 = TS13[10]
        cols[D] = profile(D, W2, rgrid)
        print(f"D = {D}: W^2 = {W2:.10f}, Y(1) = {cols[D][0]:.6f}, Y(2) = {np.interp(2.0, rgrid, cols[D]):.6f}",
              flush=True)
    limit = rgrid * np.exp(1 - rgrid**2)
    with open(os.path.join(HERE, "profiles.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["r_over_a"] + [f"Y_D{D}" for D in Ds] + ["Y_largeD_limit"])
        for i, r in enumerate(rgrid):
            w.writerow([f"{r:.10f}"] + [f"{cols[D][i]:.10e}" for D in Ds] + [f"{limit[i]:.10e}"])
    rows = [(D, TS13[D], "TS13 Table I") for D in range(4, 12)]
    for D in (20, 30, 50, 100, 200, 500, 1000, 2000):
        rows.append((D, m.W2_of(D), "modes_r.py"))
    with open(os.path.join(HERE, "growth.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["D", "W2a2", "source"])
        for D, v, src in rows:
            w.writerow([D, f"{v:.12f}", src])
    print("wrote profiles.csv, growth.csv")
