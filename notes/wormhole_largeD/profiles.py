"""Data for Paper III's figure: the growing mode's profile delta r(r) / delta r(a) at several D, and the table of growth
rates (D >= 20: the quad-double values of #195's crate, from w2_output.txt, rounded to 12 decimals; D <= 11: Torii &
Shinkai's Table I).

    python3 profiles.py                -> profiles.csv, growth.csv
    python3 profiles.py --growth-only  -> growth.csv
    python3 check_growth.py            checks growth.csv against w2_output.txt (G1-G5, R0-R2; --mutate=NAME)

The regular solution is integrated from the throat to r_m = a e^{X/s} (F = 1 - e^{-X}), and continued beyond with the
exact decaying far solution r^{(n+1)/2} K_nu(W r), nu = (n - 1)/2, scaled to match Y at r_m (mpmath for the Bessel
ratio, which overflows in double precision at large nu). The profiles' W^2 at D = 20, 100, 1000 is modes_r.py's; the
growth.csv rows do not use it.
"""
import csv
import os
import re
import sys
from decimal import ROUND_HALF_EVEN, Decimal

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


W2_OUTPUT = os.path.join(HERE, "w2_output.txt")
SRC_TS13 = "TS13 Table I"
SRC_CRATE = "w2_crate qd (#195)"
D_CRATE = (20, 30, 50, 100, 200, 500, 1000, 2000)


def crate_w2(path=W2_OUTPUT):
    """{D: string} of the crate's qd W^2 a^2 (about 20 digits) in the table "growth.csv against the crate (qd)"."""
    out, inside = {}, False
    for line in open(path):
        if "growth.csv against the crate (qd)" in line:
            inside = True
        elif inside and not line.strip() and out:
            break
        elif inside and line.split() and line.split()[0].isdigit():
            tok = line.split()
            val = [t for t in tok[1:] if re.fullmatch(r"\d+\.\d{15,}", t)]
            assert len(val) == 1, line
            out[int(tok[0])] = val[0]
    return out


def round12(text):
    """12 decimals of a decimal string, half to even, as a string (no float)."""
    return format(Decimal(text).quantize(Decimal("1e-12"), rounding=ROUND_HALF_EVEN), "f")


def write_growth():
    crate = crate_w2()
    rows = [(D, f"{TS13[D]:.12f}", SRC_TS13) for D in range(4, 12)]
    rows += [(D, round12(crate[D]), SRC_CRATE) for D in D_CRATE]
    with open(os.path.join(HERE, "growth.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["D", "W2a2", "source"])
        w.writerows(rows)
    print("wrote growth.csv")


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


def write_profiles():
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
    print("wrote profiles.csv")


if __name__ == "__main__":
    if "--growth-only" not in sys.argv[1:]:
        write_profiles()
    write_growth()
