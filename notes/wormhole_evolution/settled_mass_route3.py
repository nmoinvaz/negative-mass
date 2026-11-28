#!/usr/bin/env python3
"""A frozen refit of family G's settled-mass relation, tested blind on a third route to the member (#295). D = 4, delta = 1e-8.

#266 reads the settled mass of family G along the normal amplitude B, #267 along delta, #283 along two other routes (the normal pulse's width w = 1.0, its centre 4.0). This driver
(1) reads the fifteen members at h = 0.005 that exist (the four B of #266, the five delta of #267, the six of #283) with #283's code (settled_mass_route.py, unchanged),
(2) fits y = M / r_th - M_0 = a gamma_1 + b gamma_1^2 + c rdot_th by plain least squares with M_0 = 0.223086 (the coordinator's ruling) and prints the coefficients and the
leave-one-out errors (`refit`), which are copied into the raw README before any run of the third route (they are the constants of `CRIT3["P4"]`, and the analysis recomputes them from the
stored points and compares), and
(3) reaches gamma_1 = 0.012, 0.048, 0.075 by a third route (the normal pulse's width and centre both changed), runs the chain exactly as #283, and tests the frozen refit on the new members.

A thin driver. It imports settled_mass_route.py (#283), which imports settled_mass_b.py (#266) and settled_mass_member.py (#267), which import settled_mass.py (#41), settled_mass_g.py (#253),
threshold.py (#40), throat_locus.py (#225), dn.py (#39) and throat_drift.py (#255), all unchanged.

    python3 settled_mass_route3.py refit   [--raw=DIR]                               the fifteen points and the refit (reads the stored runs of #266, #267, #283)
    python3 settled_mass_route3.py tag|bracket|pstar|job|member ...                  as settled_mass_route.py (#283), for the third route's points
    python3 settled_mass_route3.py planted [--mutate=NAME]                           planted cases
    python3 settled_mass_route3.py analyze --raw=DIR [--mutate=NAME] [--baseline=FILE] [--rn=DIR]

Units: a = 1 is the throat radius of the initial static wormhole, 8 pi G = 1, D = 4 (n = 2), M = r/2 on the apparent horizon (r_+ = 0), t = (x+ + x-)/sqrt2.
Status words as in AGENTS.md: every mass and every member parameter is numerical evidence.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import settled_mass_route as r         # noqa: E402  (#283's driver; imports #266's and #267's drivers and everything below them unchanged)
b, m, g, th, tl, sm, dn, td = r.b, r.m, r.g, r.th, r.tl, r.sm, r.dn, r.td

SQ2 = math.sqrt(2.0)
MUT = None
ND = 2
FAM = "G"
H2 = r.H2                                       # (0.01, 0.005)
H_FIN = r.H_FIN
RAW_ROOT = os.path.expanduser("~/Source/negative-mass-raw")
I266 = m.I266
I267 = os.path.join(RAW_ROOT, "I267-mac")
I283 = os.path.join(RAW_ROOT, "I283-mac")
M0 = 0.223086                                   # B = 0's plateau at h = 0.005 (#266), the coordinator's ruling for the refit
EPS_DEP = 1e-3


def mut(name):
    return MUT == name


# --------------------------------------------------------------------------------------------------------------------
# the fifteen points
# --------------------------------------------------------------------------------------------------------------------
B_266 = (0.025, 0.05, 0.1, 0.125)               # #266's four B (B = 0 is the symmetric member, M_0)
E_267 = (7, 9, 10, 11, 12)                      # #267's five delta = 1e-e at B = 0.1 (e = 8 is #266's B = 0.1)
POINTS_283 = (("W1", (0.07905, 3.0, 1.0)), ("W2", (0.1543, 3.0, 1.0)), ("W3", (0.1909, 3.0, 1.0)),
              ("C1", (0.0393, 4.0, 0.5)), ("C2", (0.07726, 4.0, 0.5)), ("C3", (0.09576, 4.0, 0.5)))


def read_267(e, h=H_FIN):
    """#267's run at delta = 10^-e, B = 0.1 (stored, read only): plateau, member at departure (eps = 1e-3), throat speed."""
    base = os.path.join(I267, "runs", m.TAG, f"h{h:g}")
    d = dict(e=e, h=h, base=f"{base}/e{e:g}")
    d["s"] = dict(np.load(d["base"] + "/chain/series.npz"))
    d["m"] = json.load(open(d["base"] + "/member/member.json"))
    d["loc"] = dict(np.load(d["base"] + "/member/locus.npz"))
    mem = b.member_at(d, EPS_DEP, "m")
    return dict(M=g.summary_g(d["s"], h)["plat1"], g1=mem["g1"], rth=mem["r_th"], rdot=m.speed_at(d["loc"], mem["t_dep"]), t_dep=mem["t_dep"])


def points15(h=H_FIN):
    """The fifteen points at h: name, route, gamma_1, r_th, rdot_th, M (plateau). Route: B (#266), D (#267's delta scan), W, C (#283)."""
    Rb = m.load266(h)
    out = []
    for B in B_266:
        x = Rb[B]
        out.append(dict(name=f"B{B:g}", route="B", g1=x["mem"]["g1"], rth=x["mem"]["r_th"], rdot=x["speed"], M=x["plat"]))
    for e in E_267:
        x = read_267(e, h)
        out.append(dict(name=f"e{e}", route="D", g1=x["g1"], rth=x["rth"], rdot=x["rdot"], M=x["M"]))
    for nm, spec in POINTS_283:
        e = r.load_point(os.path.join(I283, "runs"), spec, h)
        x = r.read_point(e)
        out.append(dict(name=nm, route=nm[0], g1=x["mem"]["g1"], rth=x["mem"]["r_th"], rdot=x["sp"], M=x["M"]))
    return out, Rb[0.0]["plat"]


# --------------------------------------------------------------------------------------------------------------------
# the fit: plain least squares, M_0 fixed, y = M / r_th - M_0 = a g + b g^2 (+ c rdot)
# --------------------------------------------------------------------------------------------------------------------
def design(g1, rdot, form):
    g1, rdot = np.asarray(g1, float), np.asarray(rdot, float)
    cols = [g1, g1 ** 2] + ([rdot] if form == "P4" else [])
    return np.stack(cols, 1)


def fit(g1, rth, rdot, M, form="P4", m0=M0):
    """Least-squares coefficients of y = M / r_th - m0 on the form's columns (no intercept, no weights)."""
    y = np.asarray(M, float) / np.asarray(rth, float) - m0
    return np.linalg.lstsq(design(g1, rdot, form), y, rcond=None)[0]


def loo(g1, rth, rdot, M, form="P4", m0=M0):
    """Leave-one-out: fit without point k, predict point k. Returns (M - predicted) in M and (y - predicted) in y, per point."""
    g1, rth, rdot, M = (np.asarray(x, float) for x in (g1, rth, rdot, M))
    A = design(g1, rdot, form)
    y = M / rth - m0
    eM, ey = [], []
    for k in range(len(y)):
        idx = [j for j in range(len(y)) if j != k]
        co = np.linalg.lstsq(A[idx], y[idx], rcond=None)[0]
        ey.append(y[k] - A[k] @ co)
        eM.append(M[k] - rth[k] * (m0 + A[k] @ co))
    return np.array(eM), np.array(ey)


def loo_hat(g1, rth, rdot, M, form="P4", m0=M0):
    """The closed form of the leave-one-out error of a linear least-squares fit: e_k / (1 - H_kk), H the hat matrix (an independent route to loo(), in y)."""
    A = design(g1, rdot, form)
    y = np.asarray(M, float) / np.asarray(rth, float) - m0
    H = A @ np.linalg.pinv(A)
    return (y - H @ y) / (1.0 - np.diag(H))


def refit_cmd():
    pts, c0 = points15()
    print(f"== The fifteen points at h = {H_FIN:g} (B = 0's plateau {c0:.7f}; M_0 = {M0} is the ruling)")
    print("    point  route  gamma_1    r_th       rdot_th      M            y = M / r_th - M_0")
    g1, rth, rdot, M = (np.array([p[k] for p in pts]) for k in ("g1", "rth", "rdot", "M"))
    for p in pts:
        print(f"    {p['name']:6s} {p['route']:5s}  {p['g1']:.6f}  {p['rth']:.6f}  {p['rdot']:+.4e}  {p['M']:.8f}  {p['M'] / p['rth'] - M0:+.8f}")
    for form in ("P4", "P2"):
        co = fit(g1, rth, rdot, M, form)
        A = design(g1, rdot, form)
        y = M / rth - M0
        res = rth * (M0 + A @ co) - M
        eM, ey = loo(g1, rth, rdot, M, form)
        eh = loo_hat(g1, rth, rdot, M, form)
        print(f"\n   form {form}: y = " + (" + ".join(f"{c:+.10f} {n}" for c, n in zip(co, ("g", "g^2", "rdot")))))
        print(f"     in-sample M - predicted:  " + ", ".join(f"{p['name']} {-x:+.2e}" for p, x in zip(pts, res)) + f";  max {np.max(np.abs(res)):.2e}")
        print(f"     leave-one-out M - predicted: " + ", ".join(f"{p['name']} {x:+.2e}" for p, x in zip(pts, eM)) + f";  max {np.max(np.abs(eM)):.2e}")
        print(f"     leave-one-out y errors:     " + ", ".join(f"{p['name']} {x:+.2e}" for p, x in zip(pts, ey)) + f";  hat-matrix form agrees to {np.max(np.abs(ey - eh)):.1e}")
        print(f"     singular values of the design {', '.join(f'{x:.3e}' for x in np.linalg.svd(A, compute_uv=False))}")


if __name__ == "__main__":
    argv = sys.argv[1:]
    cmd = argv[0] if argv else "help"
    if cmd == "refit":
        refit_cmd()
    else:
        raise SystemExit(__doc__)
