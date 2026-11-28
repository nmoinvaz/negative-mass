#!/usr/bin/env python3
"""Post hoc description for #314 (not part of the declared checks; written after the analysis of record): how much each reading moves the prediction at the frozen coefficients, and in which direction.

Under the frozen refit of #295 (a, b, c) = P4, the prediction M_pred = r_th (M_0 + a gamma_1 + b gamma_1^2 + c rdot_th). For every reading rho and every point: M_pred(rho) - M_pred(REF) with P4 held fixed,
split into the part from r_th, from gamma_1 and from rdot_th, and the two member shifts that cancel. This says why a member reading that moves (r_th, gamma_1) by 1e-3 moves the prediction by much less.
usage: python3 settled_mass_speed_post.py            (stored runs read from ~/Source/negative-mass-raw)"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.argv = ["x"]
import settled_mass_speed as S            # noqa: E402

P = S.load_points("all")
rd = [S.readings_of(p) for p in P]
a, b, c = S.P4
M0 = S.M0


def pred(g1, rth, rdot):
    return rth * (M0 + a * g1 + b * g1 ** 2 + c * rdot)


print("== post hoc: M_pred(rho) - M_pred(REF) at the frozen P4, and its parts (from r_th, gamma_1, rdot_th), for the member readings and the speed estimators")
print("    point   | WJ: total   from r_th   from g1    | WJi: total  from r_th   from g1    | WA: total   from g1    | S2: total   S3: total   S1: total")
for p, x in zip(P, rd):
    g0, r0, v0 = x["REF"]
    row = f"    {p['name']:6s} |"
    for key in ("WJ", "WJi"):
        g1, r1, v1 = x[key]
        tot = pred(g1, r1, v1) - pred(g0, r0, v0)
        part_r = (r1 - r0) * (M0 + a * g0 + b * g0 ** 2 + c * v0)
        part_g = r0 * (a * (g1 - g0) + b * (g1 ** 2 - g0 ** 2))
        row += f" {tot:+.2e}  {part_r:+.2e}  {part_g:+.2e} |"
    g1, r1, v1 = x["WA"]
    row += f" {pred(g1, r1, v1) - pred(g0, r0, v0):+.2e}  {r0 * (a * (g1 - g0) + b * (g1 ** 2 - g0 ** 2)):+.2e} |"
    row += " " + "  ".join(f"{pred(g0, r0, x[k][2]) - pred(g0, r0, v0):+.2e}" for k in ("S2", "S3", "S1"))
    print(row)
