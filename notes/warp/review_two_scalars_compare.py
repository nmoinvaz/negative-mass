#!/usr/bin/env python3
"""Comparison of the review's results with two_scalars_output.txt (issue #79), made after the review's own runs.

Reads the printed output of two_scalars.py (not its code) and compares what it prints with what review_two_scalars_ideals.py
derives: the coefficients of 16 r^3 D at constant speed and the odd coefficients of the accelerating bubble, the sample
points of the tanh and bump profiles, the eigenvalues on the axis, the Ellis-Bronnikov b_l, and the b(rho) block.

    python3 notes/warp/review_two_scalars_compare.py                  checks, then every control
    python3 notes/warp/review_two_scalars_compare.py --no-controls    checks only
    python3 notes/warp/review_two_scalars_compare.py --mutate=NAME    one control (its target checks must fail, exit 1)
"""
import importlib.util
import os
import random
import re
import subprocess
import sys

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OPTS = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
MUT = OPTS.get("mutate", "")
CONTROLS = {
    "coef-sign": ("my c^2 coefficient with the opposite sign", ["compare_coeffs"]),
    "odd-sign": ("my c^3 coefficient with the opposite sign", ["compare_odd"]),
    "axis-eps": ("the on-axis points evaluated at c = 99/100", ["compare_points"]),
    "eig-v2": ("v doubled in my eigenvalues", ["compare_eigs"]),
    "eb-kappa": ("kappa doubled in my b_l", ["compare_eb"]),
}
if MUT and MUT not in CONTROLS:
    sys.exit(f"unknown mutation {MUT}")
RESULTS = []
mp.mp.dps = 50


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def load(name):
    saved = sys.argv
    sys.argv = [saved[0]]
    try:
        spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        sys.argv = saved
    return mod


def main():
    idl = load("review_two_scalars_ideals")
    r, c, s, w, w1, w2, aa = idl.r, idl.c, idl.s, idl.w, idl.w1, idl.w2, idl.aa
    txt = open(os.path.join(HERE, "two_scalars_output.txt"), encoding="utf-8").read().replace("\r", "")
    names = {"W0": w, "W1": w1, "W2": w2, "r": r, "R": r, "A": aa, "c": c}
    theirs = lambda e: sp.sympify(e, locals=names)

    # my 16 r^3 D at constant speed and for the accelerating bubble, as in review_two_scalars_ideals.py
    rho = sp.Symbol("rho", positive=True)
    P0, Q0, Zr0, Yr0, h0 = idl.blocks_comoving(0)
    Pm, Qm, hm = (idl.reduce_s(idl.to_slice(e, rho)) for e in (P0, Q0, h0))
    Zm = idl.reduce_s(idl.to_slice(Zr0, rho) / (r * s))
    Ym = idl.reduce_s(idl.to_slice(Yr0, rho) / (r * s))
    D16 = sp.expand(sp.cancel(16 * r**3 * idl.reduce_s(idl.D_of(Pm, Qm, Zm, Ym, hm))))
    pc = sp.Poly(D16, c)
    mine = {k: pc.coeff_monomial(c**k) for k in range(7)}
    if MUT == "coef-sign":
        mine[2] = -mine[2]
    Pa, Qa, Za, Ya, ha = idl.blocks_comoving(aa)
    Pm2, Qm2, hm2 = (idl.reduce_s(idl.to_slice(e, rho)) for e in (Pa, Qa, ha))
    Zm2 = idl.reduce_s(idl.to_slice(Za, rho) / (r * s))
    Ym2 = idl.reduce_s(idl.to_slice(Ya, rho) / (r * s))
    D16a = sp.expand(sp.cancel(16 * r**3 * idl.reduce_s(idl.D_of(Pm2, Qm2, Zm2, Ym2, hm2))))
    pa = sp.Poly(D16a, c)
    odd_mine = [pa.coeff_monomial(c**k) for k in (1, 3, 5)]
    if MUT == "odd-sign":
        odd_mine[1] = -odd_mine[1]

    # 1. the printed coefficients at constant speed
    got = {}
    for k in (0, 2, 4, 6):
        m = re.search(r"^\s+c\^%d: (.*)$" % k, txt, re.M)
        got[k] = theirs(m.group(1))
    diffs = {k: sp.expand(got[k] - mine[k]) for k in got}
    check("compare_coeffs", all(v == 0 for v in diffs.values()),
          "c^0, c^2, c^4, c^6 printed in two_scalars_output.txt section 4 minus mine: %s" % {k: (0 if v == 0 else "nonzero") for k, v in diffs.items()})

    # 2. the printed odd coefficients of the accelerating bubble
    m = re.search(r"Odd coefficients: (\[.*\])", txt)
    lst = sp.sympify(m.group(1), locals=names)
    dodd = [sp.expand(a_ - b_) for a_, b_ in zip(lst, odd_mine)]
    check("compare_odd", all(v == 0 for v in dodd), "c^1, c^3, c^5 printed in section 6 minus mine: %s"
          % ["0" if v == 0 else "nonzero" for v in dodd])

    # 3. the sample points of the tanh and bump profiles and the random-jet points
    def tanh_jets(rv, sig=8, Rb=1):
        f = lambda x: (mp.tanh(sig * (x + Rb)) - mp.tanh(sig * (x - Rb))) / (2 * mp.tanh(sig * Rb))
        return [f(rv), mp.diff(f, rv), mp.diff(f, rv, 2)]

    def bump_jets(rv, Rb=2):
        f = lambda x: (1 - x**2 / Rb**2) ** 3
        return [f(rv), mp.diff(f, rv), mp.diff(f, rv, 2)]

    def D_at(jets, v, rv, cv):
        sub = {w: v * (1 - jets[0]), w1: -v * jets[1], w2: -v * jets[2], r: rv, c: cv}
        return mp.mpf(sp.N(D16.subs(sub), 40))

    rows = []
    q = lambda a_, b_: sp.Rational(a_, b_)
    cs_tanh = [q(12, 13), q(-12, 13), 1, -1, 0]
    for cv in cs_tanh:
        for v in (sp.Integer(1), q(1, 2), sp.Integer(3)):
            rows.append(("tanh", cv, v, D_at(tanh_jets(mp.mpf(1)), v, 1, cv), abs(cv) == 1))
    for cv in [q(4, 5), q(-3, 5), 1, 0]:
        rows.append(("bump", cv, q(7, 4), D_at(bump_jets(mp.mpf(1)), q(7, 4), 1, cv), cv == 1))
    # random jets: the draws of two_scalars.py (random.seed(57), then v and three jets per point, in that order)
    random.seed(57)
    for cv, r0 in [(q(12, 13), q(13, 7)), (q(-12, 13), q(13, 7)), (q(4, 5), q(5, 7)), (q(-3, 5), q(5, 7))]:
        v = sp.Rational(random.randint(1, 9), 3)
        jets = [sp.Rational(random.randint(-9, 9), random.randint(1, 9)) for _ in range(3)]
        rows.append(("random", cv, v, mp.mpf(sp.N(D16.subs({w: v * (1 - jets[0]), w1: -v * jets[1], w2: -v * jets[2], r: r0, c: cv}), 40)), False))
    ok_pts = True
    for lab, cv, v, Dv, on_axis in rows:
        if MUT == "axis-eps" and on_axis:
            Dv = D_at(tanh_jets(mp.mpf(1)), v, 1, sp.Rational(99, 100)) if lab == "tanh" else D_at(bump_jets(mp.mpf(1)), v, 1, sp.Rational(99, 100))
        zero = abs(Dv) < mp.mpf(10) ** -30 * max(1, abs(Dv))
        if zero != on_axis:
            ok_pts = False
    n_off = sum(not x[4] for x in rows)
    check("compare_points", ok_pts, "my D at the sample points of two_scalars.py (tanh sigma = 8 R = 1, bump R_b = 2, and the random-jet draws "
          "of seed 57): D != 0 at the %d off-axis points and D = 0 on the %d on-axis points, as its section 5 reports "
          "(16 and 7 in all, with 4 random-jet points reproduced by draw order)" % (n_off, len(rows) - n_off))

    # 4. eigenvalues on the axis at x - v t = 1, tanh profile
    eigs_theirs = [sp.sympify(m_).evalf(6) for m_ in re.findall(r"eigenvalues of eta \(a a - b b\): \[([^\]]*)\]", txt)[0].split(", ")]
    ok_e = True
    for v, line in ((sp.Integer(1), 0), (q(1, 2), 1)):
        vv = 2 * v if MUT == "eig-v2" else v
        j = tanh_jets(mp.mpf(1))
        W0, W1, W2 = vv * (1 - j[0]), -vv * j[1], -vv * j[2]
        R00 = -(W0 * W2 + W1**2)  # c = 1, s = 0 in the table: R_00 = -(w w'' + w'^2), R_01 = -w'/r
        R01 = -W1
        want = (-R00, abs(R01))
        printed = re.findall(r"eigenvalues of eta \(a a - b b\): \[([^\]]*)\]", txt)[line].split(", ")
        z1 = complex(sp.sympify(printed[1]).evalf())
        ok_e &= abs(z1.real - float(want[0])) < 1e-3 and abs(abs(z1.imag) - float(want[1])) < 1e-3
    check("compare_eigs", ok_e, "eigenvalues -R_00 +- i R_01 on the axis from my table at x - v t = 1: 16 +- 4i (v = 1) and 4 +- 2i (v = 1/2), "
          "as printed")

    # 5. Ellis-Bronnikov: their printed b_l against mine, 50 digits at 6 random parameter sets
    m = re.search(r"status pass, lam = (.*?), a = (\[.*?\]), b = (\[.*\])\n", txt)
    ka, be, ll = sp.symbols("kappa beta l", positive=True)
    bl_theirs = sp.sympify(m.group(3), locals={"kappa": ka, "beta": be, "l": ll})[1]
    rng = random.Random(7918)
    ok_b = sp.sympify(m.group(1), locals={"kappa": ka}) == 2 * ka**2
    for _ in range(6):
        pv = {ka: sp.Rational(rng.randint(1, 30), 10), be: sp.Rational(rng.randint(1, 49), 100), ll: sp.Rational(rng.randint(-10, 10), 10)}
        R2 = (1 - be * sp.cos(2 * ka * ll)) / (2 * ka**2)
        kk = 2 * ka if MUT == "eb-kappa" else ka
        bl_mine = sp.sqrt(1 - be**2) / (sp.sqrt(2) * kk * R2)
        ok_b &= abs(sp.N((bl_theirs - bl_mine).subs(pv), 40)) < sp.Float(10) ** -30
    check("compare_eb", ok_b, "printed lambda = 2 kappa^2 and b_l equal mine, sqrt(1 - beta^2)/(sqrt2 kappa r^2), at 6 random (kappa, beta, l)")

    # 6. the b(rho) block
    m = re.search(r"R_3 = (\[\[.*\]\])\n", txt)
    b001, b002, rr_ = sp.symbols("b001 b002 rho")
    R3t = sp.Matrix(sp.sympify(m.group(1), locals={"b001": b001, "b002": b002, "rho": rr_}))
    hh = b001**2 / 2
    QQ = -(b002 + b001 / rr_) / 2
    R3m = sp.Matrix([[-hh, QQ, 0], [QQ, -hh, 0], [0, 0, hh]])
    check("compare_brho", sp.simplify(R3t - R3m) == sp.zeros(3), "printed R_3 for b(rho) equals [[-h, Q, 0], [Q, -h, 0], [0, 0, h]] with "
          "h = b'^2/2, Q = -(b'' + b'/rho)/2")

    fails = [n for n, ok in RESULTS if not ok]
    print("\nchecks: %d passed, %d failed %s" % (len(RESULTS) - len(fails), len(fails), fails))
    if MUT:
        fired = all(n in fails for n in CONTROLS[MUT][1])
        print("control %s: targets %s %s" % (MUT, CONTROLS[MUT][1], "FIRED" if fired else "SILENT"))
        sys.exit(1)
    status = 1 if fails else 0
    if "no-controls" not in OPTS:
        print("\n== Controls: each mutation must make its target checks fail ==")
        for name, (what, target) in CONTROLS.items():
            p = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
            failed = [ln.split("]")[1].split(":")[0].strip() for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
            fired = p.returncode == 1 and all(n in failed for n in target)
            status |= 0 if fired else 1
            print("[%s] --mutate=%s (%s): exit %d, failed %s" % ("FIRED" if fired else "SILENT", name, what, p.returncode, failed), flush=True)
        print("controls: %s" % ("all fired" if status == 0 else "SOME SILENT OR CHECK FAILED"))
    sys.exit(status)


if __name__ == "__main__":
    main()
