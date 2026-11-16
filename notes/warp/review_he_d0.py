#!/usr/bin/env python3
"""The type on D = 0 off the equator (issue #71, an extension of the review).

he_type.md lists "the type on D = 0 off the equator" as not covered. On the equator it is settled exactly (review_he_
derive.py, `Qeq0_and_edges`). Off it, D = 0 has irrational points, so this script samples the surface D = 0 at 50
digits: random rational (v, F, a, u, c, s), then b a real root of D (a polynomial in b with rational coefficients),
then the rank of M - lam on the double root lam of the cubic, by singular values.

Planted truth first: a perfect fluid (rank 1, Type I) and a Jordan block (rank 2, Type II) get their rank from the same
routine. The result is numerical evidence about a sample, not a proof.

    python3 notes/warp/review_he_d0.py                 # exit 1 if a check fails
    python3 notes/warp/review_he_d0.py --mutate=NAME   # NAME: rank-tol, off-root
"""
import argparse
import os
import random
import sys

import mpmath as mp
import sympy as sp
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import review_he_derive as rd  # noqa: E402  (this review's own code)

mp.mp.dps = 50
MUT = set()
CONTROLS = {'rank-tol': ('the rank threshold is 10^6, so every matrix has rank 0', ['planted_rank', 'D0_off_equator_type']),
            'off-root': ('b is moved by 10^-6 off the root of D', ['D0_is_zero'])}
TOL = mp.mpf(10) ** -25


def ev_ring(p, vals):
    """Value of a ring element at mp numbers vals = (v, F, a, b, u, c, s)."""
    tot = mp.mpf(0)
    for mon, co in p.terms():
        t = mp.mpf(co.numerator) / mp.mpf(co.denominator)
        for x, e in zip(vals, mon):
            if e:
                t *= x ** e
        tot += t
    return tot


def double_root(c2, c1, c0):
    return (9 * c0 - c1 * c2) / (2 * c2 ** 2 - 6 * c1)


def rank_at_double_root(M):
    """(rank of M - lam, lam, smallest and middle singular values) for a 3x3 mp matrix with a double root."""
    c2 = -(M[0, 0] + M[1, 1] + M[2, 2])
    c1 = sum(M[i, i] * M[j, j] - M[i, j] * M[j, i] for i, j in ((0, 1), (0, 2), (1, 2)))
    det = (M[0, 0] * (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]) - M[0, 1] * (M[1, 0] * M[2, 2] - M[1, 2] * M[2, 0])
           + M[0, 2] * (M[1, 0] * M[2, 1] - M[1, 1] * M[2, 0]))
    lam = double_root(c2, c1, -det)
    A = M - lam * mp.eye(3)
    sv = sorted(mp.svd_r(A, compute_uv=False), reverse=True)
    tol = mp.mpf(10) ** 6 if 'rank-tol' in MUT else TOL * max(1, sv[0])
    return sum(1 for x in sv if x > tol), lam, sv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mutate', default=None)
    ap.add_argument('--no-controls', action='store_true')
    args = ap.parse_args()
    if args.mutate:
        if args.mutate not in CONTROLS:
            sys.exit('unknown mutation')
        MUT.add(args.mutate)
    rows = []

    def check(name, ok, detail=''):
        rows.append((name, bool(ok)))
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)

    # planted truth for the rank routine: M - lam has rank 1 for a fluid and for dust, rank 2 for a Jordan block
    planted = {
        'perfect fluid, M = diag(-3/2, 1/5, 1/5)': (mp.matrix([[-1.5, 0, 0], [0, 0.2, 0], [0, 0, 0.2]]), 1),
        'dust, M = diag(-2, 0, 0)': (mp.matrix([[-2, 0, 0], [0, 0, 0], [0, 0, 0]]), 1),
        'Jordan 2x2 (eigenvalue 1) plus a simple 7': (mp.matrix([[0, -1, 0], [1, 2, 0], [0, 0, 7]]), 2),
        'Jordan 2x2 (eigenvalue 2) plus a simple -3': (mp.matrix([[1, -1, 0], [1, 3, 0], [0, 0, -3]]), 2),
    }
    bad = [nm for nm, (M, want) in planted.items() if rank_at_double_root(M)[0] != want]
    check('planted_rank', not bad, f'rank(M - lam) is 1 for a fluid and for dust and 2 for a Jordan block; wrong: {bad}')

    # random points of D = 0 off the equator
    mine = rd.ring_frame(rd.build_C())
    Rr = mine
    rng = random.Random(7)
    D = rd.red(rd.Dof(Rr))
    Bs = sp.Symbol('B')
    counts, minsep, nj0, tried = {}, mp.mpf(10) ** 9, 0, 0
    resid = []
    while sum(counts.values()) < 200 and tried < 4000:
        tried += 1
        pt = rd.randpt(rng)
        if pt[rd.c] == 0:
            continue
        vals_q = [pt[x] for x in rd.SYMS]
        # D as a polynomial in b with rational coefficients
        co = {}
        for mon, cf in D.terms():
            t = Fr(cf.numerator, cf.denominator)
            for k, (x, e) in enumerate(zip(vals_q, mon)):
                if k != 3 and e:
                    t *= x ** e
            co[mon[3]] = co.get(mon[3], 0) + t
        P = sp.Poly(sum(sp.Rational(x.numerator, x.denominator) * Bs ** e for e, x in co.items()), Bs)
        if P.degree() < 1:
            continue
        for root in P.nroots(n=50, maxsteps=500):
            if abs(sp.im(root)) > 1e-40:
                continue
            bstar = mp.mpf(str(sp.re(root)))
            if 'off-root' in MUT:
                bstar += mp.mpf(10) ** -6
            vals = [mp.mpf(x.numerator) / x.denominator for x in vals_q]
            vals[3] = bstar
            M = mp.matrix(3, 3)
            for i in range(3):
                for j in range(3):
                    M[i, j] = (-1 if i == 0 else 1) * ev_ring(Rr[i][j], vals)
            c2 = -(M[0, 0] + M[1, 1] + M[2, 2])
            c1 = sum(M[i, i] * M[j, j] - M[i, j] * M[j, i] for i, j in ((0, 1), (0, 2), (1, 2)))
            c0 = -mp.det(M)
            Dnum = 18 * c2 * c1 * c0 - 4 * c2 ** 3 * c0 + c2 ** 2 * c1 ** 2 - 4 * c1 ** 3 - 27 * c0 ** 2
            scale = max(abs(c2), abs(c1) ** mp.mpf(1) / 1, 1) ** 6
            resid.append(abs(Dnum) / scale)
            if 2 * c2 ** 2 - 6 * c1 == 0:
                continue
            jj = abs(M[0, 1]) + abs(M[0, 2])
            rk, lam, sv = rank_at_double_root(M)
            key = {1: 'I', 2: 'II', 3: '?'}.get(rk, 'I' if rk == 0 else '?')
            if jj < TOL:
                nj0 += 1
                continue
            counts[key] = counts.get(key, 0) + 1
            if len(sv) == 3 and sv[1] > 0:
                minsep = min(minsep, sv[1] / sv[0] if rk == 2 else mp.mpf(10) ** 9)
    check('D0_is_zero', resid and max(resid) < mp.mpf(10) ** -30,
          f'D (relative to its scale) at the {len(resid)} sampled roots, worst {mp.nstr(max(resid), 3)}')
    check('D0_off_equator_type', counts.get('II', 0) == sum(counts.values()) and sum(counts.values()) > 100,
          f'rank of M - lam at {sum(counts.values())} random points of D = 0 with c != 0 and j != 0: {counts}; '
          f'smallest ratio of the middle to the largest singular value {mp.nstr(minsep, 3)}; points with j = 0 '
          f'skipped: {nj0}')
    bad = sorted(n for n, ok in rows if not ok)
    print(('FAILED: ' + ','.join(bad)) if bad else f'all {len(rows)} checks passed')
    if args.mutate or args.no_controls:
        sys.exit(1 if bad else 0)
    import subprocess
    code = 1 if bad else 0
    print('\ncontrols (each is the same script with --mutate=NAME and must fail its target checks)')
    for name, (what, targets) in CONTROLS.items():
        p = subprocess.run([sys.executable, '-u', __file__, f'--mutate={name}'], capture_output=True, text=True)
        line = [ln for ln in p.stdout.splitlines() if ln.startswith('FAILED:')]
        got = set(line[-1][len('FAILED: '):].split(',')) if line else set()
        fired = p.returncode == 1 and set(targets) <= got
        print(f"[{'FIRED' if fired else 'DID NOT FIRE'}] --mutate={name}: {what}; exit {p.returncode}; "
              f"checks that fail: {', '.join(sorted(got)) or 'none'}", flush=True)
        code = code or (0 if fired else 1)
    sys.exit(code)


if __name__ == '__main__':
    main()
