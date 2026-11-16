#!/usr/bin/env python3
"""Comparison of review_he_derive.py's results with the printed output of he_type.py (issue #71).

Written and run only after review_he_derive.py had produced its result, which is the order recorded in review_he.md.
It parses the frame Ricci entries, c2, p8 and p10 that notes/warp/he_type_output.txt prints (the output of the code
under review, parsed as data, not imported), converts them from (f, f1, f2, r) to this review's (F = 1 - f, a, b,
u = 1/r), and compares them exactly with the entries and the coefficients of D derived here by Cartan's route.

    python3 notes/warp/review_he_compare.py                  # the comparisons, then every control as a subprocess
    python3 notes/warp/review_he_compare.py --no-controls    # the comparisons alone
    python3 notes/warp/review_he_compare.py --mutate=NAME    # NAME: c-no-nshift (my side), perturb-p8, perturb-R01
"""
import argparse
import os
import re
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import review_he_derive as rd  # noqa: E402  (this review's own code, not the code under review)

OUT = os.path.join(HERE, 'he_type_output.txt')
CONTROLS = {'c-no-nshift': 'my route C loses the shift in E_0 = d_t - w d_x',
            'perturb-p8': 'the parsed p8 gains the term a^2 u^2',
            'perturb-R01': 'the parsed R_01 gains the term v a u / 7'}
v, F, a, b, u, c, s = rd.SYMS
SYM_THEIRS = {n: sp.Symbol(n) for n in ('v', 'c', 's', 'f', 'f1', 'f2', 'r')}


def theirs_to_ring(text, extra=0):
    e = sp.sympify(text, locals=SYM_THEIRS)
    e = e.xreplace({SYM_THEIRS['f']: 1 - F, SYM_THEIRS['f1']: a, SYM_THEIRS['f2']: b, SYM_THEIRS['r']: 1 / u})
    return rd.to_ring(e + extra)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mutate', default=None)
    ap.add_argument('--no-controls', action='store_true')
    args = ap.parse_args()
    mutate = args.mutate
    if mutate == 'c-no-nshift':
        rd.MUT.add('c-no-nshift')
    elif mutate not in (None, 'perturb-p8', 'perturb-R01'):
        sys.exit('unknown mutation')
    text = open(OUT).read()
    mine = rd.ring_frame(rd.build_C())
    rows = []

    def check(name, ok):
        rows.append((name, bool(ok)))
        print(f"[{'PASS' if ok else 'FAIL'}] {name}", flush=True)

    # frame Ricci entries
    key = {'0': 0, '1': 1, 'rho': 2, 'varphi': 3}
    ent = dict((m.group(1) + ',' + m.group(2), m.group(3)) for m in re.finditer(
        r'^  R_\((\w+),(\w+)\) = (.*)$', text, re.M))
    assert len(ent) == 10, ent.keys()
    for kk, expr in ent.items():
        i, j = (key[x] for x in kk.split(','))
        extra = v * a * u / 7 if (mutate == 'perturb-R01' and kk == '0,1') else 0
        check(f'R_({kk}) as printed by he_type.py equals the entry derived here',
              rd.req(theirs_to_ring(expr, extra), mine[i][j]))

    M = rd.block(mine)
    c2, c1, c0 = rd.cubic(M)
    m2 = re.search(r'^  c2 = (.*)$', text, re.M).group(1)
    check('c2 as printed equals the derived c2', rd.req(theirs_to_ring(m2), rd.red(c2)))
    D = rd.red(rd.Dof(mine))
    for k in (8, 10):
        m = re.search(r'^  p%d \(coefficient of v\^%d in D\) = (.*)$' % (k, k), text, re.M).group(1)
        extra = u ** 2 * a ** 2 if (mutate == 'perturb-p8' and k == 8) else 0
        check(f'p{k} as printed equals the derived p{k} ({len(rd.coef_v(D, k))} terms)',
              rd.req(theirs_to_ring(m, extra), rd.coef_v(D, k)))
    bad = [n for n, ok in rows if not ok]
    print(('FAILED: ' + '; '.join(bad)) if bad else f'all {len(rows)} comparisons passed')
    if mutate or args.no_controls:
        sys.exit(1 if bad else 0)
    import subprocess
    code = 1 if bad else 0
    print('\ncontrols (each is the same script with --mutate=NAME and must make a comparison fail)')
    for name in CONTROLS:
        p = subprocess.run([sys.executable, '-u', __file__, f'--mutate={name}'], capture_output=True, text=True)
        nfail = sum(1 for ln in p.stdout.splitlines() if ln.startswith('[FAIL]'))
        fired = p.returncode == 1 and nfail > 0
        print(f"[{'FIRED' if fired else 'DID NOT FIRE'}] --mutate={name}: {CONTROLS[name]}; exit {p.returncode}; "
              f"{nfail} comparisons fail", flush=True)
        code = code or (0 if fired else 1)
    sys.exit(code)


if __name__ == '__main__':
    main()
