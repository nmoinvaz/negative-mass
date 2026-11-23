#!/usr/bin/env python3
"""Step 7 of #136: compares the results of review_he_type3_count.py with #133's files.

Opens #133's code and files on purpose, after the review's own runs are done; nothing here feeds back into
review_he_type3_count.py. Compares
  tables : my E1, E2, P (frame table of he_type3.md section 1) with he_type3.shadow() (the metric route of #102), up to scale;
  roots  : my refined root enclosures of P and y with the root leaves of #133's cells.tsv, and the tiling of its leaves;
  ends   : my lower bounds of G on (0, 1/4] and of S on [z_t, oo) with #133's, and both with the exact limits.

Usage (from the repository root):
  python3 notes/warp/review_he_type3_count_compare.py --cells=PATH/I133-rpi5/cells.tsv [--mutate=NAME]
Mutations: compare-drop-term (a term of my P table is dropped), compare-shift-leaf (one root leaf of #133 is moved).
"""
import os
import re
import sys
import time
import argparse
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import review_he_type3_count as R   # noqa: E402
import he_type3 as H3               # noqa: E402

MUT = None
FAILS = []
WALLS = ['le', 'rodal']


def say(ok, msg):
    print('[%s] %s' % ('PASS' if ok else 'FAIL', msg), flush=True)
    if not ok:
        FAILS.append(msg.split(':')[0])


def ratio_equal(mine, theirs):
    """mine, theirs: {key: Fraction}. True with the common ratio if they agree term by term up to one scale."""
    if set(mine) != set(theirs):
        return False, None, 'term sets differ (%d against %d terms)' % (len(mine), len(theirs))
    ratios = {theirs[k] / mine[k] for k in mine}
    if len(ratios) != 1:
        return False, None, '%d different ratios' % len(ratios)
    return True, ratios.pop(), ''


def check_tables():
    al = R.algebra()
    sh = H3.shadow()
    sp = al['sp']
    mine_P = dict(al['Ptab'])
    if MUT == 'compare-drop-term':
        del mine_P[next(iter(mine_P))]
    theirs_P = {mn: Fr(int(cf.p), int(cf.q)) for mn, cf in sp.Poly(sh['P'], H3.ph, H3.Ap, H3.Bp).terms()}
    ok, ra, why = ratio_equal(mine_P, theirs_P)
    say(ok, 'tables: P of the frame table equals he_type3.shadow()["P"] term by term, %d terms, ratio %s %s' % (len(mine_P), ra, why))
    for name, key, tab in (('E1', 0, al['E1tab']), ('E2', 1, al['E2tab'])):
        theirs = {}
        for m, co in enumerate(sh['eexact'][key]):
            for mn, p, q in co:
                theirs[(m,) + tuple(mn)] = Fr(p, q)
        ok, ra, why = ratio_equal(dict(tab), theirs)
        say(ok, 'tables: %s of the frame table equals he_type3.shadow() term by term, %d terms, ratio %s %s' % (name, len(tab), ra, why))


def check_statements():
    """Statements of he_type3_count.md that my code can test: weights, degrees and the series of Lemma 3."""
    import mpmath as mp
    mp.mp.dps = 60
    gt = R.G_table()
    wts = sorted(set(gt['wts'].values()))
    say(wts == list(range(24, 67, 2)) and len(gt['tab']) == 358,
        'statements: P-hat in (phi, psi, delta) has %d terms of weights %s..%s step 2 (he_type3_count.md Lemma 3: 358 terms, 24, 26, ..., 66)'
        % (len(gt['tab']), wts[0], wts[-1]))
    jk = sorted({m[1] + m[2] for m in R.algebra()['Ptab']})
    say(jk == list(range(9, 16)), 'statements: (A\', B\')-degrees of the terms of P are %s (he_type3_count.md Lemma 1: 9 to 15)' % jk)
    for wn in WALLS:
        W = R.walls()[wn]
        C, v = mp.cosh(R.mpf_of(2 * W.sigma * W.Rb)), R.mpf_of(W.v)
        bad = 0
        for z in (Fr(1, 5), Fr(1, 4)):
            zz = R.mpf_of(z)
            D = C + mp.cosh(zz)
            N1 = sum((mp.mpf(4) ** k + 2 * C - 2 - 2 * k * (1 + C)) / mp.factorial(2 * k) * zz ** (2 * k - 4) for k in range(2, 60))
            N2 = sum((mp.mpf(4) ** k * (2 * k + 3) / 2 - 2 * k * C) / mp.factorial(2 * k + 1) * zz ** (2 * k - 2) for k in range(1, 60))
            S2 = sum(zz ** (2 * k - 2) / mp.factorial(2 * k) for k in range(1, 60))
            theirs = (v * S2 / D, v * N1 / D ** 2, v * (1 + C) * N2 / D ** 3)       # Phi, Psi, Delta of he_type3_count.md Lemma 3
            _, Phi, Psi, Delta = R.small_z_forms(z, z, W)
            for e, t in zip((Phi, Psi, Delta), theirs):
                bad += not R._mp_enc_contains(e, t, 40)
        say(bad == 0, 'statements: %s Phi, Psi, Delta from the series S2, N1, N2 of he_type3_count.md Lemma 3 (60 terms) lie in my enclosures at z = 1/5, 1/4: %d misses' % (wn, bad))


def parse_cells(path):
    out = {}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f, wall, a, b, kind = line.rstrip('\n').split('\t')
        out.setdefault((f, wall), []).append((Fr(a), Fr(b), kind))
    for k in out:
        out[k].sort()
    return out


def check_roots(cells_path):
    cells = parse_cells(cells_path)
    first = True
    for wn, wall_name in (('le', 'Le'), ('rodal', 'Rodal')):
        if wn not in WALLS:
            continue
        W = R.walls()[wn]
        cv = R.get_cover(wn)
        for mine_key, theirs_key in (('P', 'g'), ('y', 'Yn')):
            leaves = cells[(theirs_key, wall_name)]
            roots = [lf for lf in leaves if lf[2] == 'r']
            if MUT == 'compare-shift-leaf' and first:
                a, b, k = roots[0]
                roots[0] = (a + 1, b + 1, k)
                first = False
            tl = R.tiling_ok([(a, b) + (k,) for a, b, k in leaves], Fr(1, 4), Fr(W.zt))
            enc = cv[mine_key]['enc']
            inside = [sum(1 for (a, b, _) in roots if a <= e.lo and e.hi <= b) for e in enc]
            ok = tl and len(roots) == len(enc) and inside == [1] * len(enc)
            say(ok, 'roots: %s %s: #133 has %d leaves (mine %d) tiling [1/4, %d]: %s; its %d root leaves each hold exactly one of my %d refined '
                    'root enclosures (width at most %s): %s'
                % (wall_name, mine_key, len(leaves), len(cv[mine_key]['leaves']), W.zt, tl, len(roots), len(enc),
                   R.fl(max(e.width() for e in enc), 2), inside == [1] * len(enc)))
            widths = [R.fl(Fr(b - a), 2) for a, b, _ in roots]
            print('       #133 root leaf widths: %s' % widths)


def parse_out(path):
    txt = open(path).read()
    res = {}
    m = re.search(r'\[PASS\] small_z: Le: .*?G in \[([0-9.e+-]+), ([0-9.e+-]+)\] on the first cell.*?Rodal: .*?G in \[([0-9.e+-]+), ([0-9.e+-]+)\] on the first cell', txt)
    res['G'] = {'le': (Fr(m.group(1)), Fr(m.group(2))), 'rodal': (Fr(m.group(3)), Fr(m.group(4)))}
    res['L'] = {'le': Fr(re.search(r'tail_z: Le:.*?L >= ([0-9.]+)', txt).group(1)), 'rodal': Fr(re.search(r'Rodal: P z_t = 60.*?L >= ([0-9.]+)', txt).group(1))}
    return res


def check_ends():
    theirs = parse_out(os.path.join(HERE, 'he_type3_count_output.txt'))
    for wn in WALLS:
        W = R.walls()[wn]
        # small z: the exact value G(0) = 774144 Phi0^8 Delta0^2, Phi0 = v/(2(1+C)), Delta0 = -v (C - 5)/(3 (1+C)^2)
        cells, bad = R.small_z_cover(W, Fr(1, 4))
        gl = min(c[2].lo for c in cells)
        import mpmath as mp
        mp.mp.dps = 50
        v, C = R.mpf_of(W.v), mp.cosh(R.mpf_of(2 * W.sigma * W.Rb))
        g0 = 774144 * (v / (2 * (1 + C))) ** 8 * (v * (C - 5) / (3 * (1 + C) ** 2)) ** 2
        lo2, hi2 = theirs['G'][wn]
        ok = not bad and R.mpf_of(gl) <= g0 and R.mpf_of(lo2) * (1 - mp.mpf(10) ** -5) <= g0 <= R.mpf_of(hi2) * (1 + mp.mpf(10) ** -5)
        say(ok, 'ends: %s (0, 1/4]: my lower bound of G on the whole interval %s and #133\'s enclosure [%s, %s] on [0, 1/32] both bracket G(0) = %s'
            % (wn, R.fl(gl, 6), R.fl(lo2, 6), R.fl(hi2, 6), mp.nstr(g0, 6)))
        # large z
        tb = R.tail_bounds(W, W.zt)
        smin = min(c[2].lo for c in tb['cells'])
        lim = 4 * W.v ** 3
        L = theirs['L'][wn]
        ok = L <= smin <= lim and L <= lim
        say(ok, 'ends: %s [%d, oo): my lower bound of S = P/(t z)^9 is %s, #133\'s L is %s, the limit of S as z -> oo is 4 v^3 = %s: L <= mine <= limit %s'
            % (wn, W.zt, R.fl(smin, 6), R.fl(L, 6), R.fl(lim, 6), ok))


def dump_tables(path):
    """Writes the tables of P (i, j, k, coefficient) and of E1, E2 (w-power, i, j, k, coefficient) for the Lean work of #135."""
    al = R.algebra()
    with open(path, 'w') as f:
        f.write('# table\tw-power (E1, E2 only)\tphi-power\tA\'-power\tB\'-power\tcoefficient (exact)\n')
        for name, tab in (('P', al['Ptab']), ('E1', al['E1tab']), ('E2', al['E2tab'])):
            for key in sorted(tab):
                key = tuple(key)
                wp, ex = (key[0], key[1:]) if name != 'P' else ('', key)
                f.write('%s\t%s\t%d\t%d\t%d\t%s\n' % (name, wp, ex[0], ex[1], ex[2], tab[key]))
    print('tables written: P %d terms, E1 %d, E2 %d, to %s' % (len(al['Ptab']), len(al['E1tab']), len(al['E2tab']), os.path.basename(path)), flush=True)


def main():
    global MUT
    ap = argparse.ArgumentParser()
    ap.add_argument('--cells', required=True)
    ap.add_argument('--mutate', default=None)
    ap.add_argument('--dump', default=None, help='write the tables of P, E1, E2 to this file')
    ap.add_argument('--only', default='tables,statements,ends,roots', help='comparisons to run')
    ap.add_argument('--walls', default='le,rodal')
    args = ap.parse_args()
    MUT = args.mutate
    if args.dump:
        dump_tables(args.dump)
    t0 = time.time()
    WALLS[:] = args.walls.split(',')
    only = args.only.split(',')
    if 'tables' in only:
        check_tables()
    if 'statements' in only:
        check_statements()
    if 'ends' in only:
        check_ends()
    if 'roots' in only:
        check_roots(args.cells)
    print('summary: %d comparisons failed %s, %.0f s' % (len(FAILS), FAILS, time.time() - t0))
    sys.exit(1 if FAILS else 0)


if __name__ == '__main__':
    main()
