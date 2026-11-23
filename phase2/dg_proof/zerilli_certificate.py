"""#148: a certificate for the forward direction of the Zerilli reduction, every l, without dividing by M or s.

Statement checked: Z = psi,xx - (V + s^2) psi, with psi from Dotti-Gleiser gr-qc/0604021 eq. (9) and V from the issue, equals
sum_{i,k} c_ik D^k E_i exactly, D = d/dr on jets, k <= 2, c_ik in Q(r, M, s, L), L = l(l+1) one symbol, E_i the seven radial equations of the
even sector (H0, H1, H2, K kept). The denominators of the c_ik are factored and classified. The planted case is the odd sector, Psi = f h1 / r, V_RW.

Equations: #143's library zerilli_derivation_lib.py, imported unchanged (dual-number Ricci, Legendre reduction, coefficient extraction).
Order: planted (exit 2 on failure) -> odd certificate -> even certificate. Every check prints PASS id or FAIL id; exit 0, 1 (a real check
failed) or 2 (a planted check failed).

Run: python zerilli_certificate.py [--part=planted|odd|even|all] [--L=6] [--order=2] [--attempts=N] [--dump=FILE] [--mutate=NAME[,NAME]] [--failfast]
  --L=6        pilot: L = 6 substituted before the solve; the checks that need L symbolic (S.explicit, E.L, E.M0, E.S0) are skipped
  --attempts=N the solver tries at most N column orders (default 23: O1, reversed within k, k descending, seeds 1..20) until one passes E.den
  --dump=F     writes the even certificate (str of each coefficient) as JSON
Mutations (each makes one check fail; names, targets and exit codes are in negative-mass-raw/I148-rpi5/README.md):
  leg_sign l1 den_allow_M odd_v odd_psi odd_cert odd_den expl_norm v_coef psi_coef v_s2 drop_<eq> flow_target cert_coef l_spec m0_v s0_s1
  (eq in traceless tx tr rx tt rr trace)
"""
import json
import os
import random
import resource
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zerilli_derivation_lib import *
from sympy import QQ, Derivative, Integer
from sympy.polys.fields import field

PART, LVAL, ORDER, DUMP, ATTEMPTS = 'all', None, 2, None, 23
for _a in sys.argv[1:]:
    if _a.startswith('--part='):
        PART = _a[len('--part='):]
    elif _a.startswith('--L='):
        LVAL = int(_a[len('--L='):])
    elif _a.startswith('--order='):
        ORDER = int(_a[len('--order='):])
    elif _a.startswith('--attempts='):
        ATTEMPTS = int(_a[len('--attempts='):])
    elif _a.startswith('--dump='):
        DUMP = _a[len('--dump='):]

run = Run('zerilli_certificate: certificate of the Zerilli forward direction (part=%s, L=%s, order=%d)' % (PART, 'symbolic' if LVAL is None else LVAL, ORDER))
T0 = time.time()


def lap(msg):
    print('# [%6.1f s] %s' % (time.time() - T0, msg), flush=True)


Lp = L if LVAL is None else Integer(LVAL)           # the L of the solve
f_ = 1 - 2 * M / r
EXP = sp.exp(s * t)
Px = sp.Derivative(P, x)


# ======================================================================= fields, linear forms on jets, the solver
GENS = (r, M, s) if LVAL is not None else (r, M, s, L)
FF, *_g = field(GENS, QQ)
rf = _g[0]


def spec_field(gens):
    return field(gens, QQ)[0]


def to_ff(FFx, e):
    return FFx.from_expr(sp.sympify(e))


def is_jet(sy, names):
    nm, _, n = sy.name.rpartition('_')
    return nm in names and n.isdigit()


def jkey(sy):
    nm, _, n = sy.name.rpartition('_')
    return (nm, int(n))


def to_funcs(e, names):
    """jet symbols X_n -> X(r) and its r-derivatives."""
    rep = {}
    for sy in e.free_symbols:
        if is_jet(sy, names):
            nm, n = jkey(sy)
            F = Function(nm)(r)
            rep[sy] = F if n == 0 else Derivative(F, (r, n))
    return e.subs(rep)


def form_from_expr(e, names, FFx):
    """expression linear in the jet symbols -> {(name, n): coefficient in FFx}; linearity is checked."""
    out, back = {}, 0
    for sy in e.free_symbols:
        if is_jet(sy, names):
            c = cancel(diff(e, sy))
            out[jkey(sy)] = to_ff(FFx, c)
            back += c * sy
    assert cancel(e - back) == 0, 'form is not linear in the jets'
    return {k: v for k, v in out.items() if v}


def Dform(fm):
    """d/dr of a linear form on jets: (c F_n)' = c' F_n + c F_(n+1)."""
    out = {}
    for (fn, n), c in fm.items():
        dc = c.diff(rf)
        if dc:
            out[(fn, n)] = out.get((fn, n), FF.zero) + dc
        out[(fn, n + 1)] = out.get((fn, n + 1), FF.zero) + c
    return {k: v for k, v in out.items() if v}


def solve_span(cols, Dk, Zf):
    """Gauss-Jordan elimination in the field: unknowns c_col with sum_col c_col Dk[col] = Zf, one row per jet. Returns None if inconsistent,
    else ({col: c}, rank, rows, unknowns) with the free unknowns 0."""
    jets = sorted({j for c in cols for j in Dk[c]} | set(Zf))
    nr, nc = len(jets), len(cols)
    A = [[Dk[c].get(j, FF.zero) for c in cols] + [Zf.get(j, FF.zero)] for j in jets]
    piv, row = [], 0
    for col in range(nc + 1):
        if row == nr:
            break
        cand = [i for i in range(row, nr) if A[i][col]]
        if not cand:
            continue
        i0 = min(cand, key=lambda i: len(A[i][col].numer) + len(A[i][col].denom))
        A[row], A[i0] = A[i0], A[row]
        inv = FF.one / A[row][col]
        A[row] = [a * inv for a in A[row]]
        for i in range(nr):
            if i != row and A[i][col]:
                fct = A[i][col]
                A[i] = [a - fct * b for a, b in zip(A[i], A[row])]
        piv.append(col)
        row += 1
    if nc in piv:
        return None
    return {cols[pc]: A[i][nc] for i, pc in enumerate(piv) if A[i][nc]}, len(piv), nr, nc


def spec(c, subs):
    """the coefficient c (a field element) with `subs` substituted, as a sympy expression; None if its denominator vanishes identically."""
    n_, d_ = c.numer.as_expr().subs(subs), c.denom.as_expr().subs(subs)
    if sp.expand(d_) == 0:
        return None
    return n_ / d_


def expr_terms(C, Eex, names):
    """For the sum_{(i,k)} c_ik D^k E_i: a list of (c_ik as a sympy expression, {jet: coefficient of D^k E_i}). D^k is sympy's diff on the
    equation written with Function(r) (not the solver's jet derivative)."""
    EF = {nm: to_funcs(e, names) for nm, e in Eex.items()}
    out = []
    for (nm, k), c in sorted(C.items()):
        d_ = jetify(EF[nm] if k == 0 else sp.diff(EF[nm], r, k), names)
        out.append((c, {jkey(sy): diff(d_, sy) for sy in d_.free_symbols if is_jet(sy, names)}))
    return out


def residual_field(C, Eex, Zfun, names, FFx):
    """Z - sum c D^k E as {jet: element of FFx}; the identity holds iff it is empty. Exact canonical forms."""
    res = {}
    Zj = jetify(Zfun, names)
    for sy in Zj.free_symbols:
        if is_jet(sy, names):
            res[jkey(sy)] = to_ff(FFx, diff(Zj, sy))
    for c, jc in expr_terms(C, Eex, names):
        cF = to_ff(FFx, c)
        for j, ce in jc.items():
            res[j] = res.get(j, FFx.zero) - cF * to_ff(FFx, ce)
    return {j: v for j, v in res.items() if v}


def residual_points(C, Eex, Zfun, names, points):
    """the same residual evaluated at exact rational points by sympy substitution (no field code); returns the nonzero entries."""
    Zj = jetify(Zfun, names)
    zc = {jkey(sy): diff(Zj, sy) for sy in Zj.free_symbols if is_jet(sy, names)}
    terms = expr_terms(C, Eex, names)
    bad = []
    for pt in points:
        tot = {j: ce.subs(pt) for j, ce in zc.items()}
        for c, jc in terms:
            cv = c.subs(pt)
            for j, ce in jc.items():
                tot[j] = tot.get(j, 0) - cv * ce.subs(pt)
        bad += [(j, v) for j, v in tot.items() if v != 0]
    return bad


POINTS = [{r: sp.Rational(13, 7), M: sp.Rational(-3, 4), s: sp.Rational(11, 6), L: sp.Rational(37, 5)},
          {r: sp.Rational(5, 3), M: sp.Rational(2, 9), s: sp.Rational(-7, 5), L: sp.Rational(23, 4)},
          {r: sp.Rational(29, 11), M: sp.Rational(-1, 6), s: sp.Rational(3, 8), L: sp.Rational(61, 7)}]


# ======================================================================= denominators
def is_const(e):
    return not cancel(e).free_symbols


def pure_L_ok(g):
    """g: a polynomial in L alone. Allowed if it is not L - 2 and has no rational root l(l+1) with an integer l >= 2."""
    if is_const(g / (L - 2)):
        return False, 'factor L - 2'
    for rho in sp.roots(sp.Poly(g, L), filter='Q'):
        if rho.is_integer and rho >= 6:
            q, exact = sp.integer_nthroot(int(1 + 4 * rho), 2)
            if exact and q % 2 == 1:
                return False, 'root L = %s (l = %d)' % (rho, (q - 1) // 2)
    return True, ''


def classify_den(den_expr, rule, Lv=L):
    """rule 'even': allowed factors r, (L-2) r + 6M, polynomials in L alone; rule 'odd': r, r - 2M, polynomials in L alone.
    Returns (ok, [(factor, multiplicity, verdict)])."""
    allowed = [r, (Lv - 2) * r + 6 * M] if rule == 'even' else [r, r - 2 * M]
    if mut('den_allow_M'):
        allowed = allowed + [M]
    out, ok = [], True
    for g, m in sp.factor_list(den_expr)[1]:
        verdict = None
        for a in allowed:
            if is_const(g / a):
                verdict = 'allowed'
        if verdict is None and g.free_symbols <= {L}:
            okL, why = pure_L_ok(g)
            verdict = 'allowed (L only)' if okL else 'FORBIDDEN (%s)' % why
        if verdict is None:
            verdict = 'FORBIDDEN'
        ok = ok and not verdict.startswith('FORBIDDEN')
        out.append((g, m, verdict))
    return ok, out


def den_report(C, rule, label):
    """prints each coefficient's denominator and a summary of the factors with their highest exponents. Returns (ok, summary)."""
    ok_all, summ = True, {}
    for (nm, k), c in sorted(C.items()):
        ok, fl = classify_den(c.denom.as_expr(), rule, Lp)
        ok_all = ok_all and ok
        for g, m, v in fl:
            e = summ.setdefault(str(g), [0, v])
            e[0] = max(e[0], m)
        run.note('%s c[%s, k=%d]: denominator %s (%d terms), numerator %d terms' %
                 (label, nm, k, ' '.join('(%s)^%d' % (g, m) for g, m, _ in fl) or '1', len(c.denom), len(c.numer)))
    run.note('%s denominator factors, highest exponent over all coefficients: %s' %
             (label, '; '.join('(%s)^%d %s' % (g, m, v) for g, (m, v) in sorted(summ.items())) or 'none'))
    return ok_all, summ


# ======================================================================= planted checks
g0 = sp.diag(-f_, 1 / f_, r ** 2 / (1 - x ** 2), r ** 2 * (1 - x ** 2))
R0_bg, _ = lin_ricci(g0, sp.zeros(N, N))
run.planted('PL.0c', all(tidy(R0_bg[m][n]) == 0 for m in range(N) for n in range(N)), 'Ricci of g_Schw(M) is not 0')
lap('PL.0c')
ok = True
for n_ in range(2, 7):
    for order in (2, 3, 4):
        red_ = leg_red(sp.Derivative(P, (x, order)))
        P0s, P1s = sp.symbols('P0 P1')
        lhs = red_.subs({P0s: sp.legendre(n_, x), P1s: diff(sp.legendre(n_, x), x), L: n_ * (n_ + 1)})
        ok = ok and cancel(lhs - diff(sp.legendre(n_, x), x, order)) == 0
run.planted('PL.0d', ok, 'Legendre reduction differs from the explicit derivative of P_l')
lap('PL.0d')

# PL.D.1: the classifier on 11 constructed denominators (three pass and six fail under the even rule; one passes and one fails under the odd rule)
CASES = [('even', r ** 2 * ((L - 2) * r + 6 * M) ** 2, True), ('even', r * L * (L + 1), True),
         ('even', r ** 3 * ((L - 2) * r + 6 * M) * (L - 5), True),
         ('even', r * (r - 2 * M), False), ('even', s * r, False), ('even', M * r, False), ('even', (L - 2) * r, False),
         ('even', (L - 6) * r, False), ('even', (L - 12) * r, False),
         ('odd', r ** 2 * (r - 2 * M), True), ('odd', s * r, False)]
got = [classify_den(d_, rule)[0] for rule, d_, _ in CASES]
want = [w for _, _, w in CASES]
run.planted('PL.D.1', got == want, 'classifier verdicts %s, expected %s' % (got, want))
lap('PL.D.1')
if PART == 'planted':
    run.finish()

# ======================================================================= the equations
h0f, h1f, H0f, H1f, H2f, Kf = (Function(n)(r) for n in ('h0', 'h1', 'H0', 'H1', 'H2', 'K'))


def zeros():
    return [[0] * N for _ in range(N)]


def sym_set(h, m, n, v):
    h[m][n] = v
    h[n][m] = v


def h_even_for(Pfun):
    h = zeros()
    sym_set(h, 0, 0, f_ * H0f * Pfun * EXP)
    sym_set(h, 0, 1, H1f * Pfun * EXP)
    sym_set(h, 1, 1, H2f * Pfun * EXP / f_)
    sym_set(h, 2, 2, r ** 2 * Kf * Pfun * EXP / (1 - x ** 2))
    sym_set(h, 3, 3, r ** 2 * (1 - x ** 2) * Kf * Pfun * EXP)
    return h


h_odd = zeros()
sym_set(h_odd, 0, 3, h0f * EXP * (1 - x ** 2) * Px)
sym_set(h_odd, 1, 3, h1f * EXP * (1 - x ** 2) * Px)


def red(e):
    return tidy(leg_red(sp.expand(e / EXP)))


x2 = 1 - x ** 2


def certify(Eex, names, Zfun, rule, label, order_names):
    """Solves for the certificate of Zfun in the span of D^k Eex[i], k <= ORDER. Column order O1 first; if its denominators fail the rule, the other
    orders (README) up to ATTEMPTS. Returns (kept, first): the first solution that passes the rule (or None) and the first solution (or None)."""
    forms = {nm: form_from_expr(e, names, FF) for nm, e in Eex.items()}
    Zf = form_from_expr(jetify(Zfun, names), names, FF)
    Dk = {}
    for nm, fm in forms.items():
        cur = fm
        for k in range(ORDER + 1):
            Dk[(nm, k)] = cur
            cur = Dform(cur)
    names_eq = [n_ for n_ in order_names if n_ in forms]
    cols_all = [(nm, k) for k in range(ORDER + 1) for nm in names_eq]
    attempts = [('O1', cols_all), ('reversed-within-k', [(nm, k) for k in range(ORDER + 1) for nm in reversed(names_eq)]),
                ('k-descending', [(nm, k) for k in reversed(range(ORDER + 1)) for nm in names_eq])]
    for seed in range(1, 21):
        cp = list(cols_all)
        random.Random(seed).shuffle(cp)
        attempts.append(('seed%d' % seed, cp))
    kept, first = None, None
    for tag, cols in attempts[:ATTEMPTS]:
        t1 = time.time()
        res = solve_span(cols, Dk, Zf)
        if res is None:
            run.note('%s attempt %s: no solution (%.1f s)' % (label, tag, time.time() - t1))
            break                                       # the span does not depend on the column order
        C, rank, nr, nc = res
        okd = all(classify_den(c.denom.as_expr(), rule, Lp)[0] for c in C.values())
        run.note('%s attempt %s: solution found (rank %d, %d jets, %d unknowns, %.1f s), denominators %s' %
                 (label, tag, rank, nr, nc, time.time() - t1, 'allowed' if okd else 'NOT allowed'))
        first = first or C
        if okd:
            kept = C
            break
    return kept, first


# ======================================================================= odd parity (planted)
_, R1_odd = lin_ricci(g0, sp.Matrix(h_odd))
Rr_odd = [[red(R1_odd[m][n]) if n >= m else None for n in range(N)] for m in range(N)]
ONAMES = ['h0', 'h1']
ok_o1, oeq, why_o1 = all(Rr_odd[m][n] == 0 for m, n in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2), (3, 3)]), {}, []
for nm, comp, norm, basis in (('tphi', (0, 3), 1 / x2, ['P1']), ('rphi', (1, 3), 1 / x2, ['P1']), ('xphi', (2, 3), 1, ['P0', 'xP1'])):
    d_, ok_, why_ = coef_eqs(Rr_odd[comp[0]][comp[1]], norm, basis, nm)
    ok_o1 = ok_o1 and ok_
    if not ok_:
        why_o1.append('%s: %s' % (nm, why_))
    oeq[nm] = {k: jetify(v, ONAMES) for k, v in d_.items()}
ok_o1 = ok_o1 and 'P0' in oeq['xphi'] and cancel(oeq['xphi']['P0'] + L / 2 * oeq['xphi']['xP1']) == 0
run.planted('PL.O.1', ok_o1, 'odd components do not have the expected angular form: %s' % why_o1)
lap('PL.O.1 equations')
Eodd = {'tphi': oeq['tphi']['P1'], 'rphi': oeq['rphi']['P1'], 'xphi': oeq['xphi']['xP1']}
run.note('odd radial equations (jets h_n):')
for nm_, e_ in Eodd.items():
    run.note('  %s: %s' % (nm_, sp.factor_terms(e_)))
Eodd_p = {k: v.subs(L, Lp) for k, v in Eodd.items()}


def odd_target():
    """Z_odd = Psi,xx - (V_RW + s^2) Psi as an expression in Function(r); mutation odd_den: the xphi equation scaled to leading coefficient 1 in h0."""
    if mut('odd_den'):
        e = Eodd_p['xphi']
        return to_funcs(cancel(e / diff(e, jet_sym('h0', 0))), ONAMES)
    Psi = (1 if mut('odd_psi') else f_) * h1f / r
    VRW = f_ * (Lp / r ** 2 - (5 if mut('odd_v') else 6) * M / r ** 3)
    return f_ * diff(f_ * diff(Psi, r), r) - (VRW + s ** 2) * Psi


def check_identity(C, Eex, names, Zfun, FFx, label, cid, planted, points=True):
    """exact residual in FFx (sympy-diff path) and at three exact rational points; reported under cid."""
    Cexpr = {k_: (c if not hasattr(c, 'as_expr') else c.as_expr()) for k_, c in C.items() if c}
    bad = residual_field(Cexpr, Eex, Zfun, names, FFx)
    pts = residual_points(Cexpr, Eex, Zfun, names, [{k_: v_ for k_, v_ in pt.items() if k_ in set(GENS)} for pt in POINTS]) if points else []
    chk = run.planted if planted else run.check
    chk(cid, not bad and not pts, '%s: nonzero residual jets %s; nonzero at rational points %s' % (label, sorted(bad)[:6], pts[:3]))
    return not bad and not pts


if PART in ('odd', 'all'):
    kept, first = certify(Eodd_p, ONAMES, odd_target(), 'odd', 'odd', ['xphi', 'rphi', 'tphi'])
    run.planted('PL.O.2', first is not None, 'no certificate for the odd target with k <= %d: the target is not in the span' % ORDER)
    lap('PL.O.2')
    Codd = kept if kept is not None else first
    if mut('odd_cert'):
        Codd = dict(Codd)
        Codd[sorted(Codd)[0]] = Codd[sorted(Codd)[0]] + rf
    check_identity(Codd, Eodd_p, ONAMES, odd_target(), FF, 'odd certificate', 'PL.O.3', True)
    lap('PL.O.3')
    okd, summ = den_report(Codd, 'odd', 'odd')
    run.planted('PL.O.4', okd and kept is not None, 'odd denominators: forbidden factors %s (a passing attempt exists: %s)' %
                ({k_: v_ for k_, v_ in summ.items() if v_[1].startswith('FORB')}, kept is not None))
    lap('PL.O.4')
    for k_, c_ in sorted(Codd.items()):
        run.note('odd certificate c[%s, k=%d] = %s' % (k_[0], k_[1], sp.factor(c_.as_expr())))
    run.note('odd: Psi,xx = (V_RW + s^2) Psi is a combination of the three odd radial equations and their first %d r-derivatives' % ORDER)
    if PART == 'odd':
        run.finish()

# ======================================================================= even parity
if PART in ('even', 'all'):
    # ---- A.1 P and xP' are independent for L not in {0, 2}
    c_, Ls_ = sp.symbols('c_ Ls_')
    dep = sp.solve([c_ * (c_ - 1), Ls_ - c_ * (c_ + 1)], [c_, Ls_], dict=True)
    dep_L = sorted(set(d_[Ls_] for d_ in dep))
    Lq = sp.Integer(2) if mut('l1') else L
    ok_sym = dep_L == [0, 2] and all(cancel(Lq - d_) != 0 for d_ in dep_L)
    rank_ok = True
    for lv in range(2, 31):
        pl_ = sp.Poly(sp.legendre(lv, x), x)
        q_ = sp.Poly(sp.expand(x * diff(sp.legendre(lv, x), x)), x)
        rows = sp.Matrix([[pl_.coeff_monomial(x ** k) for k in range(lv + 1)], [q_.coeff_monomial(x ** k) for k in range(lv + 1)]])
        rank_ok = rank_ok and rows.rank() == 2
    run.check('A.1', ok_sym and rank_ok, 'P and xP\' can be dependent at L in %s; this run has L = %s; ranks l = 2..30 ok: %s' % (dep_L, Lq, rank_ok))
    lap('A.1')

    # ---- E.1 the seven even radial equations
    _, R1_ev = lin_ricci(g0, sp.Matrix(h_even_for(P)))
    Rr_ev = [[red(R1_ev[m][n]) if n >= m else None for n in range(N)] for m in range(N)]
    ENAMES = ['H0', 'H1', 'H2', 'K']
    SPEC = (('tt', (0, 0), 1, ['P0']), ('tr', (0, 1), 1, ['P0']), ('rr', (1, 1), 1, ['P0']), ('tx', (0, 2), 1, ['P1']), ('rx', (1, 2), 1, ['P1']),
            ('xx', (2, 2), x2, ['P0', 'xP1']), ('pp', (3, 3), 1 / x2, ['P0', 'xP1']))
    ok_e1, why_e1, ev = all(Rr_ev[m][n] == 0 for m, n in [(0, 3), (1, 3), (2, 3)]), [], {}
    for nm, comp, norm, basis in SPEC:
        d_, ok_, why_ = coef_eqs(Rr_ev[comp[0]][comp[1]], norm, basis, nm)
        ok_e1 = ok_e1 and ok_
        if not ok_:
            why_e1.append('%s: %s' % (nm, why_))
        for k_, v_ in d_.items():
            ev[nm + '.' + k_] = jetify(v_, ENAMES)
    id1 = cancel(ev['xx.xP1'] + ev['pp.xP1'])
    id2 = cancel(ev['pp.P0'] - ev['xx.P0'] - L * ev['xx.xP1'])
    trl = cancel(ev['xx.xP1'] - (jet_sym('H0', 0) - jet_sym('H2', 0)) / 2)
    run.check('E.1', ok_e1 and id1 == 0 and id2 == 0 and trl == 0, 'even forms: %s; identities (%s, %s); traceless - (H0-H2)/2 = %s' % (why_e1, id1, id2, trl))
    lap('E.1 equations')
    Eev = {'tt': ev['tt.P0'], 'tr': ev['tr.P0'], 'rr': ev['rr.P0'], 'tx': ev['tx.P1'], 'rx': ev['rx.P1'],
           'trace': cancel(ev['xx.P0'] + ev['pp.P0']), 'traceless': ev['xx.xP1']}
    run.note('even radial equations (jets F_n = d^n F / dr^n):')
    for nm_, e_ in Eev.items():
        run.note('  %s: %s' % (nm_, sp.factor_terms(e_)))
    for nm_ in ('tt', 'tr', 'rr', 'tx', 'rx', 'trace', 'traceless'):
        if mut('drop_' + nm_):
            del Eev[nm_]
            run.note('mutation: equation %s removed from the span' % nm_)
    Eev_p = {k: v.subs(L, Lp) for k, v in Eev.items()}

    # ---- the targets, typed from the issue and DG06 (9)
    def psi9(lam_, Mm, Kx, H2x):
        den2 = lam_ if mut('psi_coef') else lam_ + 1
        return r * (r - 2 * Mm) / ((lam_ + 1) * (lam_ * r + 3 * Mm)) * (H2x - r * diff(Kx, r)) + r * Kx / den2

    def V_issue(lam_, Mm):
        nine = 8 if mut('v_coef') else 9
        return 2 * (1 - 2 * Mm / r) * (lam_ ** 2 * (lam_ + 1) * r ** 3 + 3 * lam_ ** 2 * Mm * r ** 2 + nine * lam_ * Mm ** 2 * r + 9 * Mm ** 3) / \
            (r ** 3 * (lam_ * r + 3 * Mm) ** 2)

    def Z_even(lam_, Mm, ssq):
        """psi,xx - (V + ssq) psi, d/dx = f d/dr; ssq = s^2, or 0 for the static target."""
        ff = 1 - 2 * Mm / r
        psi = psi9(lam_, Mm, Kf, H2f)
        return ff * diff(ff * diff(psi, r), r) - (V_issue(lam_, Mm) + ssq) * psi

    lam_p = (Lp - 2) / 2
    s2 = 0 if mut('v_s2') else s ** 2
    if mut('flow_target'):
        Zev = to_funcs(jet_sym('K', 1) + (r - 3 * M) / (r * (r - 2 * M)) * jet_sym('K', 0) - Lp * jet_sym('H1', 0) / (2 * s * r ** 2) - jet_sym('H2', 0) / r, ENAMES)
    else:
        Zev = Z_even(lam_p, M, s2)
    kept, first = certify(Eev_p, ENAMES, Zev, 'even', 'even', ['traceless', 'tx', 'tr', 'rx', 'tt', 'rr', 'trace'])
    lap('solve')
    run.check('E.cert', first is not None, 'the even target is not in the span of the equations and their first %d r-derivatives (no solution over Q(r, M, s%s))' %
              (ORDER, ', L' if LVAL is None else ''))
    if first is None:
        run.note('skipped E.id, E.den, E.L, E.M0, E.S0: no certificate')
        run.finish()
    Ce = kept if kept is not None else first
    if mut('cert_coef'):
        Ce = dict(Ce)
        Ce[sorted(Ce)[0]] = Ce[sorted(Ce)[0]] + rf
    check_identity(Ce, Eev_p, ENAMES, Zev, FF, 'even certificate', 'E.id', False)
    lap('E.id')
    okd, summ = den_report(Ce, 'even', 'even')
    run.check('E.den', okd and kept is not None, 'denominators: forbidden factors %s (a passing attempt exists: %s)' %
              ({k_: v_ for k_, v_ in summ.items() if v_[1].startswith('FORB')}, kept is not None))
    lap('E.den')
    run.note('even certificate: %d nonzero coefficients over %d unknowns' % (len(Ce), len(Eev) * (ORDER + 1)))
    for k_, c_ in sorted(Ce.items()):
        run.note('even certificate c[%s, k=%d] = %s' % (k_[0], k_[1], sp.factor(c_.as_expr())))
    if DUMP:
        json.dump({'%s.k%d' % k_: str(c.as_expr()) for k_, c in Ce.items()}, open(DUMP, 'w'), indent=0)

    if LVAL is not None:
        run.note('pilot (L = %d): S.explicit, E.L, E.M0, E.S0 skipped' % LVAL)
        run.note('peak rss %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
        run.finish()

    # ---- S.explicit: the coefficient forms against delta R computed with the explicit P_l (no Legendre reduction, no coefficient extraction)
    ok_ex, why_ex = True, []
    for lv in (2, 3, 4):
        Pl = sp.legendre(lv, x)
        _, R1x = lin_ricci(g0, sp.Matrix(h_even_for(Pl)))
        evl = {k_: v_.subs(L, lv * (lv + 1)) for k_, v_ in ev.items()}
        Pd = diff(Pl, x)
        rec = {(0, 0): evl['tt.P0'] * Pl, (0, 1): evl['tr.P0'] * Pl, (1, 1): evl['rr.P0'] * Pl, (0, 2): evl['tx.P1'] * Pd, (1, 2): evl['rx.P1'] * Pd,
               (2, 2): (evl['xx.P0'] * Pl + evl['xx.xP1'] * x * Pd) / (1 if mut('expl_norm') else x2),
               (3, 3): (evl['pp.P0'] * Pl + evl['pp.xP1'] * x * Pd) * x2,
               (0, 3): 0, (1, 3): 0, (2, 3): 0}
        for (m, n), want in rec.items():
            have = jetify(sp.expand(R1x[m][n] / EXP), ENAMES)
            if cancel(together(have - want)) != 0:
                ok_ex = False
                why_ex.append('l=%d component %s' % (lv, (m, n)))
        lap('S.explicit l=%d' % lv)
    run.check('S.explicit', ok_ex, 'explicit-P_l delta R differs from the coefficient forms at: %s' % why_ex[:6])

    # ---- E.L: L = 6, 12, 20; Z retyped from lambda = (l-1)(l+2)/2, the certificate and the equations at L = l(l+1)
    ok_L, why_L = True, []
    for lv in (2, 3, 4):
        Lv = lv * (lv + 1)
        FFl = spec_field((r, M, s))
        lam_v = sp.Rational((lv - 1) * (lv + (1 if mut('l_spec') else 2)), 2)
        Cl, vanish = {}, []
        for k_, c in Ce.items():
            cv = spec(c, {L: Lv})
            if cv is None:
                vanish.append(k_)
            else:
                Cl[k_] = cv
        bad = residual_field(Cl, {k_: v_.subs(L, Lv) for k_, v_ in Eev.items()}, Z_even(lam_v, M, s2), ENAMES, FFl)
        if vanish or bad:
            ok_L = False
            why_L.append('l=%d: denominators vanishing %s, residual jets %s' % (lv, vanish, sorted(bad)[:4]))
        run.note('E.L l=%d (L = %d): common denominator of the specialized coefficients %s' %
                 (lv, Lv, sp.factor(sp.lcm([c.denom.as_expr().subs(L, Lv) for c in Ce.values()]))))
        lap('E.L l=%d' % lv)
    run.check('E.L', ok_L, 'identity after L = 6, 12, 20 fails: %s' % why_L)

    # ---- E.M0: M = 0, V = L/r^2
    lam_s = (L - 2) / 2
    V_diff = cancel(V_issue(lam_s, 0) - ((L + 2) if mut('m0_v') else L) / r ** 2)
    Cm, vm = {}, []
    for k_, c in Ce.items():
        cv = spec(c, {M: 0})
        if cv is None:
            vm.append(k_)
        else:
            Cm[k_] = cv
    psi0 = psi9(lam_s, 0, Kf, H2f)
    Zm0 = diff(psi0, r, 2) - (((L + 2) if mut('m0_v') else L) / r ** 2 + s ** 2) * psi0       # f = 1 at M = 0
    bad_m = residual_field(Cm, {k_: v_.subs(M, 0) for k_, v_ in Eev.items()}, Zm0, ENAMES, spec_field((r, s, L)))
    run.check('E.M0', V_diff == 0 and not vm and not bad_m,
              'M = 0: V(issue) - V_typed = %s; vanishing denominators %s; residual jets %s' % (V_diff, vm, sorted(bad_m)[:4]))
    lap('E.M0')

    # ---- E.S0: s = 0, the static equation psi,xx = V psi
    Cs, vs = {}, []
    for k_, c in Ce.items():
        cv = spec(c, {s: 0})
        if cv is None:
            vs.append(k_)
        else:
            Cs[k_] = cv
    bad_s = residual_field(Cs, {k_: v_.subs(s, 0) for k_, v_ in Eev.items()}, Z_even(lam_s, M, 1 if mut('s0_s1') else 0), ENAMES, spec_field((r, M, L)))
    run.check('E.S0', not vs and not bad_s, 's = 0: vanishing denominators %s; residual jets %s' % (vs, sorted(bad_s)[:4]))
    lap('E.S0')
    run.note('peak rss %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))

run.finish()
