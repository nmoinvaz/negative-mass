"""#145: the static (s = 0) even-parity equation and its master variable from the linearized Einstein equations, M symbolic.

Metric g = g_Schw(M) + eps h, h in Regge-Wheeler gauge, axisymmetric, angular dependence P_l(cos theta) (l symbolic through L = l(l+1)),
NO time dependence. Library: zerilli_derivation_lib.py of #143, imported (dual-number Ricci, covariant cross-check, Legendre reduction, harness).
Order: planted checks (exit 2 on failure) -> static odd parity (planted target: Regge-Wheeler potential) -> static even parity
(H1 sector, flow, DG06 (9), V, converse) -> the s -> 0 relation to #143 -> divisions.
Sources: Dotti-Gleiser gr-qc/0604021 eqs. (2)-(4), (6), (9) (typed, see I145-rpi5/README.md); Zerilli's papers were not read.
Run: python zerilli_static.py [--part=planted|odd|h1|even|rel|all] [--mutate=NAME[,NAME]] [--failfast] [--skip=PL.G,PL.X] > zerilli_static_output.txt
Mutations (target in I145-rpi5/README.md): dual_inv ds_typo bg_nonvac leg_sign drop_trace route2_drop ang_norm_odd inj_rphi rw_typo fsign_ans cpm_typo
  cmp_blind l1 ang_norm s_leak inj_tr inj_tx inj_h1 ans_rr grr_drop tt_drop_H psi_zero mv_coef fsign_mv mv_lam v_typo res_typo dg3_typo inv_typo
  dep_typo tr_typo h1t_typo sr_typo div_set d2_zero div_s div_l2 d5_zero
"""
import os
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zerilli_derivation_lib import *
from sympy import QQ
from sympy.polys.fields import field

PART = 'all'
SKIP = []
for _a in sys.argv[1:]:
    if _a.startswith('--part='):
        PART = _a[len('--part='):]
    if _a.startswith('--skip='):
        SKIP = _a[len('--skip='):].split(',')
STAGES = ['planted', 'odd', 'h1', 'even', 'rel', 'all']


def skipped(cid):
    return any(cid.startswith(p) for p in SKIP)


run = Run('zerilli_static: static even-parity equation and master variable from the linearized Einstein equations (part=%s)' % PART)
T0 = time.time()


def lap(msg):
    print('# [%6.1f s] %s' % (time.time() - T0, msg), flush=True)


def stop_after(stage):
    if STAGES.index(PART) <= STAGES.index(stage):
        run.finish()


f = f_bg() + (M ** 2 / r ** 2 if mut('bg_nonvac') else 0)   # background metric's f (mutation bg_nonvac: not a vacuum solution)
fq = fM()                                                    # f = 1 - 2M/r for the ansatz, the master variables and the tortoise factor
g0 = sp.diag(-f, 1 / f, r ** 2 / (1 - x ** 2), r ** 2 * (1 - x ** 2))
Px = sp.Derivative(P, x)
x2 = 1 - x ** 2
E1 = sp.Integer(1)                                           # the static ansatz has no time dependence
Edep = sp.exp(s * t)                                         # used only by S.0 and S.R*


def red(e, E=E1):
    """first-order component -> divide e^{st} (E = 1: static), Legendre reduction, an expression in P0, P1 with rational coefficients."""
    return tidy(leg_red(sp.expand(e / E)))


def reduce_all(R1, E=E1):
    return [[red(R1[m][n], E) if n >= m else None for n in range(N)] for m in range(N)]


def same(e1, e2):
    d = cancel(together(e1 - e2))
    if mut('cmp_blind'):                                     # mutation cmp_blind: the comparison is evaluated at M = 0 only
        d = cancel(d.subs(M, 0))
    return d == 0


# ======================================================================= planted checks
hg = [[Function('h%d%d' % (min(i, j), max(i, j)))(t, r, x) for j in range(N)] for i in range(N)]
g_, gi_ = build_metric(g0, hg)
prod0 = sp.Matrix(N, N, lambda i, j: sum((g_[i][k] * gi_[k][j]).a0 for k in range(N)))
prod1 = sp.Matrix(N, N, lambda i, j: sum((g_[i][k] * gi_[k][j]).a1 for k in range(N)))
run.planted('PL.0a', prod0.applyfunc(tidy) == sp.eye(N) and prod1.applyfunc(tidy) == sp.zeros(N, N),
            'g g^-1 != 1 at order eps: the inverse used is not the first-order inverse')
lap('PL.0a')

a_ = Symbol('a', positive=True)
fds = 1 - r ** 2 / a_ ** 2
g_ds = sp.diag(-fds, 1 / fds, r ** 2 / (1 - x ** 2), r ** 2 * (1 - x ** 2))
R0_ds, _ = lin_ricci(g_ds, sp.zeros(N, N))
c_ds = 2 if mut('ds_typo') else 3                            # mutation ds_typo: the expected Ricci typed as (2/a^2) g
run.planted('PL.0b', all(tidy(R0_ds[m][n] - c_ds / a_ ** 2 * g_ds[m, n]) == 0 for m in range(N) for n in range(N)),
            'Ricci of static de Sitter is not (3/a^2) g')
lap('PL.0b')

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

# ----------------------------------------------------------------------- ansatz and gauge fields (all functions of r only)
r_funcs = {n: Function(n)(r) for n in ('H0', 'H1', 'H2', 'K', 'h0', 'h1', 'xa', 'xb', 'xc', 'xd')}
H0, H1, H2, K, h0, h1, xa, xb, xc, xd = (r_funcs[n] for n in ('H0', 'H1', 'H2', 'K', 'h0', 'h1', 'xa', 'xb', 'xc', 'xd'))


def zeros():
    return [[0] * N for _ in range(N)]


def sym_set(h, m, n, val):
    h[m][n] = val
    h[n][m] = val


def even_ansatz(E, gtr=1, rr_mut=False):
    h = zeros()
    sym_set(h, 0, 0, fq * H0 * P * E)
    sym_set(h, 0, 1, gtr * H1 * P * E)
    sym_set(h, 1, 1, H2 * P * E / (1 if rr_mut else fq))
    sym_set(h, 2, 2, r ** 2 * K * P * E / (1 - x ** 2))
    sym_set(h, 3, 3, r ** 2 * (1 - x ** 2) * K * P * E)
    return h


h_even = even_ansatz(E1, rr_mut=mut('ans_rr'))              # mutation ans_rr: g_rr = f^-1 + eps H2 P (no 1/f on H2)
h_odd = zeros()
sym_set(h_odd, 0, 3, h0 * (1 - x ** 2) * Px)
sym_set(h_odd, 1, 3, h1 * (1 - x ** 2) * Px)


def lie(g0m, xi):
    out = zeros()
    for m in range(N):
        for n in range(N):
            out[m][n] = sum(xi[a] * diff(g0m[m, n], COORD[a]) + g0m[a, n] * diff(xi[a], COORD[m]) + g0m[m, a] * diff(xi[a], COORD[n])
                            for a in range(N))
    return out


xi_even = [xa * P, xb * P, xc * (1 - x ** 2) * Px / r ** 2, 0]
xi_odd = [0, 0, 0, xd * Px]
h_gauge_even = lie(g0, xi_even)
h_gauge_odd = lie(g0, xi_odd)

_LIN = {}


def lin(h, key):
    """dual-number route: unreduced first-order Ricci of g0 + eps h (cached by key)."""
    if key not in _LIN:
        _, _LIN[key] = lin_ricci(g0, sp.Matrix(h))
    return _LIN[key]


if not skipped('PL.G'):
    bad = []
    for nm, hh in (('even', h_gauge_even), ('odd', h_gauge_odd)):
        Rr = reduce_all(lin(hh, 'gauge_' + nm))
        bad += [(nm, m, n) for m in range(N) for n in range(m, N) if Rr[m][n] != 0]
    run.planted('PL.G', not bad, 'pure-gauge h does not give dR = 0 at (parity, m, n) = %s' % bad)
    lap('PL.G')

if not skipped('PL.X'):
    bad = []
    for nm, hh in (('even', h_even), ('odd', h_odd), ('gauge_even', h_gauge_even)):
        R1 = lin(hh, nm if nm != 'gauge_even' else 'gauge_even')
        R1b = route2_delta_ricci(g0, hh)
        bad += [(nm, m, n) for m in range(N) for n in range(m, N) if red(R1[m][n] - R1b[m][n]) != 0]
    run.planted('PL.X', not bad, 'route 1 (dual) and route 2 (covariant) differ at (h, m, n) = %s' % bad)
    lap('PL.X')

stop_after('planted')


# ======================================================================= static odd parity (planted target: Regge-Wheeler)
def jets_of(e):
    """names of the jet symbols F_n (F a field of the ansatz) occurring in e"""
    return {sy.name for sy in sp.sympify(e).free_symbols if '_' in sy.name and sy.name.rsplit('_', 1)[1].isdigit()}


Rr_odd = reduce_all(lin(h_odd, 'odd'))
odd_names = ['h0', 'h1']
ok_o1 = all(Rr_odd[m][n] == 0 for m, n in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2), (3, 3)])
why_o1 = []
oeq = {}
NORM_T = 1 if mut('ang_norm_odd') else 1 / x2                # mutation ang_norm_odd: tphi read without 1/(1-x^2)
for nm, comp, norm, basis in (('tphi', (0, 3), NORM_T, ['P1']), ('rphi', (1, 3), 1 / x2, ['P1']), ('xphi', (2, 3), 1, ['P0', 'xP1'])):
    d_, ok_, why_ = coef_eqs(Rr_odd[comp[0]][comp[1]], norm, basis, nm)
    ok_o1 = ok_o1 and ok_
    if not ok_:
        why_o1.append('%s: %s' % (nm, why_))
    oeq[nm] = {k: jetify(v, odd_names) for k, v in d_.items()}
ok_o1 = ok_o1 and 'P0' in oeq['xphi'] and 'xP1' in oeq['xphi'] and cancel(oeq['xphi']['P0'] + L / 2 * oeq['xphi']['xP1']) == 0
run.planted('PL.O.1', ok_o1, 'static odd components do not have the expected angular form: %s' % why_o1)
Etphi, Erphi, Exphi = oeq['tphi']['P1'], oeq['rphi']['P1'], oeq['xphi']['xP1']
if mut('inj_rphi'):
    Erphi = Erphi + jet_sym('h0', 0)                         # mutation inj_rphi: an h0 term in the r-phi equation
run.note('static odd radial equations (jets h0_n, h1_n):')
for nm_, e_ in (('t-phi', Etphi), ('r-phi', Erphi), ('x-phi', Exphi)):
    run.note('  %s: %s' % (nm_, sp.factor(e_)))

h00, h01, h02 = (jet_sym('h0', n_) for n_ in range(3))
h10 = jet_sym('h1', 0)
c_r = cancel(diff(Erphi, h10)) if jets_of(Erphi) <= {'h1_0'} else 0
ok_r = jets_of(Erphi) <= {'h1_0'} and c_r != 0 and cancel(Erphi - c_r * h10) == 0
ok_x = all(not nm_.startswith('h0_') for nm_ in jets_of(Exphi)) and Exphi.subs({jet_sym('h1', n_): 0 for n_ in range(4)}) == 0
ok_t = all(not nm_.startswith('h1_') for nm_ in jets_of(Etphi)) and cancel(diff(Etphi, h02)) != 0 and jets_of(Etphi) <= {'h0_0', 'h0_1', 'h0_2'}
run.planted('PL.O.2', ok_r and ok_x and ok_t,
            'static odd reduction: r-phi is c h1 alone: %s (c = %s); x-phi has h1 jets only: %s; t-phi has h0 jets only and is second order: %s' %
            (ok_r, c_r, ok_x, ok_t))
run.note('PL.O.2: r-phi = (%s) h1, so h1 = 0 for L != 2; x-phi has h1 jets only, so it is implied; t-phi is the ODE of h0' % sp.factor(c_r))
c_h02 = cancel(diff(Etphi, h02))
F_odd = cancel(-(Etphi - c_h02 * h02) / c_h02)               # h0'' = F(h0, h0')
run.note('PL.O.2: h0\'\' = %s   (divisor: %s)' % (sp.factor(F_odd), sp.factor(c_h02)))
CPM = 3 if mut('cpm_typo') else 2                            # mutation cpm_typo: Psi = r h0' - 3 h0
psi_o = r * h01 - CPM * h00                                  # Psi = r h0' - 2 h0 (Cunningham-Price-Moncrief combination, static)


FF, rF, MF, sF, LF = field((r, M, s, L), QQ)                 # exact rational function field Q(r, M, s, L); s occurs in no static quantity


def toF(e_):
    return FF.from_expr(sp.sympify(e_))


def master2(states, flow, psi, ftort):
    """2D first-order system d states/dr = flow (linear forms), psi a linear form, d/dx = ftort d/dr. Returns the exact field data:
    (psi form, psi_x form, psi_xx form, det of (psi, psi_x) against the states, helpers). Rational function field Q(r, M, s, L)."""
    fl = [[toF(diff(flow[i], states[j])) for j in range(2)] for i in range(2)]
    ps = [toF(diff(psi, states[j])) for j in range(2)]
    ft = toF(ftort)

    def Dvec(v):
        return [v[j].diff(rF) + sum(v[i] * fl[i][j] for i in range(2)) for j in range(2)]
    psi_x = [ft * c for c in Dvec(ps)]
    psi_xx = [ft * c for c in Dvec(psi_x)]
    det = ps[0] * psi_x[1] - ps[1] * psi_x[0]
    return ps, psi_x, psi_xx, det, Dvec


def ab_from(ps, psi_x, psi_xx, det):
    a = (psi_xx[0] * psi_x[1] - psi_xx[1] * psi_x[0]) / det
    b = (-psi_xx[0] * ps[1] + psi_xx[1] * ps[0]) / det
    return cancel(a.as_expr()), cancel(b.as_expr())


ps_o, px_o, pxx_o, det_o, _ = master2([h00, h01], [h01, F_odd], psi_o, fq)
a_o, b_o = ab_from(ps_o, px_o, pxx_o, det_o)
V_RW = fq * (L / r ** 2 - (3 if mut('rw_typo') else 6) * M / r ** 3)
run.planted('PL.O.3', cancel(b_o) == 0 and same(a_o, V_RW),
            'static odd: Psi = r h0\' - 2 h0, coefficient of Psi_x = %s, V - V_RW = %s' % (b_o, cancel(a_o - V_RW)))
run.planted('PL.O.4', not same(V_RW, V_RW + M / r ** 3) and same(V_RW, sp.expand(V_RW)), 'comparison routine mis-reports a planted pair')
run.note('static odd: Psi = r h0\' - 2 h0, Psi_xx = V_RW Psi, V_RW = %s (det of (Psi, Psi_x) against (h0, h0\') = %s)' % (sp.factor(a_o), sp.factor(det_o.as_expr())))
lap('odd done')
stop_after('odd')

# ======================================================================= static even parity
Rr_ev = reduce_all(lin(h_even, 'even'))
ev_names = ['H0', 'H1', 'H2', 'K']

# ---- A.1 independence of P and x P'
c_, Ls_ = sp.symbols('c_ Ls_')
dep = sp.solve([c_ * (c_ - 1), Ls_ - c_ * (c_ + 1)], [c_, Ls_], dict=True)
dep_L = sorted(set(d_[Ls_] for d_ in dep))
Lq = sp.Integer(2) if mut('l1') else L                       # mutation l1: the check is run at L = 2 (l = 1)
ok_sym = (dep_L == [0, 2] and all(cancel(Lq - d_) != 0 for d_ in dep_L) and sp.factor(l * (l + 1) - 2) == (l - 1) * (l + 2))
rank_ok = True
for lv in range(2, 31):
    pl_ = sp.Poly(sp.legendre(lv, x), x)
    q_ = sp.Poly(sp.expand(x * diff(sp.legendre(lv, x), x)), x)
    rows = sp.Matrix([[pl_.coeff_monomial(x ** k) for k in range(lv + 1)], [q_.coeff_monomial(x ** k) for k in range(lv + 1)]])
    rank_ok = rank_ok and rows.rank() == 2
rank1 = sp.Matrix([[sp.Poly(sp.legendre(1, x), x).coeff_monomial(x ** k) for k in range(2)],
                   [sp.Poly(sp.expand(x * diff(sp.legendre(1, x), x)), x).coeff_monomial(x ** k) for k in range(2)]]).rank() == 1
run.check('A.1', ok_sym and rank_ok and rank1,
          'P and xP\' can be dependent at L in %s; this run has L = %s; ranks l=2..30 ok: %s, l=1 rank 1: %s' % (dep_L, Lq, rank_ok, rank1))
lap('A.1')

# ---- A.2 angular form of the seven components
XXNORM = 1 if mut('ang_norm') else x2                        # mutation ang_norm: wrong normalization of the xx component
SPEC = (('tt', (0, 0), 1, ['P0']), ('tr', (0, 1), 1, ['P0']), ('rr', (1, 1), 1, ['P0']),
        ('tx', (0, 2), 1, ['P1']), ('rx', (1, 2), 1, ['P1']),
        ('xx', (2, 2), XXNORM, ['P0', 'xP1']), ('pp', (3, 3), 1 / x2, ['P0', 'xP1']))


def even_equations(Rr):
    ok_ = all(Rr[m][n] == 0 for m, n in [(0, 3), (1, 3), (2, 3)])
    why_ = []
    out = {}
    for nm, comp, norm, basis in SPEC:
        d_, okk, w_ = coef_eqs(Rr[comp[0]][comp[1]], norm, basis, nm)
        ok_ = ok_ and okk
        if not okk:
            why_.append('%s: %s' % (nm, w_))
        for k_, v_ in d_.items():
            out[nm + '.' + k_] = jetify(v_, ev_names)
    ident1 = cancel(out['xx.xP1'] + out['pp.xP1'])
    ident2 = cancel(out['pp.P0'] - out['xx.P0'] - L * out['xx.xP1'])
    return out, ok_ and ident1 == 0 and ident2 == 0, why_, (ident1, ident2)


ev, ok_a2, why_a2, ids_a2 = even_equations(Rr_ev)
run.check('A.2', ok_a2, 'angular form of the even components: %s; identities %s' % (why_a2, ids_a2))
lap('A.2')

K0, K1, K2 = (jet_sym('K', n_) for n_ in range(3))
H10 = jet_sym('H1', 0)
Hh0, Hh1, Hh2 = (jet_sym('H', n_) for n_ in range(3))
H0_, H2_ = jet_sym('H0', 0), jet_sym('H2', 0)
if mut('s_leak'):
    ev['rx.P1'] = ev['rx.P1'] + s * K0                       # mutation s_leak: an s K term in the r-x equation
if mut('inj_tr'):
    ev['tr.P0'] = ev['tr.P0'] + K0                           # mutation inj_tr: a K term in the t-r equation
if mut('inj_tx'):
    ev['tx.P1'] = ev['tx.P1'] + r * K0                       # mutation inj_tx: a K term in the t-x equation
if mut('inj_h1'):
    ev['rr.P0'] = ev['rr.P0'] + H10                          # mutation inj_h1: an H1 term in the r-r equation
if mut('tt_drop_H'):
    ev['tt.P0'] = ev['tt.P0'] - diff(ev['tt.P0'], H0_) * H0_  # mutation tt_drop_H: the undifferentiated H0 term of the t-t equation dropped
run.note('static even radial equations before H0 = H2 (jets F_n = d^n F/dr^n; each is the coefficient of its angular structure, equal to 0):')
for nm_ in ('tt.P0', 'tr.P0', 'rr.P0', 'tx.P1', 'rx.P1', 'xx.P0', 'xx.xP1', 'pp.P0', 'pp.xP1'):
    run.note('  %-7s: %s' % (nm_, sp.factor_terms(ev[nm_])))

# ---- S.1 no time dependence left
run.check('S.1', all(s not in e_.free_symbols for e_ in ev.values()), 'the symbol s occurs in the static equations %s' % [k_ for k_, e_ in ev.items() if s in e_.free_symbols])
lap('S.1')

# ---- H.* the H1 sector
Etr, Etx = ev['tr.P0'], ev['tx.P1']
c_tr = cancel(diff(Etr, H10))
ok_h1 = jets_of(Etr) == {'H1_0'} and c_tr != 0 and cancel(Etr - c_tr * H10) == 0
run.check('H.1', ok_h1, 'delta R_tr is not c H1 alone: jets %s, c = %s' % (sorted(jets_of(Etr)), c_tr))
run.note('H.1: delta R_tr = (%s) H1  =>  H1 = 0 wherever that factor is nonzero (L != 0, r != 0, r != 2M)' % sp.factor(c_tr))
run.check('H.2', len(jets_of(Etx)) > 0 and all(nm_.startswith('H1_') for nm_ in jets_of(Etx)), 'delta R_tx has jets other than H1: %s' % sorted(jets_of(Etx)))
run.note('H.2: delta R_tx = %s  (H1 jets only; zero when H1 = 0)' % sp.factor_terms(Etx))
others = ('tt.P0', 'rr.P0', 'rx.P1', 'xx.P0', 'xx.xP1', 'pp.P0', 'pp.xP1')
leak = [nm_ for nm_ in others if any(j_.startswith('H1_') for j_ in jets_of(ev[nm_]))]
run.check('H.3', not leak, 'H1 jets occur in %s' % leak)
lap('H.*')
stop_after('h1')

# ---- E.1 traceless equation
Ec2 = ev['xx.xP1']
EtrP = cancel(ev['xx.P0'] + ev['pp.P0'])                    # angular trace (times 2)
c2_coef = cancel(diff(Ec2, H0_))
run.check('E.1', cancel(Ec2 - c2_coef * (H0_ - H2_)) == 0 and c2_coef != 0, 'traceless equation is not c0 (H0 - H2) with c0 != 0: %s' % sp.factor_terms(Ec2))
run.note('E.1: traceless equation = (%s) (H0 - H2), so H0 = H2 =: H (l >= 2 by A.1)' % sp.factor(c2_coef))
lap('E.1')

Hs = {}
for n_ in range(4):
    Hs[jet_sym('H0', n_)] = jet_sym('H', n_)
    Hs[jet_sym('H2', n_)] = jet_sym('H', n_)
eqsH = {nm: e_.subs(Hs, simultaneous=True) for nm, e_ in (('tt', ev['tt.P0']), ('rr', ev['rr.P0']), ('rx', ev['rx.P1']), ('trace', EtrP))}

# ---- E.2 the static flow from delta R_rx and delta G_rr
fsq = fq ** 2
Grr = eqsH['rr'] / 2 + eqsH['tt'] / (2 * fsq)
if not mut('grr_drop'):                                      # mutation grr_drop: the angular-trace term is dropped from delta G_rr
    Grr = Grr - eqsH['trace'] / (2 * fq * r ** 2)
Grr = cancel(together(Grr))
first_ok = {K0, K1, Hh0, Hh1, r, M, L}
ok_rx1 = eqsH['rx'].free_symbols <= first_ok
ok_G1 = Grr.free_symbols <= first_ok
unk = [K1, Hh1]
Mat2 = sp.Matrix([[diff(e_, v_) for v_ in unk] for e_ in (eqsH['rx'], Grr)]).applyfunc(cancel)
det2 = cancel(Mat2.det())
rhs2 = sp.Matrix([-(e_ - sum(diff(e_, v_) * v_ for v_ in unk)) for e_ in (eqsH['rx'], Grr)]).applyfunc(cancel)
if det2 != 0 and ok_rx1 and ok_G1:
    sol2 = Mat2.LUsolve(rhs2).applyfunc(cancel)
    FK, FH = sol2[0], sol2[1]
else:
    FK = FH = None
run.check('E.2', ok_rx1 and ok_G1 and det2 != 0,
          'static first-order equations: no second derivative in delta R_rx %s and in delta G_rr %s (free symbols %s), det %s' %
          (ok_rx1, ok_G1, sorted(map(str, Grr.free_symbols - first_ok)), det2))
lap('E.2')
if FK is None:
    run.note('E.2 failed: the static flow was not solved; the rest of the even part is not run')
    run.finish()
states = [K0, Hh0]
F2 = [FK, FH]
vm = {'K': K0, 'H': Hh0}
run.note('E.2: static flow, det of the coefficient matrix of (K\', H\') in (rx, G_rr) = %s' % sp.factor(det2))
for nm_, e_ in (('K\'', FK), ('H\'', FH)):
    run.note('E.2:   %s = %s' % (nm_, ' + '.join('(%s) %s' % (sp.factor(diff(e_, v_)), v_) for v_ in states)))

# ---- E.4 completeness: the other equations vanish identically on the flow
res4 = {}
for nm in ('tt', 'rr', 'trace'):
    res4[nm] = cancel(flow_sub(eqsH[nm], vm, states, F2))
    run.note('E.4: on the flow, equation %s = %s' % (nm, sp.factor(res4[nm]) if res4[nm] != 0 else 0))
run.check('E.4', all(v_ == 0 for v_ in res4.values()), 'equations not identically zero on the flow: %s' % [k_ for k_, v_ in res4.items() if v_ != 0])
lap('E.4')

# ---- E.5 DG06 (9) as the static master variable
rm2M = (r + 2 * M) if mut('fsign_mv') else (r - 2 * M)        # mutation fsign_mv: r - 2M -> r + 2M
lam1 = lam if mut('mv_coef') else lam + 1                    # mutation mv_coef: the second term has lam instead of lam + 1
lam3 = (lam * r + (2 if mut('mv_lam') else 3) * M)           # mutation mv_lam: lam r + 3M -> lam r + 2M in the prefactor
psi_e = cancel(r * rm2M / ((lam + 1) * lam3) * (Hh0 - r * FK) + r * K0 / lam1)
if mut('psi_zero'):
    psi_e = sp.Integer(0)                                    # mutation psi_zero: the master variable is 0
pcoef = (cancel(diff(psi_e, K0)), cancel(diff(psi_e, Hh0)))
ps_e, px_e, pxx_e, det_e, Dvec = master2(states, F2, psi_e, fq)
det_e_expr = cancel(det_e.as_expr())
run.check('E.5a', det_e != FF.zero, 'the map (K, H) -> (psi, psi_x) is singular (det = 0)')
run.note('E.5a: psi = (%s) K + (%s) H;  det of (psi, psi_x) against (K, H) = %s' % (sp.factor(pcoef[0]), sp.factor(pcoef[1]), sp.factor(det_e_expr)))
if det_e != FF.zero:
    a_e, b_e = ab_from(ps_e, px_e, pxx_e, det_e)
else:
    a_e, b_e = None, None
run.check('E.5b', b_e is not None and cancel(b_e) == 0, 'coefficient of psi_x in psi_xx = a psi + b psi_x: b = %s' % (sp.factor(b_e) if b_e is not None else None))
nine = 9 * lam * M ** (1 if mut('v_typo') else 2) * r        # mutation v_typo: 9 lam M^2 r -> 9 lam M r
VZ = 2 * fq * (lam ** 2 * (lam + 1) * r ** 3 + 3 * lam ** 2 * M * r ** 2 + nine + 9 * M ** 3) / (r ** 3 * (lam * r + 3 * M) ** 2)
run.check('E.5d', a_e is not None and same(a_e, VZ) and same(a_e.subs(L, l * (l + 1)), VZ.subs(L, l * (l + 1))),
          'derived V differs from the Zerilli form: difference %s' % (sp.factor(cancel(a_e - VZ)) if a_e is not None else None))
run.note('E.5: V = %s' % (sp.factor(a_e) if a_e is not None else None))
cc = L - 2                                                   # c = (l+2)(l-1) of dg_proof
lamR = cc / 2
res9 = 9 * (lamR * r - (2 if mut('res_typo') else 1))        # mutation res_typo: 9(lam r - 1) -> 9(lam r - 2)
V_res = 2 * (1 + 2 / r) * (lamR ** 2 * r ** 2 * ((lamR + 1) * r - 3) + res9) / (r ** 3 * (lamR * r - 3) ** 2)
run.check('E.5e', a_e is not None and same(a_e.subs(M, -1), V_res), 'derived V at M = -1 differs from RESULTS 4.2: %s' % (sp.factor(cancel(a_e.subs(M, -1) - V_res)) if a_e is not None else None))
lap('E.5')

# ---- E.6 the converse (DG06 (3)-(4) at s = 0) and the inversion
fF = 1 - 2 * MF / rF
Fk = toF((lam * (lam + 1) * r ** 2 + 3 * lam * M * r + 6 * M ** 2) / (r ** 2 * (lam * r + 3 * M)))                      # DG06 (4) f(r)
hk = toF((lam * r ** 2 - 3 * lam * r * M - (2 if mut('dg3_typo') else 3) * M ** 2) / ((r - 2 * M) * (lam * r + 3 * M)))  # DG06 (4) h(r)
VF = toF(VZ)
c0_ = VF / fF ** 2                                            # psi'' = c0 psi + c1 psi' from f (f psi')' = V psi
c1_ = -(2 * MF / rF ** 2) / fF


def DPair(p_):
    a_, b_ = p_
    return (a_.diff(rF) + b_ * c0_, a_ + b_.diff(rF) + b_ * c1_)


def addP(p_, q_, cp=1, cq=1):
    return (cp * p_[0] + cq * q_[0], cp * p_[1] + cq * q_[1])


Kp = (Fk, fF)                                                 # (3): K = F psi + f psi'
H1p = (FF.zero, FF.zero)                                      # (3): H1 = h psi_t + r psi_tr = 0 for a static psi
inner = (hk, rF)
finner = (fF * inner[0], fF * inner[1])
H2p = addP(DPair(finner), Kp, 1, -1)                          # (3): H2 = d/dr[f (h psi + r psi')] - K; H0 = H2
jets_dg = {}
for nm_, e0 in (('K', Kp), ('H1', H1p), ('H2', H2p), ('H0', H2p)):
    cur = e0
    for n_ in range(3):
        jets_dg[jet_sym(nm_, n_)] = cur
        cur = DPair(cur)
res6 = {}
for nm, e_ in (('tt', ev['tt.P0']), ('tr', ev['tr.P0']), ('rr', ev['rr.P0']), ('tx', ev['tx.P1']), ('rx', ev['rx.P1']), ('trace', EtrP), ('traceless', Ec2)):
    e_ = sp.expand(e_)
    tot = (FF.zero, FF.zero)
    for sy in jets_dg:
        cj = e_.coeff(sy, 1)
        if cj != 0:
            tot = addP(tot, jets_dg[sy], 1, toF(cj))
    res6[nm] = tot
run.check('E.6a', all(v_[0] == 0 and v_[1] == 0 for v_ in res6.values()),
          'DG06 (3)-(4) at s = 0 does not solve: nonzero equations %s' % [k_ for k_, v_ in res6.items() if v_[0] != 0 or v_[1] != 0])
lap('E.6a')
# (9) o (3): the typed copy of (9) (mutation inv_typo: lam + 1 -> lam + 2 in the prefactor)
pre9 = toF(r * (r - 2 * M) / (((lam + 2) if mut('inv_typo') else (lam + 1)) * (lam * r + 3 * M)))
back = addP((pre9 * (jets_dg[jet_sym('H2', 0)][0] - rF * jets_dg[jet_sym('K', 1)][0]), pre9 * (jets_dg[jet_sym('H2', 0)][1] - rF * jets_dg[jet_sym('K', 1)][1])),
            (toF(r / (lam + 1)) * jets_dg[jet_sym('K', 0)][0], toF(r / (lam + 1)) * jets_dg[jet_sym('K', 0)][1]))
ok_b1 = back[0] == 1 and back[1] == 0
# (3) o (9): linear forms in (K, H) along the static flow; psi = ps_e, psi_x = px_e
Kn = [Fk * ps_e[j] + px_e[j] for j in range(2)]
inn = [fF * hk * ps_e[j] + rF * px_e[j] for j in range(2)]
Dinn = Dvec(inn)
Hn = [Dinn[j] - Kn[j] for j in range(2)]
ok_b2 = det_e != FF.zero and Kn[0] == 1 and Kn[1] == 0 and Hn[0] == 0 and Hn[1] == 1
run.check('E.6b', ok_b1 and ok_b2, '(9) o (3) = id: %s (got %s); (3) o (9) = id on (K, H): %s (K -> %s, H -> %s)' % (ok_b1, back, ok_b2, Kn, Hn))
lap('E.6b')
stop_after('even')

# ======================================================================= the s -> 0 relation to #143
# the e^{st} equations of the same code, H1 = s h1t (h1t = H1/s), coefficient equations at s = 0
h_dep = even_ansatz(Edep, gtr=(2 if mut('dep_typo') else 1))  # mutation dep_typo: g_tr = 2 eps H1 P e^{st}
Rr_dep = reduce_all(lin(h_dep, 'even_dep'), Edep)
ev_dep, ok_dep, why_dep, ids_dep = even_equations(Rr_dep)
lap('time-dependent equations')
diffs = [k_ for k_ in ev if cancel(ev_dep[k_].subs(s, 0) - ev[k_]) != 0]
run.check('S.0', ok_dep and not diffs, 'static equations differ from the s = 0 values of the e^{st} equations at %s (angular form of the e^{st} components ok: %s %s)' % (diffs, ok_dep, why_dep))
lap('S.0')

ht = [jet_sym('h1t', n_) for n_ in range(3)]
sub_h1 = {jet_sym('H1', n_): s * ht[n_] for n_ in range(3)}
sub_all = dict(sub_h1)
sub_all.update(Hs)


def at_s0(e):
    """(e with H1 = s h1t, H0 = H2 = H) / s at s = 0"""
    return cancel(cancel(together(e.subs(sub_all, simultaneous=True)) / s).subs(s, 0))


T1 = at_s0(ev_dep['tr.P0'])
mu1 = cancel(diff(T1, K1))
typed_tr = K1 + (r - 3 * M) / (r * (r - 2 * M)) * K0 - Hh0 / r - L / ((3 if mut('tr_typo') else 2) * r ** 2) * ht[0]   # mutation tr_typo: L/(2 r^2) -> L/(3 r^2)
run.check('S.R1', mu1 != 0 and jets_of(T1) <= {'K_0', 'K_1', 'H_0', 'h1t_0'} and cancel(T1 / mu1 - typed_tr) == 0,
          '(tr/s) at s = 0 is not a multiple of the tr line of #143: T1 = %s' % sp.factor(T1))
run.note('S.R1: (delta R_tr / s)|_{s=0} = (%s) [K\' + (r-3M)/(r(r-2M)) K - H/r - L h1t/(2r^2)], h1t = H1/s' % sp.factor(mu1))
ht_sol = cancel(sp.solve(T1, ht[0])[0])                      # h1t from the s -> 0 tr equation (solved, not typed)
if mut('h1t_typo'):
    ht_sol = Q(3, 2) * ht_sol                                # mutation h1t_typo: the solved h1t multiplied by 3/2
ht_form = cancel(ht_sol.subs(K1, FK))                        # h1t as a linear form in (K, H) on the static flow
T2 = at_s0(ev_dep['tx.P1'])
T2f = cancel(flow_sub(T2, {'K': K0, 'H': Hh0, 'h1t': ht_form}, states, F2))
run.check('S.R2', T2f == 0, '(tx/s) at s = 0 does not vanish on the static flow with h1t from tr: %s' % sp.factor(T2f))
run.note('S.R2: h1t = %s;  (delta R_tx / s)|_{s=0} = %s' % (sp.factor(ht_form), sp.factor_terms(T2)))
lam3t = lam * r + (2 if mut('sr_typo') else 3) * M           # mutation sr_typo: lam r + 3M -> lam r + 2M in the typed #143 map
psi143 = cancel((r ** 2 * K0 - (r - 2 * M) * ht_form) / lam3t)
run.check('S.R3', cancel(psi143 - psi_e) == 0, 'DG06 (9) on the static flow differs from [r^2 K - (r-2M) h1t]/(lam r + 3M): %s' % sp.factor(cancel(psi143 - psi_e)))
lap('S.R')
stop_after('rel')

# ======================================================================= D.* the divisions
DIV = []


def div(step, desc, expr):
    DIV.append((step, desc, sp.factor(expr)))


div('H.1: H1 = 0', 'coefficient c_tr of H1 in delta R_tr', c_tr)
div('E.1: H0 = H2', 'coefficient of (H0 - H2) in the traceless equation', c2_coef)
det2_D5 = sp.Integer(0) if mut('d5_zero') else det2           # mutation d5_zero: the determinant of E.2 is 0
div('E.2: flow', 'det of the coefficient matrix of (K\', H\') in (rx, G_rr)', det2)
for nm_, e_ in (('E.2: flow K\'', FK), ('E.2: flow H\'', FH), ('E.5: psi coefficient of K', pcoef[0]), ('E.5: psi coefficient of H', pcoef[1]),
                ('E.5: V', a_e), ('S.R2: h1t', ht_form)):
    for fac_, mult_ in denominators_factors(e_):
        DIV.append((nm_, 'denominator factor (multiplicity %d)' % mult_, fac_))
DIV.append(('E.5: definition (9)', 'denominator (lam + 1) of both terms', sp.factor(lam + 1)))
DIV.append(('E.5: definition (9)', 'denominator (lam r + 3M)', sp.factor(lam * r + 3 * M)))
DIV.append(('E.5: definition (9)', 'denominator (r - 2M) of the inverse formula (4)', sp.factor(r - 2 * M)))
if mut('div_s'):
    DIV.append(('mutation div_s', 'a divisor s', s))          # mutation div_s: s appended to the divisors met
if mut('div_l2'):
    DIV.append(('mutation div_l2', 'a divisor L - 2', L - 2))  # mutation div_l2: L - 2 appended to the divisors met
lam_fac = (L - 2) * r + (7 if mut('div_set') else 6) * M     # lam r + 3M up to the factor 2; mutation div_set: 6M -> 7M in the stated set
ALLOWED = [s, r, r - 2 * M, L, L - 2, lam_fac, M]


def is_const(e_):
    e_ = cancel(e_)
    return e_.free_symbols == set() or e_.is_number


def classify(fac_):
    for a_ in ALLOWED:
        if is_const(fac_ / a_):
            return a_
    return None


facs_seen = {}
for step, desc, e_ in DIV:
    num_, den_ = sp.fraction(cancel(together(e_)))
    for fac_, mult_ in sp.factor_list(num_)[1] + sp.factor_list(den_)[1]:
        a_ = classify(fac_)
        facs_seen.setdefault(str(a_), []).append((step, desc))
        if a_ is None:
            facs_seen.setdefault('UNEXPECTED', []).append((step, desc, fac_))
run.note('divisors met (step | what | expression):')
for step, desc, e_ in DIV:
    run.note('  %s | %s | %s' % (step, desc, e_))
run.check('D.1', 'UNEXPECTED' not in facs_seen, 'divisor factors outside the vocabulary {s, r, r-2M, L, L-2, lam r + 3M, M}: %s' % facs_seen.get('UNEXPECTED'))
for key in sorted(k_ for k_ in facs_seen if k_ != 'UNEXPECTED'):
    run.note('  factor %s met at %d step(s)' % (key, len(facs_seen[key])))
lap('D.1')

r0_ = sp.solve(lam * r + (2 if mut('d2_zero') else 3) * M, r)[0]   # mutation d2_zero: the zero of lam r + 2M
Vnum_at = cancel(cancel(a_e * (lam * r + 3 * M) ** 2).subs(r, r0_))
ok_d2 = (cancel(r0_ + 3 * M / lam) == 0 and Vnum_at != 0 and cancel(r0_.subs(M, -1) - 6 / cc) == 0
         and r0_.subs({M: -1, L: 6}) > 0 and r0_.subs({M: 1, L: 6}) < 0 and cancel(Vnum_at - (L - 2) ** 2 * (L + 1) ** 2 / 18) == 0)
run.check('D.2', ok_d2, 'lam r + 3M: zero at %s; (lam r + 3M)^2 V there = %s; at M = -1 the zero is %s (6/c = %s)' % (r0_, Vnum_at, r0_.subs(M, -1), 6 / cc))
run.note('D.2: lam r + 3M = 0 at r = %s; (lam r + 3M)^2 V at that r = %s (nonzero: a genuine double pole); at M = -1 the zero is 6/c = r_s of dg_proof' % (r0_, sp.factor(Vnum_at)))

run.check('D.3', str(s) not in facs_seen, 'a step divides by s: %s' % facs_seen.get(str(s)))
run.note('D.3: no step of the static derivation divides by s; the symbol s occurs in none of the static equations (S.1)')

facs_f = facs_seen.get(str(r - 2 * M), [])
facs_l2 = facs_seen.get(str(L - 2), [])
facs_L = facs_seen.get(str(L), [])
run.check('D.4', len(facs_l2) == 0, 'steps dividing by L - 2 = 2 lam: %s' % facs_l2)
run.note('D.4: L = l(l+1) is a divisor at: %s' % '; '.join(sorted(set(st_ for st_, _ in facs_L))))
run.note('D.4: r - 2M (= r f) is a divisor at: %s; r - 2M > 0 for r > 0 and M < 0, and vanishes at r = 2M for M > 0' % '; '.join(sorted(set(st_ for st_, _ in facs_f))))
run.note('D.4: the restriction l >= 2 enters through A.1 and the nonvanishing of c_tr (H1 = 0) and c0 (H0 = H2), not through a division by L - 2')

run.check('D.5', det2_D5 != 0 and c_tr != 0 and det_e_expr != 0 and c2_coef != 0, 'a determinant or a coefficient is identically 0')
run.note('D.5: det of the static flow system = %s; c_tr = %s; det of (psi, psi_x) against (K, H) = %s' % (sp.factor(det2), sp.factor(c_tr), sp.factor(det_e_expr)))
lap('D')

run.finish()
