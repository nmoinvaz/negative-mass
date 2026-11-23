"""#143: the even-parity (Zerilli) potential from the linearized Einstein equations, M symbolic, then M = -1.

Metric g = g_Schw(M) + eps h, h in Regge-Wheeler gauge, axisymmetric, angular dependence P_l(cos theta) with Legendre's equation used to remove
second derivatives (l symbolic), time dependence e^{st}. Own Christoffel and Ricci code (zerilli_derivation_lib.py), first order in eps by dual numbers.
Order: planted checks (exit 2 on failure) -> odd parity (Regge-Wheeler) -> even parity (Zerilli).
Sources: Dotti-Gleiser gr-qc/0604021 eqs. (2)-(9) (data; typed, see I143-rpi5/README.md); Zerilli's papers were not read (publisher HTTP 403).
Run: python zerilli_derivation.py [--part=planted|odd|even|all] [--mutate=NAME] [--failfast] > zerilli_derivation_output.txt
Mutations: fsign_bg fsign_mv mv_coef drop_trace v_typo rw_typo l1 dual_inv leg_sign route2_drop res_typo dg3_typo skip_s
"""
import sys
import time
from zerilli_derivation_lib import *

PART = 'all'
SKIP = []
for _a in sys.argv[1:]:
    if _a.startswith('--part='):
        PART = _a[len('--part='):]
    if _a.startswith('--skip='):
        SKIP = _a[len('--skip='):].split(',')


def skipped(cid):
    return any(cid.startswith(p) for p in SKIP)

run = Run('zerilli_derivation: Zerilli potential from the linearized Einstein equations (part=%s)' % PART)
T0 = time.time()


def lap(msg):
    print('# [%6.1f s] %s' % (time.time() - T0, msg), flush=True)


E = sp.exp(s * t)
f = f_bg()                      # background metric's f (mutation fsign_bg); the rest of the code uses fM()
fq = fM()                       # f = 1 - 2M/r, always the right sign (used by the ansatz and the master variable)
g0 = sp.diag(-f, 1 / f, r ** 2 / (1 - x ** 2), r ** 2 * (1 - x ** 2))
Px = sp.Derivative(P, x)


def zero_mat(mat):
    return all(tidy(mat[m][n]) == 0 for m in range(N) for n in range(N))


def red(e):
    """first-order component -> divide e^{st}, Legendre reduction, return an expression in P0, P1 with rational coefficients."""
    e = sp.expand(e / E)
    return tidy(leg_red(e))


def reduce_all(R1):
    return [[red(R1[m][n]) if n >= m else None for n in range(N)] for m in range(N)]


# ======================================================================= planted checks
# PL.0a dual inverse
hg = [[Function('h%d%d' % (min(i, j), max(i, j)))(t, r, x) for j in range(N)] for i in range(N)]
g_, gi_ = build_metric(g0, hg)
prod0 = sp.Matrix(N, N, lambda i, j: sum((g_[i][k] * gi_[k][j]).a0 for k in range(N)))
prod1 = sp.Matrix(N, N, lambda i, j: sum((g_[i][k] * gi_[k][j]).a1 for k in range(N)))
run.planted('PL.0a', prod0.applyfunc(tidy) == sp.eye(N) and prod1.applyfunc(tidy) == sp.zeros(N, N),
            'g g^-1 != 1 at order eps: the inverse used is not the first-order inverse')
lap('PL.0a')

# PL.0b planted known Ricci: static de Sitter
a_ = Symbol('a', positive=True)
fds = 1 - r ** 2 / a_ ** 2
g_ds = sp.diag(-fds, 1 / fds, r ** 2 / (1 - x ** 2), r ** 2 * (1 - x ** 2))
R0_ds, _ = lin_ricci(g_ds, sp.zeros(N, N))
run.planted('PL.0b', all(tidy(R0_ds[m][n] - 3 / a_ ** 2 * g_ds[m, n]) == 0 for m in range(N) for n in range(N)),
            'Ricci of static de Sitter is not (3/a^2) g')
lap('PL.0b')

# PL.0c background vacuum
R0_bg, _ = lin_ricci(g0, sp.zeros(N, N))
run.planted('PL.0c', all(tidy(R0_bg[m][n]) == 0 for m in range(N) for n in range(N)), 'Ricci of g_Schw(M) is not 0')
lap('PL.0c')

# PL.0d Legendre reduction against explicit P_l
ok = True
for n_ in range(2, 7):
    for order in (2, 3, 4):
        red_ = leg_red(sp.Derivative(P, (x, order)))
        P0s, P1s = sp.symbols('P0 P1')
        lhs = red_.subs({P0s: sp.legendre(n_, x), P1s: diff(sp.legendre(n_, x), x), l: n_})
        ok = ok and cancel(lhs - diff(sp.legendre(n_, x), x, order)) == 0
run.planted('PL.0d', ok, 'Legendre reduction differs from the explicit derivative of P_l')
lap('PL.0d')

# ----------------------------------------------------------------------- ansatz and gauge fields
r_funcs = {n: Function(n)(r) for n in ('H0', 'H1', 'H2', 'K', 'h0', 'h1', 'xa', 'xb', 'xc', 'xd')}
H0, H1, H2, K, h0, h1, xa, xb, xc, xd = (r_funcs[n] for n in ('H0', 'H1', 'H2', 'K', 'h0', 'h1', 'xa', 'xb', 'xc', 'xd'))


def zeros():
    return [[0] * N for _ in range(N)]


def sym_set(h, m, n, val):
    h[m][n] = val
    h[n][m] = val


h_even = zeros()
sym_set(h_even, 0, 0, fq * H0 * P * E)
sym_set(h_even, 0, 1, H1 * P * E)
sym_set(h_even, 1, 1, H2 * P * E / fq)
sym_set(h_even, 2, 2, r ** 2 * K * P * E / (1 - x ** 2))
sym_set(h_even, 3, 3, r ** 2 * (1 - x ** 2) * K * P * E)

h_odd = zeros()
sym_set(h_odd, 0, 3, h0 * E * (1 - x ** 2) * Px)
sym_set(h_odd, 1, 3, h1 * E * (1 - x ** 2) * Px)


def lie(g0m, xi):
    out = zeros()
    for m in range(N):
        for n in range(N):
            out[m][n] = sum(xi[a] * diff(g0m[m, n], COORD[a]) + g0m[a, n] * diff(xi[a], COORD[m]) + g0m[m, a] * diff(xi[a], COORD[n])
                            for a in range(N))
    return out


xi_even = [xa * P * E, xb * P * E, xc * (1 - x ** 2) * Px * E / r ** 2, 0]
xi_odd = [0, 0, 0, xd * Px * E]
h_gauge_even = lie(g0, xi_even)
h_gauge_odd = lie(g0, xi_odd)


_LIN = {}


def lin(h, key):
    """dual-number route: unreduced first-order Ricci of g0 + eps h (cached by key)."""
    if key not in _LIN:
        _, _LIN[key] = lin_ricci(g0, sp.Matrix(h))
    return _LIN[key]


def mat_h(h):
    return sp.Matrix(h)


# PL.G pure gauge
for nm, hh in (('even', h_gauge_even), ('odd', h_gauge_odd)):
    if skipped('PL.G'):
        break
    R1 = lin(hh, 'gauge_' + nm)
    Rr = reduce_all(R1)
    run.planted('PL.G.' + nm, all(Rr[m][n] == 0 for m in range(N) for n in range(m, N)),
                'pure-gauge h (%s) does not give dR = 0: %s' % (nm, [(m, n) for m in range(N) for n in range(m, N) if Rr[m][n] != 0]))
    lap('PL.G.' + nm)

# PL.X route 2
Rx = {}
for nm, hh in (('even', h_even), ('odd', h_odd), ('gauge_even', h_gauge_even)):
    if skipped('PL.X'):
        break
    R1 = lin(hh, nm if nm != 'gauge_even' else 'gauge_even')
    R1b = route2_delta_ricci(g0, hh)
    bad = [(m, n) for m in range(N) for n in range(m, N) if red(R1[m][n] - R1b[m][n]) != 0]
    run.planted('PL.X.' + nm, not bad, 'route 1 (dual) and route 2 (covariant) differ at components %s' % bad)
    Rx[nm] = R1
    lap('PL.X.' + nm)

if PART == 'planted':
    run.finish()


# ======================================================================= odd parity (planted target: Regge-Wheeler)
x2 = 1 - x ** 2
DIV = []                         # (step, description, expression): every division made in the reduction


def div(step, desc, expr):
    DIV.append((step, desc, sp.factor(expr)))


def same(e1, e2):
    return cancel(together(e1 - e2)) == 0


Rr_odd = reduce_all(lin(h_odd, 'odd'))
odd_names = ['h0', 'h1']
ok_o1 = all(Rr_odd[m][n] == 0 for m, n in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2), (3, 3)])
why_o1 = []
oeq = {}
for nm, comp, norm, basis in (('tphi', (0, 3), 1 / x2, ['P1']), ('rphi', (1, 3), 1 / x2, ['P1']), ('xphi', (2, 3), 1, ['P0', 'xP1'])):
    d_, ok_, why_ = coef_eqs(Rr_odd[comp[0]][comp[1]], norm, basis, nm)
    ok_o1 = ok_o1 and ok_
    if not ok_:
        why_o1.append('%s: %s' % (nm, why_))
    oeq[nm] = {k: jetify(v, odd_names) for k, v in d_.items()}
ok_o1 = ok_o1 and 'P0' in oeq['xphi'] and cancel(oeq['xphi']['P0'] + L / 2 * oeq['xphi']['xP1']) == 0
run.planted('PL.O.1', ok_o1, 'odd components do not have the expected angular form: %s' % why_o1)
Etphi, Erphi, Exphi = oeq['tphi']['P1'], oeq['rphi']['P1'], oeq['xphi']['xP1']
run.note('odd radial equations (jets h0_n, h1_n):')
for nm_, e_ in (('t-phi', Etphi), ('r-phi', Erphi), ('x-phi', Exphi)):
    run.note('  %s: %s' % (nm_, sp.factor(e_)))


def Dj(e, names):
    out = diff(e, r)
    for nm_ in names:
        for n_ in range(6):
            sy = jet_sym(nm_, n_)
            if e.has(sy):
                out += diff(e, sy) * jet_sym(nm_, n_ + 1)
    return out


# h0 from the x-phi equation (algebraic in h0)
c_h0 = cancel(diff(Exphi, jet_sym('h0', 0)))
ok_h0 = cancel(Exphi - c_h0 * jet_sym('h0', 0)).free_symbols.isdisjoint({jet_sym('h0', 1), jet_sym('h0', 2)})
sol_h0 = cancel(-(Exphi - c_h0 * jet_sym('h0', 0)) / c_h0)
div('odd: h0 from x-phi', 'coefficient of h0 in the x-phi equation', c_h0)
h0e = {0: sol_h0}
for n_ in (1, 2):
    h0e[n_] = cancel(Dj(h0e[n_ - 1], ['h1']))
Er2 = cancel(Erphi.subs({jet_sym('h0', 0): h0e[0], jet_sym('h0', 1): h0e[1]}, simultaneous=True))
c_h1pp = cancel(diff(Er2, jet_sym('h1', 2)))
F2odd = cancel(-(Er2 - c_h1pp * jet_sym('h1', 2)) / c_h1pp)
div('odd: ODE of h1', 'coefficient of h1_2 after h0 is eliminated', c_h1pp)
Et2 = Etphi.subs({jet_sym('h0', 0): h0e[0], jet_sym('h0', 1): h0e[1], jet_sym('h0', 2): h0e[2]}, simultaneous=True)
Et2 = cancel(flow_sub(Et2, {'h1': jet_sym('h1', 0)}, [jet_sym('h1', 0), jet_sym('h1', 1)], [jet_sym('h1', 1), F2odd]))
run.planted('PL.O.2', ok_h0 and Et2 == 0,
            'odd reduction: h0 algebraic %s; t-phi equation on the h1 flow = %s' % (ok_h0, Et2))
u_o = [jet_sym('h1', 0), jet_sym('h1', 1)]
psi_o = fq * u_o[0] / r                                     # Psi = f h1 / r
a_o, b_o, det_o = master_reduce(u_o, [u_o[1], F2odd], psi_o, fq)
div('odd: master variable', 'det of (Psi, Psi_x) against (h1, h1_r)', det_o)
V_o = cancel(a_o - s ** 2)
V_RW = fq * (L / r ** 2 - (3 if mut('rw_typo') else 6) * M / r ** 3)
run.planted('PL.O.3', cancel(b_o) == 0 and diff(V_o, s) == 0 and same(V_o, V_RW),
            'odd: coefficient of Psi_x = %s, dV/ds = %s, V - V_RW = %s' % (b_o, diff(V_o, s), cancel(V_o - V_RW)))
run.planted('PL.O.4', not same(V_RW, V_RW + M / r ** 3) and same(V_RW, sp.expand(V_RW)), 'comparison routine mis-reports a planted pair')
run.note('odd: Psi = f h1 / r, Psi_xx = (V_RW + s^2) Psi, V_RW = %s' % sp.factor(V_o))
lap('odd done')

if PART == 'odd':
    run.finish()

# ======================================================================= even parity
Rr_ev = reduce_all(lin(h_even, 'even'))
ev_names = ['H0', 'H1', 'H2', 'K']

# ---- A.1 independence of P and x P' (needed to read the coefficients of the two structures as separate equations)
# a P + b x P' = 0 with (a, b) constants and P != 0 forces x P' = c P, so P = A x^c; the Legendre equation then gives
# (1-x^2)P'' - 2xP' + L P = A x^(c-2) [c(c-1) + (L - c(c+1)) x^2], which vanishes iff c(c-1) = 0 and L = c(c+1).
c_, Ls_ = sp.symbols('c_ Ls_')
dep = sp.solve([c_ * (c_ - 1), Ls_ - c_ * (c_ + 1)], [c_, Ls_], dict=True)
dep_L = sorted(set(d_[Ls_] for d_ in dep))                  # the values of L at which P and xP' can be dependent: [0, 2]
Lq = sp.Integer(2) if mut('l1') else L                      # mutation l1: the check is run at L = 2 (l = 1)
ok_sym = dep_L == [0, 2] and all(cancel(Lq - d_) != 0 for d_ in dep_L)
rank_ok = True
for lv in range(2, 31):
    pl_ = sp.Poly(sp.legendre(lv, x), x)
    q_ = sp.Poly(sp.expand(x * diff(sp.legendre(lv, x), x)), x)
    rows = sp.Matrix([[pl_.coeff_monomial(x ** k) for k in range(lv + 1)], [q_.coeff_monomial(x ** k) for k in range(lv + 1)]])
    rank_ok = rank_ok and rows.rank() == 2
rank1 = sp.Matrix([[0, 1], [0, 1]]).rank() == 1                # l = 1: P = x and x P' = x
run.check('A.1', ok_sym and rank_ok and rank1,
          'P and xP\' can be dependent at L in %s; this run has L = %s; ranks l=2..30 ok: %s, l=1 rank 1: %s' % (dep_L, Lq, rank_ok, rank1))
lap('A.1')

# ---- A.2 angular form of the seven components
x2 = 1 - x ** 2
SPEC = (('tt', (0, 0), 1, ['P0']), ('tr', (0, 1), 1, ['P0']), ('rr', (1, 1), 1, ['P0']),
        ('tx', (0, 2), 1, ['P1']), ('rx', (1, 2), 1, ['P1']),
        ('xx', (2, 2), x2, ['P0', 'xP1']), ('pp', (3, 3), 1 / x2, ['P0', 'xP1']))
ok_a2 = all(Rr_ev[m][n] == 0 for m, n in [(0, 3), (1, 3), (2, 3)])
why_a2 = []
ev = {}
for nm, comp, norm, basis in SPEC:
    d_, ok_, why_ = coef_eqs(Rr_ev[comp[0]][comp[1]], norm, basis, nm)
    ok_a2 = ok_a2 and ok_
    if not ok_:
        why_a2.append('%s: %s' % (nm, why_))
    for k_, v_ in d_.items():
        ev[nm + '.' + k_] = jetify(v_, ev_names)
# the tensor harmonic Z_ab = nabla_a nabla_b P + (L/2) gamma_ab P has (1-x^2) Z_xx = x P' - (L/2) P and Z_pp/(1-x^2) = -(x P' - (L/2) P):
# the traceless part gives xx.xP1 = -pp.xP1 =: c2 and pp.P0 - xx.P0 = L c2
ident1 = cancel(ev['xx.xP1'] + ev['pp.xP1'])
ident2 = cancel(ev['pp.P0'] - ev['xx.P0'] - L * ev['xx.xP1'])
run.check('A.2', ok_a2 and ident1 == 0 and ident2 == 0,
          'angular form of the even components: %s; identities (%s, %s)' % (why_a2, ident1, ident2))
lap('A.2')
Ett, Etr, Err, Etx, Erx = ev['tt.P0'], ev['tr.P0'], ev['rr.P0'], ev['tx.P1'], ev['rx.P1']
Ec2 = ev['xx.xP1']                                          # traceless angular equation
EtrP = cancel(ev['xx.P0'] + ev['pp.P0'])                    # angular trace (times 2)
run.note('even radial equations: 7 (tt, tr, rr, tx, rx, angular trace, traceless)')

# ---- E.1 traceless equation
H0_, H2_ = jet_sym('H0', 0), jet_sym('H2', 0)
c2_coef = cancel(diff(Ec2, H0_))
ok_e1 = cancel(Ec2 - c2_coef * (H0_ - H2_)) == 0 and c2_coef != 0
run.check('E.1', ok_e1, 'traceless equation is not c0 (H0 - H2) with c0 != 0: %s' % sp.factor(Ec2))
run.note('E.1: traceless equation = (%s) (H0 - H2), so H0 = H2 =: H (l >= 2 by A.1)' % c2_coef)
lap('E.1')

# H0 = H2 = H in the other six
Hs = {}
for n_ in range(4):
    Hs[jet_sym('H0', n_)] = jet_sym('H', n_)
    Hs[jet_sym('H2', n_)] = jet_sym('H', n_)
eqsH = {nm: e_.subs(Hs, simultaneous=True) for nm, e_ in (('tt', Ett), ('tr', Etr), ('rr', Err), ('tx', Etx), ('rx', Erx), ('trace', EtrP))}
K0, K1, K2 = (jet_sym('K', n_) for n_ in range(3))
H10, H11, H12 = (jet_sym('H1', n_) for n_ in range(3))
Hh0, Hh1, Hh2 = (jet_sym('H', n_) for n_ in range(3))

# ---- E.2 first-order equations and the algebraic relation
first = ('tr', 'tx', 'rx')
allowed = {K0, K1, H10, H11, Hh0, Hh1}
ok_first = all(eqsH[nm].free_symbols <= allowed | {r, s, l, M} for nm in first)
unk = [K1, H11, Hh1]
Mat3 = sp.Matrix([[diff(eqsH[nm], v) for v in unk] for nm in first]).applyfunc(cancel)
det3 = cancel(Mat3.det())
rhs3 = sp.Matrix([-(eqsH[nm] - sum(diff(eqsH[nm], v) * v for v in unk)) for nm in first]).applyfunc(cancel)
sol3 = Mat3.LUsolve(rhs3).applyfunc(cancel)
div('E.2: flow', 'det of the coefficient matrix of (K\', H1\', H\') in (tr, tx, rx)', det3)
for k_, nm_ in enumerate(('K', 'H1', 'H')):
    for d_ in sp.factor_list(sp.fraction(cancel(together(sol3[k_])))[1])[1]:
        pass
FK, FH1, FH = sol3[0], sol3[1], sol3[2]
states3 = [K0, H10, Hh0]
F3 = [FK, FH1, FH]
vm3 = {'K': K0, 'H1': H10, 'H': Hh0}
fsq = fq ** 2
Grr = Err / 2 + Ett / (2 * fsq) - EtrP / (2 * fq * r ** 2)  # delta G_rr = delta R_rr - (1/2) g_rr (g^ab delta R_ab); delta R_ab = 0 in vacuum
Grr = Grr.subs(Hs, simultaneous=True)
ok_G_first = Grr.free_symbols <= {K0, K1, K2, H10, H11, Hh0, Hh1, Hh2, r, s, l, M}
Aflow = cancel(flow_sub(Grr, vm3, states3, F3))
run.check('E.2', ok_first and det3 != 0 and Aflow.free_symbols <= set(states3) | {r, s, l, M},
          'first-order equations: no second derivatives %s, det %s, A free symbols %s' % (ok_first, det3, Aflow.free_symbols))
aK, a1, aH = (cancel(diff(Aflow, v)) for v in states3)
if mut('skip_s'):
    a1 = 2 * a1                                           # mutation skip_s: H1 coefficient of A scaled by 2 after the derivation
A = aK * K0 + a1 * H10 + aH * Hh0
run.note('E.2: algebraic relation A = %s' % sp.factor(A))
lap('E.2')

# ---- E.3 consistency of A along the flow
def D3(e):
    return diff(e, r) + sum(diff(e, states3[i]) * F3[i] for i in range(3))
dA = cancel(D3(A))
kappa = cancel(diff(dA, Hh0) / aH)
run.check('E.3', cancel(dA - kappa * A) == 0, 'dA/dr - kappa A = %s' % cancel(dA - kappa * A))
run.note('E.3: dA/dr = (%s) A' % sp.factor(kappa))
lap('E.3')

# ---- E.4 the other equations on the flow are multiples of A
res4 = {}
for nm in ('tt', 'tr', 'rr', 'tx', 'rx', 'trace'):
    ef = cancel(flow_sub(eqsH[nm], vm3, states3, F3))
    mu = cancel(diff(ef, Hh0) / aH)
    res4[nm] = (mu, cancel(ef - mu * A))
ok_e4 = all(v[1] == 0 for v in res4.values())
run.check('E.4', ok_e4, 'equations not multiples of A: %s' % [k_ for k_, v in res4.items() if v[1] != 0])
lap('E.4')
div('E.2: H from A', 'coefficient of H in A', aH)
