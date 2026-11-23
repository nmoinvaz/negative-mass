"""#146: (1) the Regge-Wheeler gauge is attainable for even-parity (l >= 2) perturbations, with e^{st} and without time dependence, exactly,
M and L symbolic; (2) #143's flow at r_s = -3M/lambda against Dotti-Gleiser's Lemma 1 condition (cz) on zeta = (2 lam r - 3 mu) psi_Z, mu = -2M.

Library: zerilli_derivation_lib.py of #143, imported unchanged. The flow of part R is derived again here with it (same code path as #143's E.1-E.2), not typed;
#143's and #145's printed flows are typed only as comparison targets (R.1, R.S.1).
Order: planted checks (exit 2 on failure) -> gauge e^{st} -> gauge static -> r_s (flow, zeta, Lemma 1, static).
Sources: Dotti-Gleiser gr-qc/0604021 (TeX in I146-rpi5/dg06_src): eqs. (3)-(4) = RW2/RW2a, (6) = poten1, (9) = RW3a, (AppB1), (RW17), (nv1), (im), (cz).
Run: python zerilli_gauge_rs.py [--part=planted|gauge|static|rs|all] [--mutate=NAME[,NAME]] [--failfast] [--skip=PL.G,PL.X] > zerilli_gauge_rs_output.txt
Mutations (targets in I146-rpi5/README.md): dual_inv ds_typo bg_nonvac leg_sign lie_drop drop_trace route2_drop ext_G solve_sign solve_drop_s l2_struct gs_a hand_typo
  gdiv_set gdiv_s l2 fam_typo st_Edep st_l2 kill_dropG flow_typo rs_typo pole_inject flow_h cz_coef cz_sign rw17_typo nv1_h1 nv1_h2 im_coef n2_sign log_inject
  rs6c_typo sflow_typo sflow_h nolog_blind ex_typed
"""
import os
import sys
import time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zerilli_derivation_lib import *

PART = 'all'
SKIP = []
for _a in sys.argv[1:]:
    if _a.startswith('--part='):
        PART = _a[len('--part='):]
    if _a.startswith('--skip='):
        SKIP = _a[len('--skip='):].split(',')
STAGES = ['planted', 'gauge', 'static', 'rs', 'all']


def skipped(cid):
    return any(cid.startswith(p) for p in SKIP)


run = Run('zerilli_gauge_rs: gauge attainability (e^{st} and static) and the flow at r_s against DG06 Lemma 1 (part=%s)' % PART)
T0 = time.time()


def lap(msg):
    print('# [%6.1f s] %s' % (time.time() - T0, msg), flush=True)


def stop_after(stage):
    if STAGES.index(PART) <= STAGES.index(stage):
        run.finish()


f = f_bg() + (M ** 2 / r ** 2 if mut('bg_nonvac') else 0)    # background f (mutation bg_nonvac: not a vacuum solution)
fq = fM()                                                     # f of the ansatz, the gauge solve and the master variables
g0 = sp.diag(-f, 1 / f, r ** 2 / (1 - x ** 2), r ** 2 * (1 - x ** 2))
Px = sp.Derivative(P, x)
x2 = 1 - x ** 2
Eest = sp.exp(s * t)                                          # e^{st}
E1 = sp.Integer(1)                                            # static
P0s, P1s = sp.symbols('P0 P1')


def red(e, E):
    """component -> divide by E, Legendre reduction: expression linear in P0, P1 with rational coefficients"""
    return tidy(leg_red(sp.expand(sp.sympify(e) / E)))


def same(e1, e2):
    return cancel(together(e1 - e2)) == 0


def zeros():
    return [[0] * N for _ in range(N)]


def sym_set(h, m, n, val):
    h[m][n] = val
    h[n][m] = val


# ======================================================================= planted checks: Ricci code (#143's)
hg = [[Function('h%d%d' % (min(i, j), max(i, j)))(t, r, x) for j in range(N)] for i in range(N)]
g_, gi_ = build_metric(g0, hg)
prod0 = sp.Matrix(N, N, lambda i, j: sum((g_[i][k] * gi_[k][j]).a0 for k in range(N)))
prod1 = sp.Matrix(N, N, lambda i, j: sum((g_[i][k] * gi_[k][j]).a1 for k in range(N)))
run.planted('PL.0a', prod0.applyfunc(tidy) == sp.eye(N) and prod1.applyfunc(tidy) == sp.zeros(N, N), 'g g^-1 != 1 at order eps')
lap('PL.0a')
a_ = Symbol('a', positive=True)
fds = 1 - r ** 2 / a_ ** 2
g_ds = sp.diag(-fds, 1 / fds, r ** 2 / (1 - x ** 2), r ** 2 * (1 - x ** 2))
R0_ds, _ = lin_ricci(g_ds, sp.zeros(N, N))
c_ds = 2 if mut('ds_typo') else 3                              # mutation ds_typo: expected Ricci typed as (2/a^2) g
run.planted('PL.0b', all(tidy(R0_ds[m][n] - c_ds / a_ ** 2 * g_ds[m, n]) == 0 for m in range(N) for n in range(N)), 'Ricci of static de Sitter is not (3/a^2) g')
lap('PL.0b')
R0_bg, _ = lin_ricci(g0, sp.zeros(N, N))
run.planted('PL.0c', all(tidy(R0_bg[m][n]) == 0 for m in range(N) for n in range(N)), 'Ricci of g_Schw(M) is not 0')
lap('PL.0c')
ok = True
for n_ in range(2, 7):
    for order in (2, 3, 4):
        red_ = leg_red(sp.Derivative(P, (x, order)))
        lhs = red_.subs({P0s: sp.legendre(n_, x), P1s: diff(sp.legendre(n_, x), x), L: n_ * (n_ + 1)})
        ok = ok and cancel(lhs - diff(sp.legendre(n_, x), x, order)) == 0
run.planted('PL.0d', ok, 'Legendre reduction differs from the explicit derivative of P_l')
lap('PL.0d')

# ----------------------------------------------------------------------- the seven functions, gauge vectors, Lie derivatives
NAMES7 = ['H0', 'H1', 'H2', 'K', 'G', 'Jt', 'Jr']
GNAMES = ['xa', 'xb', 'xc']
ENAMES = ['ea', 'eb', 'ec']
rf = {n: Function(n)(r) for n in NAMES7 + GNAMES + ENAMES}
Tstruct = x * Px - (L / 2) * P                                 # (1-x^2) T_xx; T_phiphi/(1-x^2) = -Tstruct


def build_h(F, E):
    """general even (l, m = 0) h from the seven functions F (dict name -> expression in r)"""
    h = zeros()
    sym_set(h, 0, 0, fq * F['H0'] * P * E)
    sym_set(h, 0, 1, F['H1'] * P * E)
    sym_set(h, 1, 1, F['H2'] * P * E / fq)
    sym_set(h, 0, 2, F['Jt'] * Px * E)
    sym_set(h, 1, 2, F['Jr'] * Px * E)
    sym_set(h, 2, 2, r ** 2 * E * (F['K'] * P + F['G'] * Tstruct) / x2)
    sym_set(h, 3, 3, r ** 2 * E * x2 * (F['K'] * P - F['G'] * Tstruct))
    return h


def xi_vec(a, b, c, E):
    """xi^t = a P E, xi^r = b P E, xi_x = c P' E (lower), i.e. xi^x = c (1-x^2) P' E / r^2"""
    return [a * P * E, b * P * E, c * x2 * Px * E / r ** 2, 0]


def lie_partial(g0m, xi):
    """(L_xi g)_mn = xi^a d_a g_mn + g_an d_m xi^a + g_ma d_n xi^a. Mutation lie_drop drops the last term."""
    out = zeros()
    for m in range(N):
        for n in range(N):
            tot = 0
            for a in range(N):
                tot += xi[a] * diff(g0m[m, n], COORD[a]) + g0m[a, n] * diff(xi[a], COORD[m])
                if not mut('lie_drop'):
                    tot += g0m[m, a] * diff(xi[a], COORD[n])
            out[m][n] = tot
    return out


_G0I = g0.inv().applyfunc(cancel)
_GAM = christoffel0(g0, _G0I)


def lie_cov(g0m, xi):
    """(L_xi g)_mn = nabla_m xi_n + nabla_n xi_m with xi_n = g_na xi^a and the background Christoffel symbols (independent of lie_partial)."""
    xl = [sum(g0m[a, b] * xi[b] for b in range(N)) for a in range(N)]
    Dx = [[diff(xl[n], COORD[m]) - sum(_GAM[a][m][n] * xl[a] for a in range(N)) for n in range(N)] for m in range(N)]
    return [[Dx[m][n] + Dx[n][m] for n in range(N)] for m in range(N)]


def extract(h, E):
    """the seven functions from a perturbation h (4x4 expression lists), by Legendre reduction and reading the coefficients of P, xP', P' separately.
    Returns (dict, ok, why). Inside: the odd components vanish; tt, tr, rr carry P only; tx, rx carry P' only; the angular block carries P and xP' with
    x-independent coefficients, bxx + bpp = 0 and app - axx = L bxx (the covariance identities of the traceless tensor)."""
    ok = True
    why = []
    out = {}
    for (m, n) in [(0, 3), (1, 3), (2, 3)]:
        if red(h[m][n], E) != 0:
            ok = False
            why.append('odd component (%d,%d) nonzero' % (m, n))
    for nm, (m, n), norm in (('H0', (0, 0), 1 / fq), ('H1', (0, 1), 1), ('H2', (1, 1), fq)):
        A, B, rest = ang_split(red(h[m][n], E))
        if B != 0 or rest != 0 or cancel(diff(A, x)) != 0:
            ok = False
            why.append('%s: not a pure P term with x-independent coefficient' % nm)
        out[nm] = cancel(A * norm)
    for nm, (m, n) in (('Jt', (0, 2)), ('Jr', (1, 2))):
        A, B, rest = ang_split(red(h[m][n], E))
        if A != 0 or rest != 0 or cancel(diff(B, x)) != 0:
            ok = False
            why.append('%s: not a pure P\' term with x-independent coefficient' % nm)
        out[nm] = B
    axx_, bxx_, rxx_ = ang_split(cancel(red(h[2][2], E) * x2 / r ** 2))
    app_, bpp_, rpp_ = ang_split(cancel(red(h[3][3], E) / (r ** 2 * x2)))
    bxx_, bpp_ = cancel(bxx_ / x), cancel(bpp_ / x)
    if rxx_ != 0 or rpp_ != 0 or any(cancel(diff(v, x)) != 0 for v in (axx_, bxx_, app_, bpp_)):
        ok = False
        why.append('angular block: coefficients not x-independent')
    if cancel(bxx_ + bpp_) != 0 or cancel(app_ - axx_ - L * bxx_) != 0:
        ok = False
        why.append('angular covariance identities fail')
    out['K'] = cancel((axx_ + app_) / 2)
    out['G'] = cancel(bxx_ - bpp_) if mut('ext_G') else cancel((bxx_ - bpp_) / 2)   # mutation ext_G: G twice too large
    return out, ok, '; '.join(why)


# ----------------------------------------------------------------------- jets
def jets_of_names(e, names):
    return jetify(sp.sympify(e), names)


def Dr(e, names):
    """total r-derivative of an expression in r and the jets name_n of the functions in `names`"""
    out = diff(e, r)
    for sy in e.free_symbols:
        if '_' in sy.name:
            nm, n_ = sy.name.rsplit('_', 1)
            if nm in names and n_.isdigit():
                out = out + diff(e, sy) * Symbol('%s_%d' % (nm, int(n_) + 1))
    return out


def jetD(e, k, names):
    for _ in range(k):
        e = cancel(Dr(e, names))
    return e


def jet_coeff_list(e, names):
    """expression linear in jets of `names`: dict jet symbol -> coefficient"""
    e = sp.expand(e)
    out = {}
    for sy in e.free_symbols:
        if '_' in sy.name and sy.name.rsplit('_', 1)[0] in names and sy.name.rsplit('_', 1)[1].isdigit():
            out[sy] = cancel(together(e.coeff(sy, 1)))
    return out


# ----------------------------------------------------------------------- the effect of xi on the seven functions (read off by `extract`, not typed)
_DXI = {}


def dxi(E):
    """seven functions of L_xi g for xi = (xa P E, xb P E, xc P' E): dict of jetified expressions in xa_n, xb_n, xc_n (cached by E)"""
    key = str(E)
    if key not in _DXI:
        h = lie_partial(g0, xi_vec(rf['xa'], rf['xb'], rf['xc'], E))
        d, ok_, why_ = extract(h, E)
        _DXI[key] = ({k: jets_of_names(v, GNAMES) for k, v in d.items()}, ok_, why_)
    return _DXI[key]


def apply_xi(F, D, a0, b0, c0, names):
    """F + D(xi) with xi = (a0, b0, c0) expressions in the jets of `names`; returns dict of cancelled expressions"""
    rep = {}
    for nm, e0 in (('xa', a0), ('xb', b0), ('xc', c0)):
        cur = e0
        for k in range(4):
            rep[Symbol('%s_%d' % (nm, k))] = cur
            cur = jetD(cur, 1, names)
    return {k: cancel(together(F[k] + D[k].subs(rep, simultaneous=True))) for k in NAMES7}


def solve_gauge(F, D, names, real=False):
    """(a, b, c) such that Jt, Jr, G of F + D(xi) vanish. Sequential: G fixes c, Jr fixes b (c' known), Jt fixes a. Linear in the unknown at each step."""
    xa0, xb0, xc0, xc1 = Symbol('xa_0'), Symbol('xb_0'), Symbol('xc_0'), Symbol('xc_1')
    eqG = F['G'] + D['G']
    cc = cancel(diff(eqG, xc0))
    c0 = cancel(-(eqG - cc * xc0) / cc)
    if mut('solve_sign'):
        c0 = -c0                                               # mutation solve_sign: c = +r^2 G / 2
    c1 = jetD(c0, 1, names)
    rc = {xc0: c0, xc1: c1, Symbol('xc_2'): jetD(c0, 2, names)}
    eqR = (F['Jr'] + D['Jr']).subs(rc, simultaneous=True)
    bb = cancel(diff(eqR, xb0))
    b0 = cancel(-(eqR - bb * xb0) / bb)
    eqT = F['Jt'] + D['Jt']
    if mut('solve_drop_s'):
        eqT = eqT.subs(xc0, 0)                                 # mutation solve_drop_s: the s c term of Delta Jt dropped
    eqT = eqT.subs(rc, simultaneous=True)
    aa = cancel(diff(eqT, xa0))
    a0 = cancel(-(eqT - aa * xa0) / aa)
    if real and mut('gs_a'):
        a0 = cancel(F['Jt'] / fq)                              # mutation gs_a (real data only): a = Jt/f, the s c term dropped
    return a0, b0, c0


def check_gauge_solution(F, E, names, real=False):
    D, ok_, why_ = dxi(E)
    a0, b0, c0 = solve_gauge(F, D, names, real)
    Fp = apply_xi(F, D, a0, b0, c0, names)
    return (a0, b0, c0), Fp, ok_, why_


# ----------------------------------------------------------------------- PL.L, PL.G, PL.X
xi_gen = xi_vec(rf['xa'], rf['xb'], rf['xc'], Eest)
hL_p = lie_partial(g0, xi_gen)
hL_c = lie_cov(g0, xi_gen)
bad = [(m, n) for m in range(N) for n in range(N) if cancel(together((hL_p[m][n] - hL_c[m][n]).doit())) != 0]
run.planted('PL.L', not bad, 'the partial-derivative and covariant forms of L_xi g differ at components %s' % bad)
lap('PL.L')

_LIN = {}


def lin(h, key):
    if key not in _LIN:
        _, _LIN[key] = lin_ricci(g0, sp.Matrix(h))
    return _LIN[key]


def reduce_all(R1, E):
    return [[red(R1[m][n], E) if n >= m else None for n in range(N)] for m in range(N)]


if not skipped('PL.G'):
    Rr = reduce_all(lin(hL_p, 'gauge_even'), Eest)
    run.planted('PL.G', all(Rr[m][n] == 0 for m in range(N) for n in range(m, N)),
                'pure-gauge h (even) does not give dR = 0: %s' % [(m, n) for m in range(N) for n in range(m, N) if Rr[m][n] != 0])
    lap('PL.G')

F_rw = {n: (rf[n] if n in ('H0', 'H1', 'H2', 'K') else 0) for n in NAMES7}
h_even_est = build_h(F_rw, Eest)
if not skipped('PL.X'):
    R1 = lin(h_even_est, 'even')
    R1b = route2_delta_ricci(g0, h_even_est)
    bad = [(m, n) for m in range(N) for n in range(m, N) if red(R1[m][n] - R1b[m][n], Eest) != 0]
    run.planted('PL.X', not bad, 'route 1 (dual) and route 2 (covariant) differ at components %s' % bad)
    lap('PL.X')

# ----------------------------------------------------------------------- PL.E extraction round trip
okE = True
whyE = []
for nm_, E_ in (('e^{st}', Eest), ('static', E1)):
    d_, ok_, why_ = extract(build_h({n: rf[n] for n in NAMES7}, E_), E_)
    okE = okE and ok_ and all(cancel(d_[n] - rf[n]) == 0 for n in NAMES7)
    whyE.append('%s: ok=%s %s %s' % (nm_, ok_, why_, [n for n in NAMES7 if cancel(d_[n] - rf[n]) != 0]))
run.planted('PL.E', okE, 'extraction of the seven functions from h built from them does not return them: %s' % whyE)
lap('PL.E')


# ----------------------------------------------------------------------- PL.P planted pure gauge through the same solver
def planted_pure_gauge(E):
    eta = xi_vec(rf['ea'], rf['eb'], rf['ec'], E)
    hpure = lie_cov(g0, eta)
    names = NAMES7 + ENAMES
    Fpure, ok_, why_ = extract(hpure, E)
    if not ok_:
        return False, 'extract(h_pure): ' + why_
    Fj = {k: jets_of_names(v, ENAMES) for k, v in Fpure.items()}
    (a0, b0, c0), Fp, ok2, why2 = check_gauge_solution(Fj, E, names)
    res = [cancel(a0 + Symbol('ea_0')), cancel(b0 + Symbol('eb_0')), cancel(c0 + Symbol('ec_0'))]
    allzero = all(v == 0 for v in Fp.values())
    # (ii) RW part + pure gauge
    Frw = {k: jets_of_names(Fpure[k], ENAMES) + (jets_of_names(rf[k], NAMES7) if k in ('H0', 'H1', 'H2', 'K') else 0) for k in NAMES7}
    (a1, b1, c1), Fp2, _, _ = check_gauge_solution(Frw, E, names)
    rwback = all(cancel(Fp2[k] - jets_of_names(rf[k], NAMES7)) == 0 for k in ('H0', 'H1', 'H2', 'K')) and all(Fp2[k] == 0 for k in ('G', 'Jt', 'Jr'))
    res2 = [cancel(a1 + Symbol('ea_0')), cancel(b1 + Symbol('eb_0')), cancel(c1 + Symbol('ec_0'))]
    ok_all = ok2 and allzero and all(v == 0 for v in res) and rwback and all(v == 0 for v in res2)
    return ok_all, 'xi + eta = %s, nonzero transformed functions %s, RW part returned %s, xi + eta (ii) = %s' % (res, [k for k, v in Fp.items() if v != 0], rwback, res2)


for nm_, E_ in (('e^{st}', Eest), ('static', E1)):
    ok_, msg_ = planted_pure_gauge(E_)
    run.planted('PL.P' if nm_ == 'e^{st}' else 'PL.P.static', ok_, '%s: %s' % (nm_, msg_))
    lap('PL.P ' + nm_)

# ----------------------------------------------------------------------- PL.R planted: the residue / no-logarithm routine of R.10 on flows with a known answer
def laurent_test(Rm, F1):
    """flow u' = (Rm/(r - r_s) + F1 + ...) u. Returns (eigenvalues of Rm are {0, 1}, no logarithm in the exponent-0 solution).
    Exponent-0 solution u = u0 + u1 (r - r_s) + ...: R u0 = 0 and (I - R) u1 = F1 u0, solvable iff w^T F1 u0 = 0 for w^T (I - R) = 0."""
    if mut('nolog_blind'):
        F1 = sp.zeros(2, 2)                                      # mutation nolog_blind: the next Laurent coefficient is ignored
    eig = cancel(Rm.trace() - 1) == 0 and cancel(Rm.det()) == 0 and Rm != sp.zeros(2, 2)
    if not eig:
        return False, False                                      # no resonance between the exponents 0 and 1: nothing to test
    u0_ = Rm.nullspace()[0]
    w_ = (Rm - sp.eye(2)).T.nullspace()[0]
    return eig, cancel((w_.T * F1 * u0_)[0, 0]) == 0


xx_ = sp.Symbol('xx_', positive=True)
uA, uB = sp.Function('uA')(xx_), sp.Function('uB')(xx_)
R_pl = sp.Matrix([[0, 0], [0, 1]])
eigA, nologA = laurent_test(R_pl, sp.zeros(2, 2))                 # u' = diag(0, 1) u / x: u = (c1, c2 x), analytic
eigB, nologB = laurent_test(R_pl, sp.Matrix([[0, 0], [1, 0]]))    # u1' = 0 ..., u2' = u2/x + u1: u2 = u1 x ln x + c x, a logarithm
eigC, _ = laurent_test(sp.Matrix([[0, 0], [0, 2]]), sp.zeros(2, 2))   # eigenvalues {0, 2}
solA = sp.dsolve(sp.Eq(uA.diff(xx_), uA / xx_), uA)
solB = sp.dsolve(sp.Eq(uB.diff(xx_), uB / xx_ + 1), uB)         # the logarithm, from sympy's own ODE solver
run.planted('PL.R', eigA and nologA and eigB and (not nologB) and (not eigC) and (not solA.rhs.has(sp.log)) and solB.rhs.has(sp.log),
            'planted flows: analytic flow reported eig %s nolog %s; flow with a logarithm reported nolog %s (dsolve has log: %s); eigenvalues {0,2} reported eig %s' % (
                eigA, nologA, nologB, solB.rhs.has(sp.log), eigC))
lap('PL.R')
stop_after('planted')

# ======================================================================= gauge, e^{st}
# ---- GA.1 the structure function T_l = x P_l' - (L/2) P_l
c_, Ls_ = sp.symbols('c_ Ls_')
dep = sp.solve([c_ * (c_ - 1), Ls_ - c_ * (c_ + 1)], [c_, Ls_], dict=True)
dep_L = sorted(set(d_[Ls_] for d_ in dep))                    # values of L at which P and xP' can be dependent: [0, 2]
lmin = 1 if mut('l2_struct') else 2                            # mutation l2_struct: GA.1 is run from l = 1
okT = dep_L == [0, 2]
for lv in range(lmin, 31):
    Tl = sp.Poly(sp.expand(x * diff(sp.legendre(lv, x), x) - sp.Rational(lv * (lv + 1), 2) * sp.legendre(lv, x)), x)
    kl = sp.Poly(sp.legendre(lv, x), x).LC()
    okT = okT and not Tl.is_zero and Tl.LC() == -kl * lv * (lv - 1) / 2 and lv * (lv - 1) != 0
okT = okT and all(sp.Poly(sp.expand(x * diff(sp.legendre(lv, x), x) - sp.Rational(lv * (lv + 1), 2) * sp.legendre(lv, x)), x).is_zero for lv in (0, 1))
run.check('GA.1', okT, 'T_l = xP_l\' - (L/2)P_l: dependence of P, xP\' only at L in %s; T_l nonzero with leading coefficient -k_l l(l-1)/2 for l = %d..30: %s; T_0 = T_1 = 0' % (dep_L, lmin, okT))
lap('GA.1')

# ---- GS.1 the solve, e^{st}
names_g = NAMES7
Fgen = {n: jets_of_names(rf[n], NAMES7) for n in NAMES7}
(a_s, b_s, c_s), Fp_s, ok_ex, why_ex = check_gauge_solution(Fgen, Eest, names_g, real=True)
run.check('GS.1', ok_ex and Fp_s['Jt'] == 0 and Fp_s['Jr'] == 0 and Fp_s['G'] == 0,
          'solve: extract ok = %s (%s); transformed Jt, Jr, G = %s, %s, %s' % (ok_ex, why_ex, Fp_s['Jt'], Fp_s['Jr'], Fp_s['G']))
Dest = dxi(Eest)[0]
run.note('GS.1: L_xi g, e^{st}: the effect of xi = (a, b, c) on the seven functions, read off by extract (jets xa_n = d^n a/dr^n, ...):')
for k_ in ('Jt', 'Jr', 'G', 'K', 'H0', 'H1', 'H2'):
    run.note('  Delta %s = %s' % (k_, sp.factor_terms(Dest[k_])))
run.note('GS.1: a = %s' % sp.factor_terms(a_s))
run.note('GS.1: b = %s' % sp.factor_terms(b_s))
run.note('GS.1: c = %s' % sp.factor_terms(c_s))
for k_ in ('H0', 'H1', 'H2', 'K'):
    run.note('GS.1: %s\' = %s' % (k_, sp.factor_terms(Fp_s[k_])))
lap('GS.1')

# ---- GS.2 closed forms (typed comparison targets)
G0j, G1j, Jt0, Jr0 = (Symbol('G_0'), Symbol('G_1'), Symbol('Jt_0'), Symbol('Jr_0'))
f_ = fq
fp_ = diff(fq, r)
c_t = -r ** 2 * G0j / 2
b_t = -f_ * (Jr0 - r ** 2 * G1j / 2)
a_t = (Jt0 + s * c_t) / f_
Lk = L / 4 if mut('hand_typo') else L / 2                      # mutation hand_typo: (L/2) G -> (L/4) G in the typed K'
Kp_t = Symbol('K_0') - (2 * f_ / r) * (Jr0 - r ** 2 * G1j / 2) + Lk * G0j      # K' = K + 2b/r - L c/r^2
H0p_t = Symbol('H0_0') - (fp_ / f_) * b_t - 2 * s * a_t
a1_t = jetD(a_t, 1, NAMES7)
b1_t = jetD(b_t, 1, NAMES7)
H1p_t = Symbol('H1_0') - f_ * a1_t + s * b_t / f_
H2p_t = Symbol('H2_0') - (fp_ / f_) * b_t + 2 * b1_t
okS2 = (same(a_s, a_t) and same(b_s, b_t) and same(c_s, c_t) and same(Fp_s['K'], Kp_t) and same(Fp_s['H0'], H0p_t)
        and same(Fp_s['H1'], H1p_t) and same(Fp_s['H2'], H2p_t))
run.check('GS.2', okS2, 'closed forms differ: a %s, b %s, c %s, K %s, H0 %s, H1 %s, H2 %s' % tuple(
    same(u_, v_) for u_, v_ in ((a_s, a_t), (b_s, b_t), (c_s, c_t), (Fp_s['K'], Kp_t), (Fp_s['H0'], H0p_t), (Fp_s['H1'], H1p_t), (Fp_s['H2'], H2p_t))))
lap('GS.2')


# ---- GS.3 divisors
def den_factors(e):
    """irreducible factors (with multiplicity) of the denominator of the coefficients of the jets in the expression e"""
    den = 1
    for sy, co in jet_coeff_list(e, NAMES7).items():
        den = sp.lcm(den, sp.fraction(cancel(together(co)))[1])
    return sp.factor_list(den)[1]


VOCAB = [s, r, L, L - 2, M] if mut('gdiv_set') else [s, r, r - 2 * M, L, L - 2, M]    # mutation gdiv_set: vocabulary without r - 2M


def classify(fac):
    fac = sp.expand(fac)
    for v_ in VOCAB:
        if cancel(fac - sp.expand(v_)) == 0 or cancel(fac + sp.expand(v_)) == 0:
            return str(v_)
    return 'UNEXPECTED(%s)' % fac


met = set()
for nm_, e_ in (('a', a_s), ('b', b_s), ('c', c_s), ('H0\'', Fp_s['H0']), ('H1\'', Fp_s['H1']), ('H2\'', Fp_s['H2']), ('K\'', Fp_s['K'])):
    for fac_, mult_ in den_factors(e_):
        met.add(classify(fac_))
        run.note('GS.3: %s | denominator factor %s (multiplicity %d)' % (nm_, fac_, mult_))
if mut('gdiv_s'):
    met.add('s')                                               # mutation gdiv_s: s appended to the divisors met
okS3 = met == {str(r), str(r - 2 * M)} and not any(m_.startswith('UNEXPECTED') for m_ in met)
run.check('GS.3', okS3, 'divisors met %s; expected exactly {%s, %s} (no s, no L - 2, no M)' % (sorted(met), r, r - 2 * M))
run.note('GS.3: the solve divides by %s only; constants aside' % sorted(met))
lap('GS.3')


# ---- GU.1, GU.2 explicit l: uniqueness and the L = 2 failure
def homogeneous_explicit(lv, E):
    """the homogeneous problem at explicit l (P = P_l): h = 0, L_xi g must have Jt' = Jr' = 0 and zero traceless angular part as a function of x.
    Returns (number of free functions, dict of the relations). The relations are read from the components of L_xi g, not typed."""
    Pl = sp.legendre(lv, x)
    xi = [rf['xa'] * Pl * E, rf['xb'] * Pl * E, rf['xc'] * x2 * diff(Pl, x) * E / r ** 2, 0]
    h = lie_partial(g0, xi)
    comp = lambda m, n: cancel(together(sp.expand((h[m][n] / E).doit())))
    tx = cancel(comp(0, 2) / diff(Pl, x))
    rx = cancel(comp(1, 2) / diff(Pl, x))
    Dtl = cancel(comp(2, 2) * x2 / r ** 2 - comp(3, 3) / (r ** 2 * x2))   # = 2 G T = 4 c T_l / r^2 for the pure gauge
    coeffs = sp.Poly(sp.expand(Dtl), x).all_coeffs() if Dtl != 0 else []
    forced = any(cancel(co) != 0 for co in coeffs)
    return (0 if forced else 1), {'tx': tx, 'rx': rx, 'Dtl': Dtl}


lmin_u = 1 if mut('l2') else 2                                # mutation l2: GU.1 is run from l = 1
free_u = {lv: homogeneous_explicit(lv, Eest)[0] for lv in range(lmin_u, 9)}
run.check('GU.1', all(v == 0 for v in free_u.values()), 'free functions in the homogeneous solve at l = %s: %s (0 = only xi = 0)' % (sorted(free_u), free_u))
lap('GU.1')
nfree1, rel1 = homogeneous_explicit(1, Eest)
cf = Function('cfam')(r)
a_fam = (2 if mut('fam_typo') else 1) * s * cf / fq            # mutation fam_typo: a = 2 s c / f
b_fam = -fq * (diff(cf, r) - 2 * cf / r)
okU2 = (nfree1 == 1 and rel1['Dtl'] == 0
        and cancel(rel1['tx'].subs({rf['xa']: a_fam, rf['xc']: cf, rf['xb']: b_fam}).doit()) == 0
        and cancel(rel1['rx'].subs({rf['xa']: a_fam, rf['xc']: cf, rf['xb']: b_fam}).doit()) == 0)
run.check('GU.2', okU2, 'at l = 1: free functions %d, traceless part %s, family a = s c/f, b = -f(c\' - 2c/r) solves Jt\' = Jr\' = 0: %s' % (
    nfree1, rel1['Dtl'], okU2))
run.note('GU.2: at l = 1 (L = 2) the traceless angular part of L_xi g vanishes for every xi, so Jt\' = Jr\' = 0 leaves one free function c(r): '
         'a = s c/f, b = -f (c\' - 2c/r). The RW gauge (Jt = Jr = G = 0) does not fix xi there.')
lap('GU.2')
stop_after('gauge')

# ======================================================================= gauge, static (no time dependence)
Fgen_st = {n: jets_of_names(rf[n], NAMES7) for n in NAMES7}
E_st = Eest if mut('st_Edep') else E1                          # mutation st_Edep: the "static" ansatz carries e^{st}
(a_t0, b_t0, c_t0), Fp_t0, ok_st, why_st = check_gauge_solution(Fgen_st, E_st, NAMES7)
Dst = dxi(E_st)[0]
sfree = all(not any(sy.name == 's' for sy in v.free_symbols) for v in Dst.values())
agree = same(a_t0, a_s.subs(s, 0)) and same(b_t0, b_s.subs(s, 0)) and same(c_t0, c_s.subs(s, 0))
agree_t = all(same(Fp_t0[k], Fp_s[k].subs(s, 0)) for k in ('H0', 'H1', 'H2', 'K'))
run.check('GST.1', ok_st and sfree and agree and agree_t and Fp_t0['Jt'] == 0 and Fp_t0['Jr'] == 0 and Fp_t0['G'] == 0,
          'static solve: extract ok %s; no s in the effect of xi %s; equals the s = 0 value of GS.1: xi %s, RW functions %s' % (ok_st, sfree, agree, agree_t))
run.note('GST.1: static: a = %s, b = %s, c = %s (algebraic in G, G\', Jt, Jr: no ODE in r, no integration constant)' % (
    sp.factor_terms(a_t0), sp.factor_terms(b_t0), sp.factor_terms(c_t0)))
for k_ in ('H0', 'H1', 'H2', 'K'):
    run.note('GST.1: static %s\' = %s' % (k_, sp.factor_terms(Fp_t0[k_])))
lap('GST.1')
lmin_s = 1 if mut('st_l2') else 2                              # mutation st_l2: GST.2 is run from l = 1
free_s = {lv: homogeneous_explicit(lv, E1)[0] for lv in range(lmin_s, 9)}
run.check('GST.2', all(v == 0 for v in free_s.values()), 'free functions in the static homogeneous solve at l = %s: %s' % (sorted(free_s), free_s))
lap('GST.2')
# GST.3: xi = t kappa: d_t(L_xi g) = L_kappa g must vanish; the algebraic equations of L_kappa g = 0 (Delta Jt, Delta G, Delta K) have a nonzero determinant
Dk, ok_k, why_k = dxi(E1)
alg = [Dk['Jt'], Dk['G'], Dk['K']]
if mut('kill_dropG'):
    alg = [Dk['Jt'], Dk['K'], sp.Integer(0)]                   # mutation kill_dropG: the G equation dropped
unk = [Symbol('xa_0'), Symbol('xb_0'), Symbol('xc_0')]
Malg = sp.Matrix([[diff(e_, u_) for u_ in unk] for e_ in alg]).applyfunc(cancel)
detA = cancel(Malg.det())
run.check('GST.3', ok_k and detA != 0, 'the algebraic equations of L_kappa g = 0 (Jt, G, K) for kappa = (alpha P, beta P, gamma P\') have determinant %s' % sp.factor(detA))
run.note('GST.3: det = %s, nonzero for r > 0, f != 0, and (through GA.1) l >= 2: kappa = 0, so xi = t kappa adds no residual freedom' % sp.factor(detA))
lap('GST.3')
stop_after('static')

# ======================================================================= r_s: the flow of #143, zeta, DG06's Lemma 1
# ---- the even radial equations of the e^{st} ansatz (adapted from #143's A.2, E.1, E.2: derived here, not typed)
ev_names = ['H0', 'H1', 'H2', 'K']
Rr_ev = reduce_all(lin(h_even_est, 'even'), Eest)
SPEC = (('tt', (0, 0), 1, ['P0']), ('tr', (0, 1), 1, ['P0']), ('rr', (1, 1), 1, ['P0']), ('tx', (0, 2), 1, ['P1']), ('rx', (1, 2), 1, ['P1']),
        ('xx', (2, 2), x2, ['P0', 'xP1']), ('pp', (3, 3), 1 / x2, ['P0', 'xP1']))
ev = {}
for nm, comp, norm, basis in SPEC:
    d_, ok_, why_ = coef_eqs(Rr_ev[comp[0]][comp[1]], norm, basis, nm)
    if not ok_:
        run.planted('R.0', False, 'angular form of the even component %s: %s' % (nm, why_))
    for k_, v_ in d_.items():
        ev[nm + '.' + k_] = jetify(v_, ev_names)
Ett, Etr, Err, Etx, Erx = ev['tt.P0'], ev['tr.P0'], ev['rr.P0'], ev['tx.P1'], ev['rx.P1']
Ec2 = ev['xx.xP1']
EtrP = cancel(ev['xx.P0'] + ev['pp.P0'])
H0_, H2_ = jet_sym('H0', 0), jet_sym('H2', 0)
c2_coef = cancel(diff(Ec2, H0_))
ok_trl = cancel(Ec2 - c2_coef * (H0_ - H2_)) == 0 and c2_coef != 0
Hs = {}
for n_ in range(4):
    Hs[jet_sym('H0', n_)] = jet_sym('H', n_)
    Hs[jet_sym('H2', n_)] = jet_sym('H', n_)
eqsH = {nm: e_.subs(Hs, simultaneous=True) for nm, e_ in (('tt', Ett), ('tr', Etr), ('rr', Err), ('tx', Etx), ('rx', Erx), ('trace', EtrP))}
K0, K1, K2 = (jet_sym('K', n_) for n_ in range(3))
H10, H11, H12 = (jet_sym('H1', n_) for n_ in range(3))
Hh0, Hh1, Hh2 = (jet_sym('H', n_) for n_ in range(3))
first = ('tr', 'tx', 'rx')
unk = [K1, H11, Hh1]
Mat3 = sp.Matrix([[diff(eqsH[nm], v) for v in unk] for nm in first]).applyfunc(cancel)
det3 = cancel(Mat3.det())
rhs3 = sp.Matrix([-(eqsH[nm] - sum(diff(eqsH[nm], v) * v for v in unk)) for nm in first]).applyfunc(cancel)
sol3 = Mat3.LUsolve(rhs3).applyfunc(cancel)
FK, FH1, FH = sol3[0], sol3[1], sol3[2]
states3 = [K0, H10, Hh0]
F3 = [FK, FH1, FH]
vm3 = {'K': K0, 'H1': H10, 'H': Hh0}
Grr = cancel(together((Err / 2 + Ett / (2 * fq ** 2) - EtrP / (2 * fq * r ** 2)).subs(Hs, simultaneous=True)))    # delta G_rr
Aflow = cancel(flow_sub(Grr, vm3, states3, F3))
aK, a1, aH = (cancel(diff(Aflow, v)) for v in states3)
lam_ = lam                                                      # (L - 2)/2
mu_ = -2 * M
rs = -3 * M / (lam_ + 1 if mut('rs_typo') else lam_)           # mutation rs_typo: r_s typed -3M/(lambda+1)


def D3(e, F=None):
    F = F3 if F is None else F
    return diff(e, r) + sum(diff(e, states3[i]) * F[i] for i in range(3))


# ---- R.1 the flow against #143's printed one
fK_t = -(r - (2 if mut('flow_typo') else 3) * M) / (r * (r - 2 * M))   # mutation flow_typo: typed (r - 3M) -> (r - 2M) in K'
FK_t = fK_t * K0 + L / (2 * r ** 2 * s) * H10 + Hh0 / r
FH1_t = r * s / (r - 2 * M) * K0 - 2 * M / (r * (r - 2 * M)) * H10 + r * s / (r - 2 * M) * Hh0
FH_t = -(r - 3 * M) / (r * (r - 2 * M)) * K0 + (-2 * L * M + L * r + 2 * r ** 3 * s ** 2) / (2 * r ** 2 * s * (r - 2 * M)) * H10 + (r - 4 * M) / (r * (r - 2 * M)) * Hh0
aK_t = -(-2 * L * M * r + L * r ** 2 - 6 * M ** 2 + 6 * M * r + 2 * r ** 4 * s ** 2 - 2 * r ** 2) / (2 * r ** 2 * (r - 2 * M) ** 2)
a1_t = (L * M + 2 * r ** 3 * s ** 2) / (2 * r ** 3 * s * (r - 2 * M))
aH_t = (L * r + 6 * M - 2 * r) / (2 * r ** 2 * (r - 2 * M))
det_t = -s * (r - 2 * M) / (4 * r)
okR1 = (ok_trl and same(FK, FK_t) and same(FH1, FH1_t) and same(FH, FH_t) and same(aK, aK_t) and same(a1, a1_t) and same(aH, aH_t) and same(det3, det_t))
run.check('R.1', okR1, 'derived flow differs from #143\'s printed: traceless ok %s, K\' %s, H1\' %s, H\' %s, aK %s, a1 %s, aH %s, det %s' % (
    ok_trl, same(FK, FK_t), same(FH1, FH1_t), same(FH, FH_t), same(aK, aK_t), same(a1, a1_t), same(aH, aH_t), same(det3, det_t)))
run.note('R.1: flow (K\', H1\', H\') = F (K, H1, H); A = aK K + a1 H1 + aH H; det of the coefficient matrix = %s' % sp.factor(det3))
for nm_, e_ in (('K\'', FK), ('H1\'', FH1), ('H\'', FH)):
    run.note('R.1:   %s = %s' % (nm_, ' + '.join('(%s) %s' % (sp.factor(diff(e_, v_)), v_) for v_ in states3)))
run.note('R.1: aK = %s, a1 = %s, aH = %s' % (sp.factor(aK), sp.factor(a1), sp.factor(aH)))
lap('R.1')

# ---- R.2 r_s is a regular point of the flow
VOC_R = [s, r, r - 2 * M]


def classify_r(fac, voc):
    fac = sp.expand(fac)
    for v_ in voc:
        if cancel(fac - sp.expand(v_)) == 0 or cancel(fac + sp.expand(v_)) == 0:
            return str(v_)
    return 'UNEXPECTED(%s)' % fac


met_r = set()
dens_ok = True
coefs_all = [cancel(diff(F3[i], states3[j])) for i in range(3) for j in range(3)] + [aK, a1, aH]
for co in coefs_all:
    if co == 0:
        continue
    for fac_, mult_ in denominators_factors(co):
        met_r.add(classify_r(fac_, VOC_R))
        if cancel(fac_.subs(r, rs)) == 0:
            dens_ok = False                                        # a denominator factor vanishes at r_s
rs_nonzero = cancel(rs) != 0 and cancel(rs - 2 * M) != 0
aH_rs = cancel(aH.subs(r, rs))
okR2 = (not any(m_.startswith('UNEXPECTED') for m_ in met_r)) and dens_ok and rs_nonzero and aH_rs == 0
run.check('R.2', okR2, 'denominator factors met %s (vocabulary {s, r, r-2M}); none vanishes at r_s = %s: %s; r_s, r_s - 2M nonzero: %s; aH(r_s) = %s' % (
    sorted(met_r), rs, dens_ok, rs_nonzero, aH_rs))
run.note('R.2: r_s = %s; r_s - 2M = %s; aH(r_s) = %s' % (sp.factor(rs), sp.factor(cancel(rs - 2 * M)), aH_rs))
lap('R.2')

# ---- R.3 the Taylor coefficients u^(k)(r_s) of a solution are finite (no pole of the flow at r_s)
F3j = list(F3)
if mut('pole_inject'):
    F3j[0] = FK + Hh0 / (lam_ * r + 3 * M)                       # mutation pole_inject: a pole at r_s in K' (as psi has)


def finite_at(expr, rpt):
    """True iff the denominator of every coefficient of the states in expr is nonzero at r = rpt"""
    ok = True
    for v_ in states3:
        num, den = sp.fraction(cancel(together(diff(expr, v_))))
        ok = ok and cancel(den.subs(r, rpt)) != 0
    return ok


cur = list(states3)
fin = []
for k_ in range(5):
    fin.append(all(finite_at(c_, rs) for c_ in cur))
    cur = [cancel(D3(c_, F3j)) for c_ in cur]
run.check('R.3', all(fin), 'u^(k)(r_s) finite for k = 0..4: %s' % fin)
run.note('R.3: the Taylor coefficients u^(k)(r_s), k = 0..4, are linear in (K, H1, H)(r_s) with finite coefficients; the flow has analytic coefficients at r_s, so every solution is analytic there')
lap('R.3')

# ---- R.4 zeta = (im)(K, K', H2 = H) on the flow, the matrix T, the range of T(r_s)
F3x = list(F3)
if mut('flow_h'):
    F3x[0] = FK - Hh0 / (2 * r)                                   # mutation flow_h: the coefficient of H in K' halved (after R.3)
imc_true = lam_ + 1
zeta = cancel((-2 * r ** 2 * (r + mu_) * F3x[0] + (2 * lam_ * r - 3 * mu_) * r * K0 + 2 * r * (r + mu_) * Hh0) / imc_true)   # (im), H2 = H
z1 = cancel(D3(zeta, F3x))
zeta_noH = cancel(diff(zeta, Hh0)) == 0
z1_noH = cancel(diff(z1, Hh0)) == 0
Tm = sp.Matrix([[diff(zeta, K0), diff(zeta, H10)], [diff(z1, K0), diff(z1, H10)]]).applyfunc(cancel)
detT = cancel(Tm.det())
detT_t = 4 * r * (lam_ * r + 3 * M) / s
Trs = Tm.subs(r, rs).applyfunc(cancel)
rank1 = cancel(Trs.det()) == 0 and Trs != sp.zeros(2, 2)
cstar = cancel(Trs[1, 0] / Trs[0, 0])
rowrel = [cancel(Trs[1, j] - cstar * Trs[0, j]) for j in range(2)]
cstar_cl = lam_ ** 2 / (3 * M * (2 * lam_ + 3))
okR4 = zeta_noH and z1_noH and same(detT, detT_t) and rank1 and rowrel == [0, 0] and same(cstar, cstar_cl) and same(cstar, M / (rs * (rs - 2 * M)))
run.check('R.4', okR4, 'zeta H-free %s, zeta\' H-free %s, det T = 4r(lam r+3M)/s %s, T(r_s) rank 1 %s, row2 = c* row1 %s, c* = lam^2/(3M(2 lam+3)) = M/(r_s(r_s-2M)) %s' % (
    zeta_noH, z1_noH, same(detT, detT_t), rank1, rowrel, same(cstar, cstar_cl)))
run.note('R.4: zeta = %s' % ' + '.join('(%s) %s' % (sp.factor(diff(zeta, v_)), v_) for v_ in states3 if diff(zeta, v_) != 0))
run.note('R.4: zeta\' = %s' % ' + '.join('(%s) %s' % (sp.factor(diff(z1, v_)), v_) for v_ in states3 if diff(z1, v_) != 0))
run.note('R.4: T = %s ; det T = %s' % (Tm.applyfunc(sp.factor).tolist(), sp.factor(detT)))
run.note('R.4: c* = %s (derived as T21/T11 at r_s); at M = -1: %s' % (sp.factor(cstar), sp.factor(cancel(cstar.subs(M, -1)))))
run.note('R.4: zeta\'(r_s) = c* zeta(r_s) holds for every solution of the 3-dimensional flow, whether or not A = 0')
lap('R.4')

# ---- R.5 DG06's (cz) as typed against the derived condition
cz_num = 3 * lam_ ** 2 if mut('cz_coef') else 2 * lam_ ** 2          # mutation cz_coef: 2 lam^2 -> 3 lam^2 in the typed (cz)
cz_c = (-1 if mut('cz_sign') else 1) * cz_num / (3 * mu_ * (2 * lam_ + 3))   # mutation cz_sign: typed with the derived sign
zrs = [cancel(diff(zeta, v_).subs(r, rs)) for v_ in states3]
arow = [cancel(diff(Aflow, v_).subs(r, rs)) for v_ in states3]
cuts = cancel(zrs[0] * arow[1] - zrs[1] * arow[0]) != 0
okR5 = same(cstar, -cz_c) and cancel(cz_c) != 0 and same(cstar.subs(M, -1), -cz_c.subs(M, -1)) and cuts
run.check('R.5', okR5, 'derived c* = -(typed (cz) constant): %s; at M = -1: %s; zeta(r_s) is not identically 0 on S(r_s) = {aK K + a1 H1 = 0}: %s' % (
    same(cstar, -cz_c), same(cstar.subs(M, -1), -cz_c.subs(M, -1)), cuts))
run.note('R.5: DISCREPANCY (arXiv TeX of DG06 against this derivation): (cz) as typed is zeta\' = +2 lam^2/(3 mu (2 lam+3)) zeta at r = 3 mu/(2 lam); the flow gives '
         'zeta\' = c* zeta with c* = %s, the opposite sign. On S the typed condition is zeta(r_s) = 0.' % sp.factor(cstar))
lap('R.5')

# ---- R.6 the ODE for zeta on S against (RW17), and its singular-point relation
Hexp = cancel(-(aK * K0 + a1 * H10) / aH)
Ti = Tm.inv().applyfunc(cancel)
Z0, Z1, Z2 = sp.symbols('Z0 Z1 Z2')
KH1 = (Ti * sp.Matrix([Z0, Z1])).applyfunc(cancel)
z2 = cancel(D3(z1, F3x))
z2S = cancel(z2.subs(Hh0, Hexp).subs({K0: KH1[0], H10: KH1[1]}, simultaneous=True))
alpha = cancel(diff(z2S, Z1))
beta = cancel(diff(z2S, Z0))
c6 = 5 if mut('rw17_typo') else 6                                   # mutation rw17_typo: 6 lam mu r -> 5 lam mu r in the typed (RW17)
a_O = (4 * lam_ * r ** 2 + c6 * lam_ * mu_ * r - 3 * mu_ ** 2) * (r + mu_) / (r ** 3 * (2 * lam_ * r - 3 * mu_))
b_O = (r + mu_) * (3 * mu_ ** 2 - 6 * lam_ * mu_ * r + 4 * lam_ ** 2 * r ** 2) / (r ** 4 * (2 * lam_ * r - 3 * mu_))
alpha_t = r ** 2 / (r + mu_) ** 2 * a_O
beta_t = r ** 2 / (r + mu_) ** 2 * (b_O + s ** 2)
al_res = cancel(cancel(alpha * (L * r + 6 * M - 2 * r)) / (L - 2)).subs(r, rs)
be_res = cancel(cancel(beta * (L * r + 6 * M - 2 * r)) / (L - 2)).subs(r, rs)
okR6 = same(alpha, alpha_t) and same(beta, beta_t) and same(al_res * cstar + be_res, 0)
run.check('R.6', okR6, 'zeta\'\' = alpha zeta\' + beta zeta from the flow equals (RW17) with d_t^2 -> s^2: alpha %s, beta %s; residue relation al c* + be = 0: %s' % (
    same(alpha, alpha_t), same(beta, beta_t), same(al_res * cstar + be_res, 0)))
run.note('R.6: alpha = %s' % sp.factor(alpha))
run.note('R.6: beta = %s' % sp.factor(beta))
run.note('R.6: the residues at r_s are al = %s, be = %s; al c* + be = %s' % (sp.factor(cancel(al_res)), sp.factor(cancel(be_res)), sp.factor(cancel(al_res * cstar + be_res))))
lap('R.6')

# ---- R.7 (nv1) against the flow and against (3)-(4)
z0s, z1s, z2s = sp.symbols('z0 z1 z2')
Kd = cancel(KH1[0].subs({Z0: z0s, Z1: z1s}))
H1d = cancel(KH1[1].subs({Z0: z0s, Z1: z1s}))
Hd = cancel(Hexp.subs({K0: Kd, H10: H1d}, simultaneous=True))
K_nv = ((mu_ + r) * r * z1s + (lam_ * r - mu_) * z0s) / (r ** 2 * (2 * lam_ * r - 3 * mu_))
H1_den = 2 if mut('nv1_h1') else r                                  # mutation nv1_h1: the typed H1 given the corrected denominator
H1_nv = s * (2 * r * (mu_ + r) * z1s + mu_ * z0s) / (H1_den * (r + mu_) * (2 * lam_ * r - 3 * mu_))
sg = -1 if mut('nv1_h2') else 1                                     # mutation nv1_h2: the typed H2 given the corrected sign of the zeta term
A2 = (3 * mu_ ** 2 - 6 * lam_ * mu_ * r + 4 * lam_ ** 2 * r ** 2) / (2 * r ** 2 * (2 * lam_ * r - 3 * mu_) ** 2)
H2_nv = ((r + mu_) / (2 * lam_ * r - 3 * mu_) * z2s - (4 * lam_ * r ** 2 + 6 * lam_ * mu_ * r - 3 * mu_ ** 2) / (2 * r * (2 * lam_ * r - 3 * mu_) ** 2) * z1s
         + sg * A2 * z0s)
H2_nv_onshell = H2_nv.subs(z2s, alpha * z1s + beta * z0s)
# DG06 (3)-(4) on psi = zeta/(2 lam r - 3 mu), off-shell (zeta a free function)
fD = (lam_ * (lam_ + 1) * r ** 2 + 3 * lam_ * M * r + 6 * M ** 2) / (r ** 2 * (lam_ * r + 3 * M))
hD = (lam_ * r ** 2 - 3 * lam_ * r * M - 3 * M ** 2) / ((r - 2 * M) * (lam_ * r + 3 * M))
Zf = Function('Zf')(r)
psiZ = Zf / (2 * lam_ * r - 3 * mu_)
K3 = fD * psiZ + (1 - 2 * M / r) * diff(psiZ, r)
H13 = s * (hD * psiZ + r * diff(psiZ, r))
H23 = diff((1 - 2 * M / r) * (hD * psiZ + r * diff(psiZ, r)), r) - K3


def zj(e):
    return e.doit().subs(sp.Derivative(Zf, (r, 2)), z2s).subs(sp.Derivative(Zf, r), z1s).subs(Zf, z0s)


K3, H13, H23 = zj(K3), zj(H13), zj(H23)
rel_K = same(K_nv, Kd)
rel_H1 = same(H1_nv, 2 / r * H1d)
rel_H2 = same(H2_nv_onshell - 2 * A2 * z0s, Hd)
rel3_K = same(K_nv, K3)
rel3_H1 = same(H1_nv, 2 / r * H13)
rel3_H2 = same(H2_nv - 2 * A2 * z0s, H23)
der3 = same(Kd, K3) and same(H1d, H13) and same(Hd, H23.subs(z2s, alpha * z1s + beta * z0s))
okR7 = rel_K and rel_H1 and rel_H2 and rel3_K and rel3_H1 and rel3_H2 and der3
run.check('R.7', okR7, 'typed (nv1) vs flow: K equal %s, H1 = (2/r) H1_derived %s, H2 - 2*(zeta term) = H_derived %s; vs (3)-(4): K %s, H1 %s, H2 %s; (3)-(4) = flow %s' % (
    rel_K, rel_H1, rel_H2, rel3_K, rel3_H1, rel3_H2, der3))
run.note('R.7: DISCREPANCY (arXiv TeX of DG06 against the derivation and against (3)-(4)): (nv1) K is right; (nv1) H1 is (2/r) times the derived H1 (the printed one is not dimensionless); '
         '(nv1) H2 has the opposite sign on its zeta term. Corrected: H1 = s [2r(r+mu) zeta\' + mu zeta]/(2 (r+mu)(2 lam r - 3 mu)), H2 with -(3mu^2 - 6 lam mu r + 4 lam^2 r^2)/(2 r^2 (2 lam r - 3 mu)^2) zeta.')
lap('R.7')

# ---- R.8 (im) = (9)(2 lam r - 3 mu); (9) o (3) = id
psi = Function('psi')(r)
K3p = fD * psi + (1 - 2 * M / r) * diff(psi, r)
H23p = diff((1 - 2 * M / r) * (hD * psi + r * diff(psi, r)), r) - K3p
psi9 = r * (r - 2 * M) / ((lam_ + 1) * (lam_ * r + 3 * M)) * (H23p - r * diff(K3p, r)) + r * K3p / (lam_ + 1)
imc_typed = lam_ if mut('im_coef') else lam_ + 1                  # mutation im_coef: (lam+1) -> lam in the typed (im)
im_typed = (-2 * r ** 2 * (r + mu_) * diff(K3p, r) + (2 * lam_ * r - 3 * mu_) * r * K3p + 2 * r * (r + mu_) * H23p) / imc_typed
ok9 = cancel(together((psi9 - psi).doit())) == 0
okim = cancel(together((im_typed - (2 * lam_ * r - 3 * mu_) * psi).doit())) == 0
run.check('R.8', ok9 and okim, '(9) o (3) = id on a free psi: %s; (im)((3)(psi)) = (2 lam r - 3 mu) psi: %s' % (ok9, okim))
lap('R.8')

# ---- R.9 Lemma 1 mechanics on the corrected (nv1)
nK = (mu_ + r) * r * z1s + (lam_ * r - mu_) * z0s
n1 = r * (r + mu_) * z1s + mu_ / 2 * z0s
nK_rs = cancel(nK.subs(r, rs))
n1_rs = cancel(n1.subs(r, rs))
kK = cancel(diff(nK_rs, z1s))
k1 = cancel(diff(n1_rs, z1s))
Zg = Function('Zg')(r)
N2 = (2 * r ** 2 * (r + mu_) * (2 * lam_ * r - 3 * mu_) * diff(Zg, r, 2) - r * (4 * lam_ * r ** 2 + 6 * lam_ * mu_ * r - 3 * mu_ ** 2) * diff(Zg, r)
      + (1 if mut('n2_sign') else -1) * (3 * mu_ ** 2 - 6 * lam_ * mu_ * r + 4 * lam_ ** 2 * r ** 2) * Zg)   # mutation n2_sign: the typed sign of the zeta term
w0, w1, w2, w3 = sp.symbols('w0 w1 w2 w3')


def wj(e):
    return e.doit().subs(sp.Derivative(Zg, (r, 3)), w3).subs(sp.Derivative(Zg, (r, 2)), w2).subs(sp.Derivative(Zg, r), w1).subs(Zg, w0)


N2v = cancel(wj(N2).subs(r, rs))
N2pv = cancel(wj(diff(N2, r)).subs(r, rs))
k2 = cancel(diff(N2v, w1))
okR9 = (same(nK_rs, kK * (z1s - cstar * z0s)) and kK != 0 and same(n1_rs, k1 * (z1s - cstar * z0s)) and k1 != 0
        and same(N2v, k2 * (w1 - cstar * w0)) and k2 != 0 and same(N2pv, -(L - 2) / (2 * M) * N2v) and cancel(diff(N2pv, w2)) == 0 and cancel(diff(N2pv, w3)) == 0)
run.check('R.9', okR9, 'numerators at r_s: nK = %s (z1 - c* z0), n1 = %s (...), N2 = %s (...); N2\'(r_s) = -(L-2)/(2M) N2(r_s): %s; no zeta\'\' in N2\'(r_s): %s' % (
    sp.factor(kK), sp.factor(k1), sp.factor(k2), same(N2pv, -(L - 2) / (2 * M) * N2v), cancel(diff(N2pv, w2)) == 0))
run.note('R.9: K = nK/(r^2 (2 lam r - 3 mu)), H1 = s n1/((r+mu)(2 lam r - 3 mu)), H2 = N2/(2 r^2 (2 lam r - 3 mu)^2), all single/double zeros of the numerators at r_s are the one condition '
         'zeta\' = c* zeta: K, H1, H2 smooth iff zeta smooth and zeta\'(r_s) = c* zeta(r_s) (Lemma 1 with the derived sign).')
lap('R.9')

# ---- R.10 apparent singularity of the 2-dimensional (K, H1) flow
F2 = [cancel(F3[0].subs(Hh0, Hexp)), cancel(F3[1].subs(Hh0, Hexp))]
Fm = sp.Matrix([[diff(F2[i], [K0, H10][j]) for j in range(2)] for i in range(2)]).applyfunc(cancel)
Rres = (Fm * (L * r + 6 * M - 2 * r) / (L - 2)).applyfunc(cancel).subs(r, rs).applyfunc(cancel)       # (r - r_s) F2 at r_s
F1m = (Fm - Rres / ((L * r + 6 * M - 2 * r) / (L - 2))).applyfunc(cancel).subs(r, rs).applyfunc(cancel)  # next Laurent coefficient
if mut('log_inject'):
    F1m = F1m + sp.Matrix([[0, 0], [1, 0]])                      # mutation log_inject: a perturbation of the next Laurent coefficient (after R.9)
eig_ok, nolog = laurent_test(Rres, F1m)
u0 = (Rres.nullspace()[0]).applyfunc(cancel)
wl = ((Rres - sp.eye(2)).T.nullspace()[0]).applyfunc(cancel)
run.check('R.10', eig_ok and nolog, 'residue matrix at r_s has trace 1 and det 0 (eigenvalues 0 and 1): %s; resonance compatibility w^T F1 u0 = 0 (no logarithm): %s' % (eig_ok, nolog))
run.note('R.10: residue matrix R = %s ; eigenvalues {0, 1}; w^T F1 u0 = %s' % (Rres.applyfunc(sp.factor).tolist(), sp.factor(cancel((wl.T * F1m * u0)[0, 0]))))
lap('R.10')

# ---- R.11 M = -1 and l = 2, 3, 4
r6c = (3 if mut('rs6c_typo') else 6)                                # mutation rs6c_typo: r_s at M = -1 typed 3/c
ok11 = True
for lv in (2, 3, 4):
    sub = {M: -1, L: lv * (lv + 1)}
    cc_ = (lv + 2) * (lv - 1)
    rs_n = rs.subs(sub)
    cst_n = cstar.subs(sub)
    ok11 = ok11 and cancel(rs_n - sp.Rational(r6c, cc_)) == 0 and cancel(aH.subs(sub).subs(r, rs_n)) == 0
    ok11 = ok11 and cancel(cst_n + 1 / (rs_n * (rs_n + 2))) == 0 and cancel(cst_n + cz_c.subs(sub)) == 0
    Trn = Tm.subs(sub).subs(r, rs_n).applyfunc(cancel)
    ok11 = ok11 and cancel(Trn.det()) == 0 and Trn != sp.zeros(2, 2)
    z0n = cancel(zeta.subs(sub).subs(r, rs_n))
    z2n = cancel(z2.subs(sub).subs(r, rs_n))
    aKn, a1n = cancel(aK.subs(sub).subs(r, rs_n)), cancel(a1.subs(sub).subs(r, rs_n))
    H1S = -aKn * K0 / a1n
    ratio = cancel(z2n.subs(H10, H1S) / z0n.subs(H10, H1S))
    ok11 = ok11 and cancel(ratio + 9 * s ** 2 / (lv * (lv + 1) + 1) ** 2) == 0
    run.note('R.11: l = %d (M = -1): r_s = %s = 6/c, c = %d; c* = %s = -1/(r_s (r_s+2)); typed (cz) constant = %s; zeta\'\'/zeta on S at r_s = %s' % (
        lv, rs_n, cc_, cst_n, cancel(cz_c.subs(sub)), sp.factor(ratio)))
run.check('R.11', ok11, 'M = -1, l = 2, 3, 4: r_s = 6/c, aH(r_s) = 0, T(r_s) rank 1, c* = -1/(r_s(r_s+2)) = -(typed constant), zeta\'\'/zeta = -9 s^2/(L+1)^2 on S')
lap('R.11')

# ---- R.12 a concrete smooth zeta (M = -1, l = 2): the corrected (nv1) is regular at r_s iff zeta' = c* zeta there; the typed (cz) constant gives poles
sub12 = {M: -1, L: 6}
rs12 = rs.subs(sub12)
lam12, mu12 = lam_.subs(sub12), mu_.subs(sub12)
cst12 = cstar.subs(sub12)
czt12 = cz_c.subs(sub12)


def nv1c(Zexpr):
    """corrected (nv1) K, H1, H2 for an explicit zeta(r) (M = -1, l = 2, s symbolic)"""
    Zp, Zpp = diff(Zexpr, r), diff(Zexpr, r, 2)
    K_ = ((mu12 + r) * r * Zp + (lam12 * r - mu12) * Zexpr) / (r ** 2 * (2 * lam12 * r - 3 * mu12))
    H1_ = s * (r * (r + mu12) * Zp + mu12 / 2 * Zexpr) / ((r + mu12) * (2 * lam12 * r - 3 * mu12))
    N2_ = (2 * r ** 2 * (r + mu12) * (2 * lam12 * r - 3 * mu12) * Zpp - r * (4 * lam12 * r ** 2 + 6 * lam12 * mu12 * r - 3 * mu12 ** 2) * Zp
           - (3 * mu12 ** 2 - 6 * lam12 * mu12 * r + 4 * lam12 ** 2 * r ** 2) * Zexpr)
    H2_ = N2_ / (2 * r ** 2 * (2 * lam12 * r - 3 * mu12) ** 2)
    return K_, H1_, H2_


def regular12(expr):
    num, den = sp.fraction(cancel(together(expr)))
    return cancel(den.subs(r, rs12)) != 0


c_use = czt12 if mut('ex_typed') else cst12                       # mutation ex_typed: the conforming example is built with the typed (cz) constant
Z_ok = 1 + c_use * (r - rs12) + (r - rs12) ** 2 + (r - rs12) ** 3
Z_bad = sp.Integer(1) + 0 * r                                       # zeta' = 0 != c* zeta
Z_typ = 1 + czt12 * (r - rs12)                                      # zeta' = typed constant * zeta
reg_ok = [regular12(v_) for v_ in nv1c(Z_ok)]
reg_bad = [regular12(v_) for v_ in nv1c(Z_bad)]
reg_typ = [regular12(v_) for v_ in nv1c(Z_typ)]
run.check('R.12', all(reg_ok) and not any(reg_bad) and not any(reg_typ),
          'M = -1, l = 2, r_s = %s: zeta with zeta\' = c* zeta: (K, H1, H2) regular %s; zeta = 1: %s (poles wanted); zeta with the typed (cz) constant: %s (poles wanted)' % (
              rs12, reg_ok, reg_bad, reg_typ))
run.note('R.12: l = 2, M = -1: c* = %s, typed (cz) constant = %s; a smooth zeta with zeta\'(r_s) = c* zeta(r_s) gives smooth (K, H1, H2); one with the typed constant gives poles in K, H1 and H2 at r_s' % (cst12, czt12))
lap('R.12')

# ---- R.S.1, R.S.2 the static case
st0 = {s: 0, H10: 0, H11: 0}
tr_s = cancel(eqsH['tr'].subs(s, 0))
okS1a = all(cancel(diff(tr_s, v_)) == 0 for v_ in (K0, K1, Hh0, Hh1, H11)) and cancel(diff(tr_s, H10) - L / (2 * r ** 2)) == 0
rx_s = cancel(eqsH['rx'].subs(st0))
Grr_s = cancel(Grr.subs(st0))
unk2 = [K1, Hh1]
Mat2 = sp.Matrix([[diff(rx_s, v) for v in unk2], [diff(Grr_s, v) for v in unk2]]).applyfunc(cancel)
det2 = cancel(Mat2.det())
rhs2 = sp.Matrix([-(e_ - sum(diff(e_, v) * v for v in unk2)) for e_ in (rx_s, Grr_s)]).applyfunc(cancel)
sol2 = Mat2.LUsolve(rhs2).applyfunc(cancel)
FKs, FHs = sol2[0], sol2[1]
c_K = 2 * lam_ / M if mut('sflow_typo') else lam_ / M                 # mutation sflow_typo: typed static K' coefficient lam/M -> 2 lam/M
FKs_t = c_K * K0 - (lam_ * r + 2 * M) / (M * r) * Hh0
FHs_t = lam_ / M * K0 - (lam_ * r ** 2 - 2 * (lam_ - 1) * M * r - 2 * M ** 2) / (M * r * (r - 2 * M)) * Hh0
okS1 = okS1a and same(FKs, FKs_t) and same(FHs, FHs_t) and same(det2, -M / (2 * r * (r - 2 * M)))
run.check('R.S.1', okS1, 'static tr equation is (L/(2r^2)) H1: %s; derived static flow K\' %s, H\' %s, det = -M/(2r(r-2M)) %s (#145)' % (
    okS1a, same(FKs, FKs_t), same(FHs, FHs_t), same(det2, -M / (2 * r * (r - 2 * M)))))
run.note('R.S.1: static flow K\' = %s ; H\' = %s' % (sp.factor(FKs), sp.factor(FHs)))
lap('R.S.1')
Fs = [FKs, FHs]
if mut('sflow_h'):
    Fs = [FKs - Hh0 * diff(FKs, Hh0) / 2, FHs]                       # mutation sflow_h: the coefficient of H in the static K' halved (after R.S.1)
sts = [K0, Hh0]


def D2s(e):
    return diff(e, r) + sum(diff(e, sts[i]) * Fs[i] for i in range(2))


met_s = set()
dens_s = True
for e_ in (FKs, FHs):
    for v_ in sts:
        co = cancel(diff(e_, v_))
        if co != 0:
            for fac_, mult_ in denominators_factors(co):
                met_s.add(classify_r(fac_, [M, r, r - 2 * M]))
                if cancel(fac_.subs(r, rs)) == 0:
                    dens_s = False
zeta_s = cancel((-2 * r ** 2 * (r + mu_) * Fs[0] + (2 * lam_ * r - 3 * mu_) * r * K0 + 2 * r * (r + mu_) * Hh0) / imc_true)
zs1 = cancel(D2s(zeta_s))
Ts = sp.Matrix([[diff(zeta_s, K0), diff(zeta_s, Hh0)], [diff(zs1, K0), diff(zs1, Hh0)]]).applyfunc(cancel)
detTs = cancel(Ts.det())
Tsr = Ts.subs(r, rs).applyfunc(cancel)
rank1s = cancel(Tsr.det()) == 0 and Tsr != sp.zeros(2, 2)
rowrel_s = [cancel(Tsr[1, j] - cstar * Tsr[0, j]) for j in range(2)]
okS2 = (not any(m_.startswith('UNEXPECTED') for m_ in met_s) and dens_s and same(detTs, -8 * r ** 2 * (lam_ * r + 3 * M) ** 2 / (L * M))
        and rank1s and rowrel_s == [0, 0])
run.check('R.S.2', okS2, 'static: denominator factors %s (none vanishes at r_s: %s); det T_s = -8 r^2 (lam r + 3M)^2/(L M): %s; T_s(r_s) rank 1: %s; row2 = c* row1 (c* of R.4): %s' % (
    sorted(met_s), dens_s, same(detTs, -8 * r ** 2 * (lam_ * r + 3 * M) ** 2 / (L * M)), rank1s, rowrel_s))
run.note('R.S.2: every static solution has zeta\'(r_s) = c* zeta(r_s), with the same c* = %s as for e^{st}; no relation A is needed (there is none at s = 0)' % sp.factor(cstar))
lap('R.S.2')
run.finish()
