"""#143: library of zerilli_derivation.py. Harness, dual numbers, Christoffel and Ricci code, an independent covariant route,
the Legendre reduction and the angular split. Nothing here knows the answer (V).

Coordinates (t, r, x, phi), x = cos(theta): the unit sphere is dx^2/(1-x^2) + (1-x^2) dphi^2, everything is rational in (r, x).
First order in eps is carried by dual numbers (a0 + eps a1, eps^2 = 0).

Harness: PASS <id> / FAIL <id>. Exit 0: all pass. Exit 1: a real check failed. Exit 2: a planted check failed (stop before real data).
--mutate=NAME[,NAME] switches on a deliberate error (grep "mut(" for the places). --failfast ends at the first FAIL.
"""
import sys
import sympy as sp
from sympy import Symbol, Rational as Q, cancel, together, factor, diff, Function

_MUT = set()
_FAILFAST = '--failfast' in sys.argv
for _a in sys.argv[1:]:
    if _a.startswith('--mutate='):
        _MUT |= set(_a[len('--mutate='):].split(','))


def mut(name):
    return name in _MUT


class Run:
    """planted(): a failure exits 2 at once. check(): failures counted, exit 1 at finish()."""

    def __init__(self, title):
        self.nfail = 0
        self.npass = 0
        print('#', title, flush=True)
        if _MUT:
            print('# mutations on:', ','.join(sorted(_MUT)), flush=True)

    def planted(self, cid, cond, why=''):
        if cond:
            print('PASS', cid, flush=True)
            self.npass += 1
        else:
            print('FAIL', cid + ':', why, flush=True)
            print('# planted check failed: stop before any real data', flush=True)
            sys.exit(2)

    def check(self, cid, cond, why=''):
        if cond:
            print('PASS', cid, flush=True)
            self.npass += 1
        else:
            print('FAIL', cid + ':', why, flush=True)
            self.nfail += 1
            if _FAILFAST:
                self.finish()

    def note(self, text):
        print('#', text, flush=True)

    def finish(self):
        print('# summary: %d PASS, %d FAIL' % (self.npass, self.nfail), flush=True)
        sys.exit(1 if self.nfail else 0)


# ------------------------------------------------------------------ symbols
t = Symbol('t', real=True)
r = Symbol('r', positive=True)
x = Symbol('x', real=True)
ph = Symbol('phi', real=True)
M = Symbol('M', real=True)
s = Symbol('s')
l = Symbol('l', positive=True)
L = l * (l + 1)
lam = (l - 1) * (l + 2) / 2
COORD = (t, r, x, ph)
N = 4
P = Function('P')(x)          # P_l(x), Legendre's equation (1-x^2)P'' - 2xP' + L P = 0


def fM(Mv=None):
    """f = 1 - 2M/r. Mutation fsign_bg: the wrong sign, used by the background metric."""
    Mv = M if Mv is None else Mv
    return 1 - 2 * Mv / r


def f_bg():
    return 1 + 2 * M / r if mut('fsign_bg') else 1 - 2 * M / r


# ------------------------------------------------------------------ dual numbers a0 + eps a1
class Dual:
    __slots__ = ('a0', 'a1')

    def __init__(self, a0, a1=0):
        self.a0 = sp.sympify(a0)
        self.a1 = sp.sympify(a1)

    def __add__(self, o):
        o = o if isinstance(o, Dual) else Dual(o)
        return Dual(self.a0 + o.a0, self.a1 + o.a1)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.a0, -self.a1)

    def __sub__(self, o):
        o = o if isinstance(o, Dual) else Dual(o)
        return Dual(self.a0 - o.a0, self.a1 - o.a1)

    def __mul__(self, o):
        o = o if isinstance(o, Dual) else Dual(o)
        return Dual(self.a0 * o.a0, self.a0 * o.a1 + self.a1 * o.a0)

    __rmul__ = __mul__

    def d(self, v):
        return Dual(diff(self.a0, v), diff(self.a1, v))


ZERO = Dual(0)


def dsum(it):
    out = ZERO
    for z in it:
        out = out + z
    return out


def build_metric(g0, h):
    """g = g0 + eps h, ginv = g0^-1 - eps g0^-1 h g0^-1 (exact to first order). Returns (g, ginv) as 4x4 lists of Dual."""
    G0 = sp.Matrix(g0)
    H = sp.Matrix(h)
    G0i = G0.inv().applyfunc(cancel)
    sign = 1 if mut('dual_inv') else -1                       # mutation dual_inv: +eps g^-1 h g^-1
    G1i = (sign * G0i * H * G0i).applyfunc(cancel)
    g = [[Dual(G0[i, j], H[i, j]) for j in range(N)] for i in range(N)]
    gi = [[Dual(G0i[i, j], G1i[i, j]) for j in range(N)] for i in range(N)]
    return g, gi


def christoffel(g, gi):
    """Gamma^a_bc = 1/2 g^ad (d_b g_dc + d_c g_db - d_d g_bc), first order exact."""
    Gam = [[[None] * N for _ in range(N)] for _ in range(N)]
    for a in range(N):
        for b in range(N):
            for c in range(b, N):
                acc = ZERO
                for dd in range(N):
                    acc = acc + gi[a][dd] * (g[dd][c].d(COORD[b]) + g[dd][b].d(COORD[c]) - g[b][c].d(COORD[dd]))
                acc = Dual(Q(1, 2) * acc.a0, Q(1, 2) * acc.a1)
                Gam[a][b][c] = acc
                Gam[a][c][b] = acc
    return Gam


def ricci(Gam):
    """R_mn = d_a Gam^a_mn - d_n Gam^a_ma + Gam^a_ab Gam^b_mn - Gam^a_nb Gam^b_ma. Returns the four groups per component.
    Mutation drop_trace: the first-order part of the second group is dropped (the background part is kept)."""
    R = [[None] * N for _ in range(N)]
    for m in range(N):
        for n in range(m, N):
            T1 = dsum(Gam[a][m][n].d(COORD[a]) for a in range(N))
            T2 = -dsum(Gam[a][m][a].d(COORD[n]) for a in range(N))
            if mut('drop_trace'):
                T2 = Dual(T2.a0, 0)
            T3 = dsum(Gam[a][a][b] * Gam[b][m][n] for a in range(N) for b in range(N))
            T4 = -dsum(Gam[a][n][b] * Gam[b][m][a] for a in range(N) for b in range(N))
            tot = T1 + T2 + T3 + T4
            R[m][n] = tot
            R[n][m] = tot
    return R


# ------------------------------------------------------------------ Legendre reduction
def leg_red(e):
    """Remove derivatives of P of order >= 2 with (1-x^2)P'' = 2xP' - L P. Result is linear in P, P' with rational coefficients."""
    sgn = -1 if mut('leg_sign') else 1                        # mutation leg_sign: P'' = (2xP' + L P)/(1-x^2)
    rules = {}
    P0, P1 = sp.Symbol('P0'), sp.Symbol('P1')                 # P, P'
    Pxx = (2 * x * P1 + sgn * (-L) * P0) / (1 - x ** 2)       # P'' as a function of (P, P', x)

    def d_x(expr):
        """total x-derivative of an expression in (P0, P1, x): P0' = P1, P1' = P''."""
        return diff(expr, x) + diff(expr, P0) * P1 + diff(expr, P1) * Pxx

    e = e.doit()
    # replace derivatives of P by symbols, highest first
    ders = sorted([d for d in e.atoms(sp.Derivative) if d.expr == P and d.variables and set(d.variables) == {x}],
                  key=lambda d: -len(d.variables))
    jets = {0: P0, 1: P1}
    cur = P1
    for n in range(2, max([len(d.variables) for d in ders] + [1]) + 1):
        cur = d_x(cur)
        jets[n] = cancel(cur)
    for dv in ders:
        e = e.subs(dv, jets[len(dv.variables)])
    e = e.subs(P, P0)
    return e


def ang_split(e):
    """e linear in P0, P1: returns (A, B) with e = A P0 + B P1, A, B rational in x (and the rest)."""
    P0, P1 = sp.Symbol('P0'), sp.Symbol('P1')
    e = sp.expand(e)
    A = cancel(together(e.coeff(P0, 1)))
    B = cancel(together(e.coeff(P1, 1)))
    rest = cancel(together(e - A * P0 - B * P1))
    return A, B, rest


# ------------------------------------------------------------------ the covariant route (independent of the dual-number code)
def christoffel0(g0, g0i):
    G = [[[0] * N for _ in range(N)] for _ in range(N)]
    for a in range(N):
        for b in range(N):
            for c in range(N):
                G[a][b][c] = cancel(sum(Q(1, 2) * g0i[a, dd] * (diff(g0[dd, c], COORD[b]) + diff(g0[dd, b], COORD[c]) - diff(g0[b, c], COORD[dd]))
                                        for dd in range(N)))
    return G


def route2_delta_ricci(g0, h):
    """dR_mn = nabla_a dGamma^a_mn - nabla_n dGamma^a_ma, dGamma^a_bc = 1/2 g^ad (nabla_b h_dc + nabla_c h_db - nabla_d h_bc).
    Uses the background Gamma and h only. Mutation route2_drop drops the second term."""
    g0 = sp.Matrix(g0)
    g0i = g0.inv().applyfunc(cancel)
    G = christoffel0(g0, g0i)
    # nabla_a h_bc
    Dh = [[[diff(h[b][c], COORD[a]) - sum(G[e][a][b] * h[e][c] for e in range(N)) - sum(G[e][a][c] * h[b][e] for e in range(N))
            for c in range(N)] for b in range(N)] for a in range(N)]            # Dh[a][b][c] = nabla_a h_bc
    dG = [[[sum(Q(1, 2) * g0i[a, dd] * (Dh[b][dd][c] + Dh[c][dd][b] - Dh[dd][b][c]) for dd in range(N))
            for c in range(N)] for b in range(N)] for a in range(N)]            # dG[a][b][c] = dGamma^a_bc
    # nabla_e dGamma^a_bc
    DdG = [[[[diff(dG[a][b][c], COORD[e]) + sum(G[a][e][f] * dG[f][b][c] for f in range(N))
              - sum(G[f][e][b] * dG[a][f][c] for f in range(N)) - sum(G[f][e][c] * dG[a][b][f] for f in range(N))
              for c in range(N)] for b in range(N)] for a in range(N)] for e in range(N)]  # DdG[e][a][b][c]
    R1 = [[None] * N for _ in range(N)]
    for m in range(N):
        for n in range(N):
            first = sum(DdG[a][a][m][n] for a in range(N))
            second = 0 if mut('route2_drop') else sum(DdG[n][a][m][a] for a in range(N))
            R1[m][n] = first - second
    return R1


# ------------------------------------------------------------------ driver helpers
def lin_ricci(g0, h):
    """Return (R0, R1): background Ricci and the first-order Ricci, both 4x4 lists of sympy expressions (unreduced)."""
    g, gi = build_metric(g0, h)
    Gam = christoffel(g, gi)
    R = ricci(Gam)
    return [[R[m][n].a0 for n in range(N)] for m in range(N)], [[R[m][n].a1 for n in range(N)] for m in range(N)]


def tidy(e):
    """cancel a rational expression in r, x, M, s, l (and the symbols P0, P1) component-wise."""
    return cancel(together(e))


# ------------------------------------------------------------------ jets, angular coefficient equations, 2D master reduction
def jetify(e, names):
    """Replace F(r) and its r-derivatives by symbols F_0, F_1, ... for the functions named in `names`."""
    e = e.doit()
    for nm in names:
        Fn = Function(nm)(r)
        ders = sorted([d for d in e.atoms(sp.Derivative) if d.expr == Fn], key=lambda d: -d.derivative_count)
        for dv in ders:
            e = e.subs(dv, Symbol('%s_%d' % (nm, dv.derivative_count)))
        e = e.subs(Fn, Symbol('%s_0' % nm))
    return e


def coef_eqs(e, norm, basis, cid):
    """e: reduced component (linear in P0, P1). Multiply by norm; the angular form must be A P0 + B x P1 with x-independent A and B,
    and the structures not in `basis` must vanish. Returns ({'P0': A, 'xP1': B/...}, ok, why)."""
    P0, P1 = sp.Symbol('P0'), sp.Symbol('P1')
    e = cancel(together(e * norm))
    A, B, rest = ang_split(e)
    B = cancel(B / x)
    ok = rest == 0
    why = []
    if rest != 0:
        why.append('not linear in P, P\'')
    out = {}
    B_raw = cancel(B * x)                                     # coefficient of P' itself (B above is the coefficient of x P')
    slots = {'P0': (A, 'P0' in basis), 'P1': (B_raw if 'P1' in basis else B, ('P1' in basis) or ('xP1' in basis))}
    for slot, (val, allowed) in slots.items():
        if allowed:
            if cancel(diff(val, x)) != 0:
                ok = False
                why.append('%s coefficient depends on x' % slot)
            out[slot if slot == 'P0' else ('P1' if 'P1' in basis else 'xP1')] = cancel(val)
        elif val != 0:
            ok = False
            why.append('structure %s present but not allowed' % slot)
    return out, ok, '; '.join(why)


def master_reduce(u, flow, psi, ftort):
    """Two-dimensional first-order system d u_i/dr = flow_i(u, r) (u = [u1, u2], linear forms). psi = linear form in u.
    d/dx = ftort d/dr. Returns (a, b, det T) with psi_xx = a psi + b psi_x exactly, T the matrix (psi, psi_x) = T u."""
    def D(expr):
        return diff(expr, r) + sum(diff(expr, u[i]) * flow[i] for i in range(2))

    def Dx(expr):
        return ftort * D(expr)

    psi_x = cancel(Dx(psi))
    psi_xx = cancel(Dx(psi_x))
    T = sp.Matrix([[diff(psi, u[0]), diff(psi, u[1])], [diff(psi_x, u[0]), diff(psi_x, u[1])]]).applyfunc(cancel)
    w = sp.Matrix([[diff(psi_xx, u[0]), diff(psi_xx, u[1])]]).applyfunc(cancel)
    det = cancel(T.det())
    ab = (w * T.inv()).applyfunc(cancel)
    return cancel(ab[0]), cancel(ab[1]), det


def jet_sym(nm, n):
    return Symbol('%s_%d' % (nm, n))


def flow_sub(expr, varmap, states, F, nmax=4):
    """Replace the jets nm_n of the functions in varmap ({nm: expression in the states}) by n-fold total derivatives along the
    flow d states_i / dr = F_i(states, r)."""
    def D(e):
        return diff(e, r) + sum(diff(e, states[i]) * F[i] for i in range(len(states)))
    rep = {}
    for nm, S in varmap.items():
        cur = S
        for n in range(nmax + 1):
            rep[jet_sym(nm, n)] = cancel(cur)
            cur = D(cur)
    return expr.subs(rep, simultaneous=True)


def denominators_factors(expr):
    """irreducible factors of the denominator of a rational function (list of (factor, multiplicity))."""
    num, den = sp.fraction(cancel(together(expr)))
    return sp.factor_list(den)[1]
