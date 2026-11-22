"""
Exact check that the reduced polar system  y' = J y  (polar_matrix_generated.py, from reduce_polar.py) satisfies
the separated (tt) equation and the angular trace of derive_nonradial.py (polar_eqs.pkl, H2 = H0), on solutions
of y' = J y.  verify_polar.py does this to 7e-7 with finite-difference J'; here the arithmetic is exact.
Planted controls through the same code path: (t th), (t r), (r th) and also (rr), KG_R, KG_I.

Setup (as verify_polar.py).  y = (H1/s, K, F, F', G/s, G'/s),  H1 = s y0, K = y1, F = y2, G = s y4,  H0 = H2 = h.y,
    y' = J y,   y'' = J' y + J y' = (J' + J J) y,   H0^(k) by the product rule (h', h'' = d/dr of h).
y is six free symbols: the residual of an equation is sum_i c_i y_i, and it is zero for every y iff all six c_i
are zero.  The script prints the six c_i (zero or not) and, at points, also the residual at a random rational y.

Background.  J and h are functions of (r, alpha, a, phi0, u, omega, epsilon, s, LL, Lambda, mu), with
alpha = e^{nu/2}, a = e^{lam/2}, u = phi0'/phi0.  Their r-derivatives go by the chain rule through the background,
with the first derivatives written as lincommon.bgsub writes them (lincommon.background_rules: dnu0, dlam0, ddp0):
    alpha' = alpha nu0'/2  with nu0'  = dnu0,       a' = a lam0'/2  with lam0' = dlam0,
    phi0'  = u phi0        (definition of u),        u' = phi0''/phi0 - u^2  with phi0'' = ddp0.
The equations E[...] carry the background only through nu0, lam0, phi0 (bgsub was applied by derive_nonradial.py);
exp(c nu0 + ...) becomes alpha^(2c), exp(c lam0) becomes a^(2c).

Method.  The straight-line program of J_and_h (source of polar_matrix_generated.py, read, with its two numpy
imports dropped; the printer wrote five exponents as -1.0 and one 1/3 as a Python float quotient, which are turned
into -1 and Rational(1, 3) by an AST rewrite so that no float enters) is run on truncated Taylor series in r of
order 2 whose coefficients are exact elements of a field or ring of rational functions over Q.  The order-1 coefficient of
a quantity is its total r-derivative by the chain rule through alpha, a, phi0, u (so J', h'), the order-2
coefficient is half its second derivative (h'').  pi enters as a generator PI of the number field (the identity
holds for any value of the coupling constant, so it holds for pi); epsilon is a generator too, not +-1.
Two exact number types carry the series:
    --engine=frac  sympy's rational function field: every operation cancels by gcd (canonical form).  At twelve
                   variables it did not get past J' in 5 min 48 s (a 300 s deadline, delayed by one gcd; raw I8-rpi5).
    --engine=af    gcd-free: an element is N / prod(atom_k^e_k) with N a Laurent polynomial over Q in the free
                   variables, the atoms the polynomials the program inverts (here one, an 8-term polynomial, the
                   same up to monomials and scale everywhere); sums bring both terms to the common denominator by
                   multiplying numerators by atom powers.  It is zero iff N is zero as a polynomial.
--method=points: every variable but PI is a random rational number (identity test at points, not a symbolic proof);
                 engine frac by default.
--method=symbolic: --free=all (default) keeps all twelve variables r, alpha, a, phi0, u, omega, epsilon, s, LL, Lambda,
                 mu, PI as generators; engine af by default.  A residual that is exactly the zero polynomial is an
                 exact computer-assisted identity in twelve variables (checked symbolically).  --free=r,s,... keeps
                 only those (and PI) free, the others at the base points of --seed (a ladder for the engines).
                 --deadline=<seconds> stops it with exit code 2 (not finished); a long gcd delays the stop.
Cross-checks: --selftest (series derivatives against sympy's chain rule on the symbolic J, h, at points),
--crosscheck (the symbolic objects of the af engine, all twelve variables free, evaluated at random rational points
against the field engine's exact values there: J, J', h, h', h'', every row of y, y', y'' and every coefficient of
the equations).  What the engines share: the program of J_and_h and its AST rewrite, the series class, the first
background derivatives and the extraction of the equations' coefficients from polar_eqs.pkl.

Pass criteria (written before the runs).
  - a residual vector counts as zero only if all six coefficients are the zero element;
  - points: for each of the eight equations all N points give zero (N = 20, epsilon = +1 and -1 alternating);
  - symbolic: all eight equations give zero.
  - Exit code 1 if any check fails.  Exit code 0 if all pass.
--selftest compares the series derivatives J', h', h'' with sympy's chain rule on the symbolic J, h.
Controls (--mutate), each must give a nonzero residual in (tt) or the angular trace at every point:
  --mutate=coef[:i,j]  J[i, j] is scaled by 1 + 1/7 (J' with it, being the derivative of the scaled function);
                       default J[0, 0]; --sensitivity lists every entry of J.
  --mutate=bg[:name]   one background derivative gets the wrong sign: alpha (default), a, phi0 or u.
  --mutate=eq          the coefficient of H0'' in (tt) and of K'' in the angular trace is scaled by 1 + 1/7.
A mutated run exits 1 when the mutation is detected (the check fails, as it must) and 0 when it is not detected.

Usage:  python verify_polar_exact.py [--method=points|symbolic] [--npoints=20] [--seed=8] [--mutate=...]
                                     [--sensitivity] [--selftest] [--crosscheck] [--deadline=SECONDS] [--eqs=tt,ang_tr,...]
                                     [--engine=frac|af] [--free=all|r,s,...]
"""
import os, sys, ast, time, pickle, random, signal, hashlib, argparse, subprocess
from fractions import Fraction
import sympy as sp
from sympy import QQ
from sympy.polys.fields import field
from sympy.polys.rings import ring

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.setrecursionlimit(200000)
import lincommon as lc

T0 = time.time()
NAMES = ['r', 'alpha', 'a', 'phi0', 'u', 'omega', 'epsilon', 's', 'LL', 'Lambda', 'mu', 'PI']
Q = {n: sp.Symbol(n) for n in NAMES}
BGNAMES = ('alpha', 'a', 'phi0', 'u')
EQS_PLANTED = ('tth', 'tr', 'rth', 'rr', 'KGR', 'KGI')
EQS_TARGET = ('tt', 'ang_tr')
SCALE = sp.Rational(8, 7)


def say(*a):
    print(*a, flush=True)


def tm():
    return "%.0fs" % (time.time() - T0)


# ------------------------------------------------------------------ gcd-free exact numbers
class AFField:
    """Elements N / prod(atom_k^e_k): N a Laurent polynomial over Q in the free variables, the atoms the inverted
    polynomials (monic, monomial factor removed). No gcd is ever taken; an element is zero iff N is zero."""
    def __init__(self, names):
        self.R, *self.gens = ring(','.join(names), QQ)
        self.names = names
        self.atoms, self.aid, self.pwc = [], {}, {}
        self.zero, self.one = AF(self, self.R.zero, {}), AF(self, self.R.one, {})

    def pw(self, den):
        key = tuple(sorted(den.items()))
        v = self.pwc.get(key)
        if v is None:
            v = self.R.one
            for k, e in key:
                v = v * self.atoms[k] ** e
            self.pwc[key] = v
        return v

    def poly(self, d):
        return self.R.from_dict(d)

    def from_expr(self, expr):
        out = {}
        for term in sp.Add.make_args(sp.expand(expr)):
            c, rest = term.as_coeff_Mul()
            assert c.is_Rational, term
            ev = [0] * len(self.names)
            if rest != 1:
                for b, e in rest.as_powers_dict().items():
                    assert b.is_Symbol and e.is_Integer and b.name in self.names, (term, b, e)
                    ev[self.names.index(b.name)] += int(e)
            ev = tuple(ev)
            out[ev] = out.get(ev, 0) + QQ.from_sympy(c)
        return AF(self, self.poly({k: v for k, v in out.items() if v}), {})


class AF:
    __slots__ = ('F', 'N', 'den')

    def __init__(self, F, N, den):
        self.F, self.N, self.den = F, N, (den if N else {})

    def _co(self, o):
        return o if isinstance(o, AF) else AF(self.F, self.F.R.from_dict({(0,) * len(self.F.names): QQ(int(o))}), {})

    def __bool__(self):
        return bool(self.N)

    def __neg__(self):
        return AF(self.F, -self.N, self.den)

    def __add__(self, o):
        o = self._co(o)
        if self.den == o.den:
            return AF(self.F, self.N + o.N, self.den)
        den = {k: max(self.den.get(k, 0), o.den.get(k, 0)) for k in set(self.den) | set(o.den)}
        d1 = {k: v - self.den.get(k, 0) for k, v in den.items() if v > self.den.get(k, 0)}
        d2 = {k: v - o.den.get(k, 0) for k, v in den.items() if v > o.den.get(k, 0)}
        N1 = self.N * self.F.pw(d1) if d1 else self.N
        N2 = o.N * self.F.pw(d2) if d2 else o.N
        return AF(self.F, N1 + N2, den)
    __radd__ = __add__

    def __sub__(self, o):
        return self + (-self._co(o))

    def __rsub__(self, o):
        return self._co(o) + (-self)

    def __mul__(self, o):
        o = self._co(o)
        den = dict(self.den)
        for k, v in o.den.items():
            den[k] = den.get(k, 0) + v
        return AF(self.F, self.N * o.N, den)
    __rmul__ = __mul__

    def inv(self):
        F, N = self.F, self.N
        if not N:
            raise ZeroDivisionError
        ks = list(N)
        mon = tuple(min(k[i] for k in ks) for i in range(len(ks[0])))
        P = F.R.from_dict({tuple(e - m for e, m in zip(k, mon)): c for k, c in N.items()})
        lc = P.LC
        inv_mon = F.poly({tuple(-e for e in mon): QQ(1) / lc})
        num = inv_mon * F.pw(self.den)
        if len(P) == 1:
            return AF(F, num, {})
        A = P.quo_ground(lc)
        k = F.aid.get(A)
        if k is None:
            k = F.aid[A] = len(F.atoms)
            F.atoms.append(A)
        return AF(F, num, {k: 1})


def finv(x):
    return x.inv() if isinstance(x, AF) else KF.one / x


# ------------------------------------------------------------------ truncated Taylor series
KF = None          # the field of the current context
_kc = {}


def kint(n):
    if n not in _kc:
        _kc[n] = Jet([KF.from_expr(sp.Integer(n)), KF.zero, KF.zero])
    return _kc[n]


def lift(x):
    return x if isinstance(x, Jet) else kint(int(x)) if isinstance(x, int) else NotImplemented


class Jet:
    """Taylor coefficients [f, f', f''/2] in r, each an element of KF; coefficients 1 and 2 are computed on demand."""
    __slots__ = ('c', 'g')

    def __init__(self, c, g=None):
        self.c, self.g = c, g

    def __getitem__(self, k):
        v = self.c[k]
        if v is None:
            v = self.c[k] = self.g(k)
        return v

    def __add__(self, o):
        a, b = self, lift(o)
        return Jet([a[0] + b[0], None, None], lambda k: a[k] + b[k])
    __radd__ = __add__

    def __neg__(self):
        a = self
        return Jet([-a[0], None, None], lambda k: -a[k])

    def __sub__(self, o):
        a, b = self, lift(o)
        return Jet([a[0] - b[0], None, None], lambda k: a[k] - b[k])

    def __rsub__(self, o):
        return lift(o).__sub__(self)

    def __mul__(self, o):
        a, b = self, lift(o)

        def g(k):
            t = a[0] * b[k]
            for i in range(1, k + 1):
                t = t + a[i] * b[k - i]
            return t
        return Jet([a[0] * b[0], None, None], g)
    __rmul__ = __mul__

    def inv(self):
        b = self
        i0 = finv(b[0])

        def g(k):
            if k == 1:
                return -(b[1] * i0 * i0)
            return (b[1] * b[1] - b[0] * b[2]) * i0 * i0 * i0
        return Jet([i0, None, None], g)

    def __truediv__(self, o):
        return self * lift(o).inv()

    def __rtruediv__(self, o):
        return lift(o) * self.inv()

    def __pow__(self, n):
        if not isinstance(n, int):
            raise TypeError("exponent %r" % (n,))
        if n < 0:
            return self.inv() ** (-n)
        out, base = None, self
        while n:
            if n & 1:
                out = base if out is None else out * base
            n >>= 1
            if n:
                base = base * base
        return out if out is not None else kint(1)


class Grid:
    def __init__(self):
        self.d = {}

    def __setitem__(self, ij, v):
        self.d[ij] = v

    def __getitem__(self, ij):
        return self.d.get(ij, 0)


class NPShim:
    empty = staticmethod(lambda shape, dtype=None: Grid())
    result_type = staticmethod(lambda *a: None)
    array = staticmethod(lambda x: list(x))


# ------------------------------------------------------------------ the generated program
class Fix(ast.NodeTransformer):
    """Integral floats to ints, int / int to Rational(int, int), so that nothing inexact enters."""
    def __init__(self):
        self.nfloat = self.nrat = 0

    def visit_Constant(self, node):
        if isinstance(node.value, float):
            assert node.value == int(node.value), node.value
            self.nfloat += 1
            return ast.copy_location(ast.Constant(int(node.value)), node)
        return node

    def visit_BinOp(self, node):
        self.generic_visit(node)
        if isinstance(node.op, ast.Div) and all(isinstance(x, ast.Constant) and type(x.value) is int
                                                 for x in (node.left, node.right)):
            self.nrat += 1
            return ast.copy_location(ast.Call(ast.Name('Rational', ast.Load()), [node.left, node.right], []), node)
        return node


def load_program():
    path = os.path.join(HERE, 'polar_matrix_generated.py')
    src = open(path).read()
    tree = ast.parse(src, path)
    tree.body = [n for n in tree.body if not isinstance(n, (ast.Import, ast.ImportFrom))]
    fix = Fix()
    tree = fix.visit(tree)
    ast.fix_missing_locations(tree)
    used = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    info = dict(floats=fix.nfloat, ratios=fix.nrat, exp_sqrt_log_used=sorted(used & {'exp', 'sqrt', 'log'}),
                sha256=hashlib.sha256(src.encode()).hexdigest()[:16])
    return compile(tree, path, 'exec'), info


# ------------------------------------------------------------------ the equations and the background flow
def _exp_rule(e):
    arg = e.args[0]
    c1, c2 = arg.coeff(lc.nu0), arg.coeff(lc.lam0)
    assert sp.expand(arg - c1 * lc.nu0 - c2 * lc.lam0) == 0 and (2 * c1).is_Integer and (2 * c2).is_Integer, arg
    return Q['alpha'] ** (2 * c1) * Q['a'] ** (2 * c2)


def canon(z):
    """Background functions of r to the symbols Q: phi0' = u phi0, exp(c nu0) = alpha^(2c), exp(c lam0) = a^(2c)."""
    z = z.subs(sp.diff(lc.p0, lc.r), Q['u'] * Q['phi0']).subs(lc.p0, Q['phi0'])
    z = z.replace(lambda e: isinstance(e, sp.exp), _exp_rule)
    z = z.xreplace({lc.r: Q['r'], lc.LAM: Q['Lambda'], lc.MU: Q['mu'], sp.pi: Q['PI']})
    left = [a for a in z.atoms(sp.Function, sp.Derivative)] + [a for a in z.free_symbols if a not in Q.values()]
    return z, left


FN = {n: sp.Function(n)(lc.r) for n in ('H0', 'H1', 'H2', 'K', 'F', 'G')}
JS, REP = {}, {}
for _n, _f in FN.items():
    for _k in (3, 2, 1):
        JS[(_n, _k)] = sp.Symbol(f'{_n}_d{_k}'); REP[sp.diff(_f, lc.r, _k)] = JS[(_n, _k)]
    JS[(_n, 0)] = sp.Symbol(f'{_n}_d0')
JETSYMS = set(JS.values())


def load_equations(keys):
    d = pickle.load(open(os.path.join(HERE, 'polar_eqs.pkl'), 'rb'))
    E = d['E']
    assert d['H2sol'] == FN['H0'], d['H2sol']
    eqs = {}
    for k in keys:
        z = E[k].subs(FN['H2'], FN['H0']).subs(REP).subs({FN[n]: JS[(n, 0)] for n in FN})
        z, left = canon(z)
        left = [a for a in left if a not in JETSYMS]
        assert not left, (k, left)
        used = [j for j in JETSYMS if z.has(j)]
        co = []
        for j in used:
            c = sp.diff(z, j)
            assert not any(c.has(q) for q in JETSYMS), (k, j)
            co.append((j, c))
        assert sp.expand(z.xreplace({j: 0 for j in JETSYMS})) == 0, k
        inv = {v: kk for kk, v in JS.items()}
        eqs[k] = [(inv[j], c) for j, c in co]
        assert all(nk[0] != 'H2' for nk, _ in eqs[k]), k
    return eqs


def bg_flow(mut_bg):
    """First and second derivative of alpha, a, phi0, u along the background, as expressions in Q."""
    parts = []
    for ex in (lc.dnu0, lc.dlam0, lc.ddp0):
        z, left = canon(ex)
        assert not left, left
        parts.append(z)
    dnu, dlam, ddp = parts
    u1 = sp.expand(ddp / Q['phi0']) - Q['u'] ** 2
    assert not any(p.base == Q['phi0'] and p.exp.is_negative for p in u1.atoms(sp.Pow)), "phi0 in a denominator"
    V1 = {'alpha': Q['alpha'] * dnu / 2, 'a': Q['a'] * dlam / 2, 'phi0': Q['u'] * Q['phi0'], 'u': u1}
    if mut_bg:
        V1[mut_bg] = -V1[mut_bg]

    def D(z):
        return sp.diff(z, Q['r']) + sum(sp.diff(z, Q[q]) * V1[q] for q in BGNAMES)
    V2 = {q: D(V1[q]) / 2 for q in BGNAMES}
    return V1, V2


# ------------------------------------------------------------------ the evaluation in K
class Ctx:
    def __init__(self, code, pt, flow, engine='frac', free=None):
        """pt: dict name -> Rational for the variables that are not free (None: all free). engine: 'frac' is the
        rational function field of sympy (gcd after every operation), 'af' the gcd-free AFField."""
        global KF
        self.pt, self.code, self.flow, self.engine = pt, code, flow, engine
        if free is None:
            free = NAMES if pt is None else ['PI']
        names = [n for n in NAMES if n == 'PI' or n in free]
        if engine == 'af':
            self.K = AFField(names)
            gens = [AF(self.K, g, {}) for g in self.K.gens]
        else:
            self.K, *gens = field(','.join(names), QQ)
        KF = self.K
        _kc.clear()
        self.names = names
        self.gen = dict(zip(names, gens))
        self.xmap = {Q[n]: v for n, v in (pt or {}).items() if n not in names}

    def conv(self, expr):
        return self.K.from_expr(expr.xreplace(self.xmap) if self.xmap else expr)

    def val(self, name):
        return self.gen[name] if name in self.gen else self.conv(Q[name])

    def evaluate(self, verbose=False):
        """J, J', h, h', h'' as elements of K."""
        V1, V2 = self.flow
        K = self.K
        cj = lambda name: Jet([self.val(name), K.zero, K.zero])
        inp = {'r': Jet([self.val('r'), K.one, K.zero])}
        for q in BGNAMES:
            inp[q] = Jet([self.val(q), self.conv(V1[q]), self.conv(V2[q])])
        for n in ('omega', 'epsilon', 's', 'LL', 'Lambda', 'mu'):
            inp[n] = cj(n)
        ns = dict(np=NPShim, pi=cj('PI'), exp=sp.exp, sqrt=sp.sqrt, log=sp.log,
                  Rational=lambda p, q: Jet([K.from_expr(sp.Rational(p, q)), K.zero, K.zero]))
        exec(self.code, ns)
        Jg, hg = ns['J_and_h'](*[inp[n] for n in ('r', 'alpha', 'a', 'phi0', 'u', 'omega', 'epsilon', 's', 'LL',
                                                     'Lambda', 'mu')])
        if verbose:
            say("  program run on the series (lazy): %s" % tm())
        J = [[lift(Jg[i, j])[0] for j in range(6)] for i in range(6)]
        if verbose:
            say("  J:    %s  %s" % (tm(), size(J)))
        Jp = [[lift(Jg[i, j])[1] for j in range(6)] for i in range(6)]
        if verbose:
            say("  J':   %s  %s" % (tm(), size(Jp)))
        h = [lift(x)[0] for x in hg]
        hp = [lift(x)[1] for x in hg]
        if verbose:
            say("  h, h': %s" % tm())
        hpp = [2 * lift(x)[2] for x in hg]
        if verbose:
            say("  h'':  %s  %s" % (tm(), size([hpp])))
        return dict(J=J, Jp=Jp, h=h, hp=hp, hpp=hpp)

    def make_rows(self, base, coef=None, verbose=False):
        """Rows multiplying y in each jet (field, order).  coef = (i, j): J[i, j] and J'[i, j] scaled by 8/7."""
        K = self.K
        J = [row[:] for row in base['J']]; Jp = [row[:] for row in base['Jp']]
        h, hp, hpp = base['h'], base['hp'], base['hpp']
        if coef:
            i, j = coef
            c = K.from_expr(SCALE)
            J[i][j] = J[i][j] * c; Jp[i][j] = Jp[i][j] * c
        s = self.val('s')
        rows = {}
        for n, i, sc in (('H1', 0, True), ('K', 1, False), ('F', 2, False), ('G', 4, True)):
            e = [K.one if m == i else K.zero for m in range(6)]
            f = (lambda v: [s * x for x in v]) if sc else (lambda v: v)
            rows[(n, 0)], rows[(n, 1)] = f(e), f(J[i])
            rows[(n, 2)] = f(vadd(Jp[i], vecmat(J[i], J)))
        hJ = vecmat(h, J)
        rows[('H0', 0)] = h
        rows[('H0', 1)] = vadd(hp, hJ)
        rows[('H0', 2)] = vadd(vadd(hpp, [2 * x for x in vecmat(hp, J)]), vadd(vecmat(h, Jp), vecmat(hJ, J)))
        if verbose:
            say("  rows of y, y', y'' for H0, H1, K, F, G: %s" % tm())
        return rows

    def residual(self, rows, co):
        """Six coefficients of y in the residual of one equation; co = [((name, order), coefficient)]."""
        res = [self.K.zero] * 6
        for nk, A in co:
            Ak = self.conv(A)
            for j in range(6):
                if rows[nk][j]:
                    res[j] = res[j] + Ak * rows[nk][j]
        return res


def vecmat(v, M):
    out = []
    for j in range(6):
        t = None
        for i in range(6):
            if v[i] and M[i][j]:
                p = v[i] * M[i][j]
                t = p if t is None else t + p
        out.append(KF.zero if t is None else t)
    return out


def vadd(a, b):
    return [x + y for x, y in zip(a, b)]


def sz(x):
    return (len(x.N), len(x.den)) if isinstance(x, AF) else (len(x.numer), len(x.denom))


def size(M):
    n = [sz(x) for row in M for x in row if x]
    return "terms of numerator / denominator%s: max %d / %d, total %d" % (
        " (atoms)" if n and isinstance(KF, AFField) else "", max((a for a, _ in n), default=0),
        max((b for _, b in n), default=0), sum(a + b for a, b in n))


def nnz(vec):
    return sum(1 for x in vec if x)


# ------------------------------------------------------------------ points
def make_points(n, seed):
    rng = random.Random(seed)
    q = lambda lo, hi, dlo, dhi: sp.Rational(rng.randint(lo, hi), rng.randint(dlo, dhi))
    pts = []
    for m in range(n):
        while True:
            u = q(-30, 30, 2, 13)
            if u != 0:
                break
        l = rng.randint(2, 8)
        pts.append(dict(r=q(1, 60, 3, 11), alpha=q(2, 40, 3, 13), a=q(2, 40, 3, 13), phi0=q(1, 30, 5, 31), u=u,
                        omega=q(1, 30, 5, 31), epsilon=sp.Integer(1 if m % 2 == 0 else -1),
                        s=q(1, 40, 2, 17) * rng.choice((1, -1)), LL=sp.Integer(l * (l + 1)), Lambda=q(1, 30, 3, 20),
                        mu=q(1, 30, 3, 20)))
    return pts


def run_points(a, code, eqs, order, coef, mut_bg, engine):
    flow = bg_flow(mut_bg)
    pts = make_points(a.npoints, a.seed)
    say("points (seed %d), epsilon alternates +1, -1:" % a.seed)
    for m, p in enumerate(pts):
        say("  %2d  " % (m + 1) + " ".join("%s=%s" % (k, p[k]) for k in NAMES[:-1]))
    rng = random.Random(a.seed + 1)
    stats = {k: dict(zero=0, yzero=0, nnz=[]) for k in order}
    anytgt = []
    for m, p in enumerate(pts):
        ctx = Ctx(code, p, flow, engine)
        base = ctx.evaluate()
        rows = ctx.make_rows(base, coef)
        res = {k: ctx.residual(rows, eqs[k]) for k in order}
        y = [sp.Rational(rng.randint(-9, 9), rng.randint(1, 7)) for _ in range(6)]
        line = []
        for k in order:
            c = res[k]
            ryv = sum((ctx.K.from_expr(y[j]) * c[j] for j in range(6)), ctx.K.zero)
            ryz = (not ryv) if engine == 'af' else (ryv.as_expr() == 0 and isinstance(ryv.as_expr(), sp.Rational))
            n = nnz(c)
            stats[k]['zero'] += n == 0; stats[k]['yzero'] += ryz
            stats[k]['nnz'].append(n)
            line.append("%s %s" % (k, "0" if n == 0 else "NONZERO(%d/6)" % n))
        anytgt.append(any(nnz(res[k]) for k in EQS_TARGET if k in res))
        say("point %2d  %s   [%s]" % (m + 1, "  ".join(line), tm()))
    ok = True
    say("\nsummary: residual vector (six coefficients of y) exactly zero, and residual at a random rational y exactly 0")
    for k in order:
        st = stats[k]
        z = st['zero'] == len(pts) and st['yzero'] == len(pts)
        say("  %-7s zero at %2d/%d points, residual at random y exactly 0 at %2d/%d   -> %s"
            % (k, st['zero'], len(pts), st['yzero'], len(pts), "PASS" if z else "FAIL (nonzero at %d points)"
               % (len(pts) - st['zero'])))
        ok &= z
    return ok, anytgt


def sensitivity(a, code, eqs, order, mut_bg):
    flow = bg_flow(mut_bg)
    p = make_points(1, a.seed)[0]
    ctx = Ctx(code, p, flow)
    base = ctx.evaluate()
    say("sensitivity at point 1: scale J[i, j] by 8/7 and recompute the residual")
    say("  rows i = 0..5, columns j = 0..5; per equation, '#' nonzero, '.' zero, '-' J[i, j] is 0")
    for k in order:
        say("  equation %s" % k)
        for i in range(6):
            row = ""
            for j in range(6):
                if not base['J'][i][j]:
                    row += " -"
                    continue
                res = ctx.residual(ctx.make_rows(base, (i, j)), eqs[k])
                row += " #" if nnz(res) else " ."
            say("    " + row)


ARGS = ('r', 'alpha', 'a', 'phi0', 'u', 'omega', 'epsilon', 's', 'LL', 'Lambda', 'mu')


def selftest(a, code, mut_bg):
    """The series derivatives J', h', h'' against sympy's chain rule on the symbolic J, h (same program and same
    first background derivatives, different differentiation), at two points, in K."""
    V1, _ = flow = bg_flow(mut_bg)
    ns = dict(np=NPShim, pi=Q['PI'], exp=sp.exp, sqrt=sp.sqrt, log=sp.log, Rational=sp.Rational)
    exec(code, ns)
    Jg, hg = ns['J_and_h'](*[Q[n] for n in ARGS])
    D = lambda z: sp.diff(z, Q['r']) + sum(sp.diff(z, Q[q]) * V1[q] for q in BGNAMES)
    say("selftest: symbolic J, h from the program with sympy symbols; D = d/dr + sum (q' d/dq), q' = the flow")
    ok = True
    t = time.time()
    Jd = [[D(Jg[i, j]) if not isinstance(Jg[i, j], int) else 0 for j in range(6)] for i in range(6)]
    hd = [D(x) for x in hg]; hdd = [D(x) for x in hd]
    say("  symbolic derivatives built %.0fs" % (time.time() - t))
    for m, p in enumerate(make_points(2, a.seed)):
        ctx = Ctx(code, p, flow)
        base = ctx.evaluate()
        cv = lambda z: ctx.conv(sp.sympify(z))
        cmp = {'J': (Jg, base['J'], 0), "J'": (Jd, base['Jp'], 0)}
        res = []
        for name, (S, B, _) in cmp.items():
            n = sum(cv(S[i][j] if name == "J'" else S[i, j]) == B[i][j] for i in range(6) for j in range(6))
            res.append((name, n, 36))
        for name, S, B in (('h', hg, base['h']), ("h'", hd, base['hp']), ("h''", hdd, base['hpp'])):
            res.append((name, sum(cv(S[j]) == B[j] for j in range(6)), 6))
        say("  point %d: " % (m + 1) + "  ".join("%s %d/%d equal" % x for x in res) + "   [%s]" % tm())
        ok &= all(n == tot for _, n, tot in res)
    say("RESULT selftest %s" % ("PASS: series and sympy chain rule agree exactly" if ok else "FAILED"))
    return 0 if ok else 1


def eval_af(x, F, vals):
    """Value of an AF element at rational values (list in the order of F.names), as a Fraction."""
    def poly(N):
        tot = Fraction(0)
        for ev, c in N.items():
            v = Fraction(int(c.numerator), int(c.denominator))
            for z, e in zip(vals, ev):
                if e:
                    v *= z ** e
            tot += v
        return tot
    v = poly(x.N) if x.N else Fraction(0)
    for k, e in x.den.items():
        v /= poly(F.atoms[k]) ** e
    return v


def crosscheck(a, code, eqs, order, mut_bg):
    """The symbolic objects of the gcd-free engine (all twelve variables free) evaluated at random rational points,
    against the same quantities from the field engine at those points (different arithmetic)."""
    flow = bg_flow(mut_bg)
    S = Ctx(code, None, flow, 'af')
    baseS = S.evaluate()
    rowsS = S.make_rows(baseS)
    coefS = {k: [(nk, S.conv(A)) for nk, A in eqs[k]] for k in order}
    say("crosscheck: symbolic objects built %s; atoms %d" % (tm(), len(S.K.atoms)))
    ok = True
    rng = random.Random(a.seed + 2)
    for m, pt in enumerate(make_points(a.npoints, a.seed + 2)):
        pv = sp.Rational(rng.randint(1, 60), rng.randint(1, 20))
        N = Ctx(code, pt, flow, 'frac')
        baseN = N.evaluate()
        rowsN = N.make_rows(baseN)
        vals = [Fraction(int(pt[n].p), int(pt[n].q)) if n != 'PI' else Fraction(int(pv.p), int(pv.q)) for n in NAMES]
        num = lambda x: x.as_expr().xreplace({Q['PI']: pv})
        eq = lambda xs, xn: (Fraction(int(num(xn).p), int(num(xn).q)) == eval_af(xs, S.K, vals))
        res = []
        for name in ('J', 'Jp'):
            res.append((name, sum(eq(baseS[name][i][j], baseN[name][i][j]) for i in range(6) for j in range(6)), 36))
        for name in ('h', 'hp', 'hpp'):
            res.append((name, sum(eq(baseS[name][j], baseN[name][j]) for j in range(6)), 6))
        res.append(('rows', sum(eq(rowsS[nk][j], rowsN[nk][j]) for nk in rowsS for j in range(6)), 6 * len(rowsS)))
        nc = tot = 0
        for k in order:
            for nk, A in coefS[k]:
                An = N.conv(dict(eqs[k])[nk])
                nc += eq(A, An); tot += 1
        res.append(('coefficients', nc, tot))
        say("  point %d (PI = %s): " % (m + 1, pv) + "  ".join("%s %d/%d" % x for x in res) + "   [%s]" % tm())
        ok &= all(n == tt for _, n, tt in res)
    say("RESULT crosscheck %s" % ("PASS: symbolic objects agree with the field engine at every point" if ok else "FAILED"))
    return 0 if ok else 1


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--method', default='points', choices=('points', 'symbolic'))
    ap.add_argument('--npoints', type=int, default=20)
    ap.add_argument('--seed', type=int, default=8)
    ap.add_argument('--mutate', default=None)
    ap.add_argument('--sensitivity', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--crosscheck', action='store_true')
    ap.add_argument('--deadline', type=int, default=0)
    ap.add_argument('--eqs', default=None)
    ap.add_argument('--engine', default=None, choices=('frac', 'af'))
    ap.add_argument('--free', default='all')
    a = ap.parse_args()
    order = tuple(a.eqs.split(',')) if a.eqs else EQS_PLANTED + EQS_TARGET
    coef = mut_bg = mut_eq = None
    if a.mutate:
        kind, _, arg = a.mutate.partition(':')
        if kind == 'coef':
            coef = tuple(int(x) for x in arg.split(',')) if arg else (0, 0)
        elif kind == 'bg':
            mut_bg = arg or 'alpha'
            assert mut_bg in BGNAMES
        elif kind == 'eq':
            mut_eq = True
        else:
            raise SystemExit("--mutate=coef[:i,j] or bg[:alpha|a|phi0|u] or eq")
    git = subprocess.run(['git', '-C', HERE, 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
    code, info = load_program()
    eqs = load_equations(order)
    if mut_eq:
        for key, nk in (('tt', ('H0', 2)), ('ang_tr', ('K', 2))):
            if key in eqs:
                eqs[key] = [(k, A * SCALE if k == nk else A) for k, A in eqs[key]]
    say("verify_polar_exact.py  %s" % " ".join(sys.argv[1:]))
    say("commit %s  python %s  sympy %s" % (git, sys.version.split()[0], sp.__version__))
    say("polar_matrix_generated.py sha256 %s: %d float exponents turned into ints, %d int/int quotients into Rational; "
        "exp, sqrt, log used in the program: %s" % (info['sha256'], info['floats'], info['ratios'],
                                                      info['exp_sqrt_log_used'] or 'none'))
    for k in order:
        say("  E[%-6s] jets: %s" % (k, ", ".join("%s%s" % (n, "'" * kk) for (n, kk), _ in sorted(eqs[k]))))
    engine = a.engine or ('frac' if a.method == 'points' else 'af')
    say("engine: %s" % ({'frac': 'sympy rational function field (gcd after every operation)',
                         'af': 'gcd-free: Laurent polynomial numerators, inverted polynomials as denominator atoms'}[engine]))
    if a.mutate:
        say("MUTATION %s" % a.mutate)
    if a.selftest:
        return selftest(a, code, mut_bg)
    if a.crosscheck:
        return crosscheck(a, code, eqs, order, mut_bg)
    if a.sensitivity:
        sensitivity(a, code, eqs, order, mut_bg)
        return 0
    if a.method == 'points':
        ok, anytgt = run_points(a, code, eqs, order, coef, mut_bg, engine)
        if a.mutate:
            det = all(anytgt)
            say("\nMUTATION %s: (tt) or the angular trace has a nonzero residual at %d/%d points -> %s"
                % (a.mutate, sum(anytgt), len(anytgt), "DETECTED" if det else "NOT DETECTED"))
            return 1 if det else 0
        say("\nmethod: exact evaluation at %d rational points (an identity test at points, not a symbolic proof)"
            % a.npoints)
        say("RESULT %s" % ("all eight residuals exactly 0 at all points" if ok else "FAILED"))
        return 0 if ok else 1
    # symbolic
    free = None if a.free == 'all' else a.free.split(',')
    assert free is None or all(x in NAMES for x in free), free
    full = free is None
    if a.deadline:
        def stop(*_):
            raise TimeoutError
        signal.signal(signal.SIGALRM, stop); signal.alarm(a.deadline)
    flow = bg_flow(mut_bg)
    pts = [None] if full else make_points(a.npoints, a.seed)
    say("symbolic: %s free, %s, deadline %s" % (
        "all twelve variables" if full else ", ".join(sorted(set(free) | {'PI'}, key=NAMES.index)),
        "no base point" if full else "the others at the %d points of seed %d" % (len(pts), a.seed),
        ("%d s" % a.deadline) if a.deadline else "none"))
    ok, done, detected = True, [], []
    stage = 'program'
    try:
        for m, pt in enumerate(pts):
            if not full:
                say("base point %d: %s" % (m + 1, " ".join("%s=%s" % (k, pt[k]) for k in NAMES[:-1] if k not in free)))
            ctx = Ctx(code, pt, flow, engine, free)
            stage = 'program'
            base = ctx.evaluate(verbose=True)
            if engine == 'af':
                say("  atoms (inverted polynomials): %d" % len(ctx.K.atoms))
                for k, A in enumerate(ctx.K.atoms):
                    say("    atom %d (%d terms, monic): %s" % (k, len(A), A.as_expr()))
            stage = 'rows'
            rows = ctx.make_rows(base, coef, verbose=True)
            ns = []
            for k in order:
                stage = k
                t = time.time()
                res = ctx.residual(rows, eqs[k])
                n = nnz(res)
                say("  %-7s residual %s  (%d of 6 coefficients nonzero)  %.0fs, total %s"
                    % (k, "= 0 exactly" if n == 0 else "NONZERO", n, time.time() - t, tm()))
                ns.append((k, n))
                ok &= n == 0
            done.append(ns)
            detected.append(any(n for k, n in ns if k in EQS_TARGET))
    except TimeoutError:
        say("\nNOT FINISHED: deadline of %d s reached in stage '%s'; finished: %s" % (a.deadline, stage, done))
        return 2
    if a.mutate:
        say("MUTATION %s: (tt) or the angular trace has a nonzero residual at %d/%d base points -> %s"
            % (a.mutate, sum(detected), len(detected), "DETECTED" if all(detected) else "NOT DETECTED"))
        return 1 if all(detected) else 0
    if full:
        say("\nmethod: full symbolic computation, exact, all twelve variables free (%s engine)" % engine)
        say("RESULT %s" % ("all residuals are exactly zero: checked symbolically" if ok else "FAILED"))
    else:
        say("\nmethod: exact, free variables %s, the others at the base points: an identity test in the other variables, "
            "not a full symbolic proof" % ",".join(sorted(set(free) | {'PI'}, key=NAMES.index)))
        say("RESULT %s" % ("all residuals exactly zero at the base points" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
