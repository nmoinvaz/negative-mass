"""Shared formulas and harness of #142's review of RESULTS.md (Lemma A(iii), L2-L6, the Zerilli V).

Written from phase2/dg_proof/RESULTS.md alone. What this shares with #2: the formulas typed from RESULTS.md
(c_n of Lemma A, V of 4.2, f, x, psi_0, g, V-hat, Psi_DG, k). Nothing else.

Harness: every check prints PASS <id> or FAIL <id>: <why>. Exit 0: all pass. Exit 1: a real check failed.
Exit 2: a planted check failed (the script stops before any real data). --mutate=NAME[,NAME] switches on a
deliberate error; mutations are listed next to the code they change (grep "mut(").
"""
import sys
import sympy as sp
from sympy import Rational as Q, Symbol, symbols, Poly, factorial, gamma, cancel, together, simplify, expand, factor

r = Symbol('r', positive=True)
l = Symbol('l', positive=True)

_MUT = set()
for _a in sys.argv[1:]:
    if _a.startswith('--mutate='):
        _MUT |= set(_a[len('--mutate='):].split(','))


def mut(name):
    return name in _MUT


class Run:
    """PASS/FAIL harness. planted(): a failure exits 2 at once. check(): failures counted, exit 1 at finish()."""

    def __init__(self, title):
        self.nfail = 0
        self.npass = 0
        self.ndiff = 0
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

    def published(self, cid, cond, why=''):
        """Comparison with a formula printed in a source: a difference is reported (DIFFERS), not counted as a failure."""
        if cond:
            print('PASS', cid, flush=True)
            self.npass += 1
        else:
            print('DIFFERS', cid + ':', why, flush=True)
            self.ndiff += 1

    def note(self, text):
        print('#', text, flush=True)

    def finish(self):
        print('# summary: %d PASS, %d FAIL, %d DIFFERS' % (self.npass, self.nfail, self.ndiff), flush=True)
        sys.exit(1 if self.nfail else 0)


def is_zero(expr):
    """Exactly zero after cancel (rational functions) or expand."""
    e = sp.simplify(expr) if not expr.is_rational_function() else cancel(together(expr))
    return e == 0


# ---------------------------------------------------------------- notation (RESULTS.md, line 3)
def c_of(lv):
    return (lv + 2) * (lv - 1)


def L_of(lv):
    return lv * (lv + 1)


def K_of(lv):
    return (lv - 1) * lv * (lv + 1) * (lv + 2)


def k_of(lv):
    return sp.Rational(1, 6) * c_of(lv) * (c_of(lv) + 2)   # never Python float division


def rs_of(lv):
    return Q(6) / c_of(lv)


def rf(x, m):
    """Rising factorial (x)_m for integer m >= 0; for m = -2 the convention 1/((x-1)(x-2)) (Gamma ratio)."""
    if m == -2:
        return sp.Integer(1) / ((x - 1) * (x - 2))
    out = 1
    for j in range(m):
        out = out * (x + j)
    return out


def cn(n, lv):
    """c_n of Lemma A: (n-2)((n-4)l(l+1)+n-1)(l-n+3)_{2n-4}/(2^n (n-1)!^2). lv may be a symbol (n numeric)."""
    n4 = n - 3 if mut('cn_n4to3') else n - 4   # mutation cn_n4to3: (n-4) -> (n-3)
    num = sp.sympify((n - 2) * (n4 * lv * (lv + 1) + n - 1)) * rf(lv - n + 3, 2 * n - 4)
    return sp.sympify(num) / (2 ** n * factorial(n - 1) ** 2)


def chi_poly(lv):
    """chi as a Poly in r over QQ, l an integer >= 2."""
    return Poly(sum(cn(n, lv) * r ** n for n in range(1, lv + 3)), r, domain='QQ')


def chi_over_r(lv):
    return Poly(sum(cn(n, lv) * r ** (n - 1) for n in range(1, lv + 3)), r, domain='QQ')


# ---------------------------------------------------------------- Sturm count on (0, oo)
def positive_root_count(poly):
    """Number of distinct roots in (0, oo) of a Poly with poly(0) != 0, by the Sturm sequence:
    V(0) - V(oo), V = sign changes of the sequence evaluated at 0 and at +oo (leading coefficients)."""
    p = Poly(poly, r, domain='QQ')
    assert p.eval(0) != 0
    seq = sp.sturm(p)
    at0 = [s.eval(0) for s in seq]
    atinf = [s.LC() for s in seq]
    if mut('pl_sturm_sign'):   # mutation pl_sturm_sign: evaluate "at infinity" with the sign at -infinity
        atinf = [s.LC() * (-1) ** s.degree() for s in seq]
    return sign_changes(at0) - sign_changes(atinf)


def sign_changes(vals):
    vs = [v for v in vals if v != 0]
    return sum(1 for a, b in zip(vs, vs[1:]) if (a > 0) != (b > 0))


# ---------------------------------------------------------------- f, x, derivatives in x
f_of_r = 1 + 2 / r


def dx(expr):
    """d/dx = f d/dr, with f = 1 + 2/r (RESULTS line 3: f = 1 + 2/r, x = r - 2 ln(1 + r/2), dx/dr = 1/f)."""
    return f_of_r * sp.diff(expr, r)


def V_results(lam):
    """V of RESULTS 4.2 (static Zerilli potential, M = -1): 2f[lam^2 r^2((lam+1) r - 3) + 9(lam r - 1)]/(r^3 (lam r - 3)^2).
    Mutation v_typo: 9(lam r - 2)."""
    nine = 9 * (lam * r - (2 if mut('v_typo') else 1))
    return 2 * f_of_r * (lam ** 2 * r ** 2 * ((lam + 1) * r - 3) + nine) / (r ** 3 * (lam * r - 3) ** 2)


def Vhat_of_u(u):
    """V-hat = g^2 - g_x with g = f u, u = psi_r/psi. Mutation pl_vhat: g^2 + g_x (the planted checks must fire)."""
    g = f_of_r * u
    return g ** 2 + dx(g) if mut('pl_vhat') else g ** 2 - dx(g)


def Vhat_of_psi(psi):
    """V-hat = g^2 - g_x, g = psi_x/psi (RESULTS line 4)."""
    return Vhat_of_u(sp.diff(psi, r) / psi)


# ---------------------------------------------------------------- series, exponents, planted-ODE builder
E = Symbol('E')
def lead(expr, where):
    """(power, coeff) with expr ~ coeff * r**power as r -> 0 (where = 0) or r -> oo (where = 'inf'); expr rational in r."""
    num, den = sp.fraction(cancel(together(expr)))
    pn, pd = Poly(num, r), Poly(den, r)
    if where == 0:
        mn = min(m[0] for m in pn.monoms())
        md = min(m[0] for m in pd.monoms())
        cn_ = pn.coeff_monomial(r ** mn)
        cd_ = pd.coeff_monomial(r ** md)
        return mn - md + (1 if mut('pl_order') else 0), cancel(cn_ / cd_)   # mutation pl_order: power off by one
    return pn.degree() - pd.degree(), cancel(pn.LC() / pd.LC())


def taylor(num, den, n):
    """First n Taylor coefficients at r = 0 of num/den (Polys in r), den(0) != 0."""
    nc = num.all_coeffs()[::-1] + [0] * n
    dc = den.all_coeffs()[::-1]
    assert dc[0] != 0
    out = []
    for i in range(n):
        v = nc[i] - sum(dc[j] * out[i - j] for j in range(1, min(i, len(dc) - 1) + 1))
        out.append(cancel(v / dc[0]))
    return out


def indicial(P2, P1, P0, x0):
    """Exponents at the regular singular point x0 of P2 u'' + P1 u' + P0 u = 0 (rational coefficients in r):
    s(s-1) + p0 s + q0 = 0 with p0 = lim (r-x0) P1/P2, q0 = lim (r-x0)^2 P0/P2."""
    t = r - x0
    A = cancel(together(t * P1 / P2))
    B = cancel(together(t ** 2 * P0 / P2))
    assert sp.fraction(A)[1].subs(r, x0) != 0 and sp.fraction(B)[1].subs(r, x0) != 0, 'not a regular singular point'
    p0, q0 = A.subs(r, x0), B.subs(r, x0)
    if mut('pl_indicial'):   # mutation pl_indicial: wrong sign of q0
        q0 = -q0
    ss = Symbol('ss')
    return set(sp.solve(ss * (ss - 1) + p0 * ss + q0, ss))


def wronskian_ode(y1, y2):
    """P2, P1, P0 of the ODE with solutions y1, y2: det [[u,u',u''],[y1,y1',y1''],[y2,y2',y2'']] = 0."""
    d1, d2 = sp.diff(y1, r), sp.diff(y2, r)
    e1, e2 = sp.diff(y1, r, 2), sp.diff(y2, r, 2)
    return (y1 * d2 - d1 * y2, -(y1 * e2 - e1 * y2), d1 * e2 - e1 * d2)
