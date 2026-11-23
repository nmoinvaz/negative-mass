#!/usr/bin/env python3
"""Shared code of the independent review of RESULTS.md section 4 (issue 104).

Written from RESULTS.md alone: it shares with the work under review only the formulas typed from that file.
Certified parts use mpmath `iv` (outward-rounded intervals); identities use sympy.

  rho = c r, eps = 1/c, c = (l+2)(l-1), X(rho) = sum_{n=1}^{l+2} a_n(eps) rho^n,
  a_1 = 3/2, a_2 = 0, a_n = (n-2)[(n-4)+3(n-3)eps]/(2^n (n-1)!^2) prod_{a=0}^{n-3}(1-(a-1)(a+2)eps),
  N = X - (rho-6)X', N_Z = (1+2eps)rho^3/8 - 3rho^2/4 + 9rho/2 - 9,
  Psi = rho^2 (2+eps rho) N^2 - 4 N_Z X^2,
  S = rho X [eps N + (2+eps rho) N'] - (2+eps rho) N (X + 2 rho X').

Running this file self-tests the module (exit 0).  `--mutate=NAME` breaks one ingredient (exit non-zero).
"""
import sys
from fractions import Fraction
from math import factorial

from mpmath import iv, mp

iv.prec = 192
mp.prec = 192

MUT = None  # set by the scripts from --mutate=NAME


def set_mutation(name):
    global MUT
    MUT = name


def parse_args(argv, known_mutations):
    """Return (mutation or None, remaining flags).  An unknown mutation name is an error (exit 3)."""
    mut = None
    rest = []
    for a in argv[1:]:
        if a.startswith("--mutate="):
            mut = a.split("=", 1)[1]
        else:
            rest.append(a)
    if mut is not None and mut not in known_mutations:
        print(f"unknown mutation {mut!r}; known: {', '.join(known_mutations)}")
        sys.exit(3)
    set_mutation(mut)
    return mut, rest


class Checker:
    """Prints PASS/FAIL lines.  With a mutation active it stops at the first failure (exit 1)."""

    def __init__(self, failfast=False):
        self.failfast = failfast
        self.fails = []
        self.count = 0

    def check(self, name, cond, detail="", quiet=False):
        self.count += 1
        ok = bool(cond)
        if not (quiet and ok):
            print(f"{'PASS' if ok else 'FAIL'} {name} {detail}".rstrip(), flush=True)
        if not ok:
            self.fails.append(name)
            if self.failfast:
                print(f"STOP first failure under mutation: {name}")
                sys.exit(1)
        return ok

    def finish(self):
        print(f"SUMMARY checks={self.count} failed={len(self.fails)} {self.fails if self.fails else ''}".rstrip())
        sys.exit(1 if self.fails else 0)


# ----------------------------------------------------------------------------------------------------------
# numbers
# ----------------------------------------------------------------------------------------------------------
def c_of(l):
    return (l + 2) * (l - 1)


def ivint(k):
    return iv.mpf(k)


def endpoint(num, den=1):
    """iv enclosure of the rational num/den (exact integers)."""
    if den == 1:
        return iv.mpf(num)
    return iv.mpf(num) / iv.mpf(den)


def sgn(v):
    """+1 / -1 if the enclosure excludes 0, else 0."""
    if v.a > 0:
        return 1
    if v.b < 0:
        return -1
    return 0


def hull(u, v):
    return iv.mpf([min(u.a, v.a), max(u.b, v.b)])


def lo_of(v):
    return v.a


def hi_of(v):
    return v.b


# ----------------------------------------------------------------------------------------------------------
# coefficients a_n(eps), finite l
# ----------------------------------------------------------------------------------------------------------
def q_of(a):
    """q_a = (a-1)(a+2); the mutation 'coef_q' changes it to (a-1)(a+3)."""
    if MUT == "coef_q":
        return (a - 1) * (a + 3)
    return (a - 1) * (a + 2)


def coeffs(l, K=iv.mpf):
    """[a_0 .. a_{l+2}] as K-numbers (a_0 = 0), by the product recurrence, with eps an enclosure of 1/c."""
    c = c_of(l)
    eps = K(1) / K(c)
    lin = 3 if MUT != "coef_lin" else 2
    out = [K(0), K(3) / K(2), K(0)]
    # w_n = prod_{a=0}^{n-3}(1 - q_a eps) / (2^n (n-1)!^2), starting at n = 3
    w = (1 + 2 * eps) / K(32)
    for n in range(3, l + 3):
        out.append(K(n - 2) * (K(n - 4) + lin * K(n - 3) * eps) * w)
        # w_{n+1} = w_n (1 - q_{n-2} eps) / (2 n^2)
        w = w * (1 - K(q_of(n - 2)) * eps) / K(2 * n * n)
    if MUT == "drop_last":
        out[-1] = K(0)
    return out


def a_exact(l, n):
    """Lemma A: c * c_n / c^n with c_n = (n-2)((n-4) l(l+1) + n-1) (l-n+3)_{2n-4} / (2^n (n-1)!^2), exact Fraction."""
    c = c_of(l)
    L = l * (l + 1)
    if n == 1:
        cn = Fraction(3, 2)
    else:
        rf = 1
        for j in range(2 * n - 4):
            rf *= (l - n + 3 + j)
        cn = Fraction((n - 2) * ((n - 4) * L + n - 1) * rf, 2 ** n * factorial(n - 1) ** 2)
    return cn * c / Fraction(c) ** n


def frac_in_iv(fr, v):
    """True if the exact Fraction lies in the enclosure v."""
    return to_frac(v.a) <= fr <= to_frac(v.b)


def to_frac(x):
    """Exact Fraction of the left endpoint of a point enclosure (raw mpmath tuple: sign, mantissa, exponent, bit count)."""
    sign, man, exp, _bc = x._mpi_[0]
    v = Fraction(int(man)) * (Fraction(2) ** exp)
    return -v if sign else v


# ----------------------------------------------------------------------------------------------------------
# jets of X: [X, X', X'', ...] at rho
# ----------------------------------------------------------------------------------------------------------
def jet_from_coeffs(a, rho, kmax=2):
    """sum_n a_n n^(falling k) rho^(n-k) for k = 0..kmax; rho an iv scalar (point or interval) or an mp number, rho > 0."""
    zero = rho * 0
    one = zero + 1
    S = [zero] * (kmax + 1)
    p = one  # rho^n
    for n in range(1, len(a)):
        p = p * rho
        if n == 2:
            continue
        t = a[n] * p
        S[0] = S[0] + t
        if kmax >= 1:
            u = t * n
            S[1] = S[1] + u
            if kmax >= 2:
                u = u * (n - 1)
                S[2] = S[2] + u
                if kmax >= 3:
                    u = u * (n - 2)
                    S[3] = S[3] + u
    out = []
    rk = one
    for k in range(kmax + 1):
        out.append(S[k] / rk)
        rk = rk * rho
    return out


def alpha(n, K=iv.mpf):
    """alpha_n = (n-2)(n-4)/(2^n (n-1)!^2), the eps = 0 coefficient (alpha_1 = 3/2, alpha_2 = 0).
    The mutation 'alpha_n' uses (n-2)(n-5)."""
    m = 4 if MUT != "alpha_n" else 5
    return K((n - 2) * (n - m)) / K(2 ** n * factorial(n - 1) ** 2)


def alpha_p(n, K=iv.mpf):
    """alpha'_n = 3 (n-2)(n-3)/(2^n (n-1)!^2)."""
    return K(3 * (n - 2) * (n - 3)) / K(2 ** n * factorial(n - 1) ** 2)


def falling(n, k):
    r = 1
    for j in range(k):
        r *= (n - j)
    return r


def jet_inf(rho, kmax=2, nmax=90):
    """Enclosures of X_inf^(k)(rho), k = 0..kmax, rho an iv scalar with 1 <= rho <= 12.
    Terms n <= nmax exactly; for n > nmax: |alpha_n| n^(falling k) rho^(n-k) <= n^(2+k) rho^(n-k) / (2^n (n-1)!^2) =: U(n),
    U(n+1)/U(n) = ((n+1)/n)^(2+k) rho/(2 n^2) <= r < 1 for n >= nmax, so the tail is at most U(nmax+1)/(1-r)."""
    assert rho.a >= 1 and rho.b <= 12
    out = []
    rhohi = iv.mpf(rho.b)
    for k in range(kmax + 1):
        s = iv.mpf(0)
        for n in range(1, nmax + 1):
            if n < k and falling(n, k) == 0:
                continue
            s = s + alpha(n) * falling(n, k) * rho ** (n - k)
        n0 = nmax
        U = iv.mpf((n0 + 1) ** (2 + k)) * rhohi ** (n0 + 1 - k) / iv.mpf(2 ** (n0 + 1) * factorial(n0) ** 2)
        r = (iv.mpf(n0 + 1) / n0) ** (2 + k) * rhohi / (2 * n0 * n0)
        assert r.b < 1
        T = U / (1 - r)
        out.append(s + iv.mpf([-T.b, T.b]))
    return out


def jet_inf_mp(rho, kmax=2, nmax=140):
    """Plain mp.mpf evaluation of the same series (no tail), for root guesses and reference values."""
    out = []
    for k in range(kmax + 1):
        s = mp.mpf(0)
        for n in range(1, nmax + 1):
            s += alpha(n, mp.mpf) * falling(n, k) * rho ** (n - k)
        out.append(s)
    return out


def jet_finite_mp(l, rho, kmax=2):
    a = coeffs(l, K=mp.mpf)
    S = [mp.mpf(0)] * (kmax + 1)
    for n in range(1, len(a)):
        for k in range(kmax + 1):
            S[k] += a[n] * falling(n, k) * rho ** (n - k)
    return S


# ----------------------------------------------------------------------------------------------------------
# N, Psi, S from a jet
# ----------------------------------------------------------------------------------------------------------
def NZ(rho, eps):
    """N_Z = (1+2eps) rho^3/8 - 3 rho^2/4 + 9 rho/2 - 9."""
    t = 3 if MUT != "nz_coef" else 4
    return (1 + 2 * eps) * rho ** 3 / 8 - t * rho ** 2 / 4 + 9 * rho / 2 - 9


def NZp(rho, eps):
    t = 3 if MUT != "nz_coef" else 4
    return 3 * (1 + 2 * eps) * rho ** 2 / 8 - t * rho / 2 + 4.5


def quantities(X, rho, eps):
    """X = [X, X', X'']; returns dict with N, N', Psi, Psi', S (same number type as X)."""
    N = X[0] - (rho - 6) * X[1]
    Np = -(rho - 6) * X[2]
    nz = NZ(rho, eps)
    nzp = NZp(rho, eps)
    Psi = rho ** 2 * (2 + eps * rho) * N ** 2 - 4 * nz * X[0] ** 2
    Psip = ((4 * rho + 3 * eps * rho ** 2) * N ** 2 + 2 * rho ** 2 * (2 + eps * rho) * N * Np
            - 4 * nzp * X[0] ** 2 - 8 * nz * X[0] * X[1])
    S = rho * X[0] * (eps * N + (2 + eps * rho) * Np) - (2 + eps * rho) * N * (X[0] + 2 * rho * X[1])
    if MUT == "S_flip":
        S = -S
    return dict(X=X[0], X1=X[1], X2=X[2], N=N, Np=Np, Psi=Psi, Psip=Psip, S=S)


# ----------------------------------------------------------------------------------------------------------
# bracket certification
# ----------------------------------------------------------------------------------------------------------
def sgn_checked(v):
    """Sign of an enclosure, +-1 only if it excludes 0.  Mutation 'accept_zero_enclosure' treats [0, b] as positive."""
    s = sgn(v)
    if s == 0 and MUT == "accept_zero_enclosure":
        if v.a >= 0:
            return 1
        if v.b <= 0:
            return -1
    return s


def certify_bracket(F, dF, lo, hi):
    """F: point -> iv enclosure; dF: interval -> iv enclosure of F' over it.  lo, hi: iv scalars (the two ends).
    Certified when F has definite opposite signs at the ends (enclosures exclude 0) and dF over [lo, hi] excludes 0, so that
    F has exactly one root in [lo, hi].  Returns (ok, reason)."""
    s1, s2 = sgn_checked(F(lo)), sgn_checked(F(hi))
    if s1 == 0 or s2 == 0:
        return False, "enclosure of F at an end contains 0"
    if s1 == s2:
        return False, "no sign change"
    if MUT != "no_deriv":
        d = sgn_checked(dF(hull(lo, hi)))
        if d == 0:
            return False, "derivative enclosure contains 0"
        if d != s2:
            return False, "derivative sign inconsistent with the sign change"
    return True, "certified"


# ----------------------------------------------------------------------------------------------------------
# planted problem
# ----------------------------------------------------------------------------------------------------------
def _poly_from_roots(roots):
    """Coefficients (high to low) of prod (r - root) with Fraction roots."""
    co = [Fraction(1)]
    for rt in roots:
        new = [Fraction(0)] * (len(co) + 1)
        for i, v in enumerate(co):
            new[i] += v
            new[i + 1] -= v * rt
        co = new
    return co


def _make_F(co):
    ivco = [endpoint(v.numerator, v.denominator) for v in co]
    dco = [ivco[i] * (len(co) - 1 - i) for i in range(len(co) - 1)]

    def horner(cs, x):
        s = iv.mpf(0)
        for v in cs:
            s = s * x + v
        return s

    return (lambda x: horner(ivco, x)), (lambda x: horner(dco, x))


def planted_cases():
    a = Fraction(17, 2)
    b = a + Fraction(1, 10 ** 10)
    mid = (a + b) / 2
    tiny = Fraction(1, 10 ** 14)
    base = [Fraction(3), a, b, Fraction(11)]
    P1 = _make_F(_poly_from_roots(base))
    P3 = _make_F(_poly_from_roots(base + [mid]))
    P7 = _make_F(_poly_from_roots(base + [a]))
    cases = [
        ("one root near a", P1, a - tiny, a + tiny, True, None),
        ("one root near b", P1, b - tiny, b + tiny, True, None),
        ("two roots a, b = a + 1e-10 in one bracket", P1, a - 100 * tiny, b + 100 * tiny, False, "no sign change"),
        ("three roots in one bracket", P3, a - 100 * tiny, b + 100 * tiny, False, "derivative enclosure contains 0"),
        ("endpoint exactly on a root", P1, a, a + tiny, False, "enclosure of F at an end contains 0"),
        ("no root in the bracket", P1, a + 10 ** 5 * tiny, a + 2 * 10 ** 5 * tiny, False, "no sign change"),
        ("double root at a", P7, a - tiny, a + tiny, False, "no sign change"),
    ]
    out = []
    for name, (F, dF), lo, hi, want, why in cases:
        out.append((name, F, dF, lo, hi, want, why))
    return out


def run_planted(chk, verbose=True):
    """Runs the planted problem through certify_bracket.  Returns True if every verdict and every reason is as expected."""
    allok = True
    for name, F, dF, lo, hi, want, why in planted_cases():
        elo = endpoint(lo.numerator, lo.denominator)
        ehi = endpoint(hi.numerator, hi.denominator)
        ok, reason = certify_bracket(F, dF, elo, ehi)
        good = (ok == want) and (want or reason == why)
        allok &= good
        chk.check(f"planted: {name}", good, f"verdict={'CERTIFIED' if ok else 'REFUSED'} ({reason})")
    return allok



# ----------------------------------------------------------------------------------------------------------
# helpers for the perturbation scripts (Theorems L and M)
# ----------------------------------------------------------------------------------------------------------
def tail_majorant(D, k, rho_hi, n1):
    """Upper bound (iv point) of sum_{n>n1} n^D rho_hi^(n-k) / (2^n (n-1)!^2), for rho_hi >= 1.
    U(n) = n^D rho_hi^(n-k) g_n with g_n = 1/(2^n (n-1)!^2), U(n+1)/U(n) = ((n+1)/n)^D rho_hi/(2 n^2) <= r for n >= n1, so the sum is <= U(n1+1)/(1-r)."""
    rho_hi = iv.mpf(rho_hi)
    U = iv.mpf(n1 + 1) ** D * rho_hi ** (n1 + 1 - k) / iv.mpf(2 ** (n1 + 1) * factorial(n1) ** 2)
    r = (iv.mpf(n1 + 1) / n1) ** D * rho_hi / (2 * n1 * n1)
    assert r.b < 1
    return iv.mpf((U / (1 - r)).b)


def subintervals(lo, hi, parts):
    """Cover [lo, hi] (iv scalars) by `parts` overlapping iv intervals."""
    out = []
    w = hi - lo
    for k in range(parts):
        a = lo + w * k / parts
        b = lo + w * (k + 1) / parts
        out.append(iv.mpf([a.a, b.b]))
    return out


def ub(x):
    """Upper endpoint as a point iv number."""
    return iv.mpf(x.b)


def absub(x):
    """Upper bound of |x| as a point iv number."""
    return iv.mpf(max(abs(x.a), abs(x.b)))

# ----------------------------------------------------------------------------------------------------------
# self test of this module
# ----------------------------------------------------------------------------------------------------------
KNOWN = ["coef_q", "coef_lin", "nz_coef", "accept_zero_enclosure", "no_deriv", "alpha_n", "drop_last", "S_flip"]


def main():
    mut, rest = parse_args(sys.argv, KNOWN)
    chk = Checker(failfast=mut is not None)
    print(f"review_large_l_common self test, mutation={mut}")
    # planted problem first
    if not run_planted(chk):
        print("PLANTED PROBLEM FAILED, stop")
        sys.exit(2)
    # recurrence coefficients against Lemma A's closed form (exact Fractions), l = 2..12 all n, larger l for n <= 40
    for l in list(range(2, 13)) + [30, 60, 150, 500, 3000]:
        a = coeffs(l)
        nmax = min(l + 2, 40) if l > 150 else l + 2
        ok = all(frac_in_iv(a_exact(l, n), a[n]) for n in range(1, nmax + 1))
        chk.check(f"coefficients a_n(eps) in Lemma A's closed form, l={l}, n<={nmax}", ok)
        chk.check(f"Lemma A: c_n = 0 for n = l+3..l+6 and a_{{l+2}} > 0, l={l}",
                  all(a_exact(l, n) == 0 for n in range(l + 3, l + 7)) and a[l + 2].a > 0)
    # a_n = 0 beyond l+2 (the product has the factor 1 - q_l eps = 0)
    for l in (2, 5, 10):
        c = c_of(l)
        chk.check(f"factor a = l vanishes: q_l = c, l={l}", q_of(l) == c)
    # eps = 0 limit: alpha_n = (n-2)(n-4)/(2^n (n-1)!^2)
    chk.check("alpha_1 = 3/2, alpha_2 = 0", alpha(1) == iv.mpf(3) / 2 and alpha(2) == 0)
    # N_Z(6 + tau) has positive coefficients: 18 + 9 tau + 3/2 tau^2 + tau^3/8 + 2 eps (27 + 27/2 tau + 9/4 tau^2 + tau^3/8)
    import sympy as sp
    tau, eps_s = sp.symbols("tau epsilon")
    nz = (1 + 2 * eps_s) * (6 + tau) ** 3 / 8 - sp.Rational(3, 4) * (6 + tau) ** 2 + sp.Rational(9, 2) * (6 + tau) - 9
    P = sp.Poly(sp.expand(nz), tau, eps_s)
    chk.check("N_Z(6+tau) polynomial in (tau, eps), all coefficients positive", all(co > 0 for _, co in P.terms()), str(P.terms()))
    # NZ code against sympy for a numeric point
    rr = iv.mpf(8) + iv.mpf(1) / 3
    e = iv.mpf(1) / 1234
    nzs = nz.subs({tau: sp.Rational(25, 3) - 6, eps_s: sp.Rational(1, 1234)})
    chk.check("NZ(rho, eps) matches sympy at a point", abs(NZ(rr, e).mid - iv.mpf(sp.Float(nzs, 40).__str__()).mid) < 1e-30)
    # Psi', S, N' formulas against sympy differentiation with a generic X
    rho_s, e_s = sp.symbols("rho e")
    Xf = sp.Function("X")
    Xs = Xf(rho_s)
    Ns = Xs - (rho_s - 6) * sp.diff(Xs, rho_s)
    nzs2 = (1 + 2 * e_s) * rho_s ** 3 / 8 - sp.Rational(3, 4) * rho_s ** 2 + sp.Rational(9, 2) * rho_s - 9
    Psis = rho_s ** 2 * (2 + e_s * rho_s) * Ns ** 2 - 4 * nzs2 * Xs ** 2
    d_sym = sp.diff(Psis, rho_s)
    x0, x1, x2 = sp.symbols("x0 x1 x2")
    d_sub = d_sym.subs(sp.Derivative(Xs, (rho_s, 2)), x2).subs(sp.Derivative(Xs, rho_s), x1).subs(Xs, x0)
    pt = {x0: sp.Rational(7, 3), x1: sp.Rational(-5, 7), x2: sp.Rational(11, 5), rho_s: sp.Rational(25, 3), e_s: sp.Rational(1, 100)}
    ref = d_sub.subs(pt)
    Xi = [iv.mpf(7) / 3, iv.mpf(-5) / 7, iv.mpf(11) / 5]
    qd = quantities(Xi, iv.mpf(25) / 3, iv.mpf(1) / 100)
    diff = abs(qd["Psip"].mid - iv.mpf(sp.Float(ref, 40).__str__()).mid)
    chk.check("Psi' formula equals sympy d/d rho (generic X)", diff < 1e-30, f"diff={float(diff):.1e}")
    chk.finish()


if __name__ == "__main__":
    main()
