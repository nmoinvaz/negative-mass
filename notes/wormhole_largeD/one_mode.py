#!/usr/bin/env python3
"""Checks for one_mode.md: exactly one unstable radial mode of the symmetric Ellis wormhole in every D.

    python3 one_mode.py            all checks; exit 1 if a check fails or a control passes

Notation (throat radius 1, n = D - 2, m = n - 1): r'^2 = 1 - r^(-2m), r'' = m r^(-2m-1) (l = proper distance), Paper
III's master equation (master.py, gauge delta phi = 0)
    E2[Y] = Y'' - (2 r''/r' + n r'/r) Y' + (n r'^2/r^2 + r''/r) Y - W^2 Y = 0,
s = r' r^((n-2)/2), w = s^2, and the chain Y = r Z, q = Z'/w, chi = s q.

Symbolic (sympy, n symbolic):
  S1  E2[r Z] = r (Z'' - (w'/w) Z' - W^2 Z) for any r(l): the Z term cancels. Control: Paper III's Y' coefficient with
      n -> n + 1 leaves a Z term.
  S2  with Z = w q'/W^2 the Z equation is (w q')' = W^2 w q, and with q = chi/s this is -chi'' + (s''/s) chi = -W^2 chi.
      Control: q = chi/w (the wrong weight) leaves a chi' term.
  S3  on the background, U = s''/s = (n-2)(n-4)/(4 r^2) - 3 n^2/(4 r^(2n)). Control: the same with 3 n^2 -> 3 n^2 + 1.
  S4  U is regular at the throat, U(0) = -(n + 4)(n - 1)/2, and even in l (a function of r alone).
  S5  chi_N = s H with H' = -1/s^2 solves chi'' = U chi (zero energy). Control: H' = -1/s does not.
Numerical (mpmath):
  N1  H(inf) = -B(1/2 - 1/m, -1/2)/(2m) for m = 3, 4, 6, 10, 25 (D = m + 3), from the integrated even zero-energy
      solution: H(inf) = chi_N(X)/s(X) - int_X^inf dl/s^2. Control: +B/(2m) must disagree.
  N2  B(1/2 - 1/m, -1/2) > 0 for m = 3 .. 200 (Gamma(1/2 - 1/m) > 0, Gamma(-1/2) < 0, Gamma(-1/m) < 0).
  N3  the even zero-energy solution (chi(0) = 1, chi'(0) = 0) has exactly one zero on (0, L) and the odd one (chi(0) = 0)
      none, for D = 4, 5, 6, 8, 12, 20, 50, by integrating chi'' = U(r(l)) chi.
  N4  cross-check with Paper III: the even bound state of -chi'' + U chi = E chi, found by shooting, has -E = W^2 of
      Paper III's table to 1e-6 at D = 4, 10, 20, 50. Control: the potential with 3n^2 -> 3n^2 + 1 moves it by more
      than 1e-3.
"""
import sys

import mpmath as mp
import sympy as sp

bad = []


def report(name, ok, msg, control=False):
    good = ok if not control else not ok
    print(f"{'PASS' if good else 'FAIL'}  {name}{' (control)' if control else ''}: {msg}")
    if not good:
        bad.append(name)


def symbolic():
    l, W, n = sp.symbols("l W n", positive=True)
    r = sp.Function("r")(l)
    Z = sp.Function("Z")(l)
    rp, rpp = sp.diff(r, l), sp.diff(r, l, 2)
    s = rp * r ** ((n - 2) / 2)
    w = s ** 2

    def E2(Y, shift=0):
        return (sp.diff(Y, l, 2) - (2 * rpp / rp + (n + shift) * rp / r) * sp.diff(Y, l)
                + (n * rp ** 2 / r ** 2 + rpp / r) * Y - W ** 2 * Y)

    target = r * (sp.diff(Z, l, 2) - sp.diff(w, l) / w * sp.diff(Z, l) - W ** 2 * Z)
    d = sp.simplify(sp.expand(E2(r * Z) - target))
    report("S1", d == 0, f"E2[rZ] - r(Z'' - (w'/w)Z' - W^2 Z) = {d}")
    d = sp.simplify(sp.expand(E2(r * Z, shift=1) - target))
    report("S1", d == 0, f"with n -> n + 1 in the Y' coefficient: residual {sp.factor(d)}", control=True)

    # S2: the q and chi equations
    chi = sp.Function("chi")(l)
    q = chi / s
    Zq = w * sp.diff(q, l) / W ** 2
    zeq = sp.diff(Zq, l, 2) - sp.diff(w, l) / w * sp.diff(Zq, l) - W ** 2 * Zq
    # Z = w q'/W^2 turns the Z equation into d/dl of [(w q')' - W^2 w q] / W^2 (times 1/w): check the bracket form
    qeq = sp.diff(w * sp.diff(q, l), l) - W ** 2 * w * q
    lhs = sp.simplify(sp.expand(zeq - sp.diff(qeq / w, l) * w / W ** 2))
    report("S2", lhs == 0, f"Z equation at Z = w q'/W^2 equals w (((w q')' - W^2 w q)/w)'/W^2: residual {lhs}")
    schr = sp.simplify(sp.expand(qeq - s * (sp.diff(chi, l, 2) - sp.diff(s, l, 2) / s * chi - W ** 2 * chi)))
    report("S2", schr == 0, f"(w q')' - W^2 w q = s (chi'' - (s''/s) chi - W^2 chi) with q = chi/s: residual {schr}")
    qbad = chi / w
    qeqb = sp.diff(w * sp.diff(qbad, l), l) - W ** 2 * w * qbad
    cb = sp.simplify(sp.expand(qeqb - s * (sp.diff(chi, l, 2) - sp.diff(s, l, 2) / s * chi - W ** 2 * chi) / s))
    report("S2", cb == 0, "q = chi/w instead of chi/s leaves a residual", control=True)

    # S3, S4 on the background: r'' = m r^(-2m-1), r''' = r' d/dr(r''), r'^2 = 1 - r^(-2m)
    R, P = sp.symbols("R P", positive=True)   # r and r'
    m = n - 1
    sub = {sp.Derivative(r, (l, 3)): sp.diff(m * R ** (-2 * m - 1), R) * P,
           sp.Derivative(r, (l, 2)): m * R ** (-2 * m - 1), sp.Derivative(r, l): P}
    U = (sp.diff(s, l, 2) / s)
    U = U.subs(sub).subs(r, R)
    U = sp.simplify(U.subs(P, sp.sqrt(1 - R ** (-2 * m))))
    Uform = (n - 2) * (n - 4) / (4 * R ** 2) - 3 * n ** 2 / (4 * R ** (2 * n))
    d = sp.simplify(U - Uform)
    report("S3", d == 0, f"s''/s - [(n-2)(n-4)/(4r^2) - 3n^2/(4r^(2n))] = {d}")
    d = sp.simplify(U - (Uform + 1 / (4 * R ** (2 * n))))
    report("S3", d == 0, "the formula with 3n^2 -> 3n^2 + 1 differs", control=True)
    U0 = sp.simplify(Uform.subs(R, 1) + (n + 4) * (n - 1) / 2)
    report("S4", U0 == 0, f"U(r = 1) + (n + 4)(n - 1)/2 = {U0}; U depends on l only through r(l), which is even")

    # S5: chi_N = s H, H' = -1/s^2 (as functions of l, any s)
    S = sp.Function("s")(l)
    Hf = sp.Function("H")(l)
    for pw, ctl in ((2, False), (1, True)):
        expr = sp.diff(S * Hf, l, 2) - sp.diff(S, l, 2) / S * S * Hf
        expr = expr.subs(sp.Derivative(Hf, (l, 2)), sp.diff(-1 / S ** pw, l)).subs(sp.Derivative(Hf, l), -1 / S ** pw)
        e = sp.simplify(expr)
        report("S5", e == 0, f"(sH)'' - (s''/s) sH with H' = -1/s^{pw}: {e}", control=ctl)


def numerical():
    def shoot(D, E, even, Ufun=None):
        """(r, r', chi, chi') along l from the throat: r'' = m r^(-2m-1), chi'' = (U(r) - E) chi (mpmath odefun)"""
        n = D - 2
        m = n - 1
        U = Ufun or (lambda r: (n - 2) * (n - 4) / (4 * r ** 2) - 3 * n ** 2 / (4 * r ** (2 * n)))
        y0 = [mp.mpf(1), mp.mpf(0), mp.mpf(1) if even else mp.mpf(0), mp.mpf(0) if even else mp.mpf(1)]
        return mp.odefun(lambda t, y: [y[1], m * y[0] ** (-2 * m - 1), y[3], (U(y[0]) - E) * y[2]], 0, y0)

    mp.mp.dps = 20
    # N1: H(inf) from the integrated even zero-energy solution, chi_N = chi/m with chi(0) = 1, chi'(0) = 0, as
    # chi_N(X)/s(X) - int_X^inf dl/s^2 (the tail in r, no cancellation), against -B(1/2 - 1/m, -1/2)/(2m)
    for m in (3, 4, 6, 10, 25):
        D = m + 3
        X = mp.mpf(30)
        y = shoot(D, 0, True)(X)
        r, rp = y[0], y[1]
        sX = rp * r ** (mp.mpf(m - 1) / 2)
        tail = mp.quad(lambda u: 1 / ((1 - u ** (-2 * m)) ** mp.mpf(1.5) * u ** (m - 1)), [r, 10 * r, mp.inf])
        Hinf = y[2] / m / sX - tail
        b = mp.beta(mp.mpf(1) / 2 - mp.mpf(1) / m, -mp.mpf(1) / 2) / (2 * m)
        report("N1", abs(Hinf + b) < 1e-8 * abs(b), f"m = {m}: H(inf) from the ODE {mp.nstr(Hinf, 12)}, "
                                                   f"-B/(2m) {mp.nstr(-b, 12)}")
        report("N1", abs(Hinf - b) < 1e-8 * abs(b), f"m = {m}: against +B/(2m)", control=True)
    ok = all(mp.beta(mp.mpf(1) / 2 - mp.mpf(1) / m, -mp.mpf(1) / 2) > 0 for m in range(3, 201))
    signs = (mp.gamma(mp.mpf(1) / 6) > 0, mp.gamma(-mp.mpf(1) / 2) < 0, mp.gamma(-mp.mpf(1) / 3) < 0)
    report("N2", ok and all(signs), f"B(1/2 - 1/m, -1/2) > 0 for m = 3..200; Gamma signs at m = 3: {signs}")

    mp.mp.dps = 15
    for D in (4, 5, 6, 8, 12, 20, 50):
        L = 40 if D < 8 else 15
        res = {}
        for even in (True, False):
            f = shoot(D, 0, even)
            vals = [f(mp.mpf(L) * k / 400)[2] for k in range(1, 401)]
            res[even] = sum(1 for a, b in zip(vals, vals[1:]) if a * b < 0)
        report("N3", res[True] == 1 and res[False] == 0,
               f"D = {D}: zeros on (0, {L}] of the even zero-energy solution {res[True]}, of the odd one {res[False]}")

    # N4: the even bound state of -chi'' + U chi = E chi against Paper III's W^2 (Torii-Shinkai's 1.39705243 at D 4,
    # Paper III's r-integration above); g(E) = chi' + k chi at l = L, k = sqrt(-E), vanishes at a bound state
    W2 = {4: 1.39705243, 10: 12.0442650, 20: 31.5751101276, 50: 91.3457759034}
    for D, w2 in W2.items():
        n = D - 2
        for ctl in (False, True):
            Uf = None if not ctl else (lambda r, n=n: (n - 2) * (n - 4) / (4 * r ** 2) - (3 * n ** 2 + 1) / (4 * r ** (2 * n)))
            L = 12 / mp.sqrt(w2)

            def g(E):
                y = shoot(D, E, True, Uf)(L)
                return y[3] + mp.sqrt(-E) * y[2]      # the growing solution's share at l = L (chi(0) = 1)

            try:
                E = mp.findroot(g, (-mp.mpf(w2) * (1 - 1e-4), -mp.mpf(w2) * (1 + 1e-4)), solver="anderson", tol=1e-24, verify=False)
                rel = abs(-E / w2 - 1)
            except Exception as exc:   # noqa: BLE001
                rel = mp.inf
            ok = rel < 1e-6 if not ctl else rel > 1e-3
            report("N4", ok, f"D = {D}: even bound state -E = {mp.nstr(-E, 12) if rel < mp.inf else 'none'} against "
                             f"W^2 = {w2} (relative {mp.nstr(rel, 2)})" + (" with 3n^2 -> 3n^2 + 1" if ctl else ""),
                   control=False)


if __name__ == "__main__":
    symbolic()
    numerical()
    print(f"{len(bad)} failed: {bad}" if bad else "all checks pass, all controls fire")
    sys.exit(1 if bad else 0)
