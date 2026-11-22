"""Step 3 of the independent review of #72 (issue #101): T1, T2, T3, the two identities, and cofactor certificates.

Everything here is a polynomial identity in (v, F, a, b, u, c, s, lam), reduced modulo c^2 + s^2 = 1, from my own frame
Ricci tensor (review_he_surface_lib.ricci_jet). Index convention: Sigma(lam) = R - lam eta on (n, d_x, d_rho) = indices
(0, 1, 2); m01,02 is the minor with rows (0, 1) and columns (0, 2); six distinct minors of the symmetric 3x3 Sigma.

  rotated_quantities  the entries of he_surface.md section 1 in the frame (n, r^, theta^), built from my Ricci by my own
                      rotation, for the comparison only
  identity_necessity  j_theta v c a u (1 + v^2 F^2) + m_e + c v F m_d = 0 for every lam
  D_on_T1             D(b = -a u) = v^6 a^6 u^6 (v^2 F^2 s^4 - 4 c^2) P2^2, and D = disc(chi_2) chi_2(lam*)^2
  cert_T1             the unrotated certificate  K1 (b + a u) = sum_k q_k m_k  for every lam, K1 = v^2 a u c s (1 + v^2 F^2) / 2,
                      and the a = 0 identity m01,01 = -Sigma_00^2 - Sigma_01^2
  cert_T2             the lam certificate, the P2 certificate, and the converse  m_k = (b + a u) A_k + (lam - lam*) B_k + Pu C_k
  T3_multiplicity     on T1: chi = (x - lam*) [ (x - lam*)(x - lam* - Theta) + v^2 a^2 Pu ], the eigenvector theta^, and the
                      offset Theta = v^2 a (2 a + (1 + c^2) F u)

Notation. u = 1/r, Pu = u^2 P2 = v^2 (a + F u)(a + c^2 F u) + c^2 u^2, lam* = -v^2 a s^2 (a + 2 F u)/2.
Run: python3 notes/warp/review_he_surface_typeI.py [--no-controls | --mutate=NAME | --derive]
`--derive` re-finds the T1 cofactors by linear algebra over Q (about 40 s). Written without reading he_surface.py,
he_type.py, review_he_derive.py or phantom_source.py.
"""
import sys
import time

import sympy as sp

import review_he_surface_lib as L
from review_he_surface_lib import F, a, b, c, lam, s, u, v, x

MUT = ""
T1 = {b: -a * u}


def set_mut(name):
    global MUT
    MUT = name
    L.MUT = name


def lamstar():
    return sp.expand(-sp.Rational(1, 2) * v**2 * a * s**2 * (a + 2 * F * u))


def Pu():
    return sp.expand(v**2 * (a + F * u) * (a + c**2 * F * u) + c**2 * u**2)


def Theta():
    cc = (1 + c) if MUT == "tgt-T3" else (1 + c**2)
    return sp.expand(v**2 * a * (2 * a + cc * F * u))


def setup():
    R3, _ = L.frame_block()
    S = L.sigma(R3)
    mins = [L.red(m) for m in L.minors(S)]
    return R3, S, mins


# ----------------------------------------------------------------------------------------------------------------------
def rotated(R3, lm=lam):
    """Sigma' = Q^T Sigma Q with Q = diag(1, [[c, -s], [s, c]]): the frame (n, r^, theta^), r^ = (c, s), theta^ = (-s, c)."""
    sg = -1 if MUT == "rot-sign" else 1
    Q = sp.Matrix([[1, 0, 0], [0, c, -sg * s], [0, sg * s, c]])
    S = L.sigma(R3, lm)
    return sp.Matrix(3, 3, lambda i, j: L.red((Q.T * S * Q)[i, j]))


def rotated_quantities(h):
    R3, S, mins = setup()
    Sp = rotated(R3)
    jr, jth, tau = Sp[0, 1], Sp[0, 2], Sp[1, 2]
    exp = {
        "j_r": (jr, v * c * a * u),
        "j_theta": (jth, -(v * s / 2) * (b + a * u)),
        "R'_r,theta": (tau, v**2 * F * c * s / 2 * (b + a * u)),
        "tau + c v F j_theta": (tau + c * v * F * jth, 0),
        "R'_00": (Sp[0, 0] - lam, v**2 * (F * (a * u * s**2 + b * c**2) - a**2 * (1 + c**2) / 2)),
        "R'_rr": (Sp[1, 1] + lam, v**2 * (a**2 * (1 + c**2) / 2 - F * b * c**2)),
        "R'_theta,theta": (Sp[2, 2] + lam, -sp.Rational(1, 2) * v**2 * a * s**2 * (a + 2 * F * u)),
    }
    bad = [k for k, (got, want) in exp.items() if L.red(got - want) != 0]
    h.check("rotated_quantities", not bad,
            "j_r, j_theta, R'_r,theta = -c v F j_theta, R'_00, R'_rr, R'_theta,theta agree with he_surface.md section 1 (no restriction to T1); "
            "mismatches: %s" % (bad or "none"))


def identity_necessity(h):
    R3, S, mins = setup()
    Sp = rotated(R3)
    A, Bp, Cp = Sp[0, 0], Sp[1, 1], Sp[2, 2]
    jr, jth, tau = Sp[0, 1], Sp[0, 2], Sp[1, 2]
    md = jr * tau - jth * Bp
    me = A * tau - jth * jr
    fac = (1 + v * F**2) if MUT == "tgt-ident" else (1 + v**2 * F**2)
    ident = L.red(jth * v * c * a * u * fac + me + c * v * F * md)
    # the intermediate factorizations quoted in section 2
    f1 = L.red(md + jth * (Bp + c * v * F * jr))
    f2 = L.red(me + jth * (c * v * F * A + jr))
    f3 = L.red(A + Bp - v**2 * F * a * u * s**2)
    # the a = 0 step of section 2: m_a = A B' - j_r^2 = -A^2, m_b = A C' - j_theta^2, and the displayed identity
    # j_theta^4 + m_a C'^2 + 2 A C' m_b - m_b^2 = 0 holds there for every lam
    ma, mb = A * Bp - jr**2, A * Cp - jth**2
    a0 = [L.red(e.subs(a, 0)) for e in (ma + A**2, jr, jth**4 + ma * Cp**2 + 2 * A * Cp * mb - mb**2)]
    h.check("identity_necessity", ident == 0 and f1 == 0 and f2 == 0 and f3 == 0 and all(e == 0 for e in a0),
            "j_theta v c a u (1 + v^2 F^2) + m_e + c v F m_d reduces to %s for symbolic lam; m_d = -j_theta (B' + c v F j_r): %s; "
            "m_e = -j_theta (c v F A + j_r): %s; A + B' = v^2 F a u s^2: %s; at a = 0: m_a + A^2, j_r and "
            "j_theta^4 + m_a C'^2 + 2 A C' m_b - m_b^2 reduce to %s" % (ident, f1 == 0, f2 == 0, f3 == 0, a0))


def D_on_T1(h):
    R3, S, mins = setup()
    Sx = L.sigma(R3, x)
    chi = L.red(Sx.det())
    chiT = L.red(chi.subs(T1))
    D = sp.expand(sp.discriminant(chiT, x))
    D = L.red(D)
    r = 1 / u
    P2 = v**2 * (r * a + F) * (r * a + c**2 * F) + c**2
    k4 = 3 if MUT == "tgt-DT1" else 4
    closed = v**6 * a**6 * u**6 * (v**2 * F**2 * s**4 - k4 * c**2) * P2**2
    d1 = L.red(sp.expand(D - closed))
    th = Theta()
    disc2 = L.red(sp.expand(th**2 - 4 * v**2 * a**2 * Pu() - v**2 * a**2 * u**2 * (v**2 * F**2 * s**4 - k4 * c**2)))
    d3 = L.red(sp.expand(D - (th**2 - 4 * v**2 * a**2 * Pu()) * (v**2 * a**2 * Pu()) ** 2))
    h.check("D_on_T1", d1 == 0 and disc2 == 0 and d3 == 0,
            "D(b = -a u) - v^6 a^6 u^6 (v^2 F^2 s^4 - %dc^2) P2^2 = %s; Theta^2 - 4 v^2 a^2 Pu - v^2 a^2 u^2 (v^2 F^2 s^4 - %dc^2) = %s; "
            "D - disc(chi_2) chi_2(lam*)^2 = %s" % (k4, d1, k4, disc2, d3))


# ----------------------------------------------------------------------------------------------------------------------
def k1():
    return sp.Rational(1, 2) * v**2 * a * u * c * s * (1 + v**2 * F**2)


def cof_T1():
    """cofactors of (m01,01, m01,02, m01,12, m02,02, m02,12, m12,12) in K1 (b + a u) = sum q_k m_k"""
    q = [-c * s, c**2 - s**2, v * F * c**2, c * s, v * F * c * s, 0]
    if MUT == "cof-perturb":
        q[0] = q[0] + 1
    if MUT == "cof-drop":
        q[1] = c**2
    return q


def cert_T1(h):
    R3, S, mins = setup()
    q = cof_T1()
    lhs = L.red(sum(qk * mk for qk, mk in zip(q, mins)))
    d = L.red(lhs - k1() * (b + a * u))
    # a = 0: Sigma_11 = -Sigma_00 there, so m01,01 = -Sigma_00^2 - Sigma_01^2
    sos = L.red((mins[0] + (R3[0, 0] + lam) ** 2 + R3[0, 1] ** 2).subs(a, 0))
    nz = sp.factor(k1())
    h.check("cert_T1", d == 0 and sos == 0,
            "-cs m01,01 + (c^2 - s^2) m01,02 + vF c^2 m01,12 + cs m02,02 + vF cs m02,12 - K1 (b + a u) reduces to %s for symbolic lam; "
            "at a = 0: m01,01 + Sigma_00^2 + Sigma_01^2 reduces to %s; K1 = %s" % (d, sos, nz))


def cert_T2(h):
    R3, S, mins = setup()
    ls, pu = lamstar(), Pu()
    ok = []
    # (lam) certificate: -s m01,12 + c m02,12 = v c a u (lam* - lam) + (F c s^2 v^3 / 4) (b + a u)^2
    g2 = sp.Rational(1, 4) * F * c * s**2 * v**3 * (b + a * u)
    if MUT == "cof2-drop":
        g2 = sp.Rational(1, 4) * c * s**2 * v**3 * (b + a * u)
    d2 = L.red(-s * mins[2] + c * mins[4] - (v * c * a * u * (ls - lam) + g2 * (b + a * u)))
    ok.append(("lam certificate", d2 == 0))
    # (P2) certificate: m01,01 + m02,02 + v^2 a^2 Pu in the ideal of (b + a u) and (lam - lam*), cofactors by exact division
    E = L.red(mins[0] + mins[3] + v**2 * a**2 * pu)
    E_T1 = sp.expand(E.subs(T1))
    A3 = sp.cancel((E - E_T1) / (b + a * u))
    B3 = sp.cancel((E_T1 - E_T1.subs(lam, ls)) / (lam - ls))
    d3 = L.red(E_T1.subs(lam, ls))
    gens = (v, F, a, b, u, c, s, lam)
    ok.append(("P2 certificate", d3 == 0 and A3.is_polynomial(*gens) and B3.is_polynomial(*gens)))
    # converse: m_k = (b + a u) A_k + (lam - lam*) B_k + Pu C_k with C_k typed
    Ck = [-v**2 * a**2 * c**2, -v**2 * a**2 * c * s, 0, -v**2 * a**2 * s**2, 0, 0]
    if MUT == "cof2-perturb":
        Ck[2] = Ck[2] + 1
    res = []
    for k, mk in enumerate(mins):
        mT = sp.expand(mk.subs(T1))
        Ak = sp.cancel((mk - mT) / (b + a * u))
        Bk = sp.cancel((mT - mT.subs(lam, ls)) / (lam - ls))
        poly_ok = Ak.is_polynomial(*gens) and Bk.is_polynomial(*gens)
        d = L.red(mk - (b + a * u) * Ak - (lam - ls) * Bk - pu * Ck[k])
        res.append(d == 0 and poly_ok)
    ok.append(("converse minors in <b + a u, lam - lam*, Pu>", all(res)))
    h.check("cert_T2", all(o for _, o in ok), "; ".join("%s %s" % (n, "yes" if o else "NO") for n, o in ok)
            + "; per-minor converse %s" % res)


def T3_multiplicity(h):
    R3, S, mins = setup()
    ls, pu, th = lamstar(), Pu(), Theta()
    Sx = L.sigma(R3, x)
    chiT = L.red(L.red(Sx.det()).subs(T1))
    chi2 = (x - ls) * (x - ls - th) + v**2 * a**2 * pu
    d1 = L.red(sp.expand(chiT - (x - ls) * chi2))
    # eigenvector theta^ = (0, -s, c) of M with eigenvalue lam* on T1
    ev = sp.Matrix([0, -s, c])
    Sl = L.sigma(R3, ls)
    d2 = [L.red(e) for e in (Sl * ev).subs(T1)]
    # chi_2(lam*) = v^2 a^2 Pu and chi_2'(lam*) = -Theta
    d3 = L.red(sp.expand(chi2.subs(x, ls) - v**2 * a**2 * pu))
    d4 = L.red(sp.expand(sp.diff(chi2, x).subs(x, ls) + th))
    # Theta in the rho a form of he_surface.md: B' - A = v^2 a u (2 r a + (1 + c^2) F)
    r = 1 / u
    d5 = L.red(sp.expand(th - v**2 * a * u * (2 * r * a + (1 + c**2) * F)))
    # section 2's values at lam = lam* in the rotated frame: A = -v^2 a (a + F u c^2), B' = v^2 a (a + F u), A B' - j_r^2 = -v^2 a^2 Pu
    Sp = rotated(R3, ls)
    A_, B_, jr_ = Sp[0, 0], Sp[1, 1], Sp[0, 1]
    d6 = [L.red(A_.subs(T1) + v**2 * a * (a + F * u * c**2)), L.red(B_.subs(T1) - v**2 * a * (a + F * u)),
          L.red((A_ * B_ - jr_**2).subs(T1) + v**2 * a**2 * pu), L.red((B_ - A_).subs(T1) - th)]
    h.check("T3_multiplicity", d1 == 0 and all(e == 0 for e in d2) and d3 == 0 and d4 == 0 and d5 == 0 and all(e == 0 for e in d6),
            "chi|T1 - (x - lam*)[(x - lam*)(x - lam* - Theta) + v^2 a^2 Pu] = %s; Sigma(lam*) theta^ on T1 = %s; chi_2(lam*) - v^2 a^2 Pu = %s; "
            "chi_2'(lam*) + Theta = %s; Theta - v^2 a u (2 r a + (1 + c^2) F) = %s; at lam*: A + v^2 a (a + F u c^2), B' - v^2 a (a + F u), "
            "A B' - j_r^2 + v^2 a^2 Pu, B' - A - Theta reduce to %s" % (d1, d2, d3, d4, d5, d6))


# ----------------------------------------------------------------------------------------------------------------------
def derive():
    """Find the lam-free combinations of the six minors that vanish on T1 by linear algebra over Q at random exact points."""
    import itertools
    import random
    from fractions import Fraction as Fr
    from sympy import QQ
    from sympy.polys.matrices import DomainMatrix
    t0 = time.time()
    R3, S, mins = setup()
    co = [[sp.Poly(m, lam).coeff_monomial(lam**d) for d in range(3)] for m in mins]
    terms = [[L.poly_terms(x_) for x_ in row] for row in co]
    exps = list(itertools.product(range(2), range(2), range(3), range(2)))  # v, F, c, s
    n = len(exps)
    random.seed(1)
    rows = []
    rr = lambda: Fr(random.randint(-9, 9), random.randint(1, 5))
    for _ in range(160):
        tt = Fr(random.randint(1, 12), random.randint(1, 12))
        cc, ss = (1 - tt * tt) / (1 + tt * tt), 2 * tt / (1 + tt * tt)
        vv, FF, aa, bb, uu = rr(), rr(), rr(), rr(), rr()
        for dgr in (2, 1, 0):
            b_ = bb if dgr > 0 else -aa * uu
            row = []
            for k in range(6):
                val = L.peval(terms[k][dgr], (vv, FF, aa, b_, uu, cc, ss))
                for e in exps:
                    row.append(val * vv**e[0] * FF**e[1] * cc**e[2] * ss**e[3])
            rows.append(row)
    dm = DomainMatrix([[QQ(r_.numerator, r_.denominator) for r_ in row] for row in rows], (len(rows), 6 * n), QQ)
    ns = dm.nullspace().to_Matrix()
    print("derive: %d unknowns, %d equations, nullspace dimension %d, %.0f s" % (6 * n, len(rows), ns.shape[0], time.time() - t0))
    mons = [v**e[0] * F**e[1] * c**e[2] * s**e[3] for e in exps]
    typed = sp.Matrix([[sp.Poly(L.red(qk), v, F, c, s).coeff_monomial(mo) for mo in mons] for qk in cof_T1()]).reshape(1, 6 * n)
    inspan = sp.Matrix.vstack(ns, typed).rank() == ns.rank()
    for r_ in range(ns.shape[0]):
        qs = [sum(ns[r_, k * n + i] * mons[i] for i in range(n)) for k in range(6)]
        tot = L.red(sum(qk * mk for qk, mk in zip(qs, mins)))
        print("  solution %d: lam-free %s, kappa = %s" % (r_, sp.Poly(tot, lam).degree() <= 0, sp.factor(sp.cancel(tot / (b + a * u)))))
    print("derive: the typed cofactors lie in the span of the solutions: %s" % inspan)
    return inspan


CHECKS = {
    "rotated_quantities": rotated_quantities,
    "identity_necessity": identity_necessity,
    "D_on_T1": D_on_T1,
    "cert_T1": cert_T1,
    "cert_T2": cert_T2,
    "T3_multiplicity": T3_multiplicity,
}
CONTROLS = {
    "rot-sign": (["rotated_quantities"], "the rotation is by -theta, r^ = (c, -s)"),
    "tgt-ident": (["identity_necessity"], "1 + v^2 F^2 -> 1 + v F^2 in the necessity identity"),
    "tgt-DT1": (["D_on_T1"], "4 c^2 -> 3 c^2 in the factor of D on T1"),
    "cof-perturb": (["cert_T1"], "one cofactor of the T1 certificate gets +1"),
    "cof-drop": (["cert_T1"], "one term of one cofactor of the T1 certificate is dropped"),
    "cof2-perturb": (["cert_T2"], "one cofactor of the T2 certificate gets +1"),
    "cof2-drop": (["cert_T2"], "one term of one cofactor of the T2 certificate is dropped"),
    "tgt-T3": (["T3_multiplicity"], "1 + c^2 -> 1 + c in the T3 offset"),
}

if __name__ == "__main__":
    print("review_he_surface_typeI.py: sympy %s" % sp.__version__, flush=True)
    if "--derive" in sys.argv:
        raise SystemExit(0 if derive() else 1)
    H = L.Harness(CHECKS, CONTROLS, "typeI")
    H.main(sys.argv[1:], set_mut)
