"""
Tasks 2-4 (symbolic part): EIH two-body problem with signed masses mu_A = s_A m_A.

Starting point (Landau-Lifshitz Sec. 106, eq. 106.17; Blanchet Living Rev. Relativ. 17 (2014) 2, eq. 214
at 1PN), harmonic gauge, G = c = 1:

  L = sum_A mu_A v_A^2/2 + sum_A mu_A v_A^4/8 + mu1 mu2 / r
      + mu1 mu2/(2r) [3(v1^2 + v2^2) - 7 v1.v2 - (n.v1)(n.v2)]  -  mu1 mu2 (mu1 + mu2)/(2 r^2)

(the constant -sum mu_A is dropped). The point-particle action is S = -sum mu_A int dtau_A + S_EH, i.e.
the masses enter only through mu_A, so the reduced (Fokker) action - and hence L - is the EIH one with
m_A -> mu_A (see RESULTS.md for the argument). This script:

  1. derives the Euler-Lagrange accelerations, order-reduces them to 1PN, and compares with
     Blanchet's eq. (203) harmonic-gauge acceleration (with m -> mu);
  2. builds E, P, J, G (centre-of-mass integral) and checks dE/dt, dP/dt, dJ/dt = O(c^-4),
     dG/dt - P = O(c^-4) symbolically;
(The collinear / dipole specialisation is in collinear.py.)
Bookkeeping: eps marks 1/c^2.
"""
import sympy as sp

eps = sp.Symbol('epsilon')                     # 1/c^2
mu1, mu2 = sp.symbols('mu1 mu2', real=True, nonzero=True)
x1 = sp.Matrix(sp.symbols('x1:4', real=True)); x2 = sp.Matrix(sp.symbols('y1:4', real=True))
v1 = sp.Matrix(sp.symbols('u1:4', real=True)); v2 = sp.Matrix(sp.symbols('w1:4', real=True))
a1s = sp.Matrix(sp.symbols('a1:4')); a2s = sp.Matrix(sp.symbols('b1:4'))

dx = x1 - x2
r = sp.sqrt(dx.dot(dx))
n = dx / r

L = (mu1 * v1.dot(v1) / 2 + mu2 * v2.dot(v2) / 2 + mu1 * mu2 / r
     + eps * (mu1 * v1.dot(v1)**2 / 8 + mu2 * v2.dot(v2)**2 / 8
              + mu1 * mu2 / (2 * r) * (3 * (v1.dot(v1) + v2.dot(v2)) - 7 * v1.dot(v2) - n.dot(v1) * n.dot(v2))
              - mu1 * mu2 * (mu1 + mu2) / (2 * r**2)))

def tdiff(expr):
    """Total time derivative using symbols for positions, velocities and accelerations."""
    out = 0
    for X, V, A in ((x1, v1, a1s), (x2, v2, a2s)):
        for i in range(3):
            out += sp.diff(expr, X[i]) * V[i] + sp.diff(expr, V[i]) * A[i]
    return out

# --- 1. Euler-Lagrange equations ------------------------------------------------------------
def EL(Xs, Vs):
    return sp.Matrix([tdiff(sp.diff(L, Vs[i])) - sp.diff(L, Xs[i]) for i in range(3)])

E1, E2 = EL(x1, v1), EL(x2, v2)
# EL_A = mu_A a_A + R_A, where R_A depends on accelerations only at O(eps). Order reduction:
# a_A = -R_A[a -> a_Newton]/mu_A + O(eps^2)  (standard; equivalent to solving and truncating).
aN1 = -mu2 * n / r**2
aN2 = mu1 * n / r**2
subN = {a1s[i]: aN1[i] for i in range(3)}; subN.update({a2s[i]: aN2[i] for i in range(3)})

def reduce(EA, aA, muA):
    out = []
    for i in range(3):
        R = sp.expand(EA[i] - muA * aA[i])
        assert sp.expand(R.coeff(eps, 0)).free_symbols.isdisjoint(set(a1s) | set(a2s))
        out.append(sp.expand(-(R.coeff(eps, 0) + eps * R.coeff(eps, 1).subs(subN)) / muA))
    return sp.Matrix(out)

acc1 = reduce(E1, a1s, mu1)
acc2 = reduce(E2, a2s, mu2)

# Blanchet (2014) eq. (203), harmonic gauge, 1PN part, m -> mu, n12 = (y1 - y2)/r12, v12 = v1 - v2
def blanchet_a1(m1, m2, n12, va, vb):
    return (-m2 * n12 / r**2
            + eps * ((5 * m1 * m2 / r**3 + 4 * m2**2 / r**3
                      + m2 / r**2 * (sp.Rational(3, 2) * n12.dot(vb)**2 - va.dot(va) + 4 * va.dot(vb) - 2 * vb.dot(vb))) * n12
                     + m2 / r**2 * (4 * n12.dot(va) - 3 * n12.dot(vb)) * (va - vb)))

B1 = blanchet_a1(mu1, mu2, n, v1, v2)
B2 = blanchet_a1(mu2, mu1, -n, v2, v1)

import random
import mpmath as mp
mp.mp.dps = 50
random.seed(1)
ALL = list(x1) + list(x2) + list(v1) + list(v2) + [mu1, mu2]
def rand_point():
    pt = [mp.mpf(random.randint(-9, 9)) / random.randint(1, 7) for _ in range(12)]
    pt += [mp.mpf(random.randint(1, 9)) / 5, -mp.mpf(random.randint(1, 9)) / 5]   # mu1 > 0 > mu2
    return pt

def eps_coeffs(f, pt, deg=4):
    """f(pt, eps) is a polynomial in eps of degree <= deg: recover its coefficients exactly
    (to working precision) from deg+1 evaluations."""
    es = [mp.mpf(k) for k in range(deg + 1)]
    vals = [f(*pt, e) for e in es]
    A = mp.matrix([[e**p for p in range(deg + 1)] for e in es])
    return mp.lu_solve(A, mp.matrix(vals))

def check_zero(expr, npts=5):
    f = sp.lambdify(ALL + [eps], expr, 'mpmath')
    worst = mp.mpf(0)
    for _ in range(npts):
        c = eps_coeffs(f, rand_point())
        worst = max(worst, abs(c[0]), abs(c[1]))
    return worst

print("1. max |EL(EIH) - Blanchet eq.(203)| (O(1) and O(eps) coefficients, random signed-mass points):",
      mp.nstr(max(check_zero(acc1[i] - B1[i]) for i in range(3)), 3),
      mp.nstr(max(check_zero(acc2[i] - B2[i]) for i in range(3)), 3))

# --- 2. conserved quantities -----------------------------------------------------------------
p1 = sp.Matrix([sp.diff(L, v1[i]) for i in range(3)])
p2 = sp.Matrix([sp.diff(L, v2[i]) for i in range(3)])
P = p1 + p2
E = v1.dot(p1) + v2.dot(p2) - L
J = x1.cross(p1) + x2.cross(p2)
G = (mu1 * x1 * (1 + eps * (v1.dot(v1) / 2 - mu2 / (2 * r)))
     + mu2 * x2 * (1 + eps * (v2.dot(v2) / 2 - mu1 / (2 * r))))
on_shell = {a1s[i]: acc1[i] for i in range(3)}; on_shell.update({a2s[i]: acc2[i] for i in range(3)})

def rate(Q, minus=0):
    """max over random points of the O(1) and O(eps) coefficients of dQ/dt - minus (on shell)."""
    return check_zero(tdiff(Q).subs(on_shell) - minus)

for name, Q, minus in (("E", [E], [0]), ("P", list(P), [0, 0, 0]), ("J", list(J), [0, 0, 0]),
                       ("G - P t", list(G), list(P))):
    worst = max(rate(Qi, mi) for Qi, mi in zip(Q, minus))
    print(f"2. d({name})/dt: max |O(1) and O(eps) coefficient| at random points = {mp.nstr(worst, 3)}")
# control: a wrong centre-of-mass integral (no 1PN correction) must fail
print("   control, G without 1PN terms: ", mp.nstr(max(rate(mu1 * x1[i] + mu2 * x2[i], P[i]) for i in range(3)), 3))

# Sections 3-4 (collinear pair, dipole, rigidity condition) are in collinear.py.
