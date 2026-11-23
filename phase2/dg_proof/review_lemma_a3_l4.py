"""#142: L4 of RESULTS.md (the DG mode), from RESULTS.md alone.

Psi_DG = r (r+2)^k e^{-kr/2}/(cr-6), H_Z = -d_x^2 + V (V of 4.2 at M = -1, lambda = c/2), H-hat = -d_x^2 + V-hat,
A = d_x + g, A* = -d_x + g, g = psi_0,x/psi_0, psi_0 = chi/(cr-6), f = 1 + 2/r, d_x = f d_r.
Psi_DG = E rho with E = (r+2)^k e^{-kr/2}, e = E_r/E = k/(r+2) - k/2, rho = r/(cr-6); on E-stripped functions
d_x acts as f D with D = d_r + e.
Run: python review_lemma_a3_l4.py [--mutate=NAME] > review_lemma_a3_l4_output.txt
Mutations: pl_vhat, pl_eig, l4_alg, v_typo, l4_k, cn_n4to3, l4_sign.
"""
import sys
import time
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from review_lemma_a3_common import *

run = Run('review_lemma_a3_l4: L4, the DG mode and the intertwining')
T0 = time.time()
c, k = symbols('c k', positive=True)


def Dop(rho, e):
    return sp.diff(rho, r) + e * rho


def eig_residual(V, rho, e, Eig):
    """[-Psi_xx + V Psi - Eig Psi]/E for Psi = E rho, E_r = e E: -f D(f D rho) + (V - Eig) rho.
    Mutation pl_eig: the eigenvalue term is dropped."""
    out = -f_of_r * Dop(f_of_r * Dop(rho, e), e) + V * rho
    return out if mut('pl_eig') else out - Eig * rho


# ------------------------------------------------------------------ planted problems
# PLL.1: psi = r: g = f/r, V-hat = 2(r+2)(r+3)/r^4 and V = psi_xx/psi = -2(r+2)/r^4 (by hand)
Vh = cancel(together(Vhat_of_psi(r)))
run.planted('PLL.1', cancel(Vh - 2 * (r + 2) * (r + 3) / r ** 4) == 0 and cancel(dx(dx(r)) / r + 2 * (r + 2) / r ** 4) == 0,
            'V-hat(psi=r) = %s' % sp.factor(Vh))
# PLL.2: Psi_k = (r+2)^k e^{-kr/2} = 2^k e^{-kx/2}, V = 0: eigenvalue -(k/2)^2 exactly; +1 off is -1
ek = k / (r + 2) - k / 2
res0 = cancel(eig_residual(0, sp.Integer(1), ek, -k ** 2 / 4))
res1 = cancel(eig_residual(0, sp.Integer(1), ek, -k ** 2 / 4 + 1))
run.planted('PLL.2', res0 == 0 and res1 == -1, 'residual at E = -(k/2)^2: %s; at +1: %s (want 0, -1)' % (res0, res1))

# ------------------------------------------------------------------ L4.1 intertwining algebra, generic psi(x)
X = Symbol('X')
psiX = sp.Function('psi')(X)
uX = sp.Function('u')(X)
gX = sp.diff(psiX, X) / psiX
A_ = lambda w: sp.diff(w, X) + gX * w
As_ = lambda w: -sp.diff(w, X) + (-1 if mut('l4_alg') else 1) * gX * w   # mutation l4_alg: A* = -d - g
VZ = sp.diff(psiX, X, 2) / psiX
VH_ = gX ** 2 - sp.diff(gX, X)
HZ = lambda w: -sp.diff(w, X, 2) + VZ * w
HH = lambda w: -sp.diff(w, X, 2) + VH_ * w
r1 = sp.simplify(A_(As_(uX)) - HZ(uX))
r2 = sp.simplify(As_(A_(uX)) - HH(uX))
r3 = sp.simplify(HH(As_(uX)) - As_(HZ(uX)))
r4 = sp.simplify(HZ(A_(uX)) - A_(HH(uX)))
run.check('L4.1', r1 == 0 and r2 == 0 and r3 == 0 and r4 == 0,
          'A A* - H_Z: %s; A* A - H-hat: %s; H-hat A* - A* H_Z: %s; H_Z A - A H-hat: %s' % (r1, r2, r3, r4))

# ------------------------------------------------------------------ L4.2 the DG identity with k symbolic
rho = r / (c * r - 6)
Vz = V_results(c / 2)
R = cancel(together(eig_residual(Vz, rho, ek, -k ** 2 / 4) / rho))
num = sp.numer(R)
kstar = c * (c + 2) / (5 if mut('l4_k') else 6)   # mutation l4_k: k = c(c+2)/5
zero_at = cancel(num.subs(k, kstar)) == 0
factored = sp.factor(num)
divisible = sp.fraction(cancel(num / (6 * k - c * (c + 2))))[1] == 1
nonzero_off = cancel(num.subs(k, c * (c + 2) / 6 + 1)) != 0
run.check('L4.2', zero_at and divisible and nonzero_off,
          'numerator of [-Psi_xx + V Psi + (k/2)^2 Psi]/Psi vanishes at k = %s: %s; divisible by (6k - c(c+2)): %s; nonzero at k+1: %s' % (
              kstar, zero_at, divisible, nonzero_off))
run.note('numerator factored: %s' % factored)

# ------------------------------------------------------------------ L4.3 k = K/6 integer, multiple of 4, positive
ok = all(K_of(lv) % 24 == 0 and K_of(lv) // 6 > 0 and K_of(lv) // 6 == c_of(lv) * (c_of(lv) + 2) // 6 and K_of(lv) == c_of(lv) * (c_of(lv) + 2)
         for lv in range(2, 201))
run.check('L4.3', ok, 'K = c(c+2) = 6k, k integer divisible by 4, l = 2..200')

# ------------------------------------------------------------------ L4.4, L4.5 exact, l = 2..LMAX
LMAX = 20
for lv in range(2, LMAX + 1):
    cv = c_of(lv)
    psi0 = chi_poly(lv).as_expr() / (cv * r - 6)
    res = cancel(together(dx(dx(psi0)) - V_results(Q(cv, 2)) * psi0))
    run.check('L4.4[l=%d]' % lv, res == 0, 'psi_0,xx - V psi_0 = %s' % res)

for lv in range(2, LMAX + 1):
    cv = c_of(lv)
    kv = k_of(lv)
    psi0 = chi_poly(lv).as_expr() / (cv * r - 6)
    u = sp.diff(psi0, r) / psi0
    g_ = f_of_r * u
    ev = kv / (r + 2) - kv / 2
    rho_ = r / (cv * r - 6)
    # Psi-hat/E = A* Psi_DG / E = -f D rho + g rho, from the definition (no formula for the Wronskian used)
    rhohat = cancel(together(-f_of_r * Dop(rho_, ev) + g_ * rho_))
    Vh_ = cancel(together(Vhat_of_u(u)))
    sg = -1 if mut('l4_sign') else 1   # mutation l4_sign: +(k/2)^2 as eigenvalue
    res = cancel(together(eig_residual(Vh_, rhohat, ev, sg * (-kv ** 2 / 4))))
    run.check('L4.5[l=%d]' % lv, res == 0 and rhohat != 0, 'H-hat Psi-hat + (k/2)^2 Psi-hat over E = %s; rho-hat nonzero %s' % (res, rhohat != 0))
run.note('elapsed %.1f s' % (time.time() - T0))
run.finish()
