"""#142: L2 and L3 of RESULTS.md (V-hat and the Frobenius exponents at r = 0, V-hat at infinity), from RESULTS.md alone.

Run: python review_lemma_a3_ends.py [--mutate=NAME] > review_lemma_a3_ends_output.txt
Planted PLE.1 (indicial routine), PLE.2 (V-hat routine), PLE.3 (order/leading-coefficient routine) run first.
Mutations: pl_indicial, pl_vhat, pl_order, e_const12, e_linear, e_exp, e_den, e_nn1, e_lead_inf, e_phi, cn_n4to3.
"""
import sys
import time
sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from review_lemma_a3_common import *

run = Run('review_lemma_a3_ends: L2, L3 (ends of V-hat, Frobenius exponents)')
T0 = time.time()
LMAX = 30
E = Symbol('E')


# ------------------------------------------------------------------ planted problems
P2_, P1_, P0_ = wronskian_ode(r ** 3 * (1 + r), (1 + 2 * r + 3 * r ** 2) / r)
got = indicial(P2_, P1_, P0_, 0)
run.planted('PLE.1a', got == {3, -1}, 'Wronskian-built ODE: got %s want {3, -1}' % got)
got = indicial(r ** 2, r, r ** 2 - Q(1, 4), 0)
run.planted('PLE.1b', got == {Q(1, 2), -Q(1, 2)}, 'Bessel nu=1/2: got %s' % got)
got = indicial(r ** 2, -2 * r, sp.Integer(2), 0)
run.planted('PLE.1c', got == {1, 2}, 'Euler: got %s' % got)
got = indicial((r - 1) ** 2, (r - 1), sp.Integer(-4), 1)
run.planted('PLE.1d', got == {2, -2}, 'at r=1: got %s' % got)
P2_, P1_, P0_ = wronskian_ode((r - 1) ** 2 * (1 + r), 1 / (r - 1))
got = indicial(P2_, P1_, P0_, 1)
run.planted('PLE.1e', got == {2, -1}, 'Wronskian-built ODE at r=1: got %s want {2, -1}' % got)

for n in (1, 2, 3, 5):
    V = cancel(together(Vhat_of_psi(r ** n)))
    p0_, c0_ = lead(V, 0)
    pi_, ci_ = lead(V, 'inf')
    run.planted('PLE.2[n=%d]' % n, (p0_, c0_) == (-4, 4 * n * (n + 2)) and (pi_, ci_) == (-2, n * (n + 1)),
                'n=%d: at 0 %s, at inf %s; want (-4, %d), (-2, %d)' % (n, (p0_, c0_), (pi_, ci_), 4 * n * (n + 2), n * (n + 1)))
V = cancel(together(Vhat_of_psi(r ** 3 * (1 + 3 / r))))
run.planted('PLE.2[n=3, shifted]', lead(V, 'inf') == (-2, 12), 'got %s' % (lead(V, 'inf'),))

run.planted('PLE.3a', lead((r ** 2 + 1) / (r ** 4 * (1 + r)), 0) == (-4, 1), 'got %s' % (lead((r ** 2 + 1) / (r ** 4 * (1 + r)), 0),))
run.planted('PLE.3b', lead((r ** 2 + 1) / r ** 4, 'inf') == (-2, 1), 'got %s' % (lead((r ** 2 + 1) / r ** 4, 'inf'),))
run.planted('PLE.3c', lead(r ** 3 * (2 - r), 0) == (3, 2), 'got %s' % (lead(r ** 3 * (2 - r), 0),))

# ------------------------------------------------------------------ E2.1 symbolic: r^4 V-hat = 12 + 2(c+5) r + O(r^2)
c, h2, h3 = symbols('c h2 h3')
psi_gen = r * (Q(3, 2) + h2 * r ** 2 + h3 * r ** 3) / (c * r - 6)   # chi = r h(r), h = 3/2 + 0 r + h2 r^2 + h3 r^3
Vg = cancel(together(r ** 4 * Vhat_of_psi(psi_gen)))
ng, dg = sp.fraction(Vg)
tg = taylor(Poly(ng, r), Poly(dg, r), 3)
want0 = 14 if mut('e_const12') else 12   # mutation e_const12: leading coefficient 14
want1 = 2 * (l ** 2 + l + (2 if mut('e_linear') else 3))   # mutation e_linear: l^2+l+2 instead of l^2+l+3
ok0 = cancel(tg[0] - want0) == 0
ok1 = cancel(tg[1].subs(c, c_of(l)) - want1) == 0
run.check('E2.1', ok0 and ok1, 'r^0 coefficient %s (want %s); r^1 coefficient %s (want %s)' % (tg[0], want0, tg[1], want1))
v2_sym = cancel(tg[2].subs({c: c_of(l), h2: -K_of(l) / 32}))
run.note('r^2 coefficient of r^4 V-hat, with c3 = -K/32: %s  (generic: %s)' % (sp.factor(v2_sym), sp.factor(tg[2])))

# ------------------------------------------------------------------ exact V-hat for l = 2..LMAX
VH = {}
PSI = {}
for lv in range(2, LMAX + 1):
    cv = c_of(lv)
    PSI[lv] = chi_poly(lv).as_expr() / (cv * r - 6)
    VH[lv] = cancel(together(Vhat_of_psi(PSI[lv])))
run.note('V-hat built for l = 2..%d (%.1f s)' % (LMAX, time.time() - T0))

for lv in range(2, LMAX + 1):
    num, den = sp.fraction(VH[lv])
    pn, pd = Poly(num, r), Poly(den, r)
    low = min(m[0] for m in pd.monoms())
    pole = low - min(m[0] for m in pn.monoms())
    dlow = Poly(sp.expand(pd.as_expr() / r ** low), r)
    tc = taylor(pn, dlow, 3)
    scale = tc[0] / 12
    # r^4 V-hat = (num / den_low) with den = r^low den_low; low = 4 expected
    ok = pole == 4 and low == 4 and tc[0] == 12 and tc[1] == 2 * (lv ** 2 + lv + 3) and tc[2] == v2_sym.subs(l, lv)
    run.check('E2.2[l=%d]' % lv, ok, 'pole order %s; r^4 V-hat series %s; want [12, %d, %s]' %
              (pole, tc, 2 * (lv ** 2 + lv + 3), v2_sym.subs(l, lv)))

for lv in range(2, LMAX + 1):
    num, den = sp.fraction(VH[lv])
    pd = Poly(den, r)
    pw = 3 if mut('e_den') else 4   # mutation e_den: claim r^3 instead of r^4
    target = Poly(r ** pw * chi_over_r(lv).as_expr() ** 2, r, domain='QQ')
    ok = pd.monic() == target.monic()
    run.check('E2.3[l=%d]' % lv, ok, 'den V-hat is not const * r^%d (chi/r)^2' % pw)

# ------------------------------------------------------------------ E2.4 Frobenius exponents in r at 0
v1_, v2_ = symbols('v1 v2')
n0 = 11 if mut('e_exp') else 12   # mutation e_exp: leading coefficient 11 in the symbolic ODE
Vsym = (n0 + v1_ * r + v2_ * r ** 2) / r ** 4
exps = indicial(f_of_r ** 2, f_of_r * sp.diff(f_of_r, r), E - Vsym, 0)
run.check('E2.4', exps == {3, -1}, 'symbolic V-hat = N/r^4, N(0) = %s, E symbolic: exponents %s' % (n0, exps))
for lv in range(2, LMAX + 1):
    exps = indicial(f_of_r ** 2, f_of_r * sp.diff(f_of_r, r), E - VH[lv], 0)
    run.check('E2.4[l=%d]' % lv, exps == {3, -1}, 'exponents %s' % exps)

# ------------------------------------------------------------------ E2.5 x(r) and L^2(dx)
xs = sp.series(r - 2 * sp.log(1 + r / 2), r, 0, 6).removeO()
ok_x = Poly(xs, r).all_coeffs()[::-1][:4] == [0, 0, Q(1, 4), -Q(1, 12)]
ok_dx = cancel(sp.diff(r - 2 * sp.log(1 + r / 2), r) - 1 / f_of_r) == 0
i_bad = sp.integrate(1 / (r * (r + 2)), (r, 0, 1))        # int (r^-1)^2 dx = int r^-2 dr/f
i_good = sp.integrate(r ** 7 / (r + 2), (r, 0, 1))         # int (r^3)^2 dx
run.check('E2.5', ok_x and ok_dx and i_bad == sp.oo and i_good.is_finite,
          'x series %s; dx/dr = 1/f %s; int r^-2 dx = %s; int r^6 dx = %s' % (xs, ok_dx, i_bad, i_good))

# ------------------------------------------------------------------ E2.7 phi_0 = (cr-6)/chi
phi_pow0 = 3 if mut('e_phi') else 4   # mutation e_phi: claim phi_0 ~ -3/r
for lv in range(2, 51):
    cv = c_of(lv)
    chi = chi_poly(lv)
    phi = (cv * r - 6) / chi.as_expr()
    p0_, c0_ = lead(phi, 0)
    pi_, ci_ = lead(phi, 'inf')
    cl2 = cn(lv + 2, lv)
    integrand = cancel(together(phi ** 2 / f_of_r))
    pz, cz = lead(integrand, 0)
    pinf, _ = lead(integrand, 'inf')
    rsv = Q(6, cv)
    dphi = cancel(sp.diff(phi, r).subs(r, rsv) - cv / chi.eval(rsv))
    ok = ((p0_, c0_) == (-1, -phi_pow0) and (pi_, ci_) == (-(lv + 1), cv / cl2) and cl2 > 0 and (pz, cz) == (-1, 8)
          and pinf < -1 and chi.eval(rsv) > 0 and dphi == 0)
    run.check('E2.7[l=%d]' % lv, ok, 'phi0 at 0 %s, at inf %s (c/c_{l+2} = %s); phi0^2/f at 0 %s, at inf power %s; chi(r_s) %s; dphi %s' % (
        (p0_, c0_), (pi_, ci_), cv / cl2, (pz, cz), pinf, chi.eval(rsv), dphi))

# ------------------------------------------------------------------ E2.8 phi_0,xx = V-hat phi_0
psif = sp.Function('psi')(r)
phif = 1 / psif
res = sp.simplify(dx(dx(phif)) - Vhat_of_psi(psif) * phif)
run.check('E2.8', res == 0, 'generic psi: residual %s' % res)
for lv in range(2, 21):
    cv = c_of(lv)
    phi = (cv * r - 6) / chi_poly(lv).as_expr()
    res = cancel(together(dx(dx(phi)) - VH[lv] * phi))
    run.check('E2.8[l=%d]' % lv, res == 0, 'residual %s' % res)

# ------------------------------------------------------------------ L3: V-hat at infinity
n_, a1_, a2_, a3_ = symbols('n a1 a2 a3')
S = 1 + a1_ / r + a2_ / r ** 2 + a3_ / r ** 3
u = n_ / r + sp.diff(S, r) / S      # psi_r / psi for psi = r^n S
Vinf = cancel(together(Vhat_of_u(u)))
pw, cf = lead(Vinf, 'inf')
nn = n_ * (n_ + (2 if mut('e_nn1') else 1))   # mutation e_nn1: n(n+2) instead of n(n+1)
run.check('E3.1', pw == -2 and cancel(cf - nn) == 0, 'r^2 V-hat -> %s at power %s; want n(n+1) = %s' % (cf, pw, nn))
for lv in range(2, LMAX + 1):
    pp, cp = lead(PSI[lv], 'inf')
    pv, cv_ = lead(VH[lv], 'inf')
    wl = (lv + 1) * (lv + (3 if mut('e_lead_inf') else 2))   # mutation e_lead_inf: (l+1)(l+3)
    run.check('E3.2[l=%d]' % lv, pp == lv + 1 and pv == -2 and cv_ == wl and cv_ > 0,
              'deg psi0 %s (want %d); V-hat ~ %s r^%s (want %s r^-2)' % (pp, lv + 1, cv_, pv, wl))
run.note('elapsed %.1f s' % (time.time() - T0))
run.finish()
