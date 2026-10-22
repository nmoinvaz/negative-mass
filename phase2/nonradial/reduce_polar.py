"""
Reduction of the separated polar equations (polar_eqs.pkl 'E', produced by derive_nonradial.py) to a
first-order system, done carefully at the level of jet variables.

  (t th), (t r), (r th)  -> H1', K', H0'   (first order; no F'', G'')
  (r r) after eliminating H0'', K'' (via d/dr of the first-order relations, which brings in F'', G''),
  KG_R, KG_I             -> linear system for {H0, F'', G''}, solved JOINTLY.
(Bug fixed relative to the first version, derive_nonradial.py -> gen_matrix*.py: there H0 was solved
from (rr) alone and still contained F'', G''; the code generator then silently set F'' = G'' = 0 inside
H0. polar_matrix_generated.py is regenerated here and the old matrix must not be used.)

Output: polar_matrix_generated.py  with y = (H1/s, K, F, F', G/s, G'/s), y' = J y, H0 = H2 = h.y.
Checks: J polynomial & even in s; vacuum limit = Zerilli (Moncrief map) + free KG (as gen_matrix2.py).
"""
import os, sys, pickle, time, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from lincommon import r, nu0, lam0, p0, s, w, eps, LL, LAM, MU, bgsub
T0 = time.time()
d = pickle.load(open(os.path.join(HERE, 'polar_eqs.pkl'), 'rb'))
E = d['E']
H0, H1, H2, K, F, G = [sp.Function(n)(r) for n in ('H0', 'H1', 'H2', 'K', 'F', 'G')]
D = lambda f, k=1: sp.diff(f, r, k)

def has2(z):
    return {nm: z.has(D(f, 2)) for nm, f in (('H0', H0), ('H1', H1), ('K', K), ('F', F), ('G', G))}

for k in ('tth', 'tr', 'rth', 'rr', 'KGR', 'KGI'):
    print(k, "second derivatives:", has2(E[k]))
fo = sp.solve([E['tth'], E['tr'], E['rth']], [D(H1), D(K), D(H0)], dict=True)[0]
fo = {k: sp.simplify(v) for k, v in fo.items()}
print("first-order relations: contain F'' or G''?", any(v.has(D(F, 2)) or v.has(D(G, 2)) for v in fo.values()),
      " %.0fs" % (time.time() - T0), flush=True)

def elim(z):
    """Replace H0'', K'', H1'' by d/dr of first-order relations, then H0', K', H1' -- F'', G'' kept."""
    for _ in range(3):
        z = z.subs({D(H0, 2): D(fo[D(H0)]), D(K, 2): D(fo[D(K)]), D(H1, 2): D(fo[D(H1)])})
        z = z.subs(fo)
        z = bgsub(z.doit())
    return sp.expand(z)

eqs3 = [elim(E[k]) for k in ('rr', 'KGR', 'KGI')]
for nm, z in zip(('rr', 'KGR', 'KGI'), eqs3):
    print(f"  {nm}: has H0'={z.has(D(H0))} K'={z.has(D(K))} H1'={z.has(D(H1))} F''={z.has(D(F, 2))} G''={z.has(D(G, 2))} H0={z.has(H0)}",
          flush=True)
sol = sp.solve(eqs3, [H0, D(F, 2), D(G, 2)], dict=True)
assert len(sol) == 1
sol = sol[0]
print("joint solve done %.0fs" % (time.time() - T0), flush=True)
H0sol = sol[H0]; Fpp = sol[D(F, 2)]; Gpp = sol[D(G, 2)]
for nm, z in (('H0', H0sol), ('F\'\'', Fpp), ('G\'\'', Gpp)):
    assert not any(z.has(D(f, 2)) for f in (F, G, H0, H1, K)), nm
    assert not any(z.has(D(f)) for f in (H0, H1, K)), nm
dH1 = fo[D(H1)].subs(H0, H0sol); dK = fo[D(K)].subs(H0, H0sol)
for z in (dH1, dK):
    assert not z.has(D(H0))

al, aa, P0, Uu = sp.symbols('alpha a phi0 u')
yv = sp.symbols('Ht1 K F dF Gt dGt')
def numerize(z):
    z = z.subs({D(F): yv[3], D(G): s*yv[5]})
    z = z.subs({F: yv[2], G: s*yv[4], H1: s*yv[0], K: yv[1]})
    z = z.subs(D(p0), Uu*p0)
    z = z.subs({sp.exp(nu0): al**2, sp.exp(lam0): aa**2})
    z = z.subs({nu0: 2*sp.log(al), lam0: 2*sp.log(aa)}).subs(p0, P0)
    z = sp.expand(sp.powsimp(sp.expand_power_exp(sp.expand(z)), force=True))
    assert not z.atoms(sp.Derivative) and not [a for a in z.atoms(sp.Function) if not isinstance(a, (sp.exp, sp.log))], z.atoms(sp.Function)
    return z
rows = [numerize(dH1)/s, numerize(dK), yv[3], numerize(Fpp), yv[5], numerize(Gpp)/s]
H0n = numerize(H0sol)
J = sp.zeros(6, 6)
for i, z in enumerate(rows):
    z = sp.expand(z)
    for j, v in enumerate(yv):
        J[i, j] = sp.cancel(sp.together(z.coeff(v)))
    assert sp.simplify(z - sum(J[i, j]*yv[j] for j in range(6))) == 0
hv = [sp.cancel(sp.together(sp.expand(H0n).coeff(v))) for v in yv]
print("matrix built %.0fs" % (time.time() - T0), flush=True)
bad = [(i, j) for i in range(6) for j in range(6) if sp.fraction(J[i, j])[1].has(s)]
print("entries with s in denominator:", bad, " h:", [j for j in range(6) if sp.fraction(hv[j])[1].has(s)])
odd = [(i, j) for i in range(6) for j in range(6) if sp.simplify(J[i, j].subs(s, -s) - J[i, j]) != 0]
print("entries not even in s:", odd)

# vacuum limit
M, Lm = sp.symbols('M lam')
f = 1 - 2*M/r
vac = {P0: 0, aa: 1/sp.sqrt(f), al: sp.sqrt(f)}
Jv = J.subs(vac).applyfunc(sp.simplify)
print("vacuum: metric-scalar blocks decouple:", all(Jv[i, j] == 0 for i in (0, 1) for j in (2, 3, 4, 5)),
      all(Jv[i, j] == 0 for i in (2, 3, 4, 5) for j in (0, 1)))
Z = sp.Function('Z')(r)
LLv = 2*Lm + 2
Kz = (Lm*(Lm + 1)*r**2 + 3*Lm*M*r + 6*M**2)/(r**2*(Lm*r + 3*M))*Z + f*D(Z)
Hz = (Lm*r**2 - 3*Lm*M*r - 3*M**2)/((r - 2*M)*(Lm*r + 3*M))*Z + r*D(Z)
VZ = 2*f*(Lm**2*r**2*((Lm + 1)*r + 3*M) + 9*M**2*(Lm*r + M))/(r**3*(Lm*r + 3*M)**2)
Zpp = sp.solve(sp.Eq(f*D(f*D(Z)), (s**2 + VZ)*Z), D(Z, 2))[0]
for i, lhs in enumerate((Hz, Kz)):
    res = D(lhs) - (Jv[i, 0].subs(LL, LLv)*Hz + Jv[i, 1].subs(LL, LLv)*Kz)
    print(f"vacuum Zerilli (Moncrief map) residual row {i}:", sp.simplify(res.subs(D(Z, 2), Zpp)))
Fv, Gv = sp.Function('Fv')(r), sp.Function('Gv')(r)
Zc = Fv + sp.I*s*Gv
B = Jv[2:, 2:].subs(LL, LLv)
yy = [Fv, D(Fv), Gv, D(Gv)]
rowF = sum(B[1, j]*yy[j] for j in range(4)); rowG = sum(B[3, j]*yy[j] for j in range(4))
kgc = -(s - sp.I*w)**2/f*Zc + D(r**2*f*D(Zc))/r**2 - LLv*Zc/r**2 - Zc
print("vacuum free-KG residual:", sp.simplify(sp.expand(kgc.subs({D(Fv, 2): rowF, D(Gv, 2): rowG}))))

from sympy.printing.numpy import NumPyPrinter
pr = NumPyPrinter({'fully_qualified_modules': False, 'allow_unknown_functions': True})
repl, red_ex = sp.cse(list(J) + list(hv))
with open(os.path.join(HERE, 'polar_matrix_generated.py'), 'w') as fh:
    fh.write('"""AUTO-GENERATED by reduce_polar.py. y = (H1/s, K, F, F\', G/s, G\'/s); y\' = J y; H0 = H2 = h . y.\n'
             'Background: alpha = e^{nu/2}, a = e^{lam/2}, phi0, u = phi0\'/phi0. Time dependence e^{st}; LL = l(l+1)."""\n')
    fh.write('import numpy as np\nfrom numpy import pi, exp, sqrt, log\n\n')
    fh.write('def J_and_h(r, alpha, a, phi0, u, omega, epsilon, s, LL, Lambda, mu):\n')
    for sym, ex in repl:
        fh.write(f'    {sym} = {pr.doprint(ex)}\n')
    fh.write('    J = np.empty((6, 6), dtype=np.result_type(s, float))\n')
    for k in range(36):
        fh.write(f'    J[{k//6}, {k%6}] = {pr.doprint(red_ex[k])}\n')
    fh.write('    h = np.array([' + ', '.join(pr.doprint(red_ex[36 + j]) for j in range(6)) + '])\n')
    fh.write('    return J, h\n')
pickle.dump({'J': J, 'h': hv}, open(os.path.join(HERE, 'polar_matrix_final.pkl'), 'wb'))
print("wrote polar_matrix_generated.py %.0fs" % (time.time() - T0))
