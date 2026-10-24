"""Background equations of derive_general.py in a harmonic radial coordinate y (lam0 = nu0 + n beta + c, r0 = e^beta).
Prints the static field equations as ODEs in nu0(y), beta(y), with P = phi0' and charge q."""
import sympy as sp
import derive_general as d

y = d.x
n, q = d.n, d.q
nu, beta = sp.Function("nu")(y), sp.Function("beta")(y)
c = sp.Symbol("c", real=True)
eqs, E = d.setup(perturbed=False)
sub = {d.nu0: nu, d.lam0: nu + n * beta + c, d.r0: sp.exp(beta)}
for k, v in eqs.items():
    e = sp.simplify(v.subs(sub).doit())
    print(k, "=", sp.simplify(e))
