"""Jacobi / ground-state identity used in all proofs:
   for V = phi''/phi and u = phi h:   u'^2 + V u^2 = phi^2 h'^2 + (phi phi' h^2)'   (pointwise, ' = d/dx).
Also the boundary identity for the Robin form and the zero-energy trial-function computation."""
import sympy as sp
x = sp.Symbol('x')
phi = sp.Function('phi')(x); h = sp.Function('h')(x)
u = phi*h
V = sp.diff(phi, x, 2)/phi
lhs = sp.diff(u, x)**2 + V*u**2
rhs = phi**2*sp.diff(h, x)**2 + sp.diff(phi*sp.diff(phi, x)*h**2, x)
print("u'^2 + V u^2 - [phi^2 h'^2 + (phi phi' h^2)'] =", sp.simplify(lhs - rhs))
# trial function u = phi on an interval [a,b] with V = phi''/phi: integrand phi'^2 + phi phi'' = (phi phi')'
print("phi'^2 + V phi^2 - (phi phi')' =", sp.simplify(sp.diff(phi, x)**2 + V*phi**2 - sp.diff(phi*sp.diff(phi, x), x)))
