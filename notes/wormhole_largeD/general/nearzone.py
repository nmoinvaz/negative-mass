"""Large-D limit of the general master equation in the near zone (the whole interval -pi/2 < y < pi/2).

With a = 1 (throat areal radius) and N1 = n - 1: e^c = e^{-nu_th}/(N1 cos y_th), and
  r = e^{B/N1},  B(y) = nu_th - nu(y) + ln(cos y_th / cos y)   (B = B' = 0 at the throat),
  W^2 = 2 N1 w e^{2 nu_th}  (w = W_th^2 a^2 / s, W_th the growth rate in proper time at the throat).
First-order system: Y' = nu' Y + r' L (the tx equation) and L' from the Ricci-form angular equation.
Prints L' = F(Y, L) expanded in 1/N1 at fixed y, for generic nu(y)."""
import sympy as sp

import derive_general as d

y = d.x
N1 = sp.Symbol("N1", positive=True)
w = sp.Symbol("w", positive=True)
nuth, yth = sp.symbols("nu_th y_th", real=True)
eps2, Lam2 = sp.symbols("eps2 Lam2", real=True)
nu = sp.Function("nu")(y)
Bf = nuth - nu + sp.log(sp.cos(yth) / sp.cos(y))
first = {k: d.first(v, d.setup()[1]) for k, v in d.setup()[0].items()}
Yv, Lv = sp.symbols("Yv Lv")
Lfun = d.L
eq = first["angR"]
# background
n = N1 + 1
cc = -nuth - sp.log(N1) - sp.log(sp.cos(yth))
r = sp.exp(Bf / N1)
lam = nu + n * Bf / N1 + cc
P2 = n * (1 + Lam2) / (n - 1)
q2 = n * eps2 * sp.exp(-2 * cc) / (n - 1)
W2 = 2 * N1 * w * sp.exp(2 * nuth)
# eliminate N via the scalar relation, Y'' via the tx relation
Yp = sp.diff(nu, y) * d.Y + sp.diff(d.r0, y) * Lfun
Ns = Lfun - d.n * d.Y / d.r0
e = eq.subs(d.N, Ns).doit()
e = e.subs(sp.Derivative(d.Y, (y, 2)), sp.diff(Yp, y)).subs(sp.Derivative(d.Y, y), Yp)
e = e.subs(sp.Derivative(d.Y, y), Yp)
Lp = sp.solve(e, sp.Derivative(Lfun, y))[0]
Lp = Lp.subs({d.nu0: nu, d.lam0: lam, d.r0: r}).doit()
Lp = Lp.subs({d.P**2: P2, d.q**2: q2, d.n: n, d.W: sp.sqrt(W2)})
Lp = Lp.subs({d.Y: Yv, Lfun: Lv})
Lp = sp.simplify(Lp)
# leading orders in 1/N1 at fixed y: series in x = 1/N1
xs = sp.Symbol("xs", positive=True)
ser = sp.series(Lp.subs(N1, 1 / xs), xs, 0, 2).removeO()
print("L' leading (xs = 1/N1):")
print(sp.simplify(ser.coeff(xs, -1)), " [order N1]")
print(sp.simplify(ser.coeff(xs, 0)), " [order 1]")
