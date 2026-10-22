"""r_neg(l), r_s(l), r_*(l) and min V^ for l up to 20 (exact rational V^, real-root isolation)."""
import sympy as sp, numpy as np
import common as c
r = c.rs_
print(f"{'l':>3} {'k':>6} {'r_s':>9} {'r_*':>9} {'r_neg':>9} {'l^2 r_neg':>10} {'chi>0':>6}")
for l in list(range(2, 13)) + [16, 20]:
    s = c.symbolic(l)
    num, _ = sp.fraction(sp.cancel(s['Vh']))
    P = sp.Poly(num, r)
    roots = [float(iv[0][1]) for iv in P.intervals(eps=1e-12) if iv[0][1] > 0]
    rneg = max(roots)
    numk, _ = sp.fraction(sp.cancel(s['dlogphi0_dx']))
    rst = min(float(iv[0][1]) for iv in sp.Poly(numk, r).intervals(eps=1e-12) if iv[0][0] > c.r_s(l))
    chipos = sp.Poly(sp.cancel(s['chi'] / r), r).count_roots(0, None) == 0
    print(f"{l:>3} {c.kval(l):>6} {c.r_s(l):>9.5f} {rst:>9.5f} {rneg:>9.5f} {l*l*rneg:>10.4f} {str(chipos):>6}")
