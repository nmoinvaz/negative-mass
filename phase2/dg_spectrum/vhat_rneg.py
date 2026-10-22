"""r_neg(l) = sup{r : V^(r) < 0}, recomputed here (cf. ../vhat_negativity.py) to 1e-6 by root polishing."""
import numpy as np, sympy as sp
import common as c

def _rneg(l):
    num, _ = sp.fraction(sp.cancel(c.symbolic(l)['Vh']))
    rts = [float(sp.re(z)) for z in sp.Poly(num, c.rs_).nroots(n=30) if abs(sp.im(z)) < 1e-20 and sp.re(z) > 0]
    return max(rts)

RNEG = {l: _rneg(l) for l in range(2, 9)}
if __name__ == '__main__':
    print(RNEG)
