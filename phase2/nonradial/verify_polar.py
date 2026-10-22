"""
Numerical verification that the reduced polar system  y' = J y  (polar_matrix_generated.py, from reduce_polar.py) satisfies
ALL eight separated linearized equations (tt, tr, rr, t-theta, r-theta, angular trace, KG_R, KG_I)
as derived in derive_nonradial.py (polar_eqs.pkl, H2 = H0 already imposed).

At sample radii on a real star background (eps = +1 and -1), for random state vectors y and several s:
    y' = J y,  y'' = J' y + J y'  (J' by central differences along the background, h = 1e-4),
    H1 = s y0, G = s y4, H0 = h.y (and its derivatives likewise).
Each residual is reported relative to the largest single term of that equation (|residual| / max|term|),
and compared with a control where y'' is perturbed by 1% (must give O(1e-2) residual).
This complements the exact random-rational check in derive_nonradial.py, which completed for
(t-theta, t-r, r-theta) before the run was interrupted.
"""
import os, sys, pickle, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from lincommon import r, nu0, lam0, p0, s, w, eps, LL, LAM, MU
from polar_matrix_generated import J_and_h
import background as B

d = pickle.load(open(os.path.join(HERE, 'polar_eqs.pkl'), 'rb'))
E = d['E']; H2sol = d['H2sol']
fn = {n: sp.Function(n)(r) for n in ('H0', 'H1', 'H2', 'K', 'F', 'G')}
syms = {}
rep = {}
for n, f in fn.items():
    for k in (3, 2, 1):
        sy = sp.Symbol(f'{n}_d{k}'); rep[sp.diff(f, r, k)] = sy; syms[(n, k)] = sy
    sy = sp.Symbol(f'{n}_d0'); syms[(n, 0)] = sy
NU, LA, P, DP = sp.symbols('NU LA P DP')
def prep(z):
    z = z.subs(fn['H2'], fn['H0'])
    z = z.subs(rep)
    z = z.subs({fn[n]: syms[(n, 0)] for n in fn})
    z = z.subs(sp.diff(p0, r), DP).subs(p0, P).subs({nu0: NU, lam0: LA})
    return z

def split_terms(z):
    return sp.Add.make_args(sp.expand(z))

keys = ('tt', 'tr', 'rr', 'tth', 'rth', 'ang_tr', 'KGR', 'KGI')
argsyms = [r, NU, LA, P, DP, s, w, eps, LL, LAM, MU] + [syms[(n, k)] for n in ('H0', 'H1', 'K', 'F', 'G') for k in range(3)]
funcs = {}
for k in keys:
    z = prep(E[k])
    leftover = [a for a in z.atoms(sp.Function, sp.Derivative) if not isinstance(a, sp.exp)]
    assert not leftover, (k, leftover)
    terms = split_terms(z)
    funcs[k] = (sp.lambdify(argsyms, z, 'numpy'), sp.lambdify(argsyms, sp.Matrix(terms), 'numpy'), len(terms))

def check(star, l, svals, radii, rng):
    LLv = l * (l + 1); h = 1e-4
    out = []
    for rv in radii:
        def Jh(rr, sv):
            al, a, phi, u = star.bg(np.array([rr]))
            return J_and_h(rr, al[0], a[0], phi[0], u[0], star.w, star.eps, sv, LLv, B.LAM, B.MU)
        al, a, phi, u = star.bg(np.array([rv]))
        for sv in svals:
            J0, h0 = Jh(rv, sv); Jp, hp = Jh(rv + h, sv); Jm, hm = Jh(rv - h, sv)
            dJ = (Jp - Jm) / (2 * h); dh = (hp - hm) / (2 * h); ddh = (hp - 2 * h0 + hm) / h ** 2
            y = rng.standard_normal(6)
            y1 = J0 @ y; y2 = dJ @ y + J0 @ y1
            for ctrl in (False, True):
                yy2 = y2 * (1.01 if ctrl else 1.0)
                H0 = (h0 @ y, dh @ y + h0 @ y1, ddh @ y + 2 * dh @ y1 + h0 @ yy2)
                vals = dict(H0=H0, H1=(sv * y[0], sv * y1[0], sv * yy2[0]), K=(y[1], y1[1], yy2[1]),
                            F=(y[2], y1[2], yy2[2]), G=(sv * y[4], sv * y1[4], sv * yy2[4]))
                args = [rv, np.log(al[0] ** 2), np.log(a[0] ** 2), phi[0], u[0] * phi[0], sv, star.w, star.eps, LLv,
                        B.LAM, B.MU] + [vals[n][k] for n in ('H0', 'H1', 'K', 'F', 'G') for k in range(3)]
                row = []
                for k in keys:
                    f, ft, nt = funcs[k]
                    res = f(*args); tm = np.max(np.abs(np.array(ft(*args), dtype=float)))
                    row.append(abs(res) / tm)
                out.append((rv, sv, ctrl, row))
    return out

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    for epsv, wv in ((-1, 0.7603), (+1, 0.75)):
        st = B.solve_star(wv, epsv)
        for l in (2, 3):
            res = check(st, l, [0.05, 0.7, 5.0], [0.3, 1.5, 3.0, 6.0], rng)
            good = np.array([r_[3] for r_ in res if not r_[2]]); bad = np.array([r_[3] for r_ in res if r_[2]])
            print(f"eps={epsv:+d} l={l}: max relative residual per eq (true y'')      :",
                  ' '.join(f"{k}={v:.1e}" for k, v in zip(keys, good.max(axis=0))))
            print(f"            control (y'' perturbed 1%), min per eq           :",
                  ' '.join(f"{k}={v:.1e}" for k, v in zip(keys, bad.min(axis=0))))
