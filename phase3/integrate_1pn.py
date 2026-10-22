"""
Task 3-4 (numerical): integrate the 1PN (EIH, harmonic gauge) two-body problem with signed masses.

Two integrators of the same Lagrangian L_EIH(mu_A):
  * "exact EL":       solve  W(x,v) a = F(x,v)  with W = d^2L/dv dv  (not order-reduced).  E, P, J are
                      exact Noether charges of L (time/space translations, rotations), so they must be
                      conserved to integrator tolerance. L is Poincare invariant only to O(c^-2), so
                      K = G - P t is conserved only up to O(c^-4) terms.
  * "order-reduced":  Blanchet eq.(203) acceleration (verified identical to EL at 1PN in eih.py).
                      All charges conserved up to O(c^-4) drifts.
Units G = c = 1.  Body 1: mu1 > 0 (front), body 2: mu2 < 0 (rear).
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

# ------------------------------------------------------------------ symbolic Lagrangian -> numpy
mu1, mu2 = sp.symbols('mu1 mu2')
X = sp.symbols('x1 x2 x3 y1 y2 y3'); Vv = sp.symbols('u1 u2 u3 w1 w2 w3')
x1, x2 = sp.Matrix(X[:3]), sp.Matrix(X[3:]); v1, v2 = sp.Matrix(Vv[:3]), sp.Matrix(Vv[3:])
dx = x1 - x2; r = sp.sqrt(dx.dot(dx)); n = dx / r
L = (mu1 * v1.dot(v1) / 2 + mu2 * v2.dot(v2) / 2 + mu1 * mu2 / r
     + mu1 * v1.dot(v1)**2 / 8 + mu2 * v2.dot(v2)**2 / 8
     + mu1 * mu2 / (2 * r) * (3 * (v1.dot(v1) + v2.dot(v2)) - 7 * v1.dot(v2) - n.dot(v1) * n.dot(v2))
     - mu1 * mu2 * (mu1 + mu2) / (2 * r**2))
q = sp.Matrix(X); qd = sp.Matrix(Vv)
pvec = sp.Matrix([sp.diff(L, s) for s in Vv])
W = pvec.jacobian(qd)
F = sp.Matrix([sp.diff(L, s) for s in X]) - pvec.jacobian(q) * qd
E = (qd.T * pvec)[0] - L
P = pvec[:3, 0] + pvec[3:, 0]
J = x1.cross(pvec[:3, 0]) + x2.cross(pvec[3:, 0])
G = (mu1 * x1 * (1 + v1.dot(v1) / 2 - mu2 / (2 * r)) + mu2 * x2 * (1 + v2.dot(v2) / 2 - mu1 / (2 * r)))

def blanchet(ma, mb, nab, va, vb):
    return (-mb * nab / r**2
            + ((5 * ma * mb / r**3 + 4 * mb**2 / r**3
                + mb / r**2 * (sp.Rational(3, 2) * nab.dot(vb)**2 - va.dot(va) + 4 * va.dot(vb) - 2 * vb.dot(vb))) * nab
               + mb / r**2 * (4 * nab.dot(va) - 3 * nab.dot(vb)) * (va - vb)))
Ared = sp.Matrix.vstack(blanchet(mu1, mu2, n, v1, v2), blanchet(mu2, mu1, -n, v2, v1))

args = (X, Vv, mu1, mu2)
fW, fF = sp.lambdify(args, W, 'numpy'), sp.lambdify(args, F, 'numpy')
fA = sp.lambdify(args, Ared, 'numpy')
fE, fP, fJ, fG = (sp.lambdify(args, Q, 'numpy') for Q in (E, P, J, G))
# Newtonian-order energy (for scaling the residuals)
fEN = sp.lambdify(args, mu1 * v1.dot(v1) / 2 + mu2 * v2.dot(v2) / 2 - mu1 * mu2 / r, 'numpy')

def rhs(mode, m1, m2):
    def f(t, y):
        x, v = y[:6], y[6:]
        if mode == 'exact':
            a = np.linalg.solve(np.array(fW(x, v, m1, m2), float), np.array(fF(x, v, m1, m2), float).ravel())
        else:
            a = np.array(fA(x, v, m1, m2), float).ravel()
        return np.concatenate([v, a])
    return f

def run(m1, m2, d, tmax, mode, v0=(0, 0, 0, 0, 0, 0), npts=2001):
    y0 = np.array([d, 0, 0, 0, 0, 0, *v0], float)
    ts = np.linspace(0, tmax, npts)
    sol = solve_ivp(rhs(mode, m1, m2), (0, tmax), y0, t_eval=ts, method='DOP853', rtol=1e-13, atol=1e-15)
    x, v = sol.y[:6], sol.y[6:]
    out = dict(t=sol.t, x=x, v=v)
    out['E'] = np.array([fE(x[:, i], v[:, i], m1, m2) for i in range(len(ts))], float)
    out['P'] = np.array([np.ravel(fP(x[:, i], v[:, i], m1, m2)) for i in range(len(ts))], float)
    out['J'] = np.array([np.ravel(fJ(x[:, i], v[:, i], m1, m2)) for i in range(len(ts))], float)
    out['K'] = np.array([np.ravel(fG(x[:, i], v[:, i], m1, m2)) for i in range(len(ts))], float) - out['P'] * sol.t[:, None]
    out['sep'] = np.linalg.norm(x[:3] - x[3:], axis=0)
    out['mid'] = 0.5 * (x[0] + x[3])
    out['vmid'] = 0.5 * (v[0] + v[3])
    return out

def drift(arr):
    a = np.atleast_2d(np.asarray(arr).T).T if np.ndim(arr) == 1 else np.asarray(arr)
    return np.max(np.abs(a - a[0]))

def hyperbola_centres(o, m1, m2, i):
    """Event (t0, c) about which body A's instantaneous hyperbola is centred (1D motion along x)."""
    res = []
    for A, xi, vi in ((0, 0, 0), (1, 3, 3)):
        x, v = o['x'][xi, i], o['v'][vi, i]
        a = rhs('exact', m1, m2)(0, np.concatenate([o['x'][:, i], o['v'][:, i]]))[6 + vi]
        g2 = 1 / (1 - v * v)
        alpha = a * g2**1.5                      # proper acceleration
        res.append((o['t'][i] - np.sqrt(g2) * v / alpha, x - np.sqrt(g2) / alpha))
    return res

if __name__ == "__main__":
    np.set_printoptions(precision=4)
    m, d = 0.0143, 10.0
    aN = m / d**2
    print(f"Bondi dipole, mu1=+{m}, mu2=-{m}, d0={d}  (Shirokov's campaign values);  a_N = m/d^2 = {aN:.6e}")
    print(f"  m/d = {m/d:.3e}\n")

    # instantaneous accelerations at rest (exact-EL and order-reduced)
    y0 = np.array([d, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.])
    for mode in ('exact', 'reduced'):
        a = rhs(mode, m, -m)(0, y0)[6:]
        print(f"  [{mode:7s}] at rest: a1 d^2/m = {a[0]/aN:.8f}, a2 d^2/m = {a[3]/aN:.8f}, "
              f"mid = {(a[0]+a[3])/2/aN:.8f}   (1PN prediction 1 -/+ m/d = {1-m/d:.8f}, {1+m/d:.8f}, 1)")
    for V in (0.05, 0.1, 0.2):
        yV = y0.copy(); yV[6] = yV[9] = V
        a = rhs('reduced', m, -m)(0, yV)[6:]
        print(f"  common V={V}: mid a d^2/m = {(a[0]+a[3])/2/aN:.6f}  (1 - 5V^2/2 = {1-2.5*V*V:.6f});"
              f"  proper accel * d_rest^2/m = {(a[0]+a[3])/2*(1-V*V)**-1.5*(d/np.sqrt(1-V*V))**2/m:.6f}")
    print()

    # long runs
    tmax = 2000.0
    for label, m1, m2 in (("matched |mu2| = mu1", m, -m),
                          ("rigid  |mu2| = mu1(1+mu1/d)", m, -m * (1 + m / d)),
                          ("Shirokov-like mismatch |mu2| = 0.9961 mu1", m, -0.9961 * m),
                          ("mismatch |mu2| = 1.0039 mu1", m, -1.0039 * m)):
        print(f"--- {label}   (tmax = {tmax})")
        for mode in ('exact', 'reduced'):
            o = run(m1, m2, d, tmax, mode)
            En = abs(m1 * m2) / d
            print(f"  [{mode:7s}] drift E = {drift(o['E']):.2e} (E0 = {o['E'][0]:.4e}, scale m^2/d = {En:.2e}),"
                  f" P = {drift(o['P']):.2e}, J = {drift(o['J']):.2e}, K=G-Pt = {drift(o['K']):.2e}")
        o = run(m1, m2, d, tmax, 'exact')
        for tt in (200, 400, 1000, 2000):
            i = np.argmin(abs(o['t'] - tt))
            vm = o['vmid'][i]
            dX = o['mid'][i] - o['mid'][0]
            print(f"   t={tt:5.0f}: v_mid={vm:.4f} (v^2={vm*vm:.2e})  X_mid-X0={dX:8.3f}  d={o['sep'][i]:.6f}"
                  f"  d*gamma={o['sep'][i]/np.sqrt(1-vm*vm):.6f}")
        c = [hyperbola_centres(o, m1, m2, i) for i in (0, len(o['t']) // 2, -1)]
        print("   hyperbola centres (t0, c) body1 | body2 at t=0, tmax/2, tmax:")
        for cc in c:
            print(f"      ({cc[0][0]:+.3f}, {cc[0][1]:+.3f}) | ({cc[1][0]:+.3f}, {cc[1][1]:+.3f})")
        print()

    # scaling of the O(c^-4) residuals: K drift for the exact-EL integrator vs m/d at fixed v_final
    print("--- scaling of non-conservation with m/d (d = 10, run to t = 0.3 d^2/m, i.e. v ~ 0.3)")
    for mm in (0.0143, 0.05, 0.2):
        o_e = run(mm, -mm, d, 0.3 * d * d / mm, 'exact', npts=801)
        o_r = run(mm, -mm, d, 0.3 * d * d / mm, 'reduced', npts=801)
        print(f"  m/d={mm/d:.4f}: exact-EL: dE={drift(o_e['E']):.1e} dP={drift(o_e['P']):.1e} dK/(m d)={drift(o_e['K'])/(mm*d):.2e};"
              f"  reduced: dE/(m^2/d)={drift(o_r['E'])/(mm*mm/d):.2e} dP/m={drift(o_r['P'])/mm:.2e} dK/(m d)={drift(o_r['K'])/(mm*d):.2e}")

    # off-axis run (nonzero J) for the angular-momentum check
    print("\n--- off-axis: dipole with transverse velocity 0.01 on the negative body")
    for mode in ('exact', 'reduced'):
        o = run(m, -m, d, 1000.0, mode, v0=(0, 0, 0, 0, 0.01, 0))
        print(f"  [{mode:7s}] J0 = {o['J'][0]}, drift J = {drift(o['J']):.2e}, drift E = {drift(o['E']):.2e},"
              f" drift P = {drift(o['P']):.2e}")
