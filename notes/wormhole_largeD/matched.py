"""Large-D expansion of the Ellis wormhole's growth rate by matched asymptotics.

Near the throat use R = r^s, s = 2(n - 1) = 2(D - 3), eps = 1/s, l = ln R, r = e^{eps l}; W^2 = s + Omega with
Omega = Omega0 + eps Omega1 + ...  The master equation (modes_r.py) is exactly
  Y_R = eps L r / R,
  (R - 1) L_R = (1 + eps Omega) r Y + L/(2R) + (1/2 + eps) L (R - 1)/R - (1/2 + eps) Y (R - 1)/(R r) - Y/(2 r R).
Expanding Y = sum Y_k eps^k, L = sum L_k eps^k gives at each order (R - 1) L_k' - L_k/2 = h_k(R). The solution
analytic at the throat is L_k = -2 h_k(1) + (R - 1)^{1/2} int_1^R (h_k - h_k(1)) (R' - 1)^{-3/2} dR'. At large R it
contains G_k R^{1/2} with G_k = int_1^inf (h_k - h_k(1)) (R - 1)^{-3/2} dR, and R^{1/2} = r^{n-1} is the far zone's
growing solution, so G_k = 0 at every order. G_1 = 0 identically; G_2 = 0 fixes Omega0 (Omega1 only shifts h_2 by a
constant, which drops out).
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 30
R = sp.Symbol("R", positive=True)
eps, Om0, Om1 = sp.symbols("epsilon Omega0 Omega1")
Y0, Y1, Y2, L0, L1, L2 = sp.symbols("Y0 Y1 Y2 L0 L1 L2")
lg = sp.log(R)


def h_coeffs():
    r = sp.exp(eps * lg)
    Y = Y0 + eps * Y1 + eps**2 * Y2
    L = L0 + eps * L1 + eps**2 * L2
    Om = Om0 + eps * Om1
    rhs = ((1 + eps * Om) * r * Y + L / (2 * R) + (sp.Rational(1, 2) + eps) * L * (R - 1) / R
           - (sp.Rational(1, 2) + eps) * Y * (R - 1) / (R * r) - Y / (2 * r * R))
    ser = sp.series(rhs, eps, 0, 3).removeO()
    out = []
    for k, Lk in enumerate((L0, L1, L2)):
        c = sp.expand(ser.coeff(eps, k))
        # (R-1) L_k' = c ;  c contains L_k/2 (from L/(2R) + L(R-1)/(2R)); h_k = c - L_k/2
        hk = sp.simplify(c - Lk / 2)
        assert sp.simplify(sp.diff(hk, Lk)) == 0, (k, hk)
        out.append(hk)
    # Y_k' = coefficient of eps^k in eps L r / R
    yser = sp.series(eps * L * r / R, eps, 0, 3).removeO()
    ys = [sp.expand(yser.coeff(eps, k)) for k in range(3)]
    return out, ys


if __name__ == "__main__":
    hs, ys = h_coeffs()
    print("h0 =", hs[0]); print("h1 =", hs[1]); print("h2 =", hs[2])
    print("Y0' =", ys[0], "  Y1' =", ys[1], "  Y2' =", ys[2])
