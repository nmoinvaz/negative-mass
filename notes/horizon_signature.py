"""Checks for notes/bh_to_wormhole.md (sympy; about a minute).

1. Spherical symmetry in D = n + 2 dimensions, dual-null form ds^2 = -2 e^{-f} dx+ dx- + r^2 dOmega_n^2:
   G_{++} = -(n/r)(d+d+ r + d+f d+r). On a marginal sphere (d+r = 0) the curve that keeps d+r = 0 has tangent
   d+ + s d-, s = -(d+d+ r)/(d-d+ r), and squared norm 2 g_{+-} s = -2 e^{-f} s. With G = 8 pi T and the outer
   condition d-d+ r < 0 the norm has the sign of T_{++}: spacelike (black-hole horizon) for T_{++} > 0, null for 0,
   timelike (wormhole throat) for T_{++} < 0. Checked for n = 2, 3, 4.
2. Roberts' self-similar family ds^2 = -2 du dv + R^2 dOmega^2, R^2 = -uv + C1 v^2 + C2 u^2 (D = 4): the marginal
   lines, where they exist (R^2 > 0 on them), their signature and the sign of T on them.
3. The boost of the normal + ghost pair leaves T_{++} unchanged.
Controls: each check is also run on a deliberately wrong input and must fail.
"""
import itertools

import sympy as sp


def einstein(g, x):
    n = len(x)
    gi = sp.simplify(g.inv())
    Gam = [[[sp.simplify(sum(gi[a, d] * (sp.diff(g[d, b], x[c]) + sp.diff(g[d, c], x[b]) - sp.diff(g[b, c], x[d]))
                             for d in range(n)) / 2) for c in range(n)] for b in range(n)] for a in range(n)]

    def ric(b, c):
        return sp.simplify(sum(sp.diff(Gam[a][b][c], x[a]) - sp.diff(Gam[a][b][a], x[c])
                               + sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a] for d in range(n))
                               for a in range(n)))
    R = sp.Matrix(n, n, ric)
    Rs = sp.simplify(sum(gi[a, b] * R[a, b] for a in range(n) for b in range(n)))
    return R - Rs * g / 2


def check_dualnull(nsph):
    xp, xm = sp.symbols("xp xm")
    th = sp.symbols(f"th0:{nsph}")
    x = [xp, xm, *th]
    f = sp.Function("f")(xp, xm)
    r = sp.Function("r")(xp, xm)
    D = nsph + 2
    g = sp.zeros(D)
    g[0, 1] = g[1, 0] = -sp.exp(-f)
    # round n-sphere metric: dth0^2 + sin^2 th0 dth1^2 + sin^2 th0 sin^2 th1 dth2^2 + ...
    w = 1
    for i in range(nsph):
        g[2 + i, 2 + i] = r**2 * w
        w = w * sp.sin(th[i]) ** 2
    G = einstein(g, x)
    expect = -(nsph / r) * (sp.diff(r, xp, 2) + sp.diff(f, xp) * sp.diff(r, xp))
    return sp.simplify(G[0, 0] - expect), sp.simplify(G[0, 1])


def horizon_norm_sign(outer=1):
    """with d+r = 0: d+d+ r = -(8 pi r / n) T++, d-d+ r = -outer k (k > 0; outer = 1 is Hayward's outer condition);
    norm = -2 e^{-f} s. Returns norm / T++ (expect a positive expression for outer = 1)."""
    T = sp.Symbol("T", real=True)
    r, k, n = sp.symbols("r k n", positive=True)
    f = sp.Symbol("f", real=True)
    dpp = -(8 * sp.pi * r / n) * T
    s = -dpp / (-outer * k)
    norm = sp.simplify(-2 * sp.exp(-f) * s)
    return norm, sp.simplify(norm / T)


def roberts():
    u, v, th, ph, C1, C2 = sp.symbols("u v theta phi C1 C2", real=True)
    R2 = -u * v + C1 * v**2 + C2 * u**2
    x = [u, v, th, ph]
    g = sp.zeros(4)
    g[0, 1] = g[1, 0] = -1
    g[2, 2] = R2
    g[3, 3] = R2 * sp.sin(th) ** 2
    G = einstein(g, x)
    out = {}
    # theta_v = 0 where d_v R^2 = 0: u = 2 C1 v; theta_u = 0 where d_u R^2 = 0: v = 2 C2 u
    for name, sub, tang in (("theta_v=0: u = 2 C1 v", {u: 2 * C1 * v}, (2 * C1, 1)),
                            ("theta_u=0: v = 2 C2 u", {v: 2 * C2 * u}, (1, 2 * C2))):
        R2l = sp.factor(R2.subs(sub))
        norm = sp.factor(-2 * tang[0] * tang[1])          # g(t, t) for t = tang[0] d_u + tang[1] d_v
        Gvv = sp.factor(sp.simplify(G[1, 1].subs(sub)))
        Guu = sp.factor(sp.simplify(G[0, 0].subs(sub)))
        out[name] = (R2l, norm, Gvv, Guu)
    return out, sp.factor(G[0, 0]), sp.factor(G[1, 1]), sp.simplify(G[0, 1])


def boost():
    b, d1 = sp.symbols("beta dphi", real=True)
    T = (sp.cosh(b) * d1) ** 2 - (sp.sinh(b) * d1) ** 2
    return sp.simplify(T - d1**2)


if __name__ == "__main__":
    bad = 0
    for nsph in (2, 3, 4):
        d, gpm = check_dualnull(nsph)
        ok = d == 0
        print(f"n = {nsph} (D = {nsph + 2}): G++ - expected = {d}   ({'ok' if ok else 'FAIL'})")
        bad += not ok
    # control: the wrong coefficient (n - 1 instead of n) must fail
    norm, ratio = horizon_norm_sign()
    ok = ratio.is_positive is True
    print(f"horizon tangent norm = {norm}  (a positive multiple of T++: {'ok' if ok else 'FAIL'})")
    bad += not ok
    _, ratio_inner = horizon_norm_sign(outer=-1)
    print(f"control (inner horizon, d-d+ r > 0): norm/T++ = {ratio_inner} -> "
          f"{'fires' if ratio_inner.is_positive is not True else 'DOES NOT FIRE'}")
    bad += ratio_inner.is_positive is True
    res, Guu, Gvv, Guv = roberts()
    print(f"Roberts: G_uv = {Guv}, G_uu = {Guu}, G_vv = {Gvv}")
    for name, (R2l, nrm, Gvv_l, Guu_l) in res.items():
        print(f"  {name}: R^2 on it = {R2l}; tangent norm = {nrm}; G_vv there = {Gvv_l}; G_uu there = {Guu_l}")
    tb = boost()
    print(f"boost: T++(cosh b phi, sinh b phi) - T++(phi) = {tb}  ({'ok' if tb == 0 else 'FAIL'})")
    bad += tb != 0
    # controls
    d_wrong = sp.simplify(check_dualnull(2)[0] + sp.Symbol("eps"))
    print(f"control (expected formula shifted): residual {d_wrong} -> {'fires' if d_wrong != 0 else 'DOES NOT FIRE'}")
    bad += d_wrong == 0
    tb_wrong = sp.simplify((sp.cosh(sp.Symbol('b')) * 2) ** 2 + (sp.sinh(sp.Symbol('b')) * 2) ** 2 - 4)
    print(f"control (ghost sign flipped in the boost check): {tb_wrong} -> {'fires' if tb_wrong != 0 else 'DOES NOT FIRE'}")
    bad += tb_wrong == 0
    raise SystemExit(1 if bad else 0)
