#!/usr/bin/env python3
"""Residual stage of #161's Theorem E with explicit scalars (issue #171)

planar_closure.md (#161) writes R = a a - b b for g = -2 exp(2 H(p) J(q)) dp dq + dy^2 + dz^2, with
a = (c dH - (4/c) dJ)/2 and b = (c dH + (4/c) dJ)/2, and says that the metric is a solution with one normal scalar
phi_1 ~ c H - (4/c) J and one phantom phi_2 ~ c H + (4/c) J, "by hand; the residual stage ... is not computed". This
script computes G_ab - 8 pi (T_1 - T_2)_ab for H(p), J(q) undefined and c a symbol, finds the normalizations and the
potential, and runs controls. Conventions (signature -+++, T = sum s_i (d phi_i d phi_i - g (d phi_i)^2 / 2) - g V with
(s_1, s_2) = (+1, -1)) are those of two_scalars.md, Setup.

    python3 notes/warp/planar_residual.py                 all checks, then every control (exit 1 on any failure)
    python3 notes/warp/planar_residual.py --mutate=NAME   one control (its target checks must fail, exit 1)
    python3 notes/warp/planar_residual.py --no-controls   the checks only
    python3 notes/warp/planar_residual.py --checks=a,b    only the named checks (development)

Imports nothing from the repository: the Einstein tensor and the stress tensor are written here, and share only sympy with
phantom_source.py, two_scalars.py and planar_closure.py. P1 and P2 are retyped from same_sign_scalars.md. Zero is
tested by sympy.simplify(e) == 0 and nothing weaker.
"""
import itertools
import re
import subprocess
import sys
import time

import sympy as sp

# name: (what the mutation does)
MUTATIONS = {
    "scale11-normal": "phi_1 times 11/10 in E_residual, E_rform, E_explicit",
    "scale11-phantom": "phi_2 times 11/10 in E_residual, E_rform, E_explicit",
    "flip-phantom": "the phantom's kinetic sign flipped (T_1 + T_2, and a a + b b in E_rform)",
    "swap-roles": "normal and phantom exchanged: phi_1 phantom, phi_2 normal",
    "wrong-partner": "phi_2 = alpha_2 (c H + (2/c) J)",
    "potential": "V = 1/(80 pi) in E_residual",
    "ricci-sign": "the Ricci tensor with the opposite sign in the Einstein routine",
    "conv-flip": "the sign of T flipped in the stress function",
    "p1-sign": "one kinetic sign flipped in P1",
    "p2-scale11": "one of P2's fields times 11/10",
    "flat-flip": "the phantom's sign flipped in planted_flat",
}

# name: checks that must fail (declared in the raw README after the development run of step 3)
TARGETS = {
    "scale11-normal": ["E_residual", "E_rform", "E_explicit"],
    "scale11-phantom": ["E_residual", "E_rform", "E_explicit"],
    "flip-phantom": ["E_solve", "E_residual", "E_rform", "E_explicit"],
    "swap-roles": ["E_solve", "E_residual", "E_rform", "E_explicit"],
    "wrong-partner": ["E_solve", "E_residual", "E_rform", "E_explicit"],
    "potential": ["E_residual"],
    "ricci-sign": ["planted_p1", "planted_p2", "E_solve", "E_residual", "E_rform", "E_explicit"],
    "conv-flip": ["conv_energy", "planted_p1", "planted_p2", "E_solve", "E_residual", "E_explicit"],
    "p1-sign": ["planted_p1"],
    "p2-scale11": ["planted_p2"],
    "flat-flip": ["planted_flat"],
}

RESULTS = []
MUT = None
PI = sp.pi
simp = sp.simplify


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# Einstein and stress tensors, and the residual. One function for every case.
# ---------------------------------------------------------------------------------------------------------------
def geometry(g, X):
    """Christoffels, Ricci R_bc = d_a G^a_bc - d_c G^a_ab + G^a_ad G^d_bc - G^a_cd G^d_ab, R, G_ab = R_ab - R g_ab / 2."""
    n = len(X)
    gi = g.inv().applyfunc(simp)
    Gam = [[[simp(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                      for d in range(n)) / 2) for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.zeros(n, n)
    for b in range(n):
        for c in range(b, n):
            s = sum(sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][a][b], X[c]) for a in range(n))
            s += sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][a][b] for a in range(n) for d in range(n))
            Ric[b, c] = Ric[c, b] = simp(s)
    if MUT == "ricci-sign":
        Ric = -Ric
    R = simp(sum(gi[a, b] * Ric[a, b] for a in range(n) for b in range(n)))
    G = (Ric - R * g / 2).applyfunc(simp)
    return {"g": g, "gi": gi, "X": X, "Ric": Ric, "R": R, "G": G}


def grad(geo, f):
    return sp.Matrix([sp.diff(f, x) for x in geo["X"]])


def stress(geo, fields, signs, V):
    """T_ab = sum_i s_i (d_a phi_i d_b phi_i - g_ab (d phi_i)^2 / 2) - g_ab V, as in two_scalars.md (Setup)."""
    g, gi = geo["g"], geo["gi"]
    T = sp.zeros(*g.shape)
    for s, f in zip(signs, fields):
        d = grad(geo, f)
        T += s * (d * d.T - g * (d.T * gi * d)[0] / 2)
    T = T - g * V
    return -T if MUT == "conv-flip" else T


def residual(geo, fields, signs, V):
    """G_ab - 8 pi T_ab, all components simplified."""
    return (geo["G"] - 8 * PI * stress(geo, fields, signs, V)).applyfunc(simp)


def nonzero(M):
    """The components of M that sympy.simplify does not reduce to 0, as (i, j, value)."""
    out = []
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = simp(M[i, j])
            if v != 0:
                out.append((i, j, v))
    return out


def brief(bad, k=3):
    return "; ".join(f"[{i},{j}] = {v}" for i, j, v in bad[:k]) + (" ..." if len(bad) > k else "")


# ---------------------------------------------------------------------------------------------------------------
# The Theorem E metric and its fields
# ---------------------------------------------------------------------------------------------------------------
p, q, y, z = sp.symbols("p q y z", real=True)
c = sp.Symbol("c", real=True, nonzero=True)
Hf, Jf = sp.Function("H")(p), sp.Function("J")(q)
XE = (p, q, y, z)
ALPHA = 1 / (4 * sp.sqrt(2 * PI))  # 1 / (2 sqrt(8 pi)), found by E_solve


def metric_E(H, J):
    g = sp.zeros(4, 4)
    g[0, 1] = g[1, 0] = -sp.exp(2 * H * J)
    g[2, 2] = g[3, 3] = 1
    return g


def fields_E(H, J, cc, a1, a2, scale=True):
    """phi_1 = a1 (c H - (4/c) J) normal, phi_2 = a2 (c H + (4/c) J) phantom, up to the mutation."""
    k2 = 2 if MUT == "wrong-partner" else 4
    s1 = sp.Rational(11, 10) if scale and MUT == "scale11-normal" else 1
    s2 = sp.Rational(11, 10) if scale and MUT == "scale11-phantom" else 1
    return [s1 * a1 * (cc * H - 4 * J / cc), s2 * a2 * (cc * H + k2 * J / cc)]


def signs_E():
    return {"flip-phantom": (1, 1), "swap-roles": (-1, 1)}.get(MUT, (1, -1))


_GEO_E = []


def geo_E():
    if not _GEO_E:
        _GEO_E.append(geometry(metric_E(Hf, Jf), XE))
    return _GEO_E[0]


# ---------------------------------------------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------------------------------------------
def check_conv_energy():
    t0 = time.time()
    t, x, yy, zz = sp.symbols("t x yy zz", real=True)
    w = sp.Symbol("omega", positive=True)
    geo = geometry(sp.diag(-1, 1, 1, 1), (t, x, yy, zz))
    got = {}
    for name, s in (("normal", 1), ("phantom", -1)):
        T = stress(geo, [w * t], [s], 0)
        got[name] = (simp(T[0, 0]), simp(T[1, 1]))
    ok = got["normal"] == (w**2 / 2, w**2 / 2) and got["phantom"] == (-w**2 / 2, -w**2 / 2)
    check("conv_energy", ok, f"Minkowski, phi = omega t: normal (T_tt, T_xx) = {got['normal']}, phantom "
          f"{got['phantom']}; two_scalars.md: normal +, phantom - ({time.time() - t0:.0f} s)")


def sheared(kind):
    """P1 (kind 1) and P2 (kind 2) of same_sign_scalars.md in the sheared coordinates (T, Z, R, Th)."""
    T_, Z, R, Th = sp.symbols("T Z R Th", real=True)
    X = (T_, Z, R, Th)
    tau, zz, Thh = T_ + Z / 3, Z, Th + T_ / 7 + Z / 5
    Jm = sp.Matrix([[sp.diff(f, x) for x in X] for f in (tau, zz, R, Thh)])
    d = zz if kind == 1 else tau
    g = (Jm.T * sp.diag(-1 / d**2, 1 / d**2, 1, R**2) * Jm).applyfunc(simp)
    fields = [R * sp.cos(Thh) / sp.sqrt(8 * PI), R * sp.sin(Thh) / sp.sqrt(8 * PI)]
    return geometry(g, X), fields


def check_planted(kind):
    t0 = time.time()
    name = f"planted_p{kind}"
    geo, fields = sheared(kind)
    signs = (1, 1) if kind == 1 else (-1, -1)
    V = -1 / (8 * PI) if kind == 1 else 1 / (8 * PI)
    if MUT == "p1-sign" and kind == 1:
        signs = (1, -1)
    if MUT == "p2-scale11" and kind == 2:
        fields = [fields[0], sp.Rational(11, 10) * fields[1]]
    bad = nonzero(residual(geo, fields, signs, V))
    check(name, not bad, f"G_ab - 8 pi T_ab, {'two normal scalars, V = -1/(8 pi)' if kind == 1 else 'two phantoms, V = 1/(8 pi)'}, "
          f"fields X_I / sqrt(8 pi): {'all 16 components 0' if not bad else brief(bad)} ({time.time() - t0:.0f} s)")


def check_planted_flat():
    t0 = time.time()
    t, x, yy, zz = sp.symbols("t x yy zz", real=True)
    geo = geometry(sp.diag(-1, 1, 1, 1), (t, x, yy, zz))
    f = t**2 * x + sp.sin(yy) * sp.exp(zz)
    signs = (1, 1) if MUT == "flat-flip" else (1, -1)
    bad = nonzero(residual(geo, [f, f], signs, 0))
    T1 = nonzero(stress(geo, [f], [1], 0))
    check("planted_flat", not bad and len(T1) > 0, f"Minkowski, phi_1 = phi_2 = t^2 x + sin(y) e^z, signs {signs}, V = 0: "
          f"{'all 16 components 0' if not bad else brief(bad)}; T_1 has {len(T1)} nonzero components ({time.time() - t0:.0f} s)")


def check_E_solve():
    t0 = time.time()
    geo = geo_E()
    A1, A2, Vh = sp.symbols("A1 A2 Vh", real=True)
    fl = fields_E(Hf, Jf, c, sp.sqrt(A1), sp.sqrt(A2), scale=False)
    Res = residual(geo, fl, signs_E(), Vh)
    h, j, E, P, t = sp.symbols("h j E P t", positive=True)

    def jets(e):
        e = e.subs({sp.Derivative(Hf, p): h, sp.Derivative(Jf, q): j})
        e = e.replace(lambda x: x.func == sp.exp, lambda x: E ** simp(x.args[0] / (2 * Hf * Jf)))
        e = simp(e)
        assert not e.has(Hf) and not e.has(Jf), e
        return e
    eqs, names = [], []
    for (a, b) in ((0, 0), (0, 1), (1, 1), (2, 2), (3, 3)):
        v = jets(Res[a, b])
        if v != 0:
            eqs.append(sp.numer(sp.together(v)))
            names.append("pqyz"[a] + "pqyz"[b])
    unknown = [A1, A2, Vh]
    M, _ = sp.linear_eq_to_matrix(eqs, unknown)
    Mp = M.subs(PI, P)
    minors = list({sp.expand(Mp[list(r), :].det()) for r in itertools.combinations(range(Mp.shape[0]), 3)} - {0})
    gb = sp.groebner(minors + [t * h * j * E * c * P - 1], t, h, j, E, c, P, order="grevlex")
    ok_i = list(gb.exprs) == [1]
    g3 = sp.factor(sp.gcd_list(minors)) if minors else 0
    sol = sp.linsolve(eqs, unknown)
    want = (1 / (32 * PI), 1 / (32 * PI), 0)
    pt = list(sol)[0] if len(sol) == 1 else None
    ok_ii = pt is not None and all(simp(a - b) == 0 for a, b in zip(pt, want)) and \
        all(simp(e.subs(dict(zip(unknown, want)))) == 0 for e in eqs)
    ok_iii = pt is not None and not any(sp.sympify(a).has(c) for a in pt)
    check("E_solve", ok_i and ok_ii and ok_iii,
          f"unknowns alpha_1^2, alpha_2^2, V at a point; nonzero components {names}; (i) Groebner basis of the 3x3 minors with "
          f"t h j E c P - 1 is {list(gb.exprs)} (rank 3 wherever h j E c != 0), gcd of the minors {g3}; "
          f"(ii) linsolve gives {pt if pt is not None else sol}, expected (1/(32 pi), 1/(32 pi), 0), substituted: "
          f"{'0' if ok_ii else 'not 0'}; (iii) c-free {ok_iii} ({time.time() - t0:.0f} s)")
    for nm, e in zip(names, eqs):  # information, not a check: the residual's components, = 0
        print(f"    component {nm}: {sp.factor(e)} = 0   (H' -> h, J' -> j, e^(2HJ) -> E, A1 = alpha_1^2, A2 = alpha_2^2, "
              f"Vh = V at the point; numerator)", flush=True)


def check_E_residual():
    t0 = time.time()
    geo = geo_E()
    fl = fields_E(Hf, Jf, c, ALPHA, ALPHA)
    V = 1 / (80 * PI) if MUT == "potential" else 0
    signs = signs_E()
    bad = nonzero(residual(geo, fl, signs, V))
    Gnz = [(i, j) for i in range(4) for j in range(i, 4) if simp(geo["G"][i, j]) != 0]
    T1 = nonzero(stress(geo, [fl[0]], [1], 0))
    ok = not bad and len(Gnz) > 0 and len(T1) > 0
    Kc = simp(geo["G"][2, 2])
    check("E_residual", ok, f"H(p), J(q) undefined, c a real nonzero symbol, alpha_1 = alpha_2 = 1/(4 sqrt(2 pi)), V = {V}, "
          f"signs {signs}: {'all 16 components of G_ab - 8 pi (T_1 - T_2)_ab are 0' if not bad else brief(bad)}; "
          f"G_yy = G_zz = {Kc} (nonzero components of G: {len(Gnz)}), T_1 has {len(T1)} nonzero components "
          f"({time.time() - t0:.0f} s)")
    T = (8 * PI * stress(geo, fl, signs, V)).applyfunc(simp)  # information, not a check
    nzT = [(i, j, v) for i, j, v in nonzero(T) if i <= j]
    dTc = nonzero(T.applyfunc(lambda e: sp.diff(e, c)))
    print(f"    information: components (i <= j, order p, q, y, z) of G: {[(i, j, geo['G'][i, j]) for i, j in Gnz]}; "
          f"of 8 pi (T_1 - T_2): {nzT}; components of d/dc of 8 pi (T_1 - T_2) that are not 0: {len(dTc)}", flush=True)


def check_E_rform():
    t0 = time.time()
    geo = geo_E()
    fl = fields_E(Hf, Jf, c, ALPHA, ALPHA)
    signs = signs_E()
    u = [sp.sqrt(8 * PI) * grad(geo, f) for f in fl]  # a = sqrt(8 pi) d phi_1, b = sqrt(8 pi) d phi_2
    lam = 0  # 8 pi V
    Rf = geo["Ric"] - lam * geo["g"] - sum((s * uu * uu.T for s, uu in zip(signs, u)), sp.zeros(4, 4))
    bad = nonzero(Rf)
    dH = sp.Matrix([sp.diff(Hf, p), 0, 0, 0])
    dJ = sp.Matrix([0, sp.diff(Jf, q), 0, 0])
    a_exp = (c * dH - 4 / c * dJ) / 2
    b_exp = (c * dH + 4 / c * dJ) / 2
    bad_id = nonzero(u[0] - a_exp) + nonzero(u[1] - b_exp)
    check("E_rform", not bad and not bad_id, f"R_ab - lambda g_ab - (a_a a_b - b_a b_b) with a = sqrt(8 pi) d phi_1, b = sqrt(8 pi) "
          f"d phi_2, lambda = 0: {'all 16 components 0' if not bad else brief(bad)}; a = (c dH - (4/c) dJ)/2, "
          f"b = (c dH + (4/c) dJ)/2: {'both 0' if not bad_id else brief(bad_id)} ({time.time() - t0:.0f} s)")


def check_E_explicit():
    t0 = time.time()
    msgs, ok = [], True
    cases = [(p**3 - p, q**2 + q, sp.Rational(3, 2)), (sp.sin(p), q**3 - q, sp.Rational(-2, 5))]
    for H, J, cc in cases:
        geo = geometry(metric_E(H, J), XE)
        fl = fields_E(H, J, cc, ALPHA, ALPHA)
        bad = nonzero(residual(geo, fl, signs_E(), 0))
        ok = ok and not bad
        msgs.append(f"H = {H}, J = {J}, c = {cc}: {'0' if not bad else brief(bad)}")
    sg = (sp.diff(cases[1][0], p) * sp.diff(cases[1][1], q)).subs({p: 0, q: 1}), \
        (sp.diff(cases[1][0], p) * sp.diff(cases[1][1], q)).subs({p: 0, q: 0})
    ok_sg = sg[0] > 0 and sg[1] < 0
    check("E_explicit", ok and ok_sg, f"residual with explicit functions, alpha_1 = alpha_2 = 1/(4 sqrt(2 pi)), V = 0: "
          + "; ".join(msgs) + f"; H'J' of the second at (0, 1), (0, 0) = {sg} ({time.time() - t0:.0f} s)")


CHECKS = [("conv_energy", check_conv_energy), ("planted_p1", lambda: check_planted(1)),
          ("planted_p2", lambda: check_planted(2)), ("planted_flat", check_planted_flat), ("E_solve", check_E_solve),
          ("E_residual", check_E_residual), ("E_rform", check_E_rform), ("E_explicit", check_E_explicit)]


# ---------------------------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------------------------
def run_controls(t_start):
    print("\n== Controls: the full script with --mutate=NAME, one process each ==", flush=True)
    rows, bad = [], False
    for name, what in MUTATIONS.items():
        t0 = time.time()
        r = subprocess.run([sys.executable, __file__, f"--mutate={name}"], capture_output=True, text=True)
        failed = re.findall(r"^\[FAIL\] (\w+):", r.stdout, flags=re.M)
        want = TARGETS.get(name)
        if want is None:
            outcome, extras = "undeclared", failed
        else:
            fired = all(w in failed for w in want)
            extras = [f for f in failed if f not in want]
            outcome = ("fired+extra" if extras else "fired") if fired else "SILENT"
            bad = bad or not fired
        rows.append((name, want, outcome, failed, what))
        print(f"{name}: {what}\n    target {want}; outcome {outcome}; failed checks {failed}; exit {r.returncode} "
              f"({time.time() - t0:.0f} s)", flush=True)
        show = "E_residual" if want and "E_residual" in want else (want[0] if want else None)
        for line in r.stdout.splitlines():
            if show and line.startswith(f"[FAIL] {show}:"):
                print("    " + line[:700], flush=True)
    return rows, bad


def main():
    global MUT
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    MUT = opts.get("mutate")
    if MUT is not None and MUT not in MUTATIONS:
        sys.exit(f"unknown mutation {MUT!r}. Known: {', '.join(MUTATIONS)}")
    only = opts["checks"].split(",") if "checks" in opts else None
    t_start = time.time()
    print(f"planar_residual.py (#171), sympy {sp.__version__}, python {sys.version.split()[0]}. mutation: "
          f"{MUT or 'none'}" + (f" ({MUTATIONS[MUT]})" if MUT else ""), flush=True)
    for name, fn in CHECKS:
        if only is None or name in only:
            fn()
    nfail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{len(RESULTS) - nfail} of {len(RESULTS)} checks pass ({time.time() - t_start:.0f} s)", flush=True)
    code = 1 if nfail else 0
    if MUT is None and only is None and "no-controls" not in opts:
        rows, bad = run_controls(t_start)
        code = 1 if (nfail or bad) else 0
        print(f"\n{sum(1 for r in rows if r[2] in ('fired', 'fired+extra'))} of {len(rows)} controls fire on their targets; "
              f"{sum(1 for r in rows if r[2] == 'fired+extra')} with extras; {sum(1 for r in rows if r[2] == 'SILENT')} silent; "
              f"{sum(1 for r in rows if r[2] == 'undeclared')} undeclared (total {time.time() - t_start:.0f} s)")
    sys.exit(code)


if __name__ == "__main__":
    main()
