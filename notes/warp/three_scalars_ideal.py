#!/usr/bin/env python3
"""Ideal step for three mixed-sign scalars on the constant-speed Alcubierre wall (issue #157).

At lambda = 0 three scalars of mixed signs close on an open box of the wall iff the orbit-space metric h on (t, r, theta)
is flat there (#77, Lemma C). h depends on r and theta through w(r) = v (1 - f(r)). The numerators N_ab = det(h)^2
Ric(h)_ab are polynomials in (r, c, S, W0..W4), c = cos theta, S = sin theta, W_k = w^(k)(r), linear in W4.

    reduction   N has weights r = 1, W_k = -k, so r = 1 loses nothing and u_k = r^k W_k are the jets; the saturation by
                r w' is the saturation by u1. At fixed r, N_ab = A(c) + S B(c) mod S^2 = 1 - c^2 is a trigonometric
                polynomial, so N = 0 on an interval of theta iff every coefficient of c^k and S c^k vanishes: the ideal
                J_all. N at one theta_0 is a combination of those coefficients, so each theta_0 generates a subideal.
    certificate if 1 - s u1 and the generators of a rung generate the unit ideal, Singular's `liftstd` gives cofactors
                c_i over Q with sum c_i f_i = 1. The default mode regenerates the f_i in sympy and checks the sum
                exactly, without Singular.

    python3 three_scalars_ideal.py                 verify: checks, planted cases, stored certificates (no Singular)
    python3 three_scalars_ideal.py --compute       run Singular on every rung, write inputs, outputs, certificates to RAW
    python3 three_scalars_ideal.py --mutate=NAME   one control (its target checks must fail, exit 1)
    --raw=DIR   the raw folder (default ~/Source/negative-mass-raw/I157-rpi5)
    --rungs=R1,R2,R3   which rungs (default all)

Imported read-only from three_scalars.py: `wall_frame`, `orbit_metric`, `D_jets`, `ricci_pt`, `planted_geometry`,
`closure3`, `kernel`, the symbols RW, CC, SS, WJ and the planted-block samplers. The numerators (`ric_num`), the
coefficient ideals, the Singular driver and the certificate check are new.
"""
import hashlib
import json
import math
import os
import random
import re
import subprocess
import sys
import time
from fractions import Fraction

import sympy as sp
from sympy import QQ
from sympy.polys.rings import ring

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import three_scalars as t3  # noqa: E402

MUTATIONS = {
    "num-sign": (["num_cross"], "the Gamma Gamma term of the Ricci numerator gets the opposite sign"),
    "num-det": (["num_cross"], "the det * (A d2Gamma) term loses the factor det"),
    "coef-drop": (["coef_reassemble"], "the highest-c-degree coefficients are dropped from J_all"),
    "cert-corrupt": (["cert_verify"], "1 is added to the constant term of the first nonzero cofactor"),
    "cert-satoff": (["cert_verify"], "the verification drops the generator 1 - s u1"),
    "plant-sat": (["planted_family"], "PL2 without the saturation by det h"),
    "flat-sub": (["wall_R1"], "h replaced by diag(-1, 1, 1)"),
    "rung-equator": (["wall_R3"], "the wall verdict taken from R1 alone"),
}
MUT = None
RESULTS = []
RAW = os.path.expanduser("~/Source/negative-mass-raw/I157-rpi5")
T0 = time.time()


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# Rings. Rn: numerators in (c, S, W0..W4) at r = 1. Ru: the jets only. Rs: the jets and the saturation variable s.
# ---------------------------------------------------------------------------------------------------------------
CC, SS, RW, WJ = t3.CC, t3.SS, t3.RW, t3.WJ
Rn, c, S, W0, W1, W2, W3, W4 = ring([CC, SS] + list(WJ[:5]), QQ)
JV = ["u0", "u1", "u2", "u3", "u4"]
Rs, s_, u0, u1, u2, u3, u4 = ring(["s"] + JV, QQ)  # the ring of the ideals: s first, then the jets
UV = [u0, u1, u2, u3, u4]


def from_sym(e, R):
    """A sympy expression in R's symbols as an element of R."""
    return R.from_expr(sp.expand(sp.sympify(e)))


def redS(p):
    """Reduce S^2 -> 1 - c^2 in an element of Rn."""
    one = Rn.one
    pw = {0: one}
    out = Rn.zero
    for mon, co in p.terms():
        e = mon[1]
        k = e // 2
        if k not in pw:
            pw[k] = (one - c**2) ** k
        m2 = list(mon)
        m2[1] = e % 2
        out += pw[k] * Rn.from_dict({tuple(m2): co})
    return out


# ---------------------------------------------------------------------------------------------------------------
# Numerators of det(g)^2 Ric(g) for a 3 x 3 metric at a point, from the values, first and second derivatives. The same
# function runs the wall (polynomial entries in Rn) and the planted blocks (sympy entries, with `simp`).
# ---------------------------------------------------------------------------------------------------------------
def adj3(M, simp):
    A = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            r_ = [x for x in range(3) if x != j]
            c_ = [x for x in range(3) if x != i]
            m = M[r_[0]][c_[0]] * M[r_[1]][c_[1]] - M[r_[0]][c_[1]] * M[r_[1]][c_[0]]
            A[i][j] = simp(m if (i + j) % 2 == 0 else -m)
    return A


def det3(M, simp):
    return simp(M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def ric_num(H, DH, DDH, zero, half, simp=lambda x: x):
    """Q_bc = det^2 Ric_bc and det, for the metric H (3 x 3) with DH[a][b][c] = d_a H_bc and DDH[e][a][b][c] = d_e d_a
    H_bc, all at one point. With hi = adj / det: Gam^a_bc = GN / det, d_e Gam^a_bc = (T1 + det T2) / det^2, so that
    Ric = d_a Gam^a_bc - d_c Gam^a_ab + Gam Gam - Gam Gam has the common denominator det^2. `zero` and `half` are the
    additive zero and the number 1/2 of the coefficient ring."""
    n = 3
    sm = lambda it: sum(it, zero)
    Dt = det3(H, simp)
    adj = adj3(H, simp)
    Gl = [[[simp((DH[b][d][k] + DH[k][d][b] - DH[d][b][k]) * half) for k in range(n)] for b in range(n)]
          for d in range(n)]
    dGl = [[[[simp((DDH[e][b][d][k] + DDH[e][k][d][b] - DDH[e][d][b][k]) * half) for k in range(n)] for b in range(n)]
            for d in range(n)] for e in range(n)]
    mm = lambda A, B: [[simp(sm(A[i][k] * B[k][j] for k in range(n))) for j in range(n)] for i in range(n)]
    dhi = [[[simp(-x) for x in row] for row in mm(adj, mm(DH[e], adj))] for e in range(n)]  # d_e(hi) det^2
    GN = [[[simp(sm(adj[a][d] * Gl[d][b][k] for d in range(n))) for k in range(n)] for b in range(n)]
          for a in range(n)]
    sgn = -1 if MUT == "num-sign" else 1
    Q = {}
    for b in range(n):
        for k in range(b, n):
            tot = zero
            for a in range(n):
                t1 = sm(dhi[a][a][d] * Gl[d][b][k] for d in range(n))
                t2 = sm(adj[a][d] * dGl[a][d][b][k] for d in range(n))
                tot += t1 + (t2 if MUT == "num-det" else Dt * t2)
                t1 = sm(dhi[k][a][d] * Gl[d][a][b] for d in range(n))
                t2 = sm(adj[a][d] * dGl[k][d][a][b] for d in range(n))
                tot -= t1 + (t2 if MUT == "num-det" else Dt * t2)
                for d in range(n):
                    tot += sgn * (GN[a][a][d] * GN[d][b][k] - GN[a][k][d] * GN[d][a][b])
            Q[(b, k)] = simp(tot)
    return Q, Dt


def wall_numerators():
    """N_ab = det(h)^2 Ric(h)_ab at r = 1 as elements of Rn, reduced by S^2 = 1 - c^2, and det h. h is `orbit_metric`
    of `wall_frame` (three_scalars.py), the derivatives are `D_jets`."""
    Rf = t3.wall_frame(None)
    R3 = Rf[:3, :3]
    h = sp.diag(-1, 1, 1) if MUT == "flat-sub" else t3.orbit_metric(R3, WJ[0])
    D = t3.D_jets
    dl = [lambda E: sp.Integer(0), lambda E: D(E, 1), lambda E: D(E, 2)]
    at1 = lambda e: from_sym(sp.sympify(e).subs(RW, 1), Rn)
    n = 3
    H = [[at1(h[i, j]) for j in range(n)] for i in range(n)]
    DH = [[[at1(dl[a](h[b, k])) for k in range(n)] for b in range(n)] for a in range(n)]
    DDH = [[[[at1(dl[e](dl[a](h[b, k]))) for k in range(n)] for b in range(n)] for a in range(n)] for e in range(n)]
    Q, Dt = ric_num(H, DH, DDH, Rn.zero, QQ(1, 2), redS)
    return {k: redS(v) for k, v in Q.items()}, redS(Dt)


# ---------------------------------------------------------------------------------------------------------------
# Generators of the rungs, in the ring Rs = Q[s, u0..u4]
# ---------------------------------------------------------------------------------------------------------------
def to_Rs(p):
    """An element of Rn without c, S (a polynomial in W0..W4) as an element of Rs with W_k -> u_k, s absent."""
    out = Rs.zero
    for mon, co in p.terms():
        assert mon[0] == 0 and mon[1] == 0, mon
        out += Rs.from_dict({(0,) + tuple(mon[2:]): co})
    return out


def normalize(p):
    """Integer coefficients, content 1, leading coefficient (in the ring's order) positive. Deterministic."""
    if p == 0:
        return p
    den = 1
    for co in p.coeffs():
        den = den * co.denominator // math.gcd(den, co.denominator)
    q = p * den
    g = 0
    for co in q.coeffs():
        g = math.gcd(g, int(co.numerator))
    q = q * QQ(1, g)
    return q if q.LC > 0 else -q


def eval_cS(p, cv, sv):
    """p in Rn at (c, S) = (cv, sv) rational: a polynomial in W0..W4 as an element of Rs."""
    out = Rs.zero
    for mon, co in p.terms():
        out += Rs.from_dict({(0,) + tuple(mon[2:]): co * QQ(cv) ** mon[0] * QQ(sv) ** mon[1]})
    return out


AXIS_POW = {(0, 0): 2, (0, 1): 2, (0, 2): 2, (1, 1): 2, (1, 2): 2, (2, 2): 3}


def strip_axis(Q):
    """N_ab / (1 - c^2)^k_ab, exactly (k = 2, and 3 for N_22): the numerators vanish on the axis to that order. Off the
    axis, where the box lies, 1 - c^2 = S^2 != 0. Asserts that the division is exact."""
    out = {}
    for key, p in Q.items():
        den = (1 - c**2) ** AXIS_POW[key]
        qa = Rn.zero
        for e in (0, 1):
            part = Rn.from_dict({(m[0], 0) + tuple(m[2:]): co for m, co in p.terms() if m[1] == e})
            q_, r_ = part.div(den)
            assert r_ == 0, (key, e)
            qa += q_ * (S**e)
        out[key] = qa
    return out


def coef_gens(Q, drop=False):
    """Every coefficient of c^k and S c^k of every N_ab, as (label, element of Rs). With drop the highest-c-degree
    coefficient of each N_ab is omitted (the control coef-drop)."""
    out = []
    for key in sorted(Q):
        p = Q[key]
        by = {}
        for mon, co in p.terms():
            by.setdefault((mon[0], mon[1]), {})[(0,) + tuple(mon[2:])] = co
        kmax = max((k for (k, _) in by), default=0)
        for (k, e), d in sorted(by.items()):
            if drop and k == kmax:
                continue
            out.append((f"N{key[0]}{key[1]}_c{k}_S{e}", Rs.from_dict(d)))
    return out


def rung_gens(Q, rung):
    """The generators of a rung, normalized, without duplicates: list of (label, poly in Rs). Q holds the stripped
    numerators N_ab / (1 - c^2)^k (`strip_axis`)."""
    if rung == "R1":
        raw = [(f"N{k[0]}{k[1]}@eq", eval_cS(v, 0, 1)) for k, v in sorted(Q.items())]
    elif rung == "R2":
        raw = [(f"N{k[0]}{k[1]}@eq", eval_cS(v, 0, 1)) for k, v in sorted(Q.items())]
        raw += [(f"N{k[0]}{k[1]}@35", eval_cS(v, QQ(3, 5), QQ(4, 5))) for k, v in sorted(Q.items())]
    elif rung == "R3":
        raw = coef_gens(Q, drop=(MUT == "coef-drop"))
    else:
        raise ValueError(rung)
    out, seen = [], set()
    for lab, p in raw:
        if p == 0:
            continue
        q = normalize(p)
        key = tuple(sorted(q.terms()))
        if key in seen:
            continue
        seen.add(key)
        out.append((lab, q))
    return out


SAT = 1 - s_ * u1  # the generator 1 - s u1


# ---------------------------------------------------------------------------------------------------------------
# Singular: input, run, output, parse
# ---------------------------------------------------------------------------------------------------------------
def sing_poly(p):
    """A polynomial of Rs in Singular syntax."""
    if p == 0:
        return "0"
    names = ["s"] + JV
    parts = []
    for mon, co in sorted(p.terms(), reverse=True):
        co = Fraction(int(co.numerator), int(co.denominator))
        vs = "*".join(f"{nm}^{e}" if e > 1 else nm for nm, e in zip(names, mon) if e)
        body = (str(abs(co.numerator)) + ("/" + str(co.denominator) if co.denominator != 1 else "")) if (
            abs(co) != 1 or not vs) else ""
        term = body + ("*" if body and vs else "") + vs
        parts.append(("-" if co < 0 else "+") + term)
    out = "".join(parts)
    return out[1:] if out.startswith("+") else out


def parse_poly(txt):
    """A polynomial printed by Singular (short = 0) as an element of Rs, exactly."""
    txt = txt.replace(" ", "").replace("\n", "")
    if not txt:
        return Rs.zero
    names = ["s"] + JV
    d = {}
    for m in re.finditer(r"([+-]?)([^+-]+)", txt):
        sign, body = m.group(1), m.group(2)
        co, mon = Fraction(1), [0] * 6
        for fac in body.split("*"):
            if re.fullmatch(r"\d+(/\d+)?", fac):
                co *= Fraction(fac)
            else:
                mm = re.fullmatch(r"([a-z]\d?)(\^(\d+))?", fac)
                assert mm and mm.group(1) in names, fac
                mon[names.index(mm.group(1))] += int(mm.group(3) or 1)
        if sign == "-":
            co = -co
        d[tuple(mon)] = d.get(tuple(mon), 0) + QQ(co.numerator, co.denominator)
    return Rs.from_dict({k: v for k, v in d.items() if v != 0})


def sing_header(char=0):
    return f"short = 0;\nring R = {char},(s," + ",".join(JV) + "),dp;\noption(redSB);\n"


def run_singular(case, gens, program, timeout=3600, char=0):
    """Write RAW/singular/<case>.sing, run it under /usr/bin/time -v at nice 10, keep .out/.err/.time. Returns the stdout
    (None on timeout) and the wall time."""
    d = os.path.join(RAW, "singular")
    os.makedirs(d, exist_ok=True)
    src = sing_header(char) + "ideal I = " + ",\n  ".join(sing_poly(p) for _, p in gens) + ";\n" + program + "\nquit;\n"
    base = os.path.join(d, case)
    open(base + ".sing", "w").write(src)
    t = time.time()
    try:
        p = subprocess.run(["nice", "-n", "10", "/usr/bin/time", "-v", "-o", base + ".time", "/usr/bin/Singular", "-q",
                            base + ".sing"], capture_output=True, text=True, timeout=timeout,
                           stdin=subprocess.DEVNULL)
    except subprocess.TimeoutExpired as e:
        open(base + ".out", "w").write((e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or ""))
        open(base + ".err", "w").write("TIMEOUT after %d s\n" % timeout)
        return None, time.time() - t
    open(base + ".out", "w").write(p.stdout)
    open(base + ".err", "w").write(p.stderr)
    return p.stdout, time.time() - t


# ---------------------------------------------------------------------------------------------------------------
# Certificates
# ---------------------------------------------------------------------------------------------------------------
def gens_text(gens):
    return "\n".join(f"{lab}: {sing_poly(p)}" for lab, p in gens)


def unit_program():
    """liftstd over Q: G = {1} and the cofactors T[i, 1] with sum I[i] T[i, 1] = 1, or `NOTUNIT`."""
    return """
int t0 = rtimer;
matrix T;
ideal G = liftstd(I, T);
print("LIFTSTD_SECONDS " + string(rtimer - t0));
if (G[1] != 1) { print("NOTUNIT size " + string(size(G)) + " dim " + string(dim(G))); }
else {
  print("UNIT");
  int i;
  for (i = 1; i <= ncols(I); i++) { print("COF " + string(i) + " " + string(T[i, 1])); }
}
"""


def compute_rung(rung, Q, tag=""):
    """Run Singular on rung `rung` plus the saturation generator, over Q. Writes RAW/certs/<rung><tag>.json when the
    ideal is the unit ideal. Returns (verdict, seconds)."""
    gens = rung_gens(Q, rung) + [("sat", SAT)]
    case = f"{rung}{tag}_unit"
    out, wall = run_singular(case, gens, unit_program(), timeout=3600)
    if out is None:
        print(f"{rung}: Singular did not finish in 1 h; stopped ({wall:.0f} s)", flush=True)
        return "timeout", wall
    if "NOTUNIT" in out:
        print(f"{rung}: not the unit ideal (Singular, over Q): {[l for l in out.splitlines() if 'NOTUNIT' in l]}")
        return "notunit", wall
    cof = {}
    for ln in out.splitlines():
        m = re.match(r"COF (\d+) (.*)$", ln)
        if m:
            cof[int(m.group(1))] = m.group(2)
    assert "UNIT" in out and len(cof) == len(gens), "no cofactors in Singular's output"
    cert = {"rung": rung, "ring": ["s"] + JV, "labels": [lab for lab, _ in gens], "n_gens": len(gens),
            "gen_sha256": hashlib.sha256(gens_text(gens).encode()).hexdigest(), "singular": "4.4.1",
            "singular_seconds": wall, "cofactors": [cof[i + 1] for i in range(len(gens))]}
    os.makedirs(os.path.join(RAW, "certs"), exist_ok=True)
    path = os.path.join(RAW, "certs", f"{rung}{tag}.json")
    json.dump(cert, open(path, "w"))
    print(f"{rung}: unit ideal, certificate written to {path} ({os.path.getsize(path)} bytes, {wall:.0f} s)", flush=True)
    return "unit", wall


def verify_cert(rung, Q, tag="", mut=None):
    """Regenerate the generators in sympy, parse the cofactors exactly and check sum c_i f_i = 1 in Q[s, u0..u4].
    Returns (ok, detail)."""
    path = os.path.join(RAW, "certs", f"{rung}{tag}.json")
    if not os.path.exists(path):
        return False, f"no certificate at {path}"
    cert = json.load(open(path))
    gens = rung_gens(Q, rung) + [("sat", SAT)]
    if hashlib.sha256(gens_text(gens).encode()).hexdigest() != cert["gen_sha256"]:
        return False, "the regenerated generators differ from the ones Singular was given"
    t = time.time()
    cof = [parse_poly(x) for x in cert["cofactors"]]
    nterms = sum(len(x) for x in cof)
    if mut == "cert-corrupt":
        k = next(i for i, x in enumerate(cof) if x != 0)
        cof[k] = cof[k] + 1
    use = list(range(len(gens)))
    if mut == "cert-satoff":
        use = use[:-1]
    tot = Rs.zero
    for i in use:
        tot += cof[i] * gens[i][1]
    ok = tot == 1
    return ok, (f"sum c_i f_i = 1 over {len(use)} of {len(gens)} generators, {nterms} cofactor terms, parsed and "
                f"summed in {time.time() - t:.0f} s" if ok else f"sum c_i f_i - 1 has {len(tot - 1)} terms")


def main_compute(opts):
    Q = load_numerators()
    for rung in opts.get("rungs", "R2").split(","):
        compute_rung(rung, Q)


def load_numerators():
    t = time.time()
    Q, _ = wall_numerators()
    Q = strip_axis(Q)
    print(f"numerators N_ab / (1 - c^2)^k: terms {[len(Q[k]) for k in sorted(Q)]} ({time.time() - t:.0f} s)", flush=True)
    return Q


def parse_opts(argv):
    return dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in argv if a.startswith("--"))


if __name__ == "__main__":
    opts = parse_opts(sys.argv[1:])
    MUT = opts.get("mutate")
    RAW = os.path.expanduser(opts.get("raw", RAW))
    if "compute" in opts:
        main_compute(opts)
