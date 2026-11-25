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


def peak_kb(case):
    """Peak resident size of Singular's run, from `/usr/bin/time -v`."""
    try:
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", open(os.path.join(RAW, "singular", case + ".time")).read())
        return int(m.group(1)) if m else None
    except OSError:
        return None


def store_verdict(name, rec):
    path = os.path.join(RAW, "verdicts.json")
    d = json.load(open(path)) if os.path.exists(path) else {}
    d[name] = rec
    json.dump(d, open(path, "w"), indent=1, sort_keys=True)


def load_verdicts():
    path = os.path.join(RAW, "verdicts.json")
    return json.load(open(path)) if os.path.exists(path) else {}


def compute_case(name, gens, tag=""):
    """Run Singular's liftstd over Q on gens (the last one is the saturation generator). Unit ideal: write the certificate
    RAW/certs/<name><tag>.json. Records the verdict in RAW/verdicts.json. Returns the verdict."""
    case = f"{name}{tag}_unit"
    out, wall = run_singular(case, gens, unit_program(), timeout=3600)
    if out is None:
        print(f"{name}: Singular did not finish in 1 h; stopped ({wall:.0f} s)", flush=True)
        store_verdict(name, {"verdict": "timeout", "seconds": wall, "peak_kb": peak_kb(case)})
        return "timeout"
    if "NOTUNIT" in out:
        line = [ln for ln in out.splitlines() if "NOTUNIT" in ln][0]
        print(f"{name}: not the unit ideal (Singular, over Q): {line} ({wall:.0f} s)", flush=True)
        store_verdict(name, {"verdict": "notunit", "seconds": wall, "peak_kb": peak_kb(case), "singular": line})
        return "notunit"
    cof = {}
    for ln in out.splitlines():
        m = re.match(r"COF (\d+) (.*)$", ln)
        if m:
            cof[int(m.group(1))] = m.group(2)
    assert "UNIT" in out and len(cof) == len(gens), "no cofactors in Singular's output"
    cert = {"rung": name, "ring": ["s"] + JV, "labels": [lab for lab, _ in gens], "n_gens": len(gens),
            "gen_sha256": hashlib.sha256(gens_text(gens).encode()).hexdigest(), "singular": "4.4.1",
            "singular_seconds": wall, "cofactors": [cof[i + 1] for i in range(len(gens))]}
    os.makedirs(os.path.join(RAW, "certs"), exist_ok=True)
    path = os.path.join(RAW, "certs", f"{name}{tag}.json")
    json.dump(cert, open(path, "w"))
    print(f"{name}: unit ideal, certificate {path} ({os.path.getsize(path)} bytes, Singular {wall:.0f} s)", flush=True)
    store_verdict(name, {"verdict": "unit", "seconds": wall, "peak_kb": peak_kb(case),
                         "cert_bytes": os.path.getsize(path), "cofactor_terms": sum(len(parse_poly(x)) for x in cert["cofactors"])})
    return "unit"


def compute_rung(rung, Q):
    return compute_case(rung, rung_gens(Q, rung) + [("sat", SAT)])


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


# ---------------------------------------------------------------------------------------------------------------
# Independent route to h, from the coordinate Ricci tensor of the 4D metric with a generic function w(r). It shares
# nothing with wall_frame (the frame Ricci tensor of phantom_source.py) or with ric_num: its own Christoffel code, direct
# matrix inverse, and its own 3D Ricci routine on numbers.
# ---------------------------------------------------------------------------------------------------------------
_IND = {}


def ric_coord(g, X, simp=sp.expand):
    """Ricci tensor of a metric g(X) by the textbook formula with the inverse matrix (sympy `inv`)."""
    n = len(X)
    gi = g.inv().applyfunc(sp.cancel)
    Gam = [[[simp(sum(gi[a, d] * (sp.diff(g[d, b], X[k]) + sp.diff(g[d, k], X[b]) - sp.diff(g[b, k], X[d]))
                      for d in range(n)) / 2) for k in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.zeros(n, n)
    for b in range(n):
        for k in range(b, n):
            e = sum(sp.diff(Gam[a][b][k], X[a]) - sp.diff(Gam[a][b][a], X[k]) for a in range(n))
            e += sum(Gam[a][a][d] * Gam[d][b][k] - Gam[a][k][d] * Gam[d][a][b] for a in range(n) for d in range(n))
            Ric[b, k] = Ric[k, b] = e
    return Ric


def indep_block():
    """h_{mu nu} = R_{mu nu} on (t, r, theta) of g = -dt^2 + (w dt + cos(th) dr - r sin(th) dth)^2 + (sin(th) dr + r cos(th)
    dth)^2 + r^2 sin(th)^2 dph^2, w = w(r) a sympy Function, with the generic Ricci tensor (lambda = 0, so M = R)."""
    if "h" in _IND:
        return _IND["h"]
    t, r, th, ph = sp.symbols("t r th ph", real=True)
    w = sp.Function("w")(r)
    e1 = [w, sp.cos(th), -r * sp.sin(th), 0]
    e2 = [0, sp.sin(th), r * sp.cos(th), 0]
    e0 = [1, 0, 0, 0]
    e3 = [0, 0, 0, r * sp.sin(th)]
    # coframe components over (t, r, th, ph): the metric is -e0 e0 + e1 e1 + e2 e2 + e3 e3
    g = sp.zeros(4, 4)
    for i in range(4):
        for j in range(4):
            g[i, j] = -e0[i] * e0[j] + e1[i] * e1[j] + e2[i] * e2[j] + e3[i] * e3[j]
    g = g.applyfunc(sp.simplify)
    X = [t, r, th, ph]
    Ric = ric_coord(g, X)
    _IND["h"] = (Ric[:3, :3], (r, th, w), g)
    return _IND["h"]


def indep_numerators_at(c0, s0, jets):
    """det(h)^2 Ric(h) at r = 1, (c, S) = (c0, s0), w^(k)(1) = jets[k], k = 0..6, by the independent route, as a dict
    {(a, b): Rational}."""
    h, (r, th, w), _ = indep_block()
    pt = {r: 1}
    d = {sp.Derivative(w, (r, k)): jets[k] for k in range(6, 0, -1)}

    def val(e):
        e = sp.sympify(e)
        e = e.subs(d).subs(w, jets[0]).subs({sp.cos(th): c0, sp.sin(th): s0}).subs(pt)
        e = sp.nsimplify(e) if False else sp.cancel(e)
        assert e.is_Rational, e
        return e
    n = 3
    X = [None, r, th]  # coordinate 0 (t) is cyclic
    def D(e, a):
        return sp.Integer(0) if a == 0 else sp.diff(e, X[a])
    g = sp.Matrix(n, n, lambda i, j: val(h[i, j]))
    dg = [[[val(D(h[i, j], a)) for j in range(n)] for i in range(n)] for a in range(n)]
    ddg = [[[[val(D(D(h[i, j], a), e_)) for j in range(n)] for i in range(n)] for a in range(n)] for e_ in range(n)]
    gi = g.inv()
    Gam = [[[sum(gi[a, dd] * (dg[b][dd][k] + dg[k][dd][b] - dg[dd][b][k]) for dd in range(n)) / 2 for k in range(n)]
            for b in range(n)] for a in range(n)]
    dGam = [[[[sum(-sum(gi[a, p_] * dg[e_][p_][q_] * gi[q_, dd] for p_ in range(n) for q_ in range(n)) *
                   (dg[b][dd][k] + dg[k][dd][b] - dg[dd][b][k]) / 2
                   + gi[a, dd] * (ddg[e_][b][dd][k] + ddg[e_][k][dd][b] - ddg[e_][dd][b][k]) / 2 for dd in range(n))
               for k in range(n)] for b in range(n)] for a in range(n)] for e_ in range(n)]
    out = {}
    det = g.det()
    for b in range(n):
        for k in range(b, n):
            ric = sum(dGam[a][a][b][k] - dGam[k][a][a][b] for a in range(n))
            ric += sum(Gam[a][a][dd] * Gam[dd][b][k] - Gam[a][k][dd] * Gam[dd][a][b] for a in range(n) for dd in range(n))
            out[(b, k)] = sp.cancel(det**2 * ric)
    return out


PYTH = [(sp.Rational(3, 5), sp.Rational(4, 5)), (sp.Rational(5, 13), sp.Rational(12, 13)),
        (sp.Rational(-7, 25), sp.Rational(24, 25)), (sp.Rational(20, 29), sp.Rational(21, 29)),
        (sp.Rational(-9, 41), sp.Rational(40, 41)), (sp.Rational(8, 17), sp.Rational(15, 17))]


def at_point(p, cv, sv, jets):
    """The value of an element of Rn at (c, S) = (cv, sv) and W_k = jets[k], as a sympy Rational."""
    tot = sp.Integer(0)
    for mon, co in p.terms():
        v = sp.Rational(int(co.numerator), int(co.denominator)) * cv ** mon[0] * sv ** mon[1]
        for k in range(5):
            v *= jets[k] ** mon[2 + k]
        tot += v
    return tot


def check_numerators(Qfull):
    """num_cross (ricci_pt, route E of three_scalars.py, at 6 exact points at r = 1 and 3 at r != 1), num_indep (the
    independent coordinate route, 3 points), num_weights (weights and the structure used by the coefficient step)."""
    random.seed(157)
    hsym = t3.orbit_metric(Rf_R3(), WJ[0])
    rr = lambda: sp.Rational(random.randint(-9, 9), random.randint(1, 7))
    bad, rows = [], []
    for i in range(6):
        cv, sv = PYTH[i]
        jets = [rr() for _ in range(7)]
        pt = {RW: 1, CC: cv, SS: sv}
        pt.update({WJ[k]: jets[k] for k in range(7)})
        ric = t3.ricci_pt(hsym, t3.D_jets, pt)
        det = hsym.subs(pt).det()
        for (a, b), q in Qfull.items():
            if at_point(q, cv, sv, jets) != sp.cancel(det**2 * ric[a, b]):
                bad.append((i, a, b))
    check("num_cross", not bad, f"det^2 Ric(h)_ab from `ric_num` equals the one from three_scalars.ricci_pt at 6 exact "
          f"points (r = 1, random jets W0..W6, 6 components each): {len(bad)} disagree {bad[:3]}")
    bad = []
    for i in range(3):
        cv, sv = PYTH[i + 1]
        jets = [rr() for _ in range(7)]
        ind = indep_numerators_at(cv, sv, jets)
        for (a, b), q in Qfull.items():
            if at_point(q, cv, sv, jets) != ind[(a, b)]:
                bad.append((i, a, b))
    check("num_indep", not bad, f"the same numerators agree with the independent route (coordinate Ricci tensor of the 4D "
          f"metric with generic w, own Christoffel code and 3D Ricci routine) at 3 exact points: {len(bad)} disagree "
          f"{bad[:3]}")
    # weights: N(r0; W) = r0^m N(1; r0^k W_k) with an integer m per component, at r0 != 1
    ms = {}
    bad = []
    for i in range(3):
        cv, sv = PYTH[i]
        r0 = sp.Rational(random.randint(2, 5), random.randint(2, 5))
        jets = [rr() for _ in range(7)]
        pt = {RW: r0, CC: cv, SS: sv}
        pt.update({WJ[k]: jets[k] for k in range(7)})
        ric = t3.ricci_pt(hsym, t3.D_jets, pt)
        det = hsym.subs(pt).det()
        for (a, b), q in Qfull.items():
            lhs = sp.cancel(det**2 * ric[a, b])
            rhs = at_point(q, cv, sv, [r0**k * jets[k] for k in range(5)])
            if rhs == 0 or lhs == 0:
                bad.append(("zero", i, a, b))
                continue
            m = next((m_ for m_ in range(-40, 41) if lhs == r0**m_ * rhs), None)
            if m is None:
                bad.append(("no integer weight", i, a, b))
                continue
            ms.setdefault((a, b), set()).add(sp.Integer(m))
    ok_m = all(len(v) == 1 and all(x.is_Integer for x in v) for v in ms.values()) and not bad
    struct_bad = []
    for (a, b), q in Qfull.items():
        if max(m[6] for m in q.monoms()) > 1:
            struct_bad.append(("W4 degree", a, b))
        sd = set(m[1] for m in q.monoms())
        if sd - ({1} if (a, b) in ((0, 2), (1, 2)) else {0}):
            struct_bad.append(("S parity", a, b))
    check("num_weights", ok_m and not struct_bad, f"N(r0; W) = r0^m N(1; r0^k W_k) with m = "
          f"{ {k: sorted(v) for k, v in sorted(ms.items())} } at 3 points; linear in W4 and of the S-parity stated: "
          f"{not struct_bad} {struct_bad[:2]}")


_R3 = {}


def Rf_R3():
    if "R3" not in _R3:
        _R3["R3"] = t3.wall_frame(None)[:3, :3]
    return _R3["R3"]


# ---------------------------------------------------------------------------------------------------------------
# Planted cases
# ---------------------------------------------------------------------------------------------------------------
def planted_blocks():
    """PL1, PL1b: the quotient blocks of #77's P1, P2 and D1 through `ric_num`. The block is M = R - lam g extracted on
    the coordinates other than the slice index j of `kernel`, as `closure3` does."""
    simp = lambda e: sp.cancel(sp.together(e))
    res = {}
    for tag, kind, eps in (("P1", "P1", 0), ("P2", "P2", 0), ("D1", "P1", sp.Rational(1, 10))):
        a, sg, lam0, kap = t3.warp(kind)
        g, Ric, G = t3.planted_geometry(a, eps)
        lam_expr = sp.Integer(lam0)
        M = (Ric - lam_expr * g).applyfunc(t3.simp_r)
        K, j = t3.kernel(M, t3.simp_r)
        idx = [i for i in range(4) if i != j]
        h = M.extract(idx, idx)
        Dsl = lambda E, a_: t3.Dnew(E, idx[a_])
        n = 3
        H = [[h[i, k] for k in range(n)] for i in range(n)]
        DH = [[[simp(Dsl(h[b, k], a_)) for k in range(n)] for b in range(n)] for a_ in range(n)]
        DDH = [[[[simp(Dsl(Dsl(h[b, k], a_), e)) for k in range(n)] for b in range(n)] for a_ in range(n)]
               for e in range(n)]
        Q, Dt = ric_num(H, DH, DDH, sp.Integer(0), sp.Rational(1, 2), simp)
        res[tag] = ([k for k, v in sorted(Q.items()) if simp(v) != 0], simp(Dt) != 0, K, j)
    ok = not res["P1"][0] and not res["P2"][0] and res["D1"][0] and all(res[t][1] for t in res)
    check("planted_blocks", ok, f"numerators det^2 Ric of the quotient blocks (slice index {res['P1'][3]}, "
          f"K = {list(res['P1'][2])}): P1 nonzero components {res['P1'][0]}, P2 {res['P2'][0]}, D1 {res['D1'][0]}; "
          f"all three blocks nondegenerate: {all(res[t][1] for t in res)}")


def strip_factor(p, f):
    """p divided by the largest power of f that divides it (exactly), and the power."""
    k = 0
    while True:
        q_, r_ = p.div(f)
        if r_ != 0:
            return p, k
        p, k = q_, k + 1


def planted_family_gens():
    """PL2: h = diag(-1, 1, G^2), G = r (1 + w). Returns (numerators in Rn at r = 1, det h, flat condition 2 W1 + W2,
    the generators after saturating by det h)."""
    h = sp.diag(-1, 1, RW**2 * (1 + WJ[0])**2)
    D = t3.D_jets
    dl = [lambda E: sp.Integer(0), lambda E: D(E, 1), lambda E: D(E, 2)]
    at1 = lambda e: from_sym(sp.sympify(e).subs(RW, 1), Rn)
    n = 3
    H = [[at1(h[i, j]) for j in range(n)] for i in range(n)]
    DH = [[[at1(dl[a](h[b, k])) for k in range(n)] for b in range(n)] for a in range(n)]
    DDH = [[[[at1(dl[e](dl[a](h[b, k]))) for k in range(n)] for b in range(n)] for a in range(n)] for e in range(n)]
    Q, Dt = ric_num(H, DH, DDH, Rn.zero, QQ(1, 2))
    flat = 2 * W1 + W2
    gens = [(f"N{k[0]}{k[1]}", v) for k, v in sorted(Q.items()) if v != 0]
    if MUT != "plant-sat":  # saturation by det h = -G^2: divide out powers of 1 + W0
        gens = [(lab, strip_factor(v, 1 + W0)[0]) for lab, v in gens]
    return Q, Dt, flat, gens


def planted_family():
    Q, Dt, flat, gens = planted_family_gens()
    G3 = (1 + W0) ** 3
    pred = {(1, 1): -G3 * flat, (2, 2): -G3 * (1 + W0) ** 2 * flat}
    num_ok = all((Q[k] == pred.get(k, Rn.zero)) for k in Q)
    # after the saturation every generator is (a unit times) 2 W1 + W2
    sat_ok = all(strip_factor(v, flat)[0] in (Rn.one, -Rn.one) for _, v in gens) and len(gens) >= 1
    # exact points: flat family (u0 free, u2 = -2 u1) and the spurious u0 = -1, where det h = 0
    pt_flat = {W0: sp.Rational(2, 3), W1: sp.Rational(1), W2: sp.Rational(-2)}
    spurious = {W0: sp.Rational(-1), W1: sp.Rational(1), W2: sp.Rational(5, 7)}
    ev = lambda p, pt: at_point(p, 0, 1, [pt[W0], pt[W1], pt[W2], 0, 0])
    flat_ok = all(ev(v, pt_flat) == 0 for v in Q.values())
    spur_ok = all(ev(v, spurious) == 0 for v in Q.values()) and ev(Dt, spurious) == 0 and ev(flat, spurious) != 0
    check("planted_family", num_ok and sat_ok and flat_ok and spur_ok,
          f"h = diag(-1, 1, G^2), G = r (1 + w): N_rr = -G^3 G'', N_thth = -G^5 G'' with G'' = 2 W1 + W2: {num_ok}; "
          f"generators after the saturation by det h are +-(2 W1 + W2): {sat_ok}; the flat point u = (2/3, 1, -2) "
          f"is a root: {flat_ok}; the spurious point u0 = -1 is a root of N with det h = 0 and 2 W1 + W2 != 0: {spur_ok}")


def planted_root_gens(Q):
    """PL3: the generators of R1 with constants shifted so that a random rational jet point j* (u1* != 0) is a root."""
    random.seed(1571)
    jstar = [sp.Rational(random.randint(-9, 9), random.randint(1, 5)) for _ in range(5)]
    jstar[1] = sp.Rational(random.randint(1, 9), random.randint(1, 5))
    base = rung_gens(Q, "R1")
    out = []
    for lab, p in base:
        v = sum((sp.Rational(int(co.numerator), int(co.denominator)) * sp.prod([jstar[k] ** mon[1 + k] for k in range(5)])
                 for mon, co in p.terms()), sp.Integer(0))
        q = p - Rs(QQ(int(v.p), int(v.q)))
        out.append((lab + "_shift", normalize(q)))
    return jstar, out + [("sat", SAT)]


def eval_Rs(p, sv, jets):
    """p in Rs at s = sv and u_k = jets[k], as a sympy Rational."""
    tot = sp.Integer(0)
    for mon, co in p.terms():
        v = sp.Rational(int(co.numerator), int(co.denominator)) * sv ** mon[0]
        for k in range(5):
            v *= jets[k] ** mon[1 + k]
        tot += v
    return tot


def planted_root(Q, verdicts):
    jstar, gens = planted_root_gens(Q)
    ok_pt = all(eval_Rs(p, 1 / jstar[1], jstar) == 0 for _, p in gens)
    v = verdicts.get("PL3", {}).get("verdict")
    check("planted_root", ok_pt and v == "notunit",
          f"R1's generators shifted to have the exact root u* = {jstar} (s = 1/u1*): every generator and 1 - s u1 vanish "
          f"there: {ok_pt}; Singular's verdict on the shifted system: {v} (stored)")


def check_coefficients(Q):
    """coef_reassemble: the coefficients that generate J_all reassemble each N'_ab at 7 exact points."""
    gens = coef_gens(Q, drop=(MUT == "coef-drop"))
    bad = []
    for (cv, sv) in PYTH + [(sp.Integer(0), sp.Integer(1))]:
        for key, v in Q.items():
            tot = Rs.zero
            for lab, p in gens:
                m = re.match(rf"N{key[0]}{key[1]}_c(\d+)_S(\d)$", lab)
                if m:
                    tot += p * QQ(cv.p, cv.q) ** int(m.group(1)) * QQ(sv.p, sv.q) ** int(m.group(2))
            if tot != eval_cS(v, QQ(cv.p, cv.q), QQ(sv.p, sv.q)):
                bad.append((cv, key))
    check("coef_reassemble", not bad, f"the {len(gens)} coefficients of c^k and S c^k reassemble N'_ab at 7 exact points "
          f"(c, S): {len(bad)} mismatches {bad[:2]}")


def main_compute(opts):
    Q = load_numerators()
    todo = opts.get("rungs", "R2").split(",")
    for rung in todo:
        if rung in ("R1", "R2", "R3"):
            compute_rung(rung, Q)
        elif rung == "PL2":
            _, _, _, gens = planted_family_gens()
            compute_case("PL2", [(lab, to_Rs(v)) for lab, v in gens] + [("sat", SAT)])
        elif rung == "PL3":
            compute_case("PL3", planted_root_gens(Q)[1])
        else:
            sys.exit(f"unknown rung {rung}")


def load_numerators():
    t = time.time()
    Q, _ = wall_numerators()
    Q = strip_axis(Q)
    print(f"numerators N_ab / (1 - c^2)^k: terms {[len(Q[k]) for k in sorted(Q)]} ({time.time() - t:.0f} s)", flush=True)
    return Q


def parse_opts(argv):
    return dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in argv if a.startswith("--"))


def find_point_R1(Q):
    """An exact rational point of V(R1 + sat) with u1 != 0, so R1 is not the unit ideal. At the equator N_00, N_02, N_22
    carry the factor Phi' = 4 u0^2 u1^2 - u1^4 + u1^2 + 2 u1 u2 + u2^2 (det h = 0 there) and N_11 is linear in u4: take
    u1 = 1, (1 + u2, 2 u0) a rational point of the unit circle, u3 in a short list, and solve N_11 for u4."""
    gens = rung_gens(Q, "R1")
    u4s = sp.Symbol("u4")
    syms = sp.symbols("s u0 u1 u2 u3 u4")
    for (X, Y) in PYTH + [(-x, y) for x, y in PYTH]:
        for u3v in (0, 1, -1, 2):
            pt = [Y / 2, sp.Integer(1), X - 1, sp.Integer(u3v)]
            sol = None
            for _, g in gens:
                e = sp.sympify(g.as_expr()).subs(dict(zip(syms[1:5], pt)))
                P = sp.Poly(e, syms[5])
                if P.degree() == 1:
                    sol = -P.coeff_monomial(1) / P.coeff_monomial(syms[5])
                    break
            if sol is None or not sp.sympify(sol).is_Rational:
                continue
            jets = pt + [sp.Rational(sol)]
            if all(eval_Rs(g, 1, jets) == 0 for _, g in gens + [("sat", SAT)]):
                return jets
    return None


def wall_rung(rung, Q, verdicts):
    v = verdicts.get(rung, {})
    secs = f"{v.get('seconds', float('nan')):.0f} s, peak {(v.get('peak_kb') or 0) / 1024:.0f} MB"
    if rung == "R1":
        pt = find_point_R1(Q)
        check("wall_R1", pt is not None and v.get("verdict") == "notunit",
              f"theta_0 = pi/2, saturated by u1: not the unit ideal (Singular {secs}); exact point u = {pt} lies in "
              f"V(R1 + 1 - s u1) with s = 1/u1, on det h = 0 (Phi' = 0)")
        return
    if rung not in verdicts:
        print(f"[SKIP] {rung}: not computed (no verdict in {RAW})")
        return
    ok, detail = verify_cert(rung, Q, mut=MUT)
    check(f"cert_verify_{rung}", ok, detail)
    vv = verdicts["R1"] if (MUT == "rung-equator" and rung == "R3") else v
    check(f"wall_{rung}", ok and vv.get("verdict") == "unit",
          f"{rung} with 1 - s u1: the unit ideal, certified (Singular {secs}): {ok}; stored verdict {vv.get('verdict')}")


STOP_AT = {"num-sign": 1, "num-det": 1, "plant-sat": 2, "coef-drop": 3}


def main(opts):
    global T0
    T0 = time.time()
    random.seed(157)
    sec_stop = STOP_AT.get(MUT, 99)
    print(f"three_scalars_ideal.py, issue #157. sympy {sp.__version__}, python {sys.version.split()[0]}, mutation: "
          f"{MUT or 'none'}" + (f" ({MUTATIONS[MUT][1]})" if MUT else "") + f". raw: {RAW}")
    verdicts = load_verdicts()
    print("\n== 1. The numerators det(h)^2 Ric(h) of the wall, at r = 1, and their checks ==")
    Qfull, Dt = wall_numerators()
    Q = strip_axis(Qfull)
    print(f"  N_ab terms {[len(Qfull[k]) for k in sorted(Qfull)]}; after dividing by (1 - c^2)^k (k = 2, 2, 2, 2, 2, 3): "
          f"{[len(Q[k]) for k in sorted(Q)]}; det h terms {len(Dt)} ({time.time() - T0:.0f} s)")
    check_numerators(Qfull)
    print(f"  ({time.time() - T0:.0f} s)")
    if sec_stop < 2:
        return fin()
    print("\n== 2. Planted truth: flat blocks, a family with a spurious det h = 0 component, a planted root ==")
    planted_blocks()
    planted_family()
    planted_root(Q, verdicts)
    print(f"  ({time.time() - T0:.0f} s)")
    if sec_stop < 3:
        return fin()
    print("\n== 3. The coefficient ideal J_all ==")
    check_coefficients(Q)
    if sec_stop < 4:
        return fin()
    print("\n== 4. The ladder: saturation by u1 = r w', verdict per rung ==")
    for rung in opts.get("rungs", "R1,R2,R3").split(","):
        wall_rung(rung, Q, verdicts)
        print(f"  ({time.time() - T0:.0f} s)", flush=True)
    return fin()


def fin():
    failed = [n for n, ok in RESULTS if not ok]
    print(f"\nchecks: {len(RESULTS) - len(failed)} passed, {len(failed)} failed {failed}. {time.time() - T0:.0f} s")
    if MUT is not None:
        want = MUTATIONS[MUT][0]
        print(f"control {MUT}: target checks {want}, fired {[n for n in want if n in failed]}")
    return failed


def ingest(name):
    """Rebuild the verdict and certificate of a Singular run from its stored output (no Singular run)."""
    base = os.path.join(RAW, "singular", f"{name}_unit")
    out = open(base + ".out").read()
    tm = open(base + ".time").read()
    el = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.*)", tm).group(1)
    parts = [float(x) for x in el.split(":")]
    wall = sum(x * 60**i for i, x in enumerate(reversed(parts)))
    if "NOTUNIT" in out:
        line = [ln for ln in out.splitlines() if "NOTUNIT" in ln][0]
        store_verdict(name, {"verdict": "notunit", "seconds": wall, "peak_kb": peak_kb(f"{name}_unit"), "singular": line})
        print(f"{name}: notunit ({wall:.0f} s)")
        return
    path = os.path.join(RAW, "certs", f"{name}.json")
    cert = json.load(open(path))
    store_verdict(name, {"verdict": "unit", "seconds": wall, "peak_kb": peak_kb(f"{name}_unit"),
                         "cert_bytes": os.path.getsize(path),
                         "cofactor_terms": sum(len(parse_poly(x)) for x in cert["cofactors"])})
    print(f"{name}: unit ({wall:.0f} s)")


def run_controls(opts):
    print("\n== Controls: each mutation must make its target checks fail ==", flush=True)
    ok_all = True
    for name, (want, what) in MUTATIONS.items():
        p = subprocess.run([sys.executable, os.path.abspath(__file__), f"--mutate={name}", f"--raw={RAW}"],
                           capture_output=True, text=True)
        fails = [ln.split(":")[0][7:] for ln in p.stdout.splitlines() if ln.startswith("[FAIL]")]
        fired = all(w in fails for w in want) and p.returncode == 1
        ok_all &= fired
        print(f"[{'FIRED' if fired else 'SILENT'}] --mutate={name} ({what}): exit {p.returncode}, failed {fails}", flush=True)
        if p.returncode not in (0, 1):
            print(p.stderr[-2000:])
    print(f"controls: {'all fired' if ok_all else 'SOME DID NOT FIRE'}")
    return ok_all


if __name__ == "__main__":
    opts = parse_opts(sys.argv[1:])
    MUT = opts.get("mutate")
    if MUT is not None and MUT not in MUTATIONS:
        sys.exit(f"unknown mutation {MUT!r}. Known: {', '.join(MUTATIONS)}")
    RAW = os.path.expanduser(opts.get("raw", RAW))
    if "compute" in opts:
        main_compute(opts)
    elif "ingest" in opts:
        for nm in opts["ingest"].split(","):
            ingest(nm)
    else:
        failed = main(opts)
        controls_ok = True
        if MUT is None and "no-controls" not in opts:
            controls_ok = run_controls(opts)
        if MUT is None:
            sys.exit(0 if (not failed and controls_ok) else 1)
        sys.exit(1 if failed else 0)
