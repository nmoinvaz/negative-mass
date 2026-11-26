#!/usr/bin/env python3
"""The static one-sign class with u_x rho != 0, decided by prolongation (issue #165)

#158 (one_sign_classify.md, "General static shifts with u_x rho != 0") reduced the class to an integral surface of the
rank-2 distribution X1 = D_x, X2 = D_rho on the jets (rho, u, a, beta, s), [X1, X2] = C d_s. A solution lies on C = 0, then
on X_i C = 0, and so on. This script

  1. rebuilds the system in the scale-invariant jets (u, b = rho beta, v = rho v2, s) with rho = 1 (a = v w, w = u c - sin),
  2. computes the chain f0 = C (numerator), f1 = X1 f0, f2 = X2 f0, f11 = X1 f1, f12 = X2 f1 exactly (sympy),
  3. runs the planted cases (the family u = c x + C ln rho + d at s = +-1; a toy system with a known inconsistency),
  4. verifies, without Singular, the certificates stored in the raw folder (`--raw`, default ../negative-mass-raw/I165-rpi5).

    python3 notes/warp/one_sign_static.py                  verify: checks, then every control
    python3 notes/warp/one_sign_static.py --no-controls    the checks only
    python3 notes/warp/one_sign_static.py --mutate=NAME    one control (its target checks must fail, exit 1)
    python3 notes/warp/one_sign_static.py --compute=STEP   Singular steps (chain | toy | slices | fullring), written to RAW

Imported read-only from one_sign_classify.py (#158): `static_blocks`, `v_block`, `rad_member`. Everything else is new here.
Singular is trusted for nothing: a certificate is accepted only if its identity holds exactly in this script.
"""
import argparse
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import time

import sympy as sp
from sympy import QQ
from sympy.polys.fields import field
from sympy.polys.rings import ring

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import one_sign_classify as osc  # noqa: E402

DEFAULT_RAW = os.path.normpath(os.path.join(HERE, "..", "..", "..", "negative-mass-raw", "I165-rpi5"))
SINGULAR = "/usr/bin/Singular"

# name: (checks that must fail, what the mutation does)
MUTATIONS = {
    "Phi1-sign": (["blocks", "chain", "denoms", "planted_s1", "slice_modp"], "Phi1 with the opposite sign of a^2"),
    "Psi-sign": (["blocks", "chain", "slice_modp"], "Psi with the opposite sign"),
    "AB-sign": (["AB"], "B built as PQY + Z(Q^2 + Y^2)"),
    "denoms-claim": (["denoms"], "claims the denominator 4 s u Delta' in place of 8 s u Delta'"),
    "Cpoints-s1": (["Cpoints"], "the three points drawn on s = 1, where C-bar = 0"),
    "plant-s2": (["planted_s1"], "the planted test at s = 2 in place of s = 1"),
    "toy-flat": (["toy"], "the flat toy is passed to the rejection test"),
    "sconst-kappa": (["sconst"], "kap = c sg in place of 2 c sg"),
    "cert-corrupt": (["slice_modp"], "one cofactor coefficient of the first stored slice certificate changed by 1"),
    "cert-satdrop": (["slice_modp"], "the saturation generator verified with u dropped from P"),
}
MUT = None
RESULTS = []
SKIPPED = []
NOT_MET = ["D1 slice_cert: exact certificate over Q for the saturated level-2 ideal at fixed s (none obtained; see the raw README, probes 11 and 13)",
           "D2 full_cert: exact certificate over Q for all s (none obtained; the 60 min full-ring run is recorded in certs/fullring.json)"]
INFO = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def info(text):
    INFO.append(text)
    print(f"[info] {text}", flush=True)


# --------------------------------------------------------------------------------------------------------------------
# polynomial helpers
# --------------------------------------------------------------------------------------------------------------------

def rat(c):
    """A ring coefficient (QQ element) as a sympy Rational."""
    return sp.Rational(int(c.numerator), int(c.denominator))


def qv(x):
    """A sympy Rational as a QQ element."""
    x = sp.Rational(x)
    return QQ(int(x.p), int(x.q))


def ser(poly):
    """A ring element as a sorted list of [exponents, numerator, denominator]."""
    return sorted([list(m), int(rat(c).p), int(rat(c).q)] for m, c in poly.terms())


def sha(poly):
    return hashlib.sha256(json.dumps(ser(poly)).encode()).hexdigest()


def from_ser(R, lst):
    return R.from_dict({tuple(m): QQ(n, d) for m, n, d in lst})


def normalize(R, n):
    """Primitive part with integer coefficients and a positive leading coefficient."""
    if n == 0:
        return n
    den = 1
    for c in n.coeffs():
        den = sp.ilcm(den, rat(c).q)
    n = n * den
    g = 0
    for c in n.coeffs():
        g = sp.igcd(g, int(c))
    n = n.quo_ground(g)
    return -n if n.LC < 0 else n


def strip(R, num, nv):
    """Remove the nonvanishing factors nv (ring elements) from the numerator, exactly, and normalize."""
    n = num
    if n == 0:
        return n
    for q in nv:
        while True:
            qq, rr = divmod(n, q)
            if rr == 0:
                n = qq
            else:
                break
    return normalize(R, n)


class Chain:
    """Rational vector fields X1, X2 (tuples of fraction-field elements) on the generators of the field K and the prolongation
    f -> numerator of X_i f with the nonvanishing factors nv removed."""

    def __init__(self, K, X1, X2, nv):
        self.K, self.R, self.X, self.nv = K, K.ring, (X1, X2), nv
        self.gens = K.gens

    def act(self, F, i):
        """X_i F as a fraction-field element, F a ring element."""
        Fe = self.K.new(F, self.R.one)
        return sum((Fe.diff(g) * c for g, c in zip(self.gens, self.X[i - 1])), self.K.zero)

    def step(self, F, i, stripped=True):
        g = self.act(F, i)
        return strip(self.R, g.numer, self.nv) if stripped else g.numer


# --------------------------------------------------------------------------------------------------------------------
# the static system
# --------------------------------------------------------------------------------------------------------------------

def build_system():
    """Fields at rho = 1 in the old jets (u, a, b, s), with the Euler term for d_rho (weights rho: 1, a, b: -1, so
    d_rho f = W f + a f_a + b f_b for f of weight W); sigma_x, sigma_rho from E1, E2; then the change a = v w to (u, b, v, s).
    Returns a dict with the 4D field, the fields X1c, X2c, sigma, C-bar (unreduced) and the old-coordinate pieces."""
    Ko, uo, ao, bo, so = field("u,a,b,s", QQ)
    N1 = 1 + so**2
    Dd = uo * (1 - so**2) - 2 * so
    v0, v1, v2 = bo * (1 - so**2) / N1, 2 * so * bo / N1, ao * N1 / Dd
    sg1 = 1 if MUT == "Phi1-sign" else -1
    Phi1 = (v0**2 + v2**2 + sg1 * ao**2) / uo
    Psi = -2 * v0 * v2 * (-1 if MUT == "Psi-sign" else 1)
    Phi3 = -bo + 2 * v0 * v1
    d = lambda f, x: f.diff(x)  # noqa: E731

    def Dx(f, sx):
        return d(f, uo) * ao + d(f, ao) * Phi1 + d(f, bo) * Psi + d(f, so) * sx

    def Dr(f, W, sr):
        return W * f + ao * d(f, ao) + bo * d(f, bo) + d(f, uo) * bo + d(f, ao) * Psi + d(f, bo) * Phi3 + d(f, so) * sr

    # E1 = Dr(Phi1) - Dx(Psi), E2 = Dr(Psi) - Dx(Phi3), linear in (sx, sr)
    x1, xs1 = Dx(Psi, 0), d(Psi, so)
    r1, rs1 = Dr(Phi1, -2, 0), d(Phi1, so)
    x2, xs2 = Dx(Phi3, 0), d(Phi3, so)
    r2, rs2 = Dr(Psi, -2, 0), d(Psi, so)
    E1z, E2z = r1 - x1, r2 - x2
    det = (-xs1) * rs2 - rs1 * (-xs2)
    sxo = (-(E1z * rs2 - rs1 * E2z)) / det
    sro = (-((-xs1) * E2z - (-xs2) * E1z)) / det

    K, u, b, v, s = field("u,b,v,s", QQ)
    R = K.ring
    N1n = 1 + s**2
    cc, sgm = (1 - s**2) / N1n, 2 * s / N1n
    w = u * cc - sgm
    A = v * w

    def to4(fo):
        """Old fraction-field element -> (u, b, v, s) with a = v w."""
        U, Aa, Bb, S = sp.symbols("u a b s")
        W_ = sp.Symbol("u") * (1 - S**2) / (1 + S**2) - 2 * S / (1 + S**2)
        num = fo.numer.as_expr().subs(Aa, sp.Symbol("v") * W_)
        den = fo.denom.as_expr().subs(Aa, sp.Symbol("v") * W_)
        return K.from_expr(num) / K.from_expr(den)

    sx, sr = to4(sxo), to4(sro)
    Phi1n, Psin = to4(Phi1), to4(Psi)
    w_s = (-(u * sgm + cc)) * 2 / N1n  # d w / d s

    def comp(Yu, Yb, Ys, YA):
        """Components (u, b, v, s) of a field whose old components are Y(u) = Yu, Y(b) = Yb, Y(s) = Ys, Y(a) = YA."""
        Yw = Yu * cc + w_s * Ys
        return (Yu, Yb, YA / w - v * Yw / w, Ys)

    X1c = comp(A, Psin, sx, Phi1n)
    X2c = comp(b, to4(Phi3) + b, sr, A + Psin)
    return dict(K=K, R=R, gens=(u, b, v, s), X1c=X1c, X2c=X2c, sx=sx, sr=sr, w=w, A=A, N1=N1n, c=cc, sg=sgm, to4=to4,
                old=dict(K=Ko, sx=sxo, sr=sro, Phi1=Phi1, Psi=Psi, Phi3=Phi3, gens=(uo, ao, bo, so)))


def Xapply(K, gens, F, Y):
    return sum((F.diff(g) * c for g, c in zip(gens, Y)), K.zero)


def chain_generators(sys_, with_f12=True):
    """Everything the chain needs: C-bar numerator (unstripped), f0, f1, f2, f11, f12 and the nonvanishing factors."""
    K, R = sys_["K"], sys_["R"]
    u, b, v, s = sys_["gens"]
    ru, rb, rv, rs = R.gens
    X1c, X2c, sx, sr = sys_["X1c"], sys_["X2c"], sys_["sx"], sys_["sr"]
    Cb = Xapply(K, sys_["gens"], sr, X1c) - Xapply(K, sys_["gens"], sx, X2c) + sx
    Dpol = ru * (1 - rs**2) - 2 * rs
    den_x = sx.denom
    Delta, rem = divmod(den_x, 8 * rs * ru)
    if Delta.LC < 0:
        Delta = -Delta
    nv = [ru, rs, rv, rb, rs - 1, rs + 1, rs**2 + 1, Dpol, Delta]
    ch = Chain(K, X1c, X2c, nv)
    f0 = strip(R, Cb.numer, nv)
    out = dict(Cb=Cb, Cnum=Cb.numer, f0=f0, nv=nv, Delta=Delta, Delta_rem=rem, Dpol=Dpol, chain=ch)
    out["f1"], out["f2"] = ch.step(f0, 1), ch.step(f0, 2)
    out["f11"] = ch.step(out["f1"], 1)
    if with_f12:
        out["f12"] = ch.step(out["f1"], 2)
    return out


# --------------------------------------------------------------------------------------------------------------------
# the checks
# --------------------------------------------------------------------------------------------------------------------

def check_blocks_AB():
    rho = osc.rho
    us, as_, bs, ss = sp.symbols("us as_ bs ss")
    cth, sth = (1 - ss**2) / (1 + ss**2), 2 * ss / (1 + ss**2)
    w0, w1, w2 = bs * cth, bs * sth, as_ / (us * cth - sth)
    Phi1 = (w0**2 + w2**2 + (1 if MUT == "Phi1-sign" else -1) * as_**2) / us
    Psi = -2 * w0 * w2 * (-1 if MUT == "Psi-sign" else 1)
    Phi3 = -bs / rho + 2 * w0 * w1
    ent = osc.static_blocks(us, as_, bs, Phi1, Psi, Phi3, rho)
    ident = all(sp.simplify(p - q) == 0 for p, q in zip(ent, osc.v_block(w0, w1, w2)))
    check("blocks", ident, "the block built from u_xx = Phi1, u_x rho = Psi, u_rho rho = Phi3 equals Lemma R's block "
          "-v v^T + lam eta (v0 = beta cos th, v1 = beta sin th, v2 = a/(u cos th - sin th), s = tan(th/2)), identically")
    v0, v1, v2 = sp.symbols("v0 v1 v2", real=True)
    Pq, Qq, Zq, Yq, hq = sp.symbols("Pq Qq Zq Yq hq")
    gensR = [Pq + v0**2 + v2**2, Qq + v0 * v1, Zq + v0 * v2, Yq + v1 * v2, 2 * hq - v0**2 - v1**2]
    varsR = (v0, v1, v2, Pq, Qq, Zq, Yq, hq)
    A_ = Qq * (Zq**2 + Yq**2) + 2 * hq * Zq * Yq
    B_ = Pq * Qq * Yq + (1 if MUT == "AB-sign" else -1) * Zq * (Qq**2 + Yq**2)
    mA = osc.rad_member(gensR, A_, [Zq, Yq, Qq], varsR)
    mB = osc.rad_member(gensR, B_, [Zq, Yq, Qq], varsR)
    check("AB", mA and mB, f"A = Q (Z^2 + Y^2) + 2 h Z Y in the radical: {mA}; B = P Q Y - Z (Q^2 + Y^2) in the radical: {mB} "
          "(Lemma R's ideal saturated by Z Y Q, Rabinowitsch)")


def check_fields(sys_, gen):
    K = sys_["K"]
    X1c, X2c = sys_["X1c"], sys_["X2c"]
    gens = sys_["gens"]
    comm = [Xapply(K, gens, X2c[k], X1c) - Xapply(K, gens, X1c[k], X2c) + X1c[k] for k in range(4)]
    u, b, v, s = gens
    ws = (-(u * sys_["sg"] + sys_["c"])) * 2 / sys_["N1"]
    Cb = gen["Cb"]
    # [X1, X2] + X1 = C-bar * d_s at fixed (u, A, B) = C-bar (d_s - (v w_s / w) d_v) in the coordinates (u, b, v, s)
    ok_ub = comm[0] == 0 and comm[1] == 0
    ok_v = comm[2] == Cb * (-v * ws / sys_["w"])
    ok_s = comm[3] == Cb
    check("fields", ok_ub and ok_v and ok_s, "[X1, X2] + X1 = C-bar (d_s - (v w_s/w) d_v) as rational functions of (u, b, v, s): the u and b components "
          f"vanish ({ok_ub}), the v component is -C-bar v w_s/w ({ok_v}), the s component is C-bar ({ok_s})")


def five_d_C(pt5):
    """C = D_rho sigma_x - D_x sigma_rho of #158 in 5D (rho, u, a, beta, s), from partial derivatives at an exact point."""
    K5, rho, u, a, b, s = field("rho,u,a,b,s", QQ)
    N1 = 1 + s**2
    Dd = u * (1 - s**2) - 2 * s
    v0, v1, v2 = b * (1 - s**2) / N1, 2 * s * b / N1, a * N1 / Dd
    sg1 = 1 if MUT == "Phi1-sign" else -1
    Phi1 = (v0**2 + v2**2 + sg1 * a**2) / u
    Psi = -2 * v0 * v2 * (-1 if MUT == "Psi-sign" else 1)
    Phi3 = -b / rho + 2 * v0 * v1
    d = lambda f, x: f.diff(x)  # noqa: E731

    def Dx(f, sx):
        return d(f, u) * a + d(f, a) * Phi1 + d(f, b) * Psi + d(f, s) * sx

    def Dr(f, sr):
        return d(f, rho) + d(f, u) * b + d(f, a) * Psi + d(f, b) * Phi3 + d(f, s) * sr

    x1, xs1 = Dx(Psi, 0), d(Psi, s)
    r1, rs1 = Dr(Phi1, 0), d(Phi1, s)
    x2, xs2 = Dx(Phi3, 0), d(Phi3, s)
    r2, rs2 = Dr(Psi, 0), d(Psi, s)
    E1z, E2z = r1 - x1, r2 - x2
    det = (-xs1) * rs2 - rs1 * (-xs2)
    sx = (-(E1z * rs2 - rs1 * E2z)) / det
    sr = (-((-xs1) * E2z - (-xs2) * E1z)) / det
    gens = (rho, u, a, b, s)
    vals0 = list(zip(K5.ring.gens, [qv(x) for x in pt5]))

    def ev(f):  # value of a field element at the point
        return f.numer.evaluate(vals0) / f.denom.evaluate(vals0)

    def pder(f, g):  # value of d f / d g, quotient rule on numerator and denominator (no cancellation)
        n, dn = f.numer, f.denom
        gg = K5.ring.gens[list(gens).index(g)]
        nv_, dv_ = n.evaluate(vals0), dn.evaluate(vals0)
        nd, dd = n.diff(gg).evaluate(vals0), dn.diff(gg).evaluate(vals0)
        return (nd * dv_ - nv_ * dd) / (dv_ * dv_)

    vals = [ev(Phi1), ev(Psi), ev(Phi3), ev(sx), ev(sr)]
    Phi1v, Psiv, Phi3v, sxv, srv = vals
    av, bv = qv(pt5[2]), qv(pt5[3])
    Dr_sx = pder(sx, rho) + pder(sx, u) * bv + pder(sx, a) * Psiv + pder(sx, b) * Phi3v + pder(sx, s) * srv
    Dx_sr = pder(sr, u) * av + pder(sr, a) * Phi1v + pder(sr, b) * Psiv + pder(sr, s) * sxv
    return Dr_sx - Dx_sr


def check_Cpoints(sys_, gen):
    random.seed(1651)
    K, R = sys_["K"], sys_["R"]
    Cb = gen["Cb"]
    ok, vals = True, []
    sset = (sp.Integer(1),) * 3 if MUT == "Cpoints-s1" else (sp.Rational(1, 2), sp.Rational(2), sp.Rational(3))
    for s0 in sset:
        for _try in range(200):
            u0, b0, v0, rho0 = (sp.Rational(random.randint(1, 9), random.randint(1, 4)) for _ in range(4))
            w0 = u0 * (1 - s0**2) / (1 + s0**2) - 2 * s0 / (1 + s0**2)
            Dv = u0 * (1 - s0**2) - 2 * s0
            den_ok = all(x != 0 for x in (Dv, w0))
            if not den_ok:
                continue
            pt4 = list(zip(R.gens, [qv(x) for x in (u0, b0, v0, s0)]))
            dn, nn, dlt = Cb.denom.evaluate(pt4), Cb.numer.evaluate(pt4), gen["Delta"].evaluate(pt4)
            if dn != 0 and dlt != 0:
                break
        else:
            raise SystemExit("no admissible point found")
        Cbar = rat(nn) / rat(dn)
        # 5D point at rho = rho0: beta = b0/rho0, a = A/rho0, A = v w; compare C(rho0, ...) = Cbar / rho0^2
        a5, b5 = v0 * w0 / rho0, b0 / rho0
        C5 = five_d_C((rho0, u0, a5, b5, s0))
        C5 = rat(C5)
        agree = (C5 * rho0**2 == -Cbar)
        ok &= agree and Cbar != 0
        vals.append((str(s0), str(Cbar), agree))
    check("Cpoints", ok, f"at 3 exact points (u, b, v, s) with s not in {{0, +-1}}: C-bar from the 4D chain equals -rho^2 C of the 5D formula "
          f"(C = D_rho sigma_x - D_x sigma_rho from partial derivatives at rho = rho0) and C-bar != 0; (s, agree): "
          f"{[(a, c) for a, _, c in vals]}")


def check_denoms(sys_, gen):
    K, R = sys_["K"], sys_["R"]
    ru, rb, rv, rs = R.gens
    sx, sr = sys_["sx"], sys_["sr"]
    Delta = gen["Delta"]
    cx = 4 if MUT == "denoms-claim" else 8
    ok1 = sx.denom == cx * rs * ru * Delta or sx.denom == -cx * rs * ru * Delta
    ok2 = sr.denom == cx * rs * ru * (1 + rs**2) * Delta or sr.denom == -cx * rs * ru * (1 + rs**2) * Delta
    ratio = gen["Cnum"]
    f0 = gen["f0"]
    prod = rv * rb * (rs**2 - 1) * gen["Dpol"] * f0
    cst = ratio.LC / prod.LC
    ok3 = ratio == cst * prod and cst != 0 and gen["Delta_rem"] == 0
    nterms = len(Delta.terms())
    check("denoms", ok1 and ok2 and ok3 and nterms == 26 and Delta.degree(rs) == 8,
          f"denominator of sigma_x = {cx} s u Delta' ({ok1}), of sigma_rho = {cx} s u (1+s^2) Delta' ({ok2}), Delta' has {nterms} terms and "
          f"degree {Delta.degree(rs)} in s; numer(C-bar) = c v b (s^2-1) D f0 with c = {cst} ({ok3})")
    return Delta


def check_chain(gen, rawdir, with_hash=True):
    names = ["f0", "f1", "f2", "f11", "f12"]
    want_terms = [319, 1388, 1818, 3837, 4892]
    got = [len(gen[n].terms()) for n in names]
    degs = [tuple(max(m[i] for m, _ in gen[n].terms()) for i in range(4)) for n in names]
    want_deg = [(6, 4, 4, 18), (9, 6, 6, 31), (10, 7, 6, 33), (12, 8, 8, 44), (13, 9, 8, 46)]
    vev = all(all(m[2] % 2 == 0 for m, _ in gen[n].terms()) for n in names)
    shas = {n: sha(gen[n]) for n in names}
    stored = None
    fn = os.path.join(rawdir, "chain", "chain.json")
    if os.path.exists(fn):
        stored = json.load(open(fn))["sha256"]
    hash_ok = stored is not None and all(stored.get(n) == shas[n] for n in names)
    check("chain", got == want_terms and degs == want_deg and vev and (hash_ok or not with_hash),
          f"term counts {got} (want {want_terms}), degrees (u,b,v,s) {degs == want_deg}, v even in all: {vev}, SHA-256 equal to chain/chain.json: "
          f"{hash_ok if stored is not None else 'file missing'}")
    return shas


def check_planted(sys_, gen):
    K, R = sys_["K"], sys_["R"]
    ru, rb, rv, rs = R.gens
    sx, sr = sys_["sx"], sys_["sr"]
    ch = gen["chain"]
    s_test = 2 if MUT == "plant-s2" else 1
    pts = (s_test, -s_test) if s_test == 1 else (2, -2)
    # (a) sigma numerators vanish at s = +-1; denominators do not vanish identically there
    a_ok = True
    for sv in pts:
        a_ok &= bool(sx.numer.evaluate(rs, sv).is_zero) and bool(sr.numer.evaluate(rs, sv).is_zero)
        for q in (ru, 1 + rs**2, gen["Dpol"], gen["Delta"], rs):
            a_ok &= not bool(q.evaluate(rs, sv).is_zero) if not isinstance(q, int) else True
    # (b) unstripped F0, F1, F2 vanish at s = +-1
    F0 = gen["Cnum"]
    F1 = ch.step(F0, 1, stripped=False)
    F2 = ch.step(F0, 2, stripped=False)
    b_ok = all(F.evaluate(rs, sv).is_zero for F in (F0, F1, F2) for sv in pts)
    # (c) explicit jets of u = c x + C ln rho + d at exact points: Lemma R's block (u_xx = u_x rho = 0, u_rho rho = -C/rho^2) and F_k = 0
    random.seed(1652)
    c_ok = True
    for _ in range(5):
        cc, CC, rho0, u0 = (sp.Rational(random.randint(1, 9), random.randint(1, 4)) for _ in range(4))
        # s = 1: v0 = 0, v1 = beta = C/rho, v2 = -a = -c
        P_, Q_, Z_, Y_, h_ = osc.static_blocks(u0, cc, CC / rho0, 0, 0, -CC / rho0**2, rho0)
        blk = osc.v_block(0, CC / rho0, -cc)
        c_ok &= all(sp.simplify(x - y) == 0 for x, y in zip((P_, Q_, Z_, Y_, h_), blk))
        pt = [qv(x) for x in (u0, CC, -rho0 * cc)] + [QQ(1)]
        c_ok &= all(F.evaluate(list(zip(R.gens, pt))) == 0 for F in (F0, F1, F2))
    # (d) the saturated generators do not vanish identically at s = 1
    d_ok = all(not f.evaluate(rs, 1).is_zero for f in (gen["f0"], gen["f1"], gen["f2"]))
    check("planted_s1", a_ok and b_ok and c_ok and d_ok,
          f"s = {s_test}: (a) sigma_x, sigma_rho numerators vanish at s = +-{s_test} and u, 1+s^2, D, Delta', s are nonzero there: {a_ok}; (b) the "
          f"unstripped F0, F1, F2 vanish at s = +-{s_test}: {b_ok}; (c) 5 exact jets of u = c x + C ln rho + d have Lemma R's block "
          f"and F0 = F1 = F2 = 0: {c_ok}; (d) f0, f1, f2 do not vanish at s = 1: {d_ok}")


def check_sconst():
    """s = tan(theta/2) constant = s* on an open set, s* not in {0, +-1} (c = cos theta*, sg = sin theta*, both nonzero). Then
    Q = -v0 v1 gives u_rho rho + u_rho/rho = kap u_rho^2 with kap = sin(2 theta*) = 2 c sg != 0, so B = rho u_rho has B_rho = kap B^2/rho
    and u = -ln(W)/kap + d(x), W = k(x) - kap ln(rho). Z and Y give u_x rho = -2 c u_rho a/(u c - sg) (a = u_x). The residual times
    rho W^2 equals -(k'(c u - sg) + 2 c k'/kap - 2 c W d'), linear in the independent functions ln W, W, 1 of rho: so k' = 0, d' = 0, u_x = 0."""
    x, L, c, sg, LW, Wv = sp.symbols("x L c sg LW Wv")
    kap = c * sg if MUT == "sconst-kappa" else 2 * c * sg
    k, d = sp.Function("k")(x), sp.Function("d")(x)
    W = k - kap * L
    u = -sp.log(W) / kap + d
    rho = sp.exp(L)
    uL, uLL = sp.diff(u, L), sp.diff(u, L, 2)
    u_r, u_rr = uL / rho, (uLL - uL) / rho**2
    eq_Q = sp.simplify(u_rr + u_r / rho - 2 * c * sg * u_r**2) == 0
    a = sp.diff(u, x)
    u_xr = sp.diff(uL, x) / rho
    resid = sp.simplify((u_xr * (u * c - sg) + 2 * c * u_r * a) * rho * W**2)
    kp, dp_ = sp.diff(k, x), sp.diff(d, x)
    want_sym = -(kp * (c * (-LW / kap + d) - sg) + 2 * c * kp / kap - 2 * c * Wv * dp_)
    eq_res = sp.simplify(resid - want_sym.subs({LW: sp.log(W), Wv: W})) == 0
    P = sp.Poly(sp.expand(want_sym), LW, Wv)
    coef_lnW, coef_W = P.coeff_monomial(LW), P.coeff_monomial(Wv)
    ok_coef = sp.simplify(coef_lnW - c * kp / kap) == 0 and sp.simplify(coef_W - 2 * c * dp_) == 0
    check("sconst", eq_Q and eq_res and ok_coef,
          f"for s = s* not in {{0, +-1}} (kap = 2 c sg != 0): u = -ln(W)/kap + d(x), W = k(x) - kap ln(rho) solves u_rho rho + u_rho/rho = kap u_rho^2 "
          f"({eq_Q}); (Z)(Y) residual times rho W^2 = -(k'(c u - sg) + 2 c k'/kap - 2 c W d') ({eq_res}); the coefficients of ln W and W are "
          f"c k'/kap and 2 c d' ({ok_coef}), so k' = d' = 0 and u_x = 0")


# --------------------------------------------------------------------------------------------------------------------
# toy systems through the same chain code
# --------------------------------------------------------------------------------------------------------------------

def toy_chain(flat=False):
    K, p, q, r = field("p,q,r", QQ)
    g = (K.zero if flat else p * (q**2 + 1))
    ch = Chain(K, (K.one, K.zero, K.zero), (K.zero, K.one, g), [p.numer])
    C0 = (K.zero if flat else q**2 + 1)
    gens = []
    if not flat:
        f0 = strip(K.ring, C0.numer, ch.nv)
        gens = [f0, ch.step(f0, 1), ch.step(f0, 2)]
    return ch, gens


def parse_q_poly(R, txt, names):
    """Singular output over Q -> ring element."""
    txt = txt.strip().replace(" ", "")
    if txt == "0":
        return R.zero
    out = {}
    for tm in re.findall(r"[+-]?[^+-]+", txt):
        sgn = 1
        if tm[0] in "+-":
            sgn = -1 if tm[0] == "-" else 1
            tm = tm[1:]
        coef = QQ(1)
        mon = [0] * len(names)
        for f in tm.split("*"):
            if re.fullmatch(r"\d+(/\d+)?", f):
                n, _, d = f.partition("/")
                coef *= QQ(int(n), int(d or 1))
            else:
                m = re.fullmatch(r"([A-Za-z]\w*)(?:\^(\d+))?", f)
                mon[names.index(m.group(1))] += int(m.group(2) or 1)
        out[tuple(mon)] = out.get(tuple(mon), QQ(0)) + sgn * coef
    return R.from_dict(out)


def check_toy(rawdir):
    ch, gens = toy_chain(flat=(MUT == "toy-flat"))
    flat_ch, flat_gens = toy_chain(flat=True)
    R = ch.R
    fn = os.path.join(rawdir, "certs", "toy.json")
    ok_cert, detail = False, "certs/toy.json missing"
    if os.path.exists(fn):
        rec = json.load(open(fn))
        cofs = [parse_q_poly(R, c, ["p", "q", "r"]) for c in rec["cofactors"]]
        gl = [from_ser(R, s) for s in rec["gens"]]
        same = (len(gens) > 0) and all(x == y for x, y in zip(gl, gens)) and len(gl) == len(gens)
        tot = sum((c * g for c, g in zip(cofs, gl)), R.zero)
        const = len(tot.terms()) == 1 and tot.terms()[0][0] == (0, 0, 0)
        ok_cert = same and const and tot != 0
        detail = f"stored generators equal the recomputed chain {[str(g.as_expr()) for g in gens]}: {same}; sum c_i f_i = {tot.as_expr()} (nonzero constant: {const})"
    ok_flat = len(flat_gens) == 0
    check("toy", ok_cert and ok_flat, f"inconsistent toy X2 = d_q + p (q^2 + 1) d_r: {detail}; flat toy X2 = d_q has an empty chain: {ok_flat}")


# --------------------------------------------------------------------------------------------------------------------
# Singular and modular certificates
# --------------------------------------------------------------------------------------------------------------------

def sing_poly(poly, names, sub=None):
    """A ring element (rational coefficients) as Singular text; sub = {index: integer} substitutes that variable."""
    terms = {}
    for m, c in poly.terms():
        m = list(m)
        cf = rat(c)
        if sub:
            for i, val in sub.items():
                cf *= sp.Rational(val) ** m[i]
                m[i] = 0
            m = [e for i, e in enumerate(m) if i not in sub]
        terms[tuple(m)] = terms.get(tuple(m), 0) + cf
    nm = [x for i, x in enumerate(names) if not (sub and i in sub)]
    out = []
    for m, cf in terms.items():
        if cf == 0:
            continue
        mono = "*".join(f"{n}^{e}" if e > 1 else n for n, e in zip(nm, m) if e)
        c_txt = f"({cf.p}/{cf.q})" if cf.q != 1 else f"({cf.p})"
        out.append(c_txt + ("*" + mono if mono else ""))
    return "+".join(out) if out else "0"


def parse_mod_poly(txt, p, names):
    txt = txt.strip().replace(" ", "")
    if txt == "0":
        return {}
    out = {}
    for tm in re.findall(r"[+-]?[^+-]+", txt):
        sgn = 1
        if tm[0] in "+-":
            sgn = -1 if tm[0] == "-" else 1
            tm = tm[1:]
        coef = 1
        mon = [0] * len(names)
        for f in tm.split("*"):
            if f.isdigit():
                coef *= int(f)
            else:
                m = re.fullmatch(r"([A-Za-z]\w*)(?:\^(\d+))?", f)
                mon[names.index(m.group(1))] += int(m.group(2) or 1)
        key = tuple(mon)
        out[key] = (out.get(key, 0) + sgn * coef) % p
    return {k: v for k, v in out.items() if v}


def dmul(a, b, p):
    out = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = tuple(x + y for x, y in zip(m1, m2))
            out[m] = (out.get(m, 0) + c1 * c2) % p
    return {k: v for k, v in out.items() if v}


def dadd(a, b, p):
    out = dict(a)
    for m, c in b.items():
        out[m] = (out.get(m, 0) + c) % p
    return {k: v for k, v in out.items() if v}


def slice_mod(poly, s0, p):
    """f(u, b, v, s = s0) mod p as {(t, u, b, v): coeff}."""
    out = {}
    for (mu, mb, mv, ms), c in poly.terms():
        cf = rat(c)
        assert cf.q == 1
        out[(0, mu, mb, mv)] = (out.get((0, mu, mb, mv), 0) + int(cf.p) * pow(s0, ms, p)) % p
    return {k: v for k, v in out.items() if v}


def check_slice_certs(gen, rawdir):
    """D1m: for each stored modular slice certificate (p, s0): sum c_i g_i = kappa != 0 in GF(p)[t, u, b, v]. Modular evidence
    (a statement about the reduction mod p at s = s0), not a certificate over Q: D1 (exact over Q) is NOT MET."""
    files = sorted(f for f in os.listdir(os.path.join(rawdir, "certs")) if f.startswith("slice_p")) if os.path.isdir(os.path.join(rawdir, "certs")) else []
    if not files:
        check("slice_modp", False, "no stored modular slice certificates in certs/")
        return
    all_ok, notes = True, []
    first = True
    for fn in files:
        rec = json.load(open(os.path.join(rawdir, "certs", fn)))
        p, s0 = rec["p"], rec["s0"]
        gl = [slice_mod(gen[n], s0, p) for n in ("f0", "f1", "f2", "f11")]
            # D(s0) = u (1 - s0^2) - 2 s0 and Delta'(s0) as polynomials in u, b, v
        dpol = {(0, 1, 0, 0): (1 - s0 * s0) % p, (0, 0, 0, 0): (-2 * s0) % p}
        dpol = {k: x for k, x in dpol.items() if x}
        Dl = slice_mod(gen["Delta"], s0, p)
        P = {(0, 1, 0, 0): 1}
        facs = [{(0, 0, 0, 1): 1}, {(0, 0, 1, 0): 1}, dpol, Dl]
        if MUT == "cert-satdrop":
            P = {(0, 0, 0, 0): 1}
        for f_ in facs:
            P = dmul(P, f_, p)
        tP = dmul({(1, 0, 0, 0): 1}, P, p)
        g5 = dadd({(0, 0, 0, 0): 1}, {k: (-x) % p for k, x in tP.items()}, p)
        gl.append(g5)
        cofs = []
        for i, lst in enumerate(rec["cofactors"]):
            cofs.append({tuple(m): c % p for m, c in lst})
        if MUT == "cert-corrupt" and first:
            m0 = next(iter(cofs[0]))
            cofs[0][m0] = (cofs[0][m0] + 1) % p
        first = False
        tot = {}
        for c, g in zip(cofs, gl):
            tot = dadd(tot, dmul(c, g, p), p)
        const = len(tot) == 1 and (0, 0, 0, 0) in tot
        ok = const and tot[(0, 0, 0, 0)] != 0
        all_ok &= ok
        notes.append(f"p = {p}, s0 = {s0}: sum = {'const ' + str(tot[(0, 0, 0, 0)]) if const else str(len(tot)) + ' terms'} ({'ok' if ok else 'FAIL'})")
    check("slice_modp", all_ok, "modular evidence: for each stored modular certificate (g = f0, f1, f2, f11 at s = s0 reduced mod p, and 1 - t u v b D Delta' at s0), "
          "sum c_i g_i is a nonzero constant in GF(p)[t, u, b, v]: " + "; ".join(notes))


def read_dims(rawdir):
    cd = os.path.join(rawdir, "certs")
    for fn in sorted(os.listdir(cd)) if os.path.isdir(cd) else []:
        if fn.startswith("dims_p"):
            for rec in json.load(open(os.path.join(cd, fn))):
                info(f"dims (mod {rec['p']}, s = {rec['s0']}, exit {rec['exit']}): " + "; ".join(rec["rows"]))


def read_fullring(rawdir):
    fn = os.path.join(rawdir, "certs", "fullring.json")
    if os.path.exists(fn):
        rec = json.load(open(fn))
        info(f"full-ring run of record (certs/fullring.json): {rec}")
        return rec
    info("no full-ring run of record in certs/ (certs/fullring.json missing)")
    return None


# --------------------------------------------------------------------------------------------------------------------
# --compute: Singular steps (not run by the default mode)
# --------------------------------------------------------------------------------------------------------------------

def run_singular(text, tag, rawdir, timeout):
    os.makedirs(os.path.join(rawdir, "singular"), exist_ok=True)
    inp = os.path.join(rawdir, "singular", f"{tag}.sing")
    out = os.path.join(rawdir, "singular", f"{tag}.out")
    tfile = os.path.join(rawdir, "singular", f"{tag}.time")
    for fn in (inp, out):
        if os.path.exists(fn):
            raise SystemExit(f"refusing to overwrite {fn}")
    open(inp, "w").write(text)
    t0 = time.time()
    with open(inp) as fi, open(out, "w") as fo:
        rc = subprocess.run(["/usr/bin/time", "-v", "-o", tfile, "nice", "-n", "10", "timeout", str(timeout), SINGULAR, "-q"], stdin=fi, stdout=fo,
                            stderr=subprocess.STDOUT).returncode
    return rc, time.time() - t0, out


def compute_chain(rawdir):
    sys_ = build_system()
    gen = chain_generators(sys_)
    os.makedirs(os.path.join(rawdir, "chain"), exist_ok=True)
    names = ["f0", "f1", "f2", "f11", "f12"]
    rec = {"sha256": {n: sha(gen[n]) for n in names}, "gens": {n: ser(gen[n]) for n in names}, "variables": ["u", "b", "v", "s"]}
    json.dump(rec, open(os.path.join(rawdir, "chain", "chain.json"), "w"))
    open(os.path.join(rawdir, "chain", "delta.txt"), "w").write(str(gen["Delta"].as_expr()) + "\n")
    print("chain written:", {n: len(gen[n].terms()) for n in names})


def compute_toy(rawdir):
    ch, gens = toy_chain()
    R = ch.R
    names = ["p", "q", "r"]
    text = "ring r=0,(p,q,r),dp;\noption(redSB);\nshort=0;\n" + "".join(f"poly g{i}={sing_poly(g, names)};\n" for i, g in enumerate(gens)) + \
        f"ideal J={','.join(f'g{i}' for i in range(len(gens)))};\nmatrix T; ideal G=liftstd(J,T);\n\"G1\",G[1];\n" \
        f"int i; for(i=1;i<=ncols(J);i++){{ \"C\",string(T[i,1]); }}\nquit;\n"
    rc, el, out = run_singular(text, "toy", rawdir, 120)
    lines = open(out).read().splitlines()
    cof = [ln.split(" ", 1)[1].strip().strip('"') for ln in lines if ln.startswith("C ")]
    os.makedirs(os.path.join(rawdir, "certs"), exist_ok=True)
    json.dump({"gens": [ser(g) for g in gens], "cofactors": cof}, open(os.path.join(rawdir, "certs", "toy.json"), "w"))
    print("toy certificate written, singular rc", rc, cofs_str(cof))


def cofs_str(c):
    return [x[:60] for x in c]


def compute_slices(rawdir, s0s, primes):
    gen = chain_generators(build_system(), with_f12=False)
    names = ("t", "u", "b", "v")
    os.makedirs(os.path.join(rawdir, "certs"), exist_ok=True)
    for p in primes:
        for s0 in s0s:
            tag = f"slice_p{p}_s{s0}"
            texts = []
            for n in ("f0", "f1", "f2", "f11"):
                texts.append(f"poly {n}=" + sing_poly(gen[n], ["u", "b", "v", "s"], sub={3: s0}) + ";\n")
            D_ = f"(u*(1-({s0})^2)-2*({s0}))"
            Dl = sing_poly(gen["Delta"], ["u", "b", "v", "s"], sub={3: s0})
            text = f"ring r={p},(t,u,b,v),dp;\noption(redSB);\nshort=0;\n" + "".join(texts) + f"poly Dl={Dl};\npoly P=u*v*b*{D_}*Dl;\n" + \
                "ideal J=f0,f1,f2,f11,1-t*P;\nmatrix T; ideal G=liftstd(J,T);\n\"G1\",G[1];\n\"n\",size(G);\n" \
                "int i; for(i=1;i<=ncols(J);i++){ \"C\",string(T[i,1]); }\nquit;\n"
            rc, el, out = run_singular(text, tag, rawdir, 900)
            lines = open(out).read().splitlines()
            cof = []
            for ln in lines:
                if ln.startswith("C "):
                    cof.append([[list(m), c] for m, c in parse_mod_poly(ln.split(" ", 1)[1], p, names).items()])
            g1 = [ln for ln in lines if ln.startswith("G1")]
            json.dump({"p": p, "s0": s0, "G1": g1, "cofactors": cof}, open(os.path.join(rawdir, "certs", tag + ".json"), "w"))
            print(f"{tag}: rc {rc}, {el:.1f} s, G1 {g1}, cofactor terms {[len(c) for c in cof]}", flush=True)


def compute_dims(rawdir, s0s, p):
    """Dimension and multiplicity of the saturated chain ideals at fixed s (mod p, Singular std), stored as info for the note."""
    gen = chain_generators(build_system(), with_f12=True)
    names = ["u", "b", "v", "s"]
    recs = []
    os.makedirs(os.path.join(rawdir, "certs"), exist_ok=True)
    for s0 in s0s:
        D_ = f"(u*(1-({s0})^2)-2*({s0}))"
        text = f"ring r={p},(t,u,b,v),dp;\noption(redSB);\nshort=0;\n" + "".join(f"poly {n}=" + sing_poly(gen[n], names, sub={3: s0}) + ";\n" for n in ("f0", "f1", "f2", "f11", "f12")) + \
            f"poly Dl={sing_poly(gen['Delta'], names, sub={3: s0})};\npoly Dp={D_};\npoly P=u*v*b*Dp*Dl;\n" + \
            "proc rep(string nm, ideal J){ ideal G=std(J); print(\"R \"+nm+\" dim \"+string(dim(G))+\" mult \"+string(mult(G))+\" unit \"+string(G[1]==1)); }\n" + \
            "rep(\"f0 plain\", ideal(f0));\nrep(\"f0 sat\", ideal(f0,1-t*P));\nrep(\"f0 f1 f2 plain\", ideal(f0,f1,f2));\n" \
            "rep(\"f0 f1 f2 sat\", ideal(f0,f1,f2,1-t*P));\nrep(\"f0 f1 f2 f11 sat\", ideal(f0,f1,f2,f11,1-t*P));\n" \
            "rep(\"f0 f1 f2 f12 sat\", ideal(f0,f1,f2,f12,1-t*P));\nrep(\"f0 f1 f2 f11 sat by u only\", ideal(f0,f1,f2,f11,1-t*u));\n" \
            "rep(\"f0 f1 f2 f11 sat by u v b D\", ideal(f0,f1,f2,f11,1-t*u*v*b*Dp));\nquit;\n"
        rc, el, out = run_singular(text, f"dims_p{p}_s{s0}", rawdir, 600)
        rows = [ln[2:] for ln in open(out).read().splitlines() if ln.startswith("R ")]
        recs.append({"p": p, "s0": s0, "exit": rc, "rows": rows})
        print(f"s0 = {s0}: rc {rc}, {el:.1f} s", *rows, sep="\n  ", flush=True)
    json.dump(recs, open(os.path.join(rawdir, "certs", f"dims_p{p}.json"), "w"))


def compute_fullring(rawdir, p, minutes):
    gen = chain_generators(build_system(), with_f12=False)
    names = ["u", "b", "v", "s"]
    Pfull = "u*v*b*s*(s^2-1)*(s^2+1)*(u*(1-s^2)-2*s)*Dl"
    text = f"ring r={p},(t,u,b,v,s),(dp(4),dp(1));\noption(redSB);\nshort=0;\n" + \
        "".join(f"poly {n}=" + sing_poly(gen[n], names) + ";\n" for n in ("f0", "f1", "f2", "f11")) + \
        f"poly Dl={sing_poly(gen['Delta'], names)};\npoly P={Pfull};\nideal J=f0,f1,f2,f11,1-t*P;\nint t0=rtimer;\nideal G=std(J);\n" \
        "\"std time\",rtimer-t0,\"size\",size(G),\"G1\",G[1];\nquit;\n"
    rc, el, out = run_singular(text, f"fullring_p{p}", rawdir, minutes * 60)
    rec = {"p": p, "order": "(dp(4),dp(1))", "exit": rc, "seconds": round(el, 1), "finished": rc == 0, "output_head": open(out).read()[:300]}
    os.makedirs(os.path.join(rawdir, "certs"), exist_ok=True)
    json.dump(rec, open(os.path.join(rawdir, "certs", "fullring.json"), "w"))
    print(rec)


# --------------------------------------------------------------------------------------------------------------------

def main():
    global MUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=DEFAULT_RAW)
    ap.add_argument("--mutate", default=None)
    ap.add_argument("--no-controls", action="store_true")
    ap.add_argument("--compute", default=None, choices=["chain", "toy", "slices", "dims", "fullring"])
    ap.add_argument("--s0", default="2,3,5,7,11")
    ap.add_argument("--primes", default="2147483629,2147483587")
    ap.add_argument("--minutes", type=int, default=60)
    args = ap.parse_args()
    t0 = time.time()
    if args.compute:
        if args.compute == "chain":
            compute_chain(args.raw)
        elif args.compute == "toy":
            compute_toy(args.raw)
        elif args.compute == "dims":
            compute_dims(args.raw, [int(x) for x in args.s0.split(",")], int(args.primes.split(",")[0]))
        elif args.compute == "slices":
            compute_slices(args.raw, [int(x) for x in args.s0.split(",")], [int(x) for x in args.primes.split(",")])
        else:
            compute_fullring(args.raw, int(args.primes.split(",")[0]), args.minutes)
        return 0
    if args.mutate:
        MUT = args.mutate
        if MUT not in MUTATIONS:
            raise SystemExit(f"unknown mutation {MUT}")
    print(f"one_sign_static.py{' --mutate=' + MUT if MUT else ''}  (raw: {args.raw})", flush=True)
    print("\n== 1. The system: blocks, Lemma R's A and B, the fields, C-bar ==", flush=True)
    check_blocks_AB()
    sys_ = build_system()
    print(f"[info] fields built ({time.time() - t0:.0f} s)", flush=True)
    gen = chain_generators(sys_)
    print(f"[info] chain built ({time.time() - t0:.0f} s)", flush=True)
    check_fields(sys_, gen)
    check_Cpoints(sys_, gen)
    check_denoms(sys_, gen)
    print("\n== 2. Planted cases ==", flush=True)
    check_planted(sys_, gen)
    check_sconst()
    check_toy(args.raw)
    print("\n== 3. The chain and its stored certificates ==", flush=True)
    check_chain(gen, args.raw)
    check_slice_certs(gen, args.raw)
    read_dims(args.raw)
    read_fullring(args.raw)
    bad = [n for n, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS)} checks, {len(RESULTS) - len(bad)} pass" + (f", FAILED: {bad}" if bad else "") +
          (f", skipped: {SKIPPED}" if SKIPPED else "") + f" ({time.time() - t0:.0f} s)")
    for line in NOT_MET:
        print(f"[NOT MET] {line}")
    if MUT:
        want = sorted(MUTATIONS[MUT][0])
        fired = sorted(bad)
        print(f"control {MUT}: target {want}, failed {fired}: {'ON TARGET' if fired == want else 'OFF TARGET' if fired else 'SILENT'}")
        return 1 if bad else 0
    if bad:
        return 1
    if not args.no_controls:
        silent = []
        for name in MUTATIONS:
            rc = subprocess.run([sys.executable, os.path.abspath(__file__), f"--mutate={name}", f"--raw={args.raw}"], capture_output=True, text=True)
            lastline = [ln for ln in rc.stdout.splitlines() if ln.startswith("control ")]
            on = bool(lastline) and "ON TARGET" in lastline[-1]
            print(f"  {name}: exit {rc.returncode}; {lastline[-1] if lastline else rc.stdout[-200:] + rc.stderr[-200:]}", flush=True)
            if not on:
                silent.append(name)
        print("controls: " + ("all fired on target" if not silent else f"NOT on target: {silent}"))
        return 1 if silent else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
