#!/usr/bin/env python3
"""#194: the radial boson-star modes and the phantom family by Taylor integration on multishoot and multisolver
(crate taylor_radial/).

    python taylor_radial.py                     build if needed, run every check, exit 1 if any fails
    python taylor_radial.py --only=P1,A         run the named checks only
    python taylor_radial.py --mutate=NAME       run with one deliberate error (a control: its targets must fail)
    python taylor_radial.py --gen               write taylor_radial/src/rgen.rs (taylor_radial_gen.py)
    python taylor_radial.py --build             cargo build --release (nice 10, -j 2)
    python taylor_radial.py --list              the checks and the controls
    python taylor_radial.py --raw=DIR           where the stdout of the long runs is kept (default: negative-mass-raw/I194-mac)

The criteria are the ones declared in negative-mass-raw/I194-mac/README.md before the runs of record (CRIT below repeats the
numbers). A check that fails is reported as failed; nothing here loosens a criterion after data. `--mutate=NAME` runs a
check with a deliberate error and the check it targets must fail (each control's targets are in the raw README, amendment A1).
"""
import json
import math
import os
import re
import subprocess
import threading
import sys
import time
from decimal import Decimal, ROUND_HALF_EVEN, getcontext

getcontext().prec = 60
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CRATE = os.path.join(HERE, "taylor_radial")
BIN = os.path.join(CRATE, "target", "release", "taylor_radial")
RESULTS = os.path.join(HERE, "RESULTS.md")
RAW_DEFAULT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "negative-mass-raw", "I194-mac"))
I10 = os.path.abspath(os.path.join(HERE, "..", "..", "..", "negative-mass-raw", "I10-rpi5"))

MUTATE = None
RAW = RAW_DEFAULT
OUT_RESULTS = []

CRIT = dict(
    P1_TOL={"f64": 1e-10, "dd": 1e-24, "qd": 1e-50},
    P2_TOL={"f64": 1e-9, "dd": 1e-24, "qd": 1e-48},
    P3_TOL=1e-14,
    G1_TOL=1e-24, G2_TOL=1e-20,
    S1_TOL=1e-7, S2_TOL=1e-28, S3_TOL=1e-26, S4_TOL=1e-28, S5_TOL=1e-18,
    X1_TOL=1e-12, X1_STORED=3e-8,
    Z_QD=1e-20, Z_IND=1e-18,
)


def report(name, ok, text):
    """ok is True (passed), False (failed) or None (skipped: it did not run, and counts as neither)"""
    ok = None if ok is None else bool(ok)
    OUT_RESULTS.append((name, ok))
    print(f"CHECK {name}: {'SKIPPED' if ok is None else 'PASS' if ok else 'FAIL'}  {text}", flush=True)


def gen():
    subprocess.check_call([sys.executable, os.path.join(HERE, "taylor_radial_gen.py")])


def build():
    env = dict(os.environ, CARGO_NET_GIT_FETCH_WITH_CLI="true")
    subprocess.check_call(["nice", "-n", "10", "cargo", "build", "--release", "-j", "2"], cwd=CRATE, env=env)


def rust(cmd, timeout=3600, tag=None, **kw):
    """runs the crate; returns the JSON lines of its stdout. The lines are echoed to the log and appended to RAW/<tag>.jsonl
    as they arrive, and the child (its own PID) is stopped at `timeout` s of wall time, the 1 h stop condition."""
    args = [BIN, cmd] + [f"{k}={v}" for k, v in kw.items()]
    if MUTATE:
        args.append(f"mutate={MUTATE}")
    env = dict(os.environ, RAYON_NUM_THREADS="2")
    if MUTATE:
        env["TR_MAX_SECONDS"] = "300"
    p = subprocess.Popen(["nice", "-n", "10"] + args, stdout=subprocess.PIPE, text=True, env=env)  # stderr: the driver's own
    timed_out = []

    def stop():
        timed_out.append(True)
        p.terminate()

    timer = threading.Timer(timeout, stop)
    timer.start()
    out, fh = [], None
    if tag and RAW and os.path.isdir(RAW):
        fh = open(os.path.join(RAW, f"{tag}.jsonl"), "a")
    try:
        for line in p.stdout:
            if line.startswith("{"):
                out.append(json.loads(line))
                print("  > " + line.rstrip(), flush=True)
                if fh:
                    fh.write(line)
                    fh.flush()
        p.wait()
    finally:
        timer.cancel()
        if fh:
            fh.close()
    if timed_out:
        raise RuntimeError(f"{' '.join(args)} stopped at the {timeout} s limit after {len(out)} lines")
    if p.returncode:
        raise RuntimeError(f"{' '.join(args)} exited {p.returncode}")
    return out


# ---------------------------------------------------------------------------------------------------------------
# tokens of the stored tables: a token matches when the crate's value rounds to it at the token's digits
# ---------------------------------------------------------------------------------------------------------------
def quant(token, value):
    """the crate's value (str or number) formatted with the digits of `token`; also the distance of the value from the
    token in half units of the token's last digit (<= 1 means the same rounding)"""
    t = token.strip().replace("−", "-")
    v = Decimal(str(value))
    m = re.fullmatch(r"([+-]?)(\d+\.?\d*)(?:[eE]([+-]?\d+))?", t)
    if not m:
        raise ValueError(token)
    mant, ex = m.group(2), int(m.group(3) or 0)
    dec = len(mant.split(".")[1]) if "." in mant else 0
    if m.group(3) is None:
        unit = Decimal(10) ** (-dec)
        rounded = v.quantize(unit, rounding=ROUND_HALF_EVEN)
        fmt = f"{'+' if t.startswith('+') and rounded >= 0 else ''}{rounded:.{dec}f}"
        tv = Decimal(t)
    else:
        unit = Decimal(10) ** (ex - dec)
        scaled = v / Decimal(10) ** ex
        rounded = scaled.quantize(Decimal(10) ** (-dec), rounding=ROUND_HALF_EVEN)
        if abs(rounded) >= 10:
            ex += 1
            scaled = v / Decimal(10) ** ex
            rounded = scaled.quantize(Decimal(10) ** (-dec), rounding=ROUND_HALF_EVEN)
        fmt = f"{'+' if t.startswith('+') and rounded >= 0 else ''}{rounded:.{dec}f}e{ex:+d}".replace("e+", "e").replace("e-0", "e-")
        tv = Decimal(t)
        unit = Decimal(10) ** (int(m.group(3) or 0) - dec)
    half = abs((v - tv) / (unit / 2))
    return fmt, float(half)


def tok_ok(token, value):
    t = token.strip().replace("−", "-")
    fmt, half = quant(t, value)
    return half <= 1.0 + 1e-12, fmt, half


# ---------------------------------------------------------------------------------------------------------------
# stored tables
# ---------------------------------------------------------------------------------------------------------------
def stored_scan():
    """RESULTS.md section 3: [(w, M token, N token, lowest sigma^2 token or None)]"""
    text = open(RESULTS).read()
    sec = text.split("## 3. Validation on the canonical family")[1].split("Fine scan")[0]
    rows = []
    for ln in sec.splitlines():
        m = re.match(r"\|\s*([0-9.]+)\s*\|\s*([+-][0-9.]+)\s*\|\s*([0-9.]+)\s*\|\s*(.*?)\s*\|\s*$", ln)
        if not m or m.group(1) in ("ω",):
            continue
        w = m.group(1)
        c4 = m.group(4).replace("*", "").replace("−", "-")
        r = re.search(r"([+-]?\d\.\d+e[+-]\d+)", c4)
        rows.append((float(w), m.group(2), m.group(3), r.group(1) if r else None))
    return rows


def stored_crossing():
    """RESULTS.md: the eps = +1 and eps = -1 crossing tables, as strings"""
    text = open(RESULTS).read()
    out = {}
    for eps, head in ((1, "Fine scan with 17 stars"), (-1, "Fine scan (`crossing_m.txt`):")):
        sec = text.split(head)[1].split("\n\n")[1]
        vals = re.findall(r"\|\s*(dM/dω = 0|dN/dω = 0|σ² = 0)\s*\|\s*\**([0-9.]+)\**\s*\|", sec)
        out[eps] = {k: v for k, v in vals}
    return out


def stored_modes():
    """RESULTS.md section 7: {alpha0: (w token, thr token, [root tokens])}"""
    text = open(RESULTS).read()
    sec = text.split("### Radial stability along the branch")[1].split("No root is negative")[0]
    out = {}
    for ln in sec.splitlines():
        m = re.match(r"\|\s*([0-9.]+)\s*\|\s*([0-9.]+)\s*\|\s*([0-9.e+-]+)\s*\|\s*(.*?)\s*\|\s*[0-9]+ / [0-9]+\s*\|", ln)
        if m and m.group(1) != "α₀":
            roots = [] if m.group(4).strip() == "none" else [x.strip() for x in m.group(4).split(",")]
            out[float(m.group(1))] = (m.group(2), m.group(3), roots)
    return out


def stored_family():
    """phantom_family.txt: {alpha_c: dict of tokens}"""
    out = {}
    for ln in open(os.path.join(HERE, "phantom_family.txt")):
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split()
        out[float(f[0])] = dict(alpha=f[0], w=f[1], k=f[2], N=f[3], M=f[4], max2mr=f[5], rs=f[6], R=f[7])
    return out


def stored_modes_full():
    """the full-precision roots of I10-rpi5/modes.jsonl, {alpha_c: [roots]}"""
    p = os.path.join(I10, "modes.jsonl")
    if not os.path.exists(p):
        return {}
    return {round(json.loads(x)["alpha_c"], 4): json.loads(x)["roots"] for x in open(p) if x.strip()}


# ---------------------------------------------------------------------------------------------------------------
# P1: the planted eigenproblem through the same code, in each lane
# ---------------------------------------------------------------------------------------------------------------
def check_P1(lanes=("f64", "dd", "qd")):
    ok_all, parts = True, []
    for lane in lanes:
        t0 = time.time()
        lines = rust("planted", lane=lane, n=60)
        head = [l for l in lines if l["cmd"] == "planted"][0]
        roots = [l for l in lines if l["cmd"] == "planted_root"]
        worst = max([r["rel"] for r in roots], default=float("inf"))
        ok = head["nroots"] == head["nexact"] == len(roots) and worst <= CRIT["P1_TOL"][lane]
        # the planted problem keeps an arbitrary start plane at R: its orientation (the sign of det C) must not change on the
        # scanned range, or a sign change of D would not be a root
        import mpmath as mp
        mp.mp.dps = 30
        orient = {planted_orientation(mp.mpf(-4) + mp.mpf(59 - 0) * 0 + i * mp.mpf("5.9") / 59, head["rfar"]) for i in range(60)}
        ok &= len(orient) == 1
        ok_all &= ok
        parts.append(f"{lane}: {head['nroots']} of {head['nexact']} roots, worst rel {worst:.2e} (limit {CRIT['P1_TOL'][lane]:g}), "
                     f"{head['wall']:.1f} s, R = {head['rfar']:.0f}, orientation of the start plane constant on the grid: {len(orient) == 1}")
    report("P1", ok_all, "; ".join(parts))


# ---------------------------------------------------------------------------------------------------------------
# A: the canonical scan table (RESULTS.md section 3)
# ---------------------------------------------------------------------------------------------------------------
def star_line(lines):
    return [l for l in lines if l["cmd"] == "star"][0]


def check_A(omegas=None, lane="dd"):
    stored = stored_scan()
    full = {round(r["w"], 4): r for r in json.load(open(os.path.join(HERE, "scan_eps+1.json")))}
    nbad, ntok, worst_half, out = 0, 0, 0.0, []
    for w, mtok, ntok_s, stok in stored:
        if omegas and w not in omegas:
            continue
        lines = rust("row", omega=repr(w), eps=1, lane=lane)
        st = star_line(lines)
        roots = [l["s2"] for l in lines if l["cmd"] == "root"]
        done = [l for l in lines if l["cmd"] == "scan_done"][0]
        res = []
        for name, tok, val in (("M", mtok, st["M"]), ("N", ntok_s, st["N"])):
            good, fmt, half = tok_ok(tok, val)
            ntok += 1
            worst_half = max(worst_half, half)
            nbad += not good
            res.append(f"{name} {tok}{'=' if good else '!='}{fmt} ({half:.2f})")
        if stok is None:
            good = not roots
            ntok += 1
            nbad += not good
            res.append(f"sigma^2 none: crate {'none' if not roots else roots[0]}")
        else:
            good, fmt, half = tok_ok(stok, roots[0]) if roots else (False, "none", float("inf"))
            ntok += 1
            worst_half = max(worst_half, half if half != float("inf") else 0)
            nbad += not good
            sr = full[round(w, 4)]["roots"]
            rel = abs(float(roots[0]) / sr[0][0] - 1) if roots and sr else float("nan")
            res.append(f"sigma^2 {stok}{'=' if good else '!='}{fmt} ({half:.2f}), {len(roots)} roots, vs stored double rel {rel:.1e}")
        out.append(f"  w={w:.3f}: " + "; ".join(res) + f"; fnorm {st['fnorm']:.1e}, {st['nodes']} nodes, {st['table_steps']} table steps")
    print("\n".join(out), flush=True)
    report("A", nbad == 0, f"{ntok - nbad} of {ntok} tokens equal their stored print; the largest distance from a stored token is {worst_half:.3f} half units of its last digit (limit 1)")


# ---------------------------------------------------------------------------------------------------------------
# mpmath: the generated J evaluated at 60 digits (G1), the planted reference (P2), the extremum fit (P3, X1)
# ---------------------------------------------------------------------------------------------------------------
def mp_J_entries(r, alpha, a, phi0, u, omega, eps, s2, lam, mu, pi_mp=None, dps=60):
    """J(r; s2) of y' = J y from radial_rhs_generated.py evaluated with mpmath numbers: J[i][j] = row i of the function at the
    unit vector e_j (the system is linear in y), row 0 being A' = dA"""
    import mpmath as mp
    mp.mp.dps = dps
    src = open(os.path.join(HERE, "radial_rhs_generated.py")).read()
    src = src.replace("from numpy import pi, exp, sqrt, log", "from mpmath import pi, exp, sqrt, log, mpf")
    # the generated source writes rational constants as float divisions (8/3): exact rationals here
    src = re.sub(r"(?<![\w.])(\d+)/(\d+)(?![\w.])", r"(mpf(\1)/\2)", src)
    ns = {}
    exec(compile(src, "radial_rhs_generated_mp", "exec"), ns)
    if pi_mp is not None:
        ns["pi"] = pi_mp
    args = dict(r=r, alpha=alpha, a=a, phi0=phi0, u=u, omega=omega, epsilon=eps, s2=s2, Lambda=lam, mu=mu)
    J = [[mp.mpf(0)] * 4 for _ in range(4)]
    J[0][1] = mp.mpf(1)
    for j in range(4):
        y = [mp.mpf(0)] * 4
        y[j] = mp.mpf(1)
        yy = dict(A=y[0], dA=y[1], Gt=y[2], chi=y[3])
        J[1][j] = ns["ddA"](**args, **yy)
        J[2][j] = ns["dGt"](**args, **yy)
        J[3][j] = ns["dchi"](**args, **yy)
    return J


def poly_extremum(ws, vals, deg=None):
    if MUTATE == "fit_deg":
        deg = 3
    """the extremum of the interpolating polynomial through (ws, vals): mpmath, to the digits of the inputs. ws, vals are
    mpf lists; returns the root of the derivative nearest the middle of the window"""
    import mpmath as mp
    n = len(ws)
    deg = n - 1 if deg is None else deg
    mid = ws[n // 2]
    xs = [w - mid for w in ws]
    # least-squares (exact interpolation when deg = n - 1) by normal equations in high precision
    A = mp.matrix(n, deg + 1)
    for i in range(n):
        for j in range(deg + 1):
            A[i, j] = xs[i] ** j
    b = mp.matrix(vals)
    c = mp.lu_solve(A.T * A, A.T * b) if deg < n - 1 else mp.lu_solve(A, b)
    dc = [c[j] * j for j in range(1, deg + 1)]
    # Newton on the derivative from the best grid point
    f = lambda x: sum(dc[j] * x ** j for j in range(len(dc)))
    x0 = min(xs, key=lambda x: abs(f(x)))
    root = mp.findroot(f, x0, tol=mp.mpf(10) ** (-mp.mp.dps + 8), maxsteps=200)
    return mid + root, c, mid


def check_P3():
    import mpmath as mp
    mp.mp.dps = 50
    ws0 = mp.mpf("0.92584444444444444444444444444444")
    N = lambda w: mp.mpf("0.0078099") + 3 * (w - ws0) ** 2 - 40 * (w - ws0) ** 3 + 500 * (w - ws0) ** 4
    pts = [ws0 + (i - 4) * mp.mpf("2.5e-4") + mp.mpf("1.1e-4") for i in range(9)]
    deg = 3 if MUTATE == "fit_deg" else None
    if MUTATE == "fit_deg":
        # a degree-3 least-squares fit of the planted quartic cannot recover the extremum
        pass
    wx, _, _ = poly_extremum(pts, [N(w) for w in pts], deg=deg)
    err = abs(wx - ws0)
    report("P3", err <= CRIT["P3_TOL"], f"planted extremum recovered to {float(err):.2e} (limit {CRIT['P3_TOL']:g}); degree {'8' if deg is None else deg}")


# ---------------------------------------------------------------------------------------------------------------
# G1 (J against mpmath), P2 (planted D against mpmath), G2, S1, S2
# ---------------------------------------------------------------------------------------------------------------
def mp_dd(s):
    import mpmath as mp
    return mp.mpf(s)


def check_G1(lane="dd"):
    import mpmath as mp
    mp.mp.dps = 70
    pts = [  # r, alpha, a, phi0, u, omega, eps, s2: exactly representable decimals, in the range of the canonical star
        ("0.5", "0.9921875", "1.0009765625", "0.0185546875", "-0.0302734375", "0.75", "1", "0.0625"),
        ("3.5", "0.9921875", "1.0078125", "0.0078125", "-0.4375", "0.75", "1", "0.0625"),
        ("12.25", "0.9990234375", "1.0009765625", "0.0001220703125", "-0.65625", "0.9375", "1", "-0.00390625"),
        ("40.0", "1.0", "1.0004882813", "1.5e-9", "-0.375", "0.9375", "-1", "0.0009765625"),
        ("1.0", "0.984375", "1.015625", "0.03125", "-0.125", "0.8125", "1", "-0.0625"),
        ("2.25", "0.9921875", "1.0234375", "0.015625", "-0.28125", "0.875", "1", "0.015625"),
        ("6.5", "0.99609375", "1.00390625", "0.00390625", "-0.5", "0.90625", "1", "0.00390625"),
        ("9.0", "0.998046875", "1.001953125", "0.001953125", "-0.59375", "0.9375", "1", "-0.0078125"),
        ("20.0", "0.9990234375", "1.0009765625", "3.0517578125e-05", "-0.71875", "0.75", "-1", "0.03125"),
        ("0.75", "0.9453125", "1.0078125", "0.0625", "-0.0625", "0.625", "-1", "0.0078125"),
        ("30.0", "1.0", "1.0009765625", "1.9073486328125e-06", "-0.3125", "0.9375", "1", "0.0009765625"),
        ("4.75", "0.99609375", "1.0126953125", "0.0107421875", "-0.40625", "0.84375", "-1", "-0.03125"),
    ]
    worst, nent = 0.0, 0
    for (r, al, a, ph, u, om, eps, s2) in pts:
        lines = rust("jeval", lane=lane, r=r, alpha=al, a=a, phi0=ph, u=u, omega=om, eps=eps, s2=s2)
        const = [l for l in lines if l["cmd"] == "J_const"][0]
        ref = mp_J_entries(*(mp.mpf(x) for x in (r, al, a, ph, u, om, eps, s2)), lam=mp.mpf(10240), mu=mp.mpf(const["mu"]), pi_mp=mp.mpf(const["pi"]))
        for l in [x for x in lines if x["cmd"] == "J"]:
            got, want = mp.mpf(l["v"]), ref[l["i"]][l["j"]]
            if want == 0 and got == 0:
                continue
            nent += 1
            worst = max(worst, float(abs(got - want) / max(abs(want), mp.mpf(10) ** -300)))
    report("G1", worst <= CRIT["G1_TOL"], f"J of the crate against radial_rhs_generated.py in mpmath (70 digits): worst relative {worst:.2e} over {nent} entries at {len(pts)} points (limit {CRIT['G1_TOL']:g})")


def check_G2():
    lines = rust("gcheck", omega=0.94, eps=1, x0="0.37,2.1,5.5,9.9,14.2,23.7,31.0,60.0", s2=0.003)
    w = [l for l in lines if l["cmd"] == "G2_done"][0]["worst_first12"]
    report("G2", w <= CRIT["G2_TOL"], f"the table's shifted jets against the direct series of J from the star's series, orders 0 to 11, at 8 radii: worst relative {w:.2e} (limit {CRIT['G2_TOL']:g})")


def planted_reference(s, lanes_dps=70):
    """D of the planted problem by mpmath: Legendre functions for each block's regular and decaying solutions, the planes mapped
    by S(r_m), orthonormalized once, oriented as the crate's start plane at R orients it. s is an mpf."""
    import mpmath as mp
    mp.mp.dps = lanes_dps
    lam = [mp.mpf(18) / 5, mp.mpf(27) / 10]
    b = [mp.mpf(4), mp.mpf(3)]
    c = [mp.mpf(2), mp.mpf(1)]
    rm, rfar = mp.mpf(3), None

    def block(i):
        kappa = mp.sqrt(b[i] - c[i] * s)
        L = lam[i]
        P = lambda r: mp.legenp(L, -kappa, mp.tanh(r), type=2)
        dP = lambda r: mp.diff(P, r)
        dec = (P(rm), dP(rm))
        # regular: odd combination P(x) - P(-x)
        Q = lambda r: P(r) - mp.legenp(L, -kappa, -mp.tanh(r), type=2)
        dQ = lambda r: mp.diff(Q, r)
        reg = (Q(rm), dQ(rm))
        q0 = dQ(mp.mpf(0))
        if q0 < 0:
            reg = (-reg[0], -reg[1])
        if (Q(mp.mpf(0))) != 0 and abs(Q(mp.mpf(0))) > mp.mpf(10) ** (-mp.mp.dps + 10):
            raise RuntimeError("the odd combination does not vanish at 0")
        return reg, dec, kappa

    (reg1, dec1, k1), (reg2, dec2, k2) = block(0), block(1)

    def S(r):
        g1 = mp.mpf(2) / 5 / (1 + (r / 8) ** 2)
        g2 = mp.mpf(3) / 10 + r / 20
        A = mp.eye(4)
        A[1, 2] = g1
        B = mp.eye(4)
        B[2, 1] = g2
        return A * B
    Sm = S(rm)
    yv = [Sm * mp.matrix([reg1[0], reg1[1], 0, 0]), Sm * mp.matrix([0, 0, reg2[0], reg2[1]])]
    zv = [Sm * mp.matrix([dec1[0], dec1[1], 0, 0]), Sm * mp.matrix([0, 0, dec2[0], dec2[1]])]

    def gs(cols):
        q = []
        for v in cols:
            w = v.copy()
            for u in q:
                w = w - (u.T * w)[0] * u
            q.append(w / mp.norm(w))
        return q
    qy, qz = gs(yv), gs(zv)
    M = mp.matrix(4, 4)
    for j in range(2):
        for i in range(4):
            M[i, j] = qy[j][i]
            M[i, 2 + j] = qz[j][i]
    D = mp.det(M)
    return D, (k1, k2)


def planted_orientation(s, r_far):
    """sign of det C: the components of the start plane's columns (the crate's fixed generic plane at r_far, mapped to the tilde
    coordinates) along each block's decaying solution (1, -kappa)"""
    import mpmath as mp
    b = [mp.mpf(4), mp.mpf(3)]
    c = [mp.mpf(2), mp.mpf(1)]
    r = mp.mpf(r_far)
    g1 = mp.mpf(2) / 5 / (1 + (r / 8) ** 2)
    g2 = mp.mpf(3) / 10 + r / 20
    A = mp.eye(4)
    A[1, 2] = g1
    B = mp.eye(4)
    B[2, 1] = g2
    Sinv = (A * B) ** -1
    z0 = [[1.0, 0.3, 0.0, 0.2], [0.0, 0.5, 1.0, 0.1]]
    C = mp.matrix(2, 2)
    for j in range(2):
        zt = Sinv * mp.matrix([mp.mpf(x) for x in z0[j]])
        for i in range(2):
            kappa = mp.sqrt(b[i] - c[i] * s)
            a_, b_ = zt[2 * i], zt[2 * i + 1]
            C[j, i] = (a_ - b_ / kappa) / 2
    return 1 if mp.det(C) > 0 else -1


def check_P2(lanes=("f64", "dd", "qd")):
    import mpmath as mp
    svals = ["-0.9", "-0.3125", "0.0625", "0.75", "1.5"]
    parts, ok_all = [], True
    for lane in lanes:
        dps = {"f64": 40, "dd": 70, "qd": 110}[lane]
        mp.mp.dps = dps
        lines = rust("pd", lane=lane, slist=",".join(svals))
        worst = 0.0
        for l in lines:
            s = mp.mpf(l["s"])
            ref, (k1, k2) = planted_reference(s, dps)
            kmin = float(min(k1, k2))
            r_far = 3.0 + {"f64": 19.0, "dd": 40.0, "qd": 74.0}[lane] / min(math.sqrt(4 - 2 * 1.9), math.sqrt(3 - 1.9))
            ref = ref * planted_orientation(s, r_far)
            got = mp.mpf(l["d"])
            worst = max(worst, float(abs(got - ref) / abs(ref)))
        ok = worst <= CRIT["P2_TOL"][lane]
        ok_all &= ok
        parts.append(f"{lane}: worst relative {worst:.2e} (limit {CRIT['P2_TOL'][lane]:g})")
    report("P2", ok_all, "planted D at 5 non-eigenvalue s against mpmath (Legendre solutions, one orthonormalization): " + "; ".join(parts))


# ---------------------------------------------------------------------------------------------------------------
# S1 (the star against scipy's), S4/S5 (residual, first law), X1 (the crossing)
# ---------------------------------------------------------------------------------------------------------------
def check_S1():
    import numpy as np
    import warnings
    warnings.filterwarnings("ignore")
    import background as B
    radii = [0.002, 0.05, 0.5, 2.0, 5.0, 8.0, 12.0, 16.0, 20.0, 23.0]
    worst, parts = 0.0, []
    for eps, w in ((1, 0.75), (-1, 0.7603)):
        st = B.solve_star(w, eps)
        lines = rust("starcheck", omega=repr(w), eps=eps, r=",".join(map(str, radii)))
        wmax = 0.0
        for l in [x for x in lines if x["cmd"] == "starcheck"]:
            al, a, phi0, u = st.bg(l["r"])
            for name, ref, got in (("alpha", al, l["alpha"]), ("a", a, l["a"]), ("phi0", phi0, l["phi0"]), ("u", u, l["u"])):
                rel = abs(float(got) / ref - 1) if ref != 0 else abs(float(got))
                wmax = max(wmax, rel)
        worst = max(worst, wmax)
        parts.append(f"eps={eps:+d} w={w}: worst relative {wmax:.2e}")
    report("S1", worst <= CRIT["S1_TOL"], f"the crate's star against background.solve_star (scipy BVP, tol 1e-9) at {len(radii)} radii for alpha, a, phi0, u: " + "; ".join(parts) + f" (limit {CRIT['S1_TOL']:g})")


def omega_window(eps, center, h, n, **kw):
    """stars at omega = center + (i - n//2) h: [(omega mpf, M mpf, N mpf, D0 mpf)]"""
    import mpmath as mp
    out = []
    for i in range(n):
        w = center + (i - n // 2) * h
        lines = rust("row", omega=mp.nstr(w, 40), eps=eps, s2list="0", dump=1, **kw)
        st = star_line(lines)
        d0 = [l for l in lines if l["cmd"] == "D"][0]["D"]
        out.append((mp.mpf(mp.nstr(w, 40)), mp.mpf(st["M"]), mp.mpf(st["N"]), mp.mpf(d0), float(st["fnorm"]), mp.mpf(st["omega"])))
    return out


def check_S5(eps=1, center="0.93", h="1e-4", n=13):
    import mpmath as mp
    mp.mp.dps = 50
    rows = omega_window(eps, mp.mpf(center), mp.mpf(h), n)
    ws = [r[5] for r in rows]
    mid = ws[n // 2]
    xs = [w - mid for w in ws]

    def deriv(vals):
        A = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                A[i, j] = xs[i] ** j
        c = mp.lu_solve(A, mp.matrix(vals))
        return c[1]
    dM, dN = deriv([r[1] for r in rows]), deriv([r[2] for r in rows])
    rel = abs(dM - eps * mid * dN) / abs(dM)
    fmax = max(r[4] for r in rows)
    report("S4", fmax <= CRIT["S4_TOL"], f"largest scaled Newton residual of the {n} stars of the window at w = {center}, eps = {eps:+d}: {fmax:.1e} (limit {CRIT['S4_TOL']:g})")
    report("S5", rel <= CRIT["S5_TOL"], f"first law dM = eps w dN at w = {mp.nstr(mid, 12)}: |dM/dw - eps w dN/dw| / |dM/dw| = {float(rel):.2e} from {n} stars at spacing {h}, degree {n - 1} (limit {CRIT['S5_TOL']:g})")


def secant_root(f, x0, x1, tol, maxit=30):
    import mpmath as mp
    f0, f1 = f(x0), f(x1)
    for _ in range(maxit):
        if f1 == f0:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0 = x1, f1
        x1 = x2
        f1 = f(x1)
        if abs(x1 - x0) < tol:
            break
    return x1, (x0, f0, f1)


def check_S3():
    """the star's central data and M: the dd boundary-value solution against the crate's own quad-double Taylor shooting (star_ref.rs)"""
    parts, worst = [], 0.0
    for eps, om in ((1, "0.75"), (-1, "0.7603")):
        l = [x for x in rust("s3", eps=eps, omega=om, rnorm=60, tol="1e-62") if x["cmd"] == "S3"][0]
        worst = max(worst, l["worst"])
        parts.append(f"eps={eps:+d} w={om}: P_c {l['rel_pc']:.1e}, Z_c {l['rel_zc']:.1e}, M {l['rel_M']:.1e} (qd shooting {l['t_qd']:.0f} s)")
    report("S3", worst <= CRIT["S3_TOL"], "the dd boundary-value star against the qd Taylor shooting (alpha = 1 at r = 60, tol 1e-62; not a qd boundary-value solution), worst relative difference "
           f"{worst:.1e} (limit {CRIT['S3_TOL']:g}): " + "; ".join(parts))


def check_E2(arc_log=None):
    """the lane spread, dd against qd, of D on the same dd star: omega_0 of X1 (both eps) and sigma_0^2 of Z (alpha_0 = 2.0, 3.0)"""
    import mpmath as mp
    mp.mp.dps = 50
    parts, worst = [], mp.mpf(0)
    for eps in (1, -1):
        x0 = mp.mpf(stored_crossing()[eps]["σ² = 0"])

        def d0(w, lane):
            ls = rust("row", omega=mp.nstr(w, 40), eps=eps, s2list="0", dump=1, lane=lane)
            return mp.mpf([l for l in ls if l["cmd"] == "D"][0]["D"])

        w_dd, _ = secant_root(lambda w: d0(w, "dd"), x0, x0 + mp.mpf("1e-6"), mp.mpf("1e-30"))
        w_qd, _ = secant_root(lambda w: d0(w, "qd"), x0, x0 + mp.mpf("1e-6"), mp.mpf("1e-30"))
        sp = abs(w_dd - w_qd)
        worst = max(worst, sp)
        parts.append(f"omega_0(eps={eps:+d}) dd {mp.nstr(w_dd, 24)}, qd {mp.nstr(w_qd, 24)}, difference {float(sp):.1e}")
    lines = arc_lines(arc_log)
    for al in Z_TARGETS:
        key = round(al, 6)
        base = [l for l in lines if l["cmd"] == "mroots" and l["grid"] == "base" and round(l["alpha_c"], 6) == key]
        var = [l for l in lines if l["cmd"] == "zvariant" and l.get("variant") == "qd lane" and round(l["alpha_c"], 6) == key]
        if not base or not var or not var[0].get("roots"):
            worst = mp.mpf(1)
            parts.append(f"sigma_0^2(alpha_0={al}): not computed in the qd lane")
            continue
        r_dd, r_qd = mp.mpf(base[0]["roots"][0]), mp.mpf(var[0]["roots"][0])
        sp = abs(r_qd / r_dd - 1)
        worst = max(worst, sp)
        parts.append(f"sigma_0^2(alpha_0={al}) dd {mp.nstr(r_dd, 24)}, qd {mp.nstr(r_qd, 24)}, relative difference {float(sp):.1e}")
    report("E2", worst <= mp.mpf(str(CRIT["Z_QD"])), "lane spread, dd against qd, of D on the same dd star (the qd lane of the integrator only; omega_N and omega_M in qd need qd stars at each omega: not done), worst "
           f"{float(worst):.1e} (limit {CRIT['Z_QD']:g}); " + "; ".join(parts))


def check_supplement(supp_log):
    """Amendment A3: the scans at alpha_0 = 3.5 and 3.9 on the stars as Newton converged them, outside A1.8's gate; enters no verdict"""
    ls = [json.loads(x) for x in open(supp_log) if x.startswith("{")]
    out = []
    for l in ls:
        if l["cmd"] == "mroots" and l["grid"] == "base":
            out.append(f"alpha_0={l['alpha_c']}: star residual {l['star_fnorm']:.1e} (A1.8 gate 1e-28), {l['sign_changes']} sign change(s) of D on the grid, roots {[x[:14] for x in l['roots']] or 'none'}")
        if l["cmd"] == "mroots" and l["grid"] == "dense":
            out.append(f"  4x denser grid: {l['sign_changes']} sign change(s), roots {[x[:14] for x in l['roots']] or 'none'}")
        if l["cmd"] == "mroots_error":
            out.append(f"alpha_0={l['alpha_c']}: error {l['err'][:120]}")
    print("\n".join(out), flush=True)
    report("F_supp", None, "supplementary, outside A1.8's gate (Amendment A3); no verdict, F at 3.5 and 3.9 stays NOT COMPUTED; " + " | ".join(out))


def check_X1(eps, stored_key):
    import mpmath as mp
    mp.mp.dps = 50
    stored = stored_crossing()[eps]
    w_stored = {k: mp.mpf(v) for k, v in stored.items()}
    # omega0: D(sigma^2 = 0; omega) = 0 by the secant method on stars solved at each omega
    def d0(w):
        lines = rust("row", omega=mp.nstr(w, 40), eps=eps, s2list="0", dump=1)
        return mp.mpf([l for l in lines if l["cmd"] == "D"][0]["D"])
    x0 = w_stored["σ² = 0"]
    w0, (xp, fp, fl) = secant_root(d0, x0, x0 + mp.mpf("1e-6"), mp.mpf("1e-30"))
    slope = (fl - fp) / (w0 - xp) if w0 != xp else mp.mpf(1)
    # E3: the same D(0) at omega0 with the far radius x 1.5, and with order 24 at tol 1e-22
    var = {}
    for name, kw in (("R x1.5", dict(rfar_scale=1.5)), ("order 24, tol 1e-22", dict(order=24, tol="1e-22"))):
        lines = rust("row", omega=mp.nstr(w0, 40), eps=eps, s2list="0", dump=1, **kw)
        var[name] = abs(mp.mpf([l for l in lines if l["cmd"] == "D"][0]["D"]) / slope)
    # omega_N, omega_M: the extrema of the interpolating polynomials through stars around omega0
    res = {}
    for name, h, n in (("E: window 1e-3, degree 8", mp.mpf("2.5e-4"), 9), ("E1: window 2e-3, degree 10", mp.mpf("4e-4"), 11)):
        rows = omega_window(eps, w0, h, n)
        ws = [r[5] for r in rows]
        wN, _, _ = poly_extremum(ws, [r[2] for r in rows])
        wM, _, _ = poly_extremum(ws, [r[1] for r in rows])
        res[name] = (wN, wM)
    wN, wM = res["E: window 1e-3, degree 8"]
    e1 = abs(res["E1: window 2e-3, degree 10"][0] - wN)
    d_sn, d_nm = abs(w0 - wN), abs(wN - wM)
    ok1 = d_sn <= CRIT["X1_TOL"] and d_nm <= CRIT["X1_TOL"]
    devs = {k: abs(v - m) for k, (v, m) in {"sigma^2 = 0": (w0, w_stored["σ² = 0"]), "dN/dw = 0": (wN, w_stored["dN/dω = 0"]), "dM/dw = 0": (wM, w_stored["dM/dω = 0"])}.items()}
    st_ok = all(d <= mp.mpf(CRIT["X1_STORED"]) for d in devs.values())
    report(f"X1a_eps{eps:+d}", ok1,
           f"sigma^2 = 0 at w = {mp.nstr(w0, 22)}; dN/dw = 0 at {mp.nstr(wN, 22)}; dM/dw = 0 at {mp.nstr(wM, 22)}; "
           f"|w0 - wN| = {float(d_sn):.1e}, |wN - wM| = {float(d_nm):.1e} (limit {CRIT['X1_TOL']:g}); "
           f"error bars: E1 (fit window and degree) {float(e1):.1e}, E3 shift of w0 from R x1.5 {float(var['R x1.5']):.1e}, from order 24 tol 1e-22 {float(var['order 24, tol 1e-22']):.1e}")
    report(f"X1b_eps{eps:+d}", st_ok,
           f"against the stored table ({stored['σ² = 0']}, {stored['dN/dω = 0']}, {stored['dM/dω = 0']}) within {CRIT['X1_STORED']:g}: "
           + ", ".join(f"{k} off by {float(d):.2e}" for k, d in devs.items()))


# ---------------------------------------------------------------------------------------------------------------
# the phantom family: C (rows), F (modes), Z (the lowest root at alpha_0 = 2.0 and 3.0), the omega fold
# ---------------------------------------------------------------------------------------------------------------
FAMILY_TARGETS = [1.10, 1.15, 1.17, 1.20, 1.25, 1.30] + [round(1.4 + 0.1 * i, 2) for i in range(27)] + [4.015, 4.017906]
MODES_TARGETS = [1.10, 1.20, 1.50, 2.00, 2.50, 3.00, 3.50, 3.90, 4.00, 4.015]
Z_TARGETS = [2.00, 3.00]
ARC_MAX = 2.0  # Amendment A2: arclength to here; Newton stepping in alpha_c above
fmt_list = lambda xs: ",".join(repr(float(x)) for x in xs)  # not :g, which keeps 6 digits (4.017906 -> 4.01791)
# each stage is one child run of at most 3600 s: `arc` (multisolver arclength through the fold) or `step` (Newton, alpha_c stepped with halving)
STAGES = {
    "arc": dict(cmd="arc", alphas=[x for x in FAMILY_TARGETS if x <= ARC_MAX], modes=[1.10, 1.20, 1.50, 2.00], zvariants=[2.00],
                extra=dict(alpha_stop=ARC_MAX, chunk=10, ds=0.05, ds_max=5)),
    "arcm": dict(cmd="arc", alphas=[1.10, 1.15, 1.17, 1.20, 1.25, 1.30], modes=[1.10], zvariants=[],  # the family controls (A4): alpha_c <= 1.3, one scan
                 extra=dict(alpha_stop=1.30, chunk=10, ds=0.05, ds_max=5)),
    "stepA": dict(cmd="step", alphas=[x for x in FAMILY_TARGETS if x > ARC_MAX], modes=[], zvariants=[], extra=dict(hmax=0.05, dlnk=0.25)),
    "stepB": dict(cmd="step", alphas=[2.50, 3.00], modes=[2.50, 3.00], zvariants=[3.00], extra=dict(hmax=0.05, dlnk=0.25)),
    "stepS": dict(cmd="step", alphas=[3.50, 3.90], modes=[3.50, 3.90], zvariants=[], extra=dict(hmax=0.05, dlnk=0.25, ungated=1)),  # A3, supplementary
    "stepC": dict(cmd="step", alphas=[3.50, 3.90], modes=[3.50, 3.90], zvariants=[], extra=dict(hmax=0.05, dlnk=0.25)),  # 4.0, 4.015: not computed (A2, S1)
}
ARC_LINES = None
TAG = None  # --tag=NAME: the jsonl of a stage is RAW/NAME.jsonl (development runs)


def run_stage(name):
    """one stage: the child's JSON lines (also in RAW/<name>[_MUTATE].jsonl as they arrive)"""
    st = STAGES[name]
    kw = dict(eps=-1, omega=0.66, omega_to=0.45, dw=0.003, alphas=fmt_list(st["alphas"]), verbose=1, **st["extra"])
    if st["modes"]:
        kw["modes"] = fmt_list(st["modes"])
    if st["zvariants"]:
        kw["zvariants"] = fmt_list(st["zvariants"])
    lines = rust(st["cmd"], timeout=3600, tag=f"{TAG or name}{'_' + MUTATE if MUTATE else ''}", **kw)
    kind = st["cmd"]
    for l in lines:
        l["_kind"] = kind
    return lines


def arc_lines(arc_log=None):
    """the JSON lines of the family runs: from the logs (comma-separated JSONL files of the stages), or by running every stage"""
    global ARC_LINES
    if ARC_LINES is None:
        ARC_LINES = []
        if arc_log:
            for f in arc_log.split(","):
                ls = [json.loads(x) for x in open(f) if x.startswith("{")]
                kind = "arc" if any(l["cmd"] == "arc_start" for l in ls) else "step"
                for l in ls:
                    l["_kind"] = kind
                ARC_LINES += [l for l in ls if not l.get("ungated")]
        else:
            for name in STAGES:
                ARC_LINES += [l for l in run_stage(name) if not l.get("ungated")]
    return ARC_LINES


def check_C(arc_log=None, alpha_max=None):
    """alpha_max: the family controls (A4) compare only the stored rows up to it; the verdict is then `Cs`"""
    lines = arc_lines(arc_log)
    rows = {round(l["alpha_c"], 6): l for l in lines if l["cmd"] == "row"}
    stored = stored_family()
    if alpha_max:
        stored = {a: v for a, v in stored.items() if a <= alpha_max + 1e-9}
    nbad, ntok, out = 0, 0, []
    nmiss, halves = 0, []
    worst_rel = 0.0
    for al, st in sorted(stored.items()):
        r = rows.get(round(al, 6))
        if r is None:
            nbad += 3
            ntok += 3
            nmiss += 3
            stall = [l for l in lines if l["cmd"] == "step_error"]
            why = (f"A2: Newton stepping in alpha_c stalled at alpha_c {stall[0]['last_alpha']:.7f}, no other method tried" if stall and al > stall[0]["last_alpha"]
                   else "the continuation did not reach it")
            out.append(f"  alpha_c={al}: NOT COMPUTED, no row ({why})")
            continue
        res = []
        for name, tok, val in (("w", st["w"], r["omega"]), ("N", st["N"], r["N"]), ("M", st["M"], r["M"])):
            good, fmt, half = tok_ok(tok, val)
            ntok += 1
            nbad += not good
            if not good:
                halves.append(half)
            rel = abs(float(Decimal(val) / Decimal(tok) - 1))
            worst_rel = max(worst_rel, rel)
            res.append(f"{name} {tok}{'=' if good else '!='}{fmt} ({half:.2f}, rel {rel:.1e})")
        res.append(f"max|2m/r| {st['max2mr']} / {r['max2mr']:.4f}, rs {st['rs']} / {r['rs']:.2f}")
        out.append(f"  alpha_c={al}: " + "; ".join(res))
    print("\n".join(out), flush=True)
    chunks = [l for l in lines if l["cmd"] == "arc_chunk"]
    tr = [l for l in lines if l["cmd"] == "arc_track"][0]
    err = [l for l in lines if l["cmd"] in ("arc_error", "step_error")]
    name1 = "Cs" if alpha_max else "C1"
    n_arc = len({round(l["alpha_c"], 6) for l in lines if l["cmd"] == "row" and l["_kind"] == "arc"})
    n_step = len({round(l["alpha_c"], 6) for l in lines if l["cmd"] == "row" and l["_kind"] == "step"})
    lo, hi = (min(halves), max(halves)) if halves else (0.0, 0.0)
    ends = [l for l in lines if l["cmd"] == "step_end"]
    report(name1, nbad == 0 and not err, f"{ntok - nbad} of {ntok} tokens (w, N, M of {len(stored)} stored rows) equal their stored print; largest relative difference {worst_rel:.1e}; "
           f"{len(halves)} tokens differ, by {lo:.2f} to {hi:.2f} half-units of the last printed digit (the stored rows are boundary-value solutions at 1e-8), {nmiss} are NOT COMPUTED (3 rows without a result); "
           f"{n_arc} distinct rows beside the arclength path (to alpha_c = {max(c['alpha_end'] for c in chunks):.4f} in {len(chunks)} chunks, {sum(c['points'] for c in chunks)} points, "
           f"{sum(c['refused'] for c in chunks)} refused steps), {n_step} distinct rows by Newton stepping in alpha_c (over the {len(ends)} stepping stage runs: {sum(e['solves'] for e in ends)} solves, "
           f"{sum(e['refused'] for e in ends)} refused; rows 2.5, 3.0, 3.5, 3.9 are computed in two stage runs each); "
           f"omega_min on the path {tr['omega_min_point'][1]:.8f} at alpha_c {tr['omega_min_point'][0]:.5f} (path point, not the fold)" + (f"; ERROR {err[0]['err'][:200]}" if err else ""))
    reached = max(c["alpha_end"] for c in chunks)
    if alpha_max:
        return
    report("C2", reached >= 4.017906 and not any(l["cmd"] == "arc_error" for l in lines),
           f"the pseudo-arclength continuation as worded in the README (the weak-field star up to alpha_c = 4.017906): the arclength path reached alpha_c = {reached:.4f}; "
           f"Amendment A2 ends it at {ARC_MAX} (cost: dev_02, dev_05) and takes the rows above by Newton stepping; as worded: NOT MET")


def check_F(arc_log=None, alpha_max=None):
    """alpha_max: the family controls (A4) compare only the stored alpha_0 up to it; the verdict is then `Fs`"""
    lines = arc_lines(arc_log)
    stored = stored_modes()
    if alpha_max:
        stored = {a: v for a, v in stored.items() if a <= alpha_max + 1e-9}
    fam = stored_family()
    full = stored_modes_full()
    rows = {round(l["alpha_c"], 6): l for l in lines if l["cmd"] == "row"}
    base = {round(l["alpha_c"], 6): l for l in lines if l["cmd"] == "mroots" and l["grid"] == "base"}
    dense = {round(l["alpha_c"], 6): l for l in lines if l["cmd"] == "mroots" and l["grid"] == "dense"}
    nbad, ntok, out = 0, 0, []
    for al, (wtok, thrtok, rtoks) in sorted(stored.items()):
        key = round(al, 6)
        if key not in base or key not in rows:
            nbad += 1
            ntok += 1
            gate = [l for l in lines if l["cmd"] == "mroots_error" and round(l["alpha_c"], 6) == key]
            stall = [l for l in lines if l["cmd"] == "step_error"]
            if gate and key in rows:
                why = f"the star's Newton residual {rows[key]['fnorm']:.2e} is above A1.8's gate 1e-28; A3's supplementary run is outside the gate and gives no verdict"
            elif stall and al > stall[0]["last_alpha"]:
                why = f"A2: Newton stepping in alpha_c stalled at alpha_c {stall[0]['last_alpha']:.7f}"
            else:
                why = "no scan in the logs"
            out.append(f"  alpha_0={al}: NOT COMPUTED ({why})")
            continue
        got = [Decimal(x) for x in base[key]["roots"]]
        res = []
        good, fmt, half = tok_ok(wtok, rows[key]["omega"])
        ntok += 1
        nbad += not good
        res.append(f"w {wtok}{'=' if good else '!='}{fmt}")
        good, fmt, half = tok_ok(thrtok, base[key]["thr"])
        ntok += 1
        nbad += not good
        res.append(f"thr {thrtok}{'=' if good else '!='}{fmt}")
        if not rtoks:
            good = len(got) == 0
            ntok += 1
            nbad += not good
            res.append(f"no root: crate {'none' if good else [str(x)[:12] for x in got]} (sign changes {base[key]['sign_changes']}, dense grid {dense[key]['sign_changes'] if key in dense else '?'})")
        else:
            for i, tk in enumerate(rtoks):
                ntok += 1
                if i < len(got):
                    good, fmt, half = tok_ok(tk, got[i])
                else:
                    good, fmt, half = False, "missing", float("inf")
                nbad += not good
                res.append(f"root {tk}{'=' if good else '!='}{fmt} ({half:.2f})")
            extra = len(got) - len(rtoks)
            if extra > 0:
                res.append(f"{extra} further root(s) {[str(x)[:12] for x in got[len(rtoks):]]}")
            if any(x < 0 for x in got):
                nbad += 1
                res.append("NEGATIVE root")
            fr = full.get(round(al, 4))
            if fr:
                res.append("vs stored double: " + ", ".join(f"{float(g) / s - 1:.1e}" for g, s in zip(got, fr)))
        if key in dense:
            same = len(dense[key]["roots"]) == len(base[key]["roots"]) and all(abs(Decimal(a) - Decimal(b)) <= Decimal("1e-20") * abs(Decimal(a)) + Decimal("1e-25") for a, b in zip(dense[key]["roots"], base[key]["roots"]))
            res.append(f"dense grid: {'same roots' if same else 'DIFFERENT ' + str(dense[key]['roots'])}")
            nbad += not same
            ntok += 1
        out.append(f"  alpha_0={al}: " + "; ".join(res))
    print("\n".join(out), flush=True)
    report("Fs" if alpha_max else "F", nbad == 0, f"{ntok - nbad} of {ntok} tokens of the phantom modes table (w, (1-w)^2, roots on [-0.2, (1-w)^2), 'none') equal their stored print or hold on the 4x denser grid")


def check_Z(arc_log=None):
    lines = arc_lines(arc_log)
    full = stored_modes_full()
    out, ok_all = [], True
    for al in Z_TARGETS:
        key = round(al, 6)
        base = [l for l in lines if l["cmd"] == "mroots" and l["grid"] == "base" and round(l["alpha_c"], 6) == key]
        var = [l for l in lines if l["cmd"] == "zvariant" and round(l["alpha_c"], 6) == key]
        if not base:
            ok_all = False
            out.append(f"alpha_0={al}: not computed")
            continue
        r0 = Decimal(base[0]["roots"][0])
        parts = [f"alpha_0 = {al}: sigma_0^2 = {str(r0)[:34]}"]
        spreads = {}
        for v in var:
            if v["roots"]:
                spreads[v["variant"]] = abs(Decimal(v["roots"][0]) / r0 - 1) if v["variant"] != "f64 lane" else abs(float(v["roots"][0]) / float(r0) - 1)
        for k_, v in spreads.items():
            parts.append(f"{k_}: {float(v):.1e}")
        fr = full.get(round(al, 4))
        if fr:
            parts.append(f"stored double {fr[0]:.15e} (rel {float(r0) / fr[0] - 1:.1e}); phantom_fd.py's extrapolated FD was {'6.0e-3' if al == 2.0 else '1.04e-2'} off")
        out.append("; ".join(parts))
    print("\n".join(out), flush=True)
    report("Z", ok_all, "sigma_0^2 at alpha_0 = 2.0 and 3.0 in double-double, with the spread over R x1.5, order 24 at tol 1e-22, a second node layout and the f64 lane (an estimator range over choices of the method, not an enclosure); the third declared error bar, the independent multishoot eigenvalue solver, is NOT DONE (E2, dd against qd, is its own check)")


def check_fold():
    import mpmath as mp
    mp.mp.dps = 50
    lines = rust("fold", timeout=3600, eps=-1, omega=0.7603, omega_to=0.45, dw=0.003, ds=0.01, ds_max=0.04, chunk=40, max_chunks=8)
    f = [l for l in lines if l["cmd"] == "fold"]
    if not f:
        report("FOLD", False, str([l for l in lines if l["cmd"] == "fold_error"]))
        return
    f = f[0]
    det = f["detected_between_k"]
    # omega_min and its alpha_c from the family in alpha_c: the extremum of the interpolating polynomial of omega(alpha_c)
    res = {}
    for name, h, n in (("window 1.6e-2, degree 8", mp.mpf("2e-3"), 9), ("window 3e-2, degree 10", mp.mpf("3e-3"), 11)):
        alphas = [mp.mpf("1.172") + (i - n // 2) * h for i in range(n)]
        al_f = [float(a) for a in alphas]
        rows = rust("follow", timeout=3600, eps=-1, omega=0.66, omega_to=0.45, dw=0.003, h=0.01, alphas=",".join(repr(a) for a in al_f))
        got = {round(r["alpha_c"], 12): r for r in rows if r["cmd"] == "row"}
        xs = [mp.mpf(repr(a)) for a in al_f]
        ys = [mp.mpf(got[round(a, 12)]["omega"]) for a in al_f]
        astar, c, mid = poly_extremum(xs, ys)
        wmin = sum(c[j] * (astar - mid) ** j for j in range(len(c)))
        res[name] = (astar, wmin)
    (a1, w1), (a2, w2) = res["window 1.6e-2, degree 8"], res["window 3e-2, degree 10"]
    ok = det is not None and abs(float(w1) - 0.42226) < 1e-5 and abs(float(a1) - 1.172) < 1e-3
    loc = f["omega_min"] if f["omega_min"] else f"not located: {f['locate_error']}"
    report("FOLD", ok, f"multisolver's arclength in k detects the fold of k (minimum of omega) between k = {det[0]:.10f} and {det[1]:.10f} ({f['points']} path points); locate_fold: {loc}; "
           f"omega(alpha_c) polynomial: omega_min = {mp.nstr(w1, 20)} at alpha_c = {mp.nstr(a1, 14)} (window 1.6e-2, degree 8), spread to the 3e-2 window, degree 10: omega_min {float(abs(w1 - w2)):.1e}, alpha_c {float(abs(a1 - a2)):.1e}; stored: omega_min = 0.42226 at alpha_0 = 1.172 (phantom_checks.py C4)")


CHECKS = ["P1", "P2", "P3", "G1", "G2", "S1", "S2", "S3", "S5", "A", "X1", "C", "F", "Z", "E2", "FOLD"]


def main():
    global MUTATE, RAW, TAG
    argv = sys.argv[1:]
    only = None
    arc_log = None
    stage = None
    supp_log = None
    p_lanes = ("f64", "dd", "qd")
    for a in argv:
        if a.startswith("--mutate="):
            MUTATE = a.split("=", 1)[1]
        elif a.startswith("--only="):
            only = a.split("=", 1)[1].split(",")
        elif a.startswith("--raw="):
            RAW = a.split("=", 1)[1]
        elif a.startswith("--arc-log="):
            arc_log = a.split("=", 1)[1]
        elif a.startswith("--p-lanes="):  # the lanes of the planted checks P1 and P2 (the controls noqr and lowsum never finish in dd and qd)
            p_lanes = tuple(a.split("=", 1)[1].split(","))
        elif a.startswith("--supp-log="):
            supp_log = a.split("=", 1)[1]
        elif a.startswith("--stage="):
            stage = a.split("=", 1)[1]
        elif a.startswith("--tag="):
            TAG = a.split("=", 1)[1]
    if "--gen" in argv:
        gen()
        return
    if "--build" in argv:
        build()
        return
    if "--list" in argv:
        print("checks:", ", ".join(CHECKS))
        return
    if not os.path.exists(BIN):
        build()
    t0 = time.time()
    if stage:
        lines = run_stage(stage)
        bad = [l for l in lines if l["cmd"] in ("arc_error", "step_error")]
        print(f"stage {stage}: {len(lines)} lines, {sum(1 for l in lines if l['cmd'] == 'row')} rows, {len(bad)} errors; mutation: {MUTATE}; {time.time() - t0:.0f} s")
        sys.exit(1 if bad else 0)
    table = {
        "P1": lambda: check_P1(p_lanes), "P2": lambda: check_P2(p_lanes), "P3": check_P3, "G1": check_G1, "G2": check_G2, "S1": check_S1,
        "S2": lambda: report("S2", float([l for l in rust("s2", omega="0.75", eps=1, x0=4.5) if l["cmd"] == "S2"][0]["worst_rel"]) <= CRIT["S2_TOL"],
                             "multishoot's recursion against the crate's own Num-lane recursion (star_ref.rs) at the node r = 4.5 of the canonical star"),
        "S3": check_S3, "S5": check_S5, "A": check_A, "E2": lambda: check_E2(arc_log),
        "X1": lambda: (check_X1(1, "p"), check_X1(-1, "m")),
        "C": lambda: check_C(arc_log), "F": lambda: check_F(arc_log), "Cs": lambda: check_C(arc_log, alpha_max=1.30), "Fs": lambda: check_F(arc_log, alpha_max=1.10), "Z": lambda: check_Z(arc_log), "FOLD": check_fold,
    }
    for name in (only if only else CHECKS):  # `only` may name the control checks Cs and Fs, which are not in the default list
        print(f"--- {name} ({time.time() - t0:.0f} s)", flush=True)
        try:
            table[name]()
        except Exception as e:  # a check that cannot run is a failure, not a skip
            import traceback
            traceback.print_exc()
            report(name, False, f"did not run: {e}")
    if supp_log:
        print("--- supplement (A3)", flush=True)
        check_supplement(supp_log)
    n_ok = sum(1 for _, ok in OUT_RESULTS if ok is True)
    n_bad = sum(1 for _, ok in OUT_RESULTS if ok is False)
    n_skip = sum(1 for _, ok in OUT_RESULTS if ok is None)
    print(f"{n_ok} passed, {n_bad} failed, {n_skip} skipped of {len(OUT_RESULTS)} checks; mutation: {MUTATE}; {time.time() - t0:.0f} s")
    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
