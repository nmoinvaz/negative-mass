#!/usr/bin/env python3
"""Ghost-sector (test-field) decay rates mu of Choptuik's DSS solution in D dimensions, on the real axis.

A second massless scalar psi = 0 in the critical solution, ghost or normal, perturbs as a test field: box psi = 0 on
the fixed background, psi ~ e^{mu tau} x (periodic in tau). Its multiple-shooting determinant det M_psi(mu) is
choptuik-large-d's hp/modecount_cx `test_det` (notes/GHOST_ADMIXTURE.md there); this script only reads that code.

    python3 ghost_mu.py --hp=<choptuik-large-d>/hp STATE.npz [--lo=-1.0] [--hi=-0.04] [--h=0.04] [--mode=1]
                        [--plant=X] [--threads=2] [--out=ghost_mu.jsonl] [--tag=...]

mode 1 is the test field (a and c tangents held). mode 2 is choptuik-large-d's control TEST_MUT (c held, a's response
to p, m kept): a wrong operator, which must lose the exact boost zero mu = 0 and move the decaying zeros.
--plant=X multiplies det by (mu - X): the finder must return X (a control of the scan and refinement).

On the real axis det M_psi is real (arg 0 or pi). Zeros are bracketed by sign changes on the grid lo, lo + h, ..., hi
and refined by the Illinois method on sign(det) |det| e^{-k mu}, k the bracket's slope of ln|det| (det grows like
e^{~300 mu}; removing the trend keeps the bracket well scaled). Exit 1 if a refinement fails to bracket or converge.
"""
import json
import math
import os
import sys
import time


def parse(argv):
    args = [a for a in argv if not a.startswith("--")]
    opts = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in argv if a.startswith("--"))
    return args, opts


def main():
    args, opts = parse(sys.argv[1:])
    hp = opts["hp"]
    sys.path.insert(0, hp)
    from modecount_common import setup  # noqa: E402  (choptuik-large-d, read only)

    mode = int(opts.get("mode", 1))
    plant = float(opts["plant"]) if "plant" in opts else None
    lo, hi, h = float(opts.get("lo", -1.0)), float(opts.get("hi", -0.04)), float(opts.get("h", 0.04))
    out = opts.get("out", os.path.join(os.path.dirname(os.path.abspath(__file__)), "ghost_mu.jsonl"))
    tol = float(opts.get("tol", 1e-11))
    st = setup(args[0], threads=int(opts.get("threads", 2)), log=lambda m: print(m, flush=True))
    print(f"{os.path.basename(args[0])}: D {st.D}, M {st.M}, K {st.K}, |F| {st.F:.1e}, Delta {st.Delta:.10f}, "
          f"omega {st.omega:.6f}, mode {mode}, plant {plant}", flush=True)
    nev = [0]
    t0 = time.time()

    def ev(x):
        """(ln|g|, sign g, min pivot) at real mu = x."""
        la, arg, pmin, _ = st.c.test_det(st.zh, st.zl, float(x), 0.0, 0.0, 0.0, mode)
        nev[0] += 1
        sg = 1.0 if math.cos(arg) > 0 else -1.0
        if abs(math.sin(arg)) > 1e-6:
            print(f"  warning: det not real at {x}: arg {arg}", flush=True)
        if plant is not None:
            d = x - plant
            la, sg = la + math.log(abs(d)) if d != 0 else -math.inf, sg * (1.0 if d > 0 else -1.0)
        return la, sg, pmin

    # the exact boost zero: smallest pivot at mu = 0 against mu = -h
    z0, zh_ = ev(0.0), ev(-h)
    print(f"  mu = 0: ln|det| {z0[0]:.4f}, min pivot {z0[2]:.2e};  mu = {-h}: min pivot {zh_[2]:.2e}", flush=True)

    n = int(round((hi - lo) / h))
    xs = [lo + i * h for i in range(n + 1)]
    vals = []
    for x in xs:
        v = ev(x)
        vals.append(v)
        print(f"  {x:+.4f}: ln|det| {v[0]:12.4f}  sign {v[1]:+.0f}  min pivot {v[2]:.1e}", flush=True)

    zeros, bad = [], []
    for i in range(n):
        (la, sa, _), (lb, sb, _) = vals[i], vals[i + 1]
        if sa == sb:
            continue
        a, b = xs[i], xs[i + 1]
        k = (lb - la) / (b - a)

        def f(x, la_=None):
            l_, s_, _ = ev(x) if la_ is None else la_
            return s_ * math.exp(l_ - la - k * (x - a))

        fa, fb = f(a, vals[i]), f(b, vals[i + 1])
        side, hist, fc = 0, [], 1.0
        for it in range(60):
            c = (a * fb - b * fa) / (fb - fa)
            fc = f(c)
            hist.append(c)
            if fc == 0.0:
                a = b = c
                break
            if fc * fb > 0:              # c replaces b; Illinois: halve the retained end on a repeat
                b, fb = c, fc
                if side == -1:
                    fa *= 0.5
                side = -1
            else:
                a, fa = c, fc
                if side == 1:
                    fb *= 0.5
                side = 1
            if abs(b - a) < tol:
                break
        root = 0.5 * (a + b)
        ok = abs(b - a) < 10 * tol or fc == 0.0
        zeros.append(dict(bracket=[xs[i], xs[i + 1]], mu=root, width=abs(b - a), iters=len(hist), ok=ok))
        print(f"  zero in [{xs[i]:+.4f}, {xs[i + 1]:+.4f}]: mu = {root:+.12f} (width {abs(b - a):.1e}, "
              f"{len(hist)} steps)", flush=True)
        if not ok:
            bad.append(root)
    rec = dict(state=os.path.basename(args[0]), path=args[0], D=st.D, M=st.M, K=st.K, nodes=list(st.meta["nodes"]),
               F=st.F, Delta=st.Delta, omega=st.omega, mode=mode, plant=plant, lo=lo, hi=hi, h=h,
               pivot0=z0[2], pivot_h=zh_[2], zeros=zeros, samples=[[x, v[0], v[1], v[2]] for x, v in zip(xs, vals)],
               evaluations=nev[0], sec=round(time.time() - t0, 1), tag=opts.get("tag", ""),
               time=time.strftime("%Y-%m-%d %H:%M:%S"))
    with open(out, "a") as fo:
        fo.write(json.dumps(rec) + "\n")
    print(f"  {nev[0]} evaluations, {time.time() - t0:.0f} s; record appended to {out}", flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
