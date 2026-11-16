#!/usr/bin/env python3
"""Convergence of the continued phantom stars: the same star at a larger R and at a tighter BVP tolerance.

    python phantom_converge.py SNAP_DIR OUT.json [--alphas 1.5,2.5,3.5,4.0] [--mutate=matching]
    python phantom_converge.py --selftest

For each snapshot star_a<alpha_c>.npz (solved on R = rs + 30/k with tol 1e-8), the star at the same alpha_c is solved
again with (a) R x 1.5, (b) tol 1e-9, (c) tol 1e-10, and (w, N, M) are compared with the snapshot. Every relative
difference must be at most TOL = 1e-7, an order below the 8e-8 agreement the existing sigma^2 crossings reached
(a limit set before the run). The exit status is 1 when one is larger.
--mutate=matching drops the Schwarzschild matching of nu(R) in the variants, and must make the check fail (the
alpha(inf) = 1 normalization, hence w, moves by O(m/R)). --selftest solves one planted star (alpha_c = 1.3, small R, from
the same code) as planted and with the mutation.
"""
import glob
import json
import os
import re
import sys
import warnings

warnings.filterwarnings("ignore")
import numpy as np  # noqa: E402

TOL = 1e-7


def variants(ac, g, mutate=None):
    """The star and its three variants. A variant whose BVP does not converge is None (counted as a failure)."""
    import phantom_follow as F
    base = F.solve_lapse(ac, F.EPS, g, g.R)
    out = {}
    for name, kw, R in (("R x 1.5", {}, g.R * 1.5), ("tol 1e-9", {"tol": 1e-9}, g.R), ("tol 1e-10", {"tol": 1e-10}, g.R)):
        try:
            out[name] = F.solve_lapse(ac, F.EPS, base, R, mutate=mutate, **kw)
        except RuntimeError as e:
            print(f"   variant {name}: {e}")
            out[name] = None
    return base, out


def compare(base, out):
    rows = {}
    for name, st in out.items():
        rows[name] = ({k: float("inf") for k in ("w", "N", "M")} if st is None else
                      {k: abs(getattr(st, k) - getattr(base, k)) / abs(getattr(base, k)) for k in ("w", "N", "M")})
    worst = max(v for r in rows.values() for v in r.values())
    return rows, worst


def main():
    argv = sys.argv[1:]
    mutate = next((x.split("=", 1)[1] for x in argv if x.startswith("--mutate=")), None)
    import phantom_follow as F
    if "--selftest" in argv:
        st0, ac = F.start_star(1.3)
        g = F.Guess(st0.sol.x, st0.sol.y, st0.sol.yp, st0.sol.p, st0.R)
        ok = True
        for mu in (None, "matching"):
            base, out = variants(ac, g, mu)
            rows, worst = compare(base, out)
            fired = worst > TOL
            print(f"selftest alpha_c={ac:.3f} mutate={mu}: worst rel {worst:.2e}, "
                  f"{'FIRED' if fired else 'quiet'} (limit {TOL:g})")
            ok &= (fired == (mu is not None))
        sys.exit(0 if ok else 1)
    pos = [a for a in argv if not a.startswith("--")]
    alphas = None
    if "--alphas" in argv:
        alphas = [float(x) for x in argv[argv.index("--alphas") + 1].split(",")]
        pos = [p for p in pos if p != argv[argv.index("--alphas") + 1]]
    res, ok = [], True
    for f in sorted(glob.glob(os.path.join(pos[0], "star_a*.npz"))):
        a = float(re.search(r"star_a([\d.]+)\.npz", f).group(1))
        if alphas and not any(abs(a - x) < 1e-9 for x in alphas):
            continue
        g, ac = F.load_checkpoint(f)
        base, out = variants(ac, g, mutate)
        rows, worst = compare(base, out)
        ok &= worst <= TOL
        res.append(dict(alpha_c=ac, R=g.R, w=base.w, N=base.N, M=base.M, rel=rows, worst=worst))
        print(f"alpha_c={ac:.4f} R={g.R:.0f}: worst rel {worst:.2e} " + " ".join(
            f"[{n}: w {r['w']:.1e} N {r['N']:.1e} M {r['M']:.1e}]" for n, r in rows.items()), flush=True)
    json.dump(dict(tol=TOL, mutate=mutate, stars=res), open(pos[1], "w"), indent=1)
    print(("PASS" if ok and res else "FAIL") + f": {len(res)} stars, limit {TOL:g}")
    sys.exit(0 if ok and res else 1)


if __name__ == "__main__":
    main()
