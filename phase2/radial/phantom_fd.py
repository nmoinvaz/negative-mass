#!/usr/bin/env python3
"""Independent finite-difference spectrum (fd_check.fd_eigs, dense QZ) of the continued phantom stars.

    python phantom_fd.py SNAP_DIR MODES.jsonl OUT.json --alphas 1.5,2.0,3.0 [--ns 250,400,550] [--mutate=eps]
    python phantom_fd.py --selftest

It shares with phantom_modes.py only the derived equations (radial_rhs_generated.py), not the method, as fd_check.py
does for the earlier stars. For each star:
  F1  the lowest FD real sigma^2 (the exact zero of the outer tail dropped), Richardson-extrapolated from the two largest
      n as a second-order scheme, agrees with the shooting root of MODES.jsonl to REL_TOL, and the FD values converge (the
      last change is smaller than the one before). The plain FD value at n = 550 is 4 to 5 percent off on these stars,
      because the dense QZ limits n and the wall is thin, so the plain difference is recorded and not judged. This criterion
      replaced "plain difference at most 3e-2 at n = 550" after the alpha_c = 1.5 selftest showed 3.9e-2 for a star whose
      root is known; it was fixed before the alpha_c = 3.0 run finished, and the plain-difference verdicts of the first
      run stay in fd.json as F1_plain.
  F2  there is no real FD eigenvalue below -1e-8 and no complex one with |sigma^2| < 0.3 and |Im| > 1e-7 at the
      smallest and the largest n (shooting sees only the real axis).
Exit 1 if one fails. --mutate=eps solves the FD problem for the canonical sign of the source (eps = +1) on the same
star, and F1 must fail. --selftest runs the alpha_c = 1.5 star, whose lowest root is the stored 1.2704e-2 of
strong_m.txt, once as it is and once with the mutation (F1 fires).
"""
import json
import os
import re
import sys
import warnings

warnings.filterwarnings("ignore")
import numpy as np  # noqa: E402
import scipy.linalg as sl  # noqa: E402

REL_TOL = 1e-2
NS = (250, 400, 550)
STORED_15 = 1.2704e-2


def spectra(P, n):
    """Lowest real FD roots (trivial zero dropped) and every finite eigenvalue of the dense pencil."""
    import fd_check as F
    cap = {}
    orig = sl.eig

    def grab(a, b, *args, **kw):
        cap["a"], cap["b"] = a, b
        return orig(a, b, *args, **kw)
    sl.eig = grab
    try:
        real = F.fd_eigs(P, n)
    finally:
        sl.eig = orig
    allv = orig(cap["a"], cap["b"], right=False)
    allv = allv[np.isfinite(allv)]
    real = real[np.abs(real) > 1e-8]
    return real, allv


def star_of(snap_dir, alpha, eps):
    import phantom_follow as F
    from radial_modes import RadialProblem
    g, ac = F.load_checkpoint(os.path.join(snap_dir, f"star_a{alpha:.4f}.npz"))
    st = F.solve_lapse(ac, F.EPS, g, g.R)
    P = RadialProblem(st)
    if eps != F.EPS:                       # the mutation, the canonical sign of the source on the same star
        P.eps = eps
    return P


def judge(fd, shoot, ns):
    """F1 from the FD values fd[n]: Richardson (order 2) on the two largest n against the shooting root, and convergence."""
    v = [fd[n] for n in ns]
    n1, n2 = ns[-2], ns[-1]
    extrap = (n2 ** 2 * v[-1] - n1 ** 2 * v[-2]) / (n2 ** 2 - n1 ** 2)
    rel = abs(extrap - shoot) / abs(shoot)
    conv = abs(v[-1] - v[-2]) < abs(v[-2] - v[-3])
    plain = {n: abs(fd[n] - shoot) / abs(shoot) for n in ns}
    info = dict(extrap=extrap, rel_extrap=rel, converging=conv, rel_plain=plain)
    try:                       # recorded, not judged: the order p of the three values and the extrapolation with it
        from scipy.optimize import brentq
        n0, v0 = ns[-3], v[-3]
        ratio = (v0 - v[-2]) / (v[-2] - v[-1])
        p = brentq(lambda q: (n0 ** -q - n1 ** -q) / (n1 ** -q - n2 ** -q) - ratio, 0.2, 6.0)
        ext_p = v[-1] - (v[-2] - v[-1]) * n2 ** -p / (n1 ** -p - n2 ** -p)
        info.update(order=p, extrap_fitted=ext_p, rel_extrap_fitted=(ext_p - shoot) / abs(shoot))
    except ValueError:
        pass
    return rel <= REL_TOL and conv, info


def check_star(P, shoot, ns):
    out = dict(w=P.w, shoot=shoot, fd={}, complex={}, negative={})
    for n in ns:
        real, allv = spectra(P, n)
        out["fd"][n] = float(real[0]) if len(real) else None
        out["negative"][n] = [float(x) for x in real[real < -1e-8][:3]]
        cx = allv[(np.abs(allv) < 0.3) & (np.abs(allv.imag) > 1e-7)]
        out["complex"][n] = [complex(x).__repr__() for x in cx[:4]]
        print(f"   n={n}: lowest real FD {out['fd'][n]}, shooting {shoot}, negative {out['negative'][n]}, "
              f"complex {out['complex'][n]}", flush=True)
    ok1, d = judge(out["fd"], shoot, ns)
    out["F1_plain"] = d["rel_plain"][ns[-1]] <= 3e-2
    ok2 = all(not out["negative"][n] and not out["complex"][n] for n in (ns[0], ns[-1]))
    out["rel"], out["F1"], out["F2"] = d, ok1, ok2
    return out


def planted():
    """Planted FD values s (1 + c / n^2) with a known limit s: judged as planted (must pass) and with the limit moved
    by 30 percent (the eps mutation moves it that far, must fire)."""
    ns = NS
    s0 = 1.2704e-2
    fd = {n: s0 * (1 + 1.1e3 / n ** 2 * 1.0) for n in ns}
    ok_a, da = judge(fd, s0, ns)
    fd_m = {n: 0.7 * s0 * (1 + 1.1e3 / n ** 2) for n in ns}
    ok_b, db = judge(fd_m, s0, ns)
    print(f"selftest planted FD sequence: {'passes' if ok_a else 'FAILS'} (rel {da['rel_extrap']:.1e}); "
          f"limit moved by 30 percent: {'FIRED' if not ok_b else 'DID NOT FIRE'} (rel {db['rel_extrap']:.1e})")
    return ok_a and not ok_b


def rejudge_log(path):
    """Re-judge the two blocks of a selftest log (n = ..: lowest real FD x, shooting y): the first must pass, the second fire."""
    import re
    blocks, cur = [], []
    for ln in open(path):
        m = re.match(r"\s+n=(\d+): lowest real FD ([\d.e+-]+), shooting ([\d.e+-]+)", ln)
        if m:
            cur.append((int(m[1]), float(m[2]), float(m[3])))
            if len(cur) == 3:
                blocks.append(cur); cur = []
    ok = True
    for i, b in enumerate(blocks):
        fd = {n: v for n, v, _ in b}
        j, d = judge(fd, b[0][2], tuple(n for n, _, _ in b))
        print(f"rejudge {path} block {i} ({'planted star' if i == 0 else 'eps mutation'}): "
              f"{'F1 holds' if j else 'F1 FIRED'}, extrapolated {d['extrap']:.6e}, rel {d['rel_extrap']:.2e}, "
              f"plain {{{', '.join(f'{n}: {v:.3f}' for n, v in d['rel_plain'].items())}}}")
        ok &= (j == (i == 0))
    return ok and len(blocks) == 2


def rejudge_json(path):
    """Re-judge F1 and F2 of an fd.json written before the F1 criterion was changed (the numbers are the same)."""
    d = json.load(open(path))
    ok = True
    for r in d["stars"]:
        ns = tuple(int(n) for n in r["fd"])
        fd = {int(n): v for n, v in r["fd"].items()}
        j, info = judge(fd, r["shoot"], ns)
        f2 = all(not r["negative"][str(n)] and not r["complex"][str(n)] for n in (ns[0], ns[-1]))
        print(f"alpha_c={r['alpha_c']}: F1 {'PASS' if j else 'FAIL'} (extrapolated {info['extrap']:.6e} against shooting "
              f"{r['shoot']:.6e}, rel {info['rel_extrap']:.2e}, plain {{{', '.join(f'{n}: {v:.3f}' for n, v in info['rel_plain'].items())}}}), "
              f"F2 {'PASS' if f2 else 'FAIL'}")
        if "order" in info:
            print(f"    fitted order p = {info['order']:.2f}, extrapolated with it {info['extrap_fitted']:.6e}, "
                  f"relative to shooting {info['rel_extrap_fitted']:+.2e}")
        ok &= j and f2
    return ok


def selftest(snap_dir, ns):
    ok = planted()
    for mu in (None, "eps"):
        P = star_of(snap_dir, 1.5, 1 if mu else -1)
        r = check_star(P, STORED_15, ns)
        fired = not r["F1"]
        print(f"selftest alpha_c=1.5 mutate={mu}: F1 {'FIRED' if fired else 'holds'} ({r['rel']})")
        ok &= (fired == (mu is not None))
        if mu is None:
            ok &= r["F2"]
    return ok


def main():
    argv = sys.argv[1:]
    mutate = next((a.split("=", 1)[1] for a in argv if a.startswith("--mutate=")), None)
    opts = {}
    pos = []
    i = 0
    while i < len(argv):
        if argv[i] in ("--alphas", "--ns", "--snap"):
            opts[argv[i]] = argv[i + 1]; i += 2
        elif argv[i].startswith("--"):
            i += 1
        else:
            pos.append(argv[i]); i += 1
    ns = tuple(int(x) for x in opts.get("--ns", ",".join(map(str, NS))).split(","))
    if "--rejudge-json" in argv:
        sys.exit(0 if rejudge_json(pos[0]) else 1)
    if "--rejudge-log" in argv:
        sys.exit(0 if rejudge_log(pos[0]) else 1)
    if "--selftest" in argv:
        sys.exit(0 if selftest(opts["--snap"], ns) else 1)
    snap_dir, modes_path, out_path = pos[:3]
    shoot = {round(r["alpha_c"], 4): r["roots"][0] for r in map(json.loads, open(modes_path)) if r["roots"]}
    res, ok = [], True
    for a in [float(x) for x in opts["--alphas"].split(",")]:
        P = star_of(snap_dir, a, 1 if mutate == "eps" else -1)
        print(f"alpha_c={a} w={P.w:.6f} R={P.R:.0f}", flush=True)
        r = check_star(P, shoot[round(a, 4)], ns)
        r["alpha_c"] = a
        ok &= r["F1"] and r["F2"]
        print(f"   F1 {'PASS' if r['F1'] else 'FAIL'} (rel {r['rel']}), F2 {'PASS' if r['F2'] else 'FAIL'}", flush=True)
        res.append(r)
    json.dump(dict(mutate=mutate, rel_tol=REL_TOL, stars=res), open(out_path, "w"), indent=1, default=str)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
