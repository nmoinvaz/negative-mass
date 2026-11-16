#!/usr/bin/env python3
"""Why phantom_follow.py stops: the Yukawa tail of the background BVP has a branch point at w = 1.

    python phantom_stop.py FOLLOW.jsonl STOP_CKPT.npz SNAP_DIR OUT.json [--mutate=fit|wkb]
    python phantom_stop.py --selftest

phantom_follow.py solves each star on R = rs + 30/k (rs: phi = 1e-3 phi_c, k = sqrt(1 - w^2)) with the tail condition
u(R) = -k - 1/R. Its rows end with w = 0.99994, k = 0.011, R = 1919, and the next solve (R = 2933) reports a singular
collocation Jacobian (scipy status 2, not the node limit). The tests below say whether that is the star or the boundary
condition. Each must hold, and the exit status is 1 when one fails.

  S1  w(alpha_c) over the rows with k < KMAX is a quadratic to FIT_TOL, and w = 1 lies within WINDOW above the last row
      (alpha_c*). The stop is where w reaches 1, not a turning point of the family.
  S2  The last star solved on R = 600 instead of its own R has the same w, N, M to TOL (the star does not feel R, so
      R = rs + 30/k growing without bound is the rule's, not the star's).
  S3  With the Yukawa condition on that fixed R = 600 the continuation in alpha_c goes on until it fails, and fails
      only within 1e-6 of the alpha_c* of S1 or beyond it (w >= 1 has no real k), after reaching k below K_END < the
      follow run's k.
  S4  A WKB tail condition u(R) = -kappa - 1/R, kappa = a sqrt(1 - w^2/alpha^2) from the star's own m(R) and nu(R),
      which has no k and no branch point at w = 1 (alpha^2 = 1 - 2m/R > 1 for m < 0), reproduces the Yukawa star at
      k = 0.089 (the alpha_c = 4.0 snapshot) to TOL, and continues through w = 1 (to w >= W_END) with the second
      differences of N and M at most SMOOTH of the first. It is a diagnostic here, not used for any quoted table.

--mutate=fit shifts the rows' w by 0.05 (the root leaves the window, S1 must fail), --mutate=wkb puts the Yukawa
condition, with its k, back into S4 (the continuation cannot pass w = 1, S4 must fail). A wrong sign of kappa does not
change the star (the BVP core does not feel its outer condition, as S2 shows), so it is not a usable mutation. --selftest runs S1 on a planted family w(a) =
1 - 0.2 (a* - a) - 0.3 (a* - a)^2 with the planted root, as planted and shifted.
"""
import json
import math
import os
import sys
import warnings

warnings.filterwarnings("ignore")
import numpy as np  # noqa: E402

KMAX = 0.15
FIT_TOL = 1e-8
WINDOW = 1e-3
TOL = 1e-7
K_END = 4e-3
W_END = 1.002
SMOOTH = 1e-2


def fit_alpha_star(rows):
    sel = [r for r in rows if r["k"] < KMAX]
    a = np.array([r["alpha_c"] for r in sel])
    w = np.array([r["w"] for r in sel])
    x = a - a[-1]
    c = np.polyfit(x, w, 2)
    resid = float(np.max(np.abs(np.polyval(c, x) - w)))
    roots = np.roots(c - np.array([0, 0, 1.0]))
    roots = [r.real for r in roots if abs(r.imag) < 1e-12 and r.real > 0]
    return (a[-1] + min(roots) if roots else None), resid, len(sel), float(a[-1])


def s1(rows, mutate=None):
    if mutate == "fit":
        rows = [dict(r, w=r["w"] - 0.05, k=r["k"]) for r in rows]
    a_star, resid, n, a_last = fit_alpha_star(rows)
    ok = bool(a_star is not None and resid <= FIT_TOL and 0 < a_star - a_last <= WINDOW)
    return ok, dict(alpha_star=None if a_star is None else float(a_star), resid=resid, rows=n, alpha_last=a_last)


def advance(st, ac, R, tail, step, a_end, mutate=None, hmin=1e-9):
    """alpha_c-continuation at fixed R. Stops at the first failure with h < hmin or at a_end."""
    import phantom_follow as F
    rows, h = [], step
    while ac < a_end - 1e-12:
        s = min(h, a_end - ac)
        try:
            st2 = F.solve_lapse(ac + s, F.EPS, st, R, tail=tail, mutate=mutate)
            if not (st2.w < 1.0 or tail == "wkb"):
                raise RuntimeError("w >= 1")
        except RuntimeError as e:
            h = s / 2
            if h < hmin:
                return st, ac, rows, str(e)
            continue
        ac += s
        st = st2
        h = min(h * 1.3, step)
        rows.append(dict(alpha_c=ac, w=st.w, N=st.N, M=st.M, k=math.sqrt(max(0.0, 1.0 - st.w ** 2))))
    return st, ac, rows, None


def rel(a, b):
    return abs(a - b) / abs(b)


def run(follow, ckpt, snap_dir, out, mutate):
    import phantom_follow as F
    rows = [json.loads(ln) for ln in open(follow)]
    stop = [r for r in rows if r.get("stop")][0]
    rows = [r for r in rows if not r.get("stop")]
    res, ok = {}, True

    ok1, d = s1(rows, mutate if mutate == "fit" else None)
    res["S1"] = dict(ok=ok1, **d)
    print(f"S1 {'PASS' if ok1 else 'FAIL'}: w(alpha_c) quadratic over {d['rows']} rows (k < {KMAX}), residual {d['resid']:.2e}, "
          f"w = 1 at alpha_c* = {d['alpha_star']}, last row alpha_c = {d['alpha_last']:.6f}, "
          f"follow stop record: status {stop['status']}, R tried {stop['message'].split('R=')[1].split(':')[0]}", flush=True)
    ok &= ok1

    if mutate == "fit":                       # the mutation targets S1 only
        json.dump(dict(mutate=mutate, tests=res), open(out, "w"), indent=1)
        return ok
    g, ac = F.load_checkpoint(ckpt)
    st_own = F.solve_lapse(ac, F.EPS, g, g.R)
    st600 = F.solve_lapse(ac, F.EPS, g, 600.0)
    d2 = {k: rel(getattr(st600, k), getattr(st_own, k)) for k in ("w", "N", "M")}
    ok2 = bool(max(d2.values()) <= TOL)
    res["S2"] = dict(ok=ok2, rel=d2, R_own=g.R, k=math.sqrt(1 - st_own.w ** 2))
    print(f"S2 {'PASS' if ok2 else 'FAIL'}: star at alpha_c = {ac:.6f} (k = {res['S2']['k']:.4f}) on R = {g.R:.0f} and R = 600: "
          f"relative differences {d2}", flush=True)
    ok &= ok2

    if mutate == "wkb":                       # the mutation targets S4 only
        ok2 = ok3 = True
        res["S3"] = {"ok": True, "skipped": True}
        a_star = d["alpha_star"]
    else:
        st3, ac3, rows3, err = advance(st600, ac, 600.0, "yukawa", 1e-4, ac + 0.01)
        a_star = d["alpha_star"]
        k3 = math.sqrt(1 - st3.w ** 2)
        ok3 = bool(k3 <= K_END and k3 < stop["k"] and a_star is not None and abs(ac3 - a_star) <= 1e-6 + 1e-3 * (a_star - ac))
        res["S3"] = dict(ok=ok3, alpha_end=ac3, k_end=k3, w_end=st3.w, N_end=st3.N, M_end=st3.M, last_error=err,
                         steps=len(rows3), follow_k_stop=stop["k"])
        print(f"S3 {'PASS' if ok3 else 'FAIL'}: Yukawa tail on R = 600 reaches alpha_c = {ac3:.8f}, w = {st3.w:.8f}, k = {k3:.2e} "
              f"(follow stopped at k = {stop['k']:.4f}), N = {st3.N:.4f}, M = {st3.M:+.4f}; first failure: {err}", flush=True)
        ok &= ok3

    g4, a4 = F.load_checkpoint(os.path.join(snap_dir, "star_a4.0000.npz"))
    yuk = F.solve_lapse(a4, F.EPS, g4, g4.R)
    Rw = 280.0                              # rs + 100
    tail = "yukawa" if mutate == "wkb" else "wkb"       # the mutation puts k back into the tail condition
    wk = F.solve_lapse(a4, F.EPS, g4, Rw, tail=tail)
    d4 = {k: rel(getattr(wk, k), getattr(yuk, k)) for k in ("w", "N", "M")}
    ok4a = bool(max(d4.values()) <= TOL)
    print(f"S4a {'PASS' if ok4a else 'FAIL'}: WKB tail on R = {Rw:.0f} against the Yukawa star at alpha_c = 4.0 (k = 0.089): {d4}",
          flush=True)
    st5, ac5, rows5, err5 = advance(wk, a4, Rw, tail, 1e-3, 4.03, hmin=1e-6)
    all5 = [dict(alpha_c=a4, w=wk.w, N=wk.N, M=wk.M)] + rows5
    N = np.array([r["N"] for r in all5]); M = np.array([r["M"] for r in all5]); W = np.array([r["w"] for r in all5])
    sm = max(float(np.max(np.abs(np.diff(N, 2)) / np.abs(np.diff(N))[1:])), float(np.max(np.abs(np.diff(M, 2)) / np.abs(np.diff(M))[1:]))) \
        if len(all5) > 4 else 1.0
    ok4b = bool(len(all5) > 4 and W[-1] >= W_END and sm <= SMOOTH and np.all(np.diff(N) > 0))
    n1 = float(np.interp(1.0, W, N)) if len(all5) > 4 and W[-1] > 1 else None
    m1 = float(np.interp(1.0, W, M)) if len(all5) > 4 and W[-1] > 1 else None
    res["S4"] = dict(ok=bool(ok4a and ok4b), rel_vs_yukawa=d4, alpha_end=ac5, w_end=float(W[-1]), smooth=sm, N_at_w1=n1,
                     M_at_w1=m1, rows=all5, last_error=err5)
    print(f"S4b {'PASS' if ok4b else 'FAIL'}: WKB tail continues through w = 1 to alpha_c = {ac5:.4f}, w = {W[-1]:.6f} "
          f"({len(all5)} stars), N increasing, second/first differences at most {sm:.2e}; at w = 1 N = {n1}, M = {m1}", flush=True)
    ok &= ok4a and ok4b
    json.dump(dict(mutate=mutate, tests=res), open(out, "w"), indent=1)
    return ok


def selftest():
    a_star = 4.0182
    rows = []
    for i in range(30):
        a = a_star - 0.0178 + 0.0006 * i
        w = 1 - 0.2 * (a_star - a) - 0.3 * (a_star - a) ** 2
        rows.append(dict(alpha_c=a, w=w, k=math.sqrt(1 - w * w)))
    ok, d = s1(rows)
    exact = ok and abs(d["alpha_star"] - a_star) < 1e-9
    print(f"selftest planted alpha_c* = {a_star}: {'recovered ' + format(d['alpha_star'], '.10f') if exact else 'NOT recovered ' + str(d)}")
    ok2, d2 = s1(rows, "fit")
    print(f"selftest --mutate=fit: S1 {'FIRED' if not ok2 else 'DID NOT FIRE'}")
    rows2 = [dict(r) for r in rows]
    rows2[10]["w"] += 1e-6
    ok3, _ = s1(rows2)
    print(f"selftest a planted 1e-6 kink in one row: S1 {'FIRED' if not ok3 else 'DID NOT FIRE'}")
    return exact and not ok2 and not ok3


def main():
    argv = sys.argv[1:]
    mutate = next((a.split("=", 1)[1] for a in argv if a.startswith("--mutate=")), None)
    if "--selftest" in argv:
        sys.exit(0 if selftest() else 1)
    pos = [a for a in argv if not a.startswith("--")]
    sys.exit(0 if run(pos[0], pos[1], pos[2], pos[3], mutate) else 1)


if __name__ == "__main__":
    main()
