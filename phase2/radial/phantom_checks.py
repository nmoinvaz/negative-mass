#!/usr/bin/env python3
"""Checks on the continued phantom family (phantom_follow.py output) and on its radial spectra.

    python phantom_checks.py follow.jsonl [--modes modes.jsonl] [--table] [--mutate=bump|eps|fold|sigma]
    python phantom_checks.py --selftest

Claims, each of which exits 1 when it fails:
  C1  alpha_c increases along the rows and the rows reach the stop record.
  C2  N has no extremum: every difference quotient dN/dalpha_c is positive. A sign change is reported with the
      parabola vertex of the three rows around it (the turning point the stability rule needs).
  C3  The first law dM = -w dN of the signed source holds on every step, M_{i+1} - M_i = -(w_i + w_{i+1})/2 (N_{i+1} - N_i)
      to FIRST_LAW_TOL of |dM| (the trapezoid rule has an O(h^2) error, which is far smaller; the earlier fine scans held
      to about 8e-5).
  C4  w(alpha_c) has exactly one fold in the rows (the known one, between alpha_c = 1.15 and 1.20), so N(w) is
      single-valued beyond it.
  C5  (with --modes) no radial mode has sigma^2 < 0 on any listed star, with the grid and the roots recorded.

--mutate breaks the input once and the matching claim must fire: bump (a dip in N around the middle row, C2), eps (the
first law with the wrong sign of the source, C3), fold (a second dip in w, C4), sigma (a negative root added to the
first star, C5). --selftest runs planted rows with a known turning point through the same code (the vertex must be
recovered, and a monotone family must report none), then every mutation.
"""
import json
import math
import sys

FIRST_LAW_TOL = 1e-4


def load(path):
    rows, stop = [], None
    for ln in open(path):
        r = json.loads(ln)
        if r.get("stop"):
            stop = r
        else:
            rows.append(r)
    return rows, stop


def vertex(xs, ys):
    """Vertex of the parabola through three points."""
    (x0, x1, x2), (y0, y1, y2) = xs, ys
    d = (x0 - x1) * (x0 - x2) * (x1 - x2)
    a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / d
    b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / d
    return -b / (2 * a)


def extrema(rows, key="N"):
    """alpha_c of each sign change of d key / d alpha_c, from the parabola through the rows around it."""
    al = [r["alpha_c"] for r in rows]
    v = [r[key] for r in rows]
    dq = [(v[i + 1] - v[i]) / (al[i + 1] - al[i]) for i in range(len(rows) - 1)]
    out = []
    for i in range(len(dq) - 1):
        if dq[i] * dq[i + 1] < 0:
            j = min(max(i, 0), len(rows) - 3)
            out.append(vertex(al[j:j + 3], v[j:j + 3]))
    return dq, out


def claims(rows, stop, modes=None, verbose=True):
    res = []

    def add(name, ok, text):
        res.append((name, ok))
        if verbose:
            print(("PASS " if ok else "FAIL ") + name + ": " + text)

    al = [r["alpha_c"] for r in rows]
    inc = all(b > a for a, b in zip(al, al[1:]))
    add("C1 rows", inc and stop is not None,
        f"{len(rows)} rows, alpha_c {al[0]:.4f} to {al[-1]:.6f}, "
        + (f"stop at alpha_c {stop['alpha_c']:.6f}, status {stop['status']}" if stop else "no stop record"))

    dq, ext = extrema(rows, "N")
    add("C2 no extremum of N", not ext and min(dq) > 0,
        f"min dN/dalpha_c = {min(dq):.4g} at alpha_c {al[dq.index(min(dq))]:.4f}, max {max(dq):.4g}"
        + (f"; sign change near alpha_c {ext}" if ext else ""))

    worst, at = 0.0, None
    for i in range(len(rows) - 1):
        a, b = rows[i], rows[i + 1]
        dM, dN = b["M"] - a["M"], b["N"] - a["N"]
        if dM == 0:
            continue
        e = abs(dM + 0.5 * (a["w"] + b["w"]) * dN) / abs(dM)
        if e > worst:
            worst, at = e, a["alpha_c"]
    add("C3 first law dM = -w dN", worst <= FIRST_LAW_TOL, f"worst relative residual {worst:.2e} at alpha_c {at:.4f}, "
        f"limit {FIRST_LAW_TOL:g}")

    dw, folds = extrema(rows, "w")
    ok4 = len(folds) == 1 and 1.15 < folds[0] < 1.20
    add("C4 one fold of w", ok4, f"fold(s) of w at alpha_c {['%.4f' % f for f in folds]}"
        + (f", w_min about {min(r['w'] for r in rows):.5f}" if folds else ""))

    if modes is not None:
        neg = [(m["alpha_c"], x) for m in modes for x in m["roots"] if x < 0]
        add("C5 no negative sigma^2", not neg and len(modes) > 0,
            f"{len(modes)} stars, {sum(len(m['roots']) for m in modes)} roots, negative: {neg}")
    return all(ok for _, ok in res), res


def mutate_rows(rows, how):
    rows = [dict(r) for r in rows]
    mid = len(rows) // 2
    if how == "bump":                  # N dips by more than the step around the middle row
        rows[mid]["N"] -= 3 * (rows[mid + 1]["N"] - rows[mid - 1]["N"])
    if how == "eps":                   # the source's sign flipped in the first law: M of the canonical sign
        for r in rows:
            r["M"] = -r["M"]
    if how == "fold":                  # a second fold in w
        rows[mid]["w"] -= 3 * abs(rows[mid + 1]["w"] - rows[mid - 1]["w"])
    return rows


def planted(turn=True):
    """Rows of a planted family: N = N0 - c (a - a*)^2 for turn, else rising, M from dM = -w dN, w(a) with one fold."""
    a_star, rows = 2.7183, []
    N = M = 0.0
    prev = None
    for i in range(0, 241):
        a = 1.1 + 0.01 * i
        w = 0.42 + 0.3 * (a - 1.17) ** 2 / (1 + 0.5 * (a - 1.17))
        n = 100 + (-50 * (a - a_star) ** 2 if turn else 80 * (a - 1.1))
        if prev is not None:
            M -= 0.5 * (w + prev[0]) * (n - prev[1])
        prev = (w, n)
        rows.append({"alpha_c": round(a, 6), "w": w, "N": n, "M": M})
    return rows, a_star


def selftest():
    ok = True
    rows, a_star = planted(turn=True)
    _, ext = extrema(rows, "N")
    hit = len(ext) == 1 and abs(ext[0] - a_star) < 1e-9
    print(f"selftest planted turning point {a_star}: {'recovered ' + format(ext[0], '.9f') if hit else 'NOT recovered ' + str(ext)}")
    ok &= hit
    rows, _ = planted(turn=False)
    stop = {"alpha_c": rows[-1]["alpha_c"], "status": 2}
    good, _ = claims(rows, stop, verbose=False)
    print(f"selftest planted monotone family passes C1 to C4: {good}")
    ok &= good
    for how, claim in (("bump", "C2"), ("eps", "C3"), ("fold", "C4")):
        _, res = claims(mutate_rows(rows, how), stop, verbose=False)
        fired = [n for n, o in res if not o and n.startswith(claim)]
        print(f"selftest --mutate={how}: {claim} {'FIRED' if fired else 'DID NOT FIRE'}")
        ok &= bool(fired)
    modes = [{"alpha_c": 2.0, "roots": [1e-2, 5e-2]}]
    _, res = claims(rows, stop, modes, verbose=False)
    ok &= all(o for _, o in res)
    modes[0]["roots"].insert(0, -1e-3)
    _, res = claims(rows, stop, modes, verbose=False)
    fired = [n for n, o in res if not o and n.startswith("C5")]
    print(f"selftest --mutate=sigma: C5 {'FIRED' if fired else 'DID NOT FIRE'}")
    ok &= bool(fired)
    return ok


def table(rows):
    print("| alpha_c | w | N | M | max abs 2m/r | rs | R |\n|---|---|---|---|---|---|---|")
    for r in rows:
        if r.get("target"):
            print(f"| {r['alpha_c']:.4f} | {r['w']:.6f} | {r['N']:.4f} | {r['M']:+.4f} | {r['max2mr']:.3f} | "
                  f"{r['rs']:.1f} | {r['R']:.0f} |")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    mutate = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
    if "--selftest" in sys.argv:
        if not selftest():
            sys.exit(1)
        if not args:
            return
    rows, stop = load(args[0])
    modes = None
    if "--modes" in sys.argv:
        modes = [json.loads(ln) for ln in open(sys.argv[sys.argv.index("--modes") + 1])]
        if mutate == "sigma":
            modes[0]["roots"] = [-1e-3] + modes[0]["roots"]
    if mutate:
        rows = mutate_rows(rows, mutate)
    ok, _ = claims(rows, stop, modes)
    if "--table" in sys.argv:
        table(rows)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
