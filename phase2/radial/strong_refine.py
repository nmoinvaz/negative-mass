#!/usr/bin/env python3
"""Locate the canonical strong-field turning point in a strongfield.py JSON, by the sign change of the
fundamental sigma^2 and by the maxima of N and M along alpha_c.

    python strong_refine.py strongfield_eps+1.json [--selftest]

alpha_sigma interpolates sigma^2_0 = 0 linearly between the bracketing rows. alpha_N and alpha_M are the vertices of
parabolas through the largest row and its two neighbours. Exit 1 unless there is one sign change, both maxima are
interior, and |alpha_N - alpha_sigma| and |alpha_M - alpha_sigma| are at most the grid step.
--selftest runs planted rows with a known turning point, once as planted (must pass) and once with the sigma^2
crossing moved by 5 steps (must fail).
"""
import json
import sys


def vertex(xs, ys):
    i = max(range(len(ys)), key=ys.__getitem__)
    if i in (0, len(ys) - 1):
        return None
    (x0, x1, x2), (y0, y1, y2) = xs[i - 1:i + 2], ys[i - 1:i + 2]
    d = (x0 - x1) * (x0 - x2) * (x1 - x2)
    a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / d
    b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / d
    return -b / (2 * a)


def analyse(rows):
    rows = sorted(rows, key=lambda r: r["alpha_c"])
    al = [r["alpha_c"] for r in rows]
    s2 = [r["roots"][0] if r["roots"] else float("nan") for r in rows]
    step = min(b - a for a, b in zip(al, al[1:]))
    cross = [i for i in range(len(rows) - 1) if s2[i] * s2[i + 1] < 0]
    out = {"rows": len(rows), "step": step, "sign changes": len(cross)}
    if len(cross) != 1:
        return False, out
    i = cross[0]
    out["alpha_sigma"] = al[i] - s2[i] * (al[i + 1] - al[i]) / (s2[i + 1] - s2[i])
    out["alpha_N"] = vertex(al, [r["N"] for r in rows])
    out["alpha_M"] = vertex(al, [r["M"] for r in rows])
    if out["alpha_N"] is None or out["alpha_M"] is None:
        return False, out
    out["|alpha_N - alpha_sigma|"] = abs(out["alpha_N"] - out["alpha_sigma"])
    out["|alpha_M - alpha_sigma|"] = abs(out["alpha_M"] - out["alpha_sigma"])
    ok = max(out["|alpha_N - alpha_sigma|"], out["|alpha_M - alpha_sigma|"]) <= step
    return ok, out


def planted(shift):
    a0, h = 0.2734, 0.001
    rows = []
    for k in range(-9, 10):
        a = 0.271 + k * h
        rows.append({"alpha_c": a, "N": 42.7 - 50 * (a - a0) ** 2, "M": 11.9 - 14 * (a - a0) ** 2,
                     "roots": [-1e-3 * (a - a0 - shift * h)]})
    return rows, a0


def report(ok, out):
    print(("PASS" if ok else "FAIL") + ": " + ", ".join(f"{k} {v:.6g}" if isinstance(v, float) else f"{k} {v}"
                                                      for k, v in out.items()))


def main():
    if "--selftest" in sys.argv:
        rows, a0 = planted(0)
        ok, out = analyse(rows)
        exact = ok and abs(out["alpha_N"] - a0) < 1e-9 and abs(out["alpha_sigma"] - a0) < 1e-9
        print(f"selftest planted turning point {a0}: {'recovered' if exact else 'NOT recovered'}")
        report(ok, out)
        bent, _ = planted(5)
        ok2, out2 = analyse(bent)
        print(f"selftest crossing moved by 5 steps: {'FIRED' if not ok2 else 'DID NOT FIRE'}")
        report(ok2, out2)
        if not exact or ok2:
            sys.exit(1)
        args = [a for a in sys.argv[1:] if not a.startswith("--")]
        if not args:
            return
    ok, out = analyse(json.load(open([a for a in sys.argv[1:] if not a.startswith("--")][0])))
    report(ok, out)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
