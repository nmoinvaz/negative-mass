#!/usr/bin/env python3
"""The pilot's table and checks (notes/wormhole_next.md section 1.3), from ghost_mu.jsonl.

    python3 ghost_mu_table.py [--cld=/path/to/choptuik-large-d] [--raw=/path/to/choptuik-large-d-raw]

gamma(D) is read from choptuik-large-d's double-double records (hp/data/local/gamma_dd.jsonl, D >= 8.25) and the raw
repository's T79-win/gamma.jsonl (D < 8.25). Prints mu per state, the fade-out exponent 2|mu|gamma of
notes/bh_to_wormhole.md section 3 and the charge exponent delta_Q = (D - 2 - mu) gamma of [MGM26] Eq. (66).

Checks (exit 1 if one fails):
  1. D 4 reproduces choptuik-large-d's real zero -0.3651202372 (notes/GHOST_ADMIXTURE.md section 4) to 1e-8.
  2. D 4 lies within [GMG96]'s -0.362 +- 0.012.
  3. D 5 lies within [MGM26]'s window -0.450 .. -0.446 (no error bar is given there; the window is theirs).
  4. Resolution: the M 48 and M 64 zeros at D 5 and D 6 agree to 1e-6; the second node layout at D 5 to 1e-6.
  5. Planted zero: the run with --plant=-0.17 returns -0.17 to 1e-9.
  6. Control (mode 2, a's response kept): its smallest pivot at mu = 0 exceeds 1e3 times mode 1's, and none of its
     real zeros lies within 0.012 of mode 1's (so check 2 or 3 would fail with that operator).
"""
import json
import os
import sys


def parse(argv):
    return dict(a[2:].split("=", 1) for a in argv if a.startswith("--") and "=" in a)


def main():
    opts = parse(sys.argv[1:])
    here = os.path.dirname(os.path.abspath(__file__))
    cld = opts.get("cld", "/Users/nathan/Source/choptuik-large-d")
    raw = opts.get("raw", "/Users/nathan/Source/choptuik-large-d-raw")
    gam = {}
    for f in (os.path.join(raw, "T79-win", "gamma.jsonl"), os.path.join(cld, "hp", "data", "local", "gamma_dd.jsonl")):
        for line in open(f):
            r = json.loads(line)
            gam[round(r["D"], 4)] = r["gamma"]
    recs = [json.loads(l) for l in open(os.path.join(here, "ghost_mu.jsonl"))]
    by = {r["tag"]: r for r in recs}

    def mu(tag, near=None):
        zs = [z["mu"] for z in by[tag]["zeros"] if z["ok"]]
        if near is not None:
            zs = sorted(zs, key=lambda z: abs(z - near))
        return zs[0] if zs else None

    print("tag                 D        M   K   |F|      mu (real axis)     2|mu|gamma  delta_Q   zeros in window")
    rows = []
    for r in recs:
        if r["mode"] != 1 or r["plant"] is not None:
            continue
        zs = [z["mu"] for z in r["zeros"] if z["ok"]]
        lead = max(zs) if zs else None
        g = gam.get(round(r["D"], 4))
        fe = 2 * abs(lead) * g if (lead is not None and g) else None
        dq = (r["D"] - 2 - lead) * g if (lead is not None and g) else None
        rows.append((r["D"], r["M"], lead))
        print(f"{r['tag']:18s}  {r['D']:7.4f}  {r['M']:3d}  {r['K']:2d}  {r['F']:.1e}  "
              f"{lead if lead is None else f'{lead:+.10f}':17s}  "
              f"{'' if fe is None else f'{fe:.4f}':10s}  {'' if dq is None else f'{dq:.4f}':8s}  "
              f"{[round(z, 6) for z in zs]} in [{r['lo']}, {r['hi']}]")
    bad = []

    def check(name, ok, msg):
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {msg}")
        if not ok:
            bad.append(name)

    m4 = mu("D4_M64_m1", -0.365)
    check("1 D4 reproduces", abs(m4 + 0.3651202372) < 1e-8, f"{m4:+.10f} against -0.3651202372")
    check("2 D4 GMG96", abs(m4 + 0.362) <= 0.012, f"{m4:+.4f} against -0.362 +- 0.012")
    m5 = mu("D5_M48_m1", -0.448)
    check("3 D5 MGM26", -0.450 <= m5 <= -0.446, f"{m5:+.4f} against the window -0.450 .. -0.446 (off by "
          f"{m5 + 0.446:+.4f} from its edge)")
    for a, b in (("D5_M48_m1", "D5_M64_m1"), ("D6_M48_m1", "D6_M64_m1"), ("D8.25_M32_m1", "D8.25_M48_m1"),
                 ("D5_M48_m1", "D5_M48_n11_m1")):
        if a in by and b in by:
            x, y = mu(a, -0.5), mu(b, mu(a, -0.5))
            check(f"4 {a} vs {b}", y is not None and abs(x - y) < 1e-6, f"{x:+.10f} {y if y is None else f'{y:+.10f}'}")
        else:
            check(f"4 {a} vs {b}", False, "missing run")
    if "D5_M48_plant" in by:
        p = mu("D5_M48_plant", -0.17)
        check("5 planted", p is not None and abs(p + 0.17) < 1e-9, f"{p} against -0.17")
    else:
        check("5 planted", False, "missing run")
    if "D5_M48_m2" in by:
        c = by["D5_M48_m2"]
        z2 = [z["mu"] for z in c["zeros"] if z["ok"]]
        near = min((abs(z - m5) for z in z2), default=float("inf"))
        piv = c["pivot0"] / by["D5_M48_m1"]["pivot0"]
        check("6 control", piv > 1e3 and near > 0.012,
              f"pivot at 0 {c['pivot0']:.1e} ({piv:.1e} x mode 1's), zeros {[round(z, 4) for z in z2]}, nearest to "
              f"mode 1's {near:.4f} away")
    else:
        check("6 control", False, "missing run")
    # mu_inf: least-squares fits in 1/n on D >= 8.25, one value per D (the highest M). A spread over fit choices, not
    # validated by planted truth (notes/wormhole_next.md section 1.3); printed, not checked.
    import numpy as np
    pts = {}
    for D, M, lead in rows:
        if D >= 8.2 and lead is not None and (D not in pts or M > pts[D][0]):
            pts[D] = (M, lead)
    Ds = sorted(pts)
    if len(Ds) >= 4:
        n = np.array([D - 2 for D in Ds])
        y = np.array([pts[D][1] for D in Ds])
        print(f"mu_inf fits on D = {[round(D, 2) for D in Ds]}:")
        for k in (1, 2, 3):
            for drop in (0, 3):
                if len(Ds) - drop <= k + 1:
                    continue
                A = np.vander(1 / n[drop:], k + 1, increasing=True)
                c, *_ = np.linalg.lstsq(A, y[drop:], rcond=None)
                r = float(np.abs(A @ c - y[drop:]).max())
                print(f"  terms 1/n^0..{k}, from D = {Ds[drop]:.2f}: mu_inf = {c[0]:+.4f}, max residual {r:.1e}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
