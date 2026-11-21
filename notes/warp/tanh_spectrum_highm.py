#!/usr/bin/env python3
"""Issue #125: the leading eigenvalue of m = 23/2, 29/2, 35/2 and 41/2 from #75's spectral code against #120's B-splines.

#120 computed m = 23/2 to 41/2 with tanh_spectrum_check.py's B-spline basis only, so its discretisations A and B share
the script, the formulation and the solver. This compares four of those sectors at sigma R = 2 and 4 with #75's code,
tanh_spectrum.py, whose basis (Jacobi polynomials in y = r/(r + Lmap) times spin-angular functions, dense eigensolver)
shares none of that code. What the two share: the operator S^-1 M, the tanh profile and v_c's definition.

    python notes/warp/tanh_spectrum.py --part=vc,planted,tanh --sigmas=2,4 --ms=10.5,11.5,14.5,17.5,20.5 --raw=DIR
    python notes/warp/tanh_spectrum_highm.py --raw=DIR --ref=DIR120 [--mutate=plant|conv]

DIR holds tanh_spectrum.py's files tanh_sigma{S}_m{M}.npz (vc, and muA and muB, every eigenvalue of S^-1 M in its two
discretisations A and B), DIR120 #120's files s{S}_{A|B}_m{2M}.npz (negative-mass-raw I120-win/eig). A leading value is
the largest |mu| v_c in a file.

Checks (exits 1 when any fails):
  vc      the v_c of every #75 file and of its two #120 namesakes agree to 1e-6 relative
  repro   m = 21/2: #75's A leading value, rounded to six decimals, is tanh_spectrum.md's printed one (0.083018, 0.094358)
  vs120   each of the eight (sigma R, m): #75's A leading value agrees with #120's A and with #120's B to 1e-4 relative
  conv    each of the ten sectors: #75's A and B leading values agree to 1e-6 relative, #75's own `converged` threshold.
          A and B are #75's two resolutions: dK 28 and 38, p 80 and 90 (sigma R = 2) or 100 and 110 (4), Lmap 2 and 3,
          npan 160 and 240
Mutations (each must exit 1 with its check failing):
  plant   #75's A and B values at sigma R = 4, m = 29/2 times 1.001: vs120 fails there and nowhere else
  conv    #75's B values at sigma R = 2, m = 41/2 times 1 + 1e-5: conv fails there and nowhere else
"""
import argparse
import os
import re
import sys

import numpy as np

SIGMAS = (2, 4)
MS = (11.5, 14.5, 17.5, 20.5)  # #120's sectors compared
FAILED = []


def check(name, ok, msg):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {msg}")
    if not ok:
        FAILED.append(name)


def lead75(raw, sig, m, mut):
    d = np.load(os.path.join(raw, f"tanh_sigma{sig:g}_m{m:g}.npz"))
    vc = float(d["vc"])
    a, b = np.abs(d["muA"]).max() * vc, np.abs(d["muB"]).max() * vc
    if mut == "plant" and (sig, m) == (4, 14.5):
        a, b = a * 1.001, b * 1.001
    if mut == "conv" and (sig, m) == (2, 20.5):
        b = b * (1 + 1e-5)
    return a, b, vc


def lead120(ref, sig, m, disc):
    d = np.load(os.path.join(ref, f"s{sig}_{disc}_m{int(round(2 * m))}.npz"))
    vc = float(d["vc"])
    return np.abs(d["mu"]).max() * vc, vc


def printed75(sig):
    """#75's printed leading |mu| v_c of m = 21/2, from tanh_spectrum.md section 4 (sigma R = 2 table first)."""
    md = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "tanh_spectrum.md"), encoding="utf-8").read()
    rows = re.findall(r"^\| 10\.5 \| [-\d.]+ \+ [\d.]+i \| ([\d.]+) \| [\d.]+ \| [\d.]+ \|$", md, re.M)
    assert len(rows) == 2, rows
    return rows[SIGMAS.index(sig)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--mutate", default="", choices=["", "plant", "conv"])
    args = ap.parse_args()
    for sig in SIGMAS:
        print(f"== sigma R = {sig}")
        a, b, vc = lead75(args.raw, sig, 10.5, args.mutate)
        want = printed75(sig)
        check(f"repro sR={sig}", f"{a:.6f}" == want, f"m = 10.5: #75's code {a:.8f} (A), {b:.8f} (B), printed {want}")
        check(f"conv sR={sig} m=10.5", abs(a - b) < 1e-6 * a, f"A and B differ by {abs(a - b) / a:.1e} relative")
        print("   m    | #75 A      | #75 B      | A vs B  | #120 A     | #120 B     | #75 A vs #120 A / B")
        for m in MS:
            a, b, vc = lead75(args.raw, sig, m, args.mutate)
            ra, vca = lead120(args.ref, sig, m, "A")
            rb, vcb = lead120(args.ref, sig, m, "B")
            da, db = abs(a - ra) / ra, abs(a - rb) / rb
            print(f"   {m:4.1f} | {a:.8f} | {b:.8f} | {abs(a - b) / a:.1e} | {ra:.8f} | {rb:.8f} | {da:.1e} / {db:.1e}")
            check(f"vc sR={sig} m={m}", max(abs(vc - vca), abs(vc - vcb)) < 1e-6 * vc,
                  f"#75 {vc:.6f}, #120 {vca:.6f} / {vcb:.6f}")
            check(f"vs120 sR={sig} m={m}", max(da, db) < 1e-4, f"worst {max(da, db):.1e} relative")
            check(f"conv sR={sig} m={m}", abs(a - b) < 1e-6 * a, f"A and B differ by {abs(a - b) / a:.1e} relative")
    print("FAILED: " + ", ".join(FAILED) if FAILED else "all checks passed")
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
