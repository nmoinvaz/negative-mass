#!/usr/bin/env python3
"""Table and figure of the continued phantom family, from the rows phantom_follow.py wrote.

    python phantom_graph.py FOLLOW.jsonl MODES.jsonl STOP.json phantom_family.txt phantom_family.svg

phantom_family.txt: every STRIDE-th row, the target rows and the last row (alpha_c, w, k, N, M, max|2m/r|, rs, R).
phantom_family.svg: N and |M| against w along the whole branch (a), the fundamental sigma^2 and the continuum threshold
(1 - w)^2 against alpha_c (b), and the star radius rs and 2m/r against alpha_c (c). It reads only the files it is given.
"""
import json
import sys

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

STRIDE = 10


def main():
    follow, modes, stop, txt, svg = sys.argv[1:6]
    rows = [r for r in map(json.loads, open(follow)) if not r.get("stop")]
    mrows = [json.loads(ln) for ln in open(modes)]
    s = json.load(open(stop))
    a_star = s["tests"]["S1"]["alpha_star"]
    keep = [r for i, r in enumerate(rows) if i % STRIDE == 0 or r.get("target") or i == len(rows) - 1]
    with open(txt, "w") as f:
        f.write("# phantom family (eps = -1) in the central lapse alpha_c, from phantom_follow.py (rows) and\n")
        f.write(f"# phantom_stop.py (w = 1 at alpha_c* = {a_star:.7f}); every {STRIDE}th row, the target rows and the last row.\n")
        f.write("# alpha_c  w  k=sqrt(1-w^2)  N  M  max|2m/r|  rs  R\n")
        for r in keep:
            f.write(f"{r['alpha_c']:.6f} {r['w']:.8f} {r['k']:.6f} {r['N']:.5f} {r['M']:+.5f} {r['max2mr']:.4f} "
                    f"{r['rs']:.2f} {r['R']:.0f}\n")
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.4))
    w = [r["w"] for r in rows]
    ax[0].plot(w, [r["N"] for r in rows], color="#1f77b4", label="N")
    ax[0].plot(w, [-r["M"] for r in rows], color="#d62728", label="|M|")
    ax[0].axvline(1.0, color="0.6", lw=0.8, ls=":")
    ax[0].set_xlabel("ω"); ax[0].set_ylabel("N, |M|"); ax[0].legend(frameon=False)
    ax[0].set_title("(a) N and |M| along the branch")
    ax[1].semilogy([m["alpha_c"] for m in mrows if m["roots"]], [m["roots"][0] for m in mrows if m["roots"]], "o-",
                   color="#2ca02c", label="σ² of the fundamental mode")
    ax[1].semilogy([m["alpha_c"] for m in mrows], [m["thr"] for m in mrows], "s--", color="0.4", label="(1 − ω)²")
    ax[1].set_xlabel("α₀"); ax[1].legend(frameon=False, fontsize=8)
    ax[1].set_title("(b) σ² and the continuum threshold")
    ax[2].plot([r["alpha_c"] for r in rows], [r["rs"] for r in rows], color="#9467bd", label="r(φ = 10⁻³φ_c)")
    ax[2].plot([r["alpha_c"] for r in rows], [r["max2mr"] for r in rows], color="#ff7f0e", label="max |2m/r|")
    ax[2].axvline(a_star, color="0.6", lw=0.8, ls=":")
    ax[2].set_xlabel("α₀"); ax[2].legend(frameon=False, fontsize=8)
    ax[2].set_title("(c) size and field strength")
    fig.tight_layout()
    fig.savefig(svg)


if __name__ == "__main__":
    main()
