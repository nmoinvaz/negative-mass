#!/usr/bin/env python3
"""Radial spectra (sigma^2 below the continuum threshold) of the stars saved by phantom_follow.py.

    python phantom_modes.py SNAP_DIR OUT.jsonl [--alphas 1.1,2.0,...] [--nneg 24] [--npos 32] [--s2min -0.2]
    python phantom_modes.py --compare OUT.jsonl strong_m.txt [--mutate=shift]
    python phantom_modes.py --grids A.jsonl B.jsonl [--mutate=drop|shift]
    python phantom_modes.py --selftest

Each snapshot star_a<alpha_c>.npz is re-solved at its own alpha_c (the BVP restarts on its saved mesh) and handed to
radial_modes.RadialProblem with the same lowest_modes as strongfield.py. One JSON line per star with the roots and the
grid, so a root pair lost between grid points can be looked for by rerunning with a denser grid.

--grids compares two runs of the same stars on different root grids: every root of the coarser grid must be found on the
finer one to GRID_TOL relative, and the finer grid may add no root below 0 and no root below the coarser one's highest
(a lost pair would show as an extra root). --mutate=drop deletes the finer grid's lowest root, --mutate=shift scales it
by 1 + 1e-4, and each must make it fail.
--compare checks the lowest root against the stored strong_m.txt rows (strongfield.py, alpha_c <= 1.5, the same stars on
R = 60) to REL_TOL relative. --mutate=shift scales the new roots by 1 + 1e-3 and must make it fail. --selftest runs the
compare on planted rows, as planted (must pass) and with a 1e-3 shift (must fire).
"""
import glob
import json
import os
import re
import sys
import time
import warnings

warnings.filterwarnings("ignore")
import numpy as np  # noqa: E402

REL_TOL = 6e-5            # the stored roots are printed with 5 digits, a relative 5e-5 at worst


def run(snap_dir, out_path, alphas, nneg, npos, s2min):
    import phantom_follow as F
    from radial_modes import RadialProblem, lowest_modes
    files = sorted(glob.glob(os.path.join(snap_dir, "star_a*.npz")))
    out = open(out_path, "w")
    for f in files:
        a = float(re.search(r"star_a([\d.]+)\.npz", f).group(1))
        if alphas and not any(abs(a - x) < 1e-9 for x in alphas):
            continue
        g, ac = F.load_checkpoint(f)
        st = F.solve_lapse(ac, F.EPS, g, g.R)
        t0 = time.time()
        P = RadialProblem(st)
        roots, grid, vals = lowest_modes(P, s2_min=s2min, nneg=nneg, npos=npos, nmax=3)
        thr = (1.0 - st.w) ** 2
        row = dict(alpha_c=ac, w=st.w, N=st.N, M=st.M, R=st.R, rm=P.rm, thr=thr, roots=[float(x) for x, o in roots],
                   ok=[bool(o) for x, o in roots], s2_min=s2min, nneg=nneg, npos=npos, ngrid=int(len(grid)),
                   sign_changes=int(np.sum(np.sign(vals[:-1]) != np.sign(vals[1:]))), t=time.time() - t0)
        out.write(json.dumps(row) + "\n"); out.flush()
        print(f"alpha_c={ac:.4f} w={st.w:.6f} N={st.N:.3f} M={st.M:+.3f} rm={P.rm:.1f} R={st.R:.0f} thr={thr:.3e} "
              f"roots={['%.6e%s' % (x, '' if o else '?') for x, o in roots]} t={row['t']:.0f}s", flush=True)
    out.close()


STRONG = re.compile(r"eps=-1 alpha_c=([\d.]+) w=([\d.]+) N=([\d.]+) M=([+-][\d.]+) sigma2=\[([^\]]*)\]")


def stored_strong(text):
    out = {}
    for ln in text.splitlines():
        m = STRONG.match(ln)
        if m and m[5].strip():
            out[round(float(m[1]), 4)] = float(m[5].split(",")[0].strip().strip("'").rstrip("?"))
    return out


def compare(rows, stored, mutate=None):
    fails, worst, n = [], 0.0, 0
    for r in rows:
        key = round(r["alpha_c"], 4)
        if key not in stored or not r["roots"]:
            continue
        v = r["roots"][0] * (1 + 1e-3 if mutate == "shift" else 1)
        d = abs(v - stored[key]) / abs(stored[key])
        n += 1
        worst = max(worst, d)
        if d > REL_TOL:
            fails.append(f"alpha_c={key}: new {v:.6e} stored {stored[key]:.6e} rel {d:.2e}")
    if n == 0:
        fails.append("no star in common with the stored rows")
    print(("PASS" if not fails else "FAIL") + f" lowest sigma^2 against stored strong_m.txt: {n} stars, worst rel {worst:.2e}, "
          f"limit {REL_TOL:g}")
    for f in fails[:8]:
        print("   ", f)
    return not fails


GRID_TOL = 1e-8


def compare_grids(a_rows, b_rows, mutate=None):
    """a: the coarser grid, b: the finer. Roots are matched by value, a root of a must be in b."""
    fails, n = [], 0
    b_by = {round(r["alpha_c"], 4): [float(x) for x in r["roots"]] for r in b_rows}
    for r in a_rows:
        key = round(r["alpha_c"], 4)
        if key not in b_by:
            continue
        fine = list(b_by[key])
        if mutate == "drop" and fine:
            fine = fine[1:]
        if mutate == "shift" and fine:
            fine[0] *= 1 + 1e-4
        n += 1
        coarse = [float(x) for x in r["roots"]]
        for x in coarse:
            if not any(abs(x - y) <= GRID_TOL * max(abs(x), abs(y), 1e-300) + 1e-12 for y in fine):
                fails.append(f"alpha_c={key}: root {x:.9e} of the coarse grid not on the fine grid {fine}")
        extra = [y for y in fine if not any(abs(x - y) <= GRID_TOL * max(abs(x), abs(y)) + 1e-12 for x in coarse)]
        extra = [y for y in extra if y < 0 or y < (max(coarse) if coarse else float("inf"))]
        if extra:
            fails.append(f"alpha_c={key}: roots {extra} only on the fine grid")
    if n == 0:
        fails.append("no star in common")
    print(("PASS" if not fails else "FAIL") + f" root grids: {n} stars compared, tolerance {GRID_TOL:g}")
    for f in fails[:8]:
        print("   ", f)
    return not fails


def selftest():
    text = "\n".join(f"eps=-1 alpha_c={a:.4f} w=0.43 N=2.0 M=-0.9 sigma2=['{s:.4e}', '9e-2']" for a, s in
                     ((1.1, 3.594e-2), (1.15, 2.4925e-2), (1.5, 1.2704e-2)))
    stored = stored_strong(text)
    rows = [dict(alpha_c=a, roots=[s * (1 + 3e-5)]) for a, s in ((1.1, 3.594e-2), (1.15, 2.4925e-2), (1.5, 1.2704e-2))]
    a = compare(rows, stored)
    b = compare(rows, stored, "shift")
    print(f"selftest: as planted {'passes' if a else 'FAILS'}; --mutate=shift {'FIRED' if not b else 'DID NOT FIRE'}")
    coarse = [dict(alpha_c=3.0, roots=[1.8685e-2, 4.5649e-2]), dict(alpha_c=4.0, roots=[])]
    fine = [dict(alpha_c=3.0, roots=[1.8685e-2, 4.5649e-2]), dict(alpha_c=4.0, roots=[])]
    c = compare_grids(coarse, fine)
    d = [not compare_grids(coarse, fine, m) for m in ("drop", "shift")]
    extra = compare_grids(coarse, [dict(alpha_c=3.0, roots=[-1e-3, 1.8685e-2, 4.5649e-2]), fine[1]])
    print(f"selftest grids: as planted {'passes' if c else 'FAILS'}; drop {'FIRED' if d[0] else 'DID NOT FIRE'}; "
          f"shift {'FIRED' if d[1] else 'DID NOT FIRE'}; an added negative root {'FIRED' if not extra else 'DID NOT FIRE'}")
    return a and not b and c and all(d) and not extra


def main():
    argv = sys.argv[1:]
    mutate = next((x.split("=", 1)[1] for x in argv if x.startswith("--mutate=")), None)
    opt = {}
    pos = []
    i = 0
    while i < len(argv):
        if argv[i] in ("--alphas", "--nneg", "--npos", "--s2min"):
            opt[argv[i]] = argv[i + 1]; i += 2
        elif argv[i].startswith("--"):
            i += 1
        else:
            pos.append(argv[i]); i += 1
    if "--selftest" in argv:
        sys.exit(0 if selftest() else 1)
    if "--grids" in argv:
        sys.exit(0 if compare_grids([json.loads(ln) for ln in open(pos[0])], [json.loads(ln) for ln in open(pos[1])],
                                    mutate) else 1)
    if "--compare" in argv:
        rows = [json.loads(ln) for ln in open(pos[0])]
        sys.exit(0 if compare(rows, stored_strong(open(pos[1]).read()), mutate) else 1)
    alphas = [float(x) for x in opt["--alphas"].split(",")] if "--alphas" in opt else None
    run(pos[0], pos[1], alphas, int(opt.get("--nneg", 24)), int(opt.get("--npos", 32)), float(opt.get("--s2min", -0.2)))


if __name__ == "__main__":
    main()
