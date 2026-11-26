#!/usr/bin/env python3
"""Tables and checks of step 3 (notes/wormhole_next/step3.md, issue #30), from the raw folder I30-mac:
counts/modecount_test_results.jsonl, ghost_mu3.jsonl and axis_s1.jsonl, and from this folder's ghost_mu.jsonl (the stored
real-axis zeros of section 1.3) and axis_s1.jsonl (section 6.2's scans and refinements).

    python3 step3_table.py [--raw=DIR] [--mutate=NAME] [--selftest]

DIR defaults to ../../../negative-mass-raw/I30-mac. Exit 1 if a check fails. Criteria as written in DIR/README.md:
  E1  each real-axis zero re-found in this environment equals the stored one to 1e-9
  C1  every plain box (B0, U0, full) at D = 10, 12, 20.04 is resolved and winds 0
  C2  the planted run: resolved, B0 = plain + 2, U0 and full as plain, the pair inside B0 (margin >= 0.1)
  C3  section 6.2's scan has exactly one phase jump, in (-1.05, -1.00), below the real-axis zero (D = 10, 12)
  C4  the rung-1 refinement (cell 0.005 inside (-1.05, -1.00)): (a) one jump, in one cell; (b) |Im/omega - 1| < 1e-5;
      (c) the zero's Re in the cell; (d) ln|g| at the zero >= 10 below the lower end sample of the cell (D = 10, 12)
  C5  the controls at D = 16.76 and 20.04 meet (a)-(d) in (-1.04, -1.00) and agree with section 6.2's zero to 1e-8
  X1  informative: rung 0 (section 6.2's start) at D = 10, 12 and whether C4's tests fire on it
--mutate=NAME breaks one input on purpose (expect exit 1 with only its target failing); --selftest runs every mutation
and exits 0 only if each fails exactly its targets.
"""
import copy
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_RAW = os.path.normpath(os.path.join(HERE, "..", "..", "..", "negative-mass-raw", "I30-mac"))
BOXES = ("B0", "U0", "full")
# D -> count tag, state file, tag of the stored real-axis zero in ghost_mu.jsonl, tag of E1's re-find in ghost_mu3.jsonl
STATES = {
    10: dict(tag="D10.0000_M48", state="up_D10.0000.npz", mu="D10_M48_m1", mu3="D10_M48_m1_i30"),
    12: dict(tag="D12.0000_M48", state="up_D12.0000.npz", mu="D12_M48_m1", mu3="D12_M48_m1_i30"),
    20.04: dict(tag="D20.0437_M32", state="up2_D20.0437.npz", mu="D20.04_M32_m1", mu3="D20.04_M32_m1_i30"),
}
CONTROLS = {16.76: "up2_D16.7600.npz", 20.04: "up2_D20.0437.npz"}
W_LO, W_HI = -1.05, -1.00          # section 6.2's phase-jump interval at D = 10, 12
C_LO, C_HI = -1.04, -1.00          # the controls' (section 6.2's interval at D >= 16.76)
H1, H0 = 0.005, 0.05               # rung 1's and rung 0's grid
ONLINE, DEPTH, REF_TOL, MU_TOL, MARGIN = 1e-5, 10.0, 1e-8, 1e-9, 0.1


def rows(path):
    return [json.loads(l) for l in open(path) if l.strip()] if os.path.exists(path) else []


def load(raw):
    return dict(counts=rows(os.path.join(raw, "counts", "modecount_test_results.jsonl")),
                mu3=rows(os.path.join(raw, "ghost_mu3.jsonl")), mu_stored=rows(os.path.join(HERE, "ghost_mu.jsonl")),
                axis=rows(os.path.join(raw, "axis_s1.jsonl")), old_axis=rows(os.path.join(HERE, "axis_s1.jsonl")),
                c3_zero={}, ref_shift=0.0, mu_shift=0.0, raw=raw)


def last(rs, pred):
    m = [r for r in rs if pred(r)]
    return m[-1] if m else None


def top_zero(rec):
    z = [x["mu"] for x in rec["zeros"] if x["ok"]]
    return max(z) if z else None


def jumps(rec):
    """phase jumps above pi/2 between neighbouring samples: (x_a, x_b, ln|g|_a, ln|g|_b)"""
    sm = rec["samples"]
    return [(a[0], b[0], a[1], b[1]) for a, b in zip(sm, sm[1:])
            if abs(math.remainder(b[2] - a[2], 2 * math.pi)) > math.pi / 2]


def axis_rec(rs, state, lo, hi, h):
    return last(rs, lambda r: r["state"] == state and abs(r["lo"] - lo) < 1e-9 and abs(r["hi"] - hi) < 1e-9
                and abs(r["h"] - h) < 1e-12)


def trend_depth(rec, z):
    """ln|g| of the line's samples interpolated at the zero's Re, minus ln|g| at the zero: tens for a zero, about 0 for
    an iterate that merely sits on the line's trend (None outside the samples)"""
    sm = rec["samples"]
    for a, b in zip(sm, sm[1:]):
        if a[0] <= z["re"] <= b[0]:
            return a[1] + (b[1] - a[1]) * (z["re"] - a[0]) / (b[0] - a[0]) - z["ln_g_at"]
    return None


def zero_tests(rec, lo, hi, h):
    """C4's (a)-(d) on one axis record: dict letter -> (ok, text)"""
    out = {}
    js = jumps(rec)
    cell = js[0] if len(js) == 1 else None
    out["a"] = (cell is not None and lo - 1e-9 <= cell[0] and cell[1] <= hi + 1e-9 and abs(cell[1] - cell[0] - h) < 1e-9,
                f"{len(js)} jump(s) {[(round(j[0], 6), round(j[1], 6)) for j in js]}")
    zs = [z for z in rec["zeros"] if "re" in z]
    z = zs[0] if len(zs) == 1 and len(rec["zeros"]) == 1 else None
    if z is None or cell is None:
        why = f"{len(rec['zeros'])} zero entries ({len(zs)} refined)" if z is None else "no single cell"
        for k in "bcd":
            out[k] = (False, why)
        return out
    off = z["im_over_omega"] - 1
    out["b"] = (abs(off) < ONLINE, f"Im/omega - 1 = {off:+.1e}")
    out["c"] = (cell[0] <= z["re"] <= cell[1], f"Re {z['re']:.10f} in cell ({cell[0]:.6f}, {cell[1]:.6f})")
    depth = min(cell[2], cell[3]) - z["ln_g_at"]
    out["d"] = (depth >= DEPTH, f"ln|g| at the zero {depth:+.1f} below the lower end sample")
    return out


def run_checks(d):
    """list of (name, ok, text)"""
    res = []
    add = lambda n, ok, t: res.append((n, bool(ok), t))   # noqa: E731
    # ---- E1: the real-axis zeros in this environment
    for D, s in STATES.items():
        st = last(d["mu_stored"], lambda r: r["tag"] == s["mu"] and r["mode"] == 1 and r["plant"] is None)
        new = last(d["mu3"], lambda r: r["tag"] == s["mu3"] and r["mode"] == 1 and r["plant"] is None)
        if st is None or new is None:
            add(f"E1 D{D}", False, "no record" + ("" if st else " (stored)") + ("" if new else " (re-found)"))
            continue
        a, b = top_zero(st) + d["mu_shift"], top_zero(new)
        okz = [z for z in new["zeros"] if z["ok"]]
        add(f"E1 D{D}", b is not None and len(okz) == 1 and abs(a - b) < MU_TOL,
            f"stored {a:.12f}, re-found {b if b is None else format(b, '.12f')}, difference "
            f"{'-' if b is None else format(abs(a - b), '.1e')} (limit {MU_TOL:g}), {len(okz)} zero(s) in the window")
    # ---- C1, C2: the counts
    for D, s in STATES.items():
        tag = "test_" + s["tag"]
        ok_box = lambda r: set(r["boxes"]) == set(BOXES)   # noqa: E731
        plain = last(d["counts"], lambda r: r["tag"] == tag and r.get("plant") is None and ok_box(r))
        pl = last(d["counts"], lambda r: r["tag"].startswith(tag + "_plant") and ok_box(r))
        if plain is None:
            add(f"C1 D{D}", False, "no plain record")
        else:
            w = {b: plain["boxes"][b]["winding_g_int"] for b in BOXES}
            res_ = all(plain["boxes"][b]["resolved"] for b in BOXES)
            add(f"C1 D{D}", res_ and all(v == 0 for v in w.values()), f"resolved {res_}, windings {w}")
        if plain is None or pl is None:
            add(f"C2 D{D}", False, "no planted run" if pl is None else "no plain record")
            continue
        z = complex(*pl["plant"])
        members = [z, z.conjugate()]
        inside = lambda m, r: r[0] < m.real < r[1] and r[2] < m.imag < r[3]   # noqa: E731
        want = {b: plain["boxes"][b]["winding_g_int"] + sum(inside(m, pl["boxes"][b]["rect"]) for m in members)
                for b in BOXES}
        got = {b: pl["boxes"][b]["winding_g_int"] for b in BOXES}
        r0 = pl["boxes"]["B0"]["rect"]
        margin = min(z.real - r0[0], r0[1] - z.real, r0[3] - abs(z.imag))
        res_ = all(pl["boxes"][b]["resolved"] for b in BOXES)
        add(f"C2 D{D}", got == want and margin >= MARGIN and res_ and got["B0"] == plain["boxes"]["B0"]["winding_g_int"] + 2,
            f"plant {z.real:+g} +- {abs(z.imag):g}i (margin {margin:.3f} inside B0): windings {got}, expected {want}, "
            f"resolved {res_}")
    # ---- C3, C4: the other sector at D = 10, 12
    for D in (10, 12):
        s = STATES[D]
        st = last(d["mu_stored"], lambda r: r["tag"] == s["mu"] and r["mode"] == 1 and r["plant"] is None)
        real = d["c3_zero"].get(D, top_zero(st) if st else None)
        old = last(d["old_axis"], lambda r: r["state"] == s["state"] and abs(r["lo"] + 1.6) < 1e-9)
        if old is None or real is None:
            add(f"C3 D{D}", False, "no section 6.2 scan or no real-axis zero")
        else:
            js = jumps(old)
            add(f"C3 D{D}", len(js) == 1 and abs(js[0][0] - W_LO) < 1e-9 and abs(js[0][1] - W_HI) < 1e-9
                and js[0][1] < real, f"section 6.2's scan: jumps {[(round(j[0], 4), round(j[1], 4)) for j in js]}, "
                f"real-axis zero {real:.6f}")
        r1 = axis_rec(d["axis"], s["state"], W_LO, W_HI, H1)
        if r1 is None:
            add(f"C4 D{D}", False, "no rung-1 record")
            continue
        t = zero_tests(r1, W_LO, W_HI, H1)
        add(f"C4 D{D}", all(v[0] for v in t.values()), "; ".join(f"({k}) {'ok' if v[0] else 'FAIL'} {v[1]}" for k, v in t.items()))
    # ---- C5: the controls
    for D, state in CONTROLS.items():
        r = axis_rec(d["axis"], state, C_LO, C_HI, H1)
        old = last(d["old_axis"], lambda r_: r_["state"] == state and abs(r_["lo"] + 1.6) < 1e-9)
        ref = [z["re"] for z in old["zeros"] if "re" in z] if old else []
        if r is None or len(ref) != 1:
            add(f"C5 D{D}", False, "no control record or no section 6.2 reference")
            continue
        t = zero_tests(r, C_LO, C_HI, H1)
        z = [z for z in r["zeros"] if "re" in z]
        diff = abs(z[0]["re"] - (ref[0] + d["ref_shift"])) if len(z) == 1 else math.inf
        add(f"C5 D{D}", all(v[0] for v in t.values()) and diff < REF_TOL,
            "; ".join(f"({k}) {'ok' if v[0] else 'FAIL'} {v[1]}" for k, v in t.items())
            + f"; against section 6.2's {ref[0]:.10f}: {diff:.1e} (limit {REF_TOL:g})")
    return res


MUTATIONS = {   # name -> (function on the data, set of check prefixes that must fail)
    "mu": (lambda d: d.__setitem__("mu_shift", 1e-6), {"E1"}),
    "winding": (lambda d: [r["boxes"]["U0"].__setitem__("winding_g_int", r["boxes"]["U0"]["winding_g_int"] + 1)
                           for r in d["counts"] if r["tag"].startswith("test_D10.0000_M48")], {"C1"}),
    "noplant": (lambda d: [r["boxes"].update(copy.deepcopy(next(p for p in d["counts"] if p["tag"] == "test_D12.0000_M48")
                                                              ["boxes"])) for r in d["counts"]
                           if r["tag"].startswith("test_D12.0000_M48_plant")], {"C2"}),
    "realzero": (lambda d: d["c3_zero"].update({10: -1.2, 12: -1.2}), {"C3"}),
    "offline": (lambda d: [z.__setitem__("im_over_omega", z["im_over_omega"] + 1e-3) for r in d["axis"]
                           if r["state"].startswith("up_D") and abs(r["h"] - H1) < 1e-12 for z in r["zeros"] if "re" in z],
                {"C4"}),
    "trend": (lambda d: [z.__setitem__("ln_g_at", jumps(r)[0][2]) for r in d["axis"]
                         if r["state"].startswith("up_D") and abs(r["h"] - H1) < 1e-12 and jumps(r)
                         for z in r["zeros"] if "re" in z], {"C4"}),
    "old": (lambda d: [r.__setitem__("zeros", copy.deepcopy(o["zeros"])) for r in d["axis"]
                       if r["state"].startswith("up_D") and abs(r["h"] - H1) < 1e-12
                       for o in d["old_axis"] if o["state"] == r["state"] and abs(o["lo"] + 1.6) < 1e-9], {"C4"}),
    "ctrl": (lambda d: d.__setitem__("ref_shift", 1e-6), {"C5"}),
}


def show(d, res):
    print(f"raw folder {d['raw'].replace(os.path.expanduser('~'), '~')}")
    ex = {}
    for f in sorted(glob.glob(os.path.join(d["raw"], "logs", "*.log"))):
        e = [l.split()[3] for l in open(f) if l.startswith("== ") and " exit " in l]
        ex[os.path.basename(f)[:-4]] = e[-1] if e else "running"
    print(f"jobs (informative): {len(ex)} logs, exit codes {sorted(set(ex.values()))}; not exit 0: "
          f"{[k for k, v in ex.items() if v != '0'] or 'none'}")
    envs = sorted({l.strip() for f in glob.glob(os.path.join(d["raw"], "logs", "*.log")) for l in open(f)
                   if l.startswith("-- env:")})
    print("environment lines in the job logs (informative):")
    for e in envs:
        print("  " + e)
    if len(envs) > 1:
        print("  WARNING: the environment changed between jobs")
    print("\nreal-axis zeros (E1), stored (ghost_mu.jsonl, section 1.3) against re-found here (ghost_mu3.jsonl)")
    for D, s in STATES.items():
        st = last(d["mu_stored"], lambda r: r["tag"] == s["mu"] and r["mode"] == 1 and r["plant"] is None)
        new = last(d["mu3"], lambda r: r["tag"] == s["mu3"] and r["mode"] == 1 and r["plant"] is None)
        if st and new:
            ws = [z["width"] for z in new["zeros"] if z["ok"]]
            print(f"  D {D}: stored {top_zero(st):.12f} (M {st['M']}), re-found {top_zero(new):.12f} (M {new['M']}, "
                  f"width {max(ws):.1e}), {new['evaluations']} evaluations")
    print("\ncounts (aLL < Re mu < aR, |Im mu| < 1.5 omega; B0, U0 left of aL, full right of it)")
    for D, s in STATES.items():
        tag = "test_" + s["tag"]
        plain = last(d["counts"], lambda r: r["tag"] == tag and r.get("plant") is None)
        pl = last(d["counts"], lambda r: r["tag"].startswith(tag + "_plant"))
        if plain:
            w = {b: v["winding_g_int"] for b, v in plain["boxes"].items()}
            print(f"  D {plain['D']:.4f} M {plain['M']} ({plain['state']}): aLL {plain['aLL']:.4f}, aL {plain['aL']:.4f}, "
                  f"aR {plain['aR']:.4f}, y2 {plain['y2']:.3f}; windings {w}; max step "
                  f"{max(v['max_step'] for v in plain['boxes'].values()):.3f}; {plain['evaluations']} evaluations "
                  f"({plain['new_evaluations']} new, {plain['sec_per_eval'] or 0:.1f} s each)")
        rawf = os.path.join(d["raw"], "counts", f"modecount_{tag}_raw.jsonl")
        if plain and os.path.exists(rawf):
            print(f"      determinants evaluated for the plain and planted runs together: {len(rows(rawf))}")
        if pl:
            w = {b: v["winding_g_int"] for b, v in pl["boxes"].items()}
            print(f"      planted {pl['plant'][0]:+g} +- {abs(pl['plant'][1]):g}i: windings {w}; {pl['evaluations']} "
                  f"evaluations ({pl['new_evaluations']} new)")
    print("\nother parity sector, Im mu = omega")
    for D, s in STATES.items():
        if D == 20.04:
            continue
        old = last(d["old_axis"], lambda r: r["state"] == s["state"] and abs(r["lo"] + 1.6) < 1e-9)
        oz = [z for z in old["zeros"] if "re" in z] if old else []
        print(f"  D {D} ({s['state']}): section 6.2's scan (h 0.05): jumps "
              f"{[(round(j[0], 4), round(j[1], 4)) for j in jumps(old)] if old else None}; its refinement "
              f"{[(round(z['re'], 10), round(z['im_over_omega'], 6)) for z in oz]}")
        for z in oz:
            td = trend_depth(old, z)
            print(f"      section 6.2's point: ln|g| there is {td if td is None else format(td, '+.3f')} below the line's "
                  f"trend interpolated from its samples (a zero: tens; on the trend: about 0)")
        for name, h in (("rung 0, section 6.2's start (X1)", H0), ("rung 1, the record", H1), ("rung 2", 0.002)):
            r = axis_rec(d["axis"], s["state"], W_LO, W_HI, h)
            if r is None:
                print(f"      {name}: not run")
                continue
            t = zero_tests(r, W_LO, W_HI, h)
            zz = [z for z in r["zeros"] if "re" in z]
            pos = (f"zero Re {zz[0]['re']:.10f}, Im/omega - 1 {zz[0]['im_over_omega'] - 1:+.2e}, ln|g| {zz[0]['ln_g_at']:.2f}"
                   if zz else f"zeros {r['zeros']}")
            print(f"      {name}: {pos}; tests " + " ".join(f"({k}){'ok' if v[0] else 'FAIL'}" for k, v in t.items())
                  + f"; {r['sec']:.0f} s")
            if h == H0:
                fired = [k for k, v in t.items() if not v[0]]
                same = (oz and zz and abs(oz[0]['re'] - zz[0]['re']) < 1e-6 and abs(oz[0]['im_over_omega'] -
                                                                                     zz[0]['im_over_omega']) < 1e-6)
                dre = abs(oz[0]['re'] - zz[0]['re']) if oz and zz else None
                dim = abs(oz[0]['im_over_omega'] - zz[0]['im_over_omega']) if oz and zz else None
                print(f"      X1: C4's tests {'fire (' + ','.join(fired) + ')' if fired else 'are SILENT'} on the old start; "
                      f"it {'reproduces' if same else 'does not reproduce'} section 6.2's point to 1e-6 "
                      f"(differences in Re {dre if dre is None else format(dre, '.1e')}, in Im/omega "
                      f"{dim if dim is None else format(dim, '.1e')}; the window's pole set, hence g, differs from "
                      f"the full scan's)")
    print("  controls (C5):")
    for D, state in CONTROLS.items():
        r = axis_rec(d["axis"], state, C_LO, C_HI, H1)
        zz = [z for z in r["zeros"] if "re" in z] if r else []
        print(f"    D {D}: " + (f"zero Re {zz[0]['re']:.10f}, Im/omega - 1 {zz[0]['im_over_omega'] - 1:+.2e}" if zz else "not run"))
    print("\nchecks")
    for n, ok, t in res:
        print(f"{'PASS' if ok else 'FAIL'}  {n}: {t}")


def main():
    args = dict(a[2:].split("=", 1) if "=" in a else (a[2:], "1") for a in sys.argv[1:] if a.startswith("--"))
    raw = args.get("raw", DEFAULT_RAW)
    d = load(raw)
    if "selftest" in args:
        bad = 0
        for name, (fn, targets) in MUTATIONS.items():
            m = copy.deepcopy(d)
            fn(m)
            failed = {n.split()[0] for n, ok, _ in run_checks(m) if not ok}
            base = {n.split()[0] for n, ok, _ in run_checks(d) if not ok}
            fired = (failed - base) == targets and targets <= failed
            print(f"{'PASS' if fired else 'FAIL'}  --mutate={name}: fails {sorted(failed)} (unmutated fails {sorted(base)}), "
                  f"target {sorted(targets)}")
            bad += not fired
        sys.exit(1 if bad else 0)
    if "mutate" in args:
        fn, targets = MUTATIONS[args["mutate"]]
        fn(d)
        print(f"MUTATED: --mutate={args['mutate']} (target {sorted(targets)}); the exit must be 1")
    res = run_checks(d)
    show(d, res)
    sys.exit(1 if any(not ok for _, ok, _ in res) else 0)


if __name__ == "__main__":
    main()
