#!/usr/bin/env python3
"""#191: the singularity locator of multisolver (branch feat/locate-singularity) on the series #183's walk read at its last step.

#183's raw folder keeps no Taylor coefficients, so each stored point (x_last and the state of its last step in
double-double, `s_last`) is expanded again by #183's own `hp::taylor_p` in double-double to order 30, as its double-double
walks did, in the step variable s = (x - x_last)/l with l = 2^round(log2 est.r); the harness (negative-mass-raw,
I191-rpi5/sys/compare183) runs the locator on the three fields F, F', F''. Field f of F ~ (x - x*)^alpha has the exponent
alpha - f. This script selects the points, runs the harness and prints, for each point, eta* and the exponent beside #183's
`xstar`, `bar`, fitted exponent `est.alpha` and `alpha_pred` = 2/(eta*^2 + 4).

    python3 singularity_locator.py RAW_I183_DUMP HARNESS_BIN [PLANTED_LOG] > singularity_locator_output.txt

PLANTED_LOG: the output of `cargo test --release --test locate -- --nocapture` of multisolver (branch feat/locate-singularity);
its table of the planted suite is summarized first (cells that meet the criteria, the bars against the errors).

RAW_I183_DUMP: negative-mass-raw/I183-rpi5/dump_run_03a; HARNESS_BIN: the built `compare183` binary. Read only: nothing in
#183's folders or crate is changed. A comparison, reported and not gated.
"""
import glob
import json
import math
import re
import subprocess
import sys

import mpmath as mp

mp.mp.dps = 40


def load(dump):
    """the SINGULAR records, each tagged with its file and line (#183's dump holds one point several times: first pass,
    pass of record, double- and quad-double runs)"""
    pts = []
    for path in sorted(glob.glob(f"{dump}/offgrid_*.jsonl")):
        for n, line in enumerate(open(path)):
            line = line.strip()
            if line.startswith("{"):
                r = json.loads(line)
                if r.get("class") == "SINGULAR" and r.get("est"):
                    r["_src"] = f"{path.rsplit('/', 1)[-1][8:18]}:{n}"
                    pts.append(r)
    return pts


def select(pts):
    """per point (F_inf, C, kind) the pass of record (the record with the smallest est.r) and, where one exists, the first
    pass (est.r at least 100 times larger); of the points: the 20 F->0 points of largest eta* (O8's failure regime), ten
    F->0 points at eta* about 1 ... 100, six a3->0 points"""
    groups = {}
    for p in pts:
        groups.setdefault((p["fi"], p["c"][0], p["kind"]), []).append(p)
    rec = {}
    for key, v in groups.items():
        v = sorted(v, key=lambda p: p["est"]["r"])
        first = v[-1] if v[-1]["est"]["r"] >= 100 * v[0]["est"]["r"] else None
        rec[key] = (v[0], first)
    f0 = sorted((k for k in rec if k[2] == "F->0"), key=lambda k: rec[k][0]["xstar"])
    big = f0[-20:]
    mid = [k for k in f0 if 1.0 <= rec[k][0]["xstar"] <= 100.0]
    targets = [10 ** (2.0 * i / 9) for i in range(10)]
    chosen, used = [], set()
    for t in targets:
        best = min((i for i in range(len(mid)) if i not in used), key=lambda i: abs(math.log(rec[mid[i]][0]["xstar"] / t)))
        used.add(best)
        chosen.append(mid[best])
    a3 = sorted((k for k in rec if k[2] == "a3->0"), key=lambda k: rec[k][0]["xstar"])
    step = max(1, len(a3) // 6)
    sel = []
    for tag, keys in (("O8", big), ("mid", chosen), ("a3", a3[::step][:6])):
        for k in keys:
            sel.append((tag, rec[k][0]))
            if rec[k][1] is not None:
                sel.append((tag + "-first", rec[k][1]))
    return sel


def pid(tag, p):
    return f"{tag}:{p['fi']:g}:{p['c'][0]:g}:{p['_src']}"


def harness_input(sel):
    lines = []
    for tag, p in sel:
        l = 2.0 ** round(math.log2(p["est"]["r"]))
        s = p["s_last"]
        lines.append(" ".join([pid(tag, p), repr(p["x_last"]), repr(l)] + [repr(v) for f in s for v in f]))
    return "\n".join(lines) + "\n"


def fmt_e(v):
    return "n/a" if v is None else f"{v:.6e}"


def ulp(x):
    """the spacing of f64 numbers at x: #183 stores eta* as an f64"""
    return 2.0 ** (math.floor(math.log2(abs(x))) - 52)


def dd(pair):
    return mp.mpf(pair[0]) + mp.mpf(pair[1])


ROW = re.compile(
    r"^\s*(\w+)\s+(F64|Dd|Qd)\s+(Real|Pair|Unresolved)\s+(\w+)\s+first\s+(\d+)\s+deg\s+(\d+)\s+\|\s+s: err (\S+) bar (\S+) tol (\S+)"
    r"\s+\|\s+alpha: err (\S+) bar (\S+) tol (\S+)(.*)$"
)


def planted_summary(path):
    """the planted suite's table (the test's printed rows, first occurrence of each cell), summarized"""
    rows, seen = [], set()
    for line in open(path):
        m = ROW.match(line)
        if not m or (m.group(1), m.group(2)) in seen:
            continue
        seen.add((m.group(1), m.group(2)))
        g = m.groups()
        rows.append(dict(id=g[0], lane=g[1], kind=g[2], model=g[3], es=float(g[6]), bs=float(g[7]), ts=float(g[8]),
                         ea=float(g[9]), ba=float(g[10]), ta=float(g[11]), fails="FAILS" in g[12], gated=g[0] != "N3"))
    gated = [r for r in rows if r["gated"]]
    print(f"planted suite {path.rsplit('/', 1)[-1]}: {len(rows)} rows, {len(gated)} gated cells (N3 is reported, not gated)")
    bad = [f"{r['id']} {r['lane']}" for r in gated if r["fails"]]
    print(f"  cells that fail a criterion: {len(bad)}: {bad}")
    print(f"  bar covers the error (criterion 3), s: {sum(r['es'] <= r['bs'] for r in gated)} of {len(gated)}, "
          f"alpha: {sum(r['ea'] <= r['ba'] for r in gated)} of {len(gated)}; all rows incl. N3: "
          f"{sum(r['es'] <= r['bs'] and r['ea'] <= r['ba'] for r in rows)} of {len(rows)}")
    for lane in ("F64", "Dd", "Qd"):
        v = [r for r in gated if r["lane"] == lane]
        ok = [r for r in v if not r["fails"]]
        ms = max(r["es"] / r["ts"] for r in ok)
        ma = max(r["ea"] / r["ta"] for r in ok)
        rs = sorted(r["bs"] / r["es"] for r in ok if r["es"] > 0)
        ra = sorted(r["ba"] / r["ea"] for r in ok if r["ea"] > 0)
        print(f"  {lane}: {len(ok)} of {len(v)} cells meet 1 to 4; largest error / tolerance: s {ms:.2g}, alpha {ma:.2g}; "
              f"bar / error: s median {rs[len(rs) // 2]:.3g} (min {rs[0]:.3g}, max {rs[-1]:.3g}), "
              f"alpha median {ra[len(ra) // 2]:.3g} (min {ra[0]:.3g}, max {ra[-1]:.3g})")
    print("  per cell: error of s (relative), error of alpha, bar of s (relative), bar of alpha; model and kind")
    for case in dict.fromkeys(r["id"] for r in rows):
        parts = []
        for lane in ("F64", "Dd", "Qd"):
            r = next(x for x in rows if x["id"] == case and x["lane"] == lane)
            parts.append(f"{lane} {r['es']:.1e} {r['ea']:.1e} {r['bs']:.1e} {r['ba']:.1e} {r['model'][:5]}/{r['kind']}")
        print(f"    {case:>3}: " + " | ".join(parts))
    kinds = {}
    for r in rows:
        kinds.setdefault((r["model"], r["kind"]), 0)
        kinds[(r["model"], r["kind"])] += 1
    print(f"  chosen model and kind over the rows: {kinds}")


def main():
    if len(sys.argv) > 3:
        planted_summary(sys.argv[3])
    dump, binary = sys.argv[1], sys.argv[2]
    pts = load(dump)
    sel = select(pts)
    print(f"#183 dump {dump.rstrip('/').rsplit('/', 1)[-1]} (negative-mass-raw/I183-rpi5): {len(pts)} SINGULAR records; {len(sel)} records: "
          + ", ".join(f"{t} {sum(x == t for x, _ in sel)}" for t in ("O8", "O8-first", "mid", "mid-first", "a3", "a3-first")))
    out = subprocess.run([binary, "30"], input=harness_input(sel), capture_output=True, text=True, check=True).stdout
    res = {}
    for line in out.splitlines():
        r = json.loads(line)
        res.setdefault(r["id"], {})[r["field"]] = r
    print("columns: point | #183: eta*, bar, fitted alpha (field), alpha_pred | locator, best field (smallest relative bar): "
          "kind, model, eta*, bar(eta*), alpha(F) +- bar, dEta* = eta*(here) - eta*(#183), in units of #183's bar and of this bar")
    stats = {}
    for tag, p in sel:
        i = pid(tag, p)
        l = 2.0 ** round(math.log2(p["est"]["r"]))
        best = None
        for f, r in res.get(i, {}).items():
            if r["kind"] in ("none",):
                continue
            if best is None or r["rel_bar"] < best[1]["rel_bar"]:
                best = (f, r)
        head = (f"{i:>16} | #183 eta* {p['xstar']:.15g} bar {p['bar']:.1e} alpha {p['est']['alpha']:.6e} (field {p['est']['field']}) "
                f"alpha_pred {fmt_e(p['alpha_pred'])}")
        if best is None:
            print(head + " | locator: no result")
            continue
        f, r = best
        x = mp.mpf(p["x_last"])
        eta = x + mp.mpf(l) * dd(r["re"])
        im = mp.mpf(l) * dd(r["im"])
        eta183 = mp.mpf(p["xstar"])
        a_field = dd(r["alpha"])
        a_F = a_field + f
        bar = l * r["bar"]
        d = eta - eta183
        print(head + f" | field {f} {r['kind']} {r['model']} first {r['first']} deg {r['degree']}: eta* {mp.nstr(eta, 15)} "
              f"{'+ i' + mp.nstr(im, 3) if r['kind'] == 'pair' else ''} bar {bar:.1e} alpha(F) {mp.nstr(a_F, 6)} +- {r['alpha_bar']:.1e} "
              f"dEta* {mp.nstr(d, 3)} = {mp.nstr(abs(d) / ulp(float(eta183)), 3)} ulp of f64 at eta* "
              f"({mp.nstr(abs(d) / max(p['bar'], 1e-300), 3) if p['bar'] > 0 else 'n/a'} x #183 bar, "
              f"{mp.nstr(abs(d) / bar, 3)} x this bar)")
        others = ", ".join(
            f"f{g}: {rr['kind']} rel_bar {rr['rel_bar']:.1e} alpha(F) {mp.nstr(dd(rr['alpha']) + g, 5)}"
            for g, rr in sorted(res[i].items()) if rr["kind"] != "none"
        )
        print("      all fields: " + others)
        st = stats.setdefault(tag, [])
        st.append(dict(kind=r["kind"], ulp=float(abs(d) / ulp(float(eta183))),
                       d183=float(abs(d) / p["bar"]) if p["bar"] > 0 else None,
                       a_here=a_F, a_pred=p["alpha_pred"], a_183=p["est"]["alpha"], a_bar=r["alpha_bar"]))
    print("summary per group (eta* from the locator against #183's f64 xstar; the exponent against alpha_pred = 2/(eta*^2 + 4) where it exists)")
    for tag, st in stats.items():
        us = sorted(x["ulp"] for x in st)
        line = (f"  {tag}: {len(st)} records, kinds {dict((k, sum(x['kind'] == k for x in st)) for k in sorted(set(x['kind'] for x in st)))}, "
                f"|d eta*| in f64 ulps: median {us[len(us) // 2]:.3g}, max {us[-1]:.3g}")
        pr = [x for x in st if x["a_pred"]]
        if pr:
            rh = [abs(x["a_here"] - x["a_pred"]) / x["a_pred"] for x in pr]
            r183 = [x["a_183"] / x["a_pred"] for x in pr]
            line += (f"; alpha: |here - pred| / pred max {mp.nstr(max(rh), 3)}, #183's fitted / pred from {min(r183):.4f} to {max(r183):.4f}, "
                     f"largest bar {max(x['a_bar'] for x in pr):.1e}")
        print(line)


if __name__ == "__main__":
    main()
