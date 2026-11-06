#!/usr/bin/env python3
"""Tables and checks of step 2 (notes/wormhole_next.md section 6), from counts/modecount_test_results.jsonl,
axis_s1.jsonl, ghost_mu.jsonl and ghost_mu2.jsonl.

    python3 step2_table.py

Checks (exit 1 if one fails):
  C1  every plain count: each box resolved (every phase step below pi/4, no verification mismatch) and its winding of
      g is 0, so the only zero of det M_psi in the counted region is the expected mu = 0.
  C2  every planted count: the box holding the planted pair counts 2 more (B0 is symmetric about the real axis and
      holds both members), every other box the same as the plain count. A count without its planted run fails C2.
  C3  the other parity sector: every zero found on the line Im mu = omega has Re mu below the real-axis zero of the
      same state (so the real zero is the least damped on the two lines).
Printed, not checked: the real-axis zeros mu(D) of ghost_mu2.jsonl and the layout and resolution pairs.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
bad = []


def check(name, ok, msg):
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {msg}")
    if not ok:
        bad.append(name)


def rows(path):
    return [json.loads(l) for l in open(path)] if os.path.exists(path) else []


def main():
    counts = rows(os.path.join(HERE, "counts", "modecount_test_results.jsonl"))
    plain = {r["tag"]: r for r in counts if r.get("plant") is None}
    planted = {r["tag"]: r for r in counts if r.get("plant") is not None}
    print("counts (region: aLL < Re mu < aR, |Im mu| < 1.5 omega unless only B0, U0)")
    for tag, r in sorted(plain.items(), key=lambda kv: kv[1]["D"]):
        bx = r["boxes"]
        reg = f"Re ({r['aLL']:.3f}, {max(b['rect'][1] for b in bx.values()):.3f}), |Im| < {r['y2']:.3f}"
        w = {b: v["winding_g_int"] for b, v in bx.items()}
        res = all(v["resolved"] for v in bx.values())
        step = max(v["max_step"] for v in bx.values())
        print(f"  D {r['D']:8.4f} M {r['M']}: {reg}; windings {w}; max step {step:.3f}; {r['evaluations']} evaluations")
        check(f"C1 {tag}", res and all(v == 0 for v in w.values()), f"resolved {res}, windings {w}")
        pk = [k for k in planted if k.startswith(tag)]
        if not pk:
            check(f"C2 {tag}", False, "no planted run")
            continue
        p = planted[pk[0]]
        z = complex(*p["plant"])
        inside = [b for b, v in p["boxes"].items()
                  if v["rect"][0] < z.real < v["rect"][1] and v["rect"][2] < z.imag < v["rect"][3]]
        ok = all(p["boxes"][b]["winding_g_int"] == bx[b]["winding_g_int"] + (2 if b in inside and b in ("B0", "full")
                 else 1 if b in inside else 0) for b in bx)
        check(f"C2 {tag}", ok and p["boxes"][inside[0]]["resolved"] if inside else False,
              f"plant {z}: windings {({b: v['winding_g_int'] for b, v in p['boxes'].items()})}, box holding it {inside}")
    mus = {}
    for f in ("ghost_mu.jsonl", "ghost_mu2.jsonl"):
        for r in rows(os.path.join(HERE, f)):
            if r["mode"] == 1 and r.get("plant") is None:
                z = [x["mu"] for x in r["zeros"] if x["ok"]]
                if z:
                    mus.setdefault(r["state"], []).append((r["D"], r["M"], max(z), r["tag"]))
    print("other parity sector, Im mu = omega (zeros on the line, Re mu from lo to hi)")
    for r in rows(os.path.join(HERE, "axis_s1.jsonl")):
        zs = [z["re"] for z in r["zeros"] if "re" in z]
        real = [m for _, _, m, _ in mus.get(r["state"], [])]
        ref = max(real) if real else None
        print(f"  {r['state']}: D {r['D']:.4f}, Re {r['lo']} .. {r['hi']}: zeros at Re {[round(z, 6) for z in zs]}; "
              f"real-axis zero {ref}")
        if ref is not None:
            check(f"C3 {r['state']}", all(z < ref for z in zs), f"line zeros {[round(z, 4) for z in zs]} below {ref:.4f}")
    print("real-axis zeros of step 2 (ghost_mu2.jsonl)")
    for r in rows(os.path.join(HERE, "ghost_mu2.jsonl")):
        z = [x["mu"] for x in r["zeros"] if x["ok"]]
        print(f"  {r['tag']:18s} D {r['D']:8.4f} M {r['M']} K {r['K']} |F| {r['F']:.1e}: zeros {[round(x, 10) for x in z]}"
              f" in [{r['lo']}, {r['hi']}]")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
