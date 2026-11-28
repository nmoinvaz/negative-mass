#!/usr/bin/env python3
"""Every number quoted in Paper IV (main.tex), read from the outputs of record in notes/black_universe and checked to be
printed verbatim in main.tex.

    python paper_numbers.py                 exit 1 if a number is not printed as its source gives it, or a stated relation fails
    python paper_numbers.py --mutate=NAME   the same with one source number changed as it is read (sources.MUTATIONS); must exit 1

Sources (all under notes/black_universe unless named): stability_output.txt (#212), bounce_output.txt (#223),
polar_output.txt (#218), evolve_output.txt (#217), linear_output.txt (#232), horizon_sector_output.txt (#248),
horizon_lemma_output.txt (#261), stability.py (the declared search range), and the coordinator's review sections of README.md,
polar.md, evolution.md and horizon_sector.md for the independent checks (ind212, ind218, ind217, ind248), cross-read against their
outputs in ../negative-mass-raw when that checkout is present. A last pass lists every decimal number in the body of main.tex
and fails if one is not among the checked strings.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sources as S  # noqa: E402

MUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--mutate=")), None)
S.set_mutation(MUT)
HERE = os.path.dirname(os.path.abspath(__file__))
TEX = open(os.path.join(HERE, "main.tex")).read()
TEXN = re.sub(r"\s+", " ", TEX)
bad = []
CHECKED = []                                       # every snippet confirmed in main.tex, for the coverage pass


def check(cond, what):
    if not cond:
        bad.append(what)
        print("  FAIL", what)


def need(snippet, what):
    sn = re.sub(r"\s+", " ", snippet)
    ok = sn in TEXN
    check(ok, f"{what}: '{snippet}' printed in main.tex")
    if ok:
        CHECKED.append(sn)


def sci(x, d=1, sign=False):
    """x as m\\times10^{e} with d decimals; from a string the mantissa digits are kept as printed when d is None."""
    if isinstance(x, str) and d is None:
        m, e = x.lower().split("e")
        return f"{m}\\times10^{{{int(e)}}}"
    x = float(x)
    e = math.floor(math.log10(abs(x)))
    m = f"{x / 10**e:.{d}f}"
    if m.lstrip("+-").startswith("10"):
        e += 1
        m = f"{x / 10**e:.{d}f}"
    if sign and x > 0:
        m = "+" + m
    return f"{m}\\times10^{{{e}}}"


def tnum(v):
    """a table entry: negative numbers in math mode."""
    return f"${v}$" if v.startswith("-") else v


# ---------------------------------------------------------------------------------------------------- #212, the counts
st = S.stability()
rows = st["rows"]
print(f"#212: {len(rows)} members, {st['one'][0]} of {st['one'][1]} with exactly one mode")
check(len(rows) == 20 and st["one"] == ("19", "20"), "20 members, 19 with one mode")
need("Nineteen of the 20 members have exactly one growing radial mode", "count")
need("The 19 members with $c\\ne-1$ have exactly one each", "abstract count")
need(f"at 20 members, $c={rows[0]['c']}$ to ${rows[-1]['c']}$", "range of c")
check(all(r["n1"] == r["n2"] for r in rows), "R1 and R2 counts equal")
check([r["c"] for r in rows if r["n2"] == 0] == ["-1"], "c = -1 alone has no mode")
check(all(r["nodes"] == "0" for r in rows if r["W"]), "no eigenfunction has a node")
Ws = [(float(r["c"]), float(r["W"])) for r in rows if r["W"]]
check(all(a[1] < b[1] for a, b in zip(Ws, Ws[1:])), "W falls monotonically with -c")
r12 = max(float(r["r12"]) for r in rows)
need(f"agree to ${sci(r12)}$ or better", "largest R1 v R2")
for r in rows:
    if r["W"]:
        line = (f"${r['c']}$ & {tnum(r['rho_h'])} & {r['kappa']} & {r['W']} & {r['W_kappa']} & {1 / float(r['W']):.3f} & "
                f"${sci(r['r12'])}$ \\\\")
    else:
        line = f"${r['c']}$ & {r['rho_h']} & {r['kappa']} & none$^a$ & -- & -- & -- \\\\"
    need(line, f"Table 1 row c = {r['c']}")
byc = {r["c"]: r for r in rows}
need(f"${1 / float(byc['-0.01']['W']):.3f}\\,b$ near the wormhole to ${1 / float(byc['-10']['W']):.3f}\\,b$", "e-folding range")
need(f"${float(byc['-10']['W']):.4f}/b$ at $c=-10$ to ${float(byc['-0.01']['W']):.4f}/b$ at $c=-0.01$", "abstract rate range")
k1 = byc["-1"]["kappa"]
check(abs(float(k1) - 2 / math.pi) < 1e-6, "kappa at c = -1 is 2/pi")
need(f"$\\kappa=2/\\pi={k1}$", "kappa at c = -1")
need(f"against $\\kappa={k1}$", "kappa at c = -1 (Sec. 4)")
need(f"($\\kappa={k1}$)", "kappa at c = -1 (Table 3)")
need(f"$W/\\kappa={byc['-1.001']['W_kappa']}$ and ${byc['-1.0001']['W_kappa']}$ at $c=-1.001$ and $-1.0001$, and "
     f"${byc['-0.9999']['W_kappa']}$ and ${byc['-0.999']['W_kappa']}$ at $c=-0.9999$ and $-0.999$", "W/kappa near c = -1")
w2, d2 = st["planted"]["R2"]
need(f"$W^2b^2={w2}$, against ${st['ref'][:19]}$", "planted Ellis value and reference")
need(f"a difference of ${sci(d2)}$", "planted Ellis difference")
need(f"to ${sci(st['pt']['S7b'])}$ and ${sci(st['pt']['S7d'])}$ (relative)", "Poeschl-Teller limits S7b, S7d")
need(f"numerically to ${sci(st['s8'])}$ at 30 digits", "S8")
need(f"${sci(st['s2'])}$ at 40 digits", "S2")
check(re.search(r"S8 Vhat of u = Y/rho against \[BKZ12\] \(40\), \(41\) \(mpmath 30 digits", st["text"]) is not None
      and "S2 background (mpmath 40 digits" in st["text"], "S8 at 30 digits, S2 at 40 digits")
check(len(st["controls"]) == 14, "14 controls fire")
need("Fourteen mutations, each fires", "controls")
spec = st["spec"]
check(len(spec) == 40 and st["spec_fail"] == 0 and all(x[3] == x[4] for x in spec), "spectral: 40 of 40, counts equal")
need("40 of 40, and no complex eigenvalue", "spectral 40 of 40")
uB = max((s[5] for s in spec if s[0] == "S-B" and s[1] == "u" and s[5]), key=float)
yAB = max((s[5] for s in spec if s[1] == "Y" and s[5]), key=float)
need(f"At $N=801$ its rates agree with the shooting to ${sci(uB, None)}$ for $c\\le-1$", "spectral u-form S-B worst")
need(f"they agree to ${sci(yAB, None)}$", "spectral throat members worst")
check("N = 401, X = 45; S-B N = 801, X = 90" in st["text"], "spectral resolutions in the output")
need("At $N=401$ and $801$ points", "spectral resolutions")
src = open(os.path.join(S.NB, "stability.py")).read()
lo, hi = re.search(r"W_LO, W_HI = ([\d.]+), ([\d.]+)", src).groups()
n1, n2 = re.findall(r"nW=(\d+)", src)
need(f"$W\\in[{lo}, {float(hi):g}]$ at {n1} (resolution R1) and {n2} (R2) log-spaced", "search range and grid")
# ind212 (review): the README's review section and the raw output
rd = S.read("README.md")
m = re.search(r"W equal to the shooting result to ([\d.e+-]+) or better", rd)
check(m is not None, "ind212 bound in the README review")
need(f"$W$ equal to the shooting to ${sci(m[1])}$", "ind212")
need("For the nine members with $c\\le-1$", "ind212 members")
raw212 = S.raw("I212-mac-review/ind212_run2.txt")
if raw212:
    rels = [float(x) for x in re.findall(r"\(rel ([\d.e+-]+)\)", raw212)]
    check(len(rels) == 8 and f"{max(rels):.1e}" == m[1] and "c=-1: negative eigenvalues [0, 0]" in raw212, "ind212 raw output")
else:
    print("  (negative-mass-raw not present: ind212 raw output not cross-read)")

# ---------------------------------------------------------------------------------------------------- #223, other sources
bo = S.bounce()
p2 = {}
for res, c, d in bo["P2"]:
    if d != "-":
        p2.setdefault(res, []).append(float(d))
check(len(set(c for _, c, _ in bo["P2"])) == 7, "seven zero-charge members recounted")
need(f"A recount of seven members on the wider range $W\\in[10^{{-5}},8]$", "recount")
check("low stage: W in [1e-05, 0.001]" in bo["text"] and "Main stage: W in [0.001, 8]" in bo["text"], "recount range in the output")
need(f"to ${sci(max(p2['R1']))}$ (R1) and ${sci(max(p2['R2']))}$ (R2)", "recount agreement")
check(bo["summary"] == ("54", "55"), "54 of 55")
need("54 of 55 members with a magnetic field", "abstract 54 of 55")
need("54 of 55 members have a growing radial mode", "Sec. 3.5 54 of 55")
mem = bo["members"]
grow = {k: (v["count"] > 0 or v["cplx"]) for k, v in mem.items()}
slices = {"S": "S_m", "B1": "B1_q", "B2": "B2_q", "B3": "B3_q", "P": "P2_c"}
tot = [0, 0]
for name, pre in slices.items():
    ks = [k for k in mem if k.startswith(pre)]
    n, g = len(ks), sum(grow[k] for k in ks)
    tot[0] += n
    tot[1] += g
    vals = sorted(float(k[len(pre):]) for k in ks)
    if name == "S":
        desc = f"Simpson--Visser, $m=M/a$ from {vals[0]:g} to {vals[-1]:g}"
    elif name == "P":
        desc = f"zero charge, $c={vals[0]:g}$, " + ", ".join(f"${v:g}$" for v in vals[1:])
    else:
        desc = None
    rowtail = f"& {n} & {g} \\\\"
    if desc:
        need(f"{name} & {desc} {rowtail}", f"Table 2 row {name}")
    else:
        need(f"$q$ from {vals[0]:g} to {vals[-1]:g} {rowtail}", f"Table 2 row {name}")
need(f"all & & {tot[0]} & {tot[1]} \\\\", "Table 2 total")
check(abs(float(mem["B1_q0"]["m"]) - 2 / math.pi) < 1e-4 and abs(float(mem["B3_q0"]["m"]) - 0.2 / math.pi) < 1e-4, "B1, B3 masses")
check(mem["B2_q0"]["count"] == 0 and all(mem[k]["count"] == 1 for k in mem if k.startswith("B2_q") and k != "B2_q0"), "B2 counts")
need(f"$W={mem['B2_q0.1']['main']}$ at $q=0.1$ rising to $W={mem['B2_q0.81']['main']}$ at $q=0.81$", "B2 rates")
cpl = [k for k in mem if mem[k]["cplx"]]
check(len(cpl) == 3 and sum(mem[k]["count"] == 0 for k in cpl) == 2, "three complex pairs, two without a real mode")
need("Three members have a complex pair, two of them without a real mode", "complex pairs")

# ---------------------------------------------------------------------------------------------------- #218, polar
po = S.polar()
nrun = 0
for ell, d in po["counts"].items():
    check(len(d["rows"]) == 21 and all(r[1] == "0" and r[2] == "0" and r[3] == "0" for r in d["rows"]), f"polar l = {ell}: none")
    check((d["lo"], d["hi"]) == (lo, "8"), f"polar range l = {ell}")
    nrun += 2 * len(d["rows"])
need(f"real $W\\in[{lo}, 8]$: {nrun} runs", "polar runs")
dmax = max(float(r[4]) for d in po["counts"].values() for r in d["rows"])
need(f"agree between R1 and R2 to ${sci(dmax)}$", "polar R1 v R2")
for ell, re_, im_, pre, pim in po["P2"]:
    need(f"${re_}{im_}i$", f"polar QNM l = {ell}")
    need(f"${pre}{pim}i$", f"[ABCKY23] QNM l = {ell}")
pm = S.read("polar.md")
m = re.search(r"ℓ = 2 1\.7376918464 − 0\.3051377413i, ([\d.e+-]+) from P2", pm)
check(m is not None, "ind218 l = 2 in polar.md review")
need(f"gives the $\\ell=2$ mode to ${sci(m[1])}$", "ind218")

# ---------------------------------------------------------------------------------------------------- #232, linear
li = S.linear()
R1 = li["R1"]
need(f"to ${sci(R1['-1.5'][2])}$ at $c=-1.5$, ${sci(R1['-0.3'][2])}$ at $c=-0.3$, ${sci(R1['-1.01'][2])}$ at $c=-1.01$ and "
     f"${sci(R1['-0.99'][2])}$ at $c=-0.99$", "linear solver against #212")
rate, kap, rel = li["R3"]
check(kap == k1, "kappa of R3 equals Table 1")
need(f"data G give ${rate}$ against $\\kappa={kap}$", "R3")
need(f"Sec.~\\ref{{sec:cm1}}, ${rate}$ against $\\kappa={kap}$", "R3 in the caption of Fig. 1")
c3 = re.findall(r"ratios ([\d.]+), ([\d.]+)\)", li["text"].split("[pass] C3:")[1].split("\n")[0])
check(all(x == "4.00" for t in c3 for x in t), "C3 ratios 4.00")
need("errors that fall by 4.00 per halving", "C3")
c5 = [x for t in re.findall(r"ratios ([\d.]+), ([\d.]+);", li["text"].split("[pass] C5:")[1].split("\n")[0]) for x in t]
need(f"fall by {min(c5)} to {max(c5)} per halving", "C5")
DESC = {"R": "pulse on the exterior ray (an ingoing wave)", "L": "pulse on the ray that crosses the horizon", "G": "both pulses",
        "A0": "corner $\\delta r$", "A1R": "corner $\\partial_p\\delta r$", "A1L": "corner $\\partial_\\sigma\\delta r$"}
for k, (yh, r_, lateexp, r1d, ratio) in li["classes"].items():
    yhs = "0" if float(yh) == 0 else f"${sci(yh, None)}$"
    need(f"{k} & {DESC[k]} & {yhs} & {r_} & {tnum(r1d)} \\\\", f"Table 3 row {k}")
b_, lexp, lrel = li["L3a"]
check(b_ == li["bkappa"] and abs(float(b_) - math.sqrt(2) * 2 / math.pi) < 1e-6, "b = sqrt2 kappa")
need(f"equals $\\sqrt2\\kappa={b_}$ to ${sci(lrel)}$", "L3a")
red, crate, _ = li["L3b"]
need(f"by ${float(red) / 1e4:.2f}\\times10^4$, and the combination decays, at rate ${crate}$", "L3b")
red3, rate3, W3 = li["L3c"]
need(f"only by ${float(red3):.1f}$ and leaves the rate at ${rate3}$, near $W={W3}$", "L3c")
check(f"{float(byc['-1.5']['W']):.4f}" == W3, "L3c W equals Table 1")
need(f"power ${li['L4'][0]}$", "L4")
need(f"to ${sci(li['L2'])}$ of its running maximum", "L2")
need(f"fall by {li['L5'][0]}, {li['L5'][1]} and {li['L5'][2]}, {li['L5'][3]} per halving", "L5")
drift = re.findall(r"drift ([\d.e+-]+) there", li["text"])
check(len(set(drift)) == 1, "L1 drift bar")
need(f"to ${sci(li['L1'])}$ until its own drift reaches ${sci(drift[0], 0)}$", "L1")
check(re.search(r"logarithm coefficient `π w\(1 − w\)/8`", S.read("linear.md")) is not None and "pi w (1 - w)/8" in li["derive_text"],
      "log coefficient in linear_derive_output")
need("with coefficient $\\pi w(1-w)/8$", "log coefficient")

# ---------------------------------------------------------------------------------------------------- #248, horizon sector
hs = S.horizon_sector()
CDESC = {"L": "scalar pulse on $V=-0.05$, centred on $U=0$", "G": "L and a pulse on $U=-2$", "A0": "corner $\\delta r$",
         "A1R": "corner $\\partial_V\\delta r$", "A1L": "corner $\\partial_U\\delta r$", "R": "pulse on $U=-2$ only"}
for k, v in hs["classes"].items():
    c1, c2 = v[0], v[1]
    c1s = sci(c1) if abs(float(c1)) < 1e-9 else sci(c1, None)
    c2s = f"${sci(c2, 0)}$ (seed)" if abs(float(c2)) < 1e-9 else f"${sci(c2, None)}$"
    need(f"{k} & {CDESC[k]} & ${c1s}$ & {c2s} \\\\", f"Table 4 row {k}")
check(abs(float(hs["classes"]["L"][1]) - float(hs["classes"]["G"][1])) < 1e-12, "c2 of L equals c2 of G")
need(f"has $c_2={sci(hs['classes']['L'][1], None)}$", "class L c2 in text")
hdr = re.search(r"rectangle V in \[(-?[\d.]+), (\d+)\], U in \[(-?\d+), ([\d.]+)\]; Richardson of h = ([\d.]+), ([\d.]+)", hs["text"]).groups()
need(f"evolves data on the rays $V={hdr[0]}$ and $U={hdr[2]}$", "Kruskal rays")
need(f"Richardson of $h={hdr[4]}$ and ${hdr[5]}$", "Kruskal resolutions")
check("crossing at V = 0, 0.5, 1" in hs["text"], "slices in the output")
need("at the crossings $V=0$, $0.5$ and $1$", "slices")
need(f"same $c_2$ to ${sci(hs['N4'])}$", "N4")
need(f"horizon fit to ${sci(hs['N5'])}$ after Richardson", "N5")
dev, lo2, hi2, med = hs["N11"]
check("for 20 random smooth data" in hs["text"], "20 random data")
need(f"Twenty random smooth data", "N11 count")
need(f"$|c_2|$ from ${sci(lo2, None)}$ to ${sci(hi2, None)}$, median ${sci(med, None)}$", "N11 range")
seeds = hs["N7"]
check(all(3.5 < float(x) < 4.5 for x in seeds[4:]), "class R seed falls by 4")
need(f"seed of ${sci(seeds[3])}$ at $h=0.00125$ that falls by 4 per halving", "N7")
coef, ex, mar = hs["N6"][:3]
need(f"(${sci(coef, None)}$), are recovered to ${sci(ex)}$ on exact data and ${sci(mar)}$ on the evolved solution", "N6")
n9 = max(float(x) for x in re.findall(r"relative ([\d.e+-]+)", hs["text"].split("[pass] N9:")[1].split("\n")[0]))
need(f"on the same data to ${sci(n9)}$", "N9")
check(re.search(r"amplitudes 1e-3", S.read("horizon_sector.md")) is not None, "class amplitude 1e-3 in the note")
need("for data of amplitude $10^{-3}$", "class amplitude")
hm = S.read("horizon_sector.md")
m = re.search(r"c₂ = (−[\d.]+)e-04, which agrees with the 2D Kruskal solver to ([\d.e+-]+),", hm)
check(m is not None, "ind248 in the horizon_sector.md review")
need(f"gives $c_2={m[1].replace('−', '-')}\\times10^{{-4}}$, which agrees to ${sci(m[2])}$", "ind248")
raw248 = S.raw("I248-mac-review/ind248_run1.txt")
if raw248:
    mr = re.search(r"K1: c2 = (-[\d.]+) against .*?relative ([\d.e+-]+);", raw248)
    check(mr is not None and abs(float(mr[1]) - float(m[1].replace("−", "-")) * 1e-4) < 1e-15
          and abs(float(mr[2]) - float(m[2])) <= 0.051e-8, "ind248 raw output (2.25e-8 printed as 2.3e-8 in the review)")

# ---------------------------------------------------------------------------------------------------- #261, the lemma
le = S.lemma()
check(le["checks"] == ("5", "5"), "lemma 5 of 5")
need("Five symbolic checks", "lemma checks")
L = {n.split(" (")[0]: (rp, k) for n, rp, k, _ in le["members"]}
rp15, k15 = L["c = -1.5 black universe"]
check(abs(float(L["c = -1 black universe"][1]) - 2 / math.pi) < 1e-9 and float(L["c = -1 black universe"][0]) == 0, "lemma c = -1")
need(f"black universe, $c=-1$ & 0 & $2/\\pi$ &", "Table 5 c = -1")
need(f"black universe, $c=-1.5$ & {rp15} & {k15} & pure gauge \\\\", "Table 5 c = -1.5")
check(f"{float(k15):.6f}" == byc["-1.5"]["kappa"], "lemma kappa at c = -1.5 equals Table 1")
rs, ks = L["Schwarzschild"]
need(f"Schwarzschild & {float(rs):g} & {float(ks):g} & pure gauge \\\\", "Table 5 Schwarzschild")
rsv, ksv = L["Simpson-Visser a = 2m"]
need(f"Simpson--Visser, $a=2m$ & {float(rsv):g} & {float(ksv):g} & extremal", "Table 5 Simpson-Visser")
need(f"$R'(\\rho_h)={rp15}$, and Sec.", "Appendix c = -1.5")
check("kappa - 2/pi = 0.0" in le["text"] and "mp.mp.dps = 30" in open(os.path.join(S.NB, "horizon_lemma.py")).read(),
      "kappa - 2/pi = 0 at 30 digits")
need("with $\\kappa-2/\\pi=0$ at 30 digits", "Appendix kappa at 30 digits")
check("at n = 2, 3" in le["text"], "L0 at n = 2, 3")
need("(checked symbolically at $n=2$ and $3$)", "L0")

# ---------------------------------------------------------------------------------------------------- #217, evolutions
ev = S.evolve()
worst = 0
for c in (-0.3, -1.5, -1.0):
    hs_ = sorted({h for (cc, h) in ev["runs"] if cc == c}, key=float)
    rels = [abs(float(ev["runs"][(c, h)]["rel"])) for h in hs_]
    check(all(a < b for a, b in zip(rels, rels[1:])), f"|relative| falls with h at c = {c}")
    for h in ("0.005", "0.0025"):
        r = ev["runs"][(c, h)]
        check(abs(float(r["slope"]) - float(r["ref"])) <= 2 * float(r["se"]), f"within 2 SE at c = {c}, h = {h}")
        check(float(r["lo"]) <= float(r["ref"]) <= float(r["hi"]), f"reference inside the spread at c = {c}, h = {h}")
        check(r["n"] == "18", "n = 18")
        need(f"${c:g}$ & {h} & ${r['slope']}\\pm{r['se']}$ & {r['ref']} & ${sci(r['rel'], 1, True)}$ & [{r['lo']}, {r['hi']}] \\\\",
             f"Table 6 row c = {c}, h = {h}")
        if h == "0.0025":
            worst = max(worst, abs(float(r["rel"])))
check(abs(1 / float(byc["-0.3"]["W"]) - float(ev["runs"][(-0.3, "0.005")]["ref"])) < 1e-5, "1/W at c = -0.3")
check(abs(math.pi / 2 - float(ev["runs"][(-1.0, "0.005")]["ref"])) < 1e-5, "1/kappa at c = -1")
need(f"and $1/\\kappa$ at $c=-1$, to ${sci(worst)}$", "abstract slope agreement")
o = ev["orders"]
need(f"orders {o['-0.3']} ($c=-0.3$), {o['-1.5']} ($c=-1.5$) and {o['-1']} ($c=-1$)", "A* orders")
need("$|A-A_*|\\in[10^{-7},10^{-5}]$, 9 points per side", "fit window")
evmd = S.read("evolution.md")
check("primary window |A − A\\*| ∈ [1e-7, 1e-5] (9 points per side)" in evmd and "15 fits" in evmd and "|Y| ≥ 2e-5" in evmd,
      "fit window, 15 fits and the 2e-5 threshold in evolution.md")
need("over 15 fits", "15 fits")
need("$|Y|$ on the diagonal reaches $2\\times10^{-5}$", "lifetime threshold")
need(f"at the rate ${ev['E6']}$", "E6")
xc = dict((c, (a, b)) for c, a, b in ev["XC"])
need(f"at $c=-0.3$ (${xc['-0.3'][0]}$ against ${xc['-0.3'][1]}$)", "XC c = -0.3")
need(f"at $c=-1.5$ (${xc['-1.5'][0]}$ against ${xc['-1.5'][1]}$)", "XC c = -1.5")
nm = {}
for blk in ev["text"].split("## evolve.py name")[1:]:
    c = re.search(r"--c=(-[\d.]+)", blk)[1]
    tE = re.search(r"A = A\* \+ 1e-4: .*?weighted E: t = ([\d.]+)", blk)
    nm[c] = tE[1] if tE else None
need(f"at $t={float(nm['-1.5']):.1f}$ ($c=-1.5$) and $t={float(nm['-1']):.1f}$ ($c=-1$)", "late anti-trapped cells")
check(float(nm["-1"]) <= min(float(x) for x in re.findall(r"weighted E: t = ([\d.]+) l = \+", ev["text"].split("--c=-1 ")[1].split("## ")[0])),
      "t = 50.9 is the earliest at c = -1")
pl = re.search(r"1/W inside the spread in (\d+) of (\d+), A\* recovered .*?1/\(sqrt2 W\) rejected in (\d+) of (\d+)", ev["text"]).groups()
need(f"recovered in {pl[0]} of {pl[1]} seeds, and one with the time unit off by $\\sqrt2$, rejected in {pl[2]} of {pl[3]}", "PL, PX")
check(ev["mass"]["-1.0"] == f"{2 / (3 * math.pi):.5f}", "mass at c = -1 is 2/(3 pi)")
need(f"the mass is $m={ev['mass']['-1.0']}$", "mass at c = -1")
m = re.search(r"`evolve\.F_np` agrees with them to ([\d.e+-]+) at 60 states", evmd)
check(m is not None, "ind217 in the evolution.md review")
need(f"brute-force Ricci tensor to ${sci(m[1])}$", "ind217")

# #262: the slope through c = -1 and I0 (notes/black_universe/cminus1_mode_output.txt)
CM = S.read("cminus1_mode_output.txt")
i0 = re.search(r"^\s+I0\s+= ([0-9.]+)", CM, re.M)
sl = re.search(r"d\(W/kappa\)/dc = -\(pi/2\) beta\*\s+= ([0-9.]+)", CM)
f3 = re.search(r"R2: s1 = [0-9.]+, s2 = [0-9.]+ \(second order\); s1 = ([0-9.]+),", CM)
check(i0 and sl and f3, "#262's I0, slope and third-order fit in cminus1_mode_output.txt")
need(f"$I_0={float(i0[1]):.10f}$", "I0")
need(f"$\\dd(W/\\kappa)/\\dd c={float(sl[1]):.10f}$ at $c=-1$", "d(W/kappa)/dc at c = -1")
need(f"a third-order fit of the rates gives ${float(f3[1]):.8f}$", "the third-order fit")

# #260: the count of papers searched (notes/black_universe/horizon_literature.md)
HL = S.read("horizon_literature.md")
npap = re.search(r"N = (\d+) papers whose full text was searched", HL)
check(npap is not None, "#260's paper count in horizon_literature.md")
need(f"in the {npap[1]} papers we searched", "papers searched (#260)")
nine = re.search(r"repeated without a new computation in 11 later papers", HL)
check(nine is not None, "#260: 11 later papers repeat BKZ12's verdict (Bro18, BBS20 and nine further)")
need("in~\\cite{Bro18,BBS20} and in nine further papers", "two cited and nine further = 11 (#260)")

# #258: the full-slice datum (notes/black_universe/full_slice_output.txt)
FS = S.read("full_slice_output.txt")
e4 = re.search(r"E4: .*?\(slice (-[0-9.]+e-04); raw fits [^;]*; Rich relative errors [0-9.e-]+, ([0-9.e-]+),", FS)
o1 = re.search(r"O1: two-sided .*?\(c2 slice (-[0-9.]+e-04) \(amplitude 1e-3\)", FS)
check(e4 is not None and o1 is not None, "#258's E4 and O1 in full_slice_output.txt")
need(f"$c_2={float(e4[1]) * 1e4:.6f}\\times10^{{-4}}$", "c2 of the full-slice datum")
need(f"to ${float(e4[2]) * 1e8:.1f}\\times10^{{-8}}$ (relative; Richardson of $h=0.005$ and $0.0025$)", "slice vs horizon fit")
need(f"(${float(o1[1]) * 1e4:.6f}\\times10^{{-4}}$ at amplitude $10^{{-3}}$)", "two-sided c2")

# #277: the exact family as planted truth for the nonlinear code (notes/black_universe/planted_family_output.txt)
PF = S.read("planted_family_output.txt")
rats = [float(x) for line in re.findall(r"error ratios per halving of h, r: ([0-9. ]+);", PF) for x in line.split()]
verr = [[abs(float(x)) for x in m.split(", ")] for m in re.findall(r"beta = -0\.05: v\* = [0-9.]+; found .*?; errors ([-+0-9.e, ]+); error ratios", PF)]
check(len(rats) == 16 and len(verr) == 2, f"#277: 16 error ratios and two v* error lists ({len(rats)}, {len(verr)})")
need(f"falls by ${min(rats):.3f}$ to ${max(rats):.3f}$ per halving of $h$ from $0.02$ to $0.00125$", "#277 error ratios")
need(f"found to ${max(e[0] for e in verr) * 1e4:.1f}\\times10^{{-4}}$ at the coarsest step and ${max(e[-1] for e in verr) * 1e7:.1f}\\times10^{{-7}}$ at the finest", "#277 v*")
need("at $c=-1.5$ and $-1$ with $\\beta=\\pm0.05$", "#277 members")

# #259: the static coefficient c1 (notes/black_universe/static_c1_output.txt)
SC = S.read("static_c1_output.txt")
t6 = [float(x) for x in re.findall(r"(?:A0|A1R|A1L): c1 [0-9.e+-]+, Theta ([0-9.e+-]+),", SC)]
t7 = [float(x) for x in re.findall(r"\|Theta\|max ([0-9.e+-]+) \(Richardson\)", SC)]
check(len(t6) == 3 and len(t7) == 3 and "[pass] T6:" in SC and "[pass] T7:" in SC, f"#259's T6 and T7 ({t6}, {t7})")
need(f"satisfy it to ${sci(max(t6))}$", "#259: kappa c1 = -kappa Theta_gen for A0, A1R, A1L")
need(f"have $\\Theta=0$ to ${sci(max(t7), d=0)}$", "#259: Theta of the c1 = 0 classes")
need("$M_{\\rm ADM}=2L^{3/2}/(3\\pi)$", "#259: the mass of the static member (ind259 J1)")

# #301: the static mass mode at other members (notes/black_universe/static_c1_members_output.txt)
SM = S.read("static_c1_members_output.txt")
t5 = re.search(r"c=-1\.5: r_S [0-9.]+, kappa c1 vs dM ([0-9.e+-]+) \(fails\); c=-0\.3: r_S [0-9.]+, kappa c1 vs dM ([0-9.e+-]+) \(fails\)", SM)
check(t5 is not None and "[pass] T5:" in SM and "[pass] T4:" in SM, "#301's T4 and T5")
need(f"misses by ${sci(t5[1], d=2)}$ and ${sci(t5[2], d=2)}$ if $r_S$ is set to\n$1$", "#301: the wrong-r_S control at c = -1.5, -0.3")
need("$M_{\\rm ADM}=-2cL^{3/2}/(3\\pi)$", "#301: the static member's mass at every c (ind301 J1)")

# #302: the decaying partner and the log of the exact modes at S (notes/black_universe/horizon_basis_output.txt)
HB = S.read("horizon_basis_output.txt")
l1 = [float(x) for x in re.findall(r"-> ([0-9.]+) \(Y_N\(0\)", HB)]
h3 = re.search(r"\[pass\] H3: .*?\(worst ([0-9.e+-]+);", HB)
check(len(l1) == 8 and h3 is not None and "[pass] L1:" in HB and "[pass] L2:" in HB and "[pass] S1:" in HB, f"#302's S1, L1 (8 classes), L2, H3 ({len(l1)})")
need(f"holds to ${sci(max(abs(x - 1) for x in l1), d=2)}$ on eight classes", "#302: the smoothness condition on evolved data (L1, worst Richardson ratio)")
need(f"to ${sci(h3[1])}$ of the amplitude (8 classes, 8 random data)", "#302: the three modes on both horizons (H3)")

# ---------------------------------------------------------------------------------------------------- placeholders and marks
pend = re.findall(r"\\pending\{\\#(\d+)", TEX)
needs = len(re.findall(r"% NEEDS #260", TEX))
print(f"\nplaceholders in the text: {', '.join('#' + p for p in pend)}; sentences marked NEEDS #260: {needs}")
check(not pend, "no placeholder remains (#262 and #258 filled)")

# ---------------------------------------------------------------------------------------------------- coverage
body = TEX.split("\\begin{abstract}")[1].split("\\begin{thebibliography}")[0]
body = re.sub(r"(?<!\\)%.*", "", body)
body = re.sub(r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}", "", body, flags=re.S)
body = re.sub(r"\\(ref|label|cite|eqref|includegraphics)(\[[^]]*\])?\{[^}]*\}", "", body)
ALLOWED = {}                                      # value: reason, for a decimal that is a definition, not a result
joined = " ".join(CHECKED)
decs = set(re.findall(r"(?<![\w.])\d+\.\d+", body))
got = set(re.findall(r"(?<![\w.])\d+\.\d+", joined))
missing = sorted(x for x in decs if x not in got and x not in ALLOWED)
print("allowed without a check:", {x: ALLOWED[x] for x in sorted(decs - got) if x in ALLOWED})
print(f"decimal numbers in the body: {len(decs)}; not among the checked strings: {missing}")
check(not missing, f"every decimal number in main.tex is checked (missing {missing})")

print(f"\n{len(CHECKED)} printed strings confirmed; {len(bad)} failed; mutation {MUT}")
sys.exit(1 if bad else 0)
