"""Every number quoted in Paper III (paper/3-wormhole-largeD), recomputed from the repo's data and the published tables.

    python3 paper_numbers.py            (prints; exit 1 if a stated check fails)

Inputs: growth.csv (profiles.py: #195's crate, qd, all 16 rows since #200 and #221); logs_order3.txt
(Omega1); Shinkai & Torii 2017 Table I (alpha_GB = 0 rows, transcribed below); Shinkai & Hayward 2002's lifetime
exponents 0.60 (ghost pulses) and 0.61 (normal pulses).
"""
import csv
import math
import os
import re
import sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
L2 = math.log(2)
OM0 = -4 * L2
OM1_CLOSED = 2 * math.pi**2 / 3 + 8 * L2**2
bad = 0


def check(cond, what):
    global bad
    if not cond:
        bad += 1
        print("  CHECK FAILED:", what)


# Omega1 from the matching integral (finest grid in logs_order3.txt)
om1 = float(re.findall(r"Omega1 = ([0-9.]+)", open(os.path.join(HERE, "logs_order3.txt")).read())[-1])
print(f"Omega1 (computed) = {om1:.10f};  2 pi^2/3 + 8 ln^2 2 = {OM1_CLOSED:.11f};  diff {om1 - OM1_CLOSED:.1e}")
check(abs(om1 - OM1_CLOSED) < 1e-8, "Omega1 closed form")
mp.mp.dps = 16
rel = mp.pslq([mp.mpf(repr(om1)), mp.pi**2, mp.log(2)**2, mp.log(2), 1], tol=mp.mpf("1e-9"), maxcoeff=60,
              maxsteps=10**6)
print("PSLQ relation [Omega1, pi^2, ln^2 2, ln 2, 1]:", rel)
check(rel == [3, -2, -24, 0, 0], "PSLQ relation")
mp.mp.dps = 30
# J = int_1^inf ln^2 R (R-1)^(-3/2) dR; with R = 1 + u^2 and one integration by parts, J = 8 int_0^inf ln(1+u^2)/(1+u^2) du
J = mp.quad(lambda R: mp.log(R)**2 / (R - 1)**mp.mpf(1.5), [1, 1.5, 2, 5, 10, 100, 1e4, mp.inf])
Ju = 8 * mp.quad(lambda u: mp.log(1 + u * u) / (1 + u * u), [0, 1, 10, mp.inf])
print(f"J = {mp.nstr(J, 15)} (direct), {mp.nstr(Ju, 15)} (u form);  8 pi ln 2 = {mp.nstr(8 * mp.pi * mp.log(2), 15)}")
check(abs(Ju - 8 * mp.pi * mp.log(2)) < 1e-12 and abs(J - Ju) < 1e-8, "J = 8 pi ln 2")

rows = list(csv.DictReader(open(os.path.join(HERE, "growth.csv"))))
print("\n  D    W^2 a^2 (data)      3-term formula   rel. diff    (data - 3 terms) s^2   source")
for r in rows:
    D = int(r["D"])
    w = float(r["W2a2"])
    s = 2 * (D - 3)
    f3 = s + OM0 + OM1_CLOSED / s
    print(f"{D:5d}  {w:18.10f}  {f3:15.8f}  {abs(w - f3) / w:9.2e}   {(w - f3) * s * s:12.2f}       {r['source']}")
    if 7 <= D <= 11:
        check(abs(w - f3) / w < 0.011, f"3-term within 1.1% at D = {D}")
# Paper III's expansion (#26, #208, #209: Omega_1 .. Omega_3 derived) and its Table 1, checked against main.tex
mp.mp.dps = 60
q, z3, PI = mp.log(2), mp.zeta(3), mp.pi
OMS = [-4 * q, 2 * PI**2 / 3 + 8 * q**2, 80 * z3 - mp.mpf(8) / 3 * PI**2 * q - mp.mpf(32) / 3 * q**3,
       22 * PI**4 / 45 + mp.mpf(16) / 3 * PI**2 * q**2 + mp.mpf(32) / 3 * q**4 - 320 * z3 * q - 384 * z3]
z5 = mp.zeta(5)
OM4 = (-10432 * z5 - 160 * PI**2 * z3 + 640 * q**2 * z3 - mp.mpf(64) / 9 * PI**2 * q**3 - mp.mpf(88) / 45 * PI**4 * q
       - mp.mpf(128) / 15 * q**5 - mp.mpf(2432) / 45 * PI**4 + 1536 * q * z3 + 3072 * z3)   # derived (#220, omega4.md)
OMS.append(OM4)
OM4_FIT = mp.mpf("-12794.499860023")              # #195's fit (w2_output.txt), 14 digits
o4 = open(os.path.join(HERE, "order4_output.txt")).read()
o5 = open(os.path.join(HERE, "order5_output.txt")).read()
a70 = mp.mpf(re.search(r"X1: Omega2 \(A70\) = ([0-9.]+)", o4).group(1))
check(abs(a70 - OMS[2]) < mp.mpf("1e-45"), "Omega2 closed form against order4.py's A70")
m5 = re.search(r"= (-652\.[0-9]{40,})", o5)
check(m5 is not None and abs(mp.mpf(m5.group(1)) - OMS[3]) < mp.mpf("1e-38"), "Omega3 closed form against order5_output.txt")
TEX = open(os.path.join(HERE, "..", "..", "paper", "3-wormhole-largeD", "main.tex")).read()
quoted = {"Omega1": mp.nstr(OMS[1], 12), "Omega2": mp.nstr(OMS[2], 12), "Omega3": mp.nstr(OMS[3], 12),
          "Omega1/2": mp.nstr(OMS[1] / 2, 10), "Omega2/4": mp.nstr(OMS[2] / 4, 10), "Omega3/8": mp.nstr(OMS[3] / 8, 10),
          "Omega4": mp.nstr(OM4, 10), "Omega4/16": mp.nstr(OM4 / 16, 10)}
ratios = [f"{abs(float(OMS[k + 1] / OMS[k])):.1f}" for k in range(4)]
o6 = open(os.path.join(HERE, "order6_output.txt")).read()
m6 = re.search(r"Omega4 = (-12794\.[0-9]{60,})", o6)
check(m6 is not None and abs(mp.mpf(m6.group(1)) - OM4) < mp.mpf("1e-55"), "Omega4 closed form against order6_output.txt")
check(abs(OM4 - OM4_FIT) <= mp.mpf("5e-10"), "Omega4 closed form against #195's 14-digit fit")
quoted["ratios"] = ", ".join(ratios[:3]) + " and " + ratios[3]
dB = float(re.search(r"M3: Omega2 of A70 against B2: ([0-9.]+) digits", o4).group(1))
dC = re.search(r"X3: Omega3 against.*?digits: B3 ([0-9.]+), B4", o5)
d6 = re.search(r"M4: Omega4 of G1 against G2.*?: ([0-9.]+) digits", o6)
quoted["Chebyshev digits"] = f"{int(dB)}, {int(float(dC.group(1)))} and {int(float(d6.group(1)))}" if dC and d6 else "?"
crate_rows = {int(m.group(1)): (mp.mpf(m.group(3)), mp.mpf(m.group(4))) for m in re.finditer(
    r"^\s+(\d+)\s+(TS13 Table I|modes_r\.py)\s+([0-9.]+)\s+([0-9.]+)", open(os.path.join(HERE, "w2_output.txt")).read(), re.M)}
low = max(crate_rows[D][1] - crate_rows[D][0] for D in (20, 30, 50, 100, 200, 500, 1000))
quoted["r-integration low by"] = f"{float(low):.1e}".replace("e-06", "\\times10^{-6}")
ts_digits = min(-mp.log10(abs(crate_rows[D][1] - crate_rows[D][0]) / crate_rows[D][1]) for D in range(4, 12))
check(ts_digits >= 9, f"TS13 agrees with the crate to nine or more significant figures ({float(ts_digits):.2f})")
print("\n  quoted coefficients:", quoted)
for k, v in quoted.items():
    check(v in TEX, f"{k} = {v} printed in main.tex")


def sci(x):
    x = abs(float(x))
    if x >= 0.1:
        return f"${x:.2g}$"
    e = math.floor(math.log10(x))
    return f"${x / 10**e:.1f}\\times10^{{{e}}}$"


print("\n  Table 1 rows (as printed in main.tex)")
TROWS = (4, 5, 6, 7, 8, 10, 20, 50, 100, 200, 500, 1000, 2000)
rowsD = {int(r["D"]): mp.mpf(r["W2a2"]) for r in rows}
w2o = open(os.path.join(HERE, "w2_output.txt")).read()
sec = w2o[w2o.index("growth.csv against the crate"):]
full = {int(m.group(1)): mp.mpf(m.group(2)) for m in                       # the crate's qd values, about 20 digits
        re.finditer(r"^\s+(\d+)\s+(?:TS13 Table I|modes_r\.py)\s+[0-9.]+\s+([0-9.]+)", sec, re.M)}
for D in TROWS:
    w = rowsD[D]
    s_ = mp.mpf(2 * (D - 3))
    f3 = s_ + OMS[0] + OMS[1] / s_
    f5 = f3 + OMS[2] / s_**2 + OMS[3] / s_**3
    r3 = (w - f3) * s_**2
    f5full = f5                                                             # the s^4 column needs the full digits
    r5 = f"${float((full[D] - f5full) * s_**4):.0f}$" if D >= 100 else "--"
    line = (f"{D} & {mp.nstr(w, len(str(int(w))) + 10, strip_zeros=False)} & {sci((w - f3) / w)} & {sci((w - f5) / w)} & "
            f"${float(r3):.2f}$ & {r5} \\\\")
    print("   ", line)
    check(line in TEX, f"Table 1 row D = {D} printed in main.tex")
prev = None
for D in (20, 30, 50, 100, 200, 500, 1000, 2000):                  # the five-term error falls monotonically from D = 20
    s_ = mp.mpf(2 * (D - 3))
    e5 = abs(rowsD[D] - (s_ + OMS[0] + OMS[1] / s_ + OMS[2] / s_**2 + OMS[3] / s_**3)) / rowsD[D]
    check(prev is None or e5 < prev, f"five-term error falls at D = {D}")
    prev = e5
print("  (data - 3 terms) s^2 tends to Omega2 + Omega3/s: at D = 1000 that is "
      f"{float(OMS[2] + OMS[3] / 1994):.2f}; (data - 5 terms) s^4 tends to Omega4 = {float(OM4):.5f} (derived)")
s_ = mp.mpf(2 * (2000 - 3))
O5est = (full[2000] - (s_ + sum(OMS[k] / s_**k for k in range(5)))) * s_**5
o5txt = f"{float(O5est):.0e}".replace("e+04", "\\times10^4")      # one figure: an estimate
print(f"  (data - 6 terms) s^5 at D = 2000: {float(O5est):.4g}, an estimate of Omega5; printed as {o5txt}")
check(f"$\\Omega_5\\approx{o5txt}$" in TEX, f"the Omega5 estimate {o5txt} printed in main.tex")

# the ghost sector of Choptuik's solution (notes/wormhole_next, #30): mu(D), as quoted in Sec. 7
import json
mu = {}
for fn in ("ghost_mu.jsonl", "ghost_mu2.jsonl"):
    for line in open(os.path.join(HERE, "..", "wormhole_next", fn)):
        r = json.loads(line)
        if r["mode"] == 1 and r.get("plant") is None and r["pivot0"] / r["pivot_h"] < 1e-4:
            z = [x["mu"] for x in r["zeros"] if x["ok"]]
            if z:
                mu[round(r["D"], 4)] = max(z)
mq = {D: f"{mu[D]:.4f}" for D in (4.0, 5.0, 70.5)}
print("\n  mu(D) as quoted:", mq)
for D, v in mq.items():
    check(v.replace("-", "−") in TEX or v in TEX, f"mu(D = {D}) = {v} printed in main.tex")

# e-folding times and the type I lifetime slope 1/(sqrt 2 W) in the null coordinate x^- = (t - r_*)/sqrt 2
TS = {int(r["D"]): float(r["W2a2"]) for r in rows}
print("\n  D   1/W (e-folding, units of a)   predicted slope 1/(sqrt2 W)")
for D in (4, 5, 6, 100):
    W = math.sqrt(TS[D])
    print(f"{D:4d}   {1 / W:.4f}                        {1 / (math.sqrt(2) * W):.4f}")
# Shinkai & Torii 2017, Table I, alpha_GB = 0 rows: (D, Delta E / a0, x_H^- / a0), two amplitudes per D
ST17 = {4: [(0.03, 2.94), (0.10, 2.09)], 5: [(0.15, 2.26), (0.61, 1.66)], 6: [(0.38, 1.92), (1.50, 1.46)]}
print("\n  D   ST17 two-run slope  (+/- from the table's rounding: dE +/- 0.005, x +/- 0.005)   predicted")
for D, ((e1, x1), (e2, x2)) in ST17.items():
    lr = math.log(e2 / e1)
    sl = (x1 - x2) / lr
    dlr = 0.005 / e1 + 0.005 / e2
    err = 0.01 / lr + sl * dlr / lr
    pred = 1 / (math.sqrt(2) * math.sqrt(TS[D]))
    print(f"{D:4d}   {sl:.3f} +/- {err:.3f}                                                       {pred:.3f}"
          f"   ({(sl - pred) / err:+.1f} sigma)")
    check(abs(sl - pred) < 2 * err, f"ST17 slope within 2 sigma at D = {D}")
print("\n  D = 4, Shinkai & Hayward 2002 small-amplitude runs: 0.60 (ghost pulses), 0.61 (normal pulses);"
      f" predicted {1 / (math.sqrt(2) * math.sqrt(TS[4])):.3f}")

# the threshold bisections of Sec. 6 (notes/wormhole_evolution, #40): threshold_output.txt (first code) and
# threshold_check_output.txt (the coordinator's second code, with #40's points refit over [1e-9, 1e-4])
import re
EV = os.path.join(HERE, "..", "wormhole_evolution")
sl, rel = {}, {}
for line in open(os.path.join(EV, "threshold_output.txt")):
    m = re.match(r"\s+SL D=(\d) (\w) h = ([0-9.]+): .*primary slope ([0-9.]+) .*rel dev ([-+0-9.e]+)", line)
    if m:
        sl[(int(m[1]), m[2], float(m[3]))] = (m[4], float(m[5]))
k4 = [l for l in open(os.path.join(EV, "threshold_check_output.txt")) if l.startswith("[PASS] K4")]
check(len(k4) == 1, "threshold_check_output.txt has a passing K4 line")
k4 = dict((int(D), (a, b)) for D, a, b in re.findall(r"D=(\d): ind40 ([0-9.]+), #40 refit over \[1e-9, 1e-4\] ([0-9.]+)", k4[0]))
FINE = {4: 0.005, 5: 0.00176777, 6: 0.00125}
print("\n  D   1/W        first code [1e-11,1e-3]   first code [1e-9,1e-4]   second code [1e-9,1e-4]   largest |rel| [1e-9,1e-4]")
bound = {}
for D in (4, 5, 6):
    Wi = 1 / math.sqrt(TS[D])
    first, second = float(k4[D][1]), float(k4[D][0])
    worst = max(abs(first * math.sqrt(TS[D]) - 1), abs(second * math.sqrt(TS[D]) - 1))
    e = math.floor(math.log10(worst))
    bound[D] = (math.ceil(worst / 10 ** e), e)
    row = f"{D} & {Wi:.6f} & {sl[(D, 'S', FINE[D])][0]} & {k4[D][1]} & {k4[D][0]} \\\\"
    print(f"{D:4d}   {Wi:.6f}   {sl[(D, 'S', FINE[D])][0]}                  {k4[D][1]}                 {k4[D][0]}                  {worst:.1e}")
    check(row in TEX, f"Table bisect row D = {D} printed in main.tex: {row}")
txt = ", ".join(f"${m}\\times10^{{{e}}}$" for m, e in (bound[4], bound[5])) + f" and ${bound[6][0]}\\times10^{{{bound[6][1]}}}$"
check(txt in TEX, f"the bounds {txt} printed in main.tex")
check(f"within ${bound[6][0]}\\times10^{{{bound[6][1]}}}$ in $D=4$, $5$ and $6$" in TEX, "the abstract's bound printed in main.tex")
low = (float(sl[(6, 'S', 0.00125)][0]) / float(k4[6][1]) - 1) * 100
check(f"lowers the $D=6$ slope by ${abs(low):.1f}\\%$" in TEX, f"D = 6 window effect {low:.2f} % printed in main.tex")
GH = {4: 0.005, 5: 0.00176777, 6: 0.0025}
gp = [f"${round(100 * sl[(D, 'G', GH[D])][1])}\\%$" for D in (4, 5, 6)]
check(f"{gp[0]}, {gp[1]} and {gp[2]} above $1/W$" in TEX, f"G family departures {gp} printed in main.tex")

# the expanding side (notes/wormhole_evolution/expansion_output.txt, #213): k of H = k/(tau_s - tau) for the tuned members
blk, kk = None, {}
for line in open(os.path.join(EV, "expansion_output.txt")):
    m = re.match(r"\s+--- D(\d), h = ([0-9.]+), member (\S+)", line)
    if m:
        blk = (int(m[1]), float(m[2]), m[3])
    m = re.search(r"1/H linear in tau: k = ([0-9.]+)", line)
    if m and blk and blk[2] == "tuned":
        kk[blk[:2]] = m[1]
print("\n  k of the tuned members (D, h): ", kk)
k4 = {v for (D, _), v in kk.items() if D == 4}
k5 = {f"{float(v):.2f}" for (D, _), v in kk.items() if D == 5}
check(k4 == {"2.0000"} and "$k=2.0000$ in $D=4$" in TEX, f"k at D = 4 is {k4} and is printed in main.tex")
check(len(k5) == 1 and f"$k\\approx{k5.pop()}$ in $D=5$" in TEX, "k at D = 5 printed in main.tex")

# the settled threshold mass (notes/wormhole_evolution/settled_mass_output.txt, #41): plateau Richardson value and r_AH
SM = open(os.path.join(EV, "settled_mass_output.txt")).read()
rich = [float(x) for x in re.findall(r"h = \(0\.01, 0\.005, 0\.0025\): M_end \[[^]]*\], differences [^,]*, [^,]*, Richardson \(second order\) M_end\(0\) = ([0-9.]+)", SM)]
rah = [float(x) for x in re.findall(r"^\s+9\.40\s+\S+\s+\S+\s+\S+\s+([0-9.]+)\s+0\.22308", SM, re.M)]
print("\n  settled mass (Richardson, finest triple):", rich, " r_AH at the limit:", rah)
check(len(rich) == 2 and all(f"${x:.4f}a$" in TEX for x in rich), f"the settled mass {rich} printed in main.tex as 0.2231a")
check(len(rah) >= 1 and f"radius ${rah[0]:.3f}a$" in TEX, f"r_AH {rah} printed in main.tex")
check("falls from $0.494a$ at formation" in TEX and "0.49424" in SM, "the formation value 0.494 a printed and in the output")
# family G (#253): settled mass, gap to S, and the per-decade rise for p* - p < 1e-6 (notes/wormhole_evolution/settled_mass_g_output.txt)
SG = open(os.path.join(EV, "settled_mass_g_output.txt")).read()
plg = [float(x) for x in re.findall(r"G: plateau mean over the last 1\.0 of x\+ at h = 0\.005, delta = 1e-8: ([0-9.]+)", SG)]
gap = re.findall(r"delta = 1e-08, h = \(0\.02, 0\.01, 0\.005\):[^\n]*\n\s+S plateau[^\n]*G - S at the same h: [^,]*, [^,]*, \+([0-9.e-]+)", SG)
tab = SG.split("family G, h = 0.005")[1].split("\n\n")[0]
pl = {r[0]: float(r[1]) for r in re.findall(r"^\s+(1e-\d\d)\s+\S+\s+\S+\s+\S+\s+([0-9.]+)", tab, re.M)}
inc = [pl[f"1e-{k + 1:02d}"] - pl[f"1e-{k:02d}"] for k in range(6, 12)]
print("\n  family G: plateau", plg, " G - S at h = 0.005:", gap, " per-decade rise 1e-6 .. 1e-12:", [f"{x:.2e}" for x in inc])
check(len(plg) == 1 and f"${plg[0]:.4f}a$ at $p=p_*-10^{{-8}}$" in TEX, f"G's settled mass {plg} printed in main.tex")
check(len(gap) == 1 and f"${float(gap[0]):.4f}a$ above the mirrored family" in TEX, f"G - S {gap} printed in main.tex")
check(len(inc) == 6 and f"{min(inc):.0e}" == "2e-04" and f"{max(inc):.0e}" == "6e-04" and "$2\\times10^{-4}$ to $6\\times10^{-4}$ per decade" in TEX,
      f"the per-decade rise {min(inc):.1e} .. {max(inc):.1e} printed as 2e-4 to 6e-4")
# the one-sided family's locus slope and the stall member's e-fold at D = 4 (#225, throat_locus_output.txt)
TL = open(os.path.join(EV, "throat_locus_output.txt")).read()
ra = re.findall(r"RA\[D=4 G\]: [^(]*\(([0-9.]+) \(h = 0\.01\), ([0-9.]+) \(h = 0\.005\)", TL)
mem = re.findall(r"CMP1\[D=4 G h=0\.005\][^\n]*?\(t = 12\.0: 1/lambda [0-9.]+, member ([0-9.]+)", TL)
print("  locus slope at D = 4 (h = 0.01, 0.005):", ra, " member e-fold at t = 12:", mem)
check(len(ra) == 1 and all(f"{float(x):.4f}" == "0.8495" for x in ra[0]) and "slope in $D=4$ is $0.8495$" in TEX, "the locus slope 0.8495 printed in main.tex")
check(len(mem) == 1 and f"{float(mem[0]):.3f}" == "0.847" and "e-folding time $0.847$" in TEX, "the stall member's e-fold 0.847 printed in main.tex")
print("\nchecks failed:", bad)
sys.exit(1 if bad else 0)
