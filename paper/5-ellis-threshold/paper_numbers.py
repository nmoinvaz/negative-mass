#!/usr/bin/env python3
"""Every number quoted in Paper V (paper/5-ellis-threshold/main.tex), read from the outputs of record and checked
against the text.

    python3 paper_numbers.py [--raw=~/Source/negative-mass-raw] [--mutate=NAME]

Each number is read from the output file of the script that produced it, formatted as the paper prints it, and looked
up in main.tex (comments removed, whitespace collapsed). A stated agreement that does not hold is a failure too. Exit 1
on any failure. --mutate=NAME changes one source number before it is read; every mutation must make the script fail.

Sources, notes/wormhole_evolution/ (this repository):
  threshold_output.txt (#40, threshold.py), threshold_check_output.txt (the review's second code, threshold_check.py),
  throat_locus_output.txt (#225), settled_mass_output.txt (#41), settled_mass_g_output.txt (#253),
  expansion_output.txt (#213), constraint_plateau_output.txt (#226), derive_output.txt (#39), threshold.py (parameters)
  and notes/wormhole_largeD/w2_output.txt (#195's W^2 a^2).
Sources, negative-mass-raw (the raw-output repository, read only):
  I39-dh610/runs/run01_checks.log and run01_mut_*.log (#39's checks.py), I39-mac-review (ind39),
  I40-mac-review/ind40_run1_mut_time-unit.txt, I41-mac-review (ind41, the rerun), I213-mac-review (ind213),
  I225-mac-review (ind225), I226-mac-review (ind226), I226-rpi5/supplement (roundoff of r), I253-mac-review (rerun).
Transcribed from the papers themselves: [SH02]'s exponents and masses, [GGS09b]'s radius and e-fold, [ST17]'s Table I.
"""
import argparse
import glob
import math
import os
import re
import sys
from decimal import Decimal

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, "..", "..", "notes", "wormhole_evolution")
LD = os.path.join(HERE, "..", "..", "notes", "wormhole_largeD")

# name: (source key, text in the source, replacement); each must make at least one check fail
MUTATIONS = {
    "k_D4": ("expansion", "1/H linear in tau: k = 2.0000, tau_s = 17.8902", "1/H linear in tau: k = 2.0001, tau_s = 17.8902"),
    "settled_G": ("settled_g", "delta = 1e-8: 0.235548", "delta = 1e-8: 0.235648"),
    "member_M12": ("member", "      12  0.237622   0.047732", "      12  0.237722   0.047732"),
    "spectrum_zero": ("spectrum", "zero at W = -2.61874 +0.15269i", "zero at W = -1.61874 +0.15269i"),
    "slope_S6": ("check", "ind40 0.461705", "ind40 0.461805"),
}

ap = argparse.ArgumentParser()
ap.add_argument("--raw", default=os.path.expanduser("~/Source/negative-mass-raw"))
ap.add_argument("--mutate", default=None, choices=sorted(MUTATIONS))
OPT = ap.parse_args()
RAW = OPT.raw


def read(*p):
    return open(os.path.join(*p)).read()


SRC = {
    "threshold": read(EV, "threshold_output.txt"),
    "check": read(EV, "threshold_check_output.txt"),
    "locus": read(EV, "throat_locus_output.txt"),
    "settled": read(EV, "settled_mass_output.txt"),
    "settled_g": read(EV, "settled_mass_g_output.txt"),
    "drift": read(EV, "throat_drift_output.txt"),
    "scaling": read(EV, "throat_scaling_output.txt"),
    "massb": read(EV, "settled_mass_b_output.txt"),
    "member": read(EV, "settled_mass_member_output.txt"),
    "spectrum": read(EV, "ellis_spectrum_output.txt"),
    "expansion": read(EV, "expansion_output.txt"),
    "plateau": read(EV, "constraint_plateau_output.txt"),
    "derive": read(EV, "derive_output.txt"),
    "thr_py": read(EV, "threshold.py"),
    "w2": read(LD, "w2_output.txt"),
    "checks39": read(RAW, "I39-dh610", "runs", "run01_checks.log"),
    "ind39_1": read(RAW, "I39-mac-review", "ind39_run1.log"),
    "ind39_2": read(RAW, "I39-mac-review", "ind39_run2_k2b.txt"),
    "ind40_tu": read(RAW, "I40-mac-review", "ind40_run1_mut_time-unit.txt"),
    "ind41": read(RAW, "I41-mac-review", "ind41_run2_K1b.txt"),
    "rerun41": read(RAW, "I41-mac-review", "rerun_job_compare.txt"),
    "ind213": read(RAW, "I213-mac-review", "ind213_run1.txt"),
    "ind225": read(RAW, "I225-mac-review", "ind225_run1.txt"),
    "ind226": read(RAW, "I226-mac-review", "ind226_run1.txt"),
    "r_round": read(RAW, "I226-rpi5", "supplement", "r_roundoff_static_D6.log"),
    "rerun253": read(RAW, "I253-mac-review", "compare_jobs.txt"),
}
N_MUT39 = len(glob.glob(os.path.join(RAW, "I39-dh610", "runs", "run01_mut_*.log")))
if OPT.mutate:
    key, old, new = MUTATIONS[OPT.mutate]
    assert old in SRC[key], f"mutation {OPT.mutate}: '{old}' not in {key}"
    SRC[key] = SRC[key].replace(old, new)

TEX = read(HERE, "main.tex")
TEXN = re.sub(r"\s+", " ", "\n".join(re.sub(r"(?<!\\)%.*", "", line) for line in TEX.splitlines()))
bad = 0
count = 0


def check(cond, what):
    global bad, count
    count += 1
    if not cond:
        bad += 1
    print(f"  {'ok  ' if cond else 'FAIL'} {what}")


def printed(s, what=None):
    check(re.sub(r"\s+", " ", s) in TEXN, f"printed: {s}" + (f"  ({what})" if what else ""))


def find(key, pat, n=None, flags=0):
    """The groups of the first match of pat in source key (all matches when n == 'all')."""
    if n == "all":
        return re.findall(pat, SRC[key], flags)
    m = re.search(pat, SRC[key], flags)
    if m is None:
        check(False, f"pattern not found in {key}: {pat}")
        return None
    return m.groups()


def sci(x, d):
    """x as the paper prints it in scientific form with d decimals in the mantissa: 2.55\\times10^{-6}."""
    m, e = f"{x:.{d}e}".split("e")
    return f"{m}\\times10^{{{int(e)}}}"


def sci_s(s):
    """A printed e-number string (e.g. -2.55e-06) in the paper's form, digits kept."""
    m, e = s.lower().split("e")
    return f"{m}\\times10^{{{int(e)}}}"


SQ2 = math.sqrt(2)

# ---------------------------------------------------------------------------------------------------------------------
print("W and the e-folding time (#195's crate, threshold_output.txt)")
W2 = {}
for D in (4, 5, 6):
    g = find("w2", rf"^\s+{D}\s+TS13 Table I\s+[0-9.]+\s+([0-9.]+)", flags=re.M)
    W2[D] = float(g[0])
    printed(str(Decimal(g[0]).quantize(Decimal("1e-17"))), f"W^2 a^2 at D = {D}")
Wi = {}
for D in (4, 5, 6):
    s = find("threshold", rf"SL D={D} S h = [0-9.]+: .*?1/W = ([0-9.]+)")[0]
    Wi[D] = 1 / math.sqrt(W2[D])                                     # the exact e-fold; the outputs print it to 6 digits
    check(s == f"{Wi[D]:.6f}", f"1/W = {s} at D = {D} is 1/sqrt(W^2a^2)")
printed(f"$1/(Wa)={Wi[4]:.6f}$, ${Wi[5]:.6f}$ and ${Wi[6]:.6f}$")
printed(f"buys ${math.log(10) * Wi[4]:.2f}a$", "ln 10 / W at D = 4")
printed(f"every ${Wi[4]:.3f}a$")
printed(f"$0.60\\sqrt2\\,a={0.60 * SQ2:.2f}a$ against $1/W={Wi[4]:.3f}a$")

print("Parameters of the runs (threshold.py CRIT)")
for pat, txt in ((r"G_B=0\.1\b", "$B=0.1$"), (r"pulse=dict\(xc=3\.0, w=0\.5\)", "g(s)=\\exp(-(s-3a)^2/(0.5a)^2)"),
                 (r"rcut=0\.05, fcut=25\.0", "$r<0.05a$ or $e^{-f}>25$"), (r"trap_tol=1e-2, rexp=2\.0", "$r\\ge2a$"),
                 (r"bisect_tol=1e-14", "at most $10^{-14}$"), (r"windows=\[\(1e-11, 1e-3\)", "$10^{-11}\\le|p-\\pst|\\le10^{-3}$"),
                 (r"S=\{2: 28\.0, 3: 18\.0, 4: 15\.0\}.*G=\{2: 28\.0, 3: 22\.0, 4: 28\.0\}", "$X=28$, $18$ and $15$ for S and O in $D=4$, $5$, $6$, and $X=28$, $22$, $28$ for G")):
    check(re.search(pat, SRC["thr_py"]) is not None, f"threshold.py has {pat}")
    printed(txt)

# ---------------------------------------------------------------------------------------------------------------------
print("Table 1: p*_h (threshold_output.txt, PS lines)")
ps = {}
for D, fam, a, b, c, order, rich in find("threshold", r"PS D=(\d) (\w): p\*_h = ([-0-9.e]+), ([-0-9.e]+), ([-0-9.e]+); "
                                                       r"Richardson order ([0-9.]+);.*?Richardson p\* = ([-0-9.e]+)$", "all", re.M):
    ps[(int(D), fam)] = (float(a), float(b), float(c), float(rich))
    r = float(rich)
    rs = f"${sci(r, 1)}$" if abs(r) < 1e-5 else f"${sci(r, 4)}$"
    printed(f"{D} & {fam} & ${sci(float(a), 4)}$ & ${sci(float(b), 4)}$ & ${sci(float(c), 4)}$ & {order} & {rs} \\\\")
g6 = find("threshold", r"PS D=6 G: two h only: p\*_h = ([-0-9.e]+), ([-0-9.e]+), ratio ([0-9.]+)")
printed(f"6 & G & ${sci(float(g6[0]), 4)}$ & ${sci(float(g6[1]), 4)}$ & -- & -- & -- \\\\")
printed(f"whose ratio is ${g6[2]}$")
check(len(ps) == 8, f"eight three-h rows ({len(ps)})")
rmax = max(abs(ps[(D, f)][3]) for D in (4, 5, 6) for f in "SO")
printed(f"zero to within ${sci(rmax, 1)}$", "largest |Richardson p*| of S and O")
printed(f"$p_*={sci(ps[(4, 'G')][3], 4)}$ in $D=4$ and ${sci(ps[(5, 'G')][3], 4)}$ in $D=5$")
oso = max(abs(ps[(D, "O")][i] / (2 * ps[(D, "S")][i]) - 1) for D in (4, 5, 6) for i in range(3))
check(oso < 1e-3, f"p*_h(O) = 2 p*_h(S) to {oso:.1e} < 1e-3")
printed("twice S's to a relative $10^{-3}$ at every $h$")

# ---------------------------------------------------------------------------------------------------------------------
print("Table 2 and Sec. 3.2: lifetime slopes (threshold_output.txt SL, threshold_check_output.txt K2 and K4, throat_locus_output.txt SL)")
SL = {}
for D, fam, h, s, se, sp0, sp1, rel in find("threshold", r"SL D=(\d) (\w) h = ([0-9.]+): p\*_h = [-0-9.e]+; primary slope ([0-9.]+) "
                                          r"\+- ([0-9.]+) .*?spread \[([0-9.]+), ([0-9.]+)\]; 1/W = [0-9.]+, rel dev ([-+0-9.e]+)", "all"):
    SL[(int(D), fam, float(h))] = (s, se, sp0, sp1, rel)
LOC = {}
for D, h, s, se, sp0, sp1, rel in find("locus", r"SL D=(\d) G h=([0-9.]+): p\*_h = [-0-9.e]+; locus tau slope ([0-9.]+) \+- ([0-9.]+) "
                                     r"\(OLS SE\).*?spread \[([0-9.]+), ([0-9.]+)\]; 1/W = [0-9.]+, rel ([-+0-9.e]+)", "all"):
    LOC[(int(D), float(h))] = (s, se, sp0, sp1, rel)
k4 = {int(D): (a, b, r) for D, a, b, r in find("check", r"D=(\d): ind40 ([0-9.]+), #40 refit over \[1e-9, 1e-4\] ([0-9.]+) \(22 pts\), "
                                                     r"rel ([-+0-9.e]+)", "all")}
k2 = {int(D): r2 for D, r2 in find("check", r"D=(\d): [0-9.]+, [0-9.]+ against 1/W [0-9.]+ \(rel [-+0-9.e]+, ([-+0-9.e]+);", "all")}
FINE = {4: 0.005, 5: 0.00176777, 6: 0.00125}
GFINE = {4: 0.005, 5: 0.00176777, 6: 0.0025}
for D in (4, 5, 6):
    s, se, a, b, _ = SL[(D, "S", FINE[D])]
    printed(f"{D} & {Wi[D]:.6f} & S & ${s}\\pm{se}$ & $[{a}, {b}]$ & {k4[D][1]} & {k4[D][0]} \\\\")
    s, se, a, b, _ = SL[(D, "O", FINE[D])]
    printed(f" & & O & ${s}\\pm{se}$ & $[{a}, {b}]$ & & \\\\")
    s, se, a, b, _ = SL[(D, "G", GFINE[D])]
    printed(f" & & G, line & ${s}\\pm{se}$ & $[{a}, {b}]$ & & \\\\")
    if D < 6:
        s, se, a, b, _ = LOC[(D, FINE[D])]
        printed(f" & & G, locus $\\tau$ & ${s}\\pm{se}$ & $[{a}, {b}]$ & & \\\\")
    for fam in "SO":     # 1/W inside the spread at the two finest h
        hs = sorted(h for (d, f, h) in SL if d == D and f == fam)[:2]
        check(all(float(SL[(D, fam, h)][2]) <= Wi[D] <= float(SL[(D, fam, h)][3]) for h in hs), f"1/W inside the spread, D = {D} {fam}, h = {hs}")
printed(", ".join(f"${sci_s(SL[(D, 'S', FINE[D])][4].lstrip('+-'))}$" for D in (4, 5)) +
        f" and ${sci_s(SL[(6, 'S', FINE[6])][4].lstrip('+-'))}$ of $1/W$", "S, primary window, finest h")
printed(f"by ${sci_s(k2[4])}$, ${sci_s(k2[5])}$ and ${sci_s(k2[6])}$ (relative)", "second code, finest h")
hrel = {(int(D), float(h)): float(r) for D, h, r in find("check", r"D=(\d) h=([0-9.]+): p\*_h [0-9.e-]+ \(43 bisection runs\); slope [0-9.]+ "
                                                                  r"\(22 pts; rel ([-+0-9.e]+)\)", "all")}
fac = []
for D in (4, 5, 6):
    hs = sorted((h for d, h in hrel if d == D), reverse=True)
    fac += [hrel[(D, hs[i])] / hrel[(D, hs[i + 1])] for i in range(2)]
check(all(f > 1 for f in fac), f"the second code's deviation falls with h at every halving ({', '.join(f'{f:.2f}' for f in fac)})")
printed(f"by a factor of ${min(fac):.1f}$ to ${max(fac):.1f}$ per halving")
refit = [float(k4[D][1]) / Wi[D] - 1 for D in (4, 5, 6)]
printed(f"which differ from $1/W$ by ${sci(refit[0], 1)}$, ${sci(refit[1], 1)}$ and ${sci(refit[2], 1)}$")
printed(f"agree to ${sci_s(k4[4][2].lstrip('+-'))}$, ${sci_s(k4[5][2].lstrip('+-'))}$ and ${sci_s(k4[6][2].lstrip('+-'))}$")
tau = [abs(float(a) - float(b)) for D, h in ((4, "0\\.005"), (5, "0\\.00176777"), (6, "0\\.00125"))
       for a, b in find("check", rf"D={D} h={h}: .*?slope ([0-9.]+) \(22 pts.*?in tau ([0-9.]+)", "all")]
check(len(tau) == 3 and max(tau) < 1e-5, f"tau slope equals t slope to {max(tau):.1e} < 1e-5 at the finest h")
printed("equals its slope in $t$ to $10^{-5}a$")
low = (float(SL[(6, "S", 0.00125)][0]) / float(k4[6][1]) - 1) * 100
printed(f"slope lies ${abs(low):.1f}\\%$ below the refit", "D = 6 window effect")
worst = {D: max(abs(float(k4[D][1]) / Wi[D] - 1), abs(float(k4[D][0]) / Wi[D] - 1)) for D in (4, 5, 6)}
bnd = []
for D in (4, 5, 6):
    e = math.floor(math.log10(worst[D]))
    bnd.append(f"${math.ceil(worst[D] / 10 ** e)}\\times10^{{{e}}}$")
printed(f"to within a relative {bnd[0]}, {bnd[1]} and {bnd[2]} in $D=4$, $5$ and $6$", "abstract: both codes, rounded up")
ra = {(int(D), f): r for D, f, r in find("threshold", r"RA\[D=(\d) (\w)\]: .*?rel diff ([0-9.e-]+)\)", "all")}
printed(f"changes by ${sci_s(ra[(5, 'S')])}$ in $D=5$ and ${sci_s(ra[(6, 'S')])}$ in $D=6$")
check(float(ra[(4, "S")]) < 1e-9, f"D = 4 two finest S slopes coincide (rel diff {ra[(4, 'S')]})")
gp = [round(100 * float(SL[(D, "G", GFINE[D])][4])) for D in (4, 5, 6)]
printed(f"by ${gp[0]}\\%$, ${gp[1]}\\%$ and ${gp[2]}\\%$ in $D=4$, $5$, $6$", "G on the line")
# [SH02] exponents and [ST17] Table I (alpha_GB = 0 rows), transcribed; the comparison is Paper III's
printed("exponents $0.60$ (ghost pulses) and $0.61$ (normal pulses)")
pred = {D: Wi[D] / SQ2 for D in (4, 5, 6)}
printed(f"$1/(\\sqrt2W)={pred[4]:.3f}$")
printed(f"against ${pred[4]:.3f}$, ${pred[5]:.3f}$ and ${pred[6]:.3f}$")
ST17 = {4: [(0.03, 2.94), (0.10, 2.09)], 5: [(0.15, 2.26), (0.61, 1.66)], 6: [(0.38, 1.92), (1.50, 1.46)]}
st = {}
for D, ((e1, x1), (e2, x2)) in ST17.items():
    lr = math.log(e2 / e1)
    sl = (x1 - x2) / lr
    err = 0.01 / lr + sl * (0.005 / e1 + 0.005 / e2) / lr
    st[D] = (f"{sl:.2f}", f"{err:.2f}") if D < 6 else (f"{sl:.3f}", f"{err:.3f}")
printed(f"${st[4][0]}\\pm{st[4][1]}$, ${st[5][0]}\\pm{st[5][1]}$ and ${st[6][0]}\\pm{st[6][1]}$")

# ---------------------------------------------------------------------------------------------------------------------
print("Sec. 3.3: family G on the locus (throat_locus_output.txt)")
printed(", ".join(f"${LOC[(4, h)][0]}$" for h in (0.02, 0.01)) + f" and ${LOC[(4, 0.005)][0]}$ at the three $h$")
r4 = find("locus", r"RA\[D=4 G\]: .*?; rel ([-+0-9.e]+)\)")[0]
printed(f"agreeing to ${sci_s(r4.lstrip('+'))}$")
printed(f"gives ${find('locus', r'SL D=4 G h=0[.]005:.*?[\n]\s+locus t slope ([0-9.]+)', flags=re.S)[0]}$")
pv = find("locus", r"drift v = 0\.1 .*?locus tau slope ([0-9.]+) .*?rel to 1/W ([-+0-9.e]+),.*?diagonal t \(#40's events\) ([0-9.]+) .*?rel ([-+0-9.e]+),")
printed(f"reads ${pv[0]}$ on the locus (${sci_s(pv[1])}$ against $1/W$) and ${pv[2]}$ on the line (${100 * float(pv[3]):+.2f}\\%$)")
cm = find("locus", r"CMP1\[D=4 G h=0\.005\].*?\(t = 12\.0: 1/lambda ([0-9.]+), member ([0-9.]+) \(r_th ([0-9.]+), gamma_1 ([0-9.]+)\)")
printed(f"$r_{{\\rm th}}={cm[2]}a$ and asymmetry $\\gamma_1={cm[3]}$")
printed(f"$T(\\gamma_1)\\,r_{{\\rm th}}$~\\cite{{GGS2009a}} is ${cm[1]}a$")
printed(f"gives ${cm[0]}a$")
printed(f"the e-folding time ${float(cm[1]):.3f}a$ of that member", "abstract")
blk = SRC["locus"].split("MB D=4 G h=0.005:")[1].split("[pass] CMP1[D=4 G h=0.005]")[0]
rth = {float(t): r for t, r in re.findall(r"t\s+([0-9.]+) tau\s+[0-9.]+ r_th ([0-9.]+)", blk)}
check(rth.get(13.0) == rth.get(14.0) == "1.0028", f"r_th at t = 13, 14: {rth.get(13.0)}, {rth.get(14.0)}")
printed(f"$r_{{\\rm th}}={rth[13.0]}a$ at $t=13a$ and $14a$")
qmax = max(abs(float(q)) for t, q in re.findall(r"t\s+([0-9.]+) tau.*?q ratio ([-+0-9.e]+)", blk) if 11 <= float(t) <= 15)
check(qmax < 1e-5, f"|q_n^2/q_g^2| <= {qmax:.1e} < 1e-5 for t = 11 .. 15")
printed("its share of the static charge is below $10^{-5}$")
dr = find("locus", r"DR\[D=4 G h=0\.005\].*?\(r_th ([0-9.]+) -> ([0-9.]+), change ([-+0-9.]+)\)")
printed(f"changes by ${dr[2]}a$ from $t=10a$ to $15a$")
s4 = float(LOC[(4, 0.005)][0])
both = max(abs(s4 / float(cm[1]) - 1), abs(s4 / float(cm[0]) - 1))
printed(f"The series slope ${s4:.4f}a$ is within ${100 * both:.2f}\\%$ of both")
printed(f"its slope in $D=4$ is ${s4:.4f}a$", "abstract")
d5 = find("locus", r"DR\[D=5 G h=0\.00176777\].*?\(r_th ([0-9.]+) -> ([0-9.]+), change")
printed(f"grows from ${d5[0]}a$ to ${d5[1]}a$ between $t=7a$ and $13a$")
c3 = find("locus", r"CMP3\[D=5 G h=0\.00176777\].*?member - 1 in \[([-+0-9.]+), ([-+0-9.]+)\]")
printed(f"to within ${100 * max(abs(float(c3[0])), abs(float(c3[1]))):.2f}\\%$ along the way")
printed(f"The series slope, ${LOC[(5, 0.00176777)][0]}a$ at the finest $h$")

# ---------------------------------------------------------------------------------------------------------------------
print("Sec. 3.4: the settled mass of S (settled_mass_output.txt, threshold_output.txt formation table)")
form = {(int(D), f): (float(r), float(m)) for D, f, h, r, m in re.findall(
    r"^(\d) ([SO])\s+([0-9.]+)\s+-1 collapse\s+[0-9.]+\s+\S+\s+\S+\s+([0-9.]+)\s+([0-9.]+)", SRC["threshold"], re.M)
    if float(h) == FINE[int(D)]}
fs = [f"{form[(D, 'S')][1]:.3f}" for D in (4, 5, 6)]
check(all(f"{form[(D, 'O')][1]:.3f}" == fs[i] for i, D in enumerate((4, 5, 6))), "O's formation mass equals S's to three digits")
check(all(round(form[(D, f)][0], 2) == 0.99 for D in (4, 5, 6) for f in "SO"), "r = 0.99 a at the first trapped cell, S and O")
printed(f"${fs[0]}$, ${fs[1]}$ and ${fs[2]}a^{{D-3}}$ in $D=4$, $5$, $6$, at $r\\approx0.99a$")
tbl = SRC["settled"].split("delta = 1e-08, h = 0.0025: x+_f")[1].split("==")[0]
rows = {r[0]: r for r in re.findall(r"^\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)", tbl, re.M)}
printed(f"falls from ${rows['0.00'][5]}a$ at the first trapped cell, is ${rows['2.50'][5]}a$ at $x^+-x^+_f=2.5a$")
late = find("settled", r"^\s+1e-08\s+0\.0025\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.e-]+)\s+[0-9.]+\s+[-+0-9.e]+\s+[-+0-9.e]+\s+[0-9.e-]+\s+([0-9.]+)", flags=re.M)
printed(f"settles at ${late[1]}a$ ($r_\\AH={late[3]}a$) at $h=0.0025$")
pd = find("settled", r"delta = 1e-08, h = \(0\.01, 0\.005, 0\.0025\): .*?plateau\(1\.0\) \[[^]]*\], differences [-+0-9.e]+, \+?([0-9.e-]+)")
printed(f"The two finest $h$ differ by ${sci_s(pd[0])}$")
printed(f"ranges over ${sci_s(late[2])}$")
rc7 = find("settled", r"RC7: .*?\(\|dM\| [0-9.e-]+, [0-9.e-]+, [0-9.e-]+, ([0-9.e-]+) \(bar")
printed(f"differ by ${sci_s(rc7[0])}$")
mm = find("rerun41", r"at 100% of the run: .*?\|dM\| ([0-9.e-]+)")
printed(f"at $h=0.01$ to ${sci_s(mm[0])}$", "machine to machine (I41-mac-review)")
printed(f"$M_\\AH={float(late[1]):.4f}a$ ($r_\\AH={float(late[3]):.4f}a$)")
printed(f"settles at ${float(late[1]):.4f}a$ (apparent-horizon radius ${float(late[3]):.4f}a$", "abstract")
scan = re.findall(r"^\s+(1e-0[3-8]|0\.0001|0\.001)\s+\S+\s+\S+\s+\S+\s+\S+\s+([0-9.]+)\s+[0-9.]+\s+[0-9.e+]+\s+limit", SRC["settled"], re.M)
check(len(scan) == 6 and all(m == f"{float(late[1]):.4f}" for _, m in scan), f"scan delta = 1e-8 .. 1e-3 all {[m for _, m in scan]}")
printed(f"It is the same, ${float(late[1]):.4f}a$, for every $\\delta$ from $10^{{-8}}$ to $10^{{-3}}$")
F = find("settled", r"F reached \[[0-9, ]+\] of ([0-9.e+]+)\)")
check(F[0] == "1.57286e+06" and 1.5 * 2 ** 20 == 1572864, f"limit {F[0]} = 1.5 * 2^20")
printed("reach $e^{-f}=1572864$ on the horizon")
printed("ending at $e^{-f}=1.57\\times10^6$")
check("9.40" in rows, "the last row of the M_AH table is at x+ - x+_f = 9.40")
printed("$9.4a$ of $x^+$ after formation")
p1 = find("settled", r"P1 h = 0\.0025: .*?max \|M - 0\.3\| ([0-9.e-]+)")
printed(f"mass error ${sci_s(p1[0])}$ at $h=0.0025$")

print("Table 3 and Sec. 3.4: family G (settled_mass_g_output.txt)")
g0 = find("settled_g", r"delta = 1e-08, h = 0\.005: G: x\+_f.*?\n.*?\n\s+0\.00\s+0\.0000\s+[0-9.]+\s+([0-9.]+)\s+([0-9.]+)\s+\+([0-9.e-]+)", flags=re.S)
printed(f"at the first trapped cell, $\\delta=10^{{-8}}$ & {g0[1]} & {g0[0]} \\\\")
LT = {}
for blkname, part in (("G", SRC["settled_g"].split("== Late-time values")[1].split("delta = 1e-06, h = (0.02")[0]),):
    rows_ = re.findall(r"^\s+(1e-0[68])\s+([0-9.]+)\s+[0-9.]+\s+([0-9.]+)", part, re.M)
    for d, h, p in rows_:
        LT.setdefault((d, h), []).append(p)
gS = {k: v for k, v in LT.items()}
check(len(gS.get(("1e-08", "0.005"), [])) == 2, "late-time table has G then S at h = 0.005, delta = 1e-8")
Gp8, Sp8 = gS[("1e-08", "0.005")]
Gp6, Sp6 = gS[("1e-06", "0.005")]
printed(f"settled, $\\delta=10^{{-8}}$ & {Sp8} & {Gp8} \\\\")
gr = find("settled_g", r"G: plateau mean over the last 1\.0 of x\+ at h = 0\.005, delta = 1e-8: ([0-9.]+) \(r_AH = ([0-9.]+)\)")
check(gr[0] == Gp8, f"G's plateau {gr[0]} equals the late-time table's {Gp8}")
printed(f"apparent-horizon radius $r_\\AH$ & {2 * float(Sp8):.4f} & {float(gr[1]):.4f} \\\\")
dG = find("settled_g", r"delta = 1e-08, h = \(0\.02, 0\.01, 0\.005\): G M_end .*?plateau\(1\.0\) \[[^]]*\], differences [-+0-9.e]+, \+?([0-9.e-]+), Richardson ([0-9.]+)\n"
                       r"\s+S plateau\(1\.0\) \[[^]]*\], differences [-+0-9.e]+, \+?([0-9.e-]+), Richardson ([0-9.]+); G - S at the same h: ([-+0-9.e]+), ([-+0-9.e]+), ([-+0-9.e]+)")
printed(f"$|M(h=0.01)-M(h=0.005)|$ & ${sci_s(dG[2])}$ & ${sci_s(dG[0])}$ \\\\")
printed(f"Richardson, $h=0.02$, $0.01$, $0.005$ & {dG[3]} & {dG[1]} \\\\")
printed(f"settled, $\\delta=10^{{-6}}$ & {Sp6} & {Gp6} \\\\")
check(dG[4] == dG[5] == dG[6], f"G - S the same at the three h ({dG[4]}, {dG[5]}, {dG[6]})")
printed(f"at ${Gp8}a$ ($r_\\AH={float(gr[1]):.4f}a$) for $\\delta=10^{{-8}}$, ${sci_s(dG[6].lstrip('+'))}a$ above S")
printed(f"and at ${float(Gp8):.4f}a$ for the one-sided family", "abstract")
printed(f"${float(dG[6]):.4f}a$ less than for the one-sided", "introduction")
printed(f"Of it, ${sci(float(g0[2]), 1)}a$ is present at formation")
fl = find("settled_g", r"RC4: first law .*?\n((?:\s+h=.*\n){6})")
ok = all(a == b for a, b in re.findall(r"flux ([-0-9.e]+) \(ghost alone ([-0-9.e]+)\)", fl[0]))
check(ok, "the first law with both fields equals the ghost alone (the normal field adds no flux)")
sc = SRC["settled_g"].split("family G, h = 0.005")[1].split("\n\n")[0]
pl = {d: float(p) for d, p in re.findall(r"^\s+(1e-\d\d)\s+\S+\s+\S+\s+\S+\s+([0-9.]+)", sc, re.M) if float(d) <= 1e-3}
check(min(pl, key=pl.get) == "1e-05", f"the plateau's minimum is at delta = {min(pl, key=pl.get)}")
inc = [pl[f"1e-{k + 1:02d}"] - pl[f"1e-{k:02d}"] for k in range(6, 12)]
check(all(b > a for a, b in zip(inc, inc[1:])), "the rise per decade grows down to 1e-12")
printed(f"rises by ${sci(min(inc), 0)}$ to ${sci(max(inc), 0)}$ per decade as $\\delta$ falls from $10^{{-6}}$ to $10^{{-12}}$".replace(".\\times", "\\times"))
printed(f"reaching ${pl['1e-12']:.6f}a$ at $\\delta=10^{{-12}}$")
rr = find("rerun253", r"delta 1e-06: plateau .*?diff -?([0-9.e-]+);")
printed(f"gave the plateau to ${sci_s(rr[0])}$", "I253-mac-review")

print("Sec. 3.4: [SH02]'s data through the chain (settled_mass_output.txt, h = 0.005)")
sh = SRC["settled"].split("Table I data through the chain, h = 0.005")[1].split("################ sh02_uniform")[0]
shr = re.findall(r"^\s+-[0-9.]+\s+\d+\s+([0-9.e-]+)\s+([0-9.e-]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+[0-9.]+\s+\d+\s+limit\s+[0-9.]+\s+([0-9.]+)", sh, re.M)
check(len(shr) == 6, f"six [SH02] rows ({len(shr)})")
e0 = max(abs(float(a) / float(b) - 1) for a, b, *_ in shr)
xh = max(abs(float(r[2]) / float(r[3]) - 1) for r in shr)
printed(f"injected energy to within ${100 * e0:.1f}\\%$ and their horizon coordinate to within ${100 * xh:.1f}\\%$")
printed(f"a settled mass of ${shr[0][4]}a$ for their smallest pulse against their ${shr[0][5]}a$")
dxs = [float(r[6]) for r in shr]
printed(f"equals their printed values ${min(dxs):.1f}a$ to ${max(dxs):.1f}a$ of $x^+$ after formation")
for s in ("0.30a", "0.60\\ln(1/E_0)", "rate $1.1/a$", "radius $0.45a$", "tends to $0.85a$", "$0.22a$~\\cite{GGS2009b}",
          "$r_\\AH=0.45a$", "$m_\\AH\\approx0.22$"):
    printed(s, "transcribed from [SH02] / [GGS09b]")

# ---------------------------------------------------------------------------------------------------------------------
print("Table 4 and Sec. 3.5: the expanding branch (expansion_output.txt)")
EX = {}
for D, h, mem, body in re.findall(r"--- D(\d), h = ([0-9.]+), member (\S+) \(A = [^)]*\) ---\n(.*?)(?=\n  --- |\n\[pass\] MB1|\Z)", SRC["expansion"], re.S):
    k, ts = re.search(r"1/H linear in tau: k = ([0-9.]+), tau_s = ([0-9.]+)", body).groups()
    hw = re.search(r"slope test \(#211\): H/W at r/a = 10: ([0-9.]+), 25: ([0-9.]+), 50: ([0-9.]+)", body).groups()
    rho50 = re.search(r"^\s+50\.0\d\s+(?:\S+\s+){6}(\S+)", body, re.M).group(1)
    dts = re.search(r"last live cell tau = [0-9.]+, tau_s - tau = ([0-9.]+)", body).group(1)
    ex = re.search(r"d ln\(\.\)/d ln\(tau_s - tau\): r ([-0-9.]+), H ([-0-9.]+), \|rho\| ([-0-9.]+), \|R\| ([-0-9.]+), K ([-0-9.]+)", body).groups()
    tl = re.search(r"last live diagonal cell: t = ([0-9.]+)", body).group(1)
    EX[(int(D), float(h), mem)] = dict(k=k, ts=ts, hw=hw, rho50=rho50, dts=dts, ex=ex, t=tl, body=body)
MEM = {"tuned": "$\\pst+10^{-6}$", "0.01": "$p=0.01$", "0.05": "$p=0.05$", "0.2": "$p=0.2$"}
for D, h, mem in ((4, 0.01, "tuned"), (4, 0.01, "0.01"), (4, 0.01, "0.05"), (4, 0.01, "0.2"), (5, 0.00354, "tuned"), (5, 0.00354, "0.05")):
    e = EX[(D, h, mem)]
    printed(f"{D} & {MEM[mem]} & {e['k']} & {e['ts']} & {e['rho50']} & {', '.join(e['hw'])} \\\\")
t4, t5 = EX[(4, 0.01, "tuned")], EX[(5, 0.00354, "tuned")]
printed(f"k={t4['k']}\\ (D=4)")
printed(f"k={t5['k']}\\ (D=5)")
printed(f"with $k={t4['k']}$ in $D=4$ and ${float(t5['k']):.2f}$ in $D=5$", "abstract")
printed(f"are ${t4['ex'][0]}$, ${t4['ex'][1]}$, ${t4['ex'][2]}$ and ${t4['ex'][4]}$ in $D=4$", "exponents of r, H, |rho|, K")
printed(f"$\\tau_s-\\tau={t4['dts']}a$ in $D=4$ and ${t5['dts']}a$ in $D=5$")
printed(f"$H/W$ is ${t4['hw'][0]}$, ${t4['hw'][1]}$ and ${t4['hw'][2]}$ at $r/a=10$, $25$, $50$ in $D=4$")
kd = max(abs(float(EX[(4, 0.01, m)]["k"]) - 2) for m in ("tuned", "0.01", "0.05", "0.2"))
printed(f"have the same $k$ to ${sci(kd, 0)}$".replace(".\\times", "\\times"))
res = {D: max(float(x) for line in re.findall(rf"RES D{D} \S+: .*?differs by \[([^]]+)\]", SRC["expansion"]) for x in line.split(","))
       for D in (4, 5)}
printed(f"to ${sci(res[4], 1)}$ ($D=4$) and ${sci(res[5], 1)}$ ($D=5$)")
printed(f"to $t={round(float(t4['t']))}a$ in $D=4$ and $t={round(float(t5['t']))}a$ in $D=5$")
cc = {n: n * (1 / k + (n + 1) / 2) for n, k in ((2, 2), (3, 4))}
check(cc[2] == 4 and cc[3] == 27 / 4, f"c = n(1/k + (n+1)/2): {cc}")
cm4, cm5 = 2 * (1 / float(t4["k"]) + 3 / 2), 3 * (1 / float(t5["k"]) + 2)
check(abs(cm4 - float(t4["rho50"])) < 1e-3 and abs(cm5 - float(t5["rho50"])) < 1e-3, f"c(k) {cm4:.4f}, {cm5:.4f} against |rho|/H^2 at r = 50a")
printed(f"$c={cm4:.4f}$ in $D=4$ and ${cm5:.4f}$ in $D=5$, against the measured $|\\rho|/H^2={t4['rho50']}$ and ${t5['rho50']}$")
wb = find("expansion", r"D4, h = 0\.01000, member tuned.*?w_bar in \[([-0-9.]+), ([-0-9.]+)\]", flags=re.S)
check(wb == ("-0.3333", "-0.3333"), f"w_bar = -1/3 on the D = 4 throat ({wb})")
sfit = [float(x) for x in re.findall(r"r/a <= (?:1\.5|2|3|5|10): ([0-9.]+) a\^-1", t4["body"])]
printed(f"gives ${min(sfit):.2f}$ to ${max(sfit):.2f}/a$ for $r_{{\\max}}$ from $1.5$ to $10$")
kx = {(int(D), m): k for D, m, k in find("locus", r"KX D(\d) (\S+) h=[0-9.]+: class \w+, k ([0-9.]+) against", "all")}
k4g = [kx[(4, m)] for m in ("tuned", "0.05", "0.2")]
printed(f"$k={min(k4g)}$ to ${max(k4g)}$ for three members in $D=4$ and $k={kx[(5, 'tuned')]}$ for the tuned member in $D=5$")
kk = re.findall(r"^\s+0\.3151\s+282\.84\s+.*?\s([0-9.]+)\s+[0-9.]+\s*$", t4["body"], re.M)
check(kk and kk[0] == "4.0000", f"K (tau_s - tau)^2 at the last cell, D = 4: {kk}")
printed("$K(\\tau_s-\\tau)^2\\to4$ in $D=4$")

# ---------------------------------------------------------------------------------------------------------------------
print("Sec. 5: the checks of the code (#39: raw I39-dh610 run01_checks.log, derive_output.txt)")
check("11 of 11 checks pass" in SRC["checks39"], "checks.py run of record: 11 of 11")
printed("\\texttt{checks.py}, 11 of 11 pass")
check(re.search(r"25 of 25|25 checks.*0 fail|passed: 25", SRC["derive"]) is not None, "derive.py: 25 checks pass")
printed("\\texttt{derive.py}, 25 checks")
check(re.search(r"\[pass\] K: .*max diff 0\.00e\+00", SRC["checks39"]) is not None, "numba and numpy kernels bit-identical")
o4 = find("checks39", r"P1 D=4 .*?orders \[([0-9.]+), ([0-9.]+)\]")
o5 = find("checks39", r"P1 D=5 .*?orders \[([0-9.]+), ([0-9.]+)\]")
printed(f"orders ${o4[0]}$, ${o4[1]}$ ($D=4$) and ${o5[0]}$, ${o5[1]}$ ($D=5$)")
hz = [float(x) for a, b in find("checks39", r"P2 D=\d horizon x-: .*?orders \[([0-9.]+), ([0-9.]+)\]", "all") for x in (a, b)]
printed(f"horizon position at orders ${min(hz):.3f}$ to ${max(hz):.3f}$")
p3 = find("checks39", r"r_\+ = 0: ([0-9.]+) \(exact 0\.666667\); r_- = 0: ([0-9.]+) \(exact 1\.500000\); cells misclassified .*?\[0, 0, 0\]")
printed(f"slopes ${p3[0]}$ and ${p3[1]}$ against $2/3$ and $3/2$")
g1 = find("checks39", r"G1 D=4 h=0\.02000: .*?W_fit ([0-9.]+) vs W ([0-9.]+)")
g1b = find("checks39", r"G1 D=4 h=0\.01000: .*?W_fit ([0-9.]+) vs W")
printed(f"${g1[0]}$ and ${g1b[0]}$ at $h=0.02$, $0.01$ against $W={g1[1]}$")
check(N_MUT39 == 14, f"#39's mutation logs of record: {N_MUT39}")
printed("Each of 14 mutations")
print("Sec. 5: the threshold controls (threshold_output.txt) and the second code")
pl40 = find("threshold", r"PL\[D=4\].*?\((\d+) of (\d+) inside, p\* (\d+) of (\d+)\)")
check(all(re.search(rf"PL\[D={D}\].*?\(40 of 40 inside, p\* 40 of 40\)", SRC["threshold"]) for D in (4, 5, 6)), "planted series: 40 of 40 at each D")
check(all(re.search(rf"PX\[D={D}\].*?\(40 of 40\)", SRC["threshold"]) for D in (4, 5, 6)), "sqrt2 time unit rejected 40 of 40 at each D")
printed("for 40 of 40 noise seeds in each $D$")
printed("is rejected in 40 of 40")
g6c = find("threshold", r"G6 D=6: .*?W_fit ([0-9.]+), ([0-9.]+) against W = ([0-9.]+)")
printed(f"${g6c[0]}$ and ${g6c[1]}$ against $W={g6c[2]}$")
p40 = {(int(D), float(h)): float(p) for D, h, p in find("threshold", r"SL D=(\d) S h = ([0-9.]+): p\*_h = ([0-9.e-]+);", "all")}
pchk = [abs(float(p) / p40[(int(D), round(float(h), 8))] - 1) for D, h, p in
        find("check", r"D=(\d) h=([0-9.]+): p\*_h ([0-9.e-]+) \(43 bisection runs\)", "all")]
check(len(pchk) == 9, f"nine p*_h pairs ({len(pchk)})")
printed(f"agrees with the first code's to a relative ${sci(max(pchk), 1)}$")
tu = find("ind40_tu", r"D=4 h=0\.005: .*?slope ([0-9.]+) \(22 pts")
check(abs(float(tu[0]) * SQ2 / Wi[4] - 1) < 1e-5, f"time-unit mutation slope {tu[0]} within 1e-5 of 1/(sqrt2 W) = {Wi[4] / SQ2:.6f}")
printed(f"gives ${tu[0]}$ against $1/(\\sqrt2W)={Wi[4] / SQ2:.6f}$")

print("Sec. 5: the review checks (negative-mass-raw review folders)")
check("[PASS] K1" in SRC["ind39_1"] and "symbolic residuals ['0', '0', '0', '0']" in SRC["ind39_1"], "ind39 K1: sympy residual 0")
k2b = find("ind39_2", r"\[PASS\] K2b: W_fit ([0-9.]+), ([0-9.]+) against ([0-9.]+)")
printed(f"${k2b[0]}$ and ${k2b[1]}$ against ${k2b[2]}$")
i41 = find("ind41", r"\[PASS\] K1b: .*?e\^-f 4\.074: M ([0-9.]+) ")
d41 = [abs(float(x)) for x in re.findall(r"\(([-+][0-9.e-]+)\)", SRC["ind41"].split("[PASS] K1b:")[1].split("\n")[0])]
printed(f"gives $M_\\AH={i41[0]}a$ at $e^{{-f}}=4.074$")
printed(f"agrees with the chain to ${sci(max(d41), 1)}$ up to $e^{{-f}}=11.96$")
i213 = [x for x in find("ind213", r"\[PASS\] K1: k ([0-9., ]+);")[0].split(", ")]
printed(f"$k={min(i213)}$ to ${max(i213)}$ for the eight $D=4$ runs")
i225 = find("ind225", r"\[PASS\] K1: .*?e-fold T r_th = ([0-9.]+);.*?the measured ([0-9.]+)")
printed(f"gives the e-fold ${i225[0]}a$ at the stall member, against the measured ${i225[1]}a$")
print("Sec. 5: the constraint plateau (constraint_plateau_output.txt, ind226, I226-rpi5 supplement)")
r1 = find("plateau", r"static, D = 6, region.*?ratios per halving: ([0-9.]+), ([0-9.]+), ([0-9.]+) \(outside", flags=re.S)
printed(f"ratios per halving are ${r1[0]}$, ${r1[1]}$ and ${r1[2]}$")
c1 = [float(x) for x in re.findall(r"(?:D=[56]) ([0-9.]+), ([0-9.]+)", SRC["plateau"].split("[pass] C1:")[1].split("\n")[0]) for x in x]
printed(f"fall by ${min(c1):.1f}$ to ${max(c1):.1f}$ at every halving")
c2 = find("plateau", r"\[pass\] C2: .*?D=5: a=256 at the finest declared h x([0-9.]+) .*?D=6: a=256 at the finest declared h x([0-9.]+)")
printed(f"by a factor ${round(float(c2[1]))}$ ($D=6$) and ${round(float(c2[0]))}$ ($D=5$)")
k3 = find("ind226", r"\[PASS\] K3: ratios ([0-9.]+), ([0-9.]+), ([0-9.]+);")
printed(f"own march (ind226) gives ${k3[0]}$, ${k3[1]}$ and ${k3[2]}$")
rr = find("r_round", r"max \|r_float64 - r_binary128\| = ([0-9.e-]+)\n.*?max \|r_float64 - r_exact\(l\)\| = ([0-9.e-]+)", flags=re.S)
printed(f"is ${sci(float(rr[0]), 1)}$ against a truncation error of ${sci(float(rr[1]), 1)}$ ($D=6$, $h=0.00125$)")


# #255: the D = 5 drift and G's expanding member (notes/wormhole_evolution/throat_drift_output.txt, h index 1)
dr = SRC["drift"]
m = re.search(r"D=5 h index 1: p\*_h [^;]*; exponent of p\*_h ([0-9.]+) \(0\.04 to 0\.1\)", dr)
rt = re.search(r"D=5 h index 1:.*?\n\s+r_th - 1 at t = [^:]*: .*?exponents \[([0-9.]+), [0-9.]+, ([0-9.]+)\] \(0\.04 to 0\.1\)", dr, re.S)
rd = re.search(r"D=5 h index 1:.*?\n.*?\n\s+rdot_inf \+[0-9.]+ \+- [0-9.e-]+, \+([0-9.]+) \+- [0-9.e-]+, .*?exponents ([0-9.]+),", dr, re.S)
kx = re.search(r"KX2: [^\n]*?\(h = [0-9.]+: FT, k ([0-9.]+), ratio [0-9.]+; h = [0-9.]+: FT, k ([0-9.]+),", dr)
check(m and rt and rd and kx, "#255's exponents, rdot_inf and k in throat_drift_output.txt")
printed(f"$B^{{{float(m[1]):.2f}}}$", "p*_h exponent (h index 1)")
printed(f"$B^{{{float(rt[1]):.1f}}}$ to\n$B^{{{float(rt[2]):.1f}}}$", "r_th - 1 exponents")
printed(f"$\\dot r_\\infty=0.0{round(float(rd[1]) * 1000):d}a$", "rdot_inf at B = 0.1")
printed(f"$B^{{{float(rd[2]):.2f}}}$", "rdot_inf exponent")
check(kx[1] == kx[2], "k equal at both h")
printed(f"$k={kx[1]}$ at two $h$", "G's expanding member, k")


# #230: the late state (notes/wormhole_evolution/throat_scaling_output.txt)
sc = SRC["scaling"]
k6 = [float(x) for x in re.findall(r"D6_A0\.\d+_h0\.003_X200 .*?k_fit=([0-9.]+)", sc)]
c6 = re.search(r"C6: .*?D6_A0\.05_h0\.003_X200 ([0-9.]+), D6_A0\.2_h0\.003_X200 ([0-9.]+), D6_A0\.01_h0\.003_X200 ([0-9.]+)\)", sc)
t8 = re.search(r"D4_tuned_h0\.01_X800 .*?tau_s-tau_end=([0-9.]+)\s+k_fit=([0-9.]+)\s+tau_s=([0-9.]+)", sc)
k5 = re.search(r"K5X: .*?D5_tuned_h0\.00353553 ([0-9.]+),", sc)
kp = [float(x) for x in re.findall(r"^\s+[456]\s+[-+][0-9.]+\s+[0-9.]+\s+([0-9.]+)\s+[0-9.]+\s+[0-9.]+", sc, re.M)]
jt = re.search(r"JTK: .*?\(worst ([0-9.e-]+) over 11 runs", sc)
check(len(k6) == 3 and c6 and t8 and k5 and len(kp) == 11 and jt, "#230's rows in throat_scaling_output.txt")
printed(f"$k_{{\\rm fit}}={min(k6):.3f}$ to ${max(k6):.3f}$", "D = 6 k_fit range")
cv = sorted(float(x) for x in c6.groups())
printed(f"$|\\rho|/H^2={cv[0]:.3f}$ to ${cv[-1]:.3f}$ at $r=50a$", "D = 6 c range")
printed(f"$\\tau_s-\\tau={float(t8[1]):.3f}a$ gives $\\tau_s={float(t8[3]):.5f}a$ and $k={float(t8[2]):.5f}$", "D = 4 at X = 800")
printed(f"extrapolated over the same runs is ${k5[1]}$", "D = 5 extrapolated exponent")
printed(f"for $k$ from ${min(kp):.2f}$ to ${max(kp):.1f}$, to $2\\times10^{{-4}}$", "jump-test k range")
check(f"{float(jt[1]):.0e}" == "2e-04", f"jump-test agreement {jt[1]} rounds to 2e-4")


# #266: G's settled mass against B (notes/wormhole_evolution/settled_mass_b_output.txt, the M(B) table at h = 0.005)
mb = SRC["massb"]
rowsB = {m[0]: (float(m[1]), float(m[2]), float(m[3]), float(m[4])) for m in
         re.findall(r"^\s+(0\.0|0\.025|0\.05|0\.1|0\.125)\s+[0-9.e+-]+\s+([-+0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+[0-9.]+\s+[0-9.]+\s+([0-9.]+)\s+[-+0-9.]+", mb, re.M)}
check(len(rowsB) == 5, f"#266's M(B) table: {sorted(rowsB)}")
printed(f"${rowsB['0.0'][2]:.6f}a$ at $B=0$", "M at B = 0")
printed(f"${rowsB['0.025'][2]:.6f}a$, ${rowsB['0.05'][2]:.6f}a$ and ${rowsB['0.1'][2]:.6f}a$ at $B=0.025$, $0.05$ and $0.1$", "M(B)")
printed(f"$\\gamma_1={rowsB['0.025'][0]:.5f}$, ${rowsB['0.05'][0]:.5f}$ and ${rowsB['0.1'][0]:.5f}$", "gamma_1(B)")
printed(f"$r_{{\\rm th}}={rowsB['0.025'][1]:.5f}a$, ${rowsB['0.05'][1]:.5f}a$ and ${rowsB['0.1'][1]:.5f}a$", "r_th(B)")
import numpy as _np
g = _np.array([rowsB[k][0] for k in ("0.025", "0.05", "0.1", "0.125")])
dm = _np.array([rowsB[k][3] - rowsB["0.0"][3] for k in ("0.025", "0.05", "0.1", "0.125")])     # M/r_th - M_0 from the 6-decimal column
ab = _np.linalg.lstsq(_np.vstack([g, g ** 2]).T, dm, rcond=None)[0]
printed(f"$a={ab[0]:.3f}$ and $b={ab[1]:.2f}$", "the four-point fit of M/r_th - M_0 (to the digits the printed table supports)")
lo = re.search(r"\[FAIL\] REL: .*?\(max error ([0-9.e-]+) \(bar ([0-9.e-]+)\)", mb)
lop = [abs(float(x)) for x in re.findall(r"REL-ph: leave out B = [0-9.]+: .*?error ([-+0-9.e]+)$", mb, re.M)]
check(lo is not None and len(lop) == 4, "#266's leave-one-out rows")
printed(f"(${float(lo[1]) * 1e5:.1f}\\times10^{{-5}}$ against ${float(lo[2]) * 1e5:.0f}\\times10^{{-5}}$)", "the declared leave-one-out")
printed(f"holds to ${max(lop) * 1e6:.1f}\\times10^{{-6}}$", "the post-hoc four-B leave-one-out")
st = re.search(r"first B with no stall \(gain < 2\.0\): ([0-9.]+); the last with a stall: ([0-9.]+)", mb)
check(st is not None and float(st[1]) <= 0.165 and float(st[2]) >= 0.155, f"the stall ends between {st[2] if st else '?'} and {st[1] if st else '?'} at h = 0.005")

# #267: G's settled mass along the stall, delta = 1e-7 ... 1e-12 (notes/wormhole_evolution/settled_mass_member_output.txt, h = 0.005)
mo = SRC["member"]
blk = mo.split("h = 0.005: #266's relation", 1)[1]
rowsD = {int(m[0]): (float(m[1]), float(m[2]), float(m[3]), float(m[4]), float(m[5])) for m in
         re.findall(r"^\s+(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([-+0-9.e]+)\s+\|", blk.split("[FAIL] PRED", 1)[0], re.M)}
check(sorted(rowsD) == list(range(3, 13)), f"#267's M(delta) table at h = 0.005: e = {sorted(rowsD)}")
printed(f"${rowsD[8][0]:.6f}a$ at $\\delta=10^{{-8}}$ to ${rowsD[12][0]:.6f}a$ at $10^{{-12}}$", "M(delta)")
printed(f"misses by ${rowsD[12][4] * 1e4:.1f}\\times10^{{-4}}$ at $10^{{-12}}$", "M - predicted at 1e-12")
check(all(abs(rowsD[e][4]) > 6e-5 for e in (9, 10, 11, 12)), "PRED fails at every delta from 1e-9 (bar 6e-5)")
dec = re.search(r"^\s+12\s+([-+0-9.e]+)\s+([-+0-9.e]+)\s+[-+0-9.e]+\s+[-+0-9.e]+\s+[-+0-9.e]+$", mo.split("Decomposition at h = 0.005", 1)[1], re.M)
check(dec is not None, "#267's decomposition row e = 12")
printed(f"predicts $+{float(dec[2]) * 1e3:.2f}\\times10^{{-3}}$ of the $+{float(dec[1]) * 1e3:.2f}\\times10^{{-3}}$", "dM predicted / measured, 1e-8 to 1e-12")
sp = re.search(r"D = 4 late speed rdot_inf = \+([0-9.]+) .*?ratio D = 5 / D = 4: ([0-9.]+) \(window to 20.5\)", mo)
check(sp is not None, "#267's late speed line")
printed(f"tends to ${float(sp[1]) * 1e3:.2f}\\times10^{{-3}}$ per unit proper time, ${sp[2]}$ times slower than in $D=5$", "the D = 4 drift")
check(len(re.findall(re.escape(f"${sp[2]}$ times"), TEXN)) == 3, "the ratio 20.6 is printed three times (intro, body, limits)")
t14 = re.search(r"t = 14: tau\s+[0-9.]+, r_th ([0-9.]+),", mo)
check(t14 is not None and float(t14[1]) < min(float(x) for x in re.findall(r"t = (?:12|16): tau\s+[0-9.]+, r_th ([0-9.]+),", mo)), "r_th at t = 14 is below t = 12 and 16")
printed(f"falls to ${float(t14[1]):.4f}a$ at $t=14$", "the radius's minimum region")
gl = re.search(r"\(gamma_1 ([0-9.]+) to ([0-9.]+)\)", mo)
printed(f"$\\gamma_1$ falls from ${float(gl[1]):.4f}$ to ${float(gl[2]):.4f}$", "gamma_1 along the stall")
tdl = re.search(r"^\s+h = 0\.005: e = 3: .*$", mo, re.M)[0]
tdv = {int(e): float(t) for e, t in re.findall(r"e = (\d+): T [0-9./]+, t_dep ([0-9.]+)", tdl)}
td = [round(tdv[e + 1] - tdv[e], 3) for e in range(7, 12)]
check(all(f"{x:.1f}" == "2.0" for x in td), f"t_dep gains 2.0 per decade at h = 0.005 (delta = 1e-7 ... 1e-12): {td}")
printed("moves the departure by $2.0$ in $t$", "t_dep per decade")

# #272: the linear spectrum beyond the unstable mode (notes/wormhole_evolution/ellis_spectrum_output.txt)
so = SRC["spectrum"]
s2 = re.findall(r"^PASS    S2-D[45]-(?:even|odd): no zero with Re W in \[-1\.4, -0\.2\]", so, re.M)
check(len(s2) == 4, f"#272: S2 passes in four of four sectors ({len(s2)})")
zs = [(float(a), float(b)) for a, b in re.findall(r"zero at W = ([-0-9.]+) [+0-9.]+i; over ray angle and far point \(n = \d+\): spread ([0-9.e+-]+)", so)]
conf = [-a for a, sp_ in zs if sp_ <= 1e-3]
check(len(conf) == 4, f"#272: four confirmed zeros (spread <= 1e-3): {conf}")
printed(f"its nearest zeros decay at rates ${min(conf):.2f}$", "the slowest confirmed zero")
t7 = {D: [float(x) for x in re.search(rf"T7-D{D}: .*?rate x t at t = 15, 20, 25: \[([^\]]+)\]", so)[1].split()] for D in (4, 5)}
printed(f"the speed falls only as $t^{{-{sum(t7[5]) / 3:.2f}}}$", "D = 5 tail exponent, mean of rate x t at t = 15, 20, 25")
printed(f"a speed falling as $t^{{-{sum(t7[4]) / 3:.1f}}}$", "D = 4 tail exponent, mean of rate x t")

print(f"\n{count} checks, {bad} failed; mutation {OPT.mutate}")
sys.exit(1 if bad else 0)
