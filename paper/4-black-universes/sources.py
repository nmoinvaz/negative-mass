"""Parsers of the notes' outputs of record that Paper IV quotes (notes/black_universe/*_output.txt, the review sections of
its notes, and the coordinator's review files in ../negative-mass-raw when that checkout is present).

Imported by paper_numbers.py and make_figures.py. MUTATIONS holds the --mutate= controls of paper_numbers.py: each replaces one
number in one source file as it is read, so a check downstream of it must fail.
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
NB = os.path.join(REPO, "notes", "black_universe")
RAW = os.path.normpath(os.path.join(REPO, "..", "negative-mass-raw"))     # beside the checkout (AGENTS.md)

# name: (file under notes/black_universe, text as stored, text as mutated)
MUTATIONS = {
    "W": ("stability_output.txt", "0.4941711675 [0.471899]", "0.4941711775 [0.471899]"),          # W at c = -1.5
    "c2": ("horizon_sector_output.txt", "L     -9.307319e-13  -9.610862e-04", "L     -9.307319e-13  -9.610872e-04"),
    "slope": ("evolve_output.txt", "primary slope 1.57042", "primary slope 1.57052"),                 # c = -1, h = 0.0025
    "rate": ("linear_output.txt", "(0.636627 against kappa", "(0.636637 against kappa"),            # R3, c = -1
    "slope262": ("cminus1_mode_output.txt", "beta*  = 0.3140290188873", "beta*  = 0.3140290198873"),  # d(W/kappa)/dc, #262
}
_mut = None


def set_mutation(name):
    global _mut
    if name is not None and name not in MUTATIONS:
        raise SystemExit(f"unknown mutation {name}; known: {', '.join(MUTATIONS)}")
    _mut = name


def read(name, base=NB):
    text = open(os.path.join(base, name)).read()
    if _mut and base == NB and MUTATIONS[_mut][0] == name:
        old, new = MUTATIONS[_mut][1:]
        assert text.count(old) == 1, f"mutation {_mut}: '{old}' not unique in {name}"
        text = text.replace(old, new)
    return text


def stability():
    """#212: the member table of 'stability.py report', the planted case, the spectral check and the controls."""
    t = read("stability_output.txt")
    rows = []
    rx = (r"^\s*(-?[\d.]+)\s+(-?[\d.]+)\s+([\d.]+)\s+\|\s+(\d)\s+(\d)\s+\|\s+(none|[\d.]+)(?:\s+\[([\d.]+)\])?\s+(\d+|-)"
          r"\s+([\d.]+)?\s*\|\s+([\d.e+-]+) ok")
    for m in re.finditer(rx, t, re.M):
        rows.append(dict(c=m[1], rho_h=m[2], kappa=m[3], n1=int(m[4]), n2=int(m[5]), W=None if m[6] == "none" else m[6],
                         W_bkz=m[7], nodes=m[8], W_kappa=m[9], r12=m[10]))
    planted = {r: (w, d) for r, w, d in re.findall(r"P1 planted Ellis \(c = 0\) (R\d): growing modes 1 \(declared 1\); W\^2 a\^2 = "
                                                   r"([\d.]+) ; reference [\d.]+; \|difference\| = ([\d.e+-]+)", t)}
    ref = re.search(r"reference ([\d.]+);", t)[1]
    spec = re.findall(r"^PASS  S6 spectral (S-[AB]) \((u|Y)-form\) c = (-?[\d.]+): positive real eigenvalues (\d) \(shooting (\d)\)"
                      r".*?\|dW\|/W = \[(?:'([\d.e+-]+)')?\]", t, re.M)
    spec_fail = len(re.findall(r"^FAIL  S6", t, re.M))
    controls = re.findall(r"^fires   --mutate=(\w+)", t, re.M)
    pt = {k: re.search(rf"{k} [^\n]*?= ([\d.e+-]+) \(declared <= 3e-3\)", t)[1] for k in ("S7b", "S7d")}
    s8 = re.search(r"S8 Vhat .*?max abs difference ([\d.e+-]+)", t)[1]
    s3 = re.search(r"S3 horizon .*?max abs difference ([\d.e+-]+)", t)[1]
    s2 = re.search(r"S2 background .*?max residual ([\d.e+-]+)", t)[1]
    s5 = re.search(r"S5 family .*?max abs difference ([\d.e+-]+)", t)[1]
    one = re.search(r"members with exactly one growing mode: (\d+) of (\d+)", t).groups()
    return dict(rows=rows, planted=planted, ref=ref, spec=spec, spec_fail=spec_fail, controls=controls, pt=pt, s8=s8, s3=s3,
                s2=s2, s5=s5, one=one, text=t)


def linear():
    """#232: checks, Table 1 and the c = -1 tables of linear.py."""
    t = read("linear_output.txt")
    d = dict(text=t)
    d["R1"] = {c: (lin, w, rel) for c, lin, w, rel in re.findall(
        r"c = (-[\d.]+): ([\d.]+) against W = ([\d.]+) \(relative ([\d.e+-]+);", t)}
    m = re.search(r"R3: .*?\(([\d.]+) against kappa = ([\d.]+) \(relative ([\d.e+-]+);", t)
    d["R3"] = m.groups()
    m = re.search(r"L3a: .*?late exponent ([\d.]+), [\d.]+, [\d.]+ against b = ([\d.]+) \(relative ([\d.e+-]+)\)", t)
    d["L3a"] = m.groups()
    m = re.search(r"L3b: .*?\(reduction ([\d.e+-]+), rate of the combination (-?[\d.]+), of one pulse ([\d.]+)", t)
    d["L3b"] = m.groups()
    m = re.search(r"L3c: .*?\(reduction ([\d.e+-]+), rate of the combination ([\d.]+), W = ([\d.]+)\)", t)
    d["L3c"] = m.groups()
    m = re.search(r"L4: .*?\(c = -1: ([\d.]+) \(\d+ nodes\).*?c = -1.5: ([\d.]+) \(\d+ nodes\), expected ([\d.]+)\)", t)
    d["L4"] = m.groups()
    m = re.search(r"L5: .*?linear Y_h: [^(]*\(ratios ([\d.]+), ([\d.]+)\); nonlinear drift: [^(]*\(ratios ([\d.]+), ([\d.]+)\)", t)
    d["L5"] = m.groups()
    d["L1"] = max(float(x) for x in re.findall(r"\|ratio - 1\| <= ([\d.e+-]+) for", t))
    d["L2"] = max(float(x) for x in re.findall(r"c = -[\d.]+ [RLG]: ([\d.e+-]+)[;)]", t.split("[pass] L2:")[1].split("\n")[0]))
    m = re.search(r"P1a: .*?\(W\^2 = ([\d.]+), relative difference ([\d.e+-]+);", t)
    d["P1a"] = m.groups()
    d["P1b"] = re.search(r"P1b: .*?relative difference ([\d.e+-]+)\)", t)[1]
    d["P1c"] = re.search(r"P1c: .*?W\^2 shift ([\d.e+-]+)\)", t)[1]
    m = re.search(r"(\d+) checks, (\d+) pass; failed: \[\]", t)
    d["checks"] = m.groups()
    d["mut"] = re.search(r"(\d+) of (\d+) equal", t).groups()
    tab = {}
    for line in t.split("######## linear.py table")[1].split("########")[0].splitlines():
        p = line.split()
        if len(p) == 10 and re.match(r"-[\d.]+$", p[0]):
            tab[p[0]] = p
    d["table"] = tab
    cls = {}
    for m in re.finditer(r"^\s+(R|L|G|A0|A1R|A1L)\s+.*?\s+(-?[\d.]+e[+-]\d+)\s+(-?[\d.]+)\s+([\d.]+)\s+(-?[\d.]+)\s+([\d.]+)$",
                         t.split("######## linear.py c1")[1], re.M):
        cls[m[1]] = m.groups()[1:]
    d["classes"] = cls
    m = re.search(r"c1, c2 of the three: \(([^)]*)\); \(([^)]*)\); \(([^)]*)\)", t)
    d["c12_three"] = [tuple(x.split(", ")) for x in m.groups()]
    d["bkappa"] = re.search(r"b = sqrt2 kappa = ([\d.]+)", t)[1]
    der = read("linear_derive_output.txt")
    d["derive"] = re.search(r"(\d+) checks, (\d+) pass", der).groups() if re.search(r"(\d+) checks, (\d+) pass", der) else None
    d["derive_text"] = der
    return d


def horizon_sector():
    """#248: derive, check, table of horizon_sector.py."""
    t = read("horizon_sector_output.txt")
    d = dict(text=t)
    cls = {}
    for line in t.split("######## horizon_sector.py table")[1].split("samples")[0].splitlines():
        p = line.split()
        if len(p) == 7 and p[0] in ("L", "G", "A0", "A1R", "A1L", "R"):
            cls[p[0]] = p[1:]
    d["classes"] = cls
    d["g2"] = re.search(r"\(g2 = ([^)]*)\)", t)[1]
    n4 = re.findall(r"slices V = 0, 0.5, 1: ([\d.e+-]+), ([\d.e+-]+), ([\d.e+-]+)", t.split("[pass] N4:")[1].split("\n")[0])
    d["N4"] = max(float(x) for tup in n4 for x in tup)
    n5 = t.split("[pass] N5:")[1].split("\n")[0]
    dev = re.findall(r"deviate by ([\d.e+-]+), ([\d.e+-]+), ([\d.e+-]+), ([\d.e+-]+)", n5)
    d["N5"] = max(float(x) for tup in dev for x in tup)
    d["N5_c2"] = re.findall(r"s0 = ([\d.]+): c2 = (-[\d.e+-]+);", n5)
    n6 = t.split("[pass] N6:")[1].split("\n")[0]
    d["N6"] = re.search(r"c=-1.5: known dr_h coefficient (-[\d.e+-]+); exact data ([\d.e+-]+), marched \(Richardson\) ([\d.e+-]+);"
                        r".*?Schwarzschild: known dr_h coefficient (-[\d.e+-]+); exact data ([\d.e+-]+), marched \(Richardson\) "
                        r"([\d.e+-]+);", n6).groups()
    n7 = t.split("[pass] N7:")[1].split("\n")[0]
    d["N7"] = re.search(r"class R: \|c2\| at h = 0.01 \.\.\. 0.00125 = ([\d.e+-]+), ([\d.e+-]+), ([\d.e+-]+), ([\d.e+-]+) "
                        r"\(ratios ([\d.]+), ([\d.]+), ([\d.]+)\)", n7).groups()
    n11 = t.split("[pass] N11:")[1].split("\n")[0]
    d["N11"] = re.search(r"worst slice-fit deviation ([\d.e+-]+); \|c2\| from ([\d.e+-]+) to ([\d.e+-]+); median ([\d.e+-]+)",
                         n11).groups()
    n10 = t.split("[pass] N10:")[1].split("\n")[0]
    d["N10"] = re.findall(r"ratios ([\d.]+), ([\d.]+);", n10)
    d["derive"] = re.search(r"(\d+) checks, (\d+) pass; failed: \[\]; mutation None; \d+ s\n\n######## horizon_sector.py check", t).groups()
    d["check"] = re.search(r"(\d+) checks, (\d+) pass; failed: \[\]; mutation None; \d+ s\n\n######## horizon_sector.py table", t).groups()
    d["mut"] = re.search(r"(\d+) of (\d+) equal", t).groups()
    d["d6"] = re.search(r"\((\d+) equations, (\d+) determined unknowns, (\d+) forced to zero\)", t).groups()
    return d


def evolve():
    """#217: the analyze section of evolve.py (lifetime slopes, A*_h), E6, E8 orders, the one-sided pulse."""
    t = read("evolve_output.txt")
    d = dict(text=t)
    runs = {}
    for m in re.finditer(r"c = (-[\d.]+), h = ([\d.]+): A\*_h = (-[\d.e+-]+) .*?primary slope ([\d.]+) \+- ([\d.]+) \(n = (\d+)\), "
                         r"1/(W|kappa) = ([\d.]+), relative ([-+][\d.e+-]+); .*?fit spread \[([\d.]+), ([\d.]+)\]", t):
        runs[(float(m[1]), m[2])] = dict(Astar=m[3], slope=m[4], se=m[5], n=m[6], ref=m[8], rel=m[9], lo=m[10], hi=m[11])
    d["runs"] = runs
    d["orders"] = {c: o for c, o in re.findall(r"E8c?\[c = (-[\d.]+)\]: .*?orders ([\d., ]+?); Richardson", t)}
    d["E6"] = re.search(r"-1: (-[\d.]+) \(must be < 0\)", t)[1]
    d["E6_rates"] = re.findall(r"(-[\d.]+): ([\d.]+) \(\+([\d.e+-]+)\)", t.split("[pass] E6:")[1].split("\n")[0])
    d["XC"] = re.findall(r"XC\[c = (-[\d.]+)\]: .*?\(([\d.]+) against ([\d.]+)\)", t)
    d["mass"] = dict(re.findall(r"naming runs, c = (-[\d.]+), .*?ADM mass of the member m = ([\d.]+);", t))
    d["checks"] = re.search(r"(\d+) checks, (\d+) pass, 0 fail", t).groups()
    d["mut_check"] = re.findall(r"mutation (\w+): fired .*?: equal", t)
    return d


def polar():
    """#218: polar l = 1, 2 counts and the planted cases."""
    t = read("polar_output.txt")
    d = dict(text=t)
    blocks = re.split(r"^Polar l = (\d), real W in \[([\d.]+), ([\d.]+)\]", t, flags=re.M)
    counts = {}
    for i in range(1, len(blocks), 4):
        ell, lo, hi, body = blocks[i], blocks[i + 1], blocks[i + 2], blocks[i + 3].split("\n\n")[0]
        rows = re.findall(r"^\s+(-?[\d.]+)\s+\S+\s+\S+ \| (\d)\s+(\d)\s+\| (\d)\s+([\d.e+-]+)", body, re.M)
        counts[ell] = dict(lo=lo, hi=hi, rows=rows)
    d["counts"] = counts
    d["P2"] = re.findall(r"^PASS  P2 Ellis polar l = (\d) .*?omega = ([\d.]+) (-[\d.]+)i \(printed ([\d.]+) (-[\d.]+)i\)", t, re.M)
    d["P3a"] = re.search(r"P3a .*?derived ([\d.]+) (-[\d.]+)i, master .*?\|difference\| = ([\d.e+-]+)", t).groups()
    return d


def bounce():
    """#223: the report of bounce.py (55 members)."""
    t = read("bounce_output.txt")
    d = dict(text=t)
    rep = t.split("## bounce.py report")[1].split("## bounce.py spectral")[0]
    mem = {}
    for m in re.finditer(r"^(\w+_[a-z]+[-\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\S+)\s+(\S+)\s+(\S+) \| (\d)\s+(\d)\s+\| (.*?)\s*\| (\d)\s+(\d)\s+"
                         r"\| (\S+)\s+\| (\d)\s+\| (.*)$", rep, re.M):
        mem[m[1]] = dict(m=m[2], q=m[3], xh=m[4], kappa=m[5], p=m[6], main=m[9], count=int(m[13]),
                         cplx="complex pair" in m[14], low="low " in m[9])
    d["members"] = mem
    d["summary"] = re.search(r"members with at least one growing mode .*?: (\d+) of (\d+)", t).groups()
    d["P2"] = re.findall(r"P2 (R[12]) vs #212 c = (-[\d.]+) \(n \d, W (?:[\d.]+|None)\): dW/W ([\d.e+-]+|-) ok", t)
    return d


def lemma():
    t = read("horizon_lemma_output.txt")
    d = dict(text=t)
    d["members"] = re.findall(r"^  (.*?): R'\(rho_h\) = ([\d.]+), kappa = ([\d.]+): (.*)$", t, re.M)
    d["checks"] = re.search(r"(\d+) checks, (\d+) pass; failed: \[\]", t).groups()
    return d


def notes_md(name):
    return read(name)


def raw(path):
    p = os.path.join(RAW, path)
    return open(p).read() if os.path.exists(p) else None
