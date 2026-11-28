#!/usr/bin/env python
"""#274: the winding2 boxes of #126 (31) and #117 (4) recounted with the crate's complex lane (f64) as a second evaluator.

    python recount_i274.py --raw=RAW --stars                       the star of each (omega, eps) by the crate and by the scipy path -> RAW/stars.json
    python recount_i274.py --raw=RAW --list                        the box list (job, source, region, stored result)
    python recount_i274.py --raw=RAW [--boxes=a,b] [--lane=f64]    count the boxes (default: all 35; --extra adds #115's s1 at l = 5), one crate process per box -> RAW/results/record/<box>.json
    python recount_i274.py --raw=RAW --mutate=NAME[,NAME]         a control on every box (or --boxes): exit 1, every targeted check must fail
    (add --sub=dev to a counting or control run: its results and evaluations go to results_dev/, jsonl_dev/)
    python recount_i274.py --raw=RAW --collect > OUT               the tables of the records and of the controls; exit 1 on a box not OK, missing or not run

RAW is the negative-mass-raw checkout (its folders I126-rpi5 and I117-rpi5 hold jobs.txt, the logs and the caches of winding2.py's
runs at --guard=all; the folder I274-mac of the run is RAW/I274-mac, given by --out-dir, default RAW/I274-mac).
Each box is counted exactly as winding2.py counts it: its eps, omega, l, region, --y0, --rsmult, --zfac, --plant and --guard=all, the
scalar plane started at winding2's fixed R_s = r_m + rsmult / (min Re k on the contour), its phase strip applied in double precision, and the
counter of taylor_polar.py (a copy of winding2's count() with the margins), with F(s) from the crate's `dserve`.

Per box (each check can fail; --mutate shows each failing):
  INPUT   the R_s and r_m computed here equal the ones in winding2's log (1 and 2 decimals).
  STORED  the counter on the stored cache alone (no crate) reproduces the log's winding, max steps and min|D| (equal as printed), the planted
          check, and asks no point the cache lacks.
  COUNT   the counter with the crate: winding within 0.05 of the stored one, the planted difference within 0.05 of the stored expectation,
          every phase step under 0.45 rad.
  MARGIN  the smallest distance of any compared quantity from its threshold (0.4, 0.4, 0.1 rad), without the plant and with it, is at least
          10 E, E = max(3.1e-6 rad, the largest |arg(crate / stored)| over the points both have). 3.1e-6 is the largest phase error of
          the stored caches measured in #188 (l = 6). A box under this is reported.
Recorded, not checked: the points the counter asks against the stored cache's, min|D|, the walls.

Controls (--mutate=): reverse, wrongdouble, misplant, onesided (the counter's own, as taylor_polar.py --complex): COUNT and STORED must fail
(onesided is the identity in an upper box); jsign, star_eps, rs_shift, strip_off (the crate's evaluation, eight stored points per box):
MARGIN must fail; kmin_shift (R_s computed 1 percent off): INPUT must fail.
"""
import json
import os
import pickle
import re
import statistics
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np  # noqa: E402
import taylor_polar as TP  # noqa: E402

CRIT = dict(WIND=0.05, STEP=0.45, PH_FLOOR=3.1e-6, FACTOR=10.0, TIME_LIMIT=3600.0, NCTL=8, RS_TOL=0.051, RM_TOL=0.0051)
COUNTER_MUTS = ("reverse", "wrongdouble", "misplant", "onesided")
EVAL_MUTS = ("jsign", "star_eps", "rs_shift", "strip_off")
ALL_MUTS = COUNTER_MUTS + EVAL_MUTS + ("kmin_shift",)
FOLDERS = (("126", "I126-rpi5"), ("117", "I117-rpi5"))
# beyond the two lists (--extra): #115's s1 at l = 5, counted there at --guard=all (its guard-1 count lost 2 pi); the log is #115's
EXTRA = (("s1_l5", "115", "I115-rpi5", "s1_l5_guardall", "-1 0.7603 5 0.03 0.1 2.5 --plant=0.04,0.25 --guard=all"),)

RESULT_RE = re.compile(r"RESULT (\S+): zeros in Re s∈\[(\S+),(\S+)\], (?:\|Im s\|≤(\S+)|Im s∈\[(\S+),(\S+)\]): winding = ([+-]\d+\.\d+)  "
                       r"\(max phase step ([\d.]+) rad, min\|D\| (\S+)\)")
PLANT_RE = re.compile(r"PLANTED (\S+): s0 = \((\S+)\), expected (\d) planted zeros inside; count with the factor ([+-]\d+\.\d+), without ([+-]\d+\.\d+), "
                      r"difference ([+-]\d+\.\d+)  \(max phase step ([\d.]+) rad\)  (PASS|FAIL)")
HEAD_RE = re.compile(r"(\S+): r_m=([\d.]+), R_s fixed = ([\d.]+) \(min Re k on contour ([\d.]+)\)")


class BoxTimeout(Exception):
    pass


def ckey(s):
    return TP.ckey(s)


# ---------------------------------------------------------------------------------------------------------------
# the boxes: winding2's own job lists and logs
# ---------------------------------------------------------------------------------------------------------------
def parse_args(argstr):
    tok = argstr.split()
    pos = [t for t in tok if not t.startswith("--")]
    opt = dict(t[2:].split("=", 1) for t in tok if t.startswith("--") and "=" in t)
    eps, w, l = int(pos[0]), float(pos[1]), int(pos[2])
    x0, x1, y1 = (float(a) for a in pos[3:6])
    y0 = float(opt.get("y0", 0.0))
    rsmult = float(opt.get("rsmult", 40.0))
    zfac = float(opt.get("zfac", 1.0))
    plant = tuple(float(a) for a in opt["plant"].split(","))
    assert opt.get("guard") == "all", argstr
    sym = (y0 == 0.0)
    tag = (f"eps{eps:+d}_w{w}_l{l}_{x0}_{x1}_{y1}" + ("" if sym else f"_y0{y0}") + ("" if rsmult == 40.0 else f"_rs{rsmult:g}")
           + ("" if zfac == 1.0 else f"_z{zfac:g}"))
    return dict(eps=eps, w=w, l=l, x0=x0, x1=x1, y1=y1, y0=y0, rsmult=rsmult, zfac=zfac, plant=plant, sym=sym, tag=tag)


def load_boxes(raw):
    boxes = []
    for issue, folder in FOLDERS:
        for line in open(os.path.join(raw, folder, "jobs.txt")):
            f = line.rstrip("\n").split("|")
            if len(f) < 3 or not f[0].startswith("w_") or f[0].startswith("w_mut_"):
                continue
            bx = parse_args(f[2])
            bx.update(name=f[0][2:], issue=issue, folder=folder, job=f[0], argstr=f[2])
            bx["log"] = os.path.join(raw, folder, f[0] + ".log")
            bx["cache_file"] = os.path.join(raw, folder, "cache", f"winding2_cache_{bx['tag']}.pkl")
            bx["extra"] = False
            boxes.append(bx)
    for name, issue, folder, job, argstr in EXTRA:
        bx = parse_args(argstr)
        bx.update(name=name, issue=issue, folder=folder, job=job, argstr=argstr, extra=True)
        bx["log"] = os.path.join(raw, folder, job + ".log")
        bx["cache_file"] = os.path.join(raw, folder, "cache", f"winding2_cache_{bx['tag']}.pkl")
        boxes.append(bx)
    return boxes


def stored_log(bx):
    """winding2's own printout of the run: r_m, R_s, RESULT and PLANTED lines (exactly one of each)"""
    text = open(bx["log"], encoding="utf-8").read()
    head, res, pl = HEAD_RE.findall(text), RESULT_RE.findall(text), PLANT_RE.findall(text)
    assert len(head) == 1 and len(res) == 1 and len(pl) == 1, (bx["name"], len(head), len(res), len(pl))
    h, r, p = head[0], res[0], pl[0]
    assert p[7] == "PASS", (bx["name"], "the stored planted check did not pass")
    assert h[0] == bx["tag"] and r[0] == bx["tag"] and p[0] == bx["tag"], (h[0], bx["tag"])
    return dict(rm=float(h[1]), rs=float(h[2]), kmin=float(h[3]), wind=float(r[6]), step=r[7], nmin=r[8], expected=int(p[2]), plant_with=float(p[3]),
                plant_without=float(p[4]), diff=float(p[5]), plant_step=p[6], plant_pass=(p[7] == "PASS"))


def kmin_of(bx):
    """min Re k over the contour as winding2.py computes it"""
    w = bx["w"]
    x0, x1, y0, y1 = bx["x0"], bx["x1"], bx["y0"], bx["y1"]
    ks = lambda s: (np.sqrt(1 + (s - 1j * w) ** 2), np.sqrt(1 + (s + 1j * w) ** 2))
    cands = [x0 + 1j * y for y in np.linspace(y0, y1, 200)] + [x + 1j * y1 for x in np.linspace(x0, x1, 200)]
    if not bx["sym"]:
        cands += [x1 + 1j * y for y in np.linspace(y0, y1, 200)] + [x + 1j * y0 for x in np.linspace(x0, x1, 200)]
    k = float(min(min(c.real for c in ks(c0)) for c0 in cands))
    return k * 1.01 if MUT == "kmin_shift" else k


def strip(s, rs, rm, w):
    """winding2.py's removal of the trivial phases exp(-i Im[(k+ + k-)(R_s - r_m) + (k+ + k- + s) r_m])"""
    if MUT == "strip_off":
        return complex(1.0)
    kp, km = np.sqrt(1 + (s - 1j * w) ** 2), np.sqrt(1 + (s + 1j * w) ** 2)
    return complex(np.exp(-1j * ((kp + km) * (rs - rm) + (kp + km + s) * rm).imag))


# ---------------------------------------------------------------------------------------------------------------
# stars
# ---------------------------------------------------------------------------------------------------------------
def star_key(bx):
    return f"{bx['eps']:+d}_{bx['w']}"


def make_stars(raw, outdir, boxes):
    import subprocess
    import background as B
    from polar_solver import PolarProblem
    out = {}
    for bx in boxes:
        k = star_key(bx)
        if k in out:
            continue
        t = time.time()
        st = B.solve_star(bx["w"], bx["eps"])
        P = PolarProblem(st, bx["l"])
        t_sc = time.time() - t
        t = time.time()
        p = subprocess.run(["nice", "-n", "10", TP.CBIN, "star", f"omega={bx['w']!r}", f"eps={float(bx['eps'])!r}", f"lam={TP.LAM!r}", f"mu={TP.MU!r}"],
                           capture_output=True, text=True, check=True)
        o = json.loads(p.stdout.splitlines()[0])
        out[k] = dict(w=bx["w"], eps=bx["eps"], pc=o["pc"], zc=o["zc"], rm=float(P.rm), crate_phic=o["phic"], crate_alpha_c=o["alpha_c"], scipy_phic=st.phic,
                      scipy_alpha_c=st.alpha_c, scipy_M=st.M, scipy_R=st.R, scipy_r0=st.r0, t_scipy=t_sc, t_crate=time.time() - t)
        print(f"star {k}: r_m {out[k]['rm']!r}, crate phi_c {o['phic']:.12e} alpha_c {o['alpha_c']:.12e}, scipy phi_c(r0) {st.phic:.12e} alpha_c(r0) {st.alpha_c:.12e}", flush=True)
    json.dump(out, open(os.path.join(outdir, "stars.json"), "w"), indent=1)


# ---------------------------------------------------------------------------------------------------------------
# one box
# ---------------------------------------------------------------------------------------------------------------
MUT = None


def counter_mut():
    return MUT if MUT in COUNTER_MUTS else None


class Lazy:
    """the crate's `dserve` for one box, started at the first point the jsonl file does not hold (a control that only re-counts stored
    evaluations needs no process); a new evaluation is appended to the file"""

    def __init__(self, make, jsonl):
        self.make, self.jsonl, self.srv, self.n, self.mem = make, jsonl, None, 0, {}
        if jsonl and os.path.exists(jsonl):
            for line in open(jsonl):
                o = json.loads(line)
                self.mem[(o["re"], o["im"])] = o

    @property
    def setup(self):
        return self.srv.setup if self.srv else {}

    def eval(self, s):
        k = (s.real, s.imag)
        if k not in self.mem:
            if self.srv is None:
                self.srv = self.make()
            self.mem[k] = self.srv.eval(s)
            self.n += 1
        return self.mem[k]

    def close(self):
        if self.srv:
            self.srv.close()


def crate_F(bx, star, lane, jsonl, rs, t0):
    rs_int = rs + (1.0 if MUT == "rs_shift" else 0.0)
    TP.MUTATE = MUT if MUT in ("jsign", "star_eps") else None     # the crate's own mutations; rs_shift and strip_off are in this file
    mk = lambda: TP.Server("dserve", jsonl, lane=lane, l=bx["l"], rsfix=repr(rs_int), zfac=repr(bx["zfac"]), rm=repr(star["rm"]), omega=repr(bx["w"]),
                           eps=repr(float(bx["eps"])), lam=repr(TP.LAM), mu=repr(TP.MU), pc=star["pc"], zc=star["zc"])
    srv = Lazy(mk, jsonl)
    walls = {}

    def F(s):
        if time.time() - t0 > CRIT["TIME_LIMIT"]:
            raise BoxTimeout()
        s = ckey(s)
        o = srv.eval(s)
        walls[s] = o["wall"]
        return complex(o["d_f64"][0], o["d_f64"][1]) * strip(s, rs, star["rm"], bx["w"])
    F.walls = walls
    return F, srv


def pick_points(bx, cache, n):
    """n stored points spread over the cache (not on the real axis), for the controls that mutate the crate"""
    ks = sorted((k for k in cache if abs(k.imag) > 1e-9), key=lambda c: (c.real, c.imag))
    step = max(1, len(ks) // n)
    return ks[::step][:n]


def phase_diff(a, b):
    return abs(float(np.angle(a / b)))


def run_box(bx, stars, raw, lane, jdir):
    t_box = time.time()
    star = stars[star_key(bx)]
    log = stored_log(bx)
    cache = pickle.load(open(bx["cache_file"], "rb"))
    rm = star["rm"]
    kmin = kmin_of(bx)
    rs = rm + bx["rsmult"] / kmin
    checks, notes = {}, []
    res = dict(box=bx["name"], issue=bx["issue"], job=bx["job"], l=bx["l"], eps=bx["eps"], w=bx["w"], x0=bx["x0"], x1=bx["x1"], y0=bx["y0"], y1=bx["y1"],
               plant=list(bx["plant"]), rsmult=bx["rsmult"], zfac=bx["zfac"], sym=bx["sym"], rm=rm, rs=rs, kmin=kmin, log=log, n_cache=len(cache), mutate=MUT,
               lane=lane)
    checks["INPUT"] = abs(rs - log["rs"]) <= CRIT["RS_TOL"] and abs(rm - log["rm"]) <= CRIT["RM_TOL"]
    cm = counter_mut()
    args = dict(x0=bx["x0"], x1=bx["x1"], y1=bx["y1"], y0=bx["y0"], guard="all", plant=bx["plant"], mutate=cm)

    # STORED: the counter on the stored cache alone
    asked_s, missing_s = set(), set()
    rst = TP.count_winding(TP._stored_F(cache, asked_s, missing_s), **args)
    res["stored_run"] = dict(rst, asked=len(asked_s), missing=len(missing_s))
    checks["STORED"] = (abs(rst["winding"] - log["wind"]) < CRIT["WIND"] and rst["plant_ok"] and not missing_s and f"{rst['maxstep']:.2f}" == log["step"]
                        and f"{rst['plant_maxstep']:.2f}" == log["plant_step"] and f"{rst['nmin']:.2e}" == log["nmin"]
                        and abs(rst["plant_diff"] - log["diff"]) < CRIT["WIND"])
    if MUT in EVAL_MUTS or MUT == "kmin_shift":
        # the crate is evaluated at stored points only; the margin is the stored run's
        res["points"] = None
        if MUT == "kmin_shift":
            res["checks"] = checks
            res["fired"] = not checks["INPUT"]
            res["effect"] = f"R_s {rs:.1f} against the log's {log['rs']:.1f}"
            return finish(res, checks, t_box)
        keys = pick_points(bx, cache, CRIT["NCTL"])
        F, srv = crate_F(bx, star, lane, None, rs, t_box)
        try:
            ph, md = [], []
            for k in keys:
                d = F(k)
                ph.append(phase_diff(d, cache[k]))
                md.append(abs(abs(d) / abs(cache[k]) - 1))
        finally:
            srv.close()
        E = max(CRIT["PH_FLOOR"], max(ph))
        margin = min(rst["margin"], rst["plant_margin"])
        checks["MARGIN"] = margin >= CRIT["FACTOR"] * E
        res.update(agree=dict(n=len(keys), ph_max=max(ph), mod_max=max(md)), E=E, margin_used=margin, wall_crate=sum(F.walls.values()), checks=checks)
        res["fired"] = not checks["MARGIN"]
        res["effect"] = f"{len(keys)} stored points: phase {max(ph):.2e} rad (max), modulus {max(md):.2e}; margin {margin:.2e} against 10 E = {CRIT['FACTOR'] * E:.2e}"
        return finish(res, checks, t_box)

    # the crate as the evaluator (the counter's own mutations reuse the unmutated evaluations)
    jsonl = os.path.join(jdir, f"{bx['name']}_{lane}.jsonl")
    F, srv = crate_F(bx, star, lane, jsonl, rs, t_box)
    asked = set()
    Fa = lambda s: (asked.add(ckey(s)), F(s))[1]
    stopped = None
    try:
        try:
            rc = TP.count_winding(Fa, **args)
        except BoxTimeout:
            stopped, rc = "stopped at the 1 h wall limit", None
    finally:
        srv.close()
    if stopped:
        res.update(stopped=stopped, checks=checks, verdict="NOT RUN")
        return finish(res, checks, t_box, stopped=True)
    n_new = srv.n
    common = [k for k in asked if k in cache]
    ph = [phase_diff(F(k), cache[k]) for k in common]
    md = [abs(abs(F(k)) / abs(cache[k]) - 1) for k in common]
    walls = sorted(F.walls.values())
    E = max(CRIT["PH_FLOOR"], max(ph))
    margin = min(rc["margin"], rc["plant_margin"])
    checks["COUNT"] = (abs(rc["winding"] - log["wind"]) < CRIT["WIND"] and abs(rc["plant_diff"] - log["expected"]) < CRIT["WIND"] and rc["plant_ok"]
                       and rc["maxstep"] < CRIT["STEP"] and rc["plant_maxstep"] < CRIT["STEP"])
    checks["MARGIN"] = margin >= CRIT["FACTOR"] * E
    res.update(crate=dict(rc, asked=len(asked), n_evaluated_now=n_new, wall_sum=sum(walls), wall_median=statistics.median(walls), wall_max=walls[-1],
                          t_star=srv.setup.get("t_star", 0.0), t_table=srv.setup.get("t_table", 0.0)),
               points=dict(asked=len(asked), stored=len(cache), asked_not_stored=len(asked - set(cache)), stored_not_asked=len(set(cache) - asked)),
               agree=dict(n=len(common), ph_max=max(ph), ph_median=statistics.median(ph), mod_max=max(md)), E=E, margin_used=margin, checks=checks)
    if cm:
        res["fired"] = (not checks["COUNT"]) and (not checks["STORED"])
        res["effect"] = (f"plant difference {rc['plant_diff']:+.4f} (crate), {rst['plant_diff']:+.4f} (stored) against {log['expected']}; winding {rc['winding']:+.4f}"
                         if not (cm == "onesided" and bx["sym"] is False) else "identity in an upper box")
        if cm == "onesided" and not bx["sym"]:
            res["fired"] = None
    return finish(res, checks, t_box)


def finish(res, checks, t_box, stopped=False):
    res["wall_total"] = time.time() - t_box
    res["verdict"] = "NOT RUN" if stopped else ("OK" if all(checks.values()) else "REPORT " + ",".join(k for k, v in checks.items() if not v))
    return res


def box_line(r):
    if "crate" not in r:
        return f"{r['box']}: {r.get('effect', r.get('stopped', ''))}  [{r['verdict']}]"
    c, a, p = r["crate"], r["agree"], r["points"]
    return (f"{r['box']} (l={r['l']}, {r['x0']}..{r['x1']}, Im {r['y0']}..{r['y1']}, R_s {r['rs']:.1f}): winding {c['winding']:+.4f}, max step {c['maxstep']:.2f}, "
            f"min|D| {c['nmin']:.2e}, plant {c['plant_diff']:+.4f} of {r['log']['expected']} (step {c['plant_maxstep']:.2f}); margin {c['margin']:.2e}, with plant "
            f"{c['plant_margin']:.2e}; points asked {p['asked']}, stored {p['stored']}, asked not stored {p['asked_not_stored']}, stored not asked {p['stored_not_asked']}; "
            f"crate against stored over {a['n']} points: phase {a['ph_max']:.2e} rad, modulus {a['mod_max']:.2e}; {c['n_evaluated_now']} evaluated now, "
            f"crate wall {c['wall_sum']:.1f} s; checks {' '.join(k + '=' + ('ok' if v else 'FAIL') for k, v in r['checks'].items())}  [{r['verdict']}]")


# ---------------------------------------------------------------------------------------------------------------
# the collector
# ---------------------------------------------------------------------------------------------------------------
def collect(raw, outdir, boxes):
    rec = {}
    d = os.path.join(outdir, "results", "record")
    for bx in boxes:
        f = os.path.join(d, bx["name"] + ".json")
        if os.path.exists(f) and not (MUT == "missing" and bx["name"] == boxes[0]["name"]):
            rec[bx["name"]] = json.load(open(f))
    print("# #274: the winding2 boxes of #126 (31) and #117 (4) counted with the crate's complex lane (f64) as the evaluator, --guard=all")
    lists = [b for b in boxes if not b["extra"]]
    print(f"{sum(b['name'] in rec for b in lists)} of {len(lists)} boxes of the two lists have a record"
          + (f"; {sum(b['name'] in rec for b in boxes if b['extra'])} of {sum(b['extra'] for b in boxes)} extra box outside them (#115's s1 at l = 5, last rows)" if len(lists) < len(boxes) else ""))
    print("| box | source | l | eps, omega | region | R_s | stored winding | crate winding | planted (stored) | planted (crate) | max step (stored, crate) | min abs D (stored, crate) | verdict |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    bad, bad_ctl = 0, []
    for bx in boxes:
        r = rec.get(bx["name"])
        if r is None:
            print(f"| {bx['name']} | #{bx['issue']} | {bx['l']} | | | | | | | | | | **not run: no record** |")
            bad += 1
            continue
        lg, c = r["log"], r.get("crate")
        reg = f"Re [{r['x0']}, {r['x1']}], " + (f"abs(Im) ≤ {r['y1']}" if r["sym"] else f"Im [{r['y0']}, {r['y1']}]")
        if c is None:
            print(f"| {r['box']} | #{r['issue']} | {r['l']} | {r['eps']:+d}, {r['w']} | {reg} | {r['rs']:.1f} | {lg['wind']:+.4f} | | | | | | **{r['verdict']}** {r.get('stopped', '')} |")
            bad += 1
            continue
        bad += r["verdict"] != "OK"
        print(f"| {r['box']} | #{r['issue']} | {r['l']} | {r['eps']:+d}, {r['w']} | {reg} | {r['rs']:.1f} | {lg['wind']:+.4f} | {c['winding']:+.4f} | {lg['diff']:+.4f} of {lg['expected']} | "
              f"{c['plant_diff']:+.4f} of {lg['expected']} | {lg['step']}, {c['maxstep']:.2f} | {lg['nmin']}, {c['nmin']:.2e} | {r['verdict']} |")
    print()
    print("| box | margin | margin with plant | E (phase error used) | crate against stored: phase max | modulus max | margin / 10 E | asked | stored | asked not stored | stored not asked | evaluations | crate wall (s) | median / max per D (s) | setup (s) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for bx in boxes:
        r = rec.get(bx["name"])
        if r is None or "crate" not in r:
            continue
        c, a, p = r["crate"], r["agree"], r["points"]
        print(f"| {r['box']} | {c['margin']:.2e} | {c['plant_margin']:.2e} | {r['E']:.2e} | {a['ph_max']:.2e} | {a['mod_max']:.2e} | {r['margin_used'] / (CRIT['FACTOR'] * r['E']):.1f} | "
              f"{p['asked']} | {p['stored']} | {p['asked_not_stored']} | {p['stored_not_asked']} | {c['n_evaluated_now']} | {c['wall_sum']:.1f} | {c['wall_median']:.3f} / {c['wall_max']:.3f} | "
              f"{c['t_star'] + c['t_table']:.1f} |")
    print()
    print("| box | STORED (counter on the stored cache) | checks |")
    print("|---|---|---|")
    for bx in boxes:
        r = rec.get(bx["name"])
        if r is None:
            continue
        s = r["stored_run"]
        print(f"| {r['box']} | winding {s['winding']:+.4f}, max step {s['maxstep']:.2f}, min abs D {s['nmin']:.2e}, plant {s['plant_diff']:+.4f}, asked {s['asked']}, missing {s['missing']}, "
              f"margin {s['margin']:.2e} ({s['plant_margin']:.2e} with plant) | {' '.join(k + '=' + ('ok' if v else 'FAIL') for k, v in r['checks'].items())} |")
    # the controls
    print()
    ctl = {}
    for m in ALL_MUTS:
        dm = os.path.join(outdir, "results", m)
        if os.path.isdir(dm):
            ctl[m] = {n[:-5]: json.load(open(os.path.join(dm, n))) for n in os.listdir(dm) if n.endswith(".json")}
    print("| control | target | boxes fired | not fired | effect over the boxes (smallest … largest) |")
    print("|---|---|---|---|---|")
    tgt = {m: "COUNT, STORED" for m in COUNTER_MUTS}
    tgt.update({m: "MARGIN" for m in EVAL_MUTS})
    tgt["kmin_shift"] = "INPUT"

    def effect(m, rs_):
        if m in COUNTER_MUTS:   # distance of the crate's planted difference from the expectation (the check's tolerance is 0.05)
            v = [abs(r["crate"]["plant_diff"] - r["log"]["expected"]) for r in rs_.values() if "crate" in r and r.get("fired") is not None]
            return f"planted difference off the expectation by {min(v):.4f} … {max(v):.4f} (tolerance {CRIT['WIND']})"
        if m == "kmin_shift":
            v = [abs(r["rs"] - r["log"]["rs"]) for r in rs_.values() if "log" in r]
            return f"R_s off the log's by {min(v):.1f} … {max(v):.1f} (tolerance {CRIT['RS_TOL']})"
        v = [r["agree"]["ph_max"] for r in rs_.values() if "agree" in r]
        q = [r["margin_used"] / (CRIT["FACTOR"] * r["E"]) for r in rs_.values() if "agree" in r]
        stopped = sum(1 for r in rs_.values() if "agree" not in r)
        return (f"phase against the stored values {min(v):.2e} … {max(v):.2e} rad; margin / 10 E {min(q):.2g} … {max(q):.2g} (fires below 1)"
                + (f"; the mutated crate stopped in {stopped}" if stopped else ""))
    for m, rs_ in ctl.items():
        fired = [n for n, r in rs_.items() if r.get("fired") is True]
        na = [n for n, r in rs_.items() if r.get("fired") is None]
        nf = [n for n, r in rs_.items() if r.get("fired") is False] + [b["name"] for b in lists if b["name"] not in rs_]
        bad_ctl += [f"{m} on {n}" for n in nf]
        print(f"| `{m}` | {tgt[m]} | {len(fired)} of {len(lists) - len(na)}{f' ({len(na)} n/a)' if na else ''} | {', '.join(nf) if nf else 'none'} | {effect(m, rs_)} |")
    print()
    print(f"RESULT: {len(boxes) - bad} of {len(boxes)} box rows OK" + ("" if not bad_ctl else f"; controls NOT fired: {', '.join(bad_ctl)}")
          + ("; every control fired on every box it applies to" if not bad_ctl else ""))
    return 1 if (bad or bad_ctl) else 0


# ---------------------------------------------------------------------------------------------------------------
def main():
    global MUT
    args = sys.argv[1:]
    opt = dict(a[2:].split("=", 1) for a in args if a.startswith("--") and "=" in a)
    raw = opt.get("raw", os.path.join(HERE, "..", "..", "..", "negative-mass-raw"))
    outdir = opt.get("out-dir", os.path.join(raw, "I274-mac"))
    lane = opt.get("lane", "f64")
    MUT = opt.get("mutate")
    muts = MUT.split(",") if MUT else [None]
    assert all(m in ALL_MUTS + ("missing", None) for m in muts), muts
    boxes = load_boxes(raw)
    if "boxes" in opt:
        want = opt["boxes"].split(",")
        assert all(w in [b["name"] for b in boxes] for w in want), want
        boxes = [b for b in boxes if b["name"] in want]
    elif "--extra" not in args:
        boxes = [b for b in boxes if not b["extra"]]
    if "--list" in args:
        for b in boxes:
            lg = stored_log(b)
            print(f"{b['name']}|#{b['issue']}|{b['argstr']}|winding {lg['wind']:+.4f}, max step {lg['step']}, min|D| {lg['nmin']}, planted {lg['diff']:+.4f} of {lg['expected']} {'PASS' if lg['plant_pass'] else 'FAIL'}|R_s {lg['rs']}")
        print(f"{len(boxes)} boxes")
        return 0
    if "--stars" in args:
        make_stars(raw, outdir, boxes)
        return 0
    if "--collect" in args:
        return collect(raw, outdir, [b for b in load_boxes(raw) if not b["extra"] or "--extra" in args])
    if not os.path.exists(TP.CBIN):
        sys.exit("build first: python taylor_polar.py --complex --build")
    stars = json.load(open(os.path.join(outdir, "stars.json")))
    sub = opt.get("sub")      # development runs: --sub=dev keeps their results and evaluations apart from the runs of record
    jdir = os.path.join(outdir, "jsonl_" + (sub or "rec"))
    os.makedirs(jdir, exist_ok=True)
    print(f"# recount_i274.py: start {time.strftime('%Y-%m-%d %H:%M:%S')}, load average {os.getloadavg()[0]:.2f}, lane {lane}, boxes {len(boxes)}, mutate = {MUT}", flush=True)
    bad, t0 = 0, time.time()
    for m in muts:
        MUT = m
        rdir = os.path.join(outdir, "results" + (f"_{sub}" if sub else ""), m or "record")
        os.makedirs(rdir, exist_ok=True)
        for bx in boxes:
            try:
                r = run_box(bx, stars, raw, lane, jdir)
            except Exception as e:
                if m in EVAL_MUTS:   # a mutated crate that cannot produce D (its star does not solve) has failed the check, as taylor_polar.py counts a crash
                    r = dict(box=bx["name"], verdict=f"REPORT MARGIN (the mutated crate stopped: {type(e).__name__})", checks={"MARGIN": False}, mutate=m,
                             fired=True, effect=f"the mutated crate stopped ({type(e).__name__}), no D")
                else:
                    r = dict(box=bx["name"], verdict=f"NOT RUN: {type(e).__name__}: {str(e)[:200]}", checks={"RUN": False}, mutate=m)
            json.dump(r, open(os.path.join(rdir, bx["name"] + ".json"), "w"), indent=1)
            print(box_line(r) + (f"  FIRED={r.get('fired')}  {r.get('effect', '')}" if m and "fired" in r else ""), flush=True)
            if m:
                bad += r.get("fired") is False or "RUN" in r.get("checks", {})
            else:
                bad += r["verdict"] != "OK"
    if muts != [None]:
        print(f"CONTROL {','.join(muts)}: " + ("fired on every box it applies to, as it must" if not bad else f"NOT FIRED or not run on {bad} box rows"))
        print(f"SUMMARY: {len(boxes)} boxes; {time.time() - t0:.0f} s")
        return 0 if bad else 1    # a control that works exits 1, as its targeted checks fail
    print(f"SUMMARY: {len(boxes) - bad} of {len(boxes)} boxes OK, {bad} not OK or not run; {time.time() - t0:.0f} s")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
