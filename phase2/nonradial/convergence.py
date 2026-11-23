"""
Robustness of the matching determinant D(gamma) (real growth rates) for the phantom campaign star:
integrator tolerance (1e-10 vs 1e-12), matching radius (phi_cut 1e-9 -> r_m ~ 24, 1e-12 -> r_m ~ 34),
background R (60 vs 90), and inner start r0 (1e-3 vs 3e-3).  Reports D for each variant (sign and relative
spread), plus the interior / exterior Zerilli log-derivatives q_int, q_ext at r_m whose equality is the
mode condition (q_int - q_ext is a resolution-independent, well-conditioned diagnostic; see kappa_eff.py).
Usage: python convergence.py l [w] [--mutate=nocap]
GAMS (a comma list) replaces the default gammas, NPROC (default 1) runs the (variant, gamma) evaluations in that many
worker processes, and each result is appended to convergence_l{l}_w{w}.partial.csv in PARTIAL_DIR (default this
directory) as it finishes, so a stopped run resumes. The partial file also records, as JSON objects, the task and pid of
each evaluation when it starts, finishes (with its wall and cpu seconds) or is capped.
CAP_S (seconds, default no cap) caps each evaluation's wall time, counted from its start in a worker that has already
built the problems. A capped evaluation is stopped with its worker, printed as not completed, and its row's signs and
spread are over the variants that finished. On resume it is rerun only when CAP_S is unset or larger than its cap.
--mutate=nocap ignores CAP_S, the control of the cap (#114). The script exits 0 when every row is printed.
"""
import sys, os, json, queue, numpy as np, time
import multiprocessing as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import background as B
from polar_solver import PolarProblem

ARGS = [a for a in sys.argv[1:] if not a.startswith('--')]
MUT = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--mutate=')), None)
assert MUT in (None, 'nocap')
l = int(ARGS[0]) if len(ARGS) > 0 else 2
w = float(ARGS[1]) if len(ARGS) > 1 else 0.7603
gams = [1e-3, 0.03, 0.3, 1.0, 3.0, 10.0, 40.0, 137.0, 300.0] + ([1000.0] if l >= 3 else [])
if os.environ.get('GAMS'):
    gams = [float(x) for x in os.environ['GAMS'].split(',')]
CAP = float(os.environ['CAP_S']) if os.environ.get('CAP_S') and MUT != 'nocap' else None
BASE = 'base (rtol1e-10, phi_cut1e-9, R60)'
PHI = 'phi_cut 1e-12 (larger r_m)'
VARIANTS = {
    BASE: dict(star='60'),
    'rtol 1e-12': dict(star='60', rtol=1e-12),
    PHI: dict(star='60', phi_cut=1e-12),
    'background R=90': dict(star='90'),
    'r0 = 3e-3': dict(star='60', r0=3e-3),
}
probs = {}


def init():
    st60 = B.solve_star(w, -1)
    st90 = B.solve_star(w, -1, guess=st60, R=90.0)
    for k, v in VARIANTS.items():
        v = dict(v)
        st = st60 if v.pop('star') == '60' else st90
        probs[k] = PolarProblem(st, l, **v)


def task(kg):
    k, g = kg
    T = time.time()
    if k == 'q':
        base = probs[BASE]
        zz, _ = base.interior_zerilli(g)
        f = 1 - 2 * base.Mext / base.rm
        return k, g, [f * zz[1] / zz[0], base.zerilli_q(g)], time.time() - T
    return k, g, probs[k].D(g), time.time() - T


def worker(tq, rq):
    """Build the problems once, then evaluate tasks until a None arrives, reporting each start and result."""
    init()
    while True:
        kg = tq.get()
        if kg is None:
            return
        rq.put(('start', kg, os.getpid(), time.time()))
        c0 = time.process_time()
        k, g, v, dt = task(kg)
        rq.put(('done', kg, os.getpid(), (v, dt, time.process_time() - c0)))


def evaluate(tasks, res, fp, nproc):
    """Run tasks in nproc workers; returns {task: (wall s, cap s or None, reason)} for those not completed."""
    stopped = {}
    if not tasks:
        return stopped
    tq, rq = mp.Queue(), mp.Queue()
    for kg in tasks:
        tq.put(kg)
    workers, running, left = {}, {}, len(tasks)

    def log(obj):
        fp.write(json.dumps(obj) + "\n"); fp.flush()

    def spawn():
        p = mp.Process(target=worker, args=(tq, rq), daemon=True)
        p.start()
        workers[p.pid] = p
    for _ in range(min(nproc, len(tasks))):
        spawn()
    while left:
        try:
            ev, kg, pid, pay = rq.get(timeout=1.0)
            kg = tuple(kg)
            if ev == 'start':
                running[pid] = (kg, pay)
                log({"start": list(kg), "pid": pid, "time": pay})
            else:
                v, dt, cpu = pay
                res[kg] = (v, dt)
                log([kg[0], kg[1], v, dt])
                log({"done": list(kg), "pid": pid, "wall": dt, "cpu": cpu})
                running.pop(pid, None)
                left -= 1
        except queue.Empty:
            pass
        now = time.time()
        for pid, (kg, t0) in list(running.items()):
            dead = not workers[pid].is_alive()
            if dead or (CAP is not None and now - t0 > CAP):
                if not dead:
                    workers[pid].terminate()
                workers[pid].join()
                code = workers.pop(pid).exitcode
                running.pop(pid)
                reason = f"its worker exited with code {code}" if dead else "cap"
                stopped[kg] = (now - t0, CAP, reason)
                log({"capped" if not dead else "died": list(kg), "pid": pid, "wall": now - t0, "cap": CAP,
                     "exitcode": code})
                left -= 1
                if left > len(running):
                    spawn()
    for _ in workers:
        tq.put(None)
    for p in workers.values():
        p.join(timeout=30)
    return stopped


def main():
    init()
    for k in VARIANTS:
        print(f"{k:38s}: r_m={probs[k].rm:.2f}")
    print()
    part = os.path.join(os.environ.get('PARTIAL_DIR', HERE), f"convergence_l{l}_w{w}.partial.csv")
    res, prev = {}, {}
    if os.path.exists(part):
        for line in open(part):
            rec = json.loads(line)
            if isinstance(rec, list):
                k, g, v, dt = rec
                res[(k, g)] = (v, dt)
            elif 'capped' in rec:
                prev[tuple(rec['capped'])] = (rec['wall'], rec['cap'], 'cap')
    keep = {kg: s for kg, s in prev.items() if kg not in res and CAP is not None and s[1] is not None and CAP <= s[1]}
    tasks = [(k, g) for g in sorted(gams, reverse=True) for k in list(VARIANTS) + ['q']
             if (k, g) not in res and (k, g) not in keep]
    with open(part, 'a') as fp:
        stopped = evaluate(tasks, res, fp, int(os.environ.get('NPROC', '1')))
    stopped.update(keep)
    for g in gams:
        Ds = {k: res[(k, g)][0] for k in VARIANTS if (k, g) in res}
        miss = [k for k in list(VARIANTS) + ['q'] if (k, g) not in res]
        T = sum(res[(k, g)][1] for k in list(VARIANTS) + ['q'] if (k, g) in res)
        same_rm = [Ds[k] for k in Ds if 'phi_cut' not in k]
        signs = sorted(set(np.sign(list(Ds.values()))))
        spread = f"{(max(same_rm) - min(same_rm)) / abs(np.mean(same_rm)):.1e}" if same_rm else "n/a"
        db = f"{Ds[BASE]:+.6e}" if BASE in Ds else "n/a"
        dp = f"{Ds[PHI]:+.3e}" if PHI in Ds else "n/a"
        qq = (f"q_int={res[('q', g)][0][0]:+.6e} q_ext={res[('q', g)][0][1]:+.6e}" if ('q', g) in res
              else "q_int=n/a q_ext=n/a")
        row = (f"gamma={g:8.3g}: D_base={db}  signs={signs}  rel.spread (same r_m)={spread}  D(r_m~34)={dp}  "
               f"{qq}  ({T:.0f}s)")
        if miss:
            why = []
            for k in miss:
                wall, cap, reason = stopped.get((k, g), (None, None, 'not run'))
                if reason == 'cap':
                    why.append(f"{k} not completed, stopped after {wall:.0f}s of wall time at the cap (CAP_S = {cap:.0f}s)")
                elif wall is not None:
                    why.append(f"{k} not completed, {reason} after {wall:.0f}s")
                else:
                    why.append(f"{k} not completed, not run")
            row += f"  [{'; '.join(why)}; signs and spread are over {len(Ds)} variants]"
        print(row, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
