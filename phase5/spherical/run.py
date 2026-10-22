"""Run one evolution with checkpoint/resume.  Output: runs/<name>.csv (time series) and
runs/<name>.ckpt.npz (state).  Usage:
  python3 run.py NAME STAR h R T [kind=exact] [delta=0] [dtout=1] [sko=0.2] [seamM=]
"""
import sys, os, time, json
import numpy as np
from evolve import Sim, load_star, NG

ROBS = [8.0, 12.0, 20.0, 30.0, 50.0]

def main():
    name, tag, h, R, T = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
    kw = dict(a.split('=') for a in sys.argv[6:])
    kind = kw.get('kind', 'exact'); delta = float(kw.get('delta', 0)); dtout = float(kw.get('dtout', 1.0))
    sko = float(kw.get('sko', 0.2)); seamM = (float(kw['seamM']), float(kw['jump'])) if 'seamM' in kw else None
    wall = float(kw.get('wall', 1e9))
    os.makedirs('runs', exist_ok=True)
    st = load_star(tag); eps = int(st['eps'])
    sim = Sim(st, h, R, eps, kind=kind, delta=delta, sko=sko, seam_M=seamM)
    csv = f'runs/{name}.csv'; ck = f'runs/{name}.ckpt.npz'
    if os.path.exists(ck):
        d = np.load(ck)
        sim.Phi[:] = d['Phi']; sim.Pi[:] = d['Pi']; sim.t = float(d['t']); sim.advance(0)
        print(f'resumed {name} at t={sim.t}', flush=True)
        f = open(csv, 'a')
    else:
        f = open(csv, 'w')
        f.write('# ' + json.dumps(dict(name=name, star=tag, h=h, R=R, T=T, kind=kind, delta=delta, sko=sko,
                                       seamM=seamM, eps=eps, w=float(st['w']), M_bvp=float(st['M']),
                                       N_bvp=float(st['N']), robs=ROBS)) + '\n')
        f.write('t,amp,phase,alpha_c,' + ','.join(f'm{int(r)}' for r in ROBS) + ',mR,' +
                ','.join(f'N{int(r)}' for r in ROBS) + ',maxphi,momres,content20\n')
    nout = max(1, int(round(dtout / sim.dt)))
    t0 = time.time(); lastck = t0
    def write():
        d = sim.diagnostics(ROBS)
        n = sim.n
        mr = float(np.max(np.abs(sim.mom_residual()))) if (round(sim.t / dtout) % 50 == 0) else float('nan')
        f.write(f"{d['t']:.6f},{d['amp']:.12e},{d['phase']:.12f},{d['alc']:.12f}," +
                ','.join(f'{x:.12e}' for x in d['m']) + f",{d['mR']:.12e}," +
                ','.join(f'{x:.12e}' for x in d['N']) + f",{d['maxphi']:.10e},{mr:.3e},{d['content']:.12e}\n")
    if sim.t == 0.0:
        write()
    while sim.t < T - 1e-9 and time.time() - t0 < wall:
        sim.advance(nout)
        write()
        if time.time() - lastck > 60:
            f.flush(); np.savez(ck, Phi=sim.Phi, Pi=sim.Pi, t=sim.t); lastck = time.time()
    f.flush(); f.close()
    np.savez(ck, Phi=sim.Phi, Pi=sim.Pi, t=sim.t)
    # final profile
    n = sim.n
    np.savez(f'runs/{name}.final.npz', r=sim.r, Phi=sim.Phi[NG:NG + n], Pi=sim.Pi[NG:NG + n], m=sim.m,
             lnal=sim.lnal, t=sim.t)
    print(f'{name} done t={sim.t:.2f} wall={time.time()-t0:.0f}s', flush=True)

if __name__ == '__main__':
    main()
