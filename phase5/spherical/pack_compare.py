"""Nonlinear version of family_diag: static canonical stars at the pack's end-of-run frequencies.
If the pack's lone/pair stars slid quasi-statically along the family from omega_Q1 to omega_Q4, their
peak modulus, content and N would change by these amounts.  Compare with the reported changes."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'phase2', 'radial'))
import background as bg
def diag(st):
    r = np.linspace(st.r0, st.R, 200001); y = st.sol.sol(r)
    phi = np.exp(y[0]); al = np.exp(y[3] / 2); w = st.w
    mod = phi * np.sqrt(1 + (w / al)**2)
    return dict(M=st.M, N=st.N, phic=st.phic, peak=float(mod.max()), content=float(np.trapezoid(4 * np.pi * r**2 * mod, r)))
pack = {'N128': (0.7451, 0.7148, 0.210, -0.0036), 'N192': (0.7466, 0.7218, 0.161, -0.0034),
        'N256': (0.7479, 0.7283, 0.124, -0.0025), 'lone N128': (0.7433, 0.7135, 0.200, -0.0048)}
cache = {}
def get(w):
    if w not in cache:
        cache[w] = diag(bg.solve_star(w, 1))
    return cache[w]
print('| grid | ω Q1→Q4 | static ΔM/M | static ΔN/N | static Δpeak | reported Δpeak | static Δcontent | reported Δcontent (t 0→200) |')
print('|---|---|---|---|---|---|---|---|')
for g, (w0, w1, cont, peak) in pack.items():
    a, b = get(w0), get(w1)
    rel = lambda k: b[k] / a[k] - 1
    print(f"| {g} | {w0}→{w1} | {rel('M'):+.1%} | {rel('N'):+.1%} | {rel('peak'):+.2%} | {peak:+.2%} | {rel('content'):+.1%} | {cont:+.1%} |", flush=True)
