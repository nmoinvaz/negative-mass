"""Family slopes of Shirokov's pack diagnostics, to test whether his reported drifts are a quasi-static
slide along the static family.  Pack diagnostics (review/acceleration_drift): peak = max sqrt(phi^2+Pi^2)
and content = int sqrt(phi^2 + Pi^2) dV with |Pi| = w phi / alpha (stationary star, real-field
components; phi here = |Phi|).  Coordinate volume approximated by areal 4 pi r^2 dr (isotropic vs areal
differs by O(M/r) ~ 3e-3 and cancels in log-slopes to that order)."""
import numpy as np, json
def diag(t):
    d = np.load(f'stars/{t}.npz'); r, phi, nu = d['r'], d['phi'], d['nu']; w = float(d['w'])
    al = np.exp(nu / 2)
    mod = phi * np.sqrt(1 + (w / al)**2)
    content = np.trapezoid(4 * np.pi * r**2 * mod, r)
    return w, float(d['M']), float(d['N']), float(phi[0]), float(mod.max()), float(content)
res = {}
for name, lo, mid, hi in [('c0750', 'c0740', 'c0750', 'c0760'), ('p07603', 'p07503', 'p07603', 'p07703')]:
    a, b, c = diag(lo), diag(mid), diag(hi)
    dw = c[0] - a[0]
    sl = {k: (c[i] - a[i]) / dw / b[i] for i, k in enumerate(['w', 'M', 'N', 'phic', 'peak', 'content']) if i > 0}
    res[name] = dict(values=dict(zip(['w', 'M', 'N', 'phic', 'peak', 'content'], b)), dln_dw=sl)
    print(name, json.dumps(res[name]))
# Shirokov / pack lone canonical + pair drifts (review/acceleration_drift table): Q1->Q4 omega, content t=0->200, peak
pack = {'N128': (0.7451, 0.7148, 0.210, -0.0036), 'N192': (0.7466, 0.7218, 0.161, -0.0034),
        'N256': (0.7479, 0.7283, 0.124, -0.0025), 'lone N128': (0.7433, 0.7135, 0.200, -0.0048)}
sl = res['c0750']['dln_dw']
print('\nquasi-static prediction for the pack canonical drifts (linearised in d omega):')
print('| grid | Δω | predicted Δpeak | reported Δpeak | predicted Δcontent | reported Δcontent | predicted ΔM/M |')
for g, (w0, w1, cont, peak) in pack.items():
    d = w1 - w0
    print(f"| {g} | {d:+.4f} | {sl['peak']*d:+.2%} | {peak:+.2%} | {sl['content']*d:+.2%} | {cont:+.1%} | {sl['M']*d:+.1%} |")
json.dump(res, open('family_diag.json', 'w'), indent=1)
