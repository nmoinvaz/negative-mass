"""dM/domega, dN/domega at the campaign stars from the cached BVP neighbours (centred differences)."""
import numpy as np, json
def L(t): d = np.load(f'stars/{t}.npz'); return float(d['w']), float(d['M']), float(d['N'])
out = {}
for name, lo, mid, hi in [('c0750', 'c0740', 'c0750', 'c0760'), ('p07603', 'p07503', 'p07603', 'p07703')]:
    (w0, M0, N0), (w1, M1, N1), (w2, M2, N2) = L(lo), L(mid), L(hi)
    dM = (M2 - M0) / (w2 - w0); dN = (N2 - N0) / (w2 - w0)
    out[name] = dict(w=w1, M=M1, N=N1, dMdw=dM, dNdw=dN, dlnMdw=dM / M1, dlnNdw=dN / N1,
                     first_law_ratio=dM / (w1 * dN))
    print(name, json.dumps(out[name]))
json.dump(out, open('family_slopes.json', 'w'), indent=1)
