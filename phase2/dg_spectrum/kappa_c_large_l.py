"""kappa_c(R, l) = d ln phi0/dx |_R for l up to 20 (exact formula; see truncated.py docstring).
Also checks the large-R asymptote kappa_c -> -(l+1)/R and the Robin lower bound E >= -kappa^2
valid whenever V^ >= 0 on the exterior (R > r_neg(l))."""
import numpy as np
import common as c
Ls = list(range(2, 13)) + [16, 20]
Rs = (0.05, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0)
print("kappa_c(R, l)  [unstable mode exists iff kappa < kappa_c, for R > r_s(l)]")
print(f"{'l':>3} " + " ".join(f"{'R='+str(R):>12}" for R in Rs))
rows = []
for l in Ls:
    vals = [float(c.numeric(l)['kappa_c'](R)) for R in Rs]
    rows.append([l] + vals)
    print(f"{l:>3} " + " ".join(f"{v:>12.4g}" for v in vals) + ("   (R=0.05 < r_s)" if 0.05 < c.r_s(l) else ""))
print("\nkappa_c * R  (large-R limit -(l+1)):")
for row in rows:
    print(f"{row[0]:>3} " + " ".join(f"{v*R:>12.4g}" for v, R in zip(row[1:], Rs)))
import csv, os
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kappa_c_large_l.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['l'] + [f'kappa_c_R={R}' for R in Rs]); w.writerows(rows)
