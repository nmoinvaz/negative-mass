"""
Localisation of the Gleiser-Dotti unstable mode near the negative-mass singularity.

The metric amplitude of the mode (K, H1, H2) carries the radial factor
    (r - 2M)^k exp(-k r / (2|M|)),   M = -|M|,   k = (l-1) l (l+1)(l+2) / 6,
which is maximal at r = 0 and falls off on a scale ~ 2 sqrt(2/k) |M|
(from ln(1+x) - x ~ -x^2/2 with x = r / 2|M|). The growth rate is k / (2|M|).
Units: |M| = 1.
"""
import numpy as np

r = np.linspace(1e-6, 20, 200001)
print(f"{'l':>2} {'k':>5} {'growth*|M|':>11} {'1/e radius/|M|':>15} {'2*sqrt(2/k)':>12}")
for l in range(2, 8):
    k = (l - 1) * l * (l + 1) * (l + 2) / 6
    ln_amp = k * np.log1p(r / 2) - k * r / 2  # ln[(r+2)^k e^{-kr/2}] - ln[2^k]
    w = r[np.argmax(ln_amp < -1)]
    print(f"{l:>2} {k:>5.0f} {k / 2:>11.1f} {w:>15.4f} {2 * np.sqrt(2 / k):>12.4f}")
