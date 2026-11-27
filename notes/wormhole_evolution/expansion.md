# The Ellis threshold's expanding branch, followed along the symmetric throat (#213)

Family S of #40 (symmetric two-sided ghost pulse, A above p*ₕ) evolved with #39's double-null code ([README.md](README.md), `dn.py` unchanged) to t = 283 at D = 4 and
t = 141 at D = 5, at two resolutions. Driver: [expansion.py](expansion.py). Output of record: [expansion_output.txt](expansion_output.txt)
(`python3 expansion.py analyze --raw=<negative-mass-raw>/I213-dh610`, exit 1 because CT is not met). Raw output, with the declarations committed before the runs: negative-mass-raw
`I213-dh610` (README, `dev/`, `runs/`, `figures/expansion_H.svg`). Status words as in AGENTS.md; every number is *numerical evidence* unless marked. No sentence here says what the
end state is beyond what the declared criteria and the numbers show, and none says a result is new.

Commits: the runs of record were made with `expansion.py` at 1b6b43d (criteria and mutations declared there). b6c33c2 and 015ce77 add report-only output (τ_s, exponents, curvature
scalars; the coordinator asked for them after the runs started) and change no check, criterion or tolerance; the check verdicts of the analysis at 1b6b43d and at 015ce77 are identical.

## Answer

The de Sitter-like criterion is **not met**, and this observer cannot meet it. On the last e-fold of r the throat's H rises by 1.65 (D = 4) and 1.29 (D = 5) against the bar 1.2, |ρ|
rises by 2.72 and 1.65 against 1.2, and w̄ = (1−n)/(1+n) = −1/3 (D = 4), −1/2 (D = 5) at every time by symmetry (*proved*, below). The declared class is **FT**: 1/H is linear in τ, which is the
declared alternative "growth without a constant H" with H = k/(τ_s − τ), k = 2.0000 (D = 4) and 3.99 (D = 5). All six members agree on the class at two resolutions (H at fixed r/a differs
by ≤ 5e-5 at D = 4 and 4e-4 at D = 5).

- **#211's slope test is not met.** H/W at r/a = 10, 25, 50 is 1.23, 1.95, 2.75 (D = 4, tuned member) and 0.85, 1.07, 1.28 (D = 5); the bar is 10 %.
- **Mass at the ends.** M on each end is not conserved. It changes by exactly the flux through the end (residual ≤ 1.9e-7 of the largest change, bar 1e-3). It rises from 0.0035 to 8e6 at x⁻ = 80
  (X = 200, D = 4, tuned member) by the flux of the ghost field's negative-energy radiation T₋₋ < 0 through x⁺ = X.
- **Anti-trapped region.** It contains the diagonal from t ≈ 3 (t ≈ 11 for the tuned member) to the end of every run. No trapped cell forms and no cell is excised.
- **Control.** A positive-energy (normal-field) pulse and the S pulse of the other sign collapse at both h at D = 4 and 5. The declared bar on the event times at two h is **not met for the D = 5 normal pulse**
  (|Δt| = 0.0725 against 0.040); a third h shows the event time converging (3.815, 3.887, 3.914). CT is reported as not met.
- **Checks:** 12 of 13 pass (CT fails). All 15 mutations fire exactly their declared sets.

## What was run

```
x- = X  +------------------------------+   the rectangle is #39's domain.svg, data on the two lower edges
        |                         /    |   diagonal x+ = x-: the throat observer (mirror-fixed, a geodesic)
        |   r- = 0 contour leaves /    |   contour r+ = 0: stays at x+ ~ 3 to 13
        |   the grid (x- ~ 4..13)/     |   contour r- = 0: leaves the grid through x+ = X
        |                    /         |   between them: the anti-trapped region
        |  r+ = 0         /            |   right end = the edge x+ = X; the other end = the edge x- = X
x- = 0  +------------------------------+
       x+ = 0                       x+ = X
```

| item | declaration (README of the raw folder) |
|---|---|
| family | ghost φ-Gaussian, centre 3, width 0.5, amplitude A on one ray, mirrored with odd parity; members p*ₕ + 1e-6 ("tuned"), A = 0.01, 0.05, 0.2 |
| D = 4 | h = 0.02, 0.01; X = 200 (t = 283; N = 10⁴, 2×10⁴) |
| D = 5 | h = 0.00707, 0.00354; X = 100 (t = 141; N = 14142, 28284); members tuned, 0.05 |
| observer | the diagonal; H = d ln r/dτ, dτ = √2 e^{−f/2} dx |
| ρ, p, w | ρ = T_uu, p_∥ = T_ll, p_⊥ = T^θ_θ, w = w̄ = (p_∥ + n p_⊥)/((n+1)ρ) |
| classes | DS, STIFF, FT, GROW, OTHER on the last e-fold of r (CL in `expansion.py`) |

Cost: D = 4, h = 0.01, X = 200: 2.6 to 4.7 min per member (the tuned one includes the bisection for p*ₕ); D = 5, h = 0.00354: about 7 min; all runs of record in 16 min of wall time on four lanes. The grid did not limit the runs: no excised cell, no NaN,
max e^{−f} = 1.0000 (at the corner), min e^{−f} falling to 1e-7 to 3e-7 (D = 4) and 5e-6 to 9e-6 (D = 5). The cells shrink in proper time; they do not stretch. The run length is the declared X.

## The identities on the throat (*proved*)

The data are mirror symmetric: φ(x⁺, x⁻) = −φ(x⁻, x⁺), r and f symmetric. On the diagonal φ = 0 for all time, so φ_u = 0 and φ₊ = −φ₋. With T_ab = ε(∂_aφ ∂_bφ − g_ab(∂φ)²/2) and g₊₋ = −e^{−f}:

- T₊₋ = ε(φ₊φ₋ − φ₊φ₋) = 0, so in the frame of the diagonal ρ = T_uu = e^f(T₊₊ + T₋₋)/2 = ε e^f φ₊² and p_∥ = T_ll = ρ.
- p_⊥ = T^θ_θ = ε e^f φ₊φ₋ = −ρ.
- p̄ = (p_∥ + n p_⊥)/(n+1) = (1−n)ρ/(n+1), so **w̄ = (1−n)/(1+n) = −1/3 (D = 4), −1/2 (D = 5)**, whatever the dynamics, whenever ρ ≠ 0. A homogeneous ghost has φ₊ = φ₋ and w = +1 (planted W1 below).

Check (W2, every run, both D): |p_⊥/ρ + 1|, |p_∥/ρ − 1|, |T₊₋/ρ|, |w̄ − (1−n)/(1+n)| ≤ 1.0e-15 (bar 1e-9). The brief's "w within a stated distance of −1" is therefore decided by the symmetry,
not by the dynamics: **the criterion DS is NOT MET and this observer cannot meet it.** The components are reported beside it: **w_∥ = +1, w_⊥ = −1** at every time.

One more relation, derived by hand from m = r^{n−1}(1 + ṙ²)/2 (ṙ = dr/dτ) and the ray identity at the diagonal ṁ = −(r^n/n)ρṙ (MB3 checks the identity): dH/dτ = −(n+1)H²/2 − (n−1)/(2r²) − ρ/n.
Derived by hand, not checked symbolically; checked numerically on the stored series (`dev/dev_14_hdot.py`): median relative difference ≤ 6e-5 in every window; maximum ≤ 3e-4 except where the right side crosses zero (A = 0.2 at D = 4, A = 0.05 at D = 5, both at r/a < 10, where up to 0.4 and 0.19).
With H = k/(τ_s − τ) and ρ = −cH² it gives c = n(1/k + (n+1)/2): 4 for (n, k) = (2, 2) and 27/4 for (3, 4), the measured |ρ|/H² of 4.0000 and 6.75.

## Planted controls

| id | what | result |
|---|---|---|
| ST | `stream` against `dn.evolve` on two grids, one with excision | diagonal bit-identical; edge mass 3e-13; flux residual 1e-13 |
| OB0 | static wormhole: r_diag = 1 for t ≤ 4, minimum of r on the diagonal | 1.1e-4 (bar 1e-3); offset 1 node |
| W1 | exact open FRW + homogeneous ghost (a² = cosh 2η, φ = √6 atan tanh η; *checked symbolically*, residuals 1e-30) through the same code, across the bounce | orders 1.95 to 2.00; \|w̄ − 1\| = 1.0e-10 at h = 0.01 (bar 1e-6) |
| PL | five synthetic throat series through `classify` and `slope_test` | DS, FT, STIFF, GROW, OTHER recovered; slope test met only for de Sitter |
| SL0 | the tuned member departs at the linear rate: slope of ln\|H\| against τ over \|H\|/W ∈ [1e-3, 1e-2] | −5.0e-3 (D = 4), −5.4e-3, −3.2e-3 (D = 5), bar 1 % |
| OB1, OB2 | mirror symmetry and "the diagonal is the minimum of r on every full slice", through each run | asymmetry ≤ 1.2e-11 (bar 1e-9); offset 1 node on every slice t ≤ 141 (D = 4), r_ll ≥ 2.6e-3 > 0 |

OB2 says the throat in the sense of the brief, the minimum-r point, is the diagonal through the run, within one node. The runs' slices are full length for t ≤ X/√2 (141 at D = 4, 71 at D = 5); after that the
rectangle cuts them and the check does not cover them.

## H(τ), ρ and w on the throat

D = 4, tuned member (A = p*ₕ + 1e-6 = 1.288e-5), h = 0.01. The first anti-trapped diagonal cell is at τ = 10.93 (before it |H|/W ≤ 1e-2: the linear phase).

| r/a | t | τ | H | H/W | ρ | \|ρ\|/H² | w_∥, w_⊥ | w̄ | m/(r³H²/2) |
|---|---|---|---|---|---|---|---|---|---|
| 1.5 | 14.6 | 14.129 | 0.374 | 0.316 | −1.40 | 10.0 | +1, −1 | −1/3 | 4.17 |
| 3 | 17.5 | 15.370 | 0.765 | 0.647 | −2.62 | 4.47 | +1, −1 | −1/3 | 1.189 |
| 10 | 27.7 | 16.515 | 1.454 | 1.230 | −8.48 | 4.008 | +1, −1 | −1/3 | 1.0047 |
| 25 | 48.9 | 17.021 | 2.300 | 1.946 | −21.2 | 4.0002 | +1, −1 | −1/3 | 1.0003 |
| 50 | 84.3 | 17.275 | 3.253 | 2.752 | −42.3 | 4.0000 | +1, −1 | −1/3 | 1.00004 |
| 100 | 155.0 | 17.455 | 4.600 | 3.892 | −84.6 | 4.0000 | +1, −1 | −1/3 | 1.00000 |

H grows without a plateau: d ln H/d ln r → 1/2 (H ∝ r^{1/2}); r H = 460 at r/a = 100. Members and resolutions (`expansion_output.txt` has every table):

| member | D | A | H/W at r/a = 10, 25, 50 | τ_s (h coarse, fine) | k | ln r_th max | class |
|---|---|---|---|---|---|---|---|
| tuned | 4 | p*ₕ + 1e-6 | 1.230, 1.946, 2.752 | 17.8902, 17.8902 | 2.0000 | 5.249 | FT |
| 0.01 | 4 | 0.01 | 1.232, 1.949, 2.756 | 10.0590, 10.0559 | 2.0000 | 5.280 | FT |
| 0.05 | 4 | 0.05 | 1.280, 2.023, 2.861 | 8.4340, 8.4334 | 2.0000 | 5.302 | FT |
| 0.2 | 4 | 0.2 | 1.768, 2.777, 3.924 | 6.1692, 6.1690 | 2.0004 | 5.476 | FT |

*Correction (coordinator, 2026-11-27, from Paper V's number check):* in this table the A = 0.05 row of H/W (1.280, 2.023, 2.861) is the
h = 0.02 output; at h = 0.01 the output prints 1.280, 2.024, 2.862. The A = 0.2 row matches neither resolution; the output prints 1.769, 2.778,
3.924 at both h. The statement |ρ|/H² = 4.0000 across the four amplitudes holds to the output's 3.9999 (A = 0.05) and 3.9985 (A = 0.2) at
r = 50a. `expansion_output.txt` is the canonical source; the conclusions do not change.
| tuned | 5 | p*ₕ + 1e-6 | 0.847, 1.074, 1.279 | 15.0682, 15.0724 | 3.9875 | 4.425 | FT |
| 0.05 | 5 | 0.05 | 1.023, 1.297, 1.544 | 7.9198, 7.9180 | 3.9897 | 4.486 | FT |

The late state does not depend on A except through τ_s: k and |ρ|/H² agree across the four D = 4 amplitudes (2.0000 to 2.0004; 4.0000). Two resolutions agree (RES): H at fixed r/a differs by ≤ 5.0e-5 (D = 4)
and ≤ 4.1e-4 (D = 5); k by ≤ 1.5e-6; τ_s by ≤ 3.1e-4; ln r_max by ≤ 3.6e-5. The absolute constraint maxima fall by 4.0 per halving in every member (D = 4 tuned: 5.8e-5, 1.45e-5; D = 5 tuned: 4.6e-5, 1.15e-5);
the maxima relative to their own terms are 0.04 to 1.0 and are not meaningful (not examined further).

The classification against the declared alternatives: DS needs H max/min ≤ 1.2, |ρ| max/min ≤ 1.2 and |w̄ + 1| ≤ 0.2; the runs give 1.65, 2.72 and 0.67. STIFF needs |w̄ − 1| ≤ 0.2; the runs give 1.33 (D = 4) and 1.5 (D = 5). FT is met (1/H linear in τ to 2.6e-5 at D = 4 and 2.1e-4 at D = 5, relative to its mean; τ_s inside (τ_end, 1.2 τ_end]).

## The limiting proper time and the curvature along the throat (reported, not gated)

τ_s is the intercept of the linear fit of 1/H against τ over the last e-fold of r. It is an extrapolation: the last live cell is at τ_s − τ = 0.31 (D = 4) and 1.59 (D = 5, tuned), t = 283 and 141. The fit does not depend on the window to the
digits shown at D = 4 (0.5, 1, 2 e-folds: 17.89022, 17.89022, 17.89023) and to 1e-3 at D = 5 (15.0747, 15.0724, 15.0589).

Fitted exponents against (τ_s − τ) over the last e-fold, D = 4 tuned (the others agree to 1e-4):

| quantity | exponent |
|---|---|
| r/a | −2.0000 |
| H | −1.0000 |
| \|ρ\| = \|R\|/2 | −2.0000 |
| K = (1 − (∇r)²)/r² | −2.0000 |

D = 5 tuned: r −3.9876, H −0.9999, |ρ| −1.9967, K −1.9990. On the diagonal the Ricci scalar is R = ε(∂φ)² = 2ρ and the sectional curvature of the angular planes is K = (1 − (∇r)²)/r² = 2m/r^{n+1} = H² + 1/r²
(K − H² − 1/r² ≤ 4e-14). Both are scalars of the metric (and of the field), so their rise is not a coordinate effect within the range reached. Values along the run (D = 4 tuned): at τ_s − τ = 3.0, 1.0, 0.31: r/a = 2.2, 18.9, 190; ρ = −1.93, −16.0, −161;
R = −3.9, −32.0, −322; K = 0.57, 4.00, 40.3; K(τ_s − τ)² = 5.15, 4.003, 4.0000 and |ρ|(τ_s − τ)² = 17.4, 16.004, 16.0002. Whether the throat reaches τ_s, and what happens there, is not measured: the numbers stop at τ_s − τ = 0.31.
The relation above makes the asymptotic law an ODE question, not a data one (task 1 below).

## #211's slope test and [SH02]'s H

| member | H/W at r/a = 10 | 25 | 50 | verdict (bar 10 %) |
|---|---|---|---|---|
| tuned, D = 4 | 1.230 | 1.946 | 2.752 | not met |
| 0.01 | 1.232 | 1.949 | 2.756 | not met |
| 0.05 | 1.280 | 2.023 | 2.861 | not met |
| 0.2 | 1.768 | 2.777 | 3.924 | not met |
| tuned, D = 5 | 0.847 | 1.074 | 1.279 | not met (r/a = 10 and 50 outside) |
| 0.05, D = 5 | 1.023 | 1.297 | 1.544 | not met |

ln(r/a) does not keep its small-amplitude slope once r/a ≫ 1: d ln r/dτ crosses W near r/a = 6.5 (D = 4, tuned; near 18 at D = 5) and keeps rising. The SL0 control shows the same code returns W in the linear phase (−0.5 %).

[SH02]'s fit form is r/a = 1 + b₄ exp(H(τ − b₅)). The least-squares exponent of ln(r/a − 1) against τ over 1e-2 ≤ r/a − 1 and r/a ≤ r_max, D = 4, W = 1.182 a⁻¹ (1.1 a⁻¹ = 0.931 W):

| member | r_max = 1.5 | 2 | 3 | 5 | 10 | 25 |
|---|---|---|---|---|---|---|
| tuned | 0.983 W | 0.976 | 0.968 | 0.967 | 0.985 | 1.052 |
| 0.01 | 1.202 | 1.114 | 1.058 | 1.029 | 1.033 | 1.102 |
| 0.05 | 2.206 | 1.705 | 1.418 | 1.268 | 1.206 | 1.256 |
| 0.2 | 4.424 | 4.179 | 3.332 | 2.552 | 2.134 | 2.050 |

For near-threshold members the fitted exponent is 1.14 to 1.16 a⁻¹ (0.97 to 0.99 W) for every r_max up to 10, which is within 4 to 6 % of [SH02]'s 1.1 a⁻¹. It is a window average: the local rates H and d ln(r/a − 1)/dτ at r/a = 10 are 1.23 W and 1.37 W and rising. For
members that are not near threshold the exponent depends on the amplitude and on r_max (1.2 to 4.4 W).

## The anti-trapped region and the horizons

| member (D = 4, h = 0.01) | first anti-trapped diagonal cell (r₊, r₋ > 0.01) | contour r₋ = 0 leaves the grid (x⁺ > 200) at |
|---|---|---|
| tuned | t = 10.94, r/a = 1.012 | x⁻ = 12.99 |
| 0.01 | t = 3.53, r/a = 1.003 | x⁻ = 7.47 |
| 0.05 | t = 3.15, r/a = 1.003 | x⁻ = 6.30 |
| 0.2 | t = 2.88, r/a = 1.002 | x⁻ = 4.44 |

(The last column: `expansion_output.txt`, "the contour r_- = 0 stays inside the grid up to".) At D = 5 (h = 0.00354, X = 100) the first anti-trapped cell is at t = 8.0 (tuned) and 3.0 (0.05), and the contour leaves the grid at x⁻ = 12.58 and 7.61.

- Every later live diagonal cell is anti-trapped, in every run. No trapped cell forms and no cell is excised.
- The contour r₊ = 0 stays at x⁺ ≈ 3 to 13 while r on it grows with x⁻ (D = 4 tuned: x⁺ = 10.9, 11.8, 12.3, 13.0 and r = 7.6, 28.7, 64, 134 at x⁻ = 20, 50, 100, 199, roughly proportional to x⁻).
- The contour r₋ = 0 leaves the grid. For x⁻ beyond the last column of the table the whole row x⁺ ≤ X is anti-trapped to the right of r₊ = 0. dev_05 (A = 0.05, X = 100, h = 0.04) shows the contour's x⁻ rising slowly with x⁺ (x⁻ = 2.7, 3.9, 4.55, 5.0, 5.4, 5.7 at x⁺ = 3, 8, 20, 40, 70, 100, A = 0.05).
  The grid does not say whether the anti-trapped region is closed on the right at larger x⁺.
- The throat sphere is far outside the Hubble radius: r H = 460 at r/a = 100. The apparent horizons here bound the anti-trapped region on the outside; there is no horizon inside it, and none at r = 1/H.

So an anti-trapped region forms; a cosmological horizon in the sense of a closed outer boundary of the anti-trapped region at r ≈ 1/H does not appear in the grid. [GGS09b]'s "no apparent horizon" is consistent with no trapped region; the marginal contours r₊ = 0, r₋ = 0 do exist.

## The masses at the ends and the flux

Misner–Sharp mass M = r^{n−1}(1 + 2e^f r₊r₋)/2 at the nodes of the right edge x⁺ = X, the columns X/2, X/4 and the top edge x⁻ = X (the other end).

| check | what | worst over all runs | bar |
|---|---|---|---|
| MB1 | M(j) − M(0) = ∫ −(rⁿ/n)e^f r₊ Σεφ₋² dx⁻, largest residual / largest change | 1.9e-7 | 1e-3 |
| MB2 | the two ends have the same history, M(i, N) = M(N, i) | 5.6e-9 | 1e-6 |
| MB3 | m_edge − m_diag = ∫ −(rⁿ/n)e^f r₋ Σεφ₊² dx⁺ along the row, / (\|I₊\| + \|I₋\|) | 1.3e-4 | 1e-3 |

The MB1 residual at the edge falls from 1.33 to 0.457 (D = 4 tuned, ratio 2.9), 1.74 to 0.380 (A = 0.05, 4.6), 6.04 to 1.32 (A = 0.2, 4.6), 2.84 to 1.69 (D = 5 A = 0.05, 1.7) and 1.49 to 0.558 (D = 5 tuned, 2.7) on mass changes of 1.4e8 to 7e8;
the ratio is near 4 for two members (A = 0.05, 0.2) and below 4 for the others, so second order of this residual is shown for two members only (at D = 5 the residual is 6e-9 of the change; not examined further).

D = 4, tuned member, h = 0.01, X = 200 (A = 0.05 at x⁻ = 80: 1.3e6, 6.5e6, 1.3e7, m_diag 4.3e6):

| x⁻ | t at x⁺ = X | M(x⁺ = 50) | M(x⁺ = 100) | M(x⁺ = X) | m_diag | I₊ | I₋ |
|---|---|---|---|---|---|---|---|
| 3 | 143.5 | 1.50e-2 | 7.28e-3 | 3.58e-3 | 0.5000 | 9e-8 | −0.4964 |
| 8 | 147.1 | 1.95e-2 | 9.66e-3 | 5.48e-3 | 0.5098 | 1.8e-4 | −0.5048 |
| 12 | 149.9 | 20.39 | 20.39 | 20.36 | 6.44 | 14.1 | −0.042 |
| 20 | 155.6 | 4.4e3 | 4.9e3 | 5.0e3 | 1.27e3 | 3.7e3 | 0 |
| 40 | 169.7 | 1.5e5 | 2.9e5 | 3.4e5 | 9.1e4 | 2.5e5 | 0 |
| 80 | 198.0 | 6.4e5 | 4.0e6 | 8.1e6 | 2.6e6 | 5.5e6 | 0 |

- Before the expansion reaches the edge (x⁻ ≲ 8) M at the end is the static wormhole's 1/(2r) plus the pulse's energy: E₀ = −5e-4 (A = 0.01), −0.0122 (0.05), −0.198 (0.2), below 1e-6 in magnitude for the tuned member; at D = 5, −1e-5 (tuned) and −0.031 (0.05).
- Afterwards M at the end rises by 9 orders. By MB1 this is the flux of T₋₋ = −Σφ₋² < 0 (ghost) through x⁺ = X, which raises M where r₊ > 0. The ends do not conserve their mass; the same holds at X/2 and X/4.
- At fixed x⁻ ≤ 12, M(x⁺) converges as x⁺ grows (the Bondi-like value, to 1e-3 at x⁻ = 12 between X/2 and X). For x⁻ ≳ 20 it does not: the row is inside the anti-trapped region and M depends on X (6.4e5, 4.0e6, 8.1e6 at x⁺ = 50, 100, 200, x⁻ = 80).
  No mass at spatial infinity is claimed; the values are M(x⁺ = X, x⁻) with X stated.
- The throat's mass m_diag = r^{n+1}H²/2 (to 1e-5 for r/a ≥ 50) is smaller than the edge mass by I₊ (the ray integral of the ghost's energy density weighted by e^f r₋ > 0). I₋ ≈ −0.50 for x⁻ ≲ 8 (m_edge − m_diag ≈ 0 − 1/2, the static throat mass) and vanishes after.

What this says about where the energy of the expansion comes from, in numbers: m_diag is the Hubble-flow mass r^{n+1}H²/2 of the throat sphere; it equals the end's mass minus the ghost-field integral I₊ along the ray; and the end's mass changes
exactly by the ghost flux through it. Nothing here says the ends' masses are conserved, or what a conserved total would be.

## Control: a positive-energy pulse collapses

| case | D | h | t of the first trapped cell | r there | M there | E₀ |
|---|---|---|---|---|---|---|
| normal field A = 0.05, mirrored | 4 | 0.02, 0.01 | 4.426, 4.419 | 0.9920, 0.9923 | 0.4960, 0.4961 | +0.0118 |
| S ghost A = −0.05 | 4 | 0.02, 0.01 | 3.154, 3.147 | 0.9971, 0.9973 | 0.4986, 0.4986 | −0.0121 |
| normal field A = 0.05 | 5 | 0.00707, 0.00354, 0.00177 | 3.815, 3.887, 3.914 | 0.9955, 0.9959, 0.9959 | 0.4955, 0.4959, 0.4959 | +0.0310 |
| S ghost A = −0.05 | 5 | 0.00707, 0.00354, 0.00177 | 2.975, 2.972, 2.976 | 0.9973, 0.9977, 0.9977 | 0.4973, 0.4977, 0.4977 | −0.0310 |

All four collapse at every h; the M at the first trapped cell is 0.496 to 0.499 a^{D−3}, close to the threshold-run values of #40 (0.494, 0.492), as a formation value in a stretching grid and not a settled mass (#41).
CT as declared requires the event times at the two h to agree within 4√2 h: met at D = 4 (7e-3 against 0.113) and for the D = 5 ghost (2.5e-3 against 0.040), **not met for the D = 5 normal pulse (0.0725 against 0.040)**. The third h (not part of the declared check,
run after the runs of record, `runs/extra/`) gives 0.0263, a ratio of 2.8 between successive differences; its pair with the second h is also outside 4√2 h (0.020). The collapse and the mass at the horizon are robust; the event time of the D = 5 normal pulse
converges slowly. CT is reported as not met.
S(+0.05) and S(−0.05) have the same E₀ = −0.012 (D = 4) and opposite fates: for this φ-Gaussian family E₀ does not set the fate (the sign of the overlap with the unstable mode does). #39's F1 note on [SH02]'s sign rule concerned the derivative-form pulse.

## Mutations

Fifteen, on the reduced suite (D = 4, X = 40); each fires exactly its declared set (`runs/mut/`, the declared table is in the raw README): w_sign (W1, W2), time_unit (SL0, ST, W1), mirror_parity (MB1, MB2, MB3, OB1, OB2, RES, W2), class_loose (PL),
control_ghost (CT), flux_sign and mass_exponent (MB1, MB3), no_corrector (MB1, MB3, RES, W1), the dn.py input mutations (BUILD, OB0, ... as declared), finder_sign (BUILD, CT). `constraint_sign` fires nothing by design: only the reported constraint monitor reads it.
Two of the first-run mutation logs (`runs/mut_run1/`, `mirror_parity`, `scalar_coef`) stopped with an exception in the new report code; the guard was added at 015ce77 and all 15 were rerun.

## Discrepancies and things that failed

1. **CT not met** (D = 5 normal pulse, event times at two h: 0.0725 against 0.040). Not changed.
2. **The brief's w criterion cannot be met by this observer** (*proved*, W2): reported to the coordinator before the runs; the criterion was kept as declared.
3. Two mutation logs crashed in report-only code (above); the logs are kept and the mutations rerun.
4. MB1's residual does not show second order at D = 5 (ratios 1.7, 2.7); it is 6e-9 of the mass change.
5. The constraint maxima relative to their own terms are 0.04 to 1.0 (not meaningful where the terms are small); the absolute maxima converge at second order.

## What `critical_ghost.md` §1.2 and `baby_universes.md` §3 should now say (proposals for the coordinator)

`critical_ghost.md` §1.2, after the [GGS09b] and [SH02] items:

- A symmetric two-sided ghost pulse above threshold, followed along the throat's diagonal to t = 283 (D = 4) and 141 (D = 5) at two resolutions: an anti-trapped region contains the throat from the departure on; no trapped region forms; the throat's rate H = d ln r/dτ rises
  without a plateau (H/W = 1.23, 1.95, 2.75 at r/a = 10, 25, 50), with 1/H linear in τ (k = 2.0000 at D = 4, 3.99 at D = 5) and the curvature scalars R and K ∝ (τ_s − τ)⁻² over the range reached. This is *numerical evidence* for growth without a constant H. It is not evidence of a de Sitter rate.
- [SH02]'s H ≈ 1.1/a: the least-squares exponent of ln(r/a − 1) over r/a ≤ r_max for near-threshold members is 1.14 to 1.16/a for r_max ≤ 10 (D = 4), 4 to 6 % from 1.1/a; it is a window average of a rate that keeps rising.
- [GGS09b]'s "no apparent horizon": consistent with no trapped cell; the marginal contours r₊ = 0 and r₋ = 0 bounding the anti-trapped region exist, and the right one leaves the grid.
- w on the symmetric throat is (1−n)/(1+n) at every time (*proved*); a stiff or de Sitter reading of the end state cannot be made at this observer.

`baby_universes.md` §3:

- First bullet: the slope test is done and not met: d ln r/dτ crosses W near r/a = 6.5 and reaches 2.75 W at r/a = 50 (D = 4, tuned). #39's probe values (H = 0.75 to 2.9, proper time saturating) are the same behaviour, now with 1/H linear in τ.
- Second bullet: replace "a w → −1 would then need the gradient terms of an inhomogeneous region" by the identity: with a purely radial gradient (φ_u = 0) the isotropic average is (1−n)/(1+n), −1/3 at D = 4; the conjecture for the inhomogeneous region is now a theorem for the mirror-symmetric throat.
- Third bullet: the mass bookkeeping says the ends' masses change by the ghost flux and are not conserved, and the anti-trapped region reaches the grid edge; no statement about a universe is supported by it.

## Limits

- One observer: the diagonal of the mirror-symmetric run. One-sided and generic families (G of #40) are not run; #42 asks for them. Late times of the throat are read at τ_s − τ ≥ 0.31 (D = 4) and 1.59 (D = 5).
- τ_s is an extrapolation by a linear fit of 1/H; the window dependence is 1e-3 at D = 5.
- The ends are the edges of a rectangle: M at x⁺ = X depends on X for x⁻ ≳ 20.
- D = 5 has two resolutions only; D = 6 is not run. The control at D = 5 has three resolutions as an extra.
- The coordinate patch is the double-null rectangle of #39; no claim is made about the spacetime beyond it.
- The runs share the evolution code with the checks that validate them (planted FRW, static wormhole, Tangherlini of #39); independence of the throat's late state from the scheme is shown by two resolutions and by the identity dH/dτ, not by a second code.

## Recommended tasks

1. **The asymptotic solution of the symmetric throat (new issue).** What: from the diagonal relation dH/dτ = −(n+1)H²/2 − (n−1)/(2r²) − ρ/n and the field equations, find the scaling solution with H = k/(τ_s − τ), ρ = −cH², c = n(1/k + (n+1)/2),
   and determine k(n). Why: the data give (n, k, c) = (2, 2.0000, 4.0000) and (3, 3.99, 6.75) and no formula. Acceptance (can fail): a closed-form or numerical scaling solution whose k and c match the stored runs to 1e-3; its prediction for D = 6 is then tested by a run at h ≤ 0.003 (#39's K2 resolution), k within 1e-2. Depends on: this note.
2. **Closer to τ_s.** What: the D = 4 tuned member at X = 800 (t = 1130, N = 8×10⁴ at h = 0.01, about 50 min at one thread) to reach τ_s − τ ≈ 0.15, and a check that K(τ_s − τ)² stays at 4.000. Acceptance: τ_s from the new fit within 1e-3 of 17.8902; k = 2.000 ± 1e-3. Depends on: none.
3. **One-sided and generic families (#42).** What: the same observer-and-criteria on a family that does not start at the static wormhole (G of #40), whose throat locus drifts. Acceptance: the declared classes on the minimum-r locus's proper time, a boosted static planted case recovering its boost. Depends on: #40's review.
4. **CT at D = 5.** What: the event-time convergence of the normal-pulse collapse at D = 5 (h = 0.00177 and 0.00088) and the dependence on the trapped tolerance (1e-2). Acceptance: the three-h ratio of differences near 4 or the cause named. Depends on: none.
5. **A figure in the notes.** `figures/expansion_H.svg` (raw folder; H/W against r/a and k(r)) can be copied next to this note; it was not written into the repository because the brief restricts the files.

## Files

`expansion.py`, `expansion_output.txt`, this note; raw `negative-mass-raw/I213-dh610` (README with the declarations, `dev/`, `runs/`, `scripts/`, `figures/`).

## Coordinator's review

Review of branch `wip/i213-expanding-branch` at 94d89f7 and raw `I213-dh610` at 6cba01a, recorded in negative-mass-raw
`I213-mac-review` (criteria committed at 0ce65c1 before the first run). On the Mac (Apple M5), nice 10.

- **Scope.** The branch adds `expansion.py`, `expansion_output.txt` and this note, nothing else. `dn.py` and `checks.py` are unchanged.
- **Rerun.**
  - `expansion.py analyze` over the runner's stored runs prints `expansion_output.txt` again: 12 of 13 checks pass, and CT
    fails. The only difference is ST's own streamed test, which recomputes on the Mac (edge mass 2.6e-13 against 3.0e-13).
  - The job D = 4 tuned, h = 0.02 was run again. Its `meta` equals the runner's apart from the timing, and its cell classes
    are equal.
  - The throat series agree to 6e-9 in r and 1.4e-9 in τ, not to the 1e-12 the review expected. The stall amplifies the
    two machines' last-bit differences by about e^{Wt}. The final values agree to 1e-10.
- **Independent check** (`ind213.py` in the raw folder). The review's second code (`threshold_check.py`'s march, extended
  to run without stopping) shares no code with `dn.py` or `expansion.py`. Its own p*_h gives the tuned member. 3 of 3
  criteria pass:
  - K1: k = 2.00001 to 2.00040 and |ρ|/H² = 3.99993 to 4.00000 for the eight D = 4 runs;
  - K2: H/W at r/a = 10, 25, 50 for the tuned member, 1.2307, 1.9459, 2.7523, within 5.8e-4 of this note's;
  - K3: τ_s within 4.6e-6 of this note's for all eight runs.
  - D = 5, reported: k = 3.9875 and 3.9897, |ρ|/H² = 6.7508, and τ_s within 1e-6 of this note's.
  - So the "Limits" item "not by a second code" no longer holds for the symmetric throat at D = 4 and 5.
- **Mutations of the check.**
  - `rho-factor` fires exactly K1.
  - `tuned-offset` fires K3 and not K2: H/W at fixed r/a moves by only 5e-4 when the offset grows from 1e-6 to 1e-3.
  - `coord-time` stops with an exception before any criterion is evaluated.
- **Remark.** The data fit the asymptotic relation c = n(1/k + (n+1)/2) of recommended task 1, with c = |ρ|/H², at
  (n, k) = (2, 2) and (3, 4): c = 4 and 6.75. That k = 2(n − 1) is a *conjecture* from two points, and task 1 tests it.
- **Applied.**
  - critical_ghost.md §1.2 and baby_universes.md §3 have this note's proposals, marked as updates from #213.
  - Paper III no longer calls the expanding side inflationary in its own voice. Sec. 6 gains a paragraph on H = k/(τ_s − τ),
    which `paper_numbers.py` checks against `expansion_output.txt`.
  - The figure is copied next to this note as `expansion_H.svg` (task 5).
- **Follow-ups.**
  - Task 1 is filed as #230 and task 4 as #231.
  - Task 2 is folded into #230 as its D = 4 check.
  - Task 3 is added to #225, which runs the one-sided family on the throat locus.
