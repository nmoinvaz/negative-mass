# Radial (l = 0) stability of the signed-source sextic boson stars

Model as in `../boson_star.py` (Shirokov arXiv:2608.24577): λ = 10240, μ = λ²/4.8, m = 1, G_ab = 8πε T_ab,
ε = +1 canonical, ε = −1 phantom. Perturbations ∝ e^{iσt}; **σ² < 0 is unstable**.
Throughout, **ω is the physical frequency** (α(∞) = 1), and N is the Noether particle number
N = ∫ 4π r² ω φ² (a/α) dr, which does not depend on ε.

## 1. Equations (`derive_radial.py`, output in `derive_radial_output.txt`)

Metric e^{ν}, e^{λ} perturbed by n(t,r), l(t,r). Field Φ = [φ₀ + F + iG] e^{−iωt}, and everything ∝ e^{iσt}
(F → A, G → B, n → N_, l → L). sympy builds G_ab − 8πεT_ab and the KG equation from the metric, linearizes
them, and reduces as follows:

- background: reproduces the `boson_star.py` equations exactly (residual 0). The θθ equation follows from Bianchi.
- (tr) constraint: L = 8πε [r φ₀' A − ω (a/α) Ψ / r]
- Ψ (the perturbation of the enclosed charge) is defined by iσΨ = r²(α/a)(φ₀B' − φ₀'B). With the
  gauge invariant G̃ = iσB + ωφ₀N_/2, the imaginary KG equation (current conservation) integrates to
  Ψ' = r²(a/α) φ₀ (G̃ − ωφ₀L/2 − 2ωA)
- (rr) constraint: N_' is algebraic in (A, A', G̃, Ψ), with no undifferentiated N_ (gauge).
- the Ψ definition becomes −σ²Ψ = r²(α/a)(φ₀G̃' − φ₀'G̃ − ωφ₀²N_'/2)
- the real KG equation gives A'' = (…) − σ² (a/α)² A + …, linear in σ².

The result is a closed 4th-order system for (A, A', G̃, χ = Ψ/φ₀) that contains σ only as σ², with ε kept
symbolic. **Check:** with the reduced system substituted back, the residuals of all six linearized equations
(tt, tr, rr, θθ, KG_R, KG_I) simplify to exactly 0. The full expressions are printed in
`derive_radial_output.txt`. `radial_rhs_generated.py` is generated from them, so no hand transcription is
involved. The structure (one Noether-charge function, one gauge-invariant phase, metric removed by the two
constraints) matches Gleiser–Watkins/Jetzer/Kain. I did not compare against their printed equations term by
term. The symbolic residual check is the validation of the derivation.

## 2. Numerical method

**Background (`background.py`).** The shooting stars of `boson_star.py` are trustworthy only to r ≈ 23. I
re-solve each star as a two-point BVP (`scipy.solve_bvp`, tol 1e-9, R = max(60, 30/√(1−ω²))) in the variables
ψ = ln φ, u = φ'/φ, m, ν:
- regular series at r₀ = 1e-3;
- decaying Yukawa tail u(R) = −√(1−ω²) − 1/R;
- Schwarzschild matching for ν(R), which normalizes α(∞) = 1.

The log variable keeps the exponentially small tail well conditioned (φ(R)/φ_c ~ 1e-20). Families are
continued in ω, in N with ω free (this passes ω-folds), or in α(0) with ω free (this passes N-extrema). At
the campaign stars the BVP reproduces `boson_star.py`: M = +0.014626 at ω = 0.75 and M = −0.014580 at
ω = 0.7603.

**Eigenproblem (`radial_modes.py`).** I used shooting with the compound-matrix (bivector) method:
- the 2-plane of regular solutions is integrated from r₀ (A₀, G̃₀ free; χ = O(r³));
- the 2-plane of decaying solutions is integrated inward from R. It starts on the spectral projector of the
  far-field matrix in the variables (A, A', G̃, χ/r²), whose eigenvalues reproduce ±√(1−(ω±σ)²) (checked);
- the two planes are matched at r_m (where φ = 10⁻³ φ_c);
- D(σ²) is the normalized wedge product, and eigenvalues are its zeros.

Why this method: plain 2-column shooting failed. The two decaying solutions have different rates, so their
span collapses numerically and D jumped by orders of magnitude. The bivector method keeps the plane
orthonormal. Bound modes exist only for σ² < (1−ω)², the threshold of the continuum.

**Independent check (`fd_check.py`, `complex_check.py`).** This is a second-order box-scheme finite-difference
discretization of the same first-order system on a stretched grid. It gives a generalized pencil
(K₀ + σ²K₁)y = 0, which I solve by dense QZ (cond K₀ ~ 1e22 defeats sparse shift-invert). It shares only the
derived equations with the shooting code. It also returns **complex** σ², which shooting along the real
σ²-axis cannot see. One spurious exact σ² = 0 eigenvector appears; it is localized in the outer tail, where
φ₀ underflows, and is discarded.

## 3. Validation on the canonical family (ε = +1)

Scan (`scan_canonical.py`, `scan_eps+1.json`):

| ω | M | N | lowest σ² |
|---|---|---|---|
| 0.750 | +0.0146260 | 0.0159284 | none in [−4, (1−ω)²) → stable |
| 0.800 | +0.0114777 | 0.0118555 | none (stable) |
| 0.850 | +0.0094304 | 0.0093676 | none (stable) |
| 0.880 | +0.0086225 | 0.0084324 | none (stable) |
| 0.900 | +0.0082617 | 0.0080267 | +9.462e-3 (just below threshold 1.0e-2) |
| 0.910 | +0.0081437 | 0.0078962 | +6.205e-3 |
| 0.920 | +0.0080761 | 0.0078223 | +2.376e-3 |
| 0.925 | +0.0080647 | 0.0078099 | +3.448e-4 |
| 0.930 | +0.0080709 | 0.0078166 | **−1.685e-3** |
| 0.935 | +0.0080974 | 0.0078450 | −3.638e-3 |
| 0.940 | +0.0081478 | 0.0078988 | −5.431e-3 |
| 0.950 | +0.0083413 | 0.0081034 | −8.177e-3 |
| 0.960 | +0.0087168 | 0.0084964 | −9.236e-3 |

Fine scan with 17 stars on ω ∈ [0.923, 0.931] (`crossing.py`, `crossing_p.txt`), using polynomial fits:

| | ω |
|---|---|
| dM/dω = 0 | 0.92584446 |
| dN/dω = 0 | 0.92584446 |
| σ² = 0 | **0.92584444** |

**The zero crossing of the fundamental mode coincides with the extremum of M and N to 2×10⁻⁸ in ω.** The
extremum is a minimum of M. The stable side is ω < 0.9258, where dM/dω < 0 and dN/dω < 0, as expected.
The first law dM = ω dN holds to ~5×10⁻⁵ by finite differences, which independently checks the background
and the M/N extraction.

Convergence (`convergence.py`, `convergence_output.txt`). The eigenvalue changes by ≤ 2×10⁻¹⁰ (relative
≤ 3×10⁻⁷) under each of the following:
- integrator rtol 1e-10 → 1e-12;
- r_m × 0.7 or × 1.4;
- eigen-solver R − 15;
- background tol 1e-9 → 1e-10;
- background R + 30.

For example, ω = 0.94 gives σ² = −5.4314738e-3 in every variant. The FD check gives σ² = −5.449e-3 at
n = 300, consistent with O(h²) (more resolutions are in `fd_check_output.txt`, see §6). FD finds no complex
eigenvalue with |σ²| < 0.3 at ω = 0.75 (n = 300 and 600) or at ω = 0.94 (n = 300).

Canonical low-ω / strong-field part (`lowscan.py`, `strongfield.py`; `low_p.txt`, `strong_p.txt`).
- Along the canonical family **N and M keep growing down to ω ≈ 0.200, where N peaks at ≈ 42.73
  (M ≈ 11.93, max 2m/r ≈ 0.66)**. The value "ω_min = 0.316" corresponds to N ≈ 9.2, not to a turning point.
- The fundamental σ² stays positive and drops toward 0 as that maximum is approached:

| α(0) | ω | N | M | σ²₀ |
|---|---|---|---|---|
| 0.9 | 0.3812 | 2.017 | 0.872 | 2.25e-2 |
| 0.8 | 0.3286 | 7.048 | 2.632 | 7.57e-3 |
| 0.7 | 0.2949 | 13.95 | 4.776 | 3.48e-3 |
| 0.6 | 0.2671 | 21.97 | 7.025 | 1.73e-3 |
| 0.5 | 0.2426 | 30.30 | 9.146 | 8.49e-4 |
| 0.4 | 0.2210 | 37.79 | 10.88 | 3.59e-4 |
| 0.35 | 0.2116 | 40.638 | 11.4994 | 1.96e-4 |
| 0.30 | 0.20353 | 42.4079 | 11.86739 | 6.59e-5 |
| 0.28 | 0.200834 | 42.69442 | 11.925379 | +2.00e-5 |
| 0.27 | 0.199631 | 42.72676 | 11.931862 | **−2.01e-6** |
| 0.26 | 0.198538 | 42.67740 | 11.922043 | **−2.36e-5** (past the maximum, N falling) |

**The second turning point also reproduces the criterion.** σ₀² changes sign between α₀ = 0.28 and 0.27,
and the maximum of N (and M) lies in the same interval: N = 42.6944 → 42.72676 → 42.6774, peak within
≈ 0.27 ± 0.005. Beyond it the star is unstable. The resolution is set by the α₀ step (0.01); I did not
refine the coincidence here to the 1e-8 level achieved at the high-ω turning point.
(The value "ω_min = 0.316" corresponds to N ≈ 9.2, not to a turning point; the family continues to ω ≈ 0.199.)

## 4. Phantom family (ε = −1)

High-ω scan (`scan_eps-1.json`):

| ω | M | N | lowest σ² |
|---|---|---|---|
| 0.7500 | −0.0154193 | 0.0168547 | none (stable) |
| **0.7603 (campaign)** | **−0.0145805** | **0.0157438** | **none in [−4, 0.0575) → stable** |
| 0.8000 | −0.0119681 | 0.0123897 | none (stable) |
| 0.8500 | −0.0097452 | 0.0096881 | none (stable) |
| 0.9000 | −0.0084615 | 0.0082156 | none (stable) |
| 0.9200 | −0.0082351 | 0.0079664 | +3.597e-3 |
| 0.9250 | −0.0082129 | 0.0079424 | +1.603e-3 |
| 0.9275 | −0.0082081 | 0.0079372 | +5.801e-4 |
| 0.9300 | −0.0082078 | 0.0079369 | **−4.501e-4** |
| 0.9350 | −0.0082222 | 0.0079523 | −2.494e-3 |
| 0.9400 | −0.0082593 | 0.0079919 | −4.446e-3 |
| 0.9500 | −0.0084205 | 0.0081623 | −7.705e-3 |
| 0.9600 | −0.0087488 | 0.0085059 | −9.475e-3 |

Fine scan (`crossing_m.txt`):

| | ω |
|---|---|
| dM/dω = 0 | 0.92890838 |
| dN/dω = 0 | 0.92890838 |
| σ² = 0 | **0.92890846** |

The crossing coincides with the extremum to 8×10⁻⁸. The phantom first law is **dM = εω dN = −ω dN**, which
holds to ~8×10⁻⁵. The extremum is a maximum of M, i.e. a minimum of |M|. Because dM = εω dN, M and N always
have their extrema at the same ω, for either sign of ε.

For the campaign star there is no sign change of D on σ² ∈ [−4, (1−ω)²) at two tolerances
(`deep_negative_output.txt`, `convergence_output.txt`). That excludes real growth rates up to |σ| = 2 in code
units (e-folding time ≥ 0.5) and beyond. The FD complex-spectrum check for ε = −1 finds no complex eigenvalue; see §6.

Low-ω / strong-field phantom (`lowscan.py` via N-continuation, `strongfield.py` via α(0); `strong_m.txt`):
- The phantom family has a **fold in ω (ω_min ≈ 0.42226 at N ≈ 5.6)**. N and |M| are monotone through it,
  and ω turns back up.
- The fundamental σ² stays positive and smooth across the fold. No mode crosses there.

| α(0) / N | ω | N | M | σ²₀ |
|---|---|---|---|---|
| N = 0.3 | 0.4947 | 0.30 | −0.171 | 1.02e-1 |
| N = 1 | 0.4477 | 1.00 | −0.496 | 5.23e-2 |
| α₀ = 1.10 | 0.4313 | 2.04 | −0.950 | 3.59e-2 |
| α₀ = 1.15 | 0.4228 | 4.34 | −1.932 | 2.49e-2 |
| α₀ = 1.17 (fold) | 0.42226 | 5.46 | −2.404 | 2.25e-2 |
| α₀ = 1.20 | 0.4230 | 7.33 | −3.193 | 1.99e-2 |
| α₀ = 1.25 | 0.4269 | 10.92 | −4.718 | 1.72e-2 |
| α₀ = 1.30 | 0.4328 | 15.07 | −6.502 | 1.55e-2 |
| α₀ = 1.40 | 0.4479 | 24.91 | −10.84 | 1.36e-2 |
| α₀ = 1.50 | 0.4653 | 36.65 | −16.20 | 1.27e-2 |

Over the whole range reached (max |2m/r| ≈ 0.4 at α₀ = 1.5), no extremum of N or M was found on the
phantom low-ω side, and no unstable radial mode.

## 5. Conclusions

1. **Canonical:** reproduced the textbook result at both turning points of the family. The fundamental
   mode crosses σ² = 0 at the weak-field extremum of M and N (ω = 0.925844), to 2×10⁻⁸. It crosses again
   at the strong-field maximum of N and M (ω ≈ 0.1996, N ≈ 42.727), to the scan resolution. The stable
   side is dN/dω < 0 between them.
2. **Phantom:** the campaign star ω = 0.7603 is **radially (l = 0) linearly stable**.
   - The family is stable for ω ≲ 0.9289 and unstable above it. The instability sets in exactly at the
     common extremum of M and N (minimum of |M|).
   - **For this family the turning-point criterion carries over to ε = −1, if it is stated in N:**
     instability onset at dN/dω = 0, stable side dN/dω < 0. Stated in M, the sign of the stable-side
     inequality flips, since dM/dω = −ω dN/dω > 0 there. The "dM/dω < 0" rule is not ε-invariant; the
     N-rule is.
   - Plausible reason, not proven here: at dN/dω = 0, ∂_ω(background) is a σ = 0 solution with δN = 0. That
     argument is kinematic and does not use the sign of the energy functional. Which side is stable is not
     given by that argument. Here it was found numerically.
   - The ω-fold of the phantom family (ω ≈ 0.4223) is **not** a stability change, as expected: neither M nor
     N is extremal there.
3. **Physical caveat on the high-ω result.** At ω ≈ 0.93 the stars have |M|/R ~ 10⁻³. The instability there
   is essentially the flat-space Vakhitov–Kolokolov (Q-ball) instability, dN/dω > 0, and gravity is a
   perturbation. That is why the two signs give almost the same numbers: σ² crosses at 0.92584 vs 0.92891.
   It is a weak test of how ε enters. The strong-field phantom branch (|2m/r| up to ~0.4) is the more
   meaningful test, and it is stable with no turning point within the range reached.

## 6. Caveats / open items

- **Coarse root grids on the low-ω tables.** `lowest_modes` brackets roots on a fixed σ² grid. A re-check of
  canonical N = 3 with a 5× denser grid found the fundamental at **1.637e-2**. The coarse run listed
  1.749e-1 there: it skipped a close pair. Rows from `low_p.txt` / the 16-point grid may therefore have
  missed the lowest positive root. The `strongfield.py` rows use a 64-point grid on [−0.2, 0.1], and the
  high-ω tables are bracketed densely near 0, which is where the sign question lives. **Negative σ² can
  still be missed only if an even number of roots sits between grid points.** The stronger statement for
  the campaign stars, "no root on [−4, 0)", comes from 40 log-spaced points plus the FD spectrum.
- **Strong-field canonical maximum (ω ≈ 0.2)** is resolved only to the α₀ step 0.01 (§3). Not refined.
- **FD and complex-spectrum checks (completed before a container restart cut off the remaining cells):**
  - FD lowest σ² converge under n = 250/400/550 and agree with shooting: ε=+1 ω=0.94 −5.437e-3 (shooting
    −5.431e-3); ε=−1 ω=0.9275 +5.816e-4 (shooting +5.801e-4); ε=−1 ω=0.7603 no bound mode (lowest FD value
    0.05844 sits above the continuum threshold 0.05746).
  - Complex spectrum (|σ²| < 0.3): none for ε=+1 ω=0.75, 0.94 and ε=−1 ω=0.7603, 0.85 at n = 300 and 600;
    ε=−1 ω=0.94 at n = 300 (real −4.456e-3, matching shooting −4.446e-3). Remaining planned cells were lost.
- The FD pencil is extremely ill-conditioned because coefficients carry φ₀ ~ e^{−kr}. Eigenvalues are
  stable under refinement where checked, but the dense solve limits n to ≈ 1200.
- Masses here are true ADM masses; the +2 % difference from the paper is the paper's 5 %-radius mass
  convention (see ../mass_check/RESULTS.md). Does not affect the stability statements.
- Only l = 0. The phantom star could still be unstable to l ≥ 1 or to non-linear perturbations.
- Not checked: the full (unfolded) phantom family beyond α₀ = 1.5, and whether it eventually develops an
  N-extremum.

## Files
| file | role |
|---|---|
| `derive_radial.py` → `derive_radial_output.txt`, `radial_eqs.txt`, `radial_rhs_generated.py` | sympy derivation + verification + codegen (~45 s) |
| `background.py` | BVP stars; continuation in ω, N, α(0) |
| `radial_modes.py` | compound-matrix shooting eigen-solver |
| `fd_check.py`, `complex_check.py` | independent FD discretization, complex spectrum |
| `scan_canonical.py` (both ε) → `scan_eps±1.json`, `fine_*.txt` | high-ω scans |
| `crossing.py` → `crossing_eps±1.json`, `crossing_*.txt` | σ² = 0 vs dM/dω = 0 vs dN/dω = 0 |
| `convergence.py` → `convergence_output.txt` | tolerance/radius convergence |
| `deep_negative.py` → `deep_negative_output.txt` | no real growth rate up to \|σ\| = 2 at the campaign stars |
| `fold_scan.py` → `fold_p.txt` | N-continuation of the families (folds) |
| `lowscan.py`, `strongfield.py` → `low_p.txt`, `strong_*.txt`, `low_recheck_*.txt` | low-ω / strong-field spectra |
