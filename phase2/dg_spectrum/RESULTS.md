# Spectrum of the Dotti–Gleiser operator: full problem and truncated "star" problem

Units |M| = 1 (M = −1), f = 1 + 2/r, x = r − 2 ln(1 + r/2) (x = 0 at r = 0).
Ĥ = −∂ₓ² + V̂, V̂ = g² − g′, g = ψ₀′/ψ₀, ψ₀ = χ/(6M + (l+2)(l−1) r), χ from DG eq. (55)
(same code as `../vhat_negativity.py`).

Scripts (run from this directory; outputs in `*_output.txt`, data in `*.csv`, plots `*.png` — PNGs
are git-ignored, so every plot has a CSV):

| file | what |
|---|---|
| `common.py` | exact sympy construction of V̂, regularised coefficients, φ₀ = 1/ψ₀, Chebyshev matrices |
| `full_spectrum.py` | Task 1: full problem on (0, ∞), convergence, Sturm count, l = 2 eigenfunction |
| `truncated.py` | Tasks 2–3: Robin/Dirichlet problem on x ≥ x(R), κ_c(R), grid scan, growth-rate tables |
| `vhat_rneg.py`, `rneg_large_l.py` | r_neg(l) = sup{r: V̂ < 0} (exact root isolation), also r_s, r_* up to l = 20 |
| `kappa_c_large_l.py` | κ_c(R, l) for l ≤ 20 |

## 0. A structural fact that organises everything

Because ψ₀ is a zero-energy solution of the Zerilli operator, V = g² + g′ and Ĥ = (−∂ₓ + g)(∂ₓ + g).
Hence **φ₀ = 1/ψ₀ = ((l+2)(l−1) r − 6)/χ(r) is an exact zero-energy solution of Ĥ** (checked
symbolically: residual 0 for l = 2..8). Properties (checked for l ≤ 20 by exact root counting):
χ > 0 on r > 0, so φ₀ is smooth on r > 0; φ₀ ≃ −4/r ∝ x^(−1/2) at r → 0 (the *non-admissible*
branch); φ₀ ∝ r^(−(l+1)) at ∞ (decaying); exactly one zero, at **r_s(l) = 6/((l+2)(l−1))**
(= −3M/λ, the location of the Zerilli pole); exactly one critical point beyond it, r_*(l).

This gives, via Sturm oscillation theory (regular left endpoint, limit-point at ∞), exact predictions
that the numerics below confirm:

- full problem: #(E < 0) = #zeros of the admissible E = 0 solution;
- truncated, Dirichlet at R: #(E < 0) = #zeros of φ₀ in (R, ∞) = [R < r_s];
- truncated, Robin Ψ̂′ = κΨ̂ at x_R: #(E < 0) = [R < r_s] + [κ < κ_c(R)], with
  **κ_c(R) = (d ln φ₀/dx)(R)** in closed form (κ = κ_c makes φ₀ itself the zero-energy state).

## 1. Full problem on x ∈ (0, ∞)

**Method.** Substitute Ψ̂ = r³ u. Since r ≃ 2x^{1/2}, r³ ≃ 8x^{3/2} is exactly the admissible
(Friedrichs / unique self-adjoint) branch; the other branch x^{−1/2} becomes u ∼ r⁻⁴, which a
polynomial basis cannot represent, so the boundary condition at the singularity is imposed by
construction and u is smooth at r = 0. The u-equation
−r u″ + (−6 + 2/(r+2)) u′ + C(r) u = E r³/(r+2)² u has rational C (built exactly with sympy, regular
at 0). It is discretised by Chebyshev collocation on r = a(1+t)/(1−t), N+1 Gauss–Lobatto points,
r = ∞ node dropped, dense generalised eigensolve. The map scale a sets how resolution is shared
between the core (width ∼ l⁻²) and the far zone (effective domain size); we vary both N and a.
Independent check: Sturm count with DOP853 (rtol 1e-12) for the E = 0 admissible solution.

**Results** (`full_spectrum_output.txt`, `full_convergence.csv`, `full_spectrum.png`):

| l | k | E_exact = −(k/2)² | E₀ (N=240, a=1) | rel. err | #neg (all 20 (N,a) runs) |
|---|---|---|---|---|---|
| 2 | 4 | −4 | −4.00000000000129 | 3e-13 | 1 |
| 3 | 20 | −100 | −100.00000000001231 | 1e-13 | 1 |
| 4 | 60 | −900 | −899.99999999994975 | 6e-14 | 1 |
| 5 | 140 | −4900 | −4900.00000000013915 | 3e-14 | 1 |
| 6 | 280 | −19600 | −19600.00000000148793 | 8e-14 | 1 |
| 7 | 504 | −63504 | −63504.00000000234286 | 4e-14 | 1 |
| 8 | 840 | −176400 | −176399.99999997066334 | 2e-13 | 1 |

- **Count:** in all 140 runs (l = 2..8 × N ∈ {60,100,160,240} × a ∈ {0.1,0.3,1,3,10}) there is
  **exactly one negative eigenvalue**, no complex eigenvalues, and the next eigenvalue is a positive
  discretised-continuum value that → 0⁺ as N or a grows (e.g. l = 2: 1.6e-8 at N=240, a=1;
  1.6e-10 at a=10) — i.e. no second bound state is hiding near 0.
- **Spectral convergence:** error drops from 1e-5…3e-2 (N=60) to ≤ 1e-12 (N≥160–240) for a ∈ [0.1, 1];
  large a under-resolves high-l cores (expected, l=8, a=10 still at 1e-3 at N=240).
- **Sturm check:** the admissible E = 0 solution has exactly one zero (l=2: r=2.783; l=8: r=0.144),
  so exactly one eigenvalue below 0; at E = 0.999E₀ one zero, at E = 1.001E₀ none, so E₀ is the
  ground state. This is a second, ODE-based, discretisation-independent confirmation of
  "exactly one negative eigenvalue" for l = 2..8 (numerical, not a proof for all l).
- **l = 2 eigenfunction:** the numerical eigenvector matches
  Ψ̂ ∝ e^{−2r} r³(r+2)⁴ / (r³ − 3r² + 6) to 1.6e-11 (max relative error at nodes r < 50; N=200).
  I also verified symbolically that Ψ̂ = Â Ψ_DG = (∂ₓ − g)[r(r+2)^k e^{−kr/2}/(2λr − 6)] satisfies
  ĤΨ̂ = −(k/2)²Ψ̂ for l = 2..5.
  **Note on the formula in the brief:** with M = −1 the quoted denominator 4r³ + 12Mr² + 24M³
  = 4r³ − 12r² − 24 has a root at r = 3.49 and does *not* give an eigenfunction (ĤΨ/Ψ is not
  constant). The correct one is r³ + 3Mr² − 6M³ ∝ 4r³ + 12Mr² − 24M³ (sign of the M³ term). Probably
  a transcription slip; I did not re-check the paper's typesetting.

## 2. Truncated problem: star of areal radius R, Robin condition Ψ̂′(x_R) = κΨ̂(x_R)

**Method.** Same Ψ̂ = r³u form (much better conditioned than the bare form for R ≪ 1, where
V̂ ≃ 12/r⁴; the first bare-Ψ̂ attempt lost accuracy at R ≲ 0.01), Chebyshev on
r = R + a(1+t)/(1−t), a = max(R,1), Robin row f u_r + (3f/R − κ)u = 0, Dirichlet by node removal.

**Grid scan** (`truncated_scan.csv`): l = 2..5, 25 log-spaced R ∈ [0.05, 10],
κx_R ∈ {−10, −3, −1, −0.3, 0, 0.3, 1, 3, 10} and Dirichlet — 1000 cases at N = 160 and 240.
**The number of negative eigenvalues agrees with the Sturm prediction [R<r_s] + [κ<κ_c] in all
1000 cases at both resolutions**; lowest-eigenvalue difference N=160 vs 240 ≤ 7e-8.

**κ_c(R)** (`kappa_c.csv`, `kappa_c_check.csv`, `kappa_c_large_l.csv`, left panel of
`truncated.png`): the closed form d ln φ₀/dx agrees with an independent numerical bisection on κ
to ≤ 1e-7 relative at 36 (l, R) points.

| l | R=0.05 | 0.1 | 0.3 | 1 | 3 | 10 | 100 | large-R |
|---|---|---|---|---|---|---|---|---|
| 2 | −846 | −223 | −29.9 | −6.75 | −1.94 | −0.387 | −0.0308 | → −3/R |
| 3 | −885 | −242 | −40.3 | −8.36 | −2.18 | −0.469 | −0.0407 | → −4/R |
| 4 | −934 | −269 | −223 | −13.4 | −2.47 | −0.571 | −0.0507 | → −5/R |
| 5 | −998 | −319 | +19.0 | −13.8 | −2.84 | −0.678 | −0.0608 | → −6/R |
| 8 | −1507 | +1490 | −124 | −17.8 | −4.07 | −1.00 | −0.0911 | → −9/R |
| 20 | −4140 | −1217 | −214 | −37.9 | −9.18 | −2.31 | −0.212 | → −21/R |

Reading the table (κ_c is the boundary: unstable iff κ < κ_c, plus one extra unconditional mode
when R < r_s):

| l | r_s (Dirichlet unstable iff R < r_s) | r_* (Neumann unstable iff R < r_*) | r_neg (V̂ ≥ 0 beyond) |
|---|---|---|---|
| 2 | 1.5000 | 2.2788 | 2.6741 |
| 3 | 0.6000 | 0.8789 | 1.0144 |
| 4 | 0.3333 | 0.4814 | 0.5520 |
| 5 | 0.2143 | 0.3073 | 0.3512 |
| 8 | 0.0857 | 0.1219 | 0.1388 |
| 20 | 0.01435 | 0.02032 | 0.02309 |

(r_neg l² ≈ 8.8–10.7 for l = 2..20, i.e. r_neg ≈ 9/l².)

- **Large R (R > r_*(l)):** κ_c < 0, so Dirichlet, Neumann and every κ ≥ 0 are stable; an unstable
  mode needs κ < κ_c(R) ≈ −(l+1)/R for R ≫ 1 (the static exterior solution φ₀ ∝ r^{−(l+1)}).
  In x-units, κ_c x_R ≈ −2.5 (l=2) to −4.3 (l=5) at R = 10, drifting to −(l+1) as R → ∞.
  So the matter must supply a boundary condition "softer" than the vacuum static tail by an O(1)
  factor on the scale 1/R.
- **Consistency with the analytic bound:** 282 cases with R > r_neg(l) and κ ≥ 0 (or Dirichlet):
  0 negative eigenvalues. Since r_* < r_neg, the true threshold for κ ≥ 0 stability is R > r_*(l)
  (slightly sharper than the r_neg bound); for r_* < R < r_neg, V̂ is still negative somewhere
  but not enough to bind.
- **Small R (R ≪ r_s):** κ_c x_R → −1/2 (because φ₀ ∼ x^{−1/2}). For any κ x_R > −1/2 — including
  Dirichlet and Neumann — there is exactly one negative eigenvalue and it converges to the DG value:

| l | R | E_D | E_N | E_DG | rel D | rel N |
|---|---|---|---|---|---|---|
| 2 | 0.001 | −4.0000000000 | −4.0000000000 | −4 | 9e-13 | 6e-13 |
| 2 | 0.1 | −3.99997 | −4.00009 | −4 | 8e-6 | 2e-5 |
| 2 | 1.0 | −3.432 | −5.014 | −4 | 0.14 | 0.25 |
| 5 | 0.001 | −4899.9999995 | −4900.0000014 | −4900 | 1e-10 | 3e-10 |
| 5 | 0.05 | −4895.79 | −4910.75 | −4900 | 9e-4 | 2e-3 |
| 5 | 0.2 | −1509 | −9427 | −4900 | 0.69 | 0.92 |

  Error scales roughly like R^{3}–R^{4} (`small_R.csv`): the x^{3/2} vs x^{−1/2} branch selection
  is insensitive to the boundary condition once the cut is deep inside the mode. This confirms the
  DG self-adjoint extension is the R → 0 limit of a regular cut-off with any "non-pathological" κ.
  For κ x_R < −1/2 a second negative eigenvalue appears (deeply bound, set by κ), consistent with
  the count formula.

## 3. Growth rates √(−E₀) vs l (`growth_rates.csv`, right panel of `truncated.png`)

| l | full DG k/2 | R=0.05 Dir | R=0.05 Neu | R=0.3 Dir | R=0.3 Neu | R=1 Dir | R=1 Neu | R=1 κ=−10 | R=3 κ=−10 | R=10 κ=−1 | R=10 κ=−10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2 | 2.000 | 2.000 | 1.999 | 2.002 | 1.853 | 2.239 | 8.68 | 9.79 | 0.905 | 9.990 |
| 3 | 10 | 10.000 | 10.000 | 9.846 | 10.32 | 0 | 0 | 5.72 | 9.69 | 0.865 | 9.986 |
| 4 | 30 | 29.998 | 30.006 | 19.63 | 40.28 | 0 | 0 | 0 | 9.62 | 0.802 | 9.980 |
| 5 | 70 | 69.97 | 70.08 | 0 | 43.64 | 0 | 0 | 0 | 9.52 | 0.717 | 9.973 |
| 6 | 140 | 139.7 | 140.7 | 0 | 0 | 0 | 0 | 0 | 9.38 | 0.602 | 9.964 |
| 7 | 252 | 250.0 | 256.3 | 0 | 0 | 0 | 0 | 0 | 9.23 | 0.436 | 9.955 |
| 8 | 420 | 408.3 | 442.5 | 0 | 0 | 0 | 0 | 0 | 9.05 | 0 | 9.944 |

- **Point singularity:** growth = k/2 = (l−1)l(l+1)(l+2)/12 ∼ l⁴/12, unbounded in l (confirmed
  numerically through l = 8, exactly one mode per l).
- **Any finite R, Dirichlet or Neumann (or any κ ≥ 0):** only l with R < r_*(l) ≈ 1.3·r_s(l) are
  unstable — finitely many l (roughly l ≲ 3/√R) — so the growth rate is capped, at about the
  DG value of the largest unstable l (Neumann can overshoot k/2 near threshold, e.g. l=4 at R=0.3:
  40 vs 30).
- **Any finite R, fixed κ < 0:** for every l with R > r_neg(l) (all l ≳ 3/√R) V̂ ≥ 0 on the
  exterior and the quadratic form gives E ≥ −κ², i.e. growth ≤ |κ|. Checked on 188 scan cases:
  max growth/|κ| = 0.986. So the l → ∞ growth is capped by |κ| (the R=3, R=10 columns tend to |κ|).
  Overall sup_l growth ≤ max(|κ|, growth of the finitely many l with r_neg(l) > R) < ∞.
- Hence: the unbounded-in-l growth (∼ l⁴) is a property of the point singularity; a regular cut-off
  at any R > 0 with an l-independent κ caps it. **Caveat:** if the physical interface condition made
  κ depend on l (e.g. κ ∼ −l^p/R), the cap would scale accordingly; this computation does not
  determine κ, which must come from the interior matter model.

## Caveats

- "Exactly one negative eigenvalue" for the full problem is numerical (spectral eigensolver + ODE
  Sturm count), l = 2..8 (the φ₀ facts behind the truncated count are exact for l ≤ 20). Not a proof
  for all l.
- The Robin condition is a stand-in for the matter; real stars give an energy-dependent (and
  l-dependent) interface condition, not a constant κ. All statements about stars are conditional on
  that model.
- r_s(l) coincides with the location of the Zerilli pole −3M/λ; a truncation exactly at R = r_s was
  excluded from scans (κ_c has a pole there).
- Eigenvalues are of Ĥ; growth rate in t is √(−E) (units 1/|M|).
