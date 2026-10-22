# Phase 2 results (in progress)

## 1. Where the vacuum exterior can drive an instability — `vhat_negativity.py`

After Dotti–Gleiser's intertwining transformation, the even-parity problem on negative-mass
Schwarzschild is a regular Schrödinger operator Ĥ = −∂ₓ² + V̂. A negative eigenvalue needs V̂ < 0
somewhere. Using their closed-form seed (χ polynomial, eq. 55; static residual verified 0 for l ≤ 5;
χ > 0 on r > 0 for l ≤ 8), V̂ is negative only on a thin shell near the singularity:

| l | k | V̂ < 0 on r ∈ (units \|M\|) | min V̂ | \|M\|/R bound |
|---|---|---|---|---|
| 2 | 4 | [1.545, 2.674] | −19.8 | 0.374 |
| 3 | 20 | [0.618, 1.014] | −500 | 0.986 |
| 4 | 60 | [0.343, 0.552] | −4.5e3 | 1.81 |
| 5 | 140 | [0.221, 0.351] | −2.5e4 | 2.85 |
| 6 | 280 | [0.154, 0.244] | −9.8e4 | 4.09 |
| 7 | 504 | [0.114, 0.180] | −3.2e5 | 5.54 |
| 8 | 840 | [0.088, 0.139] | −8.8e5 | 7.20 |

**Statement.** If a negative-mass star has areal surface radius R with |M|/R < 0.374 (checked l ≤ 8;
the shell moves inward with l), then V̂ ≥ 0 everywhere outside the star. For any interface condition
of Robin type Ψ̂′ = κΨ̂ with κ ≥ 0 the exterior quadratic form is non-negative, so the vacuum exterior
cannot host an unstable mode by itself: any even-parity instability must be supplied by the matter
through the interface. (The odd-parity sector is already stable, GHI.) Caveat: κ is set by the matter;
the sign of κ is not established here.

## 2. Signed-source boson stars — `boson_star.py`

Shirokov's model (arXiv:2608.24577): sextic complex scalar, λ = 10240, λ²/μ = 4.8, units m = 1,
G_ab = 8πε T_ab with ε = −1 for the phantom (field equation unchanged ⇒ positive inertial/passive mass,
negative active/ADM mass). Static shooting in Schwarzschild-like coordinates, ω_phys = ω/α(∞).

| star | ω_phys | M_ADM (ours) | M_ADM (paper) | rms radius (\|φ\|-weighted) | paper | \|M\|/R |
|---|---|---|---|---|---|---|
| canonical | 0.7500 | +0.014626 | +0.014350 | 4.344 | 4.34 | 3.4e-3 |
| phantom | 0.7603 | −0.014580 | −0.014295 | 4.402 | — | 3.3e-3 |

- Radius and compactness reproduce the paper. Our masses are 1.9–2.0 % above the paper's.
  **Resolved (see `mass_check/RESULTS.md`):** the paper's solver reports the mass enclosed at the radius
  where φ falls to 5 % of φ_c (r ≈ 5.2, still inside the matter tail), not the ADM mass at infinity. Our
  enclosed mass at that radius reproduces the paper to 3×10⁻⁴; profiles agree to 1.4×10⁻⁵. Our values are
  the true ADM masses. The bias depends on frequency (0.977–0.985), so it is not a fixed factor.
- At equal frequency the phantom is heavier in magnitude (paper's Fig. 2): canonical ≈ +0.01424 vs
  phantom −0.01482 at ω ≈ 0.757. Both families lighten toward ω → 0.9 (≈ 0.57 of the campaign mass at
  ω ≈ 0.895; paper reports a floor ≈ 0.54 near ω ≈ 0.93).

| ω_phys (canonical) | M | ω_phys (phantom) | M |
|---|---|---|---|
| 0.6916 | +0.020506 | 0.7081 | −0.019763 |
| 0.7426 | +0.015218 | 0.7573 | −0.014820 |
| 0.7934 | +0.011823 | 0.8065 | −0.011620 |
| 0.8441 | +0.009624 | 0.8558 | −0.009551 |
| 0.8949 | +0.008339 | 0.9050 | −0.008387 |

**Consequence for §1.** These stars have |M|/R ≈ 3×10⁻³, two orders below 0.374. In code units the
l = 2 negative shell would sit at r ≈ 0.02–0.04, deep inside a core of rms radius ≈ 4.4. Whatever
the stability of these stars is, it is a property of the phantom matter, not of the negative-mass vacuum.

## 3. Dotti–Gleiser spectrum reproduced; exact stability thresholds — `dg_spectrum/RESULTS.md`

- Full problem, l = 2..8: E₀ = −(k/2)² to ≤3×10⁻¹³ (Chebyshev, Ψ̂ = r³u); exactly **one** negative eigenvalue in
  all 140 resolution/domain runs, confirmed by Sturm counts. (Numerical evidence for l ≤ 8, not a proof for all l.)
- φ₀ = 1/ψ₀ is an exact zero mode of Ĥ = (−∂ₓ+g)(∂ₓ+g), with a single node at r_s = 6/((l+2)(l−1)). This gives an
  exact count for the truncated (star-like) problem with Robin condition Ψ̂′ = κΨ̂ at areal radius R:
  #negative = [R < r_s] + [κ < κ_c(R)], κ_c = d ln φ₀/dx at R (closed form; matches bisection to 1e-7; 1000-case
  scan agrees).
- Thresholds (units |M|): Dirichlet unstable iff R < r_s (l=2: 1.5); Neumann/any κ ≥ 0 unstable iff R < r_*
  (l=2: 2.279, l=3: 0.879, l=4: 0.481, l=5: 0.307) — sharper than the §1 bound r_neg (2.674, 1.014, 0.552, 0.351).
  So the κ ≥ 0 statement of §1 holds with |M|/R < 1/2.279 = 0.439.
- Large R: an instability needs κ < κ_c → −(l+1)/R, i.e. a surface condition softer than the static vacuum tail by an
  O(1) factor. Small R: truncation reproduces the DG eigenvalue with error ∝ R⁴.
- Growth rate: the ~l⁴/12 growth belongs to the point singularity; any finite R with l-independent κ caps it
  (κ ≥ 0: only l ≲ 3/√R unstable; fixed κ < 0: growth ≤ |κ|). Caveat: a real star's interface condition may depend
  on l and frequency — not determined here.
- Note: DG08 eq. (27) as rendered on ar5iv has a sign slip in the l = 2 eigenfunction denominator; the correct form is
  r³(r−2M)⁴ e^{2(t−r)/|M|}/(4r³ + 12Mr² − 24M³).

## 4. Radial (l = 0) stability — `radial/RESULTS.md`

Derived with sympy (all six linearized equations residual 0), solved by compound-matrix shooting, cross-checked
by an independent FD pencil including complex eigenvalues.
- **Canonical validation:** the fundamental mode crosses σ² = 0 at the extremum of M and N (ω = 0.925844) to 2×10⁻⁸,
  and again at the strong-field maximum of N (ω ≈ 0.1996) within scan resolution.
- **Phantom (ε = −1):** the campaign star ω = 0.7603 is **radially stable** (no real growth rate up to |σ| = 2; no
  complex eigenvalue). The family is stable for ω ≲ 0.92891 and unstable above, exactly at dN/dω = 0 (to 8×10⁻⁸).
  The turning-point criterion carries over to ε = −1 **stated in N** (first law becomes dM = εω dN, so the M-form
  of the criterion flips sign). The phantom ω-fold at ω ≈ 0.4223 is not a stability change; the strong-field
  phantom branch (|2m/r| up to 0.4) is stable with no turning point found.
- Caveat: near ω ≈ 0.93 the instability is essentially the flat-space Vakhitov–Kolokolov Q-ball instability;
  the strong-field branch is the more meaningful test of how ε enters.

## Next
Non-radial (l ≥ 2) perturbations of the phantom star, with the effective surface condition κ_eff compared to κ_c.
