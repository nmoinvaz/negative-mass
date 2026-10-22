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

## Next
Linear perturbations of the phantom star: radial (l = 0) first, then even-parity l ≥ 2 with the
interface/interior treated explicitly.
