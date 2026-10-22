# Why our boson-star ADM masses are 2% above Shirokov's

**Short answer:** the paper does not report the ADM mass of its stars. Its solver reads the
mass at a finite radius r_m, the first grid point where φ₀ < 5 % of φ_c. That radius is about
5.2, still inside the star's matter tail, so the result is the mass enclosed there. It is not
the total mass at infinity. Our M_ADM = +0.014626 / −0.014580 is correct. The paper's
+0.014350 / −0.014295 is low by the matter outside r ≈ 5.2, about 1.9–2.0 %. The same
convention applied to our solution gives the paper's numbers to ≤ 3×10⁻⁴ (relative). The
couplings, normalisation and frequency definitions all agree.

## Source of the convention (paper repo, read as data only, nothing executed)

Repo: `Nikchik-coder/GRTeclyn@develop`, cloned into scratch.

- `grteclyn-wrapper/src/grteclyn_wrapper/grtresna/profiles/boson_star_ode.py`, `_build_profile`:
  - match index: `descended = (phi < 0.05*phi_c) & (dphi < 0)`, take the first such point;
  - `psi_inf = psi + r psi'`, `alpha_inf = alpha + r alpha'` at r_m (an A + B/r fit there);
  - `omega_phys = omega_int/alpha_inf`, `r_phys = psi_inf² r`;
  - **`adm_raw = 2 r_m,phys (psi(r_m)/psi_inf − 1)`**, and that value is what gets stored as `adm_mass`;
  - outside r_m it attaches an analytic Yukawa φ tail and an isotropic-Schwarzschild ψ, α built
    from this M. The table's α column beyond the seam encodes M = 0.01434996 (canonical), which
    I confirmed by inverting the tail.
- The model matches ours exactly. U = ½φ² − ¼λφ⁴ + ⅙μφ⁶ with ½-normalised kinetic terms,
  ρ = ½ω²φ²/α² + ½ψ⁻⁴φ′² + U, ρ + S = 2ω²φ²/α² − 2U. The phantom flips only the metric sources.
  λ = 10240, μ = **21845333.0** (literal in `analysis/star_family_scan.py` and in the table
  headers).
- The paper's `stars/star_family.csv`, `star_radius.csv` and both `stars/*.dat` headers all
  carry this finite-radius `adm_mass`.

## Numerical demonstration

`iso_star.py` is an independent re-implementation of the paper's isotropic, maximally sliced
ODEs, written from their stated equations. For each solution it extracts the mass two ways:

- **paper-style:** the paper's rule, applied at the radius where φ = 5 % φ_c, interpolated;
- **true:** M = 2AB from ψ = A + B/r far outside the matter, at r ≈ 25–32 where φ ≈ 10⁻⁹–10⁻¹¹.

| star | φ_c (ours, iso) | φ_c (paper table) | M paper-style (ours) | M paper | M true (iso) | M (our Schw. code) |
|---|---|---|---|---|---|---|
| canonical ω=0.750 | 0.0197204282 | 0.0197204696 | +0.0143460 | +0.0143500 | +0.0146265 | +0.014626 |
| phantom ω=0.7603 | 0.0196642262 | 0.0196641883 | −0.0142937 | −0.0142947 | −0.0145803 | −0.014580 |

- **Profiles agree.** On the paper's own isotropic grid, out to the seam, our φ₀(r) matches
  theirs to max |Δφ|/φ_c = 1.4×10⁻⁵ (canonical) and 1.3×10⁻⁵ (phantom). The lapse matches to
  3×10⁻⁷, and φ_c to 2×10⁻⁶. They are the same star.
- **Paper-style extraction reproduces the paper.** The relative difference is −3×10⁻⁴
  (canonical) and −7×10⁻⁵ (phantom). The remaining difference comes from their unrefined
  grid-point seam: their seam sits at φ/φ_c = 0.0496–0.0498 rather than exactly 0.05.
- **The true mass is reproduced in two different gauges.** The isotropic far-field 2AB and our
  Schwarzschild-gauge m(∞) agree to 5 digits.
- **Reverse check.** Our Schwarzschild code's Misner–Sharp mass m(R) at the areal radius where
  φ = 5 % φ_c is +0.014345 (R = 5.192) and −0.014294 (R = 5.222). The paper/ours ratios are
  0.99967 and 0.99991 (`quick_m_of_r.py`).
- **The bias depends on frequency, so it is not one constant factor.** Two more points from
  `star_family.csv` (`family_check.py`):

  | star | M paper | M paper-style (ours) | paper / true |
  |---|---|---|---|
  | canonical ω=0.80 | +0.011209 | +0.011204 (−4.5×10⁻⁴) | 0.9766 |
  | phantom ω=0.70 | −0.020521 | −0.020515 (−2.7×10⁻⁴) | 0.9854 |

  The common 0.981 at the two campaign stars only happens because they have nearly the same
  ω and profile.

## Candidates ruled out

| candidate | test | effect |
|---|---|---|
| μ = 21845333 vs λ²/4.8 = 21845333.33 | `mu_check.py` | ΔM < 10⁻⁸ (relative 10⁻⁸ in μ) |
| ω_phys definition (α_inf from A + B/r at r_m vs at ∞) | `iso_star.py`, both targets | Δω ≈ 7×10⁻⁶, ΔM ≈ 5×10⁻⁷ (3×10⁻⁵ relative) |
| kinetic normalisation / potential | read from code, and profile match to 10⁻⁵ | identical |
| isotropic vs Schwarzschild gauge | true mass computed in both | agree to 5 digits |

Komar mass and time units were not tested separately. They are unnecessary: the finite-radius
read already closes the gap to < 0.05 %.

## Consequences

- The evolved spacetimes carry the **true** masses (≈ +0.014626 / −0.014580). GRTresna
  re-solves the constraints using the full painted profile, and the extraction convention only
  changes the reported label. The paper's quoted star masses, and any mass-scaled quantities,
  are therefore low by 1.9–2.3 %.
- Equal-mass matching: the paper quotes |M₋|/M₊ = 0.014295/0.014350 = 0.9962. With true masses
  it is 0.014580/0.014626 = 0.9968. The mismatch is slightly smaller than quoted, and the
  conclusion is unchanged.
- Our phase-2 §2 numbers stand as correct ADM masses. Shirokov's numbers should be read as
  "mass inside r(φ = 5 % φ_c) ≈ 5.2".

## What was verified and what was inferred

- **Verified numerically:** profile agreement; the paper-style extraction reproducing the
  paper's masses to ≤ 4.5×10⁻⁴; the true mass in two gauges; that μ rounding and the ω
  definition do not matter.
- **Read from source and not executed:** the extraction rule itself (`_build_profile`).
- **Inferred, not checked against the paper PDF or the campaign mass diagnostics:** that the
  paper quotes these `adm_mass` values as "ADM mass", and that the evolved spacetime carries the
  true mass.

## Files

- `iso_star.py` — independent isotropic solver, both mass extractions (`iso_star_output.txt`)
- `compare_tables.py` — profile and table comparison (`compare_tables_output.txt`). Run with the
  paper's `stars/` directory as the argument.
- `quick_m_of_r.py` — Misner–Sharp m(R at 5 %) from our Schwarzschild code
- `mu_check.py`, `family_check.py` — the μ rounding test and the extra family points
  (`*_output.txt`)
