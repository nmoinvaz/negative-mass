# Review: Shirokov GRTeclyn — evolution code (C++) and phantom sector

Repo: https://github.com/Nikchik-coder/GRTeclyn (branch `develop`, commit 09fbcd6). Paper: arXiv:2608.24577.
Read-only; nothing from the repo was executed. Independent check scripts in this directory:
`emtensor_check.py`, `boost_sign_check.py`, `sponge_ko_attenuation.py`, `drift_compare.py`, `momentum_vs_velocity.py`.

All packed campaign cells: matter `grtresna_bicomplex_scalar`, CCZ4, pump off (`rl_pump_stop_time = 0`),
reservoir mode 0, sponge r = 24–32 strength 4, `max_level = 0` (except amrcheck), treadmill off (except chase03c),
pair lumps painted at rest.

**Bottom line: no bug found in the evolution physics behind the runaway.**

## Findings (ranked)

### F1 — bug, does not affect paper results (verified; re-checked by coordinator)
Boosted-lump painter stores Π with the wrong sign. `grteclyn-wrapper/src/grteclyn_wrapper/grtresna/fields/boson_star.py:237-238`:
`Pi_re = -dt_phi_re * inv_lapse; Pi_im = -dt_phi_im * inv_lapse`, while the evolution uses
`rhs[c_phi_x] = vars.lapse() * phi_dot_Pi + advec_phi` (`Source/Matter/BiComplexScalarField.impl.hpp:132`), i.e.
Π = +∂ₜφ/α (zero shift). The evolution therefore receives time-reversed data: a star painted with +v moves at −v
(1D test: v = +0.3 → centroid −12.0 by t = 40 vs +12.0 with the correct sign). At rest the error only maps ω → −ω
(a charge-conjugated, still stationary star). All pair cells are painted at rest; nonzero velocities appear only in
single-star cells (`num_lumps = 1`). Any future boosted-star run would be affected.

### F2 — bug, negligible for the paper's wave-zone claims (verified; re-checked by coordinator)
Matter term of Ψ4 always zero. `Source/Matter/Weyl4WithMatter.impl.hpp:123-125`:
`amrex::Real G_Newton = 0; pp.queryAdd("G_newton", G_Newton);` — the key is set nowhere (code defaults or campaign
params), so the −4πG S^TF contribution to E_ij is multiplied by 0. The repo's own `Tests/Weyl4WithMatterTest`
defaults it to 1.0. At the R = 8 extraction sphere (~3 units from the stars) Ψ4 is wrong; at R ≥ 16 the matter is
~1e-6 of central and the effect is negligible.

### F3 — risk (verified by stencil analysis)
The "absorbing sponge" (`Source/Grids/SpongeZone.hpp:87-93`) is only extra 6th-order Kreiss–Oliger dissipation,
ramped 4((r−24)/8)^4; damping rate (σ/dx) sin⁶(k dx/2). Amplitude kept after one radial pass:

| wavelength | 128³ (dx=0.5) | 256³ (dx=0.25) |
|---|---|---|
| 8 | 0.9993 | 0.99998 |
| 4 | 0.96 | 0.9986 |
| 2 | 0.20 | 0.92 |

Resolved-mode damping scales as dx⁵ (~30× weaker at 256³), so the convergence ladder also changes the boundary
treatment. "Ejecta drain through the sponge" is not supported. Acts identically on both sectors (mirror cell
antisymmetric to 4–5 digits), so it cannot create the signed drift, but it is an unquantified, resolution-dependent
systematic. Minor: ignores `prob_lo` (0 here). `boundary_flux.dat` not checked.

### F4 — interpretation risk (verified from packed data)
Gap estimators disagree on sign at equal mass (d = 10):

| grid | barycentre | core | proper |
|---|---|---|---|
| 256³ | 10 → 9.915 | 10 → 10.013 | 10 → 10.024 |
| 192³ | 10 → 9.930 | 10 → 10.029 | 10 → 10.040 |
| 128³ | 10 → 10.003 | 10 → 10.107 | 10 → 10.006 |

"Flat at equal mass" uses the barycentre (unsigned field activity over the whole box incl. halo and sponge). The
sign flip across mass ratios 0.597/0.799/1.333 holds for every estimator. Headline acceleration is estimator-independent
(1.5491e-4 barycentre vs 1.5527e-4 core at 256³).

### F5 — minor risk (suspected, not quantified)
The repaint sets Π_im with the solved global lapse (includes the partner's potential), while GRTresna's source used
each star's own tabulated α(r) (`patches/grtresna-matter.diff`). ~1e-3 difference at d = 10, opposite sign at the two
stars — expected to cause slight breathing rather than a push. Gridinit files not in the pack.

### F6 — style / minor
Shift sampled at the integer peak cell (`sector_dynamics.py`); √γ = χ^(−3/2) assumes det h = 1 (second-order
error); collapse diagnostics at `RadialRecipeLevel.cpp:1014-1018` cover only the canonical real part (unused by the
paper); NEC monitor samples 24 directions × 3 speeds (not rigorous) and acts on the net T where sectors overlap.

## Checked and correct
- ρ, S_i, S_ij, trS (`BiComplexScalarField.impl.hpp:13-86`) match the 3+1 projection of T = Σ ε_s T[Φ_s] to ~1e-18
  on a random metric with χ, α ≠ 1, β ≠ 0 and the sextic potential; ½ kinetic normalisation consistent; trS not
  double-counted; each sector contributes exactly ε_s × the single-field expression.
- CCZ4 matter coupling (G = 1): K += 4πα(S − 3ρ); Θ += −8παρ; A_ij += −8παχ S^TF; Γ^i and B^i += −16πα h^ij S_j
  (B correct since its RHS is built from Γ's before matter); `chi_rhs_floor = 0`.
- Klein–Gordon identical and canonical for both sectors; ε never enters; potentials depend only on own |Φ|² (no
  cross-coupling); pump inactive. dV/dφ_a = 2φ_a dV/d|Φ|² matches the initial-data ODE and a finite difference;
  couplings consistent across evolution params, solver params and star tables.
- Variable layout slots 25–28 Φ+, 29–32 Φ− consistent across enum, painter, gridinit map and extractors; parities
  even, asymptotics 0; matter RHS zeroed; dissipation and sponge on all variables; t = 0 peaks match tables to 0.25%.
- Momentum: P_x[s] = ∫ −(Π₁∂ₓφ₁ + Π₂∂ₓφ₂)√γ, same convention both sectors; signed total = px_c − px_p (matter, not
  ADM, momentum); px/(|M|v) 1.03–1.17 at 256³, same for both; px_c/px_p 1.071 (128³) → 1.007 (256³).
- No sector asymmetry other than ε in RHS, constraints, energy conditions, trackers, momentum; mirror reproduces to
  4–5 digits. Constraint and energy-condition monitors include μ and use G = 1.

## Not reviewed
CCZ4 vacuum RHS, gauge and boundary internals (assumed upstream); GRTresna beyond the packed patch; treadmill
(RecentringBox), RLPumpForce (inactive), AMR constraint norms; Python Ψ4 mode decomposition; `article_figures.py`
beyond `midpoint`/`fit_a`; `boundary_flux.dat`; gridinit files (not in the pack).
