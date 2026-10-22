# Non-radial (l ≥ 2) linear stability of the signed-source boson stars

Model as in `../boson_star.py` / `../radial`: sextic complex scalar, λ = 10240, μ = λ²/4.8, m = 1,
G_ab = 8πε T_ab (ε = −1 phantom, ε = +1 canonical), field equation unchanged. Background: the BVP
solver `background.py` (a frozen copy of `../radial/background.py`), so α(∞) = 1 and ω is the physical
frequency. Perturbations go as e^{st}; **a growing mode is s = γ > 0** (or Re s > 0).

**Status:** the main results are complete. The extra phantom frequencies ω = 0.70 and 0.85 (l = 2, real γ)
are (pending).

## 1. Equations

### Polar (even), Regge–Wheeler gauge — `derive_nonradial.py`, `reduce_polar.py`
Setup: h_tt = e^ν H0 Y, h_tr = H1 Y, h_rr = e^λ H2 Y, h_AB = r² K Y γ_AB, δΦ = (F + iG) Y e^{−iωt}.
sympy linearises G_ab − 8πεT_ab and the KG equation from the metric, with ε symbolic, and separates the
angles (m = 0, Legendre). The result:
- the traceless angular equation gives **H2 = H0**;
- (tθ), (tr), (rθ) give first-order equations for H1′, K′, H0′;
- (rr), after eliminating H0″ and K″, still contains G″. It is solved **jointly** with KG_R and KG_I for
  {H0, F″, G″}. H0 is then algebraic.

The closed system is 6th order: y = (H1/s, K, F, F′, G/s, G′/s), y′ = J(r; s²) y. J is **polynomial and
even in s**, so the static limit is regular and J is real for real s. It is code-generated in
`polar_matrix_generated.py`; there is no hand transcription.

Bug caught during this work: a first version solved H0 from (rr) alone and the code generator silently
dropped the G″ that it contained. The resulting matrix failed the residual check below (O(1) residuals)
and was deleted.

Checks:
- **All 8 separated equations** (tt, tr, rr, tθ, rθ, angular trace, KG_R, KG_I) are satisfied on solutions
  of the generated system, on real ε = ±1 backgrounds, l = 2, 3, s ∈ {0.05, 0.7, 5}, 4 radii, random y
  (`verify_polar.py`, `verify_polar_output.txt`). The relative residual is ≤ 7e-7, limited by the
  finite-difference J′. A control with y″ perturbed by 1% gives 1e-5…1e-2. The exact rational check in
  `derive_nonradial.py` gives 0 for (tθ), (tr), (rθ).
- **Vacuum limit** (φ0 = 0, Schwarzschild M of either sign), symbolically: the metric and scalar blocks
  decouple. The metric block is exactly the **Zerilli equation** f(fZ′)′ = (s² + V_Z)Z under the Moncrief
  map
  K = [λ(λ+1)r² + 3λMr + 6M²]/[r²(λr + 3M)] Z + f Z′,
  H1/s = [λr² − 3λMr − 3M²]/[(r − 2M)(λr + 3M)] Z + r Z′.
  The scalar block is the free KG equation. Residuals are 0.
- **ε = +1 against a published QNM** (`mini_bs.py`, `mini_bs_output.txt`): Macedo–Pani–Cardoso–Crispino
  2013 "Mini BS I" (ω = 0.853087, our M = 0.633001 vs the paper's 0.63300). Fundamental l = 2 polar mode:
  **σ = 0.119519 − 5.35×10⁻⁵ i** (paper: 0.1195 − 5×10⁻⁵ i). Same code with Λ = μ = 0 and the outgoing
  Zerilli exterior.

### Axial (odd) — `derive_axial.py`, `axial.py`
h_tφ = h0 sinθ ∂θY, h_rφ = h1 sinθ ∂θY, δΦ = 0. With Q = e^{(ν−λ)/2} h1/r:
  d²Q/dr*² = (s² + V_ax) Q,   **V_ax = e^ν [l(l+1)/r² − 6m/r³ + 4πε(ρ − p_r)]**,  ρ − p_r = 2V(φ0²).
The remaining equation holds identically. At ε = +1 this is the Kojima / Macedo et al. boson-star
Regge–Wheeler potential (their (κ/2)(ρ − p_rad) term with κ = 8π).

## 2. Method for growing modes

- Shooting with a real determinant D(γ) at r_m, the radius where φ0 = 10⁻⁹ (r_m ≈ 24).
- **Interior:** the 3-plane of regular solutions (r^l data at r0 = 10⁻³), integrated with continuous QR
  re-orthonormalisation.
- **Exterior:** vacuum. The metric part is the decaying Zerilli solution, Riccati-integrated inward from
  far out; it matches the flat Bessel K solution to O(M), and the Moncrief map takes it to (H1/s, K). The
  scalar part is the 2-plane of decaying KG solutions (spectral projector), integrated inward with QR.
- D = det[Y | Z] × (phases of the QR diagonals). D is real and continuous for real s, and for complex s
  its argument is that of a holomorphic function, which makes argument-principle counts possible.

## 3. Results: real growth rates (`scan_growing.py`; `scan_*.csv`, `scan_*.log`)

100 log-spaced points, with brentq refinement at any sign change:

| star | l | γ range | sign changes of D |
|---|---|---|---|
| phantom ω = 0.7603, M = −0.014580 | 2 | [1e-3, 300] | **0** (D < 0 throughout) |
| phantom ω = 0.7603 | 3 | [1e-3, 1000] | **0** |
| canonical ω = 0.75, M = +0.014626 | 2 | [1e-3, 300] | **0** |
| canonical ω = 0.75 | 3 | [1e-3, 1000] (30 pts) | **0** |

The DG point-mass rate k/(2|M|) is 137 (l = 2) and 686 (l = 3); both lie inside the scanned ranges.

**Convergence** (`convergence.py`, `convergence_l2_output.txt`, phantom l = 2, γ ∈ {1e-3 … 137}):
- rtol 1e-12, background R = 90 and r0 = 3e-3 change D by ≤ 2×10⁻⁴ relative, with the sign unchanged;
- moving r_m to 34 (φ0 = 10⁻¹²) rescales D (it is a normalised determinant), and the sign is unchanged.

The mode condition q_int = q_ext (interior vs exterior Zerilli log-derivative at r_m) is far from
satisfied at every γ: q_int ≈ +√(γ² + …) and q_ext ≈ −√(γ² + …). The interior solution is the
exterior's growing branch.

## 4. Effective surface condition κ_eff (`kappa_eff.py`, `kappa_eff_l2_output.txt`, `love.py`, `love_output.txt`)

Map the regular interior solution, at r_m, to the DG variable Ψ̂ = (∂x − g)Z, with g = ∂x ln ψ0, units
|M| = 1, R = r_m/|M| = 1654. Exact identity (checked with sympy):
  **κ = Ψ̂_x/Ψ̂ = −g + s²/(q − g)**, and **κ_c(R) = −g(R)** (dg_spectrum).
So **κ_eff − κ_c = s²/(q_int − g)**, for any star and any s.

- **s → 0:** κ_eff → κ_c. The intertwiner maps every static Zerilli solution onto φ0, so the static
  limit is degenerate. The sign of the approach is sign(q_int − g), where q_int − g is the decaying
  r^{−l} admixture in the star's static tidal solution, a Love-number-type coefficient c_l. It is
  R-independent to 10⁻⁴…10⁻³ between r_m = 24 and 34, and stable under rtol.

| star | c_2 | c_3 | sign(κ_eff − κ_c) as s → 0 |
|---|---|---|---|
| phantom ω = 0.7603 | −9.27 | −53.9 | **+** |
| phantom ω = 0.85 | −8.30 | −78.9 | **+** |
| canonical ω = 0.75 | +8.64 | +48.4 | − |
| canonical ω = 0.85 | +7.86 | +75.4 | − |

  The phantom stars have κ_eff slightly **above** κ_c (stiffer than critical). The canonical star's sign
  is opposite, but its exterior has M > 0, where κ_c is not an instability threshold (no DG mode exists).
- **All γ > 0 (phantom, l = 2):** κ_eff − κ_c > 0 throughout γ ∈ [1e-4, 150].
  - κ_eff crosses 0 at γ ≈ 2.4×10⁻⁴ (s|M| ≈ 3.5×10⁻⁶) and then grows like +γ|M|.
  - The exterior's own decaying solution has κ_ext ≈ −γ|M| < κ_c.
  - The matter therefore supplies **κ_eff > κ_c (and κ_eff > 0 for γ ≳ 2.5×10⁻⁴)**. By the dg_spectrum
    criterion (unstable iff κ < κ_c for R > r_s), this is the stable side.
  - The large-R threshold κ_c ≈ −(l+1)/R = −1.81×10⁻³ (in 1/|M|) is reproduced: κ_c(R) = −1.817×10⁻³.

## 5. Matter sector

The polar system above already contains the full scalar (matter) sector: F and G couple to the metric,
and a growing matter mode is a zero of the same D. The real-γ scan therefore covers matter modes with
real growth rates as well.

Complex growth rates (oscillatory instabilities, the typical Q-ball-type case) are checked with an
argument-principle count of zeros of D in a rectangle of the right half-plane (`winding.py`).
- D(s̄) = D̄(s), so only the upper half contour is integrated. A zero-free factor e^{sL} that comes from
  normalising at different radii is stripped.
- Validation (`control_winding.py`, `control_winding_output.txt`): around the Mini BS I QNM the winding is
  **+1.0000**; displaced by 0.004i it is **0.0000**.

Production version `winding2.py`:
- the scalar plane is started at an s-independent radius;
- the trivial zero-free phases of the inward/outward exponentials are removed analytically;
- it is restart-safe (per-evaluation cache).

Control (`control_winding2.py`, appended to `control_winding_output.txt`): winding +1.0000 around the
Mini BS I QNM, 0.0000 displaced.

| star | l | box | winding (# growing modes) |
|---|---|---|---|
| phantom ω = 0.7603 | 2 | Re s ∈ [0.1, 2], \|Im s\| ≤ 2.5 | **0.0000** |
| phantom ω = 0.7603 | 3 | same | **0.0000** |
| canonical ω = 0.75 | 2 | same | **0.0000** |

(`winding2_output.txt`; max phase step 0.4 rad; D(x0), D(x1) have the same sign.) No oscillatory
instability with e-folding time ≤ 10 (code units) and |Im s| ≤ 2.5. This range covers the scalar
frequencies ω ± Im s up to the continuum and beyond.

**Not covered:** 0 < Re s < 0.1 off the real axis. That strip passes the scalar continuum thresholds
k± → 0 at Im s = ±(1−ω), which makes the integration slow (R_s ∝ 1/Re k). The real axis down to
γ = 10⁻³ is covered by the scan.

## 6. Conclusions

1. **Polar, real growth rates.** For the phantom campaign star (ω = 0.7603, M = −0.01458) there is
   **no unstable polar mode with γ ∈ [10⁻³, 300] (l = 2) or γ ∈ [10⁻³, 1000] (l = 3)**. The search used
   100 log-spaced points, and D keeps the same sign under all convergence variants. The point-mass DG
   rates (137 and 686) lie inside these ranges. The canonical star (ω = 0.75) gives the same result
   (expected: it is on the stable branch).
2. **Polar, complex growth rates.** Argument-principle count (validated on a known QNM): **0 zeros** in
   Re s ∈ [0.1, 2], |Im s| ≤ 2.5, for phantom l = 2, phantom l = 3 and canonical l = 2. This covers
   matter-sector (Q-ball-type) oscillatory instabilities in that box.
3. **Axial: stable for all l ≥ 2.** V_ax > 0 everywhere; positivity is exact once V_ax > 0 is checked on
   the grid (r²V_ax/(α² l(l+1)) ≥ 1 − 10⁻⁶).
4. **κ_eff.** κ_eff − κ_c = s²/(q_int − g) exactly, and for the phantom star q_int > g, so the matter
   supplies **κ_eff > κ_c at every γ > 0** (l = 2 tabulated for γ ∈ [10⁻⁴, 150]; l = 2, 3 at γ → 0 via
   c_l). The static limit is degenerate (κ_eff → κ_c as s²). The sign is fixed by a Love-type
   coefficient, c_l < 0 for the phantom stars (c_l > 0 canonical). The large-R requirement κ < κ_c ≈
   −(l+1)/R is **not** met. The phantom matter acts as a stiff (κ > 0 for γ ≳ 2.5×10⁻⁴) surface, so the
   DG negative-mass-vacuum instability is not triggered.

## 7. Caveats

- Real-γ scans resolve sign changes on a 100-point log grid (30 points for canonical l = 3). A pair of
  nearby roots (an even number between grid points) would be missed. The κ_eff picture (q_int − q_ext ≈
  2γ, far from 0) makes this implausible, but it is not a proof.
- Complex s: the strip 0 < Re s < 0.1 off the real axis, and |Im s| > 2.5, were not counted.
- Only l = 2, 3. The ω = 0.7603 star was done fully; c_l and axial positivity also at ω = 0.85 (and axial
  at 0.70, 0.90).
- The equations for (tt) and the angular trace are verified numerically (relative residual ≤ 7×10⁻⁷),
  not symbolically. (tθ), (tr), (rθ) are verified exactly.
- Weak-field stars only (|M|/R ~ 3×10⁻³). Nothing here says what happens for compact phantom stars.
- Linear theory, Regge–Wheeler gauge, m = 0 (sufficient by spherical symmetry).

## 8. Files

| file | role |
|---|---|
| `lincommon.py` | shared sympy linearisation of G_ab − 8πεT_ab and KG (ε symbolic) |
| `derive_nonradial.py` (+`_output.txt`) | polar linearisation and separation → `polar_eqs.pkl`; exact check of 3 eqs |
| `reduce_polar.py` (+`_output.txt`) | joint reduction → `polar_matrix_generated.py`; vacuum Zerilli/KG checks |
| `verify_polar.py` (+`_output.txt`) | numerical check of all 8 equations |
| `derive_axial.py`, `axial.py` (+outputs) | axial master equation, positivity of V_ax |
| `background.py` | frozen copy of `../radial/background.py` |
| `polar_solver.py` | determinant D(s), exterior Zerilli/KG, κ_eff helper |
| `scan_growing.py`, `scan_*.csv/log` | real-γ scans |
| `convergence.py` (+`convergence_l2_output.txt`) | robustness of D |
| `winding2.py`, `winding2_output.txt`, `control_winding*.py` | complex-s argument principle + controls (`winding.py` = first version) |
| `kappa_eff.py`, `kappa_eff_l2_*`, `love.py`, `love_output.txt` | κ_eff vs κ_c; static Love-type coefficient |
| `mini_bs.py`, `mini_bs_output.txt` | validation against Macedo et al. 2013 Mini BS I polar QNM |
| `check_bg.py`, `matter_sv.py` | background / denominator sanity checks; scalar-sector SVD diagnostic (not run) |
