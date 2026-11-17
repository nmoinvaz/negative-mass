# Non-radial (l ≥ 2) linear stability of the signed-source boson stars

Model as in `../boson_star.py` / `../radial`: sextic complex scalar, λ = 10240, μ = λ²/4.8, m = 1,
G_ab = 8πε T_ab (ε = −1 phantom, ε = +1 canonical), field equation unchanged. Background: the BVP
solver `background.py` (a frozen copy of `../radial/background.py`), so α(∞) = 1 and ω is the physical
frequency. Perturbations go as e^{st}; **a growing mode is s = γ > 0** (or Re s > 0).

**Status:** complete for the scope below (see Caveats).

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
| phantom ω = 0.70, M = −0.020824 | 2 | [1e-3, 300] | **0** |
| phantom ω = 0.70 | 3 | [1e-3, 1000] | **0** |
| phantom ω = 0.85, M = −0.009745 | 2 | [1e-3, 300] | **0** |
| phantom ω = 0.85 | 3 | [1e-3, 2000] | **0** |

The DG point-mass rate k/(2|M|) is 137 (l = 2) and 686 (l = 3) for the campaign star; both lie inside the scanned
ranges. It scales as 1/|M|, so the ω = 0.85, l = 3 scan goes to 2000.

The ω = 0.70 and 0.85 rows are from #7 (`TAG=_n100`, raw negative-mass-raw I7-win2 at 60553bc). The ω = 0.85, l = 2
scan crashed twice on win2 (access violations, at a different offset each time) and was completed on the Mac (I7-mac
at d766289). Its CSV agrees with win2's serial evaluation of D on the same grid to a relative 7×10⁻⁹ at worst
(`cmp_scan_diag.py`, pass criterion 10⁻⁸). The earlier 30-point l = 2 scans of these stars, `scan_eps-1_w0.7_l2.csv` and
`scan_eps-1_w0.85_l2.csv`, also had no sign change.

**Convergence** (`convergence.py`, `convergence_l2_output.txt`, phantom l = 2, γ ∈ {1e-3 … 137}):
- rtol 1e-12, background R = 90 and r0 = 3e-3 change D by ≤ 2×10⁻⁴ relative, with the sign unchanged;
- moving r_m to 34 (φ0 = 10⁻¹²) rescales D (it is a normalised determinant), and the sign is unchanged.

#7 repeats this for ω = 0.70 and 0.85, l = 2 and 3 (`convergence.py l w`). D is negative in all five variants at every
γ, with a relative spread at the same r_m of at most 1.4×10⁻⁴ (ω = 0.70) and 2.6×10⁻⁴ (ω = 0.85) (I7-win2 at 60553bc).
The ω = 0.85 jobs crashed on win2 before finishing and were completed on the Mac (I7-mac at d766289), with every row
both machines printed identical to the printed digits.

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
| phantom ω = 0.70 | 2 | same | **0.0000** |
| phantom ω = 0.70 | 3 | same | **0.0000** |
| phantom ω = 0.85 | 2 | same | **0.0000** |
| phantom ω = 0.85 | 3 | same | **−2.0000** |

(`winding2_output.txt` for the first three rows, #7's I7-win2 at 60553bc for the rest; max phase step 0.4 rad; D(x0),
D(x1) have the same sign.)

**Discrepancy (#7).** A count of zeros cannot be negative, so for ω = 0.85, l = 3 D has a net pole in the box. A
hypothesis, tested in #54: the Zerilli column (1, q/f) of `polar_solver.decaying_plane` has a pole wherever the
decaying Zerilli solution vanishes at r_m. If so, each winding counts zeros minus poles, and a 0 alone does not exclude a
zero. Until #54 settles it, the 0 counts are conditional on D having no pole in the box.

No oscillatory
instability with e-folding time ≤ 10 (code units) and |Im s| ≤ 2.5. This range covers the scalar
frequencies ω ± Im s up to the continuum and beyond.

**Not covered in the table above:** 0 < Re s < 0.1 off the real axis. That strip passes the scalar continuum
thresholds k± → 0 at Im s = ±(1−ω) and ±(1+ω), which makes the integration slow (R_s ∝ 1/Re k). The real axis down to
γ = 10⁻³ is covered by the scan. For the campaign star the strip, the band above |Im s| = 2.5 and the real axis are
counted in §5.1 (#5).

### 5.1 The strip, the band above |Im s| = 2.5 and the real axis (#5)

Counted for the phantom campaign star (ω = 0.7603, M = −0.014580), l = 2 and 3, with the D of §2 on rpi5
(`winding2.py` at rtol 10⁻⁸, `refine_gamma.py` at 10⁻¹⁰). Every number below is in `i5_boxes_output.txt`, which
`collect_i5.py` builds from the raw logs (negative-mass-raw, folder I5-rpi5, at RAWCOMMIT).

#### A defect in the counter

The planted zero of the first strip box was not returned. The pair planted at 0.04 ± 0.25i in Re s ∈ [0.03, 0.1],
|Im s| ≤ 2.5 gave a count of 0 instead of 2. On the left edge, one step from Im s = 0.25 to 0 had a true phase change
of 2π − 0.37 for D times the planted factor, because D alone turns by 4.36 and the factor by 1.55. The 0.4 rad rule
saw only the measured angle, −0.37, and accepted the step, losing 2π. Accepting a step only when its two half steps add
up to it did not help, since both half steps wrapped (measured −2.40 and +2.03, sum −0.37). `winding2.py` now accepts a
step of the first two bisection levels only when both half steps are at most 0.4 rad themselves and add up to it, and
`control_winding2.py` does the same. Every count below uses the new rule. The stored box (Re s ∈ [0.1, 2], |Im s| ≤ 2.5,
l = 2) still gives winding +0.0000, with the D(x0), D(x1), max phase step 0.39 and min|D| 6.70×10⁻⁴ of
`winding2_output.txt` (rerun from scratch with the stored script at 478cbd9, then with the new rule).

One such loss shifts the count of a symmetric box by ±2, the size of #7's −2 for ω = 0.85, l = 3. The windings of §5
were computed with the first rule and were not recounted here. Whether the −2 is this or a pole is for #54.

#### Boxes

```
Im s
 10 ┌──────────┬──────────┬──────────────────────┐
    │          │          │ h2 (l = 2)           │
  6 │ not      ├──────────┼──────────────────────┤
    │ counted  │ s3       │ h1 (l = 2, 3)        │
    │          │ (l = 2,3)│                      │
2.5 ├──────────┼──────────┼──────────────────────┤
    │ s2       │ s1       │ stored box           │
    │ (l = 2)  │ (l = 2,3)│ (§5, l = 2, 3)       │
0.23├──────────┴──────────┤                      │
    │ n1 (l = 2, 3)       │                      │
  0 └─────────────────────┴──────────────────────┘
    0.01      0.03       0.1                      2   Re s
```

Widths are not to scale. s1 and s2 also contain |Im s| ≤ 0.23, where n1 overlaps them, and n1 reaches Re s = 10⁻⁴
there. n3 crosses the axis (Re s ∈ [−10⁻³, 10⁻⁴], Im s ∈ [0.01, 0.23]). The threshold boxes t1 and t2 lie inside s1 and
s2. Each row is one contour, a symmetric box (|Im s| ≤ y1, doubled from the upper half) or an upper box
(Im s ∈ [y0, y1], its mirror has the same count). s3 was added after the first plan, since the strip above |Im s| = 2.5
belongs to both gaps of the issue and to none of the first boxes.

BOXTABLE

Winding is zeros minus poles of D in the box (§5, #54), so each 0 means no zero only if D has no pole there.

#### Controls

- **Planted zero.** Each row's last columns give the zero planted in the same cached evaluations, as the factor
  (s − s0)(s − s̄0) in a symmetric box and (s − s0) in an upper box, and what the count returned. The stored box has four
  mutations of the counter, each of which must make the check fail and does: reverse the contour (−2), apply the wrong
  doubling (+1), move the factor's zero out of the box (0) and plant s0 alone in the symmetric box (+1.72).
- **Known QNM, strip geometry.** On the Mini BS I star (l = 2) boxes Re s ∈ [−2×10⁻⁴, 0.1] below the threshold
  Im s = 0.1469 count the zeros of D that `polish_zero.py --mini` finds by a secant on the same D at rtol 10⁻⁸ and 10⁻¹⁰:
  −5.35350×10⁻⁵ + 0.11951908i (the published mode, spread of Re s over rtol 2×10⁻¹⁵) and three more,
  −2.62464×10⁻⁵ + 0.131605i, −1.35809×10⁻⁵ + 0.137308i and −7.7455×10⁻⁶ + 0.140372i. They are the zeros the winding
  returned: Im s ∈ [0.09, 0.125] +1, [0.09, 0.14] +3, [0.1225, 0.14] +2, [0.1225, 0.1285] 0. The same counts come at a
  doubled outer radius (R_s 419, 717 and 454 against 227, 376 and 245). The `reverse` and `shift` mutations make it fail.
  I first expected 1 and 0 in the second and third boxes, since the published mode was the only one I knew, and
  changed the expectation after the run returned +3 and +2. Only the first box and the published value are independent
  of the added zeros, the others compare the counter with a root finder that shares D.
- **Two outer radii.** Near each scalar threshold (t1 around Im s = 1 − ω = 0.2397, t2 around 1 + ω = 1.7603) the count
  is made at R_s = r_m + 40/min Re k and at 80/min Re k with the Zerilli start radius doubled. The counts and the
  min|D| agree (rows t1 and t2).
- **Axis scan, not used.** `axis_scan.py` looked for a dip of |D| on Re s = 0.005, Im s ∈ [0.005, 0.235]. Its control on
  the Mini BS fails (the known zero shows no dip, |D| varies smoothly across it), so the absence of a dip for the
  campaign star says nothing and is not used.

#### Real axis

`refine_gamma.py` starts from the stored 100-point scan and adds the log-midpoints of each level, with D from the same
`PolarProblem.D`. Three stored points recomputed agree to a relative 3×10⁻⁹ at worst.

REFTABLE

The planted pair (γ = 0.6035 and 0.6408 for l = 2, 1.115 and 1.184 for l = 3, 6 % apart) lies inside one level-0
interval, is missed there, and is found in the two intervals that contain it at level 1 and level 2. With `--mutate=coarse`
it is missed, with `--mutate=pair` an unplanted pair makes the count change between levels, and with `--mutate=stale`
the recomputed values differ from the stored ones by 1 %, each failing its check as designed. The ratio of |D_i| to
the geometric mean of its neighbours stays within 1.000 to 1.003, so ln |D| is smooth at every grid point. A pair of
roots closer than about 3 % in γ is therefore excluded only by that smoothness, not by a sign change.

#### Not covered

- NOTCOVERED

Status: numerical evidence, conditional on D having no pole in each box (#54).

## 6. Conclusions

1. **Polar, real growth rates.** For the phantom campaign star (ω = 0.7603, M = −0.01458) there is
   **no unstable polar mode with γ ∈ [10⁻³, 300] (l = 2) or γ ∈ [10⁻³, 1000] (l = 3)**. The search used
   100 log-spaced points, and D keeps the same sign under all convergence variants. The point-mass DG
   rates (137 and 686) lie inside these ranges. The canonical star (ω = 0.75) gives the same result
   (expected: it is on the stable branch). The phantom stars ω = 0.70 and 0.85 give the same result for
   l = 2 and 3 (#7).
2. **Polar, complex growth rates.** Argument-principle count (validated on a known QNM): winding **0** in
   Re s ∈ [0.1, 2], |Im s| ≤ 2.5, for phantom l = 2, phantom l = 3 and canonical l = 2, and for ω = 0.70
   (l = 2, 3) and ω = 0.85 (l = 2). This means 0 zeros only if D has no pole in the box. ω = 0.85, l = 3
   gives −2, a net pole, so that condition is open (§5, #54).
3. **Axial: stable for all l ≥ 2.** V_ax > 0 everywhere; positivity is exact once V_ax > 0 is checked on
   the grid (r²V_ax/(α² l(l+1)) ≥ 1 − 10⁻⁶).
4. **κ_eff.** κ_eff − κ_c = s²/(q_int − g) exactly, and for the phantom star q_int > g, so the matter
   supplies **κ_eff > κ_c at every γ > 0** (l = 2 tabulated for γ ∈ [10⁻⁴, 150]; l = 2, 3 at γ → 0 via
   c_l). The static limit is degenerate (κ_eff → κ_c as s²). The sign is fixed by a Love-type
   coefficient, c_l < 0 for the phantom stars (c_l > 0 canonical). The large-R requirement κ < κ_c ≈
   −(l+1)/R is **not** met. The phantom matter acts as a stiff (κ > 0 for γ ≳ 2.5×10⁻⁴) surface, so the
   DG negative-mass-vacuum instability is not triggered.

## 7. Caveats

- The real-γ scans of the other stars resolve sign changes on a 100-point log grid (30 points for canonical l = 3). A
  pair of nearby roots (an even number between grid points) would be missed there. The κ_eff picture (q_int − q_ext ≈
  2γ, far from 0) makes this implausible, but it is not a proof. For the campaign star the grid was refined to 397
  points (step ratio 1.032 for l = 2, 1.036 for l = 3) with no change in the sign count, so a pair closer than about 3 %
  in γ is excluded there only by the smoothness of ln |D| (§5.1).
- Complex s, campaign star: counted are Re s ∈ [0.1, 2] for |Im s| ≤ 6 (HMAX for l = 2), the strip down to Re s = 0.03
  (SMIN for l = 2) for |Im s| ≤ 2.5, and Re s ∈ [−10⁻³, 0.1] for |Im s| ≤ 0.23 (§5.1). Not counted: Re s below that
  edge for 0.24 < |Im s| ≤ 2.5, |Im s| above that bound, Re s > 2, and |Im s| < 0.01 left of Re s = 10⁻⁴. The other stars
  have Re s ≥ 0.1, |Im s| ≤ 2.5 only.
- The windings of §5 were computed with a step rule that can lose 2π in one step (§5.1) and were not recounted. One
  loss shifts a symmetric-box count by ±2.
- Only l = 2, 3. The ω = 0.7603, 0.70 and 0.85 stars have the real-γ scans, convergence variants and windings (0.70
  and 0.85 in #7). c_l and axial positivity are also checked at ω = 0.85, and axial positivity at 0.70 and 0.90.
- The windings count zeros minus poles of D if D has poles in the box, and ω = 0.85, l = 3 shows a net pole (§5, #54).
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
| `winding2.py`, `winding2_output.txt`, `control_winding*.py` | complex-s argument principle + controls (`winding.py` = first version). `winding2.py` takes `--y0`, `--rsmult`, `--zfac`, `--nseg`, `--cachedir`, `--plant` and `--mutate` (§5.1) |
| `refine_gamma.py` | real-γ refinement by halving the scan step, with a planted pair and mutations (§5.1) |
| `polish_zero.py` | secant polish of a zero of D near the axis, with the Mini BS I control |
| `axis_scan.py` | \|D\| on a vertical line; its Mini BS control fails, not used (§5.1) |
| `collect_i5.py`, `i5_boxes_output.txt` | the table of §5.1 from the raw logs (negative-mass-raw, I5-rpi5) |
| `kappa_eff.py`, `kappa_eff_l2_*`, `love.py`, `love_output.txt` | κ_eff vs κ_c; static Love-type coefficient |
| `mini_bs.py`, `mini_bs_output.txt` | validation against Macedo et al. 2013 Mini BS I polar QNM |
| `check_bg.py`, `matter_sv.py` | background / denominator sanity checks; scalar-sector SVD diagnostic (not run) |
