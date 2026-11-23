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
  finite-difference J′. A control with y″ perturbed by 1% gives 1e-5…1e-2. The exact rational check of (tθ), (tr), (rθ)
  that this text attributed to `derive_nonradial.py` is not in that script (the run was interrupted; #8).
- **The same 8 equations, exactly** (`verify_polar_exact.py`, `verify_polar_exact_output.txt`, #8). On y′ = J y with y″
  = (J′ + J²) y and H0 = h·y, J′, h′, h″ by the chain rule through the background (written as `lincommon.bgsub` writes
  it), every residual is **exactly 0**, for all six components of y. *Checked symbolically*: as a polynomial identity in
  r, α, a, φ0, u, ω, ε, s, l(l+1), Λ, μ and π as twelve indeterminates (ε generic, not only ±1), by exact gcd-free
  arithmetic (4 s; the generated J has one denominator, an 8-term polynomial printed in the output). It is also 0 at 20
  random rational points with ε = ±1 in sympy's rational function field in π, and, in that field, with only r, or r and
  s, or r, α and a free. Controls that fire: J[0, 0] scaled by 1 + 1/7, the sign of each of α′, a′, φ0′, u′ reversed,
  and the coefficients of H0″ in (tt) and of K″ in the angular trace scaled by 1 + 1/7 each leave a nonzero residual in
  (tt) and the angular trace, symbolically and at 20 of 20 points. The first attempt, in the rational function field
  with a gcd after every operation, did not get past J′ in 5 min 48 s (negative-mass-raw `I8-rpi5/probe`). What it
  checks: the reduced system implies the separated equations of `polar_eqs.pkl`. It does not re-derive them from G_ab −
  8πεT_ab. *Coordinator's review (#8's closing comment):* an independent derivation, with its own linearization of G_ab
  − 8πεT_ab and its own background equations, shares only `J_and_h` with this check. On y′ = J y it gives exactly 0 for
  (tt) and the angular trace at 6 random rational points, with l = 2 and 3, ε = ±1 and π replaced by a rational. Its
  controls pass: the separation, and the background (θθ) equation. Its mutations (a J coefficient scaled by 8/7, ε
  flipped in T_ab, J′ dropped) each fire. That check re-derives the separated equations, which this one takes as given;
  it is an identity test at points, not a symbolic proof (negative-mass-raw `I8-mac`).
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

## 3a. Polar modes at l = 4, 5, 6 (#4; `scan_growing.py`, `convergence.py`, `kappa_eff.py`)

Phantom campaign star (ω = 0.7603, M = −0.014580, r_m = 24.11), the same code and settings as §3, run in worker pools
(`NPROC`) with resumable partial files. Raw logs, partial files and launchers: negative-mass-raw `I4-win` at b156730.
Code commits on `wip/i4-polar-l456`: scans 05e6394 (l = 4), a1f8eb1 (l = 5), 215936b (l = 6); convergence 73a9e0b,
a1f54d1, 7682c8e; κ_eff 73a9e0b (l = 4), f2ecad0 (l = 5, 6). Every status below is *numerical evidence*.

**Control.** `mini_bs.py` rerun on this machine reproduces σ = 0.119519 − 5.35×10⁻⁵ i. Its output differs from
`mini_bs_output.txt` in 14 lines: the four pre-polish shooting rows differ in ω at about 2×10⁻⁶ and the last
|D| is 1.35×10⁻¹⁵ against 1.45×10⁻¹⁶, with the secant iterates and the target star equal to the printed digits
(`mini_bs_vs_stored.diff` in the raw folder). The coordinator took this as platform noise; the control passes.

**Real-γ scans** (100 log-spaced points on [10⁻³, γ_max], `scan_eps-1_w0.7603_l{4,5,6}.csv`):

| l | k/(2\|M\|) | γ range | sign changes of D | D at γ = 10⁻³ | D at γ_max |
|---|---|---|---|---|---|
| 4 | 2058 | [1e-3, 3000] | **0** (D < 0 at all 100 points) | −9.993e-2 | −4.595e-12 |
| 5 | 4801 | [1e-3, 6000] | **0** | −1.261e-1 | −8.936e-13 |
| 6 | 9602 | [1e-3, 12000] | **0** | −1.542e-1 | −1.596e-13 |

*Numerical evidence*: the scans reach 1.46, 1.25 and 1.25 times the point-mass rate k/(2|M|), as the l = 2, 3 scans
do (§3), and have no sign change. The largest step of ln|D| between neighbouring points is 0.45, 0.47 and 0.49, so
no jump is hidden on the grid (the 100-point caveat of §7 still applies).

`l456_graph.py` draws the scans, the convergence base values and the point-mass rates from these files.

![−D(γ) at l = 4, 5, 6 with the convergence base values and the point-mass rates](l456_mode_condition.svg)

**Convergence** (`convergence_l{4,5,6}_output.txt`): the five variants of §3 (base, rtol 1e-12, φ_cut 1e-12 at
r_m = 34.21, background R = 90, r0 = 3×10⁻³) and the q task at γ ∈ {1e-3, 0.03, 0.3, 1, 3, 10, 40, 137, 300, 1000,
γ_max}.

| l | evaluations done | sign of D | max spread, same r_m | D(r_m = 34)/D_base | rtol 1e-12 at γ_max |
|---|---|---|---|---|---|
| 4 | 66 of 66 | < 0 in all five variants at all 11 γ | 2.1e-4 | 0.67–0.68 at γ ≤ 0.03, 0.50 at γ ≥ 1 | done, 17168 s |
| 5 | 65 of 66 | < 0 in every variant that finished | 2.1e-4 | 0.66–0.67 at γ ≤ 0.03, 0.50 at γ ≥ 1 | **not completed** |
| 6 | 65 of 66 | < 0 in every variant that finished | 2.1e-4 | 0.66 at γ ≤ 0.03, 0.50–0.51 at γ ≥ 1 | **not completed** |

- *Numerical evidence*: D keeps its sign under the same variants as at l = 2, 3, and the spread over the variants at
  the same r_m is at most 2.1×10⁻⁴ (l = 2: 2×10⁻⁴), at each of the 11 γ.
- **Capped rows.** The rtol 1e-12 evaluation at the top γ was stopped, per the 5 h rule of #4, at l = 5 (γ = 6000, after
  18129 s of worker time) and at l = 6 (γ = 12000, after 22599 s, 1 h 16 min past the cap because the worker could not
  be told from the uncapped φ_cut one earlier). At those two γ the verdict rests on base, φ_cut 1e-12, R = 90, r0 and
  the q task, all with D < 0; the integrator-tolerance variant is missing there. It was tested at the ten lower γ at every
  l, and at γ = 3000 for l = 4, where D = −4.59536726×10⁻¹² against the base −4.59536727×10⁻¹². The capped rows are
  marked in the output files and were added with `cap_row.py` from the raw folder, not a repository script.
- The scan's D at γ_max and the convergence base differ by a relative 6×10⁻⁹ (l = 4), 1.8×10⁻⁸ (l = 5) and
  1.2×10⁻⁷ (l = 6), with γ differing by about 10⁻¹⁵ relative. That is evaluation noise, far below the 2×10⁻⁴ spread.
- The mode condition q_int = q_ext stays far from satisfied, q_int ≈ +γ and q_ext ≈ −γ at every γ, as at l = 2, 3.

**κ_eff** (`kappa_eff_l{4,5,6}_output.txt`, `kappa_eff_l{4,5,6}_w0.7603.csv`; 31 log-spaced γ ∈ [10⁻⁴, 150],
R = r_m/|M| = 1653.57, units |M| = 1):

| l | κ_c(R) | −(l+1)/R | min(κ_eff − κ_c), at γ = 10⁻⁴ | min κ_eff | κ_eff at γ = 150 | max(κ_ext − κ_c) |
|---|---|---|---|---|---|---|
| 4 | −3.0264e-3 | −3.0238e-3 | +6.6174e-3 | +3.5910e-3 | +2.18708 | −3.9e-10 |
| 5 | −3.6314e-3 | −3.6285e-3 | +7.8644e-3 | +4.2330e-3 | +2.18708 | −3.2e-10 |
| 6 | −4.2365e-3 | −4.2333e-3 | +9.0782e-3 | +4.8417e-3 | +2.18708 | −2.7e-10 |

- *Numerical evidence*: **κ_eff − κ_c > 0 at all 31 γ for each of l = 4, 5, 6**, so there is no γ and l among those
  tabulated where it fails. It is smallest at the smallest γ and grows with γ. The exterior's own decaying solution has
  κ_ext < κ_c at all 31 γ, as at l = 2.
- Unlike l = 2 (κ_eff crosses 0 at γ ≈ 2.4×10⁻⁴), κ_eff is positive at every tabulated γ for l = 4, 5, 6. The large-R
  threshold κ_c ≈ −(l+1)/R is reproduced to 0.1 %. By the criterion of §4 (unstable iff κ < κ_c for R > r_s) this is
  the stable side.
- Not done for l ≥ 4: the static coefficient c_l of §4 (`love.py`), so the γ → 0 limit is covered only down to
  γ = 10⁻⁴.

### Complex s at l = 4, 5, 6 (#115)

The real-γ scans above leave complex s open. #115 counts the zeros of D (zeros minus poles, see below) in three of the
boxes of §5.1, at l = 4, 5 and 6, with the current `winding2.py` (the step rule of §5.1, guard 1) at rtol 10⁻⁸ for the
phantom campaign star (ω = 0.7603, M = −0.014580). One box was counted a second time with `--guard=all` (below). The
boxes are the stored box of #7 (Re s ∈ [0.1, 2], |Im s| ≤ 2.5), #5's strip
s1 (Re s ∈ [0.03, 0.1], |Im s| ≤ 2.5) and #5's near-axis box n1 (Re s ∈ [10⁻⁴, 0.1], |Im s| ≤ 0.23). Each has a pair
planted in the same cached evaluations, at 0.3 ± 1.76i, 0.04 ± 0.25i and 0.0005 ± 0.12i, so a box passes only if the
count with the pair is 2 higher. Run on rpi5 (Raspberry Pi 5, aarch64, 4 cores shared with another project, `nice -n 10`,
one or two processes, peak rss 90 MB). Every number below is in `i115_boxes_output.txt`, which `collect_i115.py` builds
from the raw logs (negative-mass-raw, folder I115-rpi5 at 832a9ad; the commands are the lines of `jobs.txt` there). From
`phase2/nonradial`, with the raw repository cloned beside the checkout:

```
python collect_i115.py ../../../negative-mass-raw/I115-rpi5 > i115_boxes_output.txt
```

It exits 1, because s1 at l = 5 failed its planted check at the default guard (below). Its guarded rerun is a row of its
own and passes.

**Control.** `control_winding2.py --cachedir=…` on rpi5 prints PASS on both lines: +1.0000 around the Mini BS I QNM
(max step 0.33 rad) and +0.0000 displaced (0.03 rad).

| l | box | guard | winding | max step | min abs(D) | planted pair found | check | wall |
|---|---|---|---|---|---|---|---|---|
| 4 | stored | 1 | -0.0000 | 0.40 | 3.02e-03 | +2.0000 of 2 | PASS | 42 min |
| 4 | s1 | 1 | +0.0000 | 0.40 | 4.43e-03 | +2.0000 of 2 | PASS | 1 h 50 min |
| 4 | n1 | 1 | +0.0000 | 0.40 | 4.17e-02 | +2.0000 of 2 | PASS | 10 min |
| 5 | stored | 1 | -0.0000 | 0.40 | 4.70e-03 | +2.0000 of 2 | PASS | 46 min |
| 5 | s1 | 1 | -0.0000 | 0.40 | 7.21e-03 | **+0.0000 of 2** | **FAIL** | 1 h 44 min |
| 5 | n1 | 1 | -0.0000 | 0.40 | 7.17e-02 | +2.0000 of 2 | PASS | 11 min |
| 6 | stored | 1 | -0.0000 | 0.40 | 6.71e-03 | +2.0000 of 2 | PASS | 52 min |
| 6 | s1 | 1 | +0.0000 | 0.40 | 8.69e-03 | +2.0000 of 2 | PASS | 1 h 58 min |
| 6 | n1 | 1 | +0.0000 | 0.39 | 1.35e-01 | +2.0000 of 2 | PASS | 9 min |
| 5 | s1, rerun | all | -0.0000 | 0.40 | 7.21e-03 | +2.0000 of 2 | PASS | 1 h 8 min |

(The winding is printed as logged, so a rounded zero has either sign. `max step` is the largest phase step in radians and
`min abs(D)` the smallest |D| on the contour, phase-normalised, so only comparable within a row. Wall is from the first
start to the last end in `queue.log`, with one or two processes running, so it is not a clean cost.) The s1 box at l = 6
was the pilot, the hardest case: 241 evaluations, 1 h 58 min, under the 3 h stop condition. **The other eight boxes, and
the counts of §5 and §5.1, were counted with guard 1** (the half-step test at the first two bisection levels, the rule of
#54); only the last row used `--guard=all`. Recounting them with `--guard=all` is not part of #115.

**Acceptance.**

| criterion | verdict |
|---|---|
| the control passes on rpi5 | **Met.** Both lines PASS |
| each box returns its planted zero | **Not met at the default guard:** 8 of 9, s1 at l = 5 returned 0 of 2. **Met by s1 at l = 5 with `--guard=all`:** +2.0000 of 2, max step 0.40. The other eight were not recounted with `--guard=all` |
| each winding reported with max step and min abs(D), or the box where it is not 0 named | **Met.** All nine report 0 (s1 at l = 5 at both guards), max step 0.39 to 0.40, no box has a non-zero winding |
| a count is a count of zeros only where a pole exclusion covers it, and the summary says which boxes | **None covered.** Poles are not excluded in any box (below) |

**The failed check (s1 at l = 5, guard 1).** `winding2.py` exited 1: the count with the planted pair is −0.0000 instead
of +2.0000. **Cause, in one line:** the lost step is 0.1/4 of an edge wide, so it was accepted at bisection depth 2, where
the half-step test of guard 1 does not run. No threshold was changed. The raw folder's `diagnose_plant.py` replays the
counter on the 212 cached points of that box, with no new evaluation of D (`diagnose_s1.log`, `diagnose_all.log` there;
the script is not in this repository).
- It reproduces the logged cumulative phases of D on all three edges and the planted count of −0.0000.
- On the left edge, Re s = 0.03, one accepted step of D times the planted factor, from Im s = 0.25 to 0.1875, measured
  −0.3826 rad. The 14 cached steps of that stretch (15 points, 13 inside it) sum to +5.9006 rad (largest step
  0.73 rad), so the step lost one 2π and the count lost 2. The planted zero is 0.01 from that edge, so its factor turns
  by about π within one step.
  Summed over the whole contour, the dense chain of the cached points gives +6.2832 rad, a count of +2.0000.
- D alone loses nothing in any of the nine boxes. The accepted steps and the dense chain of all cached points on each
  edge give the same phase to 4 decimals (−22.8537, −0.5455, +23.3992 on the three edges for s1 at l = 5, summing to 0),
  and the largest step of the dense chain is 0.40 rad. The only step of the nine boxes where the two differ is the one
  above, in D times the factor.
This is the defect of §5.1 in a new place: one accepted step loses 2π because the 0.4 rad rule and its half-step test
see only principal values, and the test stops at depth 1. The dense chain shares D and its cache with the count, so it
tests the step rule and not D.

**The guarded rerun (s1 at l = 5, `--guard=all`).** `winding2.py` takes `--guard=all` (706b4c4), which applies the
half-step test at every bisection level; the default 1 is the rule of #54. The rerun is `winding2.py -1 0.7603 5 0.03 0.1
2.5 --plant=0.04,0.25 --guard=all` on the same cache and the same planted zero, which was not moved. It was fixed in
advance that the PLANTED line must read +2.0000 within 0.05 with a max step under 0.45 rad, and that the winding of D is
reported whatever it is.
- PLANTED: count with the factor +2.0000, without −0.0000, difference **+2.0000** of 2, max step 0.40 rad, **PASS**. Exit
  0, 1 h 8 min, 363 evaluations at the last edge (the cache goes from 212 to 387 points, none of the 212 changed).
- RESULT: winding −0.0000, max step 0.40 rad, min abs(D) 7.21e-03. The cumulative phases of D on the three edges are
  the same to 4 decimals as at guard 1 (−22.8537, −23.3992, −0.0000), so the guard changed the planted count and not D's.
- The default-guard run `s1_l5.log` stays as it is. It is the mutation that shows the check can fail: the same cache
  and plant give a count of 0 at guard 1 and +2 at guard `all`.
So s1 at l = 5 passes its planted check with `--guard=all`, and the guard leaves D's count there unchanged. The other
eight boxes passed their check with guard 1 and are not recounted here; §5.4 (#126) recounts them with `--guard=all`,
and none changes.

**What the planted check does not test** (coordinator's review of #115). The planted check compares two counts made on
the same accepted steps. Suppose a step loses 2π in D alone, and the planted factor turns little over that step. Then
the step loses 2π in D times the factor as well, and the difference is unchanged. The review removed the cached points
inside the lost stretch from the 212-point copy. D alone then lost 2π there too, and the difference stayed 2. So a
passing planted check shows that the counter finds a pair at s0. It does not show that D's own steps lost nothing; only
the step rule protects D's count. §5.4 (#126) recounts the guard-1 boxes with `--guard=all` and puts that replay in
the repository.

**Controls of the collector.** On a synthetic raw folder built from #5's logs (tags renamed to l = 4, 5, 6) the collector
exits 0, and each of `--mutate=noplant`, `wrongexpect`, `wrongbox`, `shift` and `nocontrol` exits 1. On the real logs
each of the first four lists all ten rows as FAILED, `nocontrol` lists the control and the one real failure, and
`--mutate=wrongguard` (expects guard all everywhere) lists the nine default rows.

**Pole caveat.** Each winding is zeros minus poles of D. D's known pole source is the Zerilli column (1, q/f). Its q has
a pole at each zero of the Zerilli solution at r_m (#54, #121); §5.5 (#122) checks D's other factors one by one. #121
counted those zeros in the 21 boxes of §5 and §5.1, and #124 in these nine (§5.3). #124 finds no zero of Z(r_m; s) in
any of the nine, so the Zerilli column gives D no pole there. Each 0 above then counts zeros of D, conditional on that
column being D's only pole source (§5.5).

**Not covered.** Only the stored box, s1 and n1 were counted. The other boxes of §5.1 (n3, h1, h2, s2, s3 and the
threshold boxes t1 and t2 at two outer radii) are not run for l = 4, 5, 6; Re s < 0.03 for 0.23 < |Im s| ≤ 2.5 and
|Im s| > 2.5 are open, as in §5.1. The real-γ axis at l = 4, 5, 6 is the 100-point scan above, not refined.

Status: *numerical evidence*, conditional on the Zerilli column being D's only pole source (§5.3, §5.5). Eight boxes
passed their planted check at guard 1; s1 at l = 5 failed it at guard 1 and passed it with `--guard=all`. §5.4 (#126)
recounts the eight, and #124's counts of Z, with `--guard=all`, and none changes. A 2π lost inside a step whose two
halves agree is excluded by neither guard (§5.4).

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
| phantom ω = 0.85 | 3 | same | **0.0000** (§5.2; #7's −2.0000 was the first step rule) |

(`winding2_output.txt` for the first three rows, #7's I7-win2 at 60553bc for the rest; max phase step 0.4 rad; D(x0),
D(x1) have the same sign. §5.2 recounts every row with the corrected step rule, and only the last row changes.)

**Discrepancy (#7), resolved by #54.** #7 counted −2 for ω = 0.85, l = 3. A count of zeros cannot be negative, so this
was read as a net pole of D in the box, with a hypothesis for #54 to test: the Zerilli column (1, q/f) of
`polar_solver.decaying_plane` has a pole wherever the decaying Zerilli solution vanishes at r_m. The recount of §5.2 gives
0, and the −2 was the first step rule losing one 2π on the left edge. The column has such a pole by construction
(q = f Z′/Z), and §5.3 (#121) finds no zero of Z(r_m; s) in any of these boxes, so it gives D no pole there. Each
winding counts zeros minus poles, so the 0 counts are counts of zeros of D, conditional on the Zerilli column being D's
only pole source (checked factor by factor in §5.5, #122).

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
`collect_i5.py` builds from the raw logs (negative-mass-raw, folder I5-rpi5, at e2ebcdd).

The commands of every run are the lines of `jobs.txt` there. From `phase2/nonradial` with the raw repository cloned beside
the checkout, `python collect_i5.py ../../../negative-mass-raw/I5-rpi5 > i5_boxes_output.txt` rebuilds the file and exits
non-zero if a check failed or a run is missing.

Result: all 17 boxes count 0 (zeros minus poles), each with its planted zero returned by the same code, and the real-γ
axis shows no sign change at 397 points (step ratio 1.032 for l = 2, 1.036 for l = 3). Where the counts stop is under
"Not covered" below.

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
were computed with the first rule and are recounted in §5.2 (#54), where the −2 turns out to be one such loss.

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

`boxes_graph.py` draws the same boxes at true scale from `i5_boxes_output.txt`, with the planted zeros.

![The counted boxes of the growing half plane at true scale, l = 2 and 3](complex_s_boxes.svg)

| job | what | Re s | Im s | winding | max step | min abs(D) | planted s0 | planted found | check |
|---|---|---|---|---|---|---|---|---|---|
| rep_l2 | stored box, l = 2 (#7 box of the first scans) | [0.1, 2.0] | abs(Im s) ≤ 2.5 | +0.0000 | 0.39 | 6.70e-04 | 0.3+1.76j | +2.0000 of 2 | PASS |
| s1_l2 | strip, l = 2 | [0.03, 0.1] | abs(Im s) ≤ 2.5 | +0.0000 | 0.40 | 1.00e-03 | 0.04+0.25j | +2.0000 of 2 | PASS |
| s1_l3 | strip, l = 3 | [0.03, 0.1] | abs(Im s) ≤ 2.5 | +0.0000 | 0.40 | 2.49e-03 | 0.04+0.25j | +2.0000 of 2 | PASS |
| n1_l2 | below threshold, right of the axis, l = 2 | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | -0.0000 | 0.33 | 2.82e-03 | 0.0005+0.12j | +2.0000 of 2 | PASS |
| n1_l3 | below threshold, right of the axis, l = 3 | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | +0.0000 | 0.35 | 3.16e-02 | 0.0005+0.12j | +2.0000 of 2 | PASS |
| n3_l2 | below threshold, across the axis, l = 2 | [-0.001, 0.0001] | [0.01, 0.23] | +0.0000 | 0.39 | 2.75e-03 | -0.0005+0.12j | +1.0000 of 1 | PASS |
| n3_l3 | below threshold, across the axis, l = 3 | [-0.001, 0.0001] | [0.01, 0.23] | -0.0000 | 0.38 | 3.66e-02 | -0.0005+0.12j | +1.0000 of 1 | PASS |
| h1_l2 | high Im s, l = 2 | [0.1, 2.0] | [2.5, 6.0] | +0.0000 | 0.38 | 1.05e-04 | 0.3+2.7j | +1.0000 of 1 | PASS |
| h1_l3 | high Im s, l = 3 | [0.1, 2.0] | [2.5, 6.0] | -0.0000 | 0.38 | 2.63e-04 | 0.3+2.7j | +1.0000 of 1 | PASS |
| t1_r40 | threshold 1 - w, outer radius 40, l = 2 | [0.02, 0.1] | [0.15, 0.33] | +0.0000 | 0.38 | 3.73e-03 | 0.03+0.25j | +1.0000 of 1 | PASS |
| t1_r80 | threshold 1 - w, outer radius 80, l = 2 | [0.02, 0.1] | [0.15, 0.33] | -0.0000 | 0.39 | 3.73e-03 | 0.03+0.25j | +1.0000 of 1 | PASS |
| t2_r40 | threshold 1 + w, outer radius 40, l = 2 | [0.02, 0.1] | [1.68, 1.84] | +0.0000 | 0.40 | 9.30e-04 | 0.03+1.77j | +1.0000 of 1 | PASS |
| t2_r80 | threshold 1 + w, outer radius 80, l = 2 | [0.02, 0.1] | [1.68, 1.84] | +0.0000 | 0.35 | 9.30e-04 | 0.03+1.77j | +1.0000 of 1 | PASS |
| s2_l2 | strip, closer to the axis, l = 2 | [0.01, 0.03] | abs(Im s) ≤ 2.5 | +0.0000 | 0.40 | 5.30e-04 | 0.015+0.25j | +2.0000 of 2 | PASS |
| h2_l2 | high Im s, second range, l = 2 | [0.1, 2.0] | [6.0, 10.0] | -0.0000 | 0.31 | 2.57e-05 | 0.3+6.2j | +1.0000 of 1 | PASS |
| s3_l2 | strip above Im s = 2.5, l = 2 | [0.03, 0.1] | [2.5, 6.0] | -0.0000 | 0.39 | 9.90e-05 | 0.05+3j | +1.0000 of 1 | PASS |
| s3_l3 | strip above Im s = 2.5, l = 3 | [0.03, 0.1] | [2.5, 6.0] | -0.0000 | 0.38 | 2.71e-04 | 0.05+3j | +1.0000 of 1 | PASS |

`max step` is the largest phase step in radians and `min abs(D)` the smallest |D| on the contour (phase-normalised, so
only comparable within a row).

Winding is zeros minus poles of D in the box (§5, #54). §5.3 (#121) finds no pole from the Zerilli column in any of
these boxes, and #123 none in `n3_l2` and `n3_l3`, which cross the axis ("Boxes across Re s = 0" in §5.3). So each 0
means no zero, conditional on that column being D's only pole source (§5.5).

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
`PolarProblem.D`. Three stored points recomputed agree to a relative 3.4×10⁻⁹ at worst.

| l | level | step ratio | points | sign changes | ratio of abs(D_i) to the geometric mean of its neighbours |
|---|---|---|---|---|---|
| 2 | 0 | 1.136 | 100 | 0 | not defined |
| 2 | 1 | 1.066 | 199 | 0 | 1.000 to 1.002 |
| 2 | 2 | 1.032 | 397 | 0 | 1.000 to 1.001 |
| 3 | 0 | 1.150 | 100 | 0 | not defined |
| 3 | 1 | 1.072 | 199 | 0 | 1.000 to 1.003 |
| 3 | 2 | 1.036 | 397 | 0 | 1.000 to 1.001 |

Both ranges give no sign change at any level (l = 2: γ ∈ [10⁻³, 300], l = 3: γ ∈ [10⁻³, 1000]).

The planted pair (γ = 0.6035 and 0.6408 for l = 2, 1.115 and 1.184 for l = 3, 6 % apart) lies inside one level-0
interval, is missed there, and is found in the two intervals that contain it at level 1 and level 2. With `--mutate=coarse`
it is missed, with `--mutate=pair` an unplanted pair makes the count change between levels, and with `--mutate=stale`
the recomputed values differ from the stored ones by 1 %, each failing its check as designed. The ratio of |D_i| to
the geometric mean of its neighbours stays within 1.000 to 1.003, so ln |D| is smooth at every grid point. A pair of
roots closer than about 3 % in γ is therefore excluded only by that smoothness, not by a sign change.

#### Beyond Re s = 2 and |Im s| = 10 (#117)

#117 counts four boxes outside those above, for the phantom campaign star (ω = 0.7603, M = −0.014580): Re s ∈ [2, 5] with
|Im s| ≤ 10 for l = 2 and 3 (`b25_l2`, `b25_l3`), Re s ∈ [0.1, 2] with Im s ∈ [6, 10] for l = 3 (`h2_l3`; `h2_l2` above
is its l = 2 counterpart) and Re s ∈ [0.1, 2] with Im s ∈ [10, 20] for l = 2 (`h3_l2`). Each box is one `winding2.py` job
(zeros minus poles of D) and one `zerilli_zeros.py` job (zeros of Z(r_m; s), D's poles, §5.3) with the same arguments and
`--guard=all` (§5.4), each with a planted zero.

Result: **all four boxes have winding 0 and 0 zeros of Z(r_m; s), so 0 zeros of D**, and all eight runs return their
planted zero (+2.0000 of 2 in the symmetric boxes, +1.0000 of 1 in the upper boxes). No count is nonzero. Status:
numerical evidence, conditional as §5.1's counts are (that the Zerilli column is D's only pole source; §5.5 checks
D's other factors in its 61 boxes, not in these four).

**What ran.** On rpi5 on 2026-11-21, 16:10 to 23:24 (UTC−8), with `winding2.py` and `zerilli_zeros.py` of
`wip/i117-boxes` at a05cd0b, unmodified, at nice level 10 and one thread, at most two processes at a time and the second
only while the 15-minute load was below 4 (`queue.log` gives the load at each second-process start, all below 4). The
boxes and planted zeros were fixed before the runs and not moved. Each job is `python -u <script> <arguments>
--cachedir=<raw folder>/cache`, with the eight lines `job|script|arguments` that `python collect_i117.py --jobs` prints
(raw `jobs.txt`), for example `winding2.py -1 0.7603 2 0.1 2.0 20.0 --y0=10.0 --plant=0.3,10.5 --guard=all`. The pilot
`w_h3_l2` (the largest |s| and the slow left edge, stop condition 8 h) ran alone for 3 h 55 min. The coordinator started
the other seven jobs with the same launcher while the runner session waited for a usage reset. That first start of
`w_b25_l2` and `w_b25_l3` died within 20 s and left empty logs and no cached point (by the coordinator's account,
killed when the process that started them returned, with no cache file of either box at 21:35). Raw 13360d3 holds no
cache file of them, and the 20 and 18 points in raw 705847e were committed at 21:37:58, after the restart, so they are
the restart's. The two jobs were restarted with the same command at 21:36 and began from an empty cache, so their walls
below are those of complete runs. The other six jobs ran once, and all eight exited 0.
Logs, caches and the launcher are in negative-mass-raw, folder `I117-rpi5` at 28863c5. From `phase2/nonradial`,
with the raw repository cloned beside the checkout:

```
python collect_i117.py ../../../negative-mass-raw/I117-rpi5 > i117_boxes_output.txt
```

The collector exits 0. It exits 1 when `jobs.txt` differs from the job list, a log or exit status is missing or not 0, a
count is not within 0.05 of an integer, a planted check fails, or D has opposite signs at the two ends of a symmetric box.
Its `--mutate=noplant` (the count with the planted factor read as the count without it), `shift`, `nolog`, `badbox`, `rc`,
`zcheck`, `signs` and `noguard` each exit 1 (raw `collect_mutate_<name>.log`). The four boxes:

| box | l | Re s | Im s | winding of D | max step | min abs(D) | evaluations | R_s | planted s0 | planted found (D) | count of Z | planted found (Z) | zeros of D |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| b25_l2 | 2 | [2, 5] | abs(Im s) ≤ 10 | +0.0000 | 0.39 | 1.95e-05 | 165 | 44.0 | 2.5+1j | +2.0000 of 2 | -0.0000 | +2.0000 of 2 | 0 |
| b25_l3 | 3 | [2, 5] | abs(Im s) ≤ 10 | -0.0000 | 0.40 | 4.87e-05 | 231 | 44.0 | 2.5+1j | +2.0000 of 2 | +0.0000 | +2.0000 of 2 | 0 |
| h2_l3 | 3 | [0.1, 2] | [6, 10] | +0.0000 | 0.38 | 6.43e-05 | 80 | 422.4 | 0.3+6.2j | +1.0000 of 1 | -0.0000 | +1.0000 of 1 | 0 |
| h3_l2 | 2 | [0.1, 2] | [10, 20] | -0.0000 | 0.28 | 3.38e-06 | 120 | 423.6 | 0.3+10.5j | +1.0000 of 1 | -0.0000 | +1.0000 of 1 | 0 |

`Re s`, `Im s` and `abs(Im s) ≤ y` are as in §5.1 (a symmetric box is doubled from its upper half, an upper box counted
on its four edges). `max step` is the largest phase step of D in radians and `min abs(D)` the smallest |D| on the
contour, phase-normalised and comparable within a row only. `evaluations` is the cache size when the count ended, and
the planted pass adds points. The planted found columns give the difference of the counts with and without the factor,
as in §5.1. Zeros of D = winding of D + zeros of Z (§5.3). The eight jobs:

| job | script | exit status | max step | min abs | evaluations | wall | D(x0), D(x1) | checks |
|---|---|---|---|---|---|---|---|---|
| w_b25_l2 | winding2.py | 0 | 0.39 | 1.95e-05 | 165 | 0 h 18 min | -1.773e-03, -1.895e-04 | PASS |
| w_b25_l3 | winding2.py | 0 | 0.40 | 4.87e-05 | 231 | 0 h 21 min | -4.425e-03, -4.739e-04 | PASS |
| w_h2_l3 | winding2.py | 0 | 0.38 | 6.43e-05 | 80 | 1 h 30 min |  | PASS |
| w_h3_l2 | winding2.py | 0 | 0.28 | 3.38e-06 | 120 | 3 h 55 min |  | PASS |
| z_b25_l2 | zerilli_zeros.py | 0 | 0.23 | 2.99e-01 | 61 | 24 s |  | PASS; holo defect ≤ 1.7e-08, q against zerilli_q ≤ 9.2e-12 |
| z_b25_l3 | zerilli_zeros.py | 0 | 0.23 | 2.99e-01 | 61 | 26 s |  | PASS; holo defect ≤ 1.6e-08, q against zerilli_q ≤ 1.6e-12 |
| z_h2_l3 | zerilli_zeros.py | 0 | 0.02 | 3.13e-01 | 80 | 36 s |  | PASS; holo defect ≤ 2.4e-08, q against zerilli_q ≤ 4.0e-12 |
| z_h3_l2 | zerilli_zeros.py | 0 | 0.01 | 2.23e-01 | 80 | 32 s |  | PASS; holo defect ≤ 6.1e-09, q against zerilli_q ≤ 1.5e-12 |

In the two symmetric boxes D has one sign at both ends of the real segment, as the doubling needs. In the Z jobs `max
step` and `min abs` are those of Z, and the Cauchy–Riemann defect (limit 10⁻⁶) and the difference of q from `zerilli_q`
(limit 10⁻⁶, applied wherever Re s ≥ 0.03, which is all three points of every box here) are the largest of the three
points. The four winding jobs took 6.1 h of wall time in all, and the four Z jobs about 2 min. The pilot `h3_l2` took
3 h 55 min for 120 evaluations and its planted pass; `memwatch.log` in the raw folder has the load and free memory
every minute.

| acceptance | verdict |
|---|---|
| each box returns its planted zero, with its winding, max phase step and min abs(D) reported | met in 4 of 4 boxes: +2.0000 of 2 (`b25_l2`, `b25_l3`) and +1.0000 of 1 (`h2_l3`, `h3_l2`) for D, the same for Z; max step 0.28 to 0.40 rad (D) and 0.01 to 0.23 rad (Z) against the 0.45 limit of the planted check; smallest min abs(D) 3.38×10⁻⁶ (`h3_l2`) |
| any nonzero count is reported with its box | no count is nonzero: all four windings and all four counts of Z are 0, so the Zerilli column gives D no pole in these boxes and each has 0 zeros of D; D's other factors are checked in §5.5's 61 boxes, not in these four |
| the collector exits 1 on a missing log, a winding off an integer or a failed planted check, and `--mutate=noplant` exits 1 | met: exit 0 on the logs, and each of its eight mutations exits 1 |

What this does not show. The mutations of the counter itself (reverse, wrong doubling, misplaced and one-sided plant,
§5.1) were run on the stored box and not repeated here. A 2π lost inside a step whose two halves agree is excluded by
neither guard, and no dense-chain replay (§5.4, `diagnose_plant.py`) of these four caches has been made.
`boxes_graph.py` and the figure above draw §5.1's boxes only, not these four.

*Coordinator's review (#117's closing comment).* The collector, rerun on the raw logs, reproduces
`i117_boxes_output.txt` byte for byte, and its `--mutate=noplant` exits 1. `w_b25_l2` and `z_b25_l2`, rerun from an empty
cache on another machine, give the same counts and the same 195 points of D, which agree to 2.6×10⁻¹⁰ relative. An
independent counter, which shares only `background.solve_star` and `PolarProblem.D` with these scripts, walks the full
rectangle of `b25_l2` without the doubling: 530 points, every phase step at most 0.3 rad, and every segment bisected once
more. It finds winding 0, two planted zeros at +2.0000 and a planted pole at −1.0000. D(s̄) = conj D(s) holds to
6.5×10⁻¹⁴ on 235 mirrored pairs, so the doubling holds in that box. Its mutations `reverse` and `noplant` exit 1. The
script and its logs are in negative-mass-raw, folder `I117-mac`.

#### Not covered

- Re s below the lower edges for 0.23 < |Im s| ≤ 2.5. For l = 2 that is 0 < Re s < 0.01 (the near-axis box n1 reaches
  10⁻⁴ only below Im s = 0.23), for l = 3 it is 0 < Re s < 0.03. R_s grows as 1/Re s (3831 at 0.01), and so does the
  cost of a box, which took 3.7 h at 0.01. A mode that close to marginal is neither shown nor excluded by these counts.
- The strip Re s < 0.03 for 2.5 < |Im s|, and the strip Re s < 0.1 for |Im s| > 6.
- Re s ∈ [0.1, 2] above |Im s| = 20 (l = 2) and 10 (l = 3), Re s ∈ [2, 5] above |Im s| = 10, and Re s > 5 (#117's
  edges).
- The sliver |Im s| < 0.01 left of Re s = 10⁻⁴.
- Poles of D. Every count is zeros minus poles. §5.3 (#121, and #123 for `n3_l2` and `n3_l3`) excludes the Zerilli
  column's poles in every box. §5.5 (#122) checks D's other factors one by one; its test that arg D is harmonic fails at
  Re s ≥ 0.1 (#129).
- l ≥ 4, except the stored box, s1 and n1 at l = 4, 5, 6 (§3a, #115).

Status: numerical evidence, conditional on the Zerilli column being D's only pole source (§5.3, §5.5). In `n3_l2` and
`n3_l3` the winding of D used `zerilli_q` for that column on both sides of Re s = 0, and on the left it differs from the
continued solution by up to 8.9×10⁻⁵ (§5.3, #123); §5.5 (#122) finds `zerilli_q` the non-harmonic factor there, and with
the continued column the windings stay 0. Every count here used guard 1. §5.4 (#126) recounts each with `--guard=all`,
and none changes; a 2π lost inside a step whose two halves agree is excluded by neither guard.

### 5.2 Recount of §5's windings with the corrected step rule (#54)

#54 started from #7's winding of −2 for ω = 0.85, l = 3 and the hypothesis that D has a pole there. #5 had found that the
first step rule can lose 2π in one step (§5.1), and one loss shifts a symmetric box's count by ±2. So the boxes of the §5
table were recounted first, with the current `winding2.py` (the rule of §5.1, git 59d26de), each with the zero 0.3 + 1.76i
and its conjugate planted in the same cached evaluations.

Result: all seven boxes count 0. Only ω = 0.85, l = 3 changed, from −2 to 0, and the −2 was the first step rule, not a
pole. The search for a pole (steps 1 to 3 of #54) was therefore not run.

Six boxes were counted on rpi5. The seventh, phantom ω = 0.7603 with l = 2, is #5's `rep_l2` run. The ω = 0.70 and 0.85
counts reuse #7's stored evaluations of D (byte copies of I7-win2's caches) and add the ones the new rule asks for. Every
number below is in `i54_recount_output.txt`, which `collect_i54.py` builds from the raw logs (negative-mass-raw, folders
I54-rpi5 at a0cfa51, I5-rpi5 and I7-win2). From `phase2/nonradial`, with the raw repository cloned beside the checkout:

```
python collect_i54.py ../../../negative-mass-raw > i54_recount_output.txt
```

It exits non-zero if a run is missing, a planted zero is not returned or a winding is not an integer. All boxes are
Re s ∈ [0.1, 2], |Im s| ≤ 2.5. The old values are #7's logs for ω = 0.70 and 0.85 and `winding2_output.txt` for the
others, all made with the first rule.

| star | l | old winding | new winding | max step, old / new (rad) | min abs(D), old / new | planted pair found |
|---|---|---|---|---|---|---|
| phantom ω = 0.85 | 3 | **−2.0000** | **+0.0000** | 0.40 / 0.40 | 1.23e-03 / 1.23e-03 | +2.0000 of 2 |
| phantom ω = 0.85 | 2 | −0.0000 | −0.0000 | 0.40 / 0.40 | 4.93e-04 / 4.93e-04 | +2.0000 of 2 |
| phantom ω = 0.70 | 3 | −0.0000 | +0.0000 | 0.40 / 0.40 | 1.86e-03 / 1.86e-03 | +2.0000 of 2 |
| phantom ω = 0.70 | 2 | −0.0000 | +0.0000 | 0.40 / 0.40 | 7.45e-04 / 7.45e-04 | +2.0000 of 2 |
| phantom ω = 0.7603 | 3 | +0.0000 | +0.0000 | 0.40 / 0.40 | 1.68e-03 / 1.68e-03 | +2.0000 of 2 |
| canonical ω = 0.75 | 2 | +0.0000 | +0.0000 | 0.40 / 0.40 | 6.88e-04 / 6.88e-04 | +2.0000 of 2 |
| phantom ω = 0.7603 | 2 | +0.0000 | +0.0000 | 0.39 / 0.39 | 6.70e-04 / 6.70e-04 | +2.0000 of 2 |

The old and new max step and min abs(D) agree in every row. Every planted pair is returned.

#### The lost step

`replay_old_rule_i54.py` replays the first rule on the stored evaluations of #7's four boxes without evaluating D again
(`i54_replay_output.txt`, from `phase2/nonradial`: `python replay_old_rule_i54.py ../../../negative-mass-raw`).

- It reproduces #7's logged winding, max step and min abs(D) in all four boxes, and asks for exactly the 77, 63, 104 and
  67 points of #7's caches.
- The cached points on each edge, sorted along it, form a dense chain whose steps are all at most 0.40 rad. Its winding
  is 0 in all four boxes, equal to the current rule's.
- Three boxes lose nothing. In ω = 0.85, l = 3 one step lost 2π: the left edge from Im s = 0.25 to 0 at Re s = 0.1. The
  first rule measured −0.200 rad and the true change is +6.083 rad. The 27 stored points of that stretch have |D| between
  2.15×10⁻² and 3.02×10⁻², and all 26 steps between them are positive, 0.199 to 0.397 rad. The old count of −2 is the
  current 0 minus 2 for the lost turn.

The left-edge phase of the three l = 3 boxes is now +13.17 (ω = 0.85), +13.60 (ω = 0.70) and +13.17 (ω = 0.7603). For
ω = 0.85 the first rule gave +6.88. This is a plausibility check on the corrected count, not a test.

#### Controls

- **Planted zero.** Each rerun returns the pair planted at 0.3 + 1.76i (+2.0000 of 2, PASS), with the same code as §5.1.
  Its four mutations were run there on the stored box and not repeated here.
- **Same D on two machines.** #7's evaluations were made on win2 (x86, Windows) and the new ones on rpi5 (aarch64).
  `python check_cache_i54.py ../../../negative-mass-raw/I7-win2` recomputes three stored points of each of the four caches
  on rpi5 with the same D. The worst relative difference is 1.9×10⁻¹¹ (limit 10⁻⁶, `i54_cache_check_output.txt` keeps its
  summary lines). With `--mutate=scale`, which multiplies the recomputed values by 1.001, it exits 1.
- **Mutations of the new scripts.** `replay_old_rule_i54.py --mutate=thr` (threshold 0.3 instead of 0.4) asks for points
  that are not stored and exits 1. `collect_i54.py --mutate=noplant` and `--mutate=wrongexpect` exit 1.
- **What is shared.** The counts, the dense chain and the planted zero all use the same D and the same counter, so they
  test the step rule and not D. The two-machine recomputation tests D at 12 points. The Mini BS I controls of §5.1 are
  independent of the campaign stars.

#### What this settles and what it does not

- Status: numerical evidence that D has winding 0 in each of the seven boxes of §5. The −2 of ω = 0.85, l = 3 was the
  first step rule, located at one step.
- The planted-pole control of #54 (a pole-free column returning +1 for a planted zero and 0 for a planted pole) was not
  run, as it was not needed.
- The hypothesis that the Zerilli column (1, q/f) has a pole where ψ(r_m; s) vanishes is neither confirmed nor refuted.
  Its only evidence was the −2. A winding of 0 is still zeros minus poles, and nothing here searches for poles.
- The passages this recount supersedes were edited by the coordinator at the merge: the §5 table row ω = 0.85, l = 3,
  the paragraph "Discrepancy (#7)", the "Not covered" list of §5.1 and §6 item 2. The "conditional on D having no pole"
  of §5 and of the "Status" line in §5.1 still held, since poles were not searched for. §5.3 (#121) has since counted
  the zeros of Z(r_m; s), and the coordinator narrowed those passages at its merge.

### 5.3 Poles of D (#121)

`polar_solver.decaying_plane` builds the Zerilli column of D as (1, q/f) with q = f Z′/Z at r_m, so D has a pole at each
zero of the decaying Zerilli solution Z(r_m; s), and the windings of §5, §5.1 and §5.2 count the zeros of D minus these.
#121 counts the zeros of Z(r_m; s) in each of those boxes with `zerilli_zeros.py`. It integrates Z inward from a radius
R_z fixed for the box, in WKB-scaled variables with the implicit Radau method, and counts the phase of
y₁(r_m) = Z(r_m; s) exp Ψ(r_m) with the step rule of §5.1. y₁(r_m) is holomorphic in s and has the zeros of Z(r_m; s) for
Re s > 0. The script's docstring has the details.

Result: **all 21 counted boxes have 0 zeros of Z(r_m; s)**, each with its planted zero returned. So the Zerilli column
gives D no pole in any of them, and the zeros of D equal the winding of D, 0, in each. No count is nonzero, so no box
needs a corrected count. Status: numerical evidence.

**What ran.** On win (surface-laptop-4) on 2026-11-20, with the code of `wip/i121-zerilli-zeros` at d801d1e, unmodified.
Each job is `python -u zerilli_zeros.py <arguments> --cachedir=<raw folder>/cache`, with the 25 lines `job|arguments`
that `python collect_i121.py --jobs` prints (21 boxes and 4 controls; raw `jobs.txt`). The four controls ran first, then
the pilot `n1_l3` (the largest R_z, 1.0×10⁶; 225 s in all, 146 s of it the count), then the other 20 boxes, at most two
processes at a time. Each job ran once and none died. Logs, caches and launchers are in negative-mass-raw, folder
`I121-win` at bc19469. From `phase2/nonradial`, with the raw repository cloned beside the checkout:

```
python collect_i121.py ../../../negative-mass-raw/I121-win > i121_zeros_output.txt
```

On Windows with Python 3.14 the redirect needs `PYTHONIOENCODING=utf-8` (the table has ≥ and ≤, and the default code page
raises `UnicodeEncodeError`). The collector exits 0 here. It exits 1 when a log is missing, a check is not as required or a count is not an integer.
Its `--mutate=box` (it expects `s1_l2` at Re s up to 0.2) and `--mutate=control` (it expects `mut_lshift` to pass) each
exit 1, as they should.

| job | § | star, l | Re s | Im s | zeros of Z | max step | min abs(Z) | planted found | winding of D | zeros of D |
|---|---|---|---|---|---|---|---|---|---|---|
| p7603_l2 | 5 | 0.7603, 2 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | 0.36 | 5.59e-01 | +2.0000 of 2 | +0.0000 | 0 |
| p7603_l3 | 5 | 0.7603, 3 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | 0.29 | 5.59e-01 | +2.0000 of 2 | +0.0000 | 0 |
| c075_l2 | 5 | canonical 0.75, 2 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | 0.36 | 5.59e-01 | +2.0000 of 2 | +0.0000 | 0 |
| p070_l2 | 5 | 0.70, 2 | [0.1, 2] | abs(Im s) ≤ 2.5 | -0.0000 | 0.36 | 5.59e-01 | +2.0000 of 2 | +0.0000 | 0 |
| p070_l3 | 5 | 0.70, 3 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | 0.28 | 5.59e-01 | +2.0000 of 2 | +0.0000 | 0 |
| p085_l2 | 5 | 0.85, 2 | [0.1, 2] | abs(Im s) ≤ 2.5 | -0.0000 | 0.38 | 5.59e-01 | +2.0000 of 2 | -0.0000 | 0 |
| p085_l3 | 5 | 0.85, 3 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | 0.32 | 5.59e-01 | +2.0000 of 2 | +0.0000 | 0 |
| s1_l2 | 5.1 | 0.7603, 2 | [0.03, 0.1] | abs(Im s) ≤ 2.5 | -0.0000 | 0.36 | 6.32e-01 | +2.0000 of 2 | +0.0000 | 0 |
| s1_l3 | 5.1 | 0.7603, 3 | [0.03, 0.1] | abs(Im s) ≤ 2.5 | +0.0000 | 0.35 | 6.33e-01 | +2.0000 of 2 | +0.0000 | 0 |
| n1_l2 | 5.1 | 0.7603, 2 | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | -0.0000 | 0.23 | 2.04e+00 | +2.0000 of 2 | -0.0000 | 0 |
| n1_l3 | 5.1 | 0.7603, 3 | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | -0.0000 | 0.22 | 2.07e+00 | +2.0000 of 2 | +0.0000 | 0 |
| h1_l2 | 5.1 | 0.7603, 2 | [0.1, 2] | [2.5, 6] | -0.0000 | 0.04 | 3.98e-01 | +1.0000 of 1 | +0.0000 | 0 |
| h1_l3 | 5.1 | 0.7603, 3 | [0.1, 2] | [2.5, 6] | -0.0000 | 0.04 | 3.98e-01 | +1.0000 of 1 | -0.0000 | 0 |
| t1_z1 | 5.1 | 0.7603, 2 | [0.02, 0.1] | [0.15, 0.33] | +0.0000 | 0.04 | 1.73e+00 | +1.0000 of 1 | +0.0000 | 0 |
| t1_z2 | 5.1 | 0.7603, 2 | [0.02, 0.1] | [0.15, 0.33] | +0.0000 | 0.04 | 1.73e+00 | +1.0000 of 1 | -0.0000 | 0 |
| t2_z1 | 5.1 | 0.7603, 2 | [0.02, 0.1] | [1.68, 1.84] | -0.0000 | 0.00 | 7.37e-01 | +1.0000 of 1 | +0.0000 | 0 |
| t2_z2 | 5.1 | 0.7603, 2 | [0.02, 0.1] | [1.68, 1.84] | +0.0000 | 0.00 | 7.37e-01 | +1.0000 of 1 | +0.0000 | 0 |
| s2_l2 | 5.1 | 0.7603, 2 | [0.01, 0.03] | abs(Im s) ≤ 2.5 | +0.0000 | 0.31 | 6.33e-01 | +2.0000 of 2 | +0.0000 | 0 |
| h2_l2 | 5.1 | 0.7603, 2 | [0.1, 2] | [6, 10] | +0.0000 | 0.02 | 3.13e-01 | +1.0000 of 1 | -0.0000 | 0 |
| s3_l2 | 5.1 | 0.7603, 2 | [0.03, 0.1] | [2.5, 6] | +0.0000 | 0.00 | 4.08e-01 | +1.0000 of 1 | -0.0000 | 0 |
| s3_l3 | 5.1 | 0.7603, 3 | [0.03, 0.1] | [2.5, 6] | +0.0000 | 0.00 | 4.08e-01 | +1.0000 of 1 | -0.0000 | 0 |

Stars are phantom (ε = −1) unless marked. A box with abs(Im s) ≤ y is symmetric and doubled from its upper half, one with
[a, b] is an upper box (as in §5.1). The winding of D is §5.2's recount (`i54_recount_output.txt`) for the §5 rows and
§5.1's (`i5_boxes_output.txt`) for the others, and zeros of D = winding of D + zeros of Z. `i121_zeros_output.txt` also
has the evaluations, the times and the two check columns. In `t1` and `t2` the start radius is doubled in the second row
(`z2`, #5's outer radii 40 and 80), and the counts and min abs(Z) agree.

The counted function is WKB-scaled, so `min abs(Z)` is its smallest modulus on the contour and is not a distance to a
zero. In the §5 boxes it is 5.59×10⁻¹ for every star and l, which equals |2 + 2.5i|^(−1/2), the WKB amplitude k^(−1/2)
at the far corner, to three digits.
`p7603_l2` shows 67 evaluations and 0 s because `mut_reverse` ran first and shares its cache (the reversed contour of the
same box). Its own 10 s is in `mut_reverse`'s log.

**Checks in every box.** `count`: each count within 0.05 of an integer, max step at most 0.45 rad (largest 0.38), and Z
real with one sign at both ends of a symmetric box. `holo`: Cauchy–Riemann defect at most 10⁻⁶ (worst 1.3×10⁻⁷, `s3_l3`).
`column`: q against `zerilli_q` at most 10⁻⁶ where Re s ≥ 0.03 (worst 1.6×10⁻⁹, `s3_l2`). `plant`: the planted zero
returned. Where Re s < 0.03 the difference in q is reported and not checked, since `zerilli_q` starts at a per-point radius
that leaves a WKB error there: at most 7.0×10⁻⁷ (`n1_l2`, s = 10⁻⁴ + 0.056i).

#### Controls

On the stored box (phantom 0.7603, l = 2, Re s ∈ [0.1, 2], |Im s| ≤ 2.5, planted 0.3 + 1.76i), each with its planted pair:

| job | what | ends | result |
|---|---|---|---|
| `flat` | `--flat`, Mext = 0 | all checks passed | `flat` count 0, `flatq` against r k_l(s r) at most 1.6×10⁻¹², planted +2.0000 of 2 |
| `mut_reverse` | contour clockwise | FAILED: plant | counted −2.0000 of +2 |
| `mut_nonorm` | per-point start without the WKB scaling | FAILED: holo | Cauchy–Riemann defect 1.3×10⁻² to 1.6×10⁻² |
| `mut_lshift` | l + 1 in the potential | FAILED: column | q differs from `zerilli_q` by 7.0×10⁻⁴ to 1.1×10⁻² |

Each mutation fails its own check and no other, and `flat` is a planted-truth run of the same code path (flat space has no
zero of Z with Re s > 0). The collector requires these outcomes. What is shared: `column` compares two integrations of
the same equation (Radau on the linear system, DOP853 on the Riccati one) that share `solve_star` (r_m, M_ext) and the
form of the Zerilli potential (the script has its own copy). The `flat` reference, scipy's modified spherical Bessel
function, shares only r_m. The counter's step rule is a copy of §5.1's, and
the planted zero in every box tests it. The sign with which a pole enters a winding was shown in #54 (a planted pole
counted −2) and was not repeated.

#### #124: poles of D in the l = 4, 5, 6 boxes

#124 counts the zeros of Z(r_m; s) in the nine boxes of §3a (#115: the stored box, s1 and n1 at l = 4, 5 and 6;
phantom ω = 0.7603), each with #115's planted pair. On rpi5 (Raspberry Pi 5,
aarch64, 4 cores shared with another project) on 2026-11-21, with the code of `wip/i124-zeros-l456` at b4c3a0a,
`zerilli_zeros.py` and `collect_i121.py` unmodified, at the script's default step rule as in #121. Each job is
`python -u zerilli_zeros.py <arguments> --cachedir=<raw folder>/cache`, with the 13 lines `job|arguments` that
`python collect_i121.py --jobs --only=3a` prints (nine boxes and four controls; raw `jobs.txt`), at `nice -n 10`, one
thread, at most two processes at a time. The four controls ran first, then the pilot `n1_l6` (R_z = 1.6×10⁶, 197 s in
the script), then the other eight boxes. Each job ran once and none died, and the 13 took 11 min in all. Logs and caches
are in negative-mass-raw, folder `I124-rpi5` at 68d643f. From `phase2/nonradial`, with the raw repository cloned beside
the checkout:

```
python collect_i121.py ../../../negative-mass-raw/I124-rpi5 --only=3a > i124_zeros_output.txt
```

The collector exits 0. Its `--mutate=box` (it expects `s1_l4` at Re s up to 0.2) and `--mutate=control` (it expects
`mut_lshift` to pass) each exit 1, as they should. The controls end as in the table above: `flat` all checks passed
(`flatq` at most 4.6×10⁻¹³), `mut_reverse` FAILED: plant (−2.0000 of +2), `mut_nonorm` FAILED: holo (1.3×10⁻² to
1.6×10⁻²), `mut_lshift` FAILED: column (7.0×10⁻⁴ to 1.1×10⁻²).

| job | l | Re s | Im s | zeros of Z | max step | min abs(Z) | planted found | winding of D | zeros of D |
|---|---|---|---|---|---|---|---|---|---|
| rep_l4 | 4 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | 0.27 | 5.59e-01 | +2.0000 of 2 | -0.0000 | 0 |
| s1_l4 | 4 | [0.03, 0.1] | abs(Im s) ≤ 2.5 | +0.0000 | 0.27 | 6.33e-01 | +2.0000 of 2 | +0.0000 | 0 |
| n1_l4 | 4 | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | +0.0000 | 0.22 | 2.10e+00 | +2.0000 of 2 | +0.0000 | 0 |
| rep_l5 | 5 | [0.1, 2] | abs(Im s) ≤ 2.5 | -0.0000 | 0.26 | 5.59e-01 | +2.0000 of 2 | -0.0000 | 0 |
| s1_l5 | 5 | [0.03, 0.1] | abs(Im s) ≤ 2.5 | +0.0000 | 0.27 | 6.34e-01 | +2.0000 of 2 | -0.0000 | 0 |
| n1_l5 | 5 | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | -0.0000 | 0.22 | 2.09e+00 | +2.0000 of 2 | -0.0000 | 0 |
| rep_l6 | 6 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | 0.34 | 5.59e-01 | +2.0000 of 2 | -0.0000 | 0 |
| s1_l6 | 6 | [0.03, 0.1] | abs(Im s) ≤ 2.5 | +0.0000 | 0.37 | 6.34e-01 | +2.0000 of 2 | +0.0000 | 0 |
| n1_l6 | 6 | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | +0.0000 | 0.11 | 1.94e+00 | +2.0000 of 2 | +0.0000 | 0 |

`rep` is §3a's stored box, the box of §5's `p7603_l2` at l = 4, 5, 6. The winding of D is §3a's
(`i115_boxes_output.txt`): guard 1 for eight boxes and `--guard=all` for `s1_l5`. Zeros of D = winding of D + zeros of
Z. `i124_zeros_output.txt` also has the evaluations, the times and the two check columns. Planted at 0.3 + 1.76i
(`rep`), 0.04 + 0.25i (`s1`) and 0.0005 + 0.12i (`n1`), each with its conjugate.

**All nine boxes have 0 zeros of Z(r_m; s)**, each with its planted pair returned. So the Zerilli column gives D no pole
in them, and the zeros of D equal the winding of D, 0, in each. No count is nonzero, so no box needs a corrected count
of zeros of D. Status: numerical evidence, with the conditions of the 21 boxes above: that column is D's only pole
source (§5.5). The counts of Z, and the windings of D except `s1_l5`'s, used guard 1; §5.4 (#126) recounts them with
the half-step test at every level, and none changes. The planted check
compares two counts on the same steps, so it does not see a 2π lost by Z's own steps (§3a). In the three n1 boxes the q
difference at Re s = 10⁻⁴ is reported and not checked (at most 9.3×10⁻¹⁰, `n1_l4`).

`count`, `holo`, `column` and `plant` pass in all nine. Largest max step 0.37 rad (`s1_l6`) on a count and 0.39 rad with
the planted pair (limit 0.45), worst Cauchy–Riemann defect 1.1×10⁻⁸ (`n1_l6`; limit 10⁻⁶), worst `column` difference
1.0×10⁻¹² (`s1_l4`, `s1_l5`; limit 10⁻⁶). `min abs(Z)` is 5.59×10⁻¹ in the three `rep` boxes, as in §5.

| acceptance criterion | verdict |
|---|---|
| each of the nine boxes returns its planted zero, with `count`, `holo`, `column` and `plant` passing | met in all nine |
| `collect_i121.py` exits 0 with the nine rows, and its mutations exit 1 | met: exit 0; `--mutate=box` and `--mutate=control` exit 1 |
| each nonzero count with its box and winding(D) + count(Z) | met; no count is nonzero, so winding(D) + count(Z) = 0 in all nine |

#### Not counted, and what this settles

- `rep_l2` is the first box of §5, `p7603_l2`. `n3_l2` and `n3_l3` cross Re s = 0, where k = √(s² + V) has its cut, so
  y₁(r_m) is not holomorphic there and a count would not be a count of zeros. #123 counts them with `--axis` (below).
  The 14 boxes of §5.1 counted here, with §5's seven, are the 21.
- #117's four boxes are counted in §5.1's #117 subsection, and #115's nine are #124's (above). Neither are the regions
  §5.1 lists under "Not covered" (Re s below the lower edges, above the upper ones).
- In the 21 boxes the zeros of Z(r_m; s) are absent, so the pole of D that #54's hypothesis and `decaying_plane` allow is
  absent, and each 0 of §5 and §5.1 counts zeros of D. That the Zerilli column is D's only pole source is the
  coordinator's reading of `polar_solver.py` (the regular plane is integrated from the centre, the scalar plane's
  spectral projector is holomorphic away from s = ±i(1 ± ω)). This run did not test that reading, so the "conditional on
  D having no pole" passages are unchanged here. The coordinator narrowed them to that reading at the merge, and §5.5
  (#122) checks it factor by factor.
- The exterior of the phantom star is nearly flat (r_m = 24.11, M_ext = −0.014580, so 2|M_ext|/r_m = 1.2×10⁻³), and flat
  space has no zero with Re s > 0. All counts are 0, as that suggests.

| acceptance criterion | verdict |
|---|---|
| each §5 box (7) and §5.1 box (17 less 3) returns its planted zero with `count`, `holo`, `column`, `plant` passing | met in all 21 |
| `--flat` passes `flat` and `flatq`; each mutation exits 1 on the stored box | met (table above) |
| counts reported with max step, min abs(Z), star and l from the logs by a collector whose `--mutate=` fails it; each nonzero count with its box and the corrected count of zeros of D | met; no count is nonzero, so zeros of D = 0 in all 21; both mutations exit 1 |
| #115's boxes, added as a second step | moved to #124 (above), met in all nine |

#### Boxes across Re s = 0 (#123)

#121 did not count `n3_l2` and `n3_l3`, the boxes of §5.1 that cross the axis (Re s ∈ [−10⁻³, 10⁻⁴], Im s ∈ [0.01, 0.23],
l = 2 and 3, phantom ω = 0.7603). Its function y₁ = Z exp Ψ has k = √(s² + V) on the principal branch, whose cut is on
Re s = 0, so y₁ is not holomorphic across the axis. #123 counts them with a new option, `zerilli_zeros.py --axis`.

Result: **both boxes have 0 zeros of Z(r_m; s)**, each with a planted zero returned on the left and on the right of the
axis. So the Zerilli column gives D no pole in them, and the zeros of D equal the winding of D, 0, in each. Status:
numerical evidence, with the condition of the 21 boxes above (that column is D's only pole source, §5.5).

**Method.** The counted function is y₁(r_m) = exp(s r_m*) Z(r_m; s), from y = (Z, f Z′) exp(s r*) with r* the tortoise
coordinate:

```
dy₁/dr = (y₂ + s y₁)/f,   dy₂/dr = ((s² + V) y₁ + s y₂)/f
```

It has no k, and exp(s r*) is entire and zero-free. Z is the Jost solution exp(−s r*) u, u = Σ aₙ r⁻ⁿ → 1. The aₙ follow
from f (f u′)′ − 2 s f u′ = V u by a recursion in 1/r that is singular in s only at s = 0, so the start
y(R_z) = (u, f u′ − s u) is holomorphic in the box. The series is summed until its terms are below 10⁻¹⁸. R_z is #121's
rule, r_m + (40 + 20 l)/|s|_min + 50 (8074.1 for l = 2, 10074.1 for l = 3), and y₁(r_m) does not depend on it. The
Cauchy–Riemann points are at Re s = x₀/2, 0 and x₁/2, so that their stencils lie on both sides of the axis and the limit
h ≤ Re s/2 is dropped; the criterion stays 10⁻⁶.

Two choices differ from the issue's suggestion (start (1, −k_c), integration on the real axis), and I made them after
measuring (negative-mass-raw, `I123-rpi5/explore`). For Re s < 0 the other solution grows inward on the real axis by
exp(−2 Re s (R_z* − r_m*)), e¹⁶ at Re s = −0.001 for l = 2 and e³² with `--zfac=2`, so an error at R_z is amplified
without limit; a start with the principal k is on the wrong branch there (the result at −0.0005 + 0.12i is 9×10³ times
the wanted function). With the Jost start the real-axis Radau solve at −0.0005 + 0.12i agrees with the integration
below to 2.4×10⁻¹³, but the three solves made at that point (two Jost starts, one WKB start) took 229 s, against 3.7 s
for the same three at 10⁻⁴ + 0.12i. So `--axis` integrates along the ray r = r_m − i p, p from R_z − r_m
to 0. For Im s > 0 the wanted solution decays outward there as exp(−Im s p) and the other decays inward, so the
integration is stable over the whole box (1.4 to 2.7 s per point), and its result is the continuation of the decaying
solution of Re s > 0. The count, the checks and the step rule are #121's. The docstring of `zerilli_zeros.py` has the
details.

**What ran.** On rpi5 (Raspberry Pi 5, aarch64, 4 cores shared with another project) on 2026-11-21, with the code of
`wip/i123-axis` at ce624b9 (the option at c4206eb). Each job is `python -u zerilli_zeros.py <arguments> --cachedir=<raw
folder>/cache`, with the lines `job|arguments` of `python collect_i121.py --jobs` (raw `jobs.txt`; the 25 of #121 were not
rerun), at `nice -n 10`, one thread, at most two processes at a time, a second one only while the 15-minute load was
below 4. The pilot `n3_l2` ran alone (246 s in the script, exit 0), then the other 11 jobs. 12 jobs took 50.8 min of
process time in 36.7 min of wall time. Each ran once and none died; the three that exit 1 are the mutations below. Logs and
caches are in negative-mass-raw, folder `I123-rpi5` at e43001a, with #121's logs copied beside them (sources and
sha256 in `sources.sha256`). From `phase2/nonradial`, with the raw repository cloned beside the checkout:

```
python collect_i121.py ../../../negative-mass-raw/I123-rpi5 > i121_zeros_output.txt
```

The collector exits 0, and the 21 rows of #121 regenerate byte for byte. It also exits 1 on a missing log, a tag that does
not match the arguments, a missing planted zero on either side, or a zfac difference above 10⁻⁶ (its `--mutate=axis`,
`oneside` and `zfac`, with `box` and `control`, each exit 1; logs `collect_mutate_*.log`).

| job | l | R_z | planted s₀ | zeros of Z | planted found | max step | min abs(Z) | evaluations | count (s) | winding of D | zeros of D |
|---|---|---|---|---|---|---|---|---|---|---|---|
| n3_l2 | 2 | 8074.1 | −0.0005 + 0.12i | +0.0000 | +1.0000 of 1 | 0.39 | 1.05 | 90 | 131 | +0.0000 | 0 |
| n3_l2_pr | 2 | 8074.1 | 0.00005 + 0.12i | +0.0000 | +1.0000 of 1 | 0.39 | 1.05 | 148 | 0 | +0.0000 | 0 |
| n3_l2_z2 | 2 | 16074.1 | −0.0005 + 0.12i | +0.0000 | +1.0000 of 1 | 0.39 | 1.05 | 90 | 165 | +0.0000 | 0 |
| n3_l3 | 3 | 10074.1 | −0.0005 + 0.12i | −0.0000 | +1.0000 of 1 | 0.39 | 1.11 | 98 | 195 | −0.0000 | 0 |
| n3_l3_pr | 3 | 10074.1 | 0.00005 + 0.12i | −0.0000 | +1.0000 of 1 | 0.39 | 1.11 | 156 | 0 | −0.0000 | 0 |
| n3_l3_z2 | 3 | 20074.1 | −0.0005 + 0.12i | +0.0000 | +1.0000 of 1 | 0.39 | 1.11 | 98 | 195 | −0.0000 | 0 |

The upper box (Im s ∈ [0.01, 0.23]) stands for its mirror, as in §5.1. `_pr` is the planted zero on the right of the axis;
it shares its cache with the left plant, so it shows 148 evaluations (90 counted before, 58 new for the plant's own
steps) and 0 s. `_z2` is `--zfac=2`. All six run at `--guard=all`. The winding of D is §5.1's (`i5_boxes_output.txt`), which
§5.4 recounts at `--guard=all` without change, and zeros of D = winding of D + zeros of Z. `min abs(Z)` is that of
y₁(r_m) = exp(s r_m*) Z, which is not the WKB-scaled function of the rows above. The count's largest step is 0.39 rad
(0.33 to 0.40 with the planted zero, limit 0.45), and the worst Cauchy–Riemann defect is 4.1×10⁻¹⁰ (`n3_l2_z2`; limit
10⁻⁶).

**Controls.** Each fails where it must, and each passes where it must:

| job | what | ends | result |
|---|---|---|---|
| `n3_l2_flat`, `n3_l3_flat` | `--flat` in the n3 box | all checks passed | counts 0, with no root of the reverse Bessel polynomial in the box (nearest −0.0622 ± 0.0359i for l = 2, −0.0763 ± 0.0728i and −0.0963 for l = 3); `flatq` at most 8.6×10⁻¹³ at Re s = 10⁻⁴ and −0.001 |
| `ax_flatzero` | `--flat`, box Re s ∈ [−0.09, −0.03], Im s ∈ [0.01, 0.07], planted −0.045 + 0.055i | all checks passed | counts +1.0000, the root −0.0622 + 0.0359i; `flatq` at most 6.3×10⁻¹³ |
| `ax_flatzero_rev` | the same, contour clockwise | FAILED: flat, plant | counted −1.0000 of +1 |
| `n3_l2_mprincipal`, `n3_l3_mprincipal` | `--mutate=principal`, start (1, −√(s² + V(R_z))) | FAILED: holo | Cauchy–Riemann defect 1.5×10⁻⁵ (1.2×10⁻⁵ for l = 3), then 1.2 and 1.2 at the two points whose stencil crosses the axis; `count` and `plant` still pass |

`ax_flatzero` is a planted-truth run: the zero is a root of a polynomial, at Re s < 0, and the same code counts it. The
principal-branch start makes the function jump across the axis (min abs(y₁) on the contour is 2.5×10⁻¹⁰ against 1.05), and
the count and the planted zero do not see it (they pass), so `holo` is the check that detects the defect of #121's form.
#121's three mutations are unchanged and were not rerun.

**Without `--axis`** the output of `zerilli_zeros.py` is that of 77b0541 (rerun on rpi5 from empty caches for `flat` and
`h1_l2`: identical apart from the times). Against #121's logs the RESULT lines agree except the sign of an exact zero
(`h1_l2` +0.0000 here, −0.0000 there) and the times, and the Cauchy–Riemann defects and `q` differences differ at the
10⁻⁹ level, finite-difference noise on another platform.

**zfac 1 against 2.** The start radius is doubled (8074.1 → 16074.1 for l = 2, 10074.1 → 20074.1 for l = 3). The counts are
equal, and y₁(r_m) at the three contour points of the column check differs by at most:

| box | zeros of Z, zfac 1 | zeros of Z, zfac 2 | largest relative difference of y₁ | at s |
|---|---|---|---|---|
| n3_l2 | +0.0000 | +0.0000 | 1.8×10⁻¹³ | −0.001 + 0.23i |
| n3_l3 | −0.0000 | +0.0000 | 1.4×10⁻¹³ | 0.0001 + 0.05422i |

y₁(r_m) does not depend on R_z, so these are numerical differences; the collector requires at most 10⁻⁶ (a limit I set
before the run).

**`q` against `zerilli_q`** (D's own Zerilli column; not a criterion, since the script checks `column` only where
Re s ≥ 0.03 and none of the three points is there; the line prints 0 of 3). `q` = y₂/y₁ at r_m, relative difference:

| box | points | s | side | difference |
|---|---|---|---|---|
| n3_l2 | contour | 0.0001 + 0.05422i | Re s > 0 | 6.0×10⁻⁷ |
| n3_l2 | contour | −0.001 + 0.23i | Re s < 0 | 7.3×10⁻⁶ |
| n3_l2 | contour | −0.001 + 0.05312i | Re s < 0 | 1.8×10⁻⁵ |
| n3_l2 | Cauchy–Riemann | −0.0005 + 0.076i | Re s < 0 | 8.8×10⁻⁶ |
| n3_l2 | Cauchy–Riemann | 0 + 0.142i | Re s = 0 | 4.0×10⁻⁶ |
| n3_l2 | Cauchy–Riemann | 0.00005 + 0.208i | Re s > 0 | 3.3×10⁻⁶ |
| n3_l3 | contour | 0.0001 + 0.05422i | Re s > 0 | 2.8×10⁻⁸ |
| n3_l3 | contour | −0.001 + 0.23i | Re s < 0 | 1.0×10⁻⁵ |
| n3_l3 | contour | −0.001 + 0.05312i | Re s < 0 | 1.8×10⁻⁶ |
| n3_l3 | Cauchy–Riemann | −0.0005 + 0.076i | Re s < 0 | 1.4×10⁻⁶ |
| n3_l3 | Cauchy–Riemann | 0 + 0.142i | Re s = 0 | 4.0×10⁻⁶ |
| n3_l3 | Cauchy–Riemann | 0.00005 + 0.208i | Re s > 0 | 3.7×10⁻⁶ |

At nearly the same |s| = 0.054 the difference on the left of the axis is 30 times (l = 2) and 64 times (l = 3) that on
the right. `zerilli_q` starts at a per-point radius from the WKB data, so its error from the omitted terms, damped
inward for Re s > 0, is amplified inward for Re s < 0 by exp(2 |Re s| (R_start − r_m)): 22 for R_start = 1578 at
−0.001 + 0.053i. Along the two edges of n3 (12 points, `column_scan_l2.log` and `_l3.log`) the difference is
at most 3.5×10⁻⁶ at Re s = 10⁻⁴ (3.4×10⁻⁶ for l = 3) and at most 8.9×10⁻⁵ at Re s = −0.001 (l = 2, Im s = 0.01; 1.3×10⁻⁵
for l = 3), so D's column on the left is about 10⁻⁴ off at worst, of the order that factor (9×10⁶ there) times the
right-hand difference (5.1×10⁻¹¹ at 10⁻⁴ + 0.01i) gives. This bears on D's own Zerilli column on the left of the axis,
which §5.1's n3 windings used (#122); I did not change it.

**Independent check of `q` at Re s < 0.** Raw `cross_continue_i123.py` reads the `q` of an `--axis` log at its six points
(on both sides of the axis) and compares it with #121's form of the solution, which is stable for Re s > 0, continued
across the axis by a Taylor series (order 8, fitted on 12 points of a circle of radius 0.0018 centred at 0.002 + i Im s).
It shares only the potential, r_m and M_ext with `--axis`: neither the Jost start nor the rotated path. The largest
relative difference is 2.7×10⁻⁹ for `n3_l2` (at −0.001 + 0.05312i) and 2.2×10⁻¹⁰ for `n3_l3` (at −0.001 + 0.23i), against
a limit of 10⁻⁸ set before the run. On the flat log the same run gives 2.7×10⁻⁹ at the same point while the log's `q` is
within 9.4×10⁻¹³ of the exact flat value, so the 2.7×10⁻⁹ is the continuation's own error, not the integration's. With
the series evaluated 10⁻⁴ off (`--mutate=shift`) the check fails (2.7×10⁻³). The first runs of the script missed the point
on the axis (a parser error, kept as `*_first_missed_axis_point.log`) and were rerun.

**Corrected count of zeros of D** = winding(D) + zeros of Z: 0 + 0 = 0 for `n3_l2` and for `n3_l3`. No count is nonzero.

**What this settles and what it does not.** Z(r_m; s) has no zero in either box, on either side of the axis, so the
Zerilli column adds no pole of D there, and the windings of §5.1's `n3_l2` and `n3_l3` (0, with a planted zero at
−0.0005 + 0.12i) count zeros of D, under the condition that column is D's only pole source (§5.5). The count is of the
analytic continuation of the decaying solution across Re s = 0 at Im s ∈ [0.01, 0.23]; the lower mirror has the same count
by Z(conj s) = conj Z(s), the continuation passing through Re s > 0. A zero of Z with Re s < 0 outside the box is not
excluded (flat space has them at Re s ≈ −0.06). The sliver |Im s| < 0.01 is not covered (§5.1). The checks are
numerical: no interval enclosure and no proof.

| acceptance criterion | verdict |
|---|---|
| a planted zero on each side of the axis is returned in each box | **Met**: −0.0005 + 0.12i and 0.00005 + 0.12i, +1.0000 of 1 in `n3_l2` and `n3_l3` |
| Cauchy–Riemann check at points on both sides of Re s = 0 passes at 10⁻⁶; with the principal-branch k it fails | **Met**: worst 4.1×10⁻¹⁰; `--mutate=principal` fails `holo` (1.2) in both boxes |
| `--flat` counts 0 in each box | **Met**: `n3_l2_flat`, `n3_l3_flat` |
| the count in each box and the corrected count of zeros of D | **Met**: 0 in each, and winding(D) + zeros of Z = 0 |

### 5.4 Guard at every level (#126)

The half-step test of §5.1 and §5.2 runs at bisection depths 0 and 1 only (`--guard=1`, the default). A step accepted at
depth 2 or deeper can lose 2π, and the planted check does not see it when it hits D and D times the factor alike: in
#115, s1 at l = 5 lost 2π in its planted count that way (§3a). #126 counts again at `--guard=all`, which applies the
test at every depth, each box whose count the sections above report at guard 1. Each runs on a copy of its stored cache
with the planted zero of its first run. That is 31 windings of D with `winding2.py` and 30 counts of the zeros of Z(r_m;
s) with `zerilli_zeros.py`, which gained `--guard` for this: §5's 7 boxes, §5.1's other 16 (its `rep` box at l = 2 is
§5's `p7603_l2`), #115's 8 other than s1 at l = 5 (#115 counted that one at guard all), and the 21 and 9 boxes of §5.3.

On rpi5 (Raspberry Pi 5, aarch64, 4 cores shared with another project) on 2026-11-21, 00:49 to 09:35, with the code of
`wip/i126-guard-all` at 37ec4d3. Each job is `python -u <script> <arguments> --cachedir=<raw folder>/cache`, with the 66
lines `job|script|arguments|source folder` that `python collect_i126.py --jobs` prints (raw `jobs.txt`), at nice level
10, one thread, at most two processes at a time and the second only while the 15-minute load was under 4. The five
controls ran first, then the 30 Z boxes, then the pilot `w_s2_l2` (the slowest at guard 1, 3.7 h in #5; stop condition 8
h; it took 2 h 22 min), then the other 30 W boxes, longest expected first. Each job ran once and none died. The second
process started `w_s3_l3` at 02:03 with the 15-minute load at about 6, which the rule does not allow. From 02:37 every
start of a second process waited for the load (`launch_I126_gated.sh`, the waits are in `queue.log`), and the load
stayed above 4 for most of the run, so the W boxes ran mostly on one process. This changes the times and no result. Logs
and caches are in negative-mass-raw, folder `I126-rpi5` at e4323fa. From `phase2/nonradial`, with the raw repository
cloned beside the checkout:

```
python collect_i126.py ../../../negative-mass-raw/I126-rpi5 > i126_guard_output.txt
```

The collector exits 0. **No count changed: all 61 counts at `--guard=all` equal their guard-1 values, 31 windings of D
and 30 counts of zeros of Z, all 0.** Each planted zero is returned (+2.0000 of 2 in the symmetric boxes, +1.0000 of 1
in the upper boxes) and every check passes (for the Z boxes `count`, `holo`, `column` and `plant`). The guard-1 values
are read from `i54_recount_output.txt`, `i5_boxes_output.txt`, `i115_boxes_output.txt`, `i121_zeros_output.txt` and
`i124_zeros_output.txt`. Its five controls end as required: `z_flat` all checks passed, `z_mut_reverse` FAILED: plant,
`z_mut_nonorm` FAILED: holo, `z_mut_lshift` FAILED: column and `w_mut_reverse` PLANTED FAIL. Its mutations
`--mutate=noplant`, `shift` (a count of 0.5), `changed` (a guard-1 value 2 higher), `noguard` (the job list without
`--guard=all`) and `control` (`z_mut_lshift` expected to pass) each exit 1. A copy of the raw folder without
`w_s1_l2.log` also exits 1 ("no result"); that was run once by hand and is not a mutation of the script. The 31 windings
of D:

| job | section | Re s | Im s | guard 1 | guard all | max step (rad) | min abs(D) | evaluations | planted found |
|---|---|---|---|---|---|---|---|---|---|
| p7603_l2 | 5 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | +0.0000 | 0.39 | 6.70e-04 | 135 | +2.0000 of 2 |
| p7603_l3 | 5 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | +0.0000 | 0.40 | 1.68e-03 | 195 | +2.0000 of 2 |
| c075_l2 | 5 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | +0.0000 | 0.40 | 6.88e-04 | 111 | +2.0000 of 2 |
| p070_l2 | 5 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | +0.0000 | 0.40 | 7.45e-04 | 134 | +2.0000 of 2 |
| p070_l3 | 5 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | +0.0000 | 0.40 | 1.86e-03 | 209 | +2.0000 of 2 |
| p085_l2 | 5 | [0.1, 2] | abs(Im s) ≤ 2.5 | -0.0000 | +0.0000 | 0.40 | 4.93e-04 | 129 | +2.0000 of 2 |
| p085_l3 | 5 | [0.1, 2] | abs(Im s) ≤ 2.5 | +0.0000 | +0.0000 | 0.40 | 1.23e-03 | 205 | +2.0000 of 2 |
| n1_l2 | 5.1 | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | -0.0000 | +0.0000 | 0.33 | 2.82e-03 | 84 | +2.0000 of 2 |
| n1_l3 | 5.1 | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | +0.0000 | +0.0000 | 0.35 | 3.16e-02 | 112 | +2.0000 of 2 |
| n3_l2 | 5.1 | [-0.001, 0.0001] | [0.01, 0.23] | +0.0000 | +0.0000 | 0.39 | 2.75e-03 | 116 | +1.0000 of 1 |
| n3_l3 | 5.1 | [-0.001, 0.0001] | [0.01, 0.23] | -0.0000 | +0.0000 | 0.38 | 3.66e-02 | 124 | +1.0000 of 1 |
| s1_l2 | 5.1 | [0.03, 0.1] | abs(Im s) ≤ 2.5 | +0.0000 | +0.0000 | 0.40 | 1.00e-03 | 180 | +2.0000 of 2 |
| s1_l3 | 5.1 | [0.03, 0.1] | abs(Im s) ≤ 2.5 | +0.0000 | +0.0000 | 0.40 | 2.49e-03 | 233 | +2.0000 of 2 |
| h1_l2 | 5.1 | [0.1, 2] | [2.5, 6] | +0.0000 | +0.0000 | 0.38 | 1.05e-04 | 108 | +1.0000 of 1 |
| h1_l3 | 5.1 | [0.1, 2] | [2.5, 6] | -0.0000 | +0.0000 | 0.38 | 2.63e-04 | 126 | +1.0000 of 1 |
| t1_r40 | 5.1 | [0.02, 0.1] | [0.15, 0.33] | +0.0000 | +0.0000 | 0.38 | 3.73e-03 | 94 | +1.0000 of 1 |
| t1_r80 | 5.1 | [0.02, 0.1] | [0.15, 0.33] | -0.0000 | +0.0000 | 0.39 | 3.73e-03 | 94 | +1.0000 of 1 |
| t2_r40 | 5.1 | [0.02, 0.1] | [1.68, 1.84] | +0.0000 | +0.0000 | 0.40 | 9.30e-04 | 94 | +1.0000 of 1 |
| t2_r80 | 5.1 | [0.02, 0.1] | [1.68, 1.84] | +0.0000 | +0.0000 | 0.35 | 9.30e-04 | 94 | +1.0000 of 1 |
| s2_l2 | 5.1 | [0.01, 0.03] | abs(Im s) ≤ 2.5 | +0.0000 | +0.0000 | 0.40 | 5.30e-04 | 227 | +2.0000 of 2 |
| h2_l2 | 5.1 | [0.1, 2] | [6, 10] | -0.0000 | +0.0000 | 0.31 | 2.57e-05 | 92 | +1.0000 of 1 |
| s3_l2 | 5.1 | [0.03, 0.1] | [2.5, 6] | -0.0000 | +0.0000 | 0.39 | 9.90e-05 | 106 | +1.0000 of 1 |
| s3_l3 | 5.1 | [0.03, 0.1] | [2.5, 6] | -0.0000 | +0.0000 | 0.38 | 2.71e-04 | 106 | +1.0000 of 1 |
| rep_l4 | 3a | [0.1, 2] | abs(Im s) ≤ 2.5 | -0.0000 | +0.0000 | 0.40 | 3.02e-03 | 261 | +2.0000 of 2 |
| s1_l4 | 3a | [0.03, 0.1] | abs(Im s) ≤ 2.5 | +0.0000 | +0.0000 | 0.40 | 4.43e-03 | 289 | +2.0000 of 2 |
| n1_l4 | 3a | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | +0.0000 | +0.0000 | 0.40 | 4.17e-02 | 131 | +2.0000 of 2 |
| rep_l5 | 3a | [0.1, 2] | abs(Im s) ≤ 2.5 | -0.0000 | +0.0000 | 0.40 | 4.70e-03 | 327 | +2.0000 of 2 |
| n1_l5 | 3a | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | -0.0000 | +0.0000 | 0.40 | 7.17e-02 | 128 | +2.0000 of 2 |
| rep_l6 | 3a | [0.1, 2] | abs(Im s) ≤ 2.5 | -0.0000 | +0.0000 | 0.40 | 6.71e-03 | 407 | +2.0000 of 2 |
| s1_l6 | 3a | [0.03, 0.1] | abs(Im s) ≤ 2.5 | +0.0000 | +0.0000 | 0.40 | 8.69e-03 | 433 | +2.0000 of 2 |
| n1_l6 | 3a | [0.0001, 0.1] | abs(Im s) ≤ 0.23 | +0.0000 | +0.0000 | 0.39 | 1.35e-01 | 144 | +2.0000 of 2 |

The 30 Z boxes (`i126_guard_output.txt` has every row) all count 0 at guard 1 and at guard all, with max step 0.00 to
0.38 rad and min abs(Z) 3.13e-01 to 2.10e+00. The largest max step of a D count is 0.40 rad (17 boxes), of a Z count
0.38 rad (`p085_l2`), against the 0.45 limit of the planted check. The smallest min abs(D) is 2.57e-05 (`h2_l2`). The
planted zero is returned in all 30 Z boxes: 21 with +2.0000 of 2 and 9 with +1.0000 of 1.

Evaluations of D are the points each run used: 84 (`n1_l2`) to 433 (`s1_l6`) in the 31 windings. The caches grew from
6922 to 9003 points in 65 files, 2081 new (`summary.txt` in the raw folder has each). The largest increases are
`rep_l6` (205), `s1_l6` (198) and `rep_l5` (159); s1 at l = 5 has no
job here, and its 387 points are #115's. Walls, from `queue.log`: the W boxes took 11.71 h in all, the longest `s2_l2`
(2 h 22 min), `s1_l6` (1 h 41 min) and `s1_l4` (55 min); the Z boxes 17.5 min and the controls 3.2 min. All 66 jobs took
12.06 h of process time in 8.77 h of wall time.

#### The replay

`diagnose_plant.py` replays the step rule of `winding2.py` on a stored cache, without evaluating D, and compares each
accepted step with the dense chain of every cached point on the edge. It was moved here from negative-mass-raw
`I115-rpi5`, where #115 used it, and gained `--expect` and `--mutate`. From `phase2/nonradial`, with the raw repository
cloned beside the checkout:

```
R=../../../negative-mass-raw/I115-rpi5/cache; B=${R}_before_guard; T=winding2_cache_eps-1_w0.7603_l5_0.03_0.1_2.5.pkl
{ echo "# s1 at l = 5, the 212-point cache before the guard-all rerun, guard 1 (the run that failed its planted check)" &&
  python diagnose_plant.py $B/$T 0.03 0.1 2.5 0.04,0.25 --expect=-0.3826,+5.9006 &&
  echo "# s1 at l = 5, the 387-point cache after the guard-all rerun, guard all" &&
  python diagnose_plant.py $R/$T 0.03 0.1 2.5 0.04,0.25 --guard=all --expect=none; } > i115_replay_output.txt
```

- On the 212-point cache of `cache_before_guard`, at guard 1, one step differs from the chain, in D times the factor: t
  in [0.90000, 0.92500] of the edge Re s = 0.03 (s = 0.03000+0.25000i to 0.03000+0.18750i), measured −0.3826 rad against
  the chain's +5.9006 rad over 14 cached steps (largest 0.730), a loss of one turn. The planted difference is then 0
  instead of +2, the failure of #115's first run. D alone has no such step. The check passes, so the replay reproduces
  −0.3826 against +5.9006.
- On the 387-point cache after the rerun, at guard all, no step differs and the chain's largest increment is 0.211 rad;
  D times the factor counts +2.0000 on the steps and on the chain.
- Mutations, each on the 212-point cache and each exit 1 (raw `diagnose_*.log`): `--mutate=sparse` keeps only the
  accepted steps' end points in the chain, so the chain loses the turn too and no step differs; `--mutate=wrongexpect`
  shifts the expected phases by 10⁻³; `--expect=none` at guard 1 finds the one step; `--guard=all` asks for a point the
  cache lacks, (0.1+0.00390625i), and stops, as the cache is not the one of the guard-all run.
- §3a calls the stretch "the 14 cached points". The 14 are cached steps between 15 points, 13 of them inside the
  stretch, the 13 removed in the coordinator's check (#126's description). §3a is unchanged here.

**Dense-chain replay of every box (#128).** `replay_i126.py` runs this replay on the cache of each of the 61 boxes of
this section, at `--guard=all --expect=none` and with the planted zero of its first run. `diagnose_plant.py` gained
`--y0=Y0` for it: with Y0 > 0 the contour is `winding2.py`'s upper box (four edges counterclockwise from (x1, Y0), no
doubling) and the planted factor is (s − s0) alone; Y0 = 0, the default, is the symmetric box as before, and
`i115_replay_output.txt` still regenerates byte for byte with the command above. From `phase2/nonradial`, with the raw
repository cloned beside the checkout:

```
python replay_i126.py ../../../negative-mass-raw/I126-rpi5 > i126_replay_output.txt
```

It exits 0, in about 11 s on rpi5, and prints one row per box (the counted function, the accepted steps of F and of F times
the factor, the chain's points, the replayed counts, the planted zeros found, the chain's largest increments).
**In all 61 boxes (31 of D, 30 of Z) no accepted step differs from the dense chain, in F or in F times the planted factor,
and none fails.** The chains hold 69 to 445 cached points; the largest increment of any chain is 0.314 rad (`z_n1_l2`, Z
times the factor; 0.313 rad for D, `w_s2_l2`), against the limit of 1 rad. The replayed count of F is 0 in all 61 boxes.
The planted zeros are found: +2.0000 in the 41 symmetric boxes and +1.0000 in the 20 upper boxes (11 of D, 9 of Z), which
go through `--y0`. A second script, `cross_replay_i126.py` in raw `I128-rpi5`, orders the cached points of each contour by
arclength, shares only the cache files, `jobs.txt` and the contour's geometry with the replay, and reproduces every row's
chain points and largest increments to 0.0006 rad (its `--mutate=thin` exits 1).

Controls, all in negative-mass-raw `I128-rpi5` at dc61d9d (code at 716ce01):
- Before the real caches, `synthetic_y0.py` ran `--y0` on caches of known functions: the upper box counts +1 and +2 with
  the planted zero inside and +1 and +1 below it, the symmetric box +2 and +4, as the zeros and poles inside it say. A
  function that hides 40 turns per edge from every step of depth 2 or less exits 1, with the dense chain's 40 turns
  printed on all four edges. Its `--mutate=wrongcount` exits 1.
- A `--y0` box replayed without it fails: the symmetric contour starts at (x1, 0), which an upper box's cache lacks, and
  the replay stops on that point. `replay_i126.py --mutate=noy0` stops all 20 upper boxes this way and passes the 41
  symmetric ones.
- `replay_i126.py --mutate=sparse` (passes `--mutate=sparse` to `diagnose_plant.py` on `w_h1_l2`), `noy0`, `misplant` (the
  expected planted zeros of `w_h1_l2` from a zero outside the box), `nocache` and `nojobs` (jobs.txt without `w_p7603_l2`)
  each exit 1 (raw `replay_mutate_*.log`). `diagnose_plant.py --mutate=sparse --expect=none` exits 0 on a normal box: a
  chain cut to the steps' end points agrees with the steps. The replay therefore also requires each chain to hold every
  cached point on its contour, counted from the pickle, and that is what `--mutate=sparse` fires (48 chain points against
  126 on the contour).
- `collect_i126.py --mutate=nolog` reads no log for `w_s1_l2` and exits 1; the unmutated collector still regenerates
  `i126_guard_output.txt` byte for byte.

Status: numerical evidence. The chain shares D and its cache with the count: it tests the step rule on the points the
runs cached, and is no test of D.

| acceptance criterion | verdict |
|---|---|
| in each of the 61 boxes no accepted step at `--guard=all` differs from the dense chain, in F or in F times the factor; each chain's largest increment is below 1 rad | met in 61 of 61; largest increment 0.314 rad |
| the upper boxes go through `--y0`; a `--y0` box replayed without it fails | met: 20 upper boxes with `--y0`; without it each stops on the missing point (x1, 0) |
| `--mutate=sparse` on one box and `collect_i126.py --mutate=nolog` each exit 1 | met for `replay_i126.py --mutate=sparse` (through the chain-points check) and `collect_i126.py --mutate=nolog`; `diagnose_plant.py --mutate=sparse` alone exits 0 under `--expect=none` |

#### What this settles

Status: numerical evidence that the 61 counts of §5, §5.1, §3a and §5.3 do not depend on the depth at which the
half-step test stops. The test compares a step with its two halves, so a loss that survives both halves is not detected
at any depth, and the planted check counts on the same D and the same steps. The dense chain of the replay uses every
cached point and so is a test of the step rule apart from the half-step test. #126 ran it on s1 at l = 5, and #128 on
all 61 boxes (above), where no step differs from it. It sees a loss only between cached points, so a turn hidden between
two neighbouring cached points stays undetected. The Z boxes use the same rule. The other conditions of §5.3 stand as
they are (that the Zerilli column is D's only pole source, §5.5).

| acceptance criterion | verdict |
|---|---|
| at `--guard=all` each box returns its planted zero, with max step, min abs(D), evaluations; a winding that differs from guard 1 is named with its box | met in all 61 (31 D, 30 Z); none differs, so none is named |
| a collector exits 1 on a failed planted check, a missing log or a non-integer winding, and its mutations exit 1 | met: exit 0 on the real logs; `--mutate=noplant`, `shift`, `changed`, `noguard` and `control` each exit 1; a missing log exits 1 (by hand) |
| the replay reproduces −0.3826 against +5.9006 on `cache_before_guard/` and exits 1 under a mutation | met: PASS on the 212-point cache; `--mutate=sparse`, `--mutate=wrongexpect`, `--expect=none` and `--guard=all` on it each exit 1 |
| the Z boxes of §5.3 recounted at `--guard=all` (added by the coordinator) | met in all 30; none changed |

### 5.5 D's factors and its poles (#122)

§5.3 finds no zero of the decaying Zerilli solution in the counted boxes, so the Zerilli column gives D no pole there.
That makes each winding of 0 a count of zeros of D only if the Zerilli column is D's only pole source. #122 lists every
factor of `polar_solver.PolarProblem.D` along its code path and checks the others with the new script `d_factors.py`
(`d_factors_output.txt`). Five parts check the factors, the scalar plane's spectral projector, its eigenvector
normalisation, a pole planted in two factors, and the n3 boxes with the continued Zerilli column. A sixth part, `mean`,
was added after the first failed. Status: numerical evidence.

**Result.** The checks on the individual factors are met in the 61 boxes. The numerical check that D as computed is
holomorphic, `harm`, is not met.
- Met: the two selected eigenvalues of the scalar block stay apart from the other two, with a smallest relative
  separation of 0.122 in all 61 boxes (`s2_l2`; limit 10⁻³) and a box around s = i(1 − ω) failing at 6×10⁻⁹; D does not
  depend on the eigenvector normalisation (at most 1.0×10⁻¹³ relative, limit 10⁻¹⁰); a pole planted in the Zerilli
  column and one planted in the scalar plane's first column each change the winding of `n1_l2` by −2.0000 through D's
  code path; with q of `zerilli_zeros.py --axis` in place of `zerilli_q` the windings of `n3_l2` and `n3_l3` stay 0 and
  the planted zero is returned.
- **Not met: `harm`.** arg D is not harmonic to the criterion fixed before the runs in 5 of the 6 boxes checked. At
  winding2's rtol 10⁻⁸ the four boxes with Re s ≥ 0.1 exceed it by up to 129 times, and `n3_l2` by 840 to 1.3×10⁵. In
  `n3_l2` it is the `zerilli_q` factor (D with the continued q passes), as #123 found. In the other four the cause is
  not established; below says what was shown and what is a reading.

#### The factors

In order along `D`, with `D = det[Y, Z] · p_Y · p_Z · exp(−i Im Θ)`, Θ = (k₊ + k₋)(R_s − r_m) + (k₊ + k₋ + s) r_m, k± =
√(1 + (s ∓ iω)²) (principal branch) and p the products of the QR phases. `d_factors.py` prints this table.

| # | factor | holomorphic in a box? | zero-free? | checked by |
|---|---|---|---|---|
| 1 | interior start data Y0 (`regular_plane`) | constant, real, s-independent | 3 independent columns | structure |
| 2 | J(r; s²), interior and vacuum | polynomial in s² (entire); singularities in r only | - | `jpoly` |
| 3 | interior flow Φ(r_m, r0; s), DOP853 in chunks | entire: linear ODE with entire coefficients | Φ invertible, so the plane keeps rank 3 | `harm` (numerical) |
| 4 | QR of the regular plane: Q, phase of R's diagonal | D = det[Y, Z] / ∏ abs(R_ii), a positive divisor; arg D = arg det | abs(R_ii) > 0 | `pole` (the factor enters as c / abs(c)), structure |
| 5 | Moncrief map T(r_m) | constant in s (poles in r at 2 M_ext and −3 M_ext/λ, not at r_m) | det T ≠ 0 | `moncrief` |
| 6 | Zerilli column zv = (1, q/f), q = f Z′/Z | meromorphic; poles at the zeros of Z(r_m; s) (the pole source; none in the 61 boxes: #121, #123, #124) | first entry 1 | `axis` (D_axis against D), `harm` at n3_l2 |
| 7 | scalar block Js = Jvac(R_s, s)[2:, 2:], R_s fixed per box | entire (factor 2) | - | `jpoly` |
| 8 | eig(Js), the two eigenvalues of smallest real part | eigenvectors carry arbitrary normalisations and the choice is by real part: neither is holomorphic by itself | - | `norm`, `order` |
| 9 | projector P = V_sel V⁻¹_sel | holomorphic while the selected pair stays apart from the other two (Riesz projector); independent of the eigenvector normalisation | - | `gap`, `norm` |
| 10 | start columns P[:, 0], P[:, 2] | holomorphic | independent iff P e₀ and P e₂ are | `cols` |
| 11 | scalar plane flow R_s → r_m and its QR phases | as 3 and 4 | as 3 | `harm`, `pole` (scalar column) |
| 12 | decaying plane Z = [m, S] and its QR | m = T zv and S as above; positive divisor as 4 | m ≠ 0 (T invertible, zv[0] = 1); S rank 2 (σ_min(S) = 1.0000 at the `norm` points, reported) | `pole` (Zerilli column), `norm` |
| 13 | det[Y, Z] times the QR phases | holomorphic matching determinant divided by a positive function | zeros are the modes | `harm` |
| 14 | winding2's exp(−i Im Θ) | a pure phase, not holomorphic; its argument −Im Θ is that of exp(−Θ), holomorphic and zero-free where k± is, i.e. off the cut on Re s = 0, abs(Im s ∓ ω) ≥ 1 | never zero | `harm`, `gap` |

Rows 8 to 10 are the factors the issue names (the projector and the normalisation). Row 10 was not in its list: D
contains the columns P[:, 0] and P[:, 2], and if they became dependent D would have a zero that the decaying plane does
not.

**What makes arg D the argument of a holomorphic function.** The QRs divide each plane by a positive number, since Y =
QR gives det[Y, Z] = det[Q_Y, Q_Z] det R_Y det R_Z and D has det[Q_Y, Q_Z] times the phases of det R, so D = det[Y, Z] /
∏ abs(R_ii) for any QR. With numpy the diagonal of R is real, so these phases are signs, and a column's complex phase is
in Q, not in R (a pole planted as a column factor c enters D as c / abs(c), `pole` below). The unnormalised matrix is
holomorphic apart from the Zerilli column's poles. winding2's factor exp(−i Im Θ) is a pure phase, which no holomorphic
function except a constant can be; but D times it is (det[Y, Z] exp(−Θ)) times a positive function, because exp(−i Im Θ)
= exp(−Θ) · exp(Re Θ), so its argument is that of a holomorphic function. Where k± is holomorphic, arg D is therefore
harmonic away from zeros and poles of det[Y, Z] exp(−Θ), and `harm` tests that numerically.

#### What ran

On rpi5 (Raspberry Pi 5, aarch64, 4 cores shared with another project) on 2026-11-21, from 12:11 to 15:45, with the code
of `wip/i122-factors`: `d_factors.py` at 9f34847 until 13:12, then the working tree that became 2222a98 (it adds the
file suffix of non-default rtol, the two stencil Laplacians in the `harm` lines and the `mean` part, which those jobs
did not call), and 7ace97a for the last output run. Each job is `python -u
d_factors.py <arguments> --raw=<raw folder>` from `jobs.txt`, at `nice -n 10`, one thread, in a load-gated queue (a
second process only while the 15-minute load was below 4; it stayed above 4 for most of the run, so most jobs ran one at
a time). 38 jobs finished (17 exit 0 and 21 exit 1: 12 are mutations, which must fail, and 9 are failed checks), 4.09 h
of process time. D evaluations take 4 to 12 s (a stored box about 12 s). Logs, evaluation caches (one pickle per part
and box) and the queue scripts are in negative-mass-raw, folder `I122-rpi5` at 813692a. The box caches of `I126-rpi5`
and the `--axis` logs of `I123-rpi5` are read, never written. From `phase2/nonradial` with the raw repository cloned
beside the checkout, the output file regenerates from the caches in about 6 min:

```
python d_factors.py --part=factors,gap,norm,pole,axis,mean --raw=../../../negative-mass-raw/I122-rpi5 > d_factors_output.txt
```

It exits 1, because `harm` and `mean` fail. All criteria are in the script's docstring and were fixed before the runs,
except `mean`'s, which was written after `harm` failed and before any `mean` run.

#### `moncrief`, `jpoly`

`moncrief`: |det T| / (|col 1| |col 2|) at r_m is 0.071 or more in the 10 (star, l) of the boxes (limit 10⁻³), so T is
invertible. `jpoly`: J(r; s²) is linear in s² and even in s, to 6.8×10⁻¹⁴ relative at five radii in [r₀, r_m] and three
in vacuum for each (star, l) (limit 10⁻¹⁰); a planted J with a t² term gives 4.3×10⁻² and fails.

#### `harm`: arg D harmonic

At three interior points of six boxes, D as `winding2.py` evaluates it (rtol 10⁻⁸), the five-point Laplacian of φ = arg
D at h and h/2, Richardson-extrapolated, h = 10⁻³ of the box's longer side L. Criterion: |Δφ| L² ≤ 10⁻⁶ (1 + |∇φ| L).
The boxes are the stored box at l = 2 and 3 (phantom 0.7603), the canonical star's at l = 2, phantom 0.85 at l = 3,
`n1_l4` and the upper box `n3_l2`. In each, the planted phase of (s − a)(s − b) eˢ passes (at most 2×10⁻³ of the limit)
and that phase plus 10⁻³ (Re s)² / L² fails (2.0×10⁻³ against limits of 2×10⁻⁶ to 9×10⁻⁶).

| box | s | abs(∇φ) | abs(Δφ) L² | limit | ratio | verdict |
|---|---|---|---|---|---|---|
| `p7603_l2` | +0.575+0.75i | 3.77 | 9.06×10⁻⁵ | 1.04×10⁻⁵ | 8.71 | **exceeds** |
| `p7603_l2` | +1.05+1.5i | 2.85 | 2.24×10⁻⁶ | 8.12×10⁻⁶ | 0.276 | pass |
| `p7603_l2` | +1.525+2.25i | 2.3 | 4.20×10⁻⁵ | 6.76×10⁻⁶ | 6.21 | **exceeds** |
| `p7603_l3` | +0.575+0.75i | 5.92 | 9.02×10⁻⁵ | 1.58×10⁻⁵ | 5.71 | **exceeds** |
| `p7603_l3` | +1.05+1.5i | 4.3 | 6.10×10⁻⁶ | 1.17×10⁻⁵ | 0.521 | pass |
| `p7603_l3` | +1.525+2.25i | 3.36 | 3.63×10⁻⁵ | 9.40×10⁻⁶ | 3.86 | **exceeds** |
| `c075_l2` | +0.575+0.75i | 3.17 | 1.08×10⁻⁴ | 8.93×10⁻⁶ | 12.1 | **exceeds** |
| `c075_l2` | +1.05+1.5i | 2.19 | 1.16×10⁻⁵ | 6.48×10⁻⁶ | 1.79 | **exceeds** |
| `c075_l2` | +1.525+2.25i | 1.68 | 6.70×10⁻⁴ | 5.21×10⁻⁶ | 129 | **exceeds** |
| `p085_l3` | +0.575+0.75i | 6 | 1.28×10⁻⁵ | 1.60×10⁻⁵ | 0.8 | pass |
| `p085_l3` | +1.05+1.5i | 4.12 | 2.07×10⁻⁵ | 1.13×10⁻⁵ | 1.83 | **exceeds** |
| `p085_l3` | +1.525+2.25i | 3.25 | 1.04×10⁻³ | 9.11×10⁻⁶ | 114 | **exceeds** |
| `n1_l4` | +0.025075+0.069i | 23.7 | 4.46×10⁻⁸ | 6.45×10⁻⁶ | 0.00691 | pass |
| `n1_l4` | +0.05005+0.138i | 24.6 | 6.06×10⁻⁸ | 6.65×10⁻⁶ | 0.00911 | pass |
| `n1_l4` | +0.075025+0.207i | 34.5 | 6.48×10⁻⁷ | 8.94×10⁻⁶ | 0.0725 | pass |
| `n3_l2` | −0.0005+0.076i | 19.2 | 6.87×10⁻¹ | 5.22×10⁻⁶ | 1.32×10⁵ | **exceeds** |
| `n3_l2` | +0+0.142i | 7.96 | 2.11×10⁻¹ | 2.75×10⁻⁶ | 7.67×10⁴ | **exceeds** |
| `n3_l2` | +5e-05+0.208i | 19.8 | 4.51×10⁻³ | 5.37×10⁻⁶ | 840 | **exceeds** |

**FAIL in 5 of 6 boxes.**

- The four boxes with Re s ≥ 0.1 exceed the limit at 9 of 12 points, by up to 129 times (`c075_l2` at 1.525 + 2.25i,
  with 6.7×10⁻⁴ against 5.2×10⁻⁶). `p7603_l2` and `p7603_l3` give nearly equal values (9.1×10⁻⁵ and 9.0×10⁻⁵ at 0.575 +
  0.75i).
- In `c075_l2` at 1.525 + 2.25i the Laplacian times L² is −1.9×10⁻³ at h and +2.8×10⁻⁵ at h/2, so it does not scale as
  h².
- `n1_l4` (Re s ≥ 0.025) passes, at 7 % of the limit.
- `n3_l2` exceeds it by 840 to 1.3×10⁵, with the same value at h and h/2 (smooth, not noise). With q of `--axis` in
  place of `zerilli_q` the same three points give 4.2×10⁻⁸, 4.5×10⁻⁸ and 2.2×10⁻⁷ (`axis` below), so the non-harmonic
  factor is `zerilli_q`.

*Shown:* the failures above, and that the criterion's planted pair and its mutation `nonholo` (D times exp(i 10⁻³ (Re
s)² / L²) fails `n1_l4`, 310 times the limit) work. In the four boxes with Re s ≥ 0.1, a rerun of `p7603_l2` at rtol
10⁻¹⁰ (`harm_diag_p7603_l2_rtol1e-10.log`) gives 6.8×10⁻⁷ (passes), 1.4×10⁻⁵ and 1.7×10⁻⁵ (exceed 8.1×10⁻⁶ and 6.8×10⁻⁶)
at the three points; at rtol 10⁻¹² (`harm_diag_p7603_l2_rtol1e-12.log`) 1.8×10⁻⁸ and 5.4×10⁻⁶ (pass) and 1.0×10⁻⁵
(exceeds 6.8×10⁻⁶). The worst ratio to the limit is 8.7 at rtol 10⁻⁸, 2.6 at 10⁻¹⁰ and 1.5 at 10⁻¹².

*Reading, not tested:* that the failures at Re s ≥ 0.1 are integrator noise (of order 10⁻¹¹ in arg D at rtol 10⁻⁸, which
the stencil multiplies by 4 (L/h)² ≈ 4×10⁶ in abs(Δφ) L²), sometimes with a discrete jump in the stencil (`nchunks` of
`_integrate_plane` depends on |s|, and 1.525 + 2.25i lies within 0.0014 of a threshold). The first point's values fall
with the tolerance (9.1×10⁻⁵, 6.8×10⁻⁷ and 1.8×10⁻⁸ at rtol 10⁻⁸, 10⁻¹⁰ and 10⁻¹²), the third point's fall more slowly
(4.2×10⁻⁵, 1.7×10⁻⁵, 1.0×10⁻⁵) and stay above the limit, and the second point's do not fall monotonically (2.2×10⁻⁶,
1.4×10⁻⁵, 5.4×10⁻⁶). So rtol alone does not remove the failure at 10⁻¹². `zerilli_q` integrates at a fixed rtol of 10⁻¹¹
(`polar_solver.py`) that `rtol` does not set, which is a candidate floor: also a reading, not tested.

#### `mean`: the mean value property

Added after `harm` failed, as a better-conditioned test of the same property, with its criterion fixed before its runs.
φ at the middle interior point of each box against the mean of φ over 16 points of a circle of radius ρ = 0.1 L (5 h in
`n3_l2`, which crosses Re s = 0 and where `zerilli_q` degrades outside the box); 4 (mean − φ_centre) / ρ² stands for the
Laplacian in the same inequality, and for a harmonic φ it is exact with no 1/h² amplification of noise. The planted pair
gives at most 1.5×10⁻¹¹ (passes) and 2.0×10⁻³ (fails) in every box; `nonholo` fails `n1_l4` (300 times).

| box | s | ρ | abs(∇φ) | abs(Δφ) L² | limit | ratio | verdict |
|---|---|---|---|---|---|---|---|
| `p7603_l2` | +1.05+1.5i | 0.25 | 2.85 | 6.02×10⁻⁶ | 8.12×10⁻⁶ | 0.741 | pass |
| `p7603_l3` | +1.05+1.5i | 0.25 | 4.3 | 5.74×10⁻⁶ | 1.17×10⁻⁵ | 0.489 | pass |
| `c075_l2` | +1.05+1.5i | 0.25 | 2.19 | 4.20×10⁻⁶ | 6.48×10⁻⁶ | 0.649 | pass |
| `p085_l3` | +1.05+1.5i | 0.25 | 4.12 | 2.22×10⁻⁵ | 1.13×10⁻⁵ | 1.97 | **exceeds** |
| `n1_l4` | +0.05005+0.138i | 0.023 | 24.6 | 1.78×10⁻⁸ | 6.65×10⁻⁶ | 0.00268 | pass |
| `n3_l2` | +0+0.142i | 0.0011 | 7.96 | 2.38×10⁻¹ | 2.75×10⁻⁶ | 8.63×10⁴ | **exceeds** |

Two boxes fail, as computed: `p085_l3` (ratio 1.97) and `n3_l2` (8.6×10⁴, `zerilli_q` again). The other three pass, at
ratios of 0.49 to 0.74 that are close to 1, with ρ = 0.25 (L = 2.5). *Reading, not tested:* that a 16-point trapezoid on
a circle of radius 0.25 has an aliasing error of about (ρ/d)¹⁶ for a singularity of det[Y, Z] at distance d, which at d
≈ 0.7 would give about 4×10⁻⁵ in this quantity, as large as the values above. The test is a smaller ρ at the same point,
`--rho=`, not run before this report (the runs are in `listF_after_report.txt` in the raw folder). #129 tested it:
the failure is integrator error, not aliasing (below).

#### `gap`: the projector's eigenvalue separation

Over each box's contour (401 points per edge), the points of its cache and a 20 × 20 grid, polished by a local search
inside the box: the relative separation of the two selected eigenvalues of the vacuum scalar block at R_s from the other
two (`gap`), the gap between their real parts (`order`: the selection by real part is stable), and σ_min / σ_max of the
columns 0 and 2 of the projector (`cols`).

| check | smallest over the 61 boxes | box | limit | verdict |
|---|---|---|---|---|
| `gap` | 0.1222 at s = 0.01 + 1.7603i | `s2_l2` (Re s ∈ [0.01, 0.03]) | 10⁻³ | pass |
| `order` | 6.77×10⁻³ | `s2_l2` | 10⁻³ | pass |
| `cols` | 0.3436 | `p085_l2`, `p085_l3` | 10⁻³ | pass |

The branch points s = ±i(1 ± ω) are at ±0.2397i, ±1.7603i (ω = 0.7603), ±0.25i, ±1.75i (0.75), ±0.3i, ±1.7i (0.70) and
±0.15i, ±1.85i (0.85). The nearest a counted box comes is 0.010 (`s2_l2` and the `n1`, `n3` boxes), 0.020 (the threshold
boxes) and 0.030 (`s1`); the stored boxes are 0.100 from the axis, `s3` and `h1` 0.74. The 61 rows (the contour, cache
and grid minima, R_s, the polished minimum and where it is, and the distance to the branch points) are in
`d_factors_output.txt`. Planted: boxes of ±0.01 around s = i(1 − ω) and i(1 + ω) give polished separations of 6.4×10⁻⁹
and 1.3×10⁻⁸ (at 0.2397i and 1.7603i) and an order gap below 10⁻¹⁸, so they fail, and their 20 × 20 grids alone
(9.0×10⁻² and 3.3×10⁻²) would not: the polish is what finds the point. A projector with P e₀ = P e₂ gives a column ratio
of 0 and fails `cols`. Mutations: `branch` (every box moved onto i(1 − ω)) fails all three (1.4×10⁻⁹, 0, 6×10⁻¹⁷) and
`cols` fails `cols`.

#### `norm`: the eigenvector normalisation

D with each eigenvector of `np.linalg.eig` in `scalar_plane` multiplied by a random nonzero complex factor (modulus
log-uniform in [0.1, 10], uniform phase) at each s, at five points in each of the six boxes: the largest relative
difference from D is 1.0×10⁻¹³ (`p7603_l2`; the others 4.9×10⁻¹⁴, 6.1×10⁻¹⁴, 2.9×10⁻¹⁴, 1.5×10⁻¹⁴ and 5.5×10⁻¹⁵; limit
10⁻¹⁰), and σ_min of the scalar plane's 4 × 2 block is 1.0000. Mutation `noproj` (the selected eigenvectors themselves
as start columns, `--boxes=n1_l4,c075_l2`) gives 1.99 and fails.

#### `pole`: a pole planted through D's code path

In `n1_l2` (phantom 0.7603, l = 2, Re s ∈ [10⁻⁴, 0.1], |Im s| ≤ 0.23, 0 at `--guard=all` in #126, planted zero 0.0005 +
0.12i) the factor 1/((s − s₀)(s − s̄₀)), s₀ = 0.05 + 0.1i, is multiplied onto the Zerilli column (zv in
`decaying_plane`) and, in a second run, onto the scalar plane's first start column, in copies of those two methods on
the script's own instance of `PolarProblem` (the files are unchanged). The box is recounted with a copy of
`winding2.py`'s step rule at `--guard=all`.

| run | winding | change against #126's 0 | max step | points, time | D_pole / D against c / abs(c) |
|---|---|---|---|---|---|
| Zerilli column | −2.0000 | −2.0000 | 0.40 rad | 81, 348 s | 1.1×10⁻¹⁵ at 63 cached points |
| scalar plane, first column | −2.0000 | −2.0000 | 0.40 rad | 81, 395 s | 4.3×10⁻¹⁴ at 63 cached points |

The copied counter on #126's cache reproduces its winding (0), max step (0.33 rad), min abs(D) (2.82×10⁻³) and planted
pair (+2.0000; 71 distinct points of the 84 its cache held), and the copied methods with no factor reproduce the cached
D at three points exactly. The pole enters D as the factor's phase c / abs(c), to rounding, so the count follows.
Mutations, each exit 1: `reverse` (contour clockwise) gives +2.0000; `nophase` (the phase factor of the planted column's
QR dropped, as specified) does not lose the pole, since numpy's R has a real diagonal, but removes a sign, so D
jumps by ±1 and the identity deviates by 2.00; the count had not finished after 300 evaluations; `nogauge` (that
column's phase fixed by its largest component, which does lose the pole) gives +8.0000 with a step of 3.14 rad and a
deviation of 2.00 after 293 points.

#### `axis`: n3 with the continued Zerilli column

`n3_l2` and `n3_l3` (Re s ∈ [−10⁻³, 10⁻⁴], Im s ∈ [0.01, 0.23]) with q from `zerilli_zeros.py --axis` (its function
copied; it reproduces the six q values of the I123 logs exactly) in place of `zerilli_q`, at every cached contour point
(144 and 152 points, 294 and 310 evaluations). The dense-chain winding is that of the cached points on the four edges.

| box | edge | points | largest abs(D_axis − D) / abs(D) | largest abs(D_axis − D) | of min abs(D) |
|---|---|---|---|---|---|
| `n3_l2` | right Re s = 10⁻⁴ | 53 | 4.1×10⁻⁶ | 8.0×10⁻⁷ | 2.9×10⁻⁴ |
| `n3_l2` | top Im s = 0.23 | 21 | 2.2×10⁻⁶ | 6.1×10⁻⁹ | 2.2×10⁻⁶ |
| `n3_l2` | left Re s = −10⁻³ | 53 | 7.9×10⁻⁵ | 4.2×10⁻⁶ | 1.5×10⁻³ |
| `n3_l2` | bottom Im s = 0.01 | 21 | 7.9×10⁻⁵ | 4.2×10⁻⁶ | 1.5×10⁻³ |
| `n3_l3` | right | 57 | 3.7×10⁻⁶ | 6.8×10⁻⁷ | 1.9×10⁻⁵ |
| `n3_l3` | top | 21 | 2.4×10⁻⁶ | 1.1×10⁻⁷ | 3.1×10⁻⁶ |
| `n3_l3` | left | 57 | 2.5×10⁻⁵ | 3.2×10⁻⁶ | 8.9×10⁻⁵ |
| `n3_l3` | bottom | 21 | 1.1×10⁻⁷ | 8.3×10⁻⁹ | 2.3×10⁻⁷ |

min abs(D) on the contour is 2.75×10⁻³ (`n3_l2`) and 3.61×10⁻² (`n3_l3`). The difference is largest on the left edge,
where `zerilli_q` is off by up to 8.9×10⁻⁵ (§5.3), and the relative difference is below 10⁻⁴ at every cached point, so
the argument of D_axis / D stays within ±π/2 on the chain and D_axis and D have the same winding. Winding of the cached
D on the dense chain: 0 (largest increment 0.210 and 0.201 rad); of D_axis: 0 (0.209 and 0.133 rad); with the planted
zero −0.0005 + 0.12i, D_axis: +1.0000 in both, so the criterion (winding 0, planted zero returned) is met. Harmonic at
the three points of `harm`: D_axis passes at 4 % of the limit in `n3_l2` and 16 % in `n3_l3`, while D (with `zerilli_q`)
fails by 1.3×10⁵ and 4.1×10⁴. Mutations: `reverse` gives −1.0000 for the planted zero and fails `axis_plant`; `wrongq`
(q times 1 + 10⁻⁴) fails `axis_q` (10⁻⁴ against 10⁻¹⁰).

#### What this settles and what it does not

- *Shown (numerical evidence):* the projector's separation, the order of the selection and the independence of its start
  columns hold in the 61 boxes with margins of 7 to 340 times the limits; D does not depend on the eigenvector
  normalisation; a pole in either column of the decaying plane enters D exactly as its phase and changes the winding by
  −2; the Zerilli column's non-holomorphy in `n3` (via `zerilli_q`) does not change the windings, and the continued
  column gives windings of 0 with the planted zero returned.
- *Not shown:* that arg D is harmonic to 10⁻⁶ in the boxes with Re s ≥ 0.1 at rtol 10⁻⁸ (the check fails there), so a
  non-holomorphic part of D at that level is not excluded by this check. That such a part cannot change a count is a
  reading, not tested here: it would have to move arg D by about the 0.4 rad of a counted step, and the failures are of
  the order 10⁻⁵ to 10⁻³ in abs(Δφ) L².
- *Reading, not tested:* that those failures are integrator noise and that the `mean` failure at `p085_l3` is aliasing
  (above). #129 tested both: integrator error at all 12 points, and not aliasing at `p085_l3` (below).
- *Coordinator's review (#122's closing comment):* a fresh rerun of `c075_l2` on another machine fails `harm` again,
  with the same value at h at the worst point and a different one at h/2. arg D at rtol 10⁻⁸ differs from arg D at
  10⁻¹⁰ by far more than the 10⁻¹¹ of the reading above, so that estimate does not hold. A mean-value test on 64 points
  at the worst point fails narrowly in `c075_l2` at rtol 10⁻⁸ and passes there at 10⁻¹⁰, and passes in `p085_l3` at
  10⁻⁸: at those points the failure is the integrator's error. #129 tests every point.

**Integrator error or a non-holomorphic part (#129).** Status: numerical evidence. `harm` failed at 9 of 12 points of
the four boxes with Re s ≥ 0.1, and the paragraphs above read that as integrator error without testing it.
`harm_circles.py` tests φ = arg D by the mean value property at those 12 points: 64 points of circles of radius ρ = 0.05
L and 0.025 L (L = 2.5), 4 mean(φ_k − φ₀) / ρ² for Δφ, `harm`'s criterion |Δφ| L² ≤ 10⁻⁶ (1 + |∇φ| L) unchanged. The
mean does not amplify a point error by 1/h²: an error ε in a point moves |Δφ| L² by at most 4 ε L² / ρ², that is 1600 ε
at 0.05 L and 6400 ε at 0.025 L. D is winding2's D at three settings. S8: rtol 10⁻⁸, `nchunks` as `polar_solver` chooses
it at each point, `zerilli_q` at its own 10⁻¹¹. S10: rtol 10⁻¹⁰, `nchunks` fixed across each circle at the centre's
value, `zerilli_q` at 10⁻¹². S12: 10⁻¹², fixed, 10⁻¹³. The fixed `nchunks` and `zerilli_q`'s tolerance are set on a
subclass in the script; `polar_solver.py` is unchanged. At every point and radius a planted holomorphic phase, (s − a)(s
− c) e²ˢ, must pass at 10⁻³ of the limit, and the same times exp(i 10⁻³ (Re s)² / L²) must fail, through the same circle
code.

What ran, from `phase2/nonradial` (the raw folder is `I129-win` of negative-mass-raw, beside this checkout; each
`<setting>` run is 3 points × 129 evaluations, cached after every point):

```
python harm_circles.py --box=<box> --setting=<setting> --raw=../../../negative-mass-raw/I129-win
python harm_circles.py --box=c075_l2 --setting=S12 --raw=../../../negative-mass-raw/I129-win --mutate=nonholo
python harm_circles.py --box=c075_l2 --setting=S12 --raw=../../../negative-mass-raw/I129-win --mutate=half
python harm_circles.py --collect --raw=../../../negative-mass-raw/I129-win > harm_circles_output.txt
```

The twelve runs (`<box>` in `p7603_l2`, `p7603_l3`, `c075_l2`, `p085_l3`; `<setting>` in S12, S10, S8) exited 0 at S12
and S10 (eight runs) and 1 at S8 (four runs: a failed check, a result). Both mutations, on the cached S12 data of
`c075_l2` with no new evaluation, exit 1, and `--collect` exits 0 (`RESULT PASS`, every control passing, no cache
missing, D passing at S12). Time per evaluation on one thread of `win`: 17.6 to 45.2 s at S12 and 7.0 to 19.2 s at S8.
The ratio |Δφ| L² / limit (pass at or below 1) at both radii, from `harm_circles_output.txt` (pivoted by
`results_table.py` in the raw folder):

| box | point | S8, 0.05 L | S8, 0.025 L | S10, 0.05 L | S10, 0.025 L | S12, 0.05 L | S12, 0.025 L | first passing |
|---|---|---|---|---|---|---|---|---|
| `p7603_l2` | +0.575+0.75i | 8.47 **FAIL** | 8.35 **FAIL** | 0.0767 | 0.102 | 0.000771 | 0.00134 | S10 |
| `p7603_l2` | +1.05+1.5i | 0.799 | 0.692 | 0.0086 | 0.0091 | 5.01×10⁻⁶ | 0.000671 | S8 |
| `p7603_l2` | +1.525+2.25i | 1.43 **FAIL** | 1.32 **FAIL** | 0.0131 | 0.0139 | 0.000422 | 0.00114 | S10 |
| `p7603_l3` | +0.575+0.75i | 5.68 **FAIL** | 5.61 **FAIL** | 0.0511 | 0.0678 | 0.000193 | 0.000507 | S10 |
| `p7603_l3` | +1.05+1.5i | 0.502 | 0.405 | 0.00607 | 0.00675 | 4.32×10⁻⁵ | 0.000551 | S8 |
| `p7603_l3` | +1.525+2.25i | 1.07 **FAIL** | 1.08 **FAIL** | 0.00971 | 0.0103 | 0.000181 | 0.000683 | S10 |
| `c075_l2` | +0.575+0.75i | 11.8 **FAIL** | 11.4 **FAIL** | 0.0926 | 0.0918 | 0.000718 | 0.000446 | S10 |
| `c075_l2` | +1.05+1.5i | 0.749 | 0.823 | 0.00939 | 0.00981 | 0.00022 | 0.00154 | S8 |
| `c075_l2` | +1.525+2.25i | 1.78 **FAIL** | 1.54 **FAIL** | 0.0176 | 0.0203 | 0.000664 | 0.00253 | S10 |
| `p085_l3` | +0.575+0.75i | 1.08 **FAIL** | 1.13 **FAIL** | 0.0173 | 0.0151 | 0.000206 | 0.000175 | S10 |
| `p085_l3` | +1.05+1.5i | 2.11 **FAIL** | 2.16 **FAIL** | 0.00358 | 0.0357 | 0.000378 | 0.000761 | S10 |
| `p085_l3` | +1.525+2.25i | 0.261 | 0.0579 | 0.00405 | 0.00336 | 0.000232 | 0.000692 | S8 |

- *Every point passes at S12*, at ratios 5.0×10⁻⁶ to 2.5×10⁻³ (|Δφ| L² at most 1.4×10⁻⁸ against limits of 5.2×10⁻⁶ to
  1.6×10⁻⁵). No point fails at rtol 10⁻¹² with `nchunks` fixed, so none of the 12 is, by this test, a part of D's phase
  that the integrator does not explain.
- *Which setting removed each failure.* S8 fails 8 of the 12 points (16 of 24 rows, ratios 1.07 to 11.8): point 0 and 2
  of `p7603_l2`, `p7603_l3` and `c075_l2`, and points 0 and 1 of `p085_l3`. S10 passes all eight (3.6×10⁻³ to 0.102), so
  S10 is the first setting that passes there. The other four points pass at S8 (0.058 to 0.82) and S8 is their first.
  S10 changes rtol, the fixed `nchunks` and `zerilli_q`'s tolerance at once, so these runs do not say which of the three
  removes a failure.
- *The values fall with rtol.* The largest |Δφ| L² is 1.06×10⁻⁴, 1.07×10⁻⁶ and 1.4×10⁻⁸ at S8, S10 and S12 (99 and 77
  times smaller per two decades of rtol), and the smallest point error that explains it under the bound above is
  6.6×10⁻⁸, 5.2×10⁻¹⁰ and 5.0×10⁻¹², 6.6, 5.2 and 5.0 times rtol. A part of arg D that did not depend on rtol would not
  fall. What remains at S12 is below 1.4×10⁻⁸ in |Δφ| L² on these discs, 7×10⁻⁶ of the planted non-harmonic part's
  2.0×10⁻³.
- *Controls.* In all 72 rows (12 points, 2 radii, 3 settings) the planted holomorphic phase gives at most 4.3×10⁻¹³
  (passes) and the planted non-harmonic one 2.0×10⁻³ (fails). `nonholo` fails D at all six rows, at 224 to 384 times the
  limit, and `half` fails the planted holomorphic control (240 and 480, where at most 10⁻³ of its limit passes) and D
  (2.4×10⁷ to 6.8×10⁷ times).
- *Against the text above.* The `mean` failure at `p085_l3` (1.97, 16 points at ρ = 0.25) was read as aliasing, about
  (ρ/d)ᴺ. It is not: with 64 points at ρ = 0.125 and 0.0625 that estimate is about 10⁻⁴⁸ at d ≈ 0.7, yet `p085_l3` point
  1 fails at S8 by 2.11 and 2.16, and passes at S10 (3.6×10⁻³, 3.6×10⁻²). The review's mean-value result holds:
  `c075_l2` point 2 fails narrowly at S8 (1.78 and 1.54) and passes at S10 (0.018, 0.020), and `p085_l3` point 2 passes
  at S8 (0.26, 0.058); the other two points of `p085_l3` fail at S8. `harm` and the circles do not agree at every point:
  at point 2 of `c075_l2` and `p085_l3` `harm` gave 6.7×10⁻⁴ and 1.0×10⁻³ and the circles at S8 (0.05 L) give 9.3×10⁻⁶
  and 2.4×10⁻⁶, while at point 0 of `p7603_l2` and `c075_l2` they are alike (9.1×10⁻⁵ and 1.1×10⁻⁴ against 8.8×10⁻⁵ and
  1.1×10⁻⁴).
- *Not shown.* That `harm`'s stencil (h = 10⁻³ L) passes at S12, which was not run, so the reading that rtol alone does
  not remove its failure at 10⁻¹² is neither confirmed nor contradicted. That arg D is harmonic at any other point than
  these 12, at other radii, or in `n3_l2` (where `zerilli_q` is the cause). The `harm` rows above stay as recorded.

| #129 acceptance criterion | verdict |
|---|---|
| at every point the planted holomorphic phase passes and the planted non-harmonic one fails through the same code; `nonholo` fails | met: 72 of 72 rows (planted holomorphic at most 4.3×10⁻¹³, planted non-harmonic 2.0×10⁻³); `nonholo` and `half` exit 1 |
| each point's verdict at each rtol and radius as computed; if every point passes at S12, which setting removed each failure | met: the table; all 12 pass at S12; S10 removes the 8 failures of S8, and the runs do not separate its three changes |
| a point that still fails at rtol 10⁻¹² with fixed `nchunks` is reported as a part of D's phase that the integrator does not explain | none: 0 of 12 points fail at S12 (24 of 24 rows pass) |

- *Coordinator's review (#129's closing comment):* on the Mac, `--collect` on the `I129-win` caches gives the table above
  (one S12 entry differs in the fourth digit) and `--mutate=nonholo` exits 1. A fresh S12 run of `c075_l2` point 2 agrees
  with the laptop's arg D to 3.3×10⁻¹² and passes at 7.9×10⁻⁴ and 2.9×10⁻³ of the limit. An independent check
  (negative-mass-raw `I129-mac`, `ind129b.py`), with its own transform and a two-radius estimate that removes the ρ⁴
  term, reproduces every ratio of the table and passes at S12 (at most 3.2×10⁻³). Over the 1548 cached points arg D at
  S8 differs from S12 by at most 6.8×10⁻⁷ rad (S10: 6.3×10⁻⁹), far below the 0.4 rad of a counted step. |D| is
  phase-normalised by `PolarProblem.D`, so log|D| is not harmonic (10⁵ to 10⁶ times the limit) and only arg D is
  tested.

- The checks are numerical, one phase of each plane, and no interval enclosure or proof. They do not look beyond the 61
  boxes or outside the part of Re s = 0 that the n3 boxes cover.

| acceptance criterion | verdict |
|---|---|
| a table of D's factors, each with its argument or check of holomorphy | met (the table above, in `d_factors_output.txt`) |
| the smallest separation of the selected eigenvalues per box, a planted box around s = i(1 − ω) failing | met: 0.122 at worst, planted 6×10⁻⁹ |
| D unchanged to 10⁻¹⁰ under random eigenvector factors; a mutation skipping the normalisation failing | met: 1.0×10⁻¹³; `noproj` 1.99 |
| a pole planted in one factor changes the winding (by −2 in the symmetric box) | met: −2.0000 in the Zerilli column and in the scalar plane's first column; `reverse`, `nophase`, `nogauge` fail |
| `harm`: φ = arg D harmonic, planted pair passing and failing | **not met**: 5 of 6 boxes exceed the limit (worst ratios 5.7 to 129 in the four boxes at rtol 10⁻⁸, 1.3×10⁵ in `n3_l2`); the planted pair and `nonholo` behave |
| `axis`: winding of D_axis equals 0 and the planted zero is returned | met in `n3_l2` and `n3_l3` |

## 6. Conclusions

1. **Polar, real growth rates.** For the phantom campaign star (ω = 0.7603, M = −0.01458) there is
   **no unstable polar mode with γ ∈ [10⁻³, 300] (l = 2) or γ ∈ [10⁻³, 1000] (l = 3)**. The search used
   100 log-spaced points, and D keeps the same sign under all convergence variants. The point-mass DG
   rates (137 and 686) lie inside these ranges. The canonical star (ω = 0.75) gives the same result
   (expected: it is on the stable branch). The phantom stars ω = 0.70 and 0.85 give the same result for
   l = 2 and 3 (#7).
2. **Polar, complex growth rates.** Argument-principle count (validated on a known QNM): winding **0** in
   Re s ∈ [0.1, 2], |Im s| ≤ 2.5, for phantom l = 2, phantom l = 3 and canonical l = 2, and for ω = 0.70
   and 0.85 (l = 2, 3). #7's −2 for ω = 0.85, l = 3 was an error of the first step rule, and the recount gives 0
   (§5.2). §5.3 (#121) finds no zero of the decaying Zerilli solution at r_m in these boxes, so the 0 counts zeros
   of D, conditional on that column being D's only pole source (§5.5). At l = 4, 5, 6 the stored box, s1 and n1
   count 0 too (§3a, #115), and #124 (§5.3) finds no zero of the Zerilli solution there, under the same condition.
   Every count used the step rule at guard 1, and the recount at every level gives the same (§5.4, #126).
3. **Axial: stable for all l ≥ 2.** V_ax > 0 everywhere; positivity is exact once V_ax > 0 is checked on
   the grid (r²V_ax/(α² l(l+1)) ≥ 1 − 10⁻⁶).
4. **κ_eff.** κ_eff − κ_c = s²/(q_int − g) exactly, and for the phantom star q_int > g, so the matter
   supplies **κ_eff > κ_c at every γ > 0** (l = 2 tabulated for γ ∈ [10⁻⁴, 150]; l = 2, 3 at γ → 0 via
   c_l). The static limit is degenerate (κ_eff → κ_c as s²). The sign is fixed by a Love-type
   coefficient, c_l < 0 for the phantom stars (c_l > 0 canonical). The large-R requirement κ < κ_c ≈
   −(l+1)/R is **not** met. The phantom matter acts as a stiff (κ > 0 for γ ≳ 2.5×10⁻⁴) surface, so the
   DG negative-mass-vacuum instability is not triggered.

## 7. Caveats

- The real-γ scans of the other stars, and of the campaign star at l = 4, 5, 6, resolve sign changes on a 100-point
  log grid (30 points for canonical l = 3). A pair of nearby roots (an even number between grid points) would be missed
  there. The κ_eff picture (q_int − q_ext ≈ 2γ, far from 0) makes this implausible, but it is not a proof. For the
  campaign star at l = 2, 3 the grid was refined to 397 points (step ratio 1.032 for l = 2, 1.036 for l = 3) with no
  change in the sign count, so a pair closer than about 3 % in γ is excluded there only by the smoothness of ln |D|
  (§5.1).
- Complex s, campaign star, l = 2 (l = 3 where it differs): counted are Re s ∈ [0.1, 2] for |Im s| ≤ 20 (10 for l = 3,
  #117), Re s ∈ [2, 5] for |Im s| ≤ 10 (l = 2 and 3, #117), Re s ∈ [0.03, 0.1] for |Im s| ≤ 6, Re s ∈ [0.01, 0.03] for
  |Im s| ≤ 2.5 (l = 2 only), and Re s ∈ [−10⁻³, 0.1] for |Im s| ≤ 0.23 (§5.1). The rest is listed under "Not covered"
  there. The other stars have Re s ≥ 0.1, |Im s| ≤ 2.5 only. At l = 4, 5, 6 the campaign star has the stored box, s1 and
  n1 only (§3a).
- Every winding except the rerun of s1 at l = 5 used guard 1, the half-step test at bisection depths 0 and 1 only. In
  #115 a planted count lost 2π at depth 2. The planted check cannot see a loss that hits D and D times the factor
  alike (§3a). §5.4 (#126) recounts every guard-1 count of §3a, §5 and §5.1 to §5.3 with `--guard=all`, and none
  changes. A loss inside a step whose two halves agree is excluded by neither guard. The dense chain of every cached
  point agrees with every accepted step in all 61 boxes (§5.4, #128), so such a loss would have to fall between two
  neighbouring cached points.
- The windings of §5 were computed with a step rule that can lose 2π in one step (§5.1). One loss shifts a symmetric-box
  count by ±2. They were recounted with the current rule in §5.2. Only ω = 0.85, l = 3 changed, from −2 to 0, because the
  old rule lost one 2π there.
- l = 2, 3 for every star, and l = 4, 5, 6 only for the ω = 0.7603 phantom star (§3a, #4), with the real-γ scans,
  convergence variants and κ_eff, windings in three boxes only (the stored box, s1 and n1; eight of the nine returned
  their planted pair at guard 1, s1 at l = 5 did only with `--guard=all`; §3a, #115; no pole of the Zerilli column
  in them, §5.3, #124) and no c_l. The ω = 0.7603, 0.70 and 0.85 stars have the l = 2, 3 real-γ scans, convergence
  variants and windings (0.70 and 0.85 in #7). c_l and axial positivity are also checked at ω = 0.85, and axial
  positivity at 0.70 and 0.90. At l = 5 and 6 the rtol 1e-12 variant is missing at the top γ
  (§3a).
- The windings count zeros minus poles of D. The −2 that suggested a net pole at ω = 0.85, l = 3 was the old step rule
  (§5.2), so no box shows a net pole. §5.3 (#121) finds no pole from the Zerilli column in 21 boxes. That it is D's only
  pole source is checked factor by factor in §5.5 (#122), whose test that arg D is harmonic fails at Re s ≥ 0.1 at rtol
  10⁻⁸ (#129). #123 finds none in `n3_l2` and `n3_l3` either, across Re s = 0; there D's own column, from `zerilli_q`,
  is up to 8.9×10⁻⁵ off on the left of the axis (§5.3).
- The equations for (tt) and the angular trace were verified numerically (relative residual ≤ 7×10⁻⁷,
  `verify_polar.py`) and are now *checked symbolically* (#8, `verify_polar_exact.py`, §1 Checks): on solutions of
  y′ = J y both residuals are exactly 0, as are those of (tθ), (tr), (rθ), (rr), KG_R and KG_I. The check takes the
  separated equations of `polar_eqs.pkl` as given. The coordinator's independent re-derivation from G_ab − 8πεT_ab
  agrees exactly at 6 rational points (§1 Checks).
- Weak-field stars only (|M|/R ~ 3×10⁻³). Nothing here says what happens for compact phantom stars.
- Linear theory, Regge–Wheeler gauge, m = 0 (sufficient by spherical symmetry).

## 8. Files

| file | role |
|---|---|
| `lincommon.py` | shared sympy linearisation of G_ab − 8πεT_ab and KG (ε symbolic) |
| `derive_nonradial.py` (+`_output.txt`) | polar linearisation and separation → `polar_eqs.pkl` (the exact check of 3 equations once attributed to it is not in it; #8) |
| `reduce_polar.py` (+`_output.txt`) | joint reduction → `polar_matrix_generated.py`; vacuum Zerilli/KG checks |
| `verify_polar.py` (+`_output.txt`) | numerical check of all 8 equations |
| `verify_polar_exact.py` (+`_output.txt`) | exact check of all 8 equations on y′ = J y (#8): *checked symbolically* in twelve variables (gcd-free, 4 s) and exactly 0 at 20 rational points; `--mutate=coef\|bg\|eq`, `--selftest`, `--crosscheck`, `--sensitivity`, `--method`, `--engine`, `--free`; exits 1 when a check fails or a mutation is detected |
| `derive_axial.py`, `axial.py` (+outputs) | axial master equation, positivity of V_ax |
| `background.py` | frozen copy of `../radial/background.py` |
| `polar_solver.py` | determinant D(s), exterior Zerilli/KG, κ_eff helper |
| `scan_growing.py`, `scan_*.csv/log` | real-γ scans |
| `convergence.py` (+`convergence_l2_output.txt`) | robustness of D |
| `scan_eps-1_w0.7603_l{4,5,6}.csv` | l = 4, 5, 6 real-γ scans of #4 (§3a); raw logs in negative-mass-raw `I4-win` |
| `convergence_l{4,5,6}_output.txt` | l = 4, 5, 6 convergence variants (§3a), capped rows marked |
| `kappa_eff_l{4,5,6}_output.txt`, `_w0.7603.csv` | l = 4, 5, 6 κ_eff against κ_c (§3a) |
| `winding2.py`, `winding2_output.txt`, `control_winding*.py` | complex-s argument principle + controls (`winding.py` = first version). `winding2.py` takes `--y0`, `--rsmult`, `--zfac`, `--nseg`, `--cachedir`, `--plant`, `--guard` (§3a, #115) and `--mutate` (§5.1) |
| `refine_gamma.py` | real-γ refinement by halving the scan step, with a planted pair and mutations (§5.1) |
| `polish_zero.py` | secant polish of a zero of D near the axis, with the Mini BS I control |
| `axis_scan.py` | \|D\| on a vertical line; its Mini BS control fails, not used (§5.1) |
| `collect_i5.py`, `i5_boxes_output.txt` | the table of §5.1 from the raw logs (negative-mass-raw, I5-rpi5) |
| `collect_i117.py`, `i117_boxes_output.txt` | the table of §5.1's #117 subsection from the raw logs (negative-mass-raw, I117-rpi5): four boxes beyond Re s = 2 and abs(Im s) = 10, each counted with `winding2.py` and `zerilli_zeros.py` at `--guard=all` |
| `collect_i54.py`, `i54_recount_output.txt` | the recount table of §5.2 from the raw logs (I54-rpi5, I5-rpi5, I7-win2) |
| `collect_i115.py`, `i115_boxes_output.txt` | the table of §3a's complex-s subsection from the raw logs (I115-rpi5); exits 1 because s1, l = 5 failed at guard 1; its `--guard=all` rerun is a row of its own |
| `replay_old_rule_i54.py`, `i54_replay_output.txt` | the first step rule replayed on #7's stored boxes, and the lost step (§5.2) |
| `check_cache_i54.py`, `i54_cache_check_output.txt` | stored D values recomputed on a second machine (§5.2) |
| `zerilli_zeros.py` | zeros of the decaying Zerilli solution Z(r_m; s) in a box, the poles of D's Zerilli column (§5.3) |
| `collect_i121.py`, `i121_zeros_output.txt` | the job list and table of §5.3 from the raw logs (negative-mass-raw, I121-win) |
| `i124_zeros_output.txt` | the table of §5.3's #124 paragraph, the nine boxes of §3a, from the raw logs (I124-rpi5), built by `collect_i121.py --only=3a` |
| `collect_i126.py`, `i126_guard_output.txt` | the recount of §5.4 at `--guard=all`, 61 boxes and 5 controls, from the raw logs (negative-mass-raw, I126-rpi5) |
| `diagnose_plant.py`, `i115_replay_output.txt` | the replay of a lost step on a cached box, and its two runs on #115's s1 at l = 5 before and after the guard-all rerun (§5.4) |
| `replay_i126.py`, `i126_replay_output.txt` | the dense-chain replay of §5.4's 61 boxes at `--guard=all` on the caches of negative-mass-raw `I126-rpi5` (`diagnose_plant.py --y0` for the upper boxes); exits 1 when a box fails (§5.4, #128) |
| `d_factors.py`, `d_factors_output.txt` | the checks of D's factors and its pole sources of §5.5 (`--part=factors,gap,norm,pole,axis,mean`, each with a `--mutate=`), 61 boxes, its output from the caches of negative-mass-raw `I122-rpi5`; exits 1 because `harm` and `mean` fail (§5.5, #122) |
| `harm_circles.py`, `harm_circles_output.txt` | #129's mean value test of arg D at the 12 `harm` points of the four boxes with Re s ≥ 0.1, at S8, S10 and S12 (§5.5), with `--mutate=nonholo` and `--mutate=half`; `--collect` reads the caches of negative-mass-raw `I129-win` and exits 1 when a control fails, a cache is missing or D fails at S12 |
| `l456_graph.py` → `l456_mode_condition.svg` | the figure of §3a |
| `boxes_graph.py` → `complex_s_boxes.svg` | the figure of §5.1 |
| `data_graphs.py` | SVG drawing primitives of both figures, vendored unchanged (its docstring names the canonical copy) |
| `kappa_eff.py`, `kappa_eff_l2_*`, `love.py`, `love_output.txt` | κ_eff vs κ_c; static Love-type coefficient |
| `mini_bs.py`, `mini_bs_output.txt` | validation against Macedo et al. 2013 Mini BS I polar QNM |
| `check_bg.py` | background / denominator sanity checks |
| `matter_sv.py`, `matter_sv_output.txt` | scalar-sector SVD diagnostic at r_m over γ ∈ [10⁻³, 0.6] (25 values), run in #8 for (ε, ω, l) = (−1, 0.7603, 2), (−1, 0.7603, 3), (+1, 0.75, 2), (+1, 0.75, 3): sv2/sv1 does not approach 0. Its minimum is 0.5845, 0.5973, 0.6051, 0.6179, each at the smallest γ = 10⁻³, and it rises to 0.959–0.960 at γ = 0.6; sv3/sv1 ≤ 8×10⁻¹⁷ (the generic null vector). The γ-range and the four stars are all it covers |
