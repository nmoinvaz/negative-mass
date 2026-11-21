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
the half-step test of guard 1 does not run. No threshold was changed. The raw folder's `diagnose_plant.py` replays the counter on the 212 cached points of
that box, with no new evaluation of D (`diagnose_s1.log`, `diagnose_all.log` there; the script is not in this repository).
- It reproduces the logged cumulative phases of D on all three edges and the planted count of −0.0000.
- On the left edge, Re s = 0.03, one accepted step of D times the planted factor, from Im s = 0.25 to 0.1875, measured
  −0.3826 rad. The 14 cached points of that stretch sum to +5.9006 rad (largest step 0.73 rad), so the step lost one 2π
  and the count lost 2. The planted zero is 0.01 from that edge, so its factor turns by about π within one step.
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
So the count of D at s1, l = 5 is confirmed by its planted check, with `--guard=all`. The other eight boxes passed their
check with guard 1 and are not recounted here.

**Controls of the collector.** On a synthetic raw folder built from #5's logs (tags renamed to l = 4, 5, 6) the collector
exits 0, and each of `--mutate=noplant`, `wrongexpect`, `wrongbox`, `shift` and `nocontrol` exits 1. On the real logs
each of the first four lists all ten rows as FAILED, `nocontrol` lists the control and the one real failure, and
`--mutate=wrongguard` (expects guard all everywhere) lists the nine default rows.

**Pole caveat.** Each winding is zeros minus poles of D. D has one pole source, the Zerilli column (1, q/f), whose q has a
pole at each zero of the Zerilli solution at r_m (#54, #121). #121 counts those poles in every counted box, this task's
included. **Until then no box of this task is covered by a pole exclusion**, and a count of 0 is not a proof that there is
no zero.

**Not covered.** Only the stored box, s1 and n1 were counted. The other boxes of §5.1 (n3, h1, h2, s2, s3 and the
threshold boxes t1 and t2 at two outer radii) are not run for l = 4, 5, 6; Re s < 0.03 for 0.23 < |Im s| ≤ 2.5 and
|Im s| > 2.5 are open, as in §5.1. The real-γ axis at l = 4, 5, 6 is the 100-point scan above, not refined.

Status: *numerical evidence*, conditional on D having no pole in each box. Eight boxes passed their planted check at
guard 1; s1 at l = 5 failed it at guard 1 and passed it with `--guard=all`.

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
0, and the −2 was the first step rule losing one 2π on the left edge. The hypothesis is untested. Each winding counts
zeros minus poles, so the 0 counts remain conditional on D having no pole in the box.

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

#### Not covered

- Re s below the lower edges for 0.23 < |Im s| ≤ 2.5. For l = 2 that is 0 < Re s < 0.01 (the near-axis box n1 reaches
  10⁻⁴ only below Im s = 0.23), for l = 3 it is 0 < Re s < 0.03. R_s grows as 1/Re s (3831 at 0.01), and so does the
  cost of a box, which took 3.7 h at 0.01. A mode that close to marginal is neither shown nor excluded by these counts.
- The strip Re s < 0.03 for 2.5 < |Im s|, and the strip Re s < 0.1 for |Im s| > 6.
- Re s ∈ [0.1, 2] above |Im s| = 10 (l = 2) and 6 (l = 3), and Re s > 2.
- The sliver |Im s| < 0.01 left of Re s = 10⁻⁴.
- Poles of D. Every count is zeros minus poles, and a count of 0 is not a proof that there is no zero. #54 found no net
  pole (§5.2) but did not search for poles.
- l ≥ 4 (#4).

Status: numerical evidence, conditional on D having no pole in each box (#54).

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
  of §5 and of the "Status" line in §5.1 still holds, since poles were not searched for.

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
   (§5.2). This means 0 zeros only if D has no pole in the box, which was not searched for.
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
- Complex s, campaign star, l = 2 (l = 3 where it differs): counted are Re s ∈ [0.1, 2] for |Im s| ≤ 10 (6 for l = 3),
  Re s ∈ [0.03, 0.1] for |Im s| ≤ 6, Re s ∈ [0.01, 0.03] for |Im s| ≤ 2.5 (l = 2 only), and Re s ∈ [−10⁻³, 0.1] for
  |Im s| ≤ 0.23 (§5.1). The rest is listed under "Not covered" there. The other stars have Re s ≥ 0.1, |Im s| ≤ 2.5 only.
- The windings of §5 were computed with a step rule that can lose 2π in one step (§5.1). One loss shifts a symmetric-box
  count by ±2. They were recounted with the current rule in §5.2. Only ω = 0.85, l = 3 changed, from −2 to 0, because the
  old rule lost one 2π there.
- l = 2, 3 for every star, and l = 4, 5, 6 only for the ω = 0.7603 phantom star (§3a, #4), with the real-γ scans,
  convergence variants and κ_eff, windings in three boxes only (the stored box, s1 and n1; eight of the nine returned
  their planted pair at guard 1, s1 at l = 5 did only with `--guard=all`; §3a, #115) and no c_l. The ω = 0.7603, 0.70 and 0.85 stars have the l = 2, 3
  real-γ scans, convergence variants and windings (0.70 and 0.85 in #7). c_l and axial positivity are also checked at
  ω = 0.85, and axial positivity at 0.70 and 0.90. At l = 5 and 6 the rtol 1e-12 variant is missing at the top γ
  (§3a).
- The windings count zeros minus poles of D. The −2 that suggested a net pole at ω = 0.85, l = 3 was the old step rule
  (§5.2), so no box shows a net pole. Poles were not searched for, and the hypothesis of #54 that the Zerilli column has
  one where ψ(r_m; s) vanishes is untested.
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
| `scan_eps-1_w0.7603_l{4,5,6}.csv` | l = 4, 5, 6 real-γ scans of #4 (§3a); raw logs in negative-mass-raw `I4-win` |
| `convergence_l{4,5,6}_output.txt` | l = 4, 5, 6 convergence variants (§3a), capped rows marked |
| `kappa_eff_l{4,5,6}_output.txt`, `_w0.7603.csv` | l = 4, 5, 6 κ_eff against κ_c (§3a) |
| `winding2.py`, `winding2_output.txt`, `control_winding*.py` | complex-s argument principle + controls (`winding.py` = first version). `winding2.py` takes `--y0`, `--rsmult`, `--zfac`, `--nseg`, `--cachedir`, `--plant`, `--guard` (§3a, #115) and `--mutate` (§5.1) |
| `refine_gamma.py` | real-γ refinement by halving the scan step, with a planted pair and mutations (§5.1) |
| `polish_zero.py` | secant polish of a zero of D near the axis, with the Mini BS I control |
| `axis_scan.py` | \|D\| on a vertical line; its Mini BS control fails, not used (§5.1) |
| `collect_i5.py`, `i5_boxes_output.txt` | the table of §5.1 from the raw logs (negative-mass-raw, I5-rpi5) |
| `collect_i54.py`, `i54_recount_output.txt` | the recount table of §5.2 from the raw logs (I54-rpi5, I5-rpi5, I7-win2) |
| `collect_i115.py`, `i115_boxes_output.txt` | the table of §3a's complex-s subsection from the raw logs (I115-rpi5); exits 1 because s1, l = 5 failed at guard 1; its `--guard=all` rerun is a row of its own |
| `replay_old_rule_i54.py`, `i54_replay_output.txt` | the first step rule replayed on #7's stored boxes, and the lost step (§5.2) |
| `check_cache_i54.py`, `i54_cache_check_output.txt` | stored D values recomputed on a second machine (§5.2) |
| `l456_graph.py` → `l456_mode_condition.svg` | the figure of §3a |
| `boxes_graph.py` → `complex_s_boxes.svg` | the figure of §5.1 |
| `data_graphs.py` | SVG drawing primitives of both figures, vendored unchanged (its docstring names the canonical copy) |
| `kappa_eff.py`, `kappa_eff_l2_*`, `love.py`, `love_output.txt` | κ_eff vs κ_c; static Love-type coefficient |
| `mini_bs.py`, `mini_bs_output.txt` | validation against Macedo et al. 2013 Mini BS I polar QNM |
| `check_bg.py`, `matter_sv.py` | background / denominator sanity checks; scalar-sector SVD diagnostic (not run) |
