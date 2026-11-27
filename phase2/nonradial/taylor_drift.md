# #196: where scipy's polar D(γ) drifts from the crate as γ grows

Status: **numerical evidence**. At (l = 5, γ = 1000, 2000, 3000) the gap between scipy's D and the crate's is the rounding floor of scipy's interior plane
times D's conditioning, not a drift in an input. The declared rules **do not name a cause** ("not established": statuses below); the reading that the data support
is given after them and marked post hoc. Every number is in `taylor_drift_output.txt` (the runs of record; raw folder `negative-mass-raw/I196-mac`, README first with the
criteria, the checks and their controls, and amendments A1 to A7, each declared before the data it concerns). Two things are not met: 7 of the 24 P0 cells were not run (time box, A7), and COMPLETE fails on them.
No sentence here claims novelty.

## Answer

1. **At these γ the gap is realization noise.** It changes sign and size when only the chunk grid changes (one extra chunk of `_integrate_plane` moves D by 4e-10 to 2e-8, as much as the gap, 1.5e-11 to 1.7e-8),
   and the win laptop's stored D differs from this Mac's by 5.6e-10 at (5, 1000). It falls with rtol in the base cells from 1e-11 to 1e-13 (1.2e-9 to 4.9e-11 at 1000, 5.7e-9 to 2.8e-10 at 2000) and not in the φ_cut cells (7.8e-10 to 6.9e-10 at 1000),
   and each point is one realization. On the exact-J system (below) the scatter falls by 0.55, 0.90, 0.16 from rtol 1e-10 to 1e-13 to a floor of 1.6e-10 … 6e-10 that neither rtol nor atol lowers.
2. **The error is in scipy's interior plane at r_m.** The largest principal angle between scipy's interior plane and the crate's double-double plane is 4e-16 … 4e-14 (median 7e-15 at (1000, base), 2.2e-14 at (3000, base));
   the exterior plane agrees to 3.3e-16 and its metric direction to the last bit. D is the product of the sines of the three principal angles between the two planes; the smallest is 4.8e-5 at γ = 1000 and 1.6e-5 at 3000
   (base; half that for φ_cut), so 2e-14 in the plane is 1e-9 in D.
3. **Not the star or M_ext** (P1: the crate's star in `PolarProblem` moves D by less than three times the chunk-grid scatter in 5 of 6 cells), **not the exterior** (the scalar plane at r_m is the local attractor: star J and vacuum J give the same
   D to the last bit; D's sensitivity to q's relative error is 0.2 … 1.5, q's own error 6e-14), **not atol** (1e-14 to 1e-17: no change), **not the interpolation kinks** of the `solve_bvp` star (the scatter does not fall with the crate's smooth star),
   and **not the generated J's rounding alone** (the crate's tables, J rounded once from double-double, rms 1.1 ulp per entry against 3e7 ulp, give scatter 3e-10 … 3e-9 against 7e-10 … 2e-8).
4. **#187's step "the gap grows with γ, so it is not a static input" does not follow, and its conclusion stands by direct test.** The plane error is about the same size at every γ and grows only as γ; the smallest sine falls as 1/γ.
   Post hoc: gap ≈ ε_Y/sin θ_min with ε_Y ∝ γ and sin θ_min ∝ 1/γ is ∝ γ². Fixed at γ = 1000 (A = 1.5e-10 base, 4.4e-10 φ_cut) it gives 5.3e-9 and 1.6e-8 at (5, 6000) against the stored 5.4e-9 and **4.7e-8 (CAP5c's failing cell)**,
   and 2.1e-8 and 6.3e-8 at (6, 12000) against 1.9e-8 and 5.1e-8. CAP5c's bar 1e-8 was below scipy's own floor at that γ.

## What was run

| name | what | cases |
|---|---|---|
| reference | the crate's D (double-double, order 30, tol 1.2e-30), against order 24 at 1e-22 (≤ 4.5e-24) and, at 1000, quad-double (3.9e-27, 1.1e-26) | 6 cells |
| P0 | `PolarProblem` as `convergence.py` builds it (scipy's `solve_bvp` star), base and φ_cut, γ = 1000, 2000, 3000, rtol 1e-10 … 1e-13 | 17 of 24 |
| P1 | the crate's star and its M_ext in `PolarProblem`, rtol 1e-10; r_m is the same grid point (check RMX, last bit) | 6 |
| N0, N1 | P0 and P1 with one extra chunk: the realization noise η | 6 + 6 |
| CJE | the exact-J system: the crate's star, M_ext and J tables; 3 cells × (rtol, atol) ∈ {(1e-10, 1e-14), (1e-13, 1e-14), (1e-13, 1e-17)} × 8 chunk grids | 72 |
| CJV, CJQ, CJN | the same with the exterior scalar plane in vacuum; with q planted 1e-13, 1e-12 too large; with planted per-entry noise of 1, 10, 100 ulp | 24, 4, 12 |
| planes | scipy's two orthonormal planes at r_m against the crate's (`rows planes=1`); the number of renormalizations 1/4 … 8 times | 24 + 4 + 4 |
| PLN | scipy's integrator on #186's planted system (exact J, D ≈ γ⁻³): D error 4e-14 … 2.3e-13, plane error 4e-17 … 3e-16 | 6 |

Not run: P0 at (2000, φ, 1e-11), (3000, base and φ, 1e-11), (3000, base and φ, 1e-12), (3000, base, 1e-13: stopped by its alarm after 7200 s wall, 55 min of it paused by me) and (3000, φ, 1e-13).
The load on the Mac was 4 to 50 from other users, which doubled the times. One (5, 1000) case takes 136 s with the generated J and 11 s with the crate's tables.

## The gap against rtol, the issue's first bullet (|D_scipy / D_crate − 1|, one realization per cell)

| γ | variant | rtol 1e-10 | 1e-11 | 1e-12 | 1e-13 | η: one extra chunk, rtol 1e-10 |
|---|---|---|---|---|---|---|
| 1000 | base | 3.4e-10 | 1.2e-9 | 1.1e-10 | 4.9e-11 | 6.8e-10 |
| 1000 | φ_cut | 7.8e-10 | 5.2e-11 | 6.4e-10 | 6.9e-10 | 1.9e-9 |
| 2000 | base | 2.6e-9 | 5.7e-9 | 9.4e-10 | 2.8e-10 | 1.7e-9 |
| 2000 | φ_cut | 1.5e-11 | not run | 4.1e-9 | 1.6e-9 | 6.8e-9 |
| 3000 | base | 1.7e-8 | not run | not run | not run | 2.1e-8 |
| 3000 | φ_cut | 1.3e-9 | not run | not run | not run | 3.9e-10 |

The crate's D is the reference, so these are the gaps of the issue's text; the signs are mixed (the report prints them). The classifier of the rtol curve (gap(1e-13) / gap(1e-10) ≤ 0.1 truncation, ≥ 0.5 floor) labels the four cells with both ends:
1000 base 0.14 mixed, 1000 φ 0.87 floor, 2000 base 0.11 mixed, 2000 φ 104 floor. A single point cannot tell a trend from its own scatter, so the ensembles below carry the conclusion.

## The planted input, the issue's second bullet (the crate's star and M_ext in `PolarProblem`, rtol 1e-10)

| γ | variant | D moves | 3 η | moves D beyond 3 η | gap without → with |
|---|---|---|---|---|---|
| 1000 | base | 4.1e-10 | 2.0e-9 | no | 3.4e-10 → 6.4e-11 |
| 1000 | φ_cut | 1.2e-9 | 5.6e-9 | no | 7.8e-10 → 1.9e-9 |
| 2000 | base | 3.2e-9 | 5.2e-9 | no | 2.6e-9 → 5.8e-9 |
| 2000 | φ_cut | 4.3e-9 | 2.0e-8 | no | 1.5e-11 → 4.3e-9 |
| 3000 | base | 1.4e-8 | 6.2e-8 | no | 1.7e-8 → 3.4e-9 |
| 3000 | φ_cut | 1.1e-8 | 1.2e-9 | yes (η = 3.9e-10 is one small draw) | 1.3e-9 → 1.2e-8 |

**The planted input does not move the gap beyond the realization noise**: D moves by 4e-10 to 1.4e-8, always at the scale of the chunk-grid scatter, and the gap goes up as often as down. With the crate's star the scatter η is as large as with scipy's (η1/η0 = 1.5, 0.6, 4.9, 0.6, 0.3, 22).
"The crate's r_m" is scipy's r_m to the last bit on both stars (RMX), so substituting it changes nothing.

## The exact-J ensembles (CJE: 8 chunk grids per cell, signed D/D_crate − 1)

| cell | rtol | atol | mean | std | t = mean / (std/√8) |
|---|---|---|---|---|---|
| 1000 base | 1e-10 | 1e-14 | −1.2e-10 | 3.0e-10 | −1.1 |
| 1000 base | 1e-13 | 1e-14 | +8.2e-11 | 1.6e-10 | +1.4 |
| 1000 base | 1e-13 | 1e-17 | −3.4e-11 | 3.5e-10 | −0.3 |
| 1000 φ_cut | 1e-10 | 1e-14 | −4.7e-10 | 6.5e-10 | −2.1 |
| 1000 φ_cut | 1e-13 | 1e-14 | +3.0e-10 | 5.8e-10 | +1.5 |
| 1000 φ_cut | 1e-13 | 1e-17 | +2.3e-10 | 5.1e-10 | +1.3 |
| 3000 base | 1e-10 | 1e-14 | −6.1e-10 | 2.9e-9 | −0.6 |
| 3000 base | 1e-13 | 1e-14 | **−1.3e-9** | 4.6e-10 | **−8.2** |
| 3000 base | 1e-13 | 1e-17 | **−1.6e-9** | 6.1e-10 | **−7.2** |

At (3000, base) the 8 realizations share an offset of −1.3e-9 to −1.6e-9 (all 16 negative) that neither rtol, atol, the vacuum scalar plane (CJV: identical to the last bit) nor q (CJQ) moves.
The planes show where it comes from: ε_Y = 2.4e-14, 7.8e-15, 2.7e-14, 2.0e-14 on four chunk grids, signed gaps −1.5e-9, −4.9e-10, −1.7e-9, −1.3e-9 (sin θ_min = 1.6e-5). That the sign repeats is not explained here.

## Where the error is, and how it scales

| cell (rtol 1e-13, exact J) | ε_Y (interior), 4 chunk grids | median | ε_Z (exterior) | its metric direction | sin θ_min |
|---|---|---|---|---|---|
| 1000 base | 4.2e-15, 4.6e-15, 9.5e-15, 2.0e-14 | 7.0e-15 | ≤ 3.3e-16 | 0 | 4.8e-5 |
| 1000 φ_cut | 5.7e-15, 1.6e-14, 1.4e-14, 7.0e-15 | 1.0e-14 | ≤ 3.3e-16 | 0 | 2.4e-5 |
| 3000 base | 2.4e-14, 7.8e-15, 2.7e-14, 2.0e-14 | 2.2e-14 | 2.2e-16 | 0 | 1.6e-5 |

The median ε_Y grows 3.1 times from γ = 1000 to 3000 (γ ratio 3) and sin θ_min falls 3 times, so the gap ε_Y/sin θ_min grows 9 times. Dose-response to planted per-entry noise on the exact-J system at (1000, base):
std of the signed gap 3.0e-10 at none (the tables' own 1.1 ulp), 1.0e-9 at 1 ulp, 4.8e-9 at 10 ulp, 2.8e-8 at 100 ulp (4 seeds each; no rule was declared for this reading). Scipy's integrator on the planted exact-A system of #186, whose entries are dyadic,
has a plane error of 4e-17 … 3e-16: the same integrator, chunking and tolerances, with no rounding of the data. **Reading, post hoc:** rounding of J's entries to double (about 1 ulp per entry, which no f64 J avoids; the crate's own f64 lane
differs from its double-double lane by 5e-11 … 2.8e-9 at these cells) is amplified by the dynamics to a plane error of ~γ·ε_mach and by D's conditioning to the gap.
The number of renormalizations of the interior plane changes ε_Y little (median 1.6e-14, 9.6e-15, 7.0e-15, 5.5e-15, 6.9e-15, 1.4e-14 at 1/4, 1/2, 1, 2, 4, 8 times; H8's rule gives *undecided*).

## The hypotheses by the declared rules

| | hypothesis | read on | result | status |
|---|---|---|---|---|
| H1 | truncation by rtol | std(rtol 1e-13)/std(1e-10) on the exact-J ensembles: 0.55, 0.90, 0.16 | | undecided |
| H2 | the star and M_ext as inputs | P1 against 3 η: moves D in 1 of 6 cells | | ruled out |
| H3 | rounding in the generated J_and_h | η(crate J tables)/η(J_and_h): 0.33, 0.82, 0.62 | | undecided |
| H4 | kinks of the `solve_bvp` interpolant | η(crate star)/η(scipy star), 6 cells: 1.5, 0.6, 4.9, 0.6, 0.3, 22 | | ruled out |
| H5 | truncation by atol | std(atol 1e-17)/std(1e-14): 2.2, 0.88, 1.3 | | ruled out |
| H6 | the exterior scalar plane through the star's J | the mean gap at (3000, base): −1.33e-9 both ways | ratio 1.00 | ruled out |
| H7 | the Riccati of the exterior q at rtol 1e-11 | D's sensitivity 0.43 × q's error 6e-14 against the offset 1.3e-9: 2e-5 | | ruled out |
| H8 | roundoff accumulation of the renormalizations | median ε_Y at 8× over 1×: 2.0, at 1/4× over 1×: 2.3 | | undecided |

**Cause: not established** (no single hypothesis supported with the rest ruled out). Ruled out: H2, H4, H5, H6, H7. The post hoc reading above is numerical evidence and is the proposal to test next (task 1).

## Checks (`python taylor_drift.py check`; each control run alone and fires)

| check | result | control | effect |
|---|---|---|---|
| IDENT | pass | `ident_atol` | D differs by 1.2e-8, fail |
| CS | pass, 1.9e-16 | `cs_order` | 4.5e-6, fail |
| RMX | pass | `rm_cut` | r_m 30.82 against 34.21, fail |
| CJ | pass, 1.12 ulp rms | `cj_shift` | 1.8e15 ulp, fail |
| REP | pass, 1.5e-9 | `rep_variant` | 1.0, fail |
| REF | pass, 4.5e-24 | `ref_tol` | 1.7e-14, fail |
| CLS | pass | `cls_swap` | labels swapped, fail |
| PLN0 | pass, 1.1e-14 | `pln_sign` | 2.0, fail |
| COMPLETE | **fail**: 147 of 154 | `drop_cell` | 146 of 154 (on a complete synthetic set it passes and the control fails it, `dev_43`, `dev_44`) |

## Proposed section for `taylor_polar.md`

> ## #196: where scipy's D drifts from the crate
>
> Status: **numerical evidence**; cause by the declared rules not established. At (l = 5, γ = 1000 … 3000) the gap between `PolarProblem`'s D and the crate's is not a drift: one extra chunk moves scipy's D by 4e-10 … 2e-8,
> and the gap does not fall with rtol (1e-10 … 1e-13) or atol (1e-14 … 1e-17). The crate's star and M_ext in `PolarProblem` move D by no more than that scatter; the exterior (the scalar plane, q), the star's interpolation and the generated J's rounding
> do not make it. Scipy's interior plane at r_m differs from the crate's double-double plane by 4e-16 … 4e-14, the exterior plane by 3e-16; D is the product of the sines of the three principal angles between the planes, the smallest 4.8e-5 at γ = 1000 and
> 1.6e-5 at 3000, so the plane error is amplified 2e4 … 6e4 times. **#187's "the gap grows with γ, so it is not a static input" does not follow**: a plane error that grows as γ and a sine that falls as 1/γ give a gap ∝ γ².
> Post hoc, A = ε_Y/sin θ_min at γ = 1000 gives A (γ/1000)² = 5.3e-9 and 1.6e-8 at (5, 6000) against the stored 5.4e-9 and 4.7e-8 (CAP5c), 2.1e-8 and 6.3e-8 at (6, 12000) against 1.9e-8 and 5.1e-8: **CAP5c's bar 1e-8 was below scipy's own floor**.
> Scipy with the crate's J tables (J rounded once from double-double) runs 10 times faster (11 s at (5, 1000)). Raw `I196-mac`; script `taylor_drift.py`, note `taylor_drift.md`.

## Defects and discrepancies found

- **#187's text** ("so it is not a static input such as the star or M_ext, which would be the same at every γ", `taylor_polar.md`) rests on a step that does not follow (above). The conclusion it draws about the star and M_ext holds by direct test (P1, planes).
- `PolarProblem.zerilli_q` hard-codes `rtol = 1e-11, atol = 1e-14` and ignores `self.rtol`: the issue's rtol sweep never reaches the metric exterior. It does not matter here (CJQ), and it is worth knowing for a sweep of `PolarProblem`.
- `scalar_plane` integrates through the star's J up to r = 60 where the crate uses vacuum J. At r_m the plane is the local attractor and the two give the same D to the last bit (`dev_48`).
- **Generated J in f64 loses digits**: against the same source at 60 digits, rms 2651 ulp per entry for r ≥ 0.5 and 3e7 ulp from r = 0.002, but 1.05 ulp of its row maximum (`jround`). The integrator does not thrash on it, and D is as noisy with exact J.
- The unit of a "realization": scipy's `D` is deterministic on one machine (the dev and record values at (1000, base) are equal to the last bit) and differs between machines by 3.8e-10 … 1.5e-9; the stored comparisons of #187 are single draws.
- **Not met:** 7 of the 24 P0 cells (time box); COMPLETE fails. The case (3000, base, 1e-13) was stopped by its own alarm, which counted a pause I made. The planted-noise cases at 0.1 ulp of the row maximum were not run (0.01 ulp took 500 s a case, `dev_38`, `CJN`).
- My mistakes (each in the raw README): no pilot of the planted-noise knob before the first run (`dev_38`), `set -- $m` in zsh in the first control loop (`dev_09` … `dev_13`), a python heredoc with backticks that ran them, six amendments after reading partial data (each declared before the data it concerns).

## Recommended tasks

1. **Test the post hoc reading at (5, 6000) and (6, 12000).** What: the CJE and planes protocol with the crate's J tables at those γ (minutes each, not hours), with 8 realizations per cell, and the gap of #114's stored values against the ensemble.
   Why: CAP5c's failing cell and #198's comparison. Accept: for each cell the std of the signed gap and the median ε_Y, a prediction A (γ/1000)² made before the run, and a fail if the predicted rms is more than 3 times off or #114's stored value lies beyond 3 std from the ensemble mean. Depends: none.
2. **Set #198's bar from the floor.** What: before the l = 6 value of #114 is compared, the bar is 3 × the predicted rms (task 1), not 1e-8. Accept: the bar is committed before the value is read. Depends: 1, #114.
3. **Why the sign repeats.** What: at (3000, base) the 16 realizations of CJE share a negative offset of 1.3e-9 … 1.6e-9, with ε_Y 8e-15 … 3e-14. Candidates: a systematic part of the interior plane's rounding. Accept: a planted test that moves the offset, or the offset gone with the interior plane in double-double. Depends: none; low priority.
4. **A double-double scipy path.** What: `PolarProblem` with J and the plane arithmetic in double-double (or the crate's `d_eval` with the chunking of scipy). Why: the only direct test that the floor is f64 rounding of J. Accept: the gap against the crate falls to 1e-14 or better at (5, 3000). Depends: none.

## Coordinator's review

Review of branch `wip/i196-polar-drift` at 64b4e35 and raw `I196-mac` at edd24a2a, recorded in negative-mass-raw `I196-mac-review`.
On the Mac, nice 10, under a load of 20 to 60 from other projects.

- **Scope.** The branch adds the script, its output and this note. The crate gains commands only: `starsteps`, `jtable`, the variant
  `phio24`, and `rows planes=1`. `DOut` gains the two planes. The default rows are unchanged (REP, and #187's rows through `rows`).
- **Rerun.** `taylor_drift.py check` on the stored runs, with the crate built in the review worktree, gives the note's verdicts: 8
  pass, and COMPLETE fails on the seven cases not run. The first two rerun attempts failed on the coordinator's side, from missing
  arguments and an unbuilt crate; their logs are kept.
- **No independent code.** This is a diagnostic study. Its claims rest on controlled comparisons: the extra-chunk realizations, the
  exact-J ensembles, and the plane angles against the crate's double-double planes. REP compares this Mac with the win laptop's stored D,
  so two machines enter.
- **Status.** The cause is "not established" by the declared rules, as the note says. The γ² reading is post hoc, and it predicts
  #187's top-γ gaps to within a factor of 3.
- **Applied.** `taylor_polar.md` gains the proposed #196 section and a correction of #187's step that does not follow.
- **Follow-ups.**
  - Task 1 (test the reading at (5, 6000) and (6, 12000)) is #271.
  - Task 2 (#198's bar from the floor) is added to #198.
  - Tasks 3 and 4 are not filed: the offset's sign and a double-double scipy path matter only if #271 fails.
