# #186: the polar matching determinant D(γ) by a Taylor integrator in Rust

Status: **numerical evidence**. The crate `taylor_polar/` reproduces the stored base-variant rows of
`convergence_l{2,4,5,6}_output.txt` (the digits they print) at every γ, including l = 6 at γ = 1000 and 12000, and a fresh scipy
evaluation here to 4.5e-9 at worst. Its time per evaluation is 9 to 21 times lower than the scipy path in double-double and
114 to 750 times lower in f64, on the same Mac. Every number below is in `taylor_polar_output.txt` (the runs of record; the
raw folder is `negative-mass-raw/I186-mac`, README first, last commit when the output was assembled in its first line).
No sentence here claims novelty.

## What it is

- **Crate** `taylor_polar/` (Cargo.toml, Cargo.lock, src/). `taylor_polar.py` builds it, regenerates `src/jgen.rs`
  (`--gen`), runs the checks, and runs a control with `--mutate=NAME`. `python taylor_polar.py` takes about 13 min and exits 1 on
  any failed check.
- **Construction** as in `polar_solver.py`: interior 3-plane from r0 with flat data r^l, Riccati q = f Z′/Z from
  r_m + 40/γ + 20 l/γ + 50 and the Moncrief map, the two decaying scalar solutions from the contracting eigenspace at R_s,
  D = det of the orthonormalized planes at r_m (Gram-Schmidt with positive diagonal, one reorthogonalization after every step,
  which gives the same signed D as the LAPACK QR with its phase product).
- **Integrator.** Taylor order 30 for y′ = J y from the jets of J's entries; the step is the largest whose last two terms are
  below `tol` (multishoot's `plan_path` rule; the crate does not call multishoot's `Equations` or `Shooting`, whose 8-lane
  double-double boundary-value layout has no place for the orthonormalization between steps). Lanes f64, double-double,
  quad-double, octo-double are multidouble's scalar types. Reciprocals 1/k are formed in the lane (R1 holds, and its control
  `recip_dd` fires). In #186 multishoot 4.2.0 and multidouble main ba0ad101ce76 (series module) resolved together by `[patch.crates-io]`
  (multishoot's requirement is the caret "5.1.0"); no multishoot change was needed. Since #192 the crate takes multishoot
  5.0.0 and multidouble 6.1.0 from crates.io, with no patch. `mullow` and `Relaxed` are not used: their
  crossovers (512 in double-double, 64 in quad-double, 16 in octo-double) lie above the 30 coefficients of every product here.
- **J.** `jgen.rs` is the syntax-tree translation of `polar_matrix_generated.py`, which `reduce_polar.py` wrote (the issue's
  `gen_matrix2.py` is not in the repository). **The symbolic source is shared with `polar_solver.py`; it is not independent.**
  J = J0 + s² J1 is stored as Taylor tables along the star, one per l.
- **The star** by Newton shooting in the crate (6 secant steps, P_c by bisection on the shot's class), its own Taylor series by
  step; nothing from scipy's mesh. Against `solve_star` at 10 radii (α, a, φ0, u): 4.1e-8 for the exact regular start, and
  1.05e-10 when the shooting starts at r0 = 1e-3 with scipy's leading-order boundary condition. So the 4e-8 is scipy's
  boundary condition at r0, not shooting error; φ_c = 0.0196642435159 (exact) against 0.0196642426625 (`solve_star`). The
  central data move with the start radius as rs², ratio 1.00e-4 between rs = 1e-8 vs 1e-9 and 1e-6 vs 1e-7 (S3).

## Results of the runs of record

| Check | Result |
|---|---|
| A2, l = 2, 8 stored rows (1e-3 … 137) | worst \|D/D_stored − 1\| 4.2e-7 (the rows print 7 digits), signs all negative |
| A6, l = 6, 11 rows (… 1000, 12000) | 4.6e-7; γ = 1000: −2.7575156e-10, γ = 12000: −1.5957913e-13, stored −2.757516e-10, −1.595791e-13 |
| A4 (l = 4, …, 3000), A5 (l = 5, …, 6000) | 3.9e-7, 4.0e-7 |
| TS, a fresh scipy D here, 11 points | crate against scipy 1.0e-13 … 4.5e-9; l = 2, γ = 300 (the row the file lacks): scipy −1.0204252e-9, crate −1.0204252e-9 |
| V4, convergence.py's variants through both, 12 cases | crate against scipy 8.5e-10 worst |
| P1 planted (growth e^{γ r}, L = 24, γ ≤ 12000), planes and D | f64 8e-16, dd 5e-32, qd 3e-64; ln vol relative 2e-12, 1e-28, 2e-62 |
| P0 transient, P2 variable coefficients | dd 1.2e-32; 2.5e-30 (mpmath's own limit) |
| T1, T2, T3 (J and its table against mpmath and the direct evaluation) | 5.5e-28, 3.1e-29, 4.6e-23 |
| L1, dd against qd at (2, 300), (6, 137) | 5.7e-26, 5.4e-29 |
| V3, order 24 at 1e-22 against order 30 | 4.7e-25 |
| the 14 controls | each fails its targeted checks (output, run_03 and run_06) |

**Which lane.** f64 differs from double-double by 6.1e-10 at (2, 300), 4.6e-11 at (6, 1000), 6.9e-9 at (6, 12000): inside the
2.1e-4 bar by four orders, so f64 is enough to reproduce these rows. The error grows as 1/|D| (D falls like γ⁻³). Double-double
agrees with quad-double to 5.7e-26 and 5.4e-29 where quad-double ran, so it is converged; quad-double at γ = 12000 was not run (about
10 times the steps at 10 times the cost per step). The planted system with D ∼ γ⁻³ down to 2.5e-12 is reproduced in f64 to 8e-16.
Choptuik-large-d #1082's noise did not appear: that determinant was formed from unnormalized solutions, and here the planes are
orthonormalized after every step.

**The stored spread of 2.1e-4 is not integration error.** V4 shows each variant of `convergence.py` in scipy and in the crate at
the same (r0, r_m): rtol 1e-12 moves D by at most 2.3e-9, R = 90 by +1.13e-4 and r0 = 3e-3 by −9.9e-5, and the crate reproduces
each shift to 2e-14 … 8.5e-10. R = 90 and r0 = 3e-3 change the grid choice of r_m (24.1084, 24.1110 against 24.1098), and D, a
normalized determinant, depends on r_m at that size.

## Timings (one process, one thread, nice 10, load 1.8 to 3.0, Apple M5)

Per evaluation, base variant, excluding the one-time star shooting (2.7 s) and J table (3.2 s per l, double-double whatever the
lane) of the crate, and scipy's `solve_star`; the scipy values are one evaluation each (run_01), the crate's from run_02.

| l | γ | scipy (s) | crate dd (s) | scipy / dd | crate f64 (s) | scipy / f64 |
|---|---|---|---|---|---|---|
| 2 | 1 | 0.4 | 0.046 | 9 | 0.0031 | 138 |
| 2 | 10 | 1.5 | 0.090 | 16 | 0.0032 | 459 |
| 2 | 40 | 5.2 | 0.278 | 19 | 0.0086 | 602 |
| 2 | 137 | 17.7 | 0.888 | 20 | 0.0257 | 687 |
| 2 | 300 | 40.2 | 1.918 | 21 | 0.0544 | 739 |
| 6 | 0.001 | 0.5 | 0.054 | 10 | 0.0047 | 114 |
| 6 | 1 | 0.5 | 0.048 | 10 | 0.0032 | 157 |
| 6 | 40 | 5.2 | 0.280 | 18 | 0.0086 | 602 |
| 6 | 137 | 17.4 | 0.889 | 20 | 0.0256 | 679 |
| 6 | 300 | 38.8 | 1.923 | 20 | 0.0544 | 714 |
| 6 | 1000 | 135.6 | 6.501 | 21 | 0.1807 | 750 |

At l = 6, γ = 12000 the crate took 76.4 s (double-double) and 2.13 s (f64). The scipy
evaluation there was not run (the brief's instruction); #114's 89 h is the rtol 1e-12 variant on surface-laptop-4 and is not
comparable to these. Both costs grow linearly with γ (scipy 38.8 s at 300, 135.6 s at 1000; the crate 6.4 s at 1000, 76.4 s at 12000).
Peak memory: scipy run 96 MiB, checks 107 MiB (reported only).

## What it takes to run #114's rows on it

The rows print five D variants and the interior and exterior Zerilli log-derivatives. The crate has the base variant, r0 = 3e-3,
r_m = 34.2 (V1, V2), the order/tolerance analog of rtol (V3) and the r_m shifts of the R = 90 variant (V4), and prints the exterior
q. Missing: the interior q (the 4 × 3 null-vector step of `interior_zerilli`, about 100 lines in the lane) and a driver that prints
rows in `convergence.py`'s format. A row at l = 6, γ = 12000 is then five evaluations: about 11 s in f64, 6.5 min in double-double
on this Mac, against hours in scipy.

## Defects and discrepancies found

- **J's table steps.** The first table cut only where the star's series needed: near r = 0 the star's steps are long (its
  functions are regular), but J has 1/r and 1/r² terms, so a shift of J's series by a whole star step went past its radius
  (r ≈ r0 to 2e-3: errors up to 560 %), and the sign of D depended on the star's step structure. Table steps now stop at 0.3 of the
  width where J's own series has converged to 1e-25. T3 and the control `wide_table` catch it; the earlier agreement with stored
  rows at one start radius was partly luck. See the raw README for the logged reproduction (`dev_40`–`dev_44`).
- **`convergence_l2_output.txt` has eight rows, not nine** (the script's list has γ = 300; both entered in commit 81cc13d and the
  file has no note). The coordinator confirmed the issue's "nine" was an error. A fresh scipy D at (2, 300) is in TS.
- **`gen_matrix2.py` is not in the repository**; the source is `reduce_polar.py`'s generated file.
- Harness defects found by the checks and fixed (the raw README lists them): decimal planted inputs not exact in both languages,
  mpmath's `expm` and `eig` limits, a 40-digit print hiding quad-double, and an f64 square of s in the table check.

## Recommended tasks

1. **Rows of `convergence.py` from the crate.** What: the interior q, a `convergence_rust.py` printing rows in the stored format, and
   l = 2 … 6 rows at the stored γ, plus #114's two capped rows. Why: they take minutes. Accept: every stored D agrees to its printed
   digits and q_int, q_ext to theirs, a check fails if one does not, with a control. Depends: this issue.
2. **The star by multishoot's `Shooting`.** What: the star's boundary-value problem as an `Equations` with exact Jacobians in
   double-double, compared with the crate's shooting. Why: a second solver for the star (S1 compares with scipy's mesh and S2 with
   scipy's right-hand side; no check solves the boundary-value problem a second way). Accept: central data agree to 1e-16, a
   mutated equation fails. Depends: none.
3. **Complex s.** What: the same engine over a complex scalar, for the argument-principle counts of #115, #117, #121 to #126 that now
   take hours per box. Why: those counts call D thousands of times. Accept: D(s) on a stored box agrees with the scipy path's
   phase and modulus to 1e-9. Depends: 1.
4. **Quad-double at large γ.** What: one evaluation at (6, 12000) in quad-double (about 3 h estimated). Why: closes the lane
   question where D is smallest. Accept: agreement with double-double to 1e-20. Depends: none; low priority, since double-double
   agrees with quad-double to 1e-26 where both ran.

## Coordinator's review

Reviewed 2026-11-26 at f4edcb1 (raw `I186-mac` 2a4fbb7). The review's files are in `negative-mass-raw/I186-mac-review`
(raw 77ae052), with its criteria committed before the first run of `ind186.py` (raw f143b75).

| item | result |
|---|---|
| scope | the crate, the script, its output and this note, and the raw folder |
| rerun on the Mac | the default mode: 22 checks reported passed, 873 s. Every value equals the committed output to its printed digits; only per-γ timings differ |
| independent check `ind186.py` | out of sample, at (l, γ) = (2, 50), (3, 10) and (4, 300), which the runner did not use. The crate's dd agrees with a fresh scipy D to 1.3e-10, 1.2e-12 and 2.1e-10, and f64 agrees with dd to 3.9e-11 at worst. 2 of 2 criteria pass; `f64-claim` fires on target, and `wrong-l` fires R1 only, off target against its declared R1, R2, as its README anticipated |

Findings:
- **The result stands** as numerical evidence. The crate reproduces the scipy path's D at the stored rows and at points
  it was not tuned on. The 2.1e-4 spread of the stored rows is the shift of r_m between variants (V4), not integration
  error.
- **The speedup grows with γ.** Each call has a fixed setup of about 3 s (the star's shooting and J's tables). At small
  γ the scipy path is as fast or faster: (3, 10) takes 1.4 s in scipy against 3.1 s in the crate. At (4, 300) it is
  37.9 s against 5.1 s, and at l = 6, γ = 12000 about 6.5 min against 89 h. The note's 9–21× is per evaluation, after
  setup.
- **A reporting defect.** The default mode prints check TS as SKIPPED when no `--scipy-file` is given, and still counts
  it among the passed checks. A skipped check should be counted apart.
- **Found and fixed by the runner:** the J table's steps ran past J's 1/r² radius of convergence. T3 and `wide_table`
  now catch it.

Follow-ups: #187 (`convergence.py`'s rows and #114's capped rows from the crate), #188 (complex s for the box counts).
The note's recommended tasks 2 (the star by multishoot's `Shooting`) and 4 (quad-double at (6, 12000)) are not filed:
they are low priority beside #187.

## #187: convergence.py's rows, and #114's two capped rows, from the crate

Status: **numerical evidence**, with one criterion **not met** (CAP5c below). The crate prints `convergence.py`'s rows for l = 2 … 6 in the
stored format. Every token of the four stored files (41 rows of 7 tokens and 40 header tokens, 327 tokens) equals the crate's token as a string. The two
capped rows of #114 come out with D < 0 in all five variants. Against #114's finished l = 5 evaluation the φ_cut variant differs by 4.7e-8
where the bar, set before any run, was 1e-8. Every number below is in `taylor_rows_output.txt` (the runs of record; raw folder
`negative-mass-raw/I187-mac`, README first with the criteria, amendments A1 to A3, last commit when assembled in the file's second line).
No sentence here claims novelty.

### What was added

- **q_int** (`polar.rs: interior_q`): the null vector of [Y[2:] | −S] (4 × 5) by signed 4 × 4 cofactors in the lane, (H1/s, K) = Y[:2] c, the
  Moncrief map, q_int = f Z′/Z. The scipy path takes an SVD null vector; the symbolic source, the planes and the Moncrief map are shared.
  q_ext is #186's `q`.
- **The r_m rule of `PolarProblem`** (`rows.rs: grid_rm`): the first point with φ0 < φ_cut of numpy's `linspace(r0, R, 20000)`, φ0 from the crate's
  own star. The crate's r_m equals scipy's to the last bit for all five variants (RM); the margin of the rule (ln φ/φ_cut at r_m) is at least 7.8e-5
  (the R = 90 variant), 1.0e-3 for the base.
- **The five variants** (`rows` command; one star and one J table per l, then one `d_eval` per (γ, variant), each a flushed JSON line):

  | convergence.py | the crate's analogue |
  |---|---|
  | base | order 30, tol 1.2e-30 (dd), φ_cut 1e-9, grid [1e-3, 60], r0 = 1e-3 |
  | `rtol 1e-12` | **the crate has no scipy rtol**: a looser step control, order 30, tol 1e-12. Printed `rtol 1e-12 [crate: Taylor tol 1e-12]` |
  | `phi_cut 1e-12` | as base with φ_cut 1e-12 (r_m = 34.2051) |
  | `background R=90` | the same star (P_c, Z_c are the same 33 digits for classification radius 60 and 90) and the grid over [r0, 90] (r_m = 24.1084) |
  | `r0 = 3e-3` | interior start at 3e-3 and the grid from 3e-3 (r_m = 24.1110) |

  The spread keeps convergence.py's quirk: it is taken over the names without `phi_cut`, which leaves base out (#114 found this).
- **`python taylor_polar.py --rows`** prints the rows and runs the checks below; `--mutate=NAME` runs a control; `--jsonl-dir=DIR` keeps and resumes the
  crate's lines. #186's 22 checks still pass by default (run_18: 22 passed, 0 failed, 1 skipped of 23).
- **The reporting fix.** A skipped check reports `SKIPPED` and is counted apart: before, `--only=TS` printed `0 of 0 checks passed` and exit 0 (`dev_37`);
  now `0 passed, 0 failed, 1 skipped (TS: not run) of 1` (`dev_38`). My first version missed a numpy bool (S2's), so `run_07` summed to 22 of 23; the verdict is now
  coerced to a bool and the default mode was run again as `run_18`.

### Results of the runs of record

| Check | Result |
|---|---|
| QP, planted q_int (dyadic planes, exact sympy null space, l = 3) | well-conditioned f64 6.6e-15, dd 3.6e-31, qd 2.1e-63; near-degenerate planes (2⁻²⁷) f64 3.2e-10, dd 1.8e-26, qd 1.8e-58 |
| RP, planted r_m (φ = 0.02 e^{−0.65 r}) against numpy | 4 of 4 cases equal (index and float) |
| RM, r_m of the five variants against `PolarProblem.rm` on scipy's stars | 5 of 5 equal to the last bit |
| RA2, RA4, RA5, RA6: every token of the stored file as a string | 66 of 66, 87 of 87, 87 of 87, 87 of 87 |
| RX, l = 2, γ = 300 (no stored row) against #186's fresh scipy D | −1.0204252e-9 both |
| CAP5a, CAP6a: D < 0 in the five variants | met, both |
| CAP5b, CAP6b: the rtol analogue within 2.1e-4 of base | E1 = 1.4e-14, 1.8e-15 |
| **CAP5c**: the five crate D against #114's finished scipy D, 1e-8 each | **not met**: base 5.4e-9, rtol 3.1e-9, R = 90 8.0e-10, r0 4.2e-9, **φ_cut 4.7e-8** |
| **CAP5s**: the same with the crate's star started as scipy's (rs = 1e-3) | **not met**, the same five numbers |
| #186's default mode | 22 passed, 0 failed, 1 skipped (TS), 847 s |
| the controls | each fires on its declared targets (below) |

The closest any D or q token comes to a rounding boundary of its print is 8.3e-10 relative (l = 2, q_ext at γ = 0.001), 6.2e-9 (l = 4), 1.5e-9 (l = 5, q_int
at γ = 0.001) and 5.2e-9 (l = 6, D at γ = 137): the equality of the strings is not by luck of rounding at the 1e-8 scale, but q_ext at γ = 0.001 is within 1e-9 of flipping.

**The capped rows** (dd, one thread; D, crate):

| | (l = 5, γ = 6000) | (l = 6, γ = 12000) |
|---|---|---|
| base | −8.93608292806e-13 | −1.59579128004e-13 |
| rtol analogue (Taylor tol 1e-12) | −8.93608292806e-13 | −1.59579128004e-13 |
| φ_cut 1e-12 | −4.43895159701e-13 | −7.92691720653e-14 |
| R = 90 | −8.93709571709e-13 | −1.59597214579e-13 |
| r0 = 3e-3 | −8.93519613293e-13 | −1.59563291452e-13 |
| stored row (4 variants, scipy) | D_base −8.936083e-13, spread 2.1e-4 | D_base −1.595791e-13, spread 2.1e-4 |
| #114's scipy rtol 1e-12 (win) | −8.936082956e-13 | still running, not touched |
| row time, five evaluations | 176 s | 352 s |

Error bars, each named and none rigorous:
- **E1**, the spread over two Taylor tolerances (1e-12 and 1.2e-30, order 30, dd): 1.4e-14 and 1.8e-15. An estimator range for the looser setting.
- **E2**, the spread over order and tolerance (order 24 at 1e-22 against order 30 at 1.2e-30, dd): 1.3e-24 and 2.4e-26.
- **E3**, the spread over two lanes (f64 against dd): 1.5e-9 and 6.9e-9. This is f64's roundoff scale, not dd's, which is smaller by about 1e-16. At (5, 6000) the f64 lane differs from dd by
  1.5e-9 (base, rtol, φ_cut), 1.7e-9 (R = 90) and 7.4e-9 (r0); at (6, 12000) by 4.7e-9 to 1.7e-8 (`run_20`).
- **E4**, the spread over two solvers (the crate against #114's scipy evaluation): 5.4e-9 (base) at l = 5; it contains scipy's own error, which #114 measures as 2.33e-9
  between its rtol 1e-10 and 1e-12 evaluations here.
- The 2.1e-4 spread between the variants is not an error bar: it is the shift of r_m (24.1098, 24.1084, 24.1110) to which D, a normalized determinant, is sensitive (#186, V4); R = 90 is +1.13e-4 and
  r0 = 3e-3 is −9.9e-5 of base at both rows.
- No quad-double evaluation at these γ was made (about 1.5 to 3 h each); #186's dd against qd 5.7e-26 and 5.4e-29 at (2, 300) and (6, 137) is cited, not repeated.

**CAP5c is not met, and CAP5s says the star's start is not why.** The φ_cut variant at (5, 6000) differs from scipy's by 4.7e-8 against the bar 1e-8, set a priori as twice
#186's worst 4.5e-9; the coordinator's ruling kept it. CAP5s repeats the comparison with the crate's star shot from rs = 1e-3 with scipy's leading-order boundary condition (the start that
agrees with `solve_star` to 1e-10 in #186's S1): the five numbers are the same, and the crate's D does not change in its f64 digits between the two starts (a development value at l = 2, γ = 1 in dd: 5e-21, `dev_21`), because the star beyond
the centre is the same decaying solution (m(r_m) agrees to 5e-18). What `run_19` shows beyond that: the crate agrees with scipy's full-precision D to 8e-13 up to γ = 10 in every variant and every l (4e-12 at γ = 40), and the gap
grows with γ (up to 2e-11 at 137, 2e-10 at 300, 3e-9 at 1000, and 5e-8 to 8e-8 at the top γ of l = 4, 5, 6, in the φ_cut, R = 90 and r0 columns), so it is not a static input such as the star or M_ext, which would be the same at every γ. The crate's own convergence *Coordinator, from #196:* this step does not follow. A plane error that grows as γ and a smallest sine that falls as 1/γ give a gap ∝ γ². #196's direct test (the crate's star and M_ext in `PolarProblem`) supports the conclusion.
is far tighter (E1, E2), and its f64 lane differs from dd by about the size of the gap for the base, R = 90 and r0 columns but not for φ_cut (1.5e-9 against 4.7e-8). **The cause is not established.**
Candidates: scipy's integration error at rtol 1e-10 over the longer path to r_m = 34, and roundoff in scipy's f64. Against truncation alone: scipy's rtol 1e-12 value is still 3.1e-9 from the crate, only 1.7 times closer than
its rtol 1e-10 value (5.4e-9) for a hundred times tighter rtol.

The crate's check numbers for the stored rows alone, before any variant: the worst crate-against-scipy relative differences at the stored γ are in the RA lines of the output (base 2.1e-9, 5.4e-9, 1.9e-8 at l = 4, 5, 6).

**Timings** (one process, one thread, nice 10, load 1.9 to 3.0, Apple M5; the crate in dd; the evaluations of every (γ, variant) of the stored rows, setup 3.3 to 3.5 s per l extra):

| l | evaluations | crate (s) | stored time column, summed (s) |
|---|---|---|---|
| 2 | 45 | 16 | 462 |
| 3 | 50 | 45 | no stored file |
| 4 | 55 | 134 | 36810 |
| 5 | 55 | 221 | 51614 |
| 6 | 55 | 400 | 91606 |

The stored column is the sum of worker wall times on the win laptop (below-normal priority, other work running), and covers the evaluation of q and the capped rows differently, so the two are not like for
like; #186 measured scipy on this Mac 4 to 5 times faster than the stored times. Peak resident: 57 MiB per rows run, 114 MiB for the default mode (reported only).

**The controls** (each run alone, exit 1; the values are in `run_08` to `run_17`): `qint_cofsign` (QP 9.0, RA2 63 of 66), `qint_f` (QP 1.3e-3, RA2 58 of 66), `rm_off1` (RP 4 of 4 differ, RM 5 of 5, RA2 33 of 66), `rm_cut`
(RM φ_cut 30.82 against 34.21, RA2 57 of 66), `r0_grid` (RM r0 24.1098 against 24.1110, RA2 58 of 66), `qstart` (RA2 54 of 66), `moncrief` (RA2 35 of 66, RX −3.7e-6, QP 0.26), `tol_loose` at tol 1e-1 (CAP5b E1 = 1.4e-3; the `tol`
entry of CAP5c and of CAP5s 1.4e-3, run_15 and run_17), `jsign` (CAP5a: D > 0 in base, rtol and r0).

### Defects and discrepancies found

- **CAP5c and CAP5s not met** (above).
- **The control `qint_sign` could not fire** (a sign of S does not change the null vector); replaced by `qint_cofsign` before any run.
- **`qstart` does not fire on RX**: at γ = 300 the Riccati forgets its start; RX is not its target (`dev_31`).
- **`tol_loose` at tol 1e-3 did not fire** (E1 3.9e-6 against the bar 2.1e-4, `dev_33`); the control runs at 1e-1. A sweep run at tol = 1 did not finish in 8 minutes (the integrator at that tolerance
  had not finished) and was stopped by its PID with its script (`dev_34`).
- **My skip fix missed S2's numpy bool** (above); found by the run of record of the default mode and rerun.
- **`convergence.py`'s spread leaves out the base variant** (its name contains `phi_cut`; #114 noted it). The crate's rows keep it so that the stored spreads reproduce.
- The R = 90 grid index is the least certain of the five (ln φ/φ_cut 7.8e-5 from flipping): a φ error of that size at r = 24 would change r_m by one grid point (about 3e-4 in D at γ ≥ 1).

### Recommended tasks

1. **Compare #114's l = 6, γ = 12000 rtol 1e-12 value with the crate's five D.** What: when #114 finishes, a one-line comparison of its D with the crate's −1.59579128004e-13 (rtol analogue) and of the other four
   variants with its stored values. Why: the second source for the row that took 89 h. Accept: a bar set before the value is read (at l = 5 the same difference was 3.1e-9, and the crate's f64 lane differs from dd by 6.9e-9 here); reported as it is. Depends: #114.
2. **Where scipy's D drifts from the crate with γ.** What: scipy's φ_cut and base D at (5, 1000) to (5, 3000) at rtol 1e-10, 1e-11, 1e-12, 1e-13, and its D with the crate's M_ext and star substituted. Why: CAP5c is not met and its cause is open: truncation,
   roundoff or an input. Accept: the gap against rtol curve at each γ, with a mutated input that moves it. Depends: none.
3. **Complex s** (#188), the engine for the argument-principle counts. Accept as in #186's recommended task 3, with the rows' q_int as a second output. Depends: this issue.
4. **Quad-double at (5, 6000)** (about 1.5 h): closes the lane question for the capped rows. Accept: agreement with dd to 1e-20. Depends: none; low priority.


### Coordinator's review (#187)

Under review: `wip/i187-taylor-rows` at 2a84083, raw `I187-mac` at 89bc112. Review folder `negative-mass-raw/I187-mac-review`
(README with the criteria, committed before the first run of `ind187.py`). The changed files are those the issue names.

- **Rerun.** run_02 to run_06 and run_18 from a worktree at 2a84083, with a fresh `--jsonl-dir`, so nothing was read
  back from the runner's cache. Each log is identical to the runner's line for line once timings and load are removed, with
  the same exit code (`compare_rerun_1.txt`; run_05 exits 1 on CAP5c and CAP5s).
- **ind187, an independent test of the rows.** It does not format the crate's values as strings. It parses each stored
  token and tests whether the crate's dd value lies within half a unit of its last digit. It computes the spread and the
  signs itself.
  - **K1**: 41 rows, 266 tokens, all met. The closest is `D(r_m~34)` at l = 4, γ = 1, at 0.994 of half a unit. That is
    3e-8 absolute from the rounding boundary, against crate–scipy agreement near 1e-14 at γ ≤ 10, so the equality is
    not luck.
  - **K2**: l = 3, which had no stored file and had not been compared with scipy. At γ = 0.3 and 40, D, the φ_cut D,
    q_int and q_ext agree with scipy to 1.6e-13 … 2.5e-11.
  - **K3**: CAP5c's five numbers, recomputed from #114's file: 5.44e-9, 3.11e-9, 4.70e-8, 7.98e-10, 4.16e-9.
  - Each of the mutations `k1-ulp`, `k2-l` and `k3-pair` fails exactly its target.
- **CAP5c stays not met**, as reported. Its cause is #196. The quad-double lane at (5, 6000) is #197. The comparison with
  #114's l = 6 value, once that run ends, is #198, with its bar set there before the value is read.
- *Added by the coordinator:* the note's "closest rounding boundary" is relative to the value. Measured in half units of
  the printed digit, the closest token is the one K1 names above, not q_ext at l = 2.
- *Added by the coordinator:* #189's runner found that `plane.rs`'s doc comment says "classical Gram–Schmidt with one
  reorthogonalization", but the loop is modified Gram–Schmidt applied twice. The comment will be corrected when #189
  merges, so that #189's branch is not changed under it.
- The rerun's timings are within 4% of the runner's, at load 2.7 to 5.1.

Follow-ups: #196 (scipy's drift with γ), #197 (quad-double at (5, 6000)), #198 (#114's l = 6 value), and #188 (complex s,
already filed).

**Update from #189.** The crate can step its planes with multishoot's plane stepper (pull request #70, unmerged) behind the *Update 2026-11-28:* #70 is merged and released in multishoot 5.0.0; #81, the complex plane of #188, is merged into multishoot main (75d1e1e), not yet released.
cargo feature `lib_plane`. With it, every D and q token of #187's rows is the same; see [taylor_plane_lib.md](taylor_plane_lib.md) and its
review.

## #196: where scipy's D drifts from the crate

Status: **numerical evidence**; cause by the declared rules not established. At (l = 5, γ = 1000 … 3000) the gap between `PolarProblem`'s D and the crate's is not a drift: one extra chunk moves scipy's D by 4e-10 … 2e-8,
and the gap does not fall with rtol (1e-10 … 1e-13) or atol (1e-14 … 1e-17). The crate's star and M_ext in `PolarProblem` move D by no more than that scatter; the exterior (the scalar plane, q), the star's interpolation and the generated J's rounding
do not make it. Scipy's interior plane at r_m differs from the crate's double-double plane by 4e-16 … 4e-14, the exterior plane by 3e-16; D is the product of the sines of the three principal angles between the planes, the smallest 4.8e-5 at γ = 1000 and
1.6e-5 at 3000, so the plane error is amplified 2e4 … 6e4 times. **#187's "the gap grows with γ, so it is not a static input" does not follow**: a plane error that grows as γ and a sine that falls as 1/γ give a gap ∝ γ².
Post hoc, A = ε_Y/sin θ_min at γ = 1000 gives A (γ/1000)² = 5.3e-9 and 1.6e-8 at (5, 6000) against the stored 5.4e-9 and 4.7e-8 (CAP5c), 2.1e-8 and 6.3e-8 at (6, 12000) against 1.9e-8 and 5.1e-8: **CAP5c's bar 1e-8 was below scipy's own floor**.
Scipy with the crate's J tables (J rounded once from double-double) runs 10 times faster (11 s at (5, 1000)). Raw `I196-mac`; script `taylor_drift.py`, note `taylor_drift.md`.

## #188: complex s

Status: **numerical evidence**. The crate evaluates winding2's D (its fixed R_s and its phase strip) at a complex growth rate s on a pair lane
`Cx<T>` of f64, double-double or quad-double. At the 147 stored points of the §5 box (l = 2, Re s ∈ [0.1, 2], |Im s| ≤ 2.5) it agrees with the
scipy path at rtol 1e-13 to 1.2e-12 in modulus and 5.1e-11 rad in phase (bar 1e-9), in both f64 and double-double, and 8 points of an
l = 6 box at R_s = 1293 agree to 5.4e-12 and 8.5e-11 rad. winding2's counter with the crate as its evaluator gives the stored count
(winding 0, planted pair 2) of both boxes, asking for exactly the stored points. **One declared check, XR, is not met** (below). Every number is in
`taylor_complex_output.txt` (the runs of record; raw folder `negative-mass-raw/I188-mac`, README first, design in `design.md`, last commit when
the output was assembled in its first line); numbers named `dev_NN` come from development runs in that folder and are not in the output file. No sentence here claims novelty.

### What was added

- **Feature `complex`** (off by default; exclusive with `lib_plane`, `compile_error!`). `src/complex.rs`: `Cx<T: Num> { re, im }` implements `Num`, so the
  series, the plane stepper, `det_n`, the projector, the Riccati integration, `JTable::to_lane` and `VacSource` run on it unchanged. A `Cx` with im = 0 is
  its real part bit for bit (products with exact zeros; quotient, reciprocal, square root, logarithm take the lane's own operation); `ln`, `exp`, `to_f64`
  of a number with im ≠ 0 panic. `src/cdrive.rs`: `dc` (a list of s) and `dserve` (s on stdin, one JSON line out each). `src/cplanted.rs`: `cp1` (the planted
  complex plane) and `czserve` (the planted determinant). `table.rs` is not touched.
- **Shared code** (every change is the identity on the real lanes; R1 compares the binaries): four trait defaults in `num.rs` (`conj`, `abs2`, `mag`,
  `re_f64`); `gram_schmidt` takes `conj(y_p) · y_c` and `abs2`; nine `.to_f64().abs()` became `.mag()`; `zerilli_q` has `polar_solver`'s branch rule for
  a complex s (`if Re(q/s) > 0 { q = -q }`, never true for real γ > 0); `d_eval` is split into itself and `d_eval_s(cfg, table, mext, s, r_s, r_far)`.
  725 lines are added in 8 files and 16 changed.
- **Gauge.** Gram–Schmidt with a real positive R diagonal, the phase in Q, `ln_vol` real. LAPACK's complex QR also has a real R diagonal (`dev_16`), so
  `polar_solver`'s `ph` is a sign, and Hermitian Gram–Schmidt gives the same D (`dev_07`: 3.3e-12). D = det[Q_Y | Q_Z] · det R_Y · det R_Z with
  det R real positive, so its phase is that of the holomorphic determinant.
- **winding2's normalization**: the Python driver passes the fixed R_s (`rsfix`, from winding2's own kmin) and `zfac`, and applies the phase strip
  exp(−i Im[(k₊+k₋)(R_s − r_m) + (k₊+k₋+s) r_m]) in double precision to the crate's unstripped D.
- **`python taylor_polar.py --complex`**: the checks and controls below; the counter is a copy of winding2's `count()` and `--plant`; `--scipy-ref-make=BOX`
  makes the scipy references; `--list`, `--only=`, `--mutate=`.

### Results of the runs of record

| Check | Result |
|---|---|
| R1, the real lanes: main 0d1d376 against this branch, default and `lib_plane`, 65 records each (D at l = 2, 6, γ = 0.001, 1, 40 in f64, dd, qd; planted; star), every field but the times | 0 differ in both pairs |
| R2, the complex lane at Im s = 0 against the real lane, 19 cases | all equal as strings, Im D = 0, steps, q, r_s, r_far equal |
| P1C, planted complex plane (P1 of #189 with phases 7, −3, 5, −11), Λ = ln_vol + ln det(Q₀ᵀQ) = 384.8047 + 48i | plane f64 1.1e-14, dd 6.8e-32, qd 1.4e-64; Λ 5.5e-13, 6.2e-31, 2.5e-61; orthonormality 2e-16, 9e-33, 1e-64; every one under its bar (1e4 ε, 1e4 ε max(1, \|Λ\|), 1e2 ε) |
| X1, D at 147 stored points of the l = 2 box against scipy at rtol 1e-13 | dd and f64 worst modulus ratio 1.2e-12, worst phase 5.1e-11 rad, 0 over 1e-9 |
| X1b, 8 points of the l = 6 box Re s ∈ [0.03, 0.1], R_s = 1293.18 | worst modulus ratio 5.4e-12, phase 8.5e-11 rad |
| W0, winding2's counter ported, on the stored values alone | winding +0.0000, max step 0.39, plant +2.0000 (0.40), all 147 and 195 stored points asked, none missing: the stored results |
| CX2, the counter with the crate (dd) as the evaluator | b1 and b2: winding 0.0000, max step 0.39, min\|D\| 6.70e-04 and 1.95e-05, plant +2.0000 (0.40); the evaluated points are exactly the stored 147 and 195 |
| XRc, b1 counted with the scalar plane started at rsmult 40 and 80 (R_s 404.9, 785.6) | winding 0, plant +2.0000 at both |
| CXZ, a planted determinant with zeros at 0.75 ± 0.5i, through the stepper | D against mpmath's `expm` at 20 contour points 18 ε·(growth) (bar 1e4); winding +2.0000 in Re s ∈ [0.2, 1.2], 0.0000 in [1.3, 2.0] |
| L1C, qd against dd at 3 points | 1.5e-29 modulus, 7.6e-27 phase (bar 1e-20); f64 against dd 3.4e-14, 6.7e-14 |
| **XR**, the stripped D at R_s = 200, 300, 600 against 404.85 (8 points) | **not met**: modulus 2.2e-16, **phase 0.11 rad** (bar 1e-9) |
| the default mode (#186's checks) with this branch's default binary | 22 passed, 0 failed, 1 skipped (TS) |

**The stored values.** The caches were made at scipy's rtol 1e-8, and their phase is off: against scipy at rtol 1e-13 the stored values of the l = 2 box differ by
6.0e-10 in modulus and **9.8e-7 rad in phase** (up to 3.1e-6 at l = 6), so the issue's 1e-9 cannot hold against them for any accurate evaluator; X1 and X1b
use the rtol 1e-13 reference, the stored comparison is reported without a bar, and the crate against the stored values (9.8e-7, 3.1e-6) is the same number. The
phase error falls like rtol (`dev_14`: successive phase changes 5.8e-7, 5.6e-8, 6.0e-9, 6.5e-10, 7.2e-11 for rtol 1e-8 → 1e-13). **It cannot move a stored
count**: the smallest distance of any quantity the counter compares with its thresholds (0.4, 0.4, 0.1 rad) from them is 2.2e-3 rad on b1 (2.2e-3 with the
plant) and 5.6e-3 (2.9e-3 with the plant) on b2, against a phase error of 1e-6, and CX2 reproduces every stored count and asks for the same points.

**XR and its cause.** winding2 strips the exponential rotation e^((k₊+k₋)(R_s−r_m)); the exterior solution's tails remain: a 1/R_s term from the start
vectors (the eigenvectors of the frozen matrix) and the angular term, and a ln R_s term from the mass. At 8 points the phase relative to R_s = 1200 falls like
1/R_s with M = 0 (`mext0`: 2.6e-2, 1.7e-2, 1.2e-2, 7.1e-3, 4.0e-3, 2.4e-3, 7.9e-4 at R_s = 100 … 900 for l = 2, s = 0.1 + i) and like ln R_s with the mass
(0.16, 0.13, 0.11, 0.081, 0.056, 0.039, 0.016); the modulus agrees to 4e-14 throughout (`run_19_xr_diag`). Both tails are holomorphic and nonzero at a fixed R_s,
which is all the argument principle needs, and XRc confirms the counts. The design's XR premise (that the strip makes D independent of R_s) was wrong; winding2's
own docstring claims only holomorphy up to a zero-free factor at a fixed R_s. XR is reported as not met, not changed.

### Timings (one process, nice 10, load 5 to 7, Apple M5; the crate with one star and one J table per process, setup 5 s per l not included)

| box | R_s | scipy rtol 1e-8 (the stored) | scipy rtol 1e-13 | crate f64 | crate dd |
|---|---|---|---|---|---|
| l = 2, Re s ∈ [0.1, 2], 8 points | 404.9 | 4.59 s (2.0 to 8.4) | 13.6 s | 0.053 s (×87) | 1.16 s (×4) |
| l = 6, Re s ∈ [0.03, 0.1], 8 points | 1293.2 | 12.1 s (6.3 to 23.0) | 38.5 s | 0.163 s (×74) | 2.77 s (×4) |

f64 agrees with scipy at 1e-13 as well as dd does (X1), so f64 is the lane for counts; quad-double takes 88 to 104 s per D at the box's points (L1C). The cost
is the vacuum plane: R_s − r_m steps at about 1.3 per unit, each re-expanding J's series at the lane (`run_18_tc_big_Rs`, l = 2, s = 0.1 + i, f64 / dd per D):
0.059 / 0.744 s at R_s = 405, 0.131 / 1.72 s at 1000, 0.262 / 3.16 s at 2000, 0.560 / 6.33 s at 4000, 1.19 / 12.6 s at 8000; the interior plane takes
0.005 / 0.11 s. The design's estimate for double-double (0.3 to 0.6 s at R_s = 405) was low by 1.5 to 2: the vacuum series, not the planes, dominate.

**What #116's strips take.** The call is the same with `rsfix` = r_m + 40 / kmin; R_s is 1.3 × 10³ for the l = 6 box above and a few 10³ at Re s of order 10⁻⁴ (kmin about 0.01, an estimate), which costs 0.13 to 1.2 s per D in f64 (R_s 1000 to 8000), so a box of
about 150 evaluations is minutes, against 3.7 h on rpi5. The phase strip's argument reaches 6 × 10³ rad at R_s = 4000 (1.2 × 10⁴ at 8000), whose double-precision error is
about 1e-12 rad, below the counter's 0.4 rad steps. No check covers R_s above 1293 against scipy (the l = 6 scipy values there take 20 to 70 s each at rtol 1e-13).

### The controls (each run alone, exit 1; effects as computed in the development runs before they were declared)

| control | target | effect |
|---|---|---|
| `r1_ulp` | R1 | 1 of 65 records differs in each pair |
| `cx_recip_general` (the general complex quotient at im = 0) | R2 | 18 of 19 cases differ |
| `bilinear` (no conjugate: "drops the phase") | P1C | plane 4.9e-15, orthonormality 1.2e-16 and ln_vol's real part 5.5e-13 unchanged (f64); **the phase of Λ off by 2.27 rad** in every lane (on D, in a development run on X1's 8 points: modulus 4.6, phase 3.0) |
| `redot` (the projector's imaginary part dropped) | P1C | orthonormality 1.0 (f64), Λ off by 55; the plane stays |
| `f_scaled` (f₂ × (1 + 1e-6) in the system) | P1C | Im Λ off by 3.60e-5 in every lane |
| `cx_mul_sign` | P1C, X1 | Λ off by 3.2e2 (f64); D modulus off by 9.3 |
| `strip_off`, `conj_s`, `rs_shift` (R_s + 1, strip unchanged), `strip_sign` | X1 (and X1b for the first three) | phase off by 2.4, 2.9, 2.6, 3.1 rad at 6 of 8 points (the two real-axis points do not move) |
| `lowsum` | L1C | qd against dd 1.1e-6 modulus, 3.8e-3 phase |
| #186's `jsign`, `moncrief`, `swapcols`, `mext0`, `star_eps`, `shift_off` | X1 | 8 of 8 points over the bar: modulus 0.99, 39, 3e-13 (phase 3.1), 1.4e-2, 4.4e-3, 9.1e-2 |
| the counter's `reverse`, `wrongdouble`, `misplant`, `onesided` | W0, CX2, XRc | plant difference −2.0000, +1.0000, 0.0000, +1.7195 (b2: +1.4735) against 2 |
| `conj_s`, `reverse`, `wrongdouble` | CXZ | winding −2.0000, −2.0000, +1.0000 against +2; for `conj_s` D is also off from mpmath (6.0e31 in the check's units, bar 1e4) |

XR has no control: it fails on the unmutated crate, so none could serve.

### Defects and discrepancies found

- **XR** (above): a declared check not met because its premise was wrong.
- The stored caches' phase errors of 1e-7 to 3e-6 rad (above); the 1e-9 bar of the issue was moved to the rtol 1e-13 reference with the coordinator's approval.
- The crate's option for the fixed R_s could not be `rs`: `rs` is the star's start radius (first comparison `dev_22`, crashed; renamed `rsfix`).
- The control `bilinear` crashed on `ln` of a complex norm until `ln` took ln |z| under the mutation; the first conjugate-free probe of stage 1 (`dev_01`) shared
  a class attribute with its Hermitian twin and read the same as it.
- The first launch of the runs of record ran every command with exit 127 (zsh does not split `$TP`); kept in `failed_launch_1/`, rerun.
- The crate in double-double is 4 times faster than the stored scipy time, not the 4 to 25 times the design estimated (the vacuum plane's series per step, see the timings); f64 is 74 to 87 times faster.

### Recommended tasks

1. **Recount the winding2 boxes of #5, #54, #115, #117 and #126 with the crate (f64) as a second evaluator.** What: `taylor_polar.py --complex`'s counter on the 31 winding2 boxes of #126's list
   (the Z boxes need task 3) and #117's four. Why: it is the first independent evaluation of those counts and takes minutes per box. Accept: every
   count and planted check equals the stored one, with the thresholds' margin (as W0) printed per box and a control as CX2's. Depends: this issue.
2. **#116 on the crate.** What: the strips next to Re s = 0 at `rsfix` 10³ to 10⁴. Why: 3.7 h per box in scipy. Accept: the count equal to a scipy count on a short strip at rtol 1e-12, f64 and dd.
   Depends: this issue.
3. **Zerilli zeros on the crate.** What: the counts of zeros of Z(r_m; s) (`zerilli_zeros.py`, the poles of D) need Z itself, not q; one more complex linear integration. Why: #121 and #124 call it
   thousands of times. Accept: Z at stored points to 1e-9 against its scipy path at tight rtol. Depends: this issue.
4. **Faster vacuum jets.** What: the vacuum plane re-expands J's series per step; at R_s = 4000 it is 98% (dd) and 99% (f64) of the cost. Why: 2 to 10 times in dd. Accept: D unchanged to 1e-14 relative at 8 points, a control that
   changes the jets. Depends: none.
5. **Move the crate to the library's complex plane** once multishoot's pull request is merged and released, and `Cx<T>::Lib` can be `plane::Complex<T::Lib>` so that `lib_plane` and `complex`
   can be on together. Depends: the pull request below, #199.

### The library pull request (stage 3)

multishoot's plane module takes a complex scalar with a small change: `plane::Complex<R>` (a pair of lane values, a `Real` itself), four defaulted methods on `Real`
(`conj`, `abs2`, `mag`, `re_f64`: the identity, `z * z`, `|x|`, `x` on the real lanes), the Hermitian product in `orthonormalize`, and the control `NoConj`. Pull request https://github.com/nmoinvaz/multishoot/pull/81 (branch
`feat/complex-plane`, 54d6fca after CodeRabbit's four findings, not merged and not published): `tests/plane_complex.rs` runs P1c (this note's planted system) in f64, dd and qd against a closed form with the bars of the plane tests, the real lanes' digests
of P1, P2, P3 are those of #189's runs of record (10 of 10), the complex lane with real data has the real lane's digest, and 16 digests are the same on the default backend (aarch64 NEON) and with `--features portable`;
`cargo fmt` and `cargo clippy --all-targets --features mutate -- -D warnings` (also `portable,mutate`) are clean. `+avx2,+fma` was not run here (no x86 machine); the pull request's CI jobs for AVX2 and ARM64 NEON passed. With the pull request merged and released the crate's `Cx<T>` can
become `plane::Complex<T::Lib>` (recommended task 5).

### Coordinator's review (#188)

Review of branch `wip/i188-complex-s` at 75eb5f7 and raw `I188-mac` at 7b55688d, recorded in negative-mass-raw `I188-mac-review`. On the
Mac, nice 10, at a load of 20 to 60 from other projects.

- **Design.** The design was approved with one added check, X1b (an l = 6 box at large R_s). It was built as designed.
- **Merge.** It was merged after #196. The two conflicts were the command list in `main.rs` and the appended sections here; both are
  kept. #196's `DOut` planes and #188's `d_eval_s` auto-merged. Both builds compile: default, and `--features complex` in `target/complex`.
- **Rerun on the merged tree.**
  - `taylor_polar.py --complex` gives P1C, XRc, L1C, W0, CX2 and CXZ pass. XR fails, as reported.
  - R1, R2, TC and the X1 checks skip without the runner's references.
  - With `--ref-dir` set to the runner's `ref/`: X1 (1.2e-12 in modulus over 147 points), X1-f64 and X1b (5.4e-12 over 8 points at
    R_s = 1293) pass.
  - Harness defect, reported: `--only=X1-f64` prints its PASS line and then raises `KeyError: 'X1-f64'`. The check itself is unaffected.
- **No independent code.** The scipy path is the comparison, and it shares the problem's definition, not the integrator. CXZ's planted
  determinant uses mpmath's `expm`.
- **XR not met.** The stripped D is not independent of R_s: the phase differs by 0.11 rad. The cause is the exterior solution's tails in
  1/R_s and ln R_s, which are holomorphic and zero-free, and XRc shows the counts are unaffected. The finding stands as reported, and the
  bar is not moved.
- **The library PR.** multishoot#81 (`feat/complex-plane`) was merged into multishoot main on 2026-11-28 (75d1e1e). It is not yet in a release.
- **Follow-ups.**
  - Task 1 (recount #126's boxes with the crate as a second evaluator) is #274.
  - Task 2 is added to #116.
  - Tasks 3–5 are not filed yet. Task 5 waits for a release that contains #81.

## #192: J's table from multishoot

Status: **numerical evidence**. The crate takes multishoot 5.0.0 and multidouble 6.1.0 from crates.io, with no git pin, and gives
the same bits as before: 3532 records of three builds and every printed token of the check logs (default mode, `--rows`,
`--rows` with `lib_plane`, `--complex`) are equal, old pin against new pin, except load averages, the timing check TC and one
R1 line (all classified below). Behind the cargo feature `lib_table`, J's table is built by multishoot's `table` module with
#186's settings, and its step list is #186's bit for bit: 464 steps at l = 2 to 6 and 389 for the r0 = 3e-3 range, every
breakpoint and every stored coefficient equal. The 22 default checks and the `--rows` checks give the same verdicts and the
same tokens with the feature on. **The feature waits for a release**: the module is merged into multishoot main and is in no
release, so the optional dependency is a git revision. Raw folder `negative-mass-raw/I192-mac` (README: the criteria of this
stage first, amendments A4 to A13; last commit when assembled: 2cad9bb5). No sentence here claims novelty.

### What was added

| File | Change |
|---|---|
| `Cargo.toml`, `Cargo.lock` | `multishoot = "=5.0.0"` and `multidouble = "6.1.0"` from crates.io, `[patch.crates-io]` removed, multisolver follows to 0.3.0. **The repin forced no source change**: the default, `lib_plane` and `complex` builds compile as they were. The feature `lib_table` enables the optional dependency `multishoot_tables` (package multishoot, git, rev 75d1e1e, library feature `mutate`) |
| `src/table.rs` | the pieces both builds share are helpers (`j_about`, `entry_series`, `linear_rel`, `jsign_flip`); the own `build_table` and `width_of` are `#[cfg(not(feature = "lib_table"))]`; with the feature, `JSource` (a `SeriesSource`) and a `build_table` that calls `Table::build` |
| `src/main.rs` | the `setup` record of `d` gets `"table":"multishoot"` under the feature, nothing without it |
| `taylor_polar.py` | `--build --features=NAME` passes a cargo feature; nothing else |

`JTable`, its lookups (`jets`, with the crate's `taylor_shift` in the lane), `to_lane`, `VacSource` and the plane integrators are
the crate's own in both builds: the library builds the table in double-double and the table is copied into `JTable<DD>`.

**The pin.** multishoot#79 was merged into multishoot main by rebase, so the pull request's head 60d21ab is not on main. The
dependency is `rev = 75d1e1e` (main, where `src/table.rs` is the file of 60d21ab byte for byte; `plane.rs` differs by #81).
5.0.0 predates the module, so **the flag waits for the next multishoot release**; the default build stays on crates.io. Build:
`cargo build --release --features lib_table`, or `python taylor_polar.py --build --features=lib_table`. The builds
`lib_table,lib_plane` and `lib_table,complex` compile and were not run.

**The library call.** `nb` 30 (`NS`), `tail` 1e-25, `last` 2 passed explicitly (the library's default is 4), `fraction` 0.3, no
`hmax`, no `atol`. The source's `valid_to` is the end of the star step. The functions are J0 and J1 = J(s = 1) - J(s = 0) of
each entry of the pattern (two per entry, combined as `a0 + s² a1` before the shift, as before). The pattern is that of the
first step, and every later step is checked to hold nothing outside it (a panic otherwise). The path is
[max(first star step's x0, 0.999 rmin), x1 of the last star step with x0 <= rmax], as #186's loop.

### Results

| Criterion (raw README A4) | Result |
|---|---|
| RP1, builds; no git source in the default tree | **met**: default, `lib_plane` and `complex` build with 5.0.0; `cargo tree --locked -e normal` lists crates.io packages only |
| RP2, 1167 + 1167 + 1198 records of `d` (f64, dd, qd, l = 2, 6; D to 130 digits), `planted`, `starsteps`, `jtable` (l = 2, 6: every coefficient of J0 and J1, 20 digits), `recipcheck` (f64, dd, qd, od), `tabcheck`, `jjets`, and for `complex` `dc` and `cp1`: old pin against new pin, every field as a string except times | **met**: 0 differ in the default, `lib_plane` and `complex` pairs. Control (`mutate=jsign` on the new binary only): 949 of 1167 differ, exit 1 |
| RP3, default mode, old pin against new pin | new pin 22 passed, 1 skipped (TS); old pin 21 passed, 1 skipped, and **R1 failed there by my own `kill` of its child process**. Of the 24 lines (1401 tokens) only R1 and the SUMMARY line that counts it differ. R1 rerun alone on the old tree passes and its 45 tokens equal the new run's |
| RP3, `--rows` | 12 passed, 2 failed (CAP5c, CAP5s, as #187) in both; 1451 tokens, **2 differ: the load average printed in CAP5b and CAP6b** (`load N,`; the rule replaced `load average N` only) |
| RP3, `--rows` with the `lib_plane` binaries (#189's rows) | the same verdicts, 1451 tokens, the same 2 load tokens; the printed E2 and E3 equal (1.42e-24, 3.58e-9; 8.81e-25, 9.60e-9) |
| RP3, `--complex` (13 checks) | 12 passed, 1 failed (XR, as #188) in both; 1034 tokens, **13 differ, all in TC, the timing check** (load average, times per D and their ranges, speed-up ratios); the other 12 check lines are token for token equal |
| the refactor of `build_table` into helpers | bit-neutral: the binary with the code and the feature off against the repin-only binary, 1167 records, 0 differ |
| JT1, `lib_table`, `lib_table,lib_plane`, `lib_table,complex` | **met**: they build |
| JT2, the step list against #186's (every x0, x1 within 1e-9) | **met, with nothing to spare**: 464 against 464 steps (l = 2 to 6) and 389 against 389; all 5418 breakpoints bitwise equal; the coefficients of every step equal as 20-digit strings; head records equal. Control (`mutate=lib_last4`, the library called with `last` 4): 532 and 451 steps, exit 1 |
| JT3, with the feature: default mode, `--rows` | default mode 22 passed, 0 failed, 1 skipped; tokens equal to the repin-only build's log (1463, 0 differ). Rows 12 passed, 2 failed (CAP5c, CAP5s); 1451 tokens, the same 2 load tokens; RA2, RA4, RA5, RA6 pass |
| JT4, `--mutate=wide_table` | **met**: T3 (worst entry relative to itself 26.3, bar 1e-20) and A2 (D of the opposite sign) fail in both builds, and the two logs are equal token for token (193). The step count under it is 87, the star's, in both (464 without) |
| JT5, the table's build time | the same: `t_table` 5.38 to 5.58 s (own) against 5.27 to 5.48 s (library), l = 2 and 6, three runs each; 422 table steps and `linear_defect` (5.124e-25, 3.335e-24) equal. The star's solve and J's series dominate |

### Found, and not fixed in the crate's own `table.rs`

1. **#186's jump to the star step for `h < 1e-14 b` never fired on these inputs** (design, section 9, item 1). The library has no
   such jump and returns `Stall` where a step does not move x, and it gives the same 464 and 389 steps. The jump would hide a
   singularity on the path; multishoot#79's planted case PG shows what the library does there (`Stall` in f64, `NotFinite` in dd and qd).
2. **#186's `width_of` forms `|c_k| h^k` in plain f64** (item 2). The library scales by powers of two, and at nb 30 and r >= 5e-4
   the two give the same bits. The scaled form matters from nb 60 with a pole near 1e-3 (multishoot#79's case PS), not here.
3. **The library's default `last` is 4.** With it J's table has 532 steps instead of 464 (+14.7 %; +15.9 % for the second
   range) at the same `tail`. The crate passes 2 to reproduce #186. A symmetry of period 3 in a J series would need 3.
4. **multidouble 6's reciprocals** (#193, multidouble#188), named as a candidate cause of a difference: none was found.
   `recipcheck` is equal in f64, dd, qd and od, and so is every D to 130 digits in dd and qd.
5. The octo-double `d` record was dropped from RP2: it had not finished in 31 minutes on the old binary alone. The od lane is
   compared through `recipcheck` only.

### Not run, and what the numbers do not say

- `lib_table` with `lib_plane` or `complex` was built, not run. aarch64 NEON only; no x86 build here. The timings were taken at
  a load average of 4 to 33 from other projects.
- The runs compare two builds of one crate on one machine: the repin and the table swap change no printed bit. They do not say
  J's table is right; that is T1, T2 and T3 (against direct evaluation), which pass in both builds.
- The consumer's reach (the plane's steps past the table's width) is multishoot#79's open point, filed as #256. The crate's checks
  pass at `fraction` 0.3 as before.
- Process notes (raw README A5 to A8): the first od record, a scratch tree without `phase2/boson_star.py` (S1 to S3 failed
  there), a jjets parsing fault in my comparison script, and one `kill` of a child I did not mean to stop. None changed a bar.

### Recommended tasks

1. **Take `lib_table` to the release.** What: when multishoot releases a version that contains `table`, replace the git
   dependency by the release, then drop the own `build_table`, `width_of` and the feature. Why: one code path and the library's
   controls. Accept: `jtable` for l = 2 to 6 bitwise equal to the lists stored in the raw folder (`s3_steps_out/`, `run_52`), the 22 default
   checks pass, and `--mutate=wide_table` fails T3 and A2. Depends: the release, and #256 for the plane's cap.
2. **Look up J from the library too.** What: `Table::shifted` and `TableJets` for the lookups, with J0 + s² J1 formed before the
   shift (the library has no parameter), under `lib_plane`. Why: the crate's shift in its own lane duplicates `table::shift`.
   Accept: D equal to this section's rows to the last digit, or the difference explained. Depends: task 1.

### Coordinator's review (#192, stage 3)

Merged at bc59048. Review folder `negative-mass-raw/I192-mac-review/s3` at aca36da9, criteria committed (2eacc699) before the
first run. The changed files are those the brief names.

- **Rerun on the Mac** from a7a57dd: the default and `lib_table` builds succeed, and the default tree names multishoot 5.0.0
  and no git source. `s3_steps.py` on the two new binaries passes, and each of its six step lists is byte-equal (`cmp`) to the
  runner's `s3_steps_out/`.
- **Default mode with `lib_table`:** 21 passed, 0 failed, and every verdict equals the runner's `run_44` except R1. R1 is
  skipped here, because it compares the old-pin and new-pin binaries, which the review tree does not have. Its content is the
  runner's RP2 (3532 records, 0 differ), which was not rerun.
- **Stale text fixed:** the pin passages of this note (§ on the series, #186's patch) and of `taylor_plane_lib.md` now read as
  history, and the split sentence in `taylor_plane_lib.md`'s review is repaired.
- `lib_table`'s git dependency on multishoot main 75d1e1e stays until a release contains `table`; the default build is
  crates.io only.
