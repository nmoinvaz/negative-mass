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
  `recip_dd` fires). multishoot 4.2.0 and multidouble main ba0ad101ce76 (series module) resolve together by `[patch.crates-io]`
  (multishoot's requirement is the caret "5.1.0"); no multishoot change was needed. `mullow` and `Relaxed` are not used: their
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
