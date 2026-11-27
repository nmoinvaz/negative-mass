# #194: the radial boson-star modes and the phantom family on multishoot and multisolver (Rust, `taylor_radial`)

Status: **numerical evidence**. The crate `taylor_radial/` reproduces what it was asked to reproduce for the canonical star and for
the phantom family up to alpha_c = 3.9, and four declared criteria are **not met** or **not computed** (section 6). Every number
below is in `taylor_radial_output.txt` (run_13_final) or in the raw folder `negative-mass-raw/I194-mac` (README first; last
commit when this note was written: e7336a6). No sentence here claims novelty.

| What | Verdict |
|---|---|
| Planted truth (Poeschl-Teller blocks mixed by a variable-coefficient shear), three lanes | passed: f64 7.4e-15, dd 5.1e-31, qd 4.1e-63 on the roots |
| Canonical scan table (13 rows, lowest sigma^2, M, N) | passed: 39 of 39 tokens equal the stored print |
| Crossing, eps = +1 and eps = -1 | passed: sigma^2 = 0, dN/dw = 0, dM/dw = 0 agree to 1.6e-19 and 1.8e-19. Against the stored table: eps = +1 within 3e-8, **eps = -1 not met** (5.11e-8) |
| Phantom rows alpha_c = 1.10 to 4.017906 (C1) | **not met**: 86 of 105 tokens equal; 10 differ by 1.03 to 1.45 half-units; 9 not computed (rows 4.0, 4.015, 4.017906) |
| The family by multisolver's arclength to 4.017906 (C2) | **not met as worded**: arclength to alpha_c = 2.0 only (Amendment A2) |
| Phantom modes (F) | **not met**: 34 of 38 tokens; 3.5, 3.9 not computed under the star gate, 4.0, 4.015 not computed |
| sigma_0^2 at alpha_0 = 2.0 and 3.0 (Z) | reported to 30 digits with the spread over five variants; the independent multishoot eigenvalue solver is **not done** |

## 1. What it is

- **Crate** `taylor_radial/` (Cargo.toml, Cargo.lock, 17 files in src/, 4145 lines). Rust 1.94.0, multishoot 4.2.0, multisolver 0.2.2,
  multidouble 5.1.0. No change to any of the three libraries was needed; no library PR exists.
- **Driver** `taylor_radial.py` builds the crate, regenerates `src/rgen.rs` (`--gen`), runs the checks, runs one stage of the family
  (`--stage=arc|stepA|stepB|stepC|stepS|arcm`), and runs a control with `--mutate=NAME`. `taylor_radial_gen.py` writes `rgen.rs`.
- **Output** `taylor_radial_output.txt`, produced from `phase2/radial` by

      python taylor_radial.py --arc-log=$RAW/run_09_arc.jsonl,$RAW/run_05_stepA.jsonl,$RAW/run_10_stepB.jsonl,$RAW/run_07_stepC.jsonl \
          --supp-log=$RAW/run_11_stepS.jsonl            # RAW=../../../negative-mass-raw/I194-mac

  The stage logs are the runs of record. Without `--arc-log` the driver runs every stage itself, each child capped at 3600 s of wall
  time (about 1 h of stages and 25 min of checks). A stage alone: `python taylor_radial.py --stage=arc --tag=NAME`. The driver exits
  1 when a check fails; run_13_final: 17 passed, 4 failed, 1 skipped.
- **Machine and timing.** Apple M5, 10 cores, 32 GB, macOS; `RAYON_NUM_THREADS=2`, `nice -n 10`, one process at a time; other users'
  jobs were running (load averages 2 to 11, listed). From the `.err` and `.load` files of the raw folder:

| run | rc | wall s | CPU s | peak RSS MB | load at start / end |
|---|---|---|---|---|---|
| run_04_arc (arclength to 2.0, 4 scans) | 0 | 286 | 444 | 501 | 2.7 3.2 4.5 / 2.5 3.0 4.1 |
| run_09_arc (the same with the qd variant) | 0 | 522 | 769 | 437 | 11.0 6.0 4.2 / 4.0 5.5 5.0 |
| run_05_stepA (rows 2.1 to 4.017906, stalls) | 1 | 1189 | 2134 | 598 | 2.1 2.9 3.9 / 2.5 2.9 3.5 |
| run_06_stepB (scans 2.5, 3.0) | 0 | 229 | 374 | 877 | 2.5 2.9 3.5 / 3.3 3.3 3.5 |
| run_07_stepC (rows 3.5, 3.9) | 0 | 669 | 1211 | 352 | 3.3 3.3 3.5 / 3.3 3.7 3.7 |
| run_11_stepS (supplement) | 0 | 779 | 1372 | 290 | 3.5 4.3 4.5 / 4.1 4.2 4.4 |
| run_13_final (the driver) | 1 | 1470 | 1705 | 147 | 2.6 3.2 3.6 / 3.4 3.0 3.0 |

  A star by BVP takes 2.6 s in dd, a D evaluation 0.07 s (dd, omega = 0.94) and 4.4 s in qd (order 30, omega_0, R about 2200);
  S3's qd shooting takes about 300 s per star.

## 2. Method

- **Star.** (P, u, m, Z) = (phi_0^2, phi_0'/phi_0, m, alpha^-2) by multishoot's boundary-value layout (`Equations`, `Shooting`,
  `plan_segments`, `newton`), with the exact regular start at r_s = 1e-8 (leading-order series, error O(r_s^5)). Two layouts:
  k = sqrt(1 - w^2) global with a pin row (the canonical family, the crossing), and alpha_c global with k a constant field (the family in
  alpha_c; this is what `in_global` hands to multisolver). Outward segments at most 2.5/k, tail segments at most 10/k, r_t = r_cut + 40/k
  with r_cut where P falls to 1e-44; right end: u = -k - 1/r_t, Z = 1/(1 - 2m_R/r_t), k' = 0. Newton to a scaled residual 1e-28 with
  iterate scaling (`Scaling::Iterates`; without it Newton stalled at 2e-24 on tail unknowns of size 1e-70). The model is lambda = 10240,
  mu = lambda^2/4.8, with pi, lambda and mu as double-doubles in every lane (one model in all lanes).
- **J.** `rgen.rs` is the syntax-tree translation of `radial_rhs_generated.py` (J0 + sigma^2 J1 for y = (A, A', G~, chi)), written by
  `taylor_radial_gen.py`. **The symbolic source is shared with `radial_modes.py`; it is not independent.** J is stored as Taylor tables
  along the star's steps, order 30.
- **D(sigma^2).** The regular 2-plane (from r_s outward to r_m) and the decaying 2-plane (from r_far inward to r_m) are
  orthonormalized after every step (Gram-Schmidt, positive diagonal, twice); D = det of the two planes at r_m. The decaying plane
  starts from the spectral projector of the frozen matrix (matrix sign function in the lane): an arbitrary start plane gave a
  spurious sign change of D that moved with R (Amendment A1.1). r_far = r_m + 40/kappa_min (dd; f64 19, qd 74).
- **Lanes.** f64, double-double (the default), quad-double (the check lane); the J table carries 30 coefficients in every lane, so the
  qd integrator runs at order 30 on a table that is itself dd-accurate: **the qd lane of D is a lane check of the integrator, not a
  qd star.** Roots: sign changes on the declared grid, Illinois to 1e-9, secant polish.
- **Continuation.** `arclength::continuation` on the alpha_c-pinned layout, chunks of 10 points; after each chunk the star at its end
  is solved again on a fresh layout (new nodes and Taylor plan). Rows at the stated alpha_c by Newton from the nearest path point.

## 3. Where multishoot and multisolver fit, and where they do not

| Part | Used | Note |
|---|---|---|
| The star (both layouts) | multishoot `Equations`, `Shooting`, `newton`, `in_global` | S2: its recursion equals the crate's own `Num`-lane recursion to the 1e-28 bar |
| The family in alpha_c | multisolver `arclength::start`, `continuation` | to alpha_c = 2.0 through the omega-fold; 210 points, 2 refused steps |
| The omega-fold | multisolver detection (sign change of t_k) | `locate_fold` **did not converge** on the k-pinned family (section 6) |
| D(sigma^2) | neither | the planes need an orthonormalization between steps and the layout has no place for it (#186's finding) |
| An independent eigenvalue solver with sigma^2 global | not built | **NOT DONE** |

Nothing is wrong with the libraries that I could show. `par = false` in my own `exact_f` made the corrector serial: with `par = true`
(multishoot documents that the result does not depend on it) one chunk of 10 points gives 7 identical JSON lines, residual norms and 30-digit
rows included (dev_06, dev_07), and took 12.9 s instead of 16.5 s (one measurement, load 8 to 11).

## 4. Planted truth first, through the same code

Two Poeschl-Teller blocks (lambda = 18/5, b = 4, c = 2 and lambda = 27/10, b = 3, c = 1) mixed by a variable-coefficient shear
S = (I + g1 e12)(I + g2 e21), g1 = (2/5)/(1 + (r/8)^2), g2 = 3/10 + r/20; exact s = (b - (lambda - n)^2)/c for odd n: -1.38, 0.11, 1.82.

| Check | f64 | dd | qd | limit (f64, dd, qd) |
|---|---|---|---|---|
| P1 roots, worst relative, 3 of 3 found | 7.44e-15 | 5.14e-31 | 4.08e-63 | 1e-10, 1e-24, 1e-50 |
| P2 D against mpmath (Legendre) at 5 non-eigenvalues | 1.93e-15 | 1.07e-31 | 6.13e-61 | 1e-9, 1e-24, 1e-48 |
| P3 planted extremum (degree 8) | | 6.68e-51 | | 1e-14 |

The orientation of the start plane (the sign of det C) is constant on the grid in every lane (P1 prints it).

## 5. Results

### 5.1 Checks of the star and J

| Check | Result |
|---|---|
| G1, J against `radial_rhs_generated.py` in mpmath at 70 digits, exact rationals, 12 points | 1.61e-30 over 144 entries (limit 1e-24) |
| G2, the table's shifted jets against the direct series, orders 0 to 11, 8 radii | 6.84e-22 (limit 1e-20) |
| S1, the star against `solve_star` (scipy, tol 1e-9) at 10 radii, alpha, a, phi_0, u | 4.15e-8 (eps = +1, w = 0.75), 4.10e-8 (eps = -1, w = 0.7603); limit 1e-7. This is scipy's boundary condition at r_0, not the crate |
| S2, multishoot's recursion against the crate's own | passed (1e-28) |
| S3, the dd boundary-value star against the qd Taylor shooting (alpha = 1 at r = 60, tol 1e-62): P_c, Z_c, M | 3.0e-33, 4.8e-32, 9.2e-32 (eps = +1); 3.5e-32, 4.4e-32, 5.2e-32 (eps = -1); limit 1e-26. **Not a qd boundary-value solution** |
| S4, largest scaled Newton residual of the 13 stars of the window | 6.6e-30 (limit 1e-28) |
| S5, first law dM = eps w dN, 13 stars at spacing 1e-4 at w = 0.93, degree 12 | 1.55e-27 (limit 1e-18) |

### 5.2 The canonical scan table (check A)

All 39 tokens equal the stored print; the largest distance from a token is 0.988 half-units of its last digit (limit 1). Lowest sigma^2:

| w | 0.90 | 0.91 | 0.92 | 0.925 | 0.93 | 0.935 | 0.94 | 0.95 | 0.96 |
|---|---|---|---|---|---|---|---|---|---|
| sigma^2 (4 digits) | +9.462e-3 | +6.205e-3 | +2.376e-3 | +3.448e-4 | -1.685e-3 | -3.638e-3 | -5.431e-3 | -8.177e-3 | -9.236e-3 |

No sign change of D on the grid at w = 0.75, 0.80, 0.85, 0.88. At w = 0.94: sigma^2 = -5.43147383546923512697530554200e-3 (dd).
The stored full-precision values differ by 8.6e-12 there and by 2.1e-7 at w = 0.90 (scipy's error; reported, not judged).

### 5.3 The crossing (X1), by the secant method on D(0; w) and the extremum of the interpolating polynomial of N(w), M(w)

| eps | sigma^2 = 0 | dN/dw = 0 | dM/dw = 0 | the three agree to |
|---|---|---|---|---|
| +1 | 0.9258444557595977278577 | 0.9258444557595977280211 | 0.9258444557595977280332 | 1.6e-19, 1.2e-20 |
| -1 | 0.9289084310581665250429 | 0.9289084310581665252239 | 0.928908431058166525237 | 1.8e-19, 1.3e-20 |

The tolerance set before the run was 1e-12 on both differences: met by seven orders. Error bars, named: **E1** is a spread over the fit window
and degree, not an enclosure: 1.7e-19 (eps = +1), 1.9e-19 (eps = -1). **E3**, the shift of omega_0 under R x 1.5 and under order 24 at tol 1e-22:
1.1e-31 and 7.8e-28 (+1), 5.8e-32 and 6.4e-24 (-1). **E2**, omega_0 with D in qd on the same dd star: 2.5e-27 (+1), 1.1e-28 (-1); omega_N and omega_M in
qd would need a qd star at every omega: not done.

Against the stored table (0.92584444, 0.92584446, 0.92584446 and 0.92890846, 0.92890838, 0.92890838), within 3e-8 each:

| eps | sigma^2 = 0 off by | dN/dw = 0 off by | dM/dw = 0 off by | verdict |
|---|---|---|---|---|
| +1 | 1.58e-8 | 4.24e-9 | 4.24e-9 | met |
| -1 | 2.89e-8 | **5.11e-8** | **5.11e-8** | **not met** |

The stored values are polynomial fits to 17 stars at boundary-value tolerance 1e-9, so the 5.1e-8 is the stored fit's error. The criterion is not reworded.

### 5.4 The phantom family (check C1; stored rows are boundary-value solutions at 1e-8)

13 rows at alpha_c <= 2.0 (arclength path, 210 points, 2 refused steps; the path reaches 2.0184) and 19 distinct rows 2.1 to 3.9 by Newton stepping
in alpha_c (A2). **Rows 4.0, 4.015 and 4.017906: not computed.** The ten tokens that are not equal to the stored token:

| alpha_c | quantity | stored | crate (to the stored digits) | crate (16 digits) | half-units | relative |
|---|---|---|---|---|---|---|
| 2.0 | N | 120.14873 | 120.14872 | 120.1487245048392 | 1.10 | 4.6e-8 |
| 2.1 | N | 141.36645 | 141.36644 | 141.3664447506266 | 1.05 | 3.7e-8 |
| 2.3 | omega | 0.62757001 | 0.62757000 | 0.627570004505126 | 1.10 | 8.8e-9 |
| 2.4 | N | 213.38902 | 213.38901 | 213.3890127494153 | 1.45 | 3.4e-8 |
| 2.4 | M | -116.14250 | -116.14249 | -116.142494219367 | 1.16 | 5.0e-8 |
| 2.8 | omega | 0.73465189 | 0.73465188 | 0.734651884514128 | 1.10 | 7.5e-9 |
| 2.8 | N | 328.15125 | 328.15124 | 328.1512430989620 | 1.38 | 2.1e-8 |
| 3.5 | N | 578.11028 | 578.11027 | 578.1102733974376 | 1.32 | 1.1e-8 |
| 3.9 | omega | 0.97410304 | 0.97410303 | 0.974103033967024 | 1.21 | 6.2e-9 |
| 3.9 | N | 748.13401 | 748.13400 | 748.1340048394089 | 1.03 | 6.9e-9 |

All other 76 computed tokens are equal. The criterion is "each equal to the stored token": **not met** for these ten, each by 1.03 to 1.45 half-units of its last printed digit, and
for the nine tokens of the three rows not computed. The declared information tier (relative 2e-8) is exceeded by five tokens (N at 2.0, 2.1, 2.4, 2.8, M at 2.4).
max |2m/r| and r_s are printed beside every row in the output.

**The stall.** Newton in alpha_c with halving (first step 0.01, at most 0.05, x1.3 after a success, halved after a refusal; a step is taken when the scaled
residual is <= 1e-24 and |d ln k| <= 0.25) stalled at **alpha_c = 3.9935489** (target 4.0; k = 0.1038, omega = 0.99460): the residual stays at 2e-6 to 8e-6 for every step
down to 1e-7, with a first iterate of 1.5. The stored family is smooth there (no turning point in alpha_c). The rule of the ruling was applied: stop, report, no other method.
362 solves, 106 refused (run_05_stepA). The accepted rows' scaled residuals are not a clean function of alpha_c (2.2e-27 at 3.0, 1.5e-29 to 3.9e-29 at 3.1 to 3.4), and rise from 6.6e-29 at 3.5 through 2.7e-28, 1.3e-27, 5.5e-27 to 7.3e-27 at 3.9.

### 5.5 The omega-fold

multisolver's arclength in k detects the fold (a sign change of t_k between k = 0.9064748256 and 0.9064748213, 160 path points); the family in alpha_c crosses it in
chunk 4 (alpha_c 1.1443 to 1.1834; the path minimum of omega is 0.42226152 at alpha_c 1.17311, a path point). `locate_fold` **did not converge** at s = 1.2877e-2 (not located; Amendment A1.7). The
polynomial of omega(alpha_c): **omega_min = 0.42225985865266326181 at alpha_c = 1.1718557760503** (window 1.6e-2, degree 8; spread to the 3e-2 window, degree 10: 4.3e-18 in omega_min,
1.4e-14 in alpha_c, a spread over fit choices, not an enclosure). Stored: omega_min = 0.42226 at alpha_0 = 1.172 (a fit in `phantom_checks.py` C4).

### 5.6 Phantom modes (F) and sigma_0^2 (Z)

Roots of D on [-0.2, (1 - w)^2) (declared grid and a 4 x denser one: the same roots); the stored tables list at most three roots, the crate finds more where noted.

| alpha_0 | omega | (1 - w)^2 | roots (5 digits) | further roots (10 digits) | verdict |
|---|---|---|---|---|---|
| 1.10 | 0.431253 | 3.235e-1 | 3.5940e-2, 1.4917e-1 | | equal |
| 1.20 | 0.422999 | 3.329e-1 | 1.9926e-2, 7.6515e-2, 1.7287e-1 | 0.3096768717 | equal |
| 1.50 | 0.465340 | 2.859e-1 | 1.2704e-2, 4.2725e-2, 9.3272e-2 | 0.1647328544, 0.2573865297 | equal |
| 2.00 | 0.564546 | 1.896e-1 | 1.2843e-2, 3.7466e-2, 7.8791e-2 | 0.1369650961 | equal |
| 2.50 | 0.670180 | 1.088e-1 | 1.5289e-2, 4.0335e-2, 8.2368e-2 | | equal |
| 3.00 | 0.777898 | 4.933e-2 | 1.8685e-2, 4.5649e-2 | | equal |
| 3.50, 3.90 | | | | | **not computed** under A1.8's gate (the star's residual 1.02e-28 and 7.6e-27 against 1e-28) |
| 4.00, 4.015 | | | | | **not computed** (the stall above) |

F: 34 of 38 tokens. **Supplementary run, outside A1.8's gate (Amendment A3), no verdict:** on the stars as Newton converged them, D has no sign change on the grid or on the 4 x denser grid at 3.5 (star
residual 1.0e-28) and at 3.9 (7.6e-27). F stays not computed there.

**sigma_0^2, the lowest root, in dd** (E3 variants: R x 1.5, order 24 at tol 1e-22, layout dmax 6, f64 lane, and E2: D in qd on the same dd star; relative differences from the dd value):

| alpha_0 | sigma_0^2 | R x 1.5 | order 24, 1e-22 | layout dmax 6 | f64 lane | qd lane (E2) |
|---|---|---|---|---|---|---|
| 2.0 | 1.28430640508399647929969563322e-2 | 0 | 4.2e-24 | 6.2e-28 | 8.4e-15 | 5.9e-30 |
| 3.0 | 1.86850874779842754263187486000e-2 | 0 | 1.4e-25 | 8.0e-29 | 2.0e-14 | 1.1e-28 |

These are spreads over choices of the method, not rigorous enclosures. The qd lane agrees to 1e-28 or better (limit 1e-20), which is the dd star's accuracy: the star itself is not in qd. The stored double values
(1.284306405866834e-2 and 1.868508754613376e-2) are off by -6.1e-10 and -3.6e-9. `phantom_fd.py`'s extrapolated finite-difference values were 6.0e-3 and 1.04e-2 away from the crate's values.
The third declared error bar, the independent multishoot eigenvalue solver, is **NOT DONE**.

## 6. Not met, not computed, not done

| Item | Status | Why, with the numbers |
|---|---|---|
| X1b, eps = -1, stored dN/dw = 0 and dM/dw = 0 | not met | off by 5.11e-8 against 3e-8: the stored 17-star fits' error |
| C1, rows against the stored tokens | not met | 10 tokens differ by 1.03 to 1.45 half-units (table 5.4); the stored rows are 1e-8 boundary-value solutions |
| C1, rows 4.0, 4.015, 4.017906 | not computed | stall at alpha_c = 3.9935489, k = 0.1038 |
| C2, the arclength continuation to 4.017906 as worded | not met | arclength to 2.0184 only: the cost, 1370 points to 3.736 (dev_02: 2168 s), 3.0 s per point under load (dev_05); the projection to 4.017906 is 1900 to 3600 points, 0.9 to 3 h; run_03 was stopped at 2 h 29 min |
| F at 3.5, 3.9 | not computed | A1.8's gate (Newton to 1e-28): the stars reached 1.02e-28 and 7.6e-27 |
| F at 4.0, 4.015 | not computed | the stall |
| `locate_fold` on the k-pinned family | not located | the corrector did not converge at s = 1.2877e-2, tolerance 1e-23 (A1.7); the fold is detected and omega_min is from the polynomial |
| E2 for omega_N, omega_M, and a qd boundary-value star | not done | needs a qd star at every omega (the table carries 30 coefficients; S3 compares one star by shooting) |
| The independent multishoot eigenvalue solver (sigma^2 the global) | **NOT DONE** | not built |

**Where the work went wrong, so that it is not repeated:** run_03 asked for one arclength run to 4.017906 under a 45 to 60 min estimate that ignored dev_02's 37 min to 3.736; the driver held the child's
output in a pipe (nothing visible for 2 h 29 min) and its timeout was 3 h. G1 ran at 4 points where the README says 12 (run_08), and `:g` formatting turned the target 4.017906 into 4.01791.
My step test also required Newton's own flag (tolerance 1e-28) beyond A2's residual 1e-24 (dev_09). All four are fixed in the code and noted in the raw README.

## 7. Hypotheses, untested (no run spent on them)

1. **The arclength step length.** The step stays near 1e-3 in alpha_c per point (d ln k per point 6e-4 on the weak side, 1.6e-3 at k = 0.35). The scaled step norm
   (`Scaling::Iterates`, needed so that Newton sees the tail unknowns of size 1e-44 and below) may let the tail, ln P about -2kr over 44 decades, set the step. Not tested.
2. **The floor at k < 0.104.** The layout's outer radius R = r_t against the decay length 1/k = 9.6 at the stall, or the accuracy of the tail condition u = -k - 1/r_t (no O(1/r^2)
   or gravitational correction) at that radius. r_t of the rows is 377.5 (alpha_c = 3.7, k = 0.367), 425.0 (3.8, 0.305), 514 to 516 (3.9, 0.226); the same layout rule
   (r_cut = r_x + ln(P(r_x)/1e-44)/(2k), r_t = r_cut + 40/k) extrapolated from these three rows gives R of about 930 to 980 at the stall (an extrapolation, not a measurement), against
   R = 1919 at k = 0.0109 in `phantom_follow.py`. The rows' residuals rise from 6.6e-29 at 3.5 to 7.3e-27 at 3.9, and the first iterate is 1.5 at a step of 1e-7, so the layout rebuilt from the previous star's
   profile does not reproduce it at this k. `RESULTS.md` section 7 says the old solver's stop near omega = 1 is a limit of the boundary condition, not of the star; the crate neither confirms nor refutes this.

## 8. What `RESULTS.md` should now say

- **Section 3, the crossing (lines 83 to 91 and 158 to 166).** For eps = +1: omega_0 = 0.92584445575959772786, omega_N = 0.92584445575959772802, omega_M = 0.92584445575959772803, equal to
  1.6e-19 (E1 spread 1.7e-19); the stored fits were off by 1.6e-8, 4.2e-9, 4.2e-9, so "coincides to 2e-8" describes the stored fits, not the star. For eps = -1: 0.92890843105816652504,
  ...52522, ...52524, equal to 1.8e-19; the stored dN/dw = 0 and dM/dw = 0 (0.92890838) are off by 5.1e-8, so "coincides to 8e-8" is the stored fit's error.
- **Section 3, w = 0.94 (line 104).** The crate gives sigma^2 = -5.43147383546923512697530554200e-3, equal to every stored digit; the finite-difference value -5.449e-3 is off by 3.2e-3 relative:
  the finite-difference discretization, not the root.
- **Section 7, sigma_0^2.** At alpha_0 = 2.0 and 3.0: 1.28430640508399647929969563322e-2 and 1.86850874779842754263187486000e-2. The stored doubles are off by 6.1e-10 and 3.6e-9; the
  finite-difference extrapolations were 6.0e-3 and 1.04e-2 away. A sentence saying that the FD check confirms the root to a few percent should go.
- **Section 7, the table and the stop.** The stored rows are reproduced to 9e-9 in omega and 5e-8 in N and M (10 of 105 tokens differ in the last digit); rows 4.0, 4.015, 4.017906 and the scans at
  3.5 to 4.015 were not reproduced (the crate stalls at alpha_c = 3.9935489). The statement about the boundary condition at omega = 1 stands as the stored solver's, with the crate's stall as an
  open, different observation (hypothesis 2). The "none at alpha_0 >= 3.5" is consistent with the supplementary scans at 3.5 and 3.9 (outside the gate), not confirmed by a verdict.
- **Section 6, caveats.** The stored root lists are capped at three: the crate finds a fourth root at alpha_0 = 1.20 (0.3097), 1.50 (0.1647, 0.2574) and 2.00 (0.1370).
  The stored omega_min = 0.42226358 at alpha_0 = 1.172 is a fit: the polynomial of the crate's family gives 0.42225985865266326181 at alpha_c = 1.1718557760503.

## 9. Controls (`--mutate=NAME`; targets from development runs, Amendment A4 of the raw README)

Each control fired or stayed silent on its declared checks in the runs of record (`run_ctl_*`, `run_ctlf_*`), as in the development runs. Silent is a result.
`Cs` and `Fs` are C1 and F restricted to alpha_c <= 1.30 and alpha_0 = 1.10, on the stage `arcm` run under the mutation (the unmutated `arcm` passes both).

| control | fires on | silent on |
|---|---|---|
| `jsign` (a J entry's sign) | A, X1a, Fs | G1 (the entry is flipped in the table, not in `jeval`'s series), Cs |
| `jomega` (omega in J off by 1e-6) | A, X1a, Fs | Cs |
| `jmu` (mu in J times 1 + 1e-6) | X1a | |
| `eps` (eps = +1 in the phantom star) | S1, S5, Cs, Fs (the arclength stage never reaches 1.10: stopped at its 600 s cap) | S4 |
| `noqr` (no orthonormalization), f64 lane | P1, P2, A | the dd and qd lanes never finish (1800 s cap) |
| `lowsum` (Taylor sum cut at order 8), f64 lane | P1, P2 | G2, A |
| `shortR` (R - r_m = 3/kappa), qd lane | P1 | Cs, Fs (the scan at 1.10); Z at 2.0 and 3.0 was not run under it |
| `start` (wrong leading term of chi) | | A, X1a, X1b (eps = +1) |
| `plant_sign` | P1, P2 | |
| `fit_deg` (degree 3) | P3, X1a | X1b (eps = +1) |
| `tailk` (u(R) = -k + 1/R) | | S1, S4, S5, Cs, Fs (silent by construction: the star does not feel its far condition; replaced by `znorm`) |
| `znorm` (Z(R) = 1) | S1, S5, Cs, Fs | S4 |
| `nonewton` (one Newton step) | Cs, Fs | S4, S5 |
| `arcfold` (tangent flipped at the first point) | Cs, Fs | |
| `arbstart` (arbitrary decaying plane) | A (374.6 half-units at the worst token) | |
| `start_term` (wrong leading coefficient of u) | | S1, S2, A |
| `shift_off`, `wide_table` (table shift, step width) | G2, A | |
| `recip_dd` (reciprocals in dd inside qd), qd lane | P1 (1.14e-33 against 1e-50) | |

X1b for eps = -1 fails unmutated, so it is not read as an effect. 16 of the 19 controls fire on at least one check; `start`, `tailk` and `start_term` are silent on everything declared for them.

## 10. Recommended tasks

| Title | What | Why | Failable acceptance | Depends on |
|---|---|---|---|---|
| Phantom endpoint: the star solver at k < 0.104 | Diagnose the stall at alpha_c = 3.9935489 (outer radius, tail condition, residual floor) with a star checkpoint (dump and restart) so that an experiment starts at 3.9; then rows 4.0, 4.015, 4.017906 and the scans at 3.5 to 4.015 | The three end rows and the scans test the far boundary | rows 4.0, 4.015, 4.017906 equal the stored tokens or differ by a stated 1e-8 bound; the scans at 3.5 to 4.015 under a gate declared before the run show no sign change | this crate |
| A star gate that scales with the star | Replace "Newton to 1e-28" by a gate that holds for stars with k < 0.4 (the rows' residuals rise from 6.6e-29 at 3.5 to 7.3e-27 at 3.9), or a change of unknown (ln P) | F at 3.5 and 3.9 is not computed only because of the gate | scans at 3.5 and 3.9 inside the declared gate | this crate |
| Independent eigenvalue solver on multishoot | sigma^2 the global, the mode normalized by an extra row, on the layout | the third error bar of Z; shares J and the star, none of the planes, matching or root grid | sigma_0^2 at 2.0 and 3.0 equal the crate's to 1e-18 | this crate |
| A qd boundary-value star | the star BVP in the Quad lane (4-part data) and the J table at order 56 | E2 for omega_N, omega_M; S3 compares one shot | omega_N, omega_M dd against qd within 1e-20 | this crate |
| The step-length hypothesis | Run the arclength with the tail unknowns excluded from the step norm (or with ln P as the unknown) | cost of the family: 3.0 s per point and about 1e-3 in alpha_c per point | points to alpha_c = 3.7 fewer than 1370 at the same tolerances, rows unchanged | this crate |
| `RESULTS.md` sections 3, 6, 7 | the edits of section 8 | the stored statements of "2e-8" and "a few percent" are the stored solvers' errors | each quoted number is in `taylor_radial_output.txt` | this note |

## Files

`taylor_radial/` (crate), `taylor_radial.py` (driver), `taylor_radial_gen.py` (J generator), `taylor_radial_output.txt` (run_13_final), this note. Raw logs, stage logs
(`*.jsonl`), the amendments A1 to A4 and the index of runs: `negative-mass-raw/I194-mac/README.md`.

## Coordinator's review

Review of branch `wip/i194-radial-crate` at 2ed755e and raw `I194-mac` at e7336a6, recorded in negative-mass-raw
`I194-mac-review` (criteria 849e2b7; amendments A1 and A2 before runs 2 and 3). On the Mac, nice 10.

- **Scope.** The branch adds the crate, the driver, the generator, the output and this note under `phase2/radial/`, nothing
  else. No library change.
- **Rerun.** `taylor_radial.py` over the stage logs of record (run_13_final) gives rc 1 with 17 passed, 4 failed and 1
  skipped. **All 820 lines equal `taylor_radial_output.txt`** once wall times and loads are removed.
- **Independent check** (`ind194.py` in the raw folder). The coordinator's own shooting code, in float64 with DOP853 and a
  numba right-hand side, shares only the field equations of `boson_star.py`'s docstring.
  - **The crossing.** Over a ±4 % window (K1) its dM/dω = 0 lay 4.7e-9 above the crate's at both signs of ε. K1 is not met:
    the cause was the window's fifth-order content. Over ±0.5 % (K1b) at rtol 1e-11 it gives 0.9258444559483 and
    0.9289084310226, **+1.9e-10 and −3.6e-11 from the crate**.
  - K1b is not met as worded, because its two tolerances differ by about 2e-9: the narrow window amplifies the rtol 1e-10
    error. So the crate's crossing is confirmed to 2e-10, and the stored fits' 4.2e-9 and 5.1e-8 are the stored fits' errors.
  - **The first law.** dM = εω dN along the canonical family holds to 1.3e-12 at ω = 0.95 (K0). The points at ω = 0.90
    and for ε = −1 were not computed before the 2 h stop, so K0 is not met as worded.
  - **Not computed.** The strong-field phantom rows, thin-walled stars of radius 80 to 170, cannot be held by shooting
    from the centre in float64 (A1). So C1's ten last-digit differences are not arbitrated by the review. The crate's own
    S3 check, dd BVP against qd shooting at 9e-32, is the evidence for the crate there.
- **Applied.**
  - RESULTS.md carries updates from #194: the crossings, ω = 0.94, σ₀² at α₀ = 2.0 and 3.0, ω_min, the fourth roots, and
    the rows and the stall.
  - Paper II quotes the crossings as 0.9258444558 and 0.9289084311, coinciding with the extrema to 2×10⁻¹⁹, and names
    this code in its method.
  - The note's proposal to delete RESULTS.md's sentence on the finite-difference check was not taken. The sentence
    describes what that check can show, and it stays beside the crate's values.
- **Follow-ups.**
  - #233: the end of the phantom family, the stall and the star gate (tasks 1 and 2).
  - #234: an eigenvalue solver independent of the planes, and a qd star (tasks 3 and 4).
  - #235: the arclength step norm (task 5).
