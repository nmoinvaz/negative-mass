# A singularity locator for Taylor coefficients, in multisolver

November 2026, issue [#191](https://github.com/nmoinvaz/negative-mass/issues/191). Pull request
[nmoinvaz/multisolver#26](https://github.com/nmoinvaz/multisolver/pull/26) (module `multisolver::locate`, branch
`feat/locate-singularity`, commits a2d767e and 0db3966 on multisolver main 8e699ab). Script
[`singularity_locator.py`](singularity_locator.py), its stdout [`singularity_locator_output.txt`](singularity_locator_output.txt).
Raw folder: negative-mass-raw `I191-rpi5` (README with the criteria and amendments A0 to A6, every log). It generalizes the
locators of [`conical_global.md`](conical_global.md) (#167) and [`conical_offgrid.md`](conical_offgrid.md) (#183) and uses
#184's series in [`conical_loglog.md`](conical_loglog.md) only through #183's stored points. Status words as in AGENTS.md:
**numerical evidence** for the planted and the #183 results; nothing is proved or certified. Every number quoted below is in
the output file or in the logs named. No claim is made that the method is new; the closest prior work has not been read.

## Answer

The library routine exists, is tested on planted series in f64, double-double and quad-double, and its error bar covers the
planted error in every cell. **68 of the 69 gated cells meet all four criteria; one does not: D6 in f64.** There the locator
says `Unresolved` with a bar of 1.1 |s| (criterion 3, that the bar covers the error, holds): it gives no confident wrong answer. A
window starting near 0.8 K would resolve that regime; it is not implemented.

| question | answer |
|---|---|
| which library | multisolver, module `locate`: `singularity` already holds the estimators it generalizes, and multishoot depends on multisolver |
| method for a real singularity | the ratio method of Domb and Sykes, c_k/c_(k-1) as a polynomial in 1/k, fitted over a grid of windows and degrees (degree 1 to 9, first order 0.2 K to 0.6 K) |
| method for a pair | **a three-term recurrence with terms in 1/k** (not Mercer and Roberts), exact for a conjugate pair or two real singularities of one exponent; the exponent comes from the first-order terms |
| the error bar | **an estimator range over the fit choices** (2 times the largest disagreement of the eight neighboring cells), floored at a few roundings; **not a bound** |
| does it cover the planted error | yes in 69 of 69 gated cells and 72 of 72 rows, for s and for alpha; the bar is 20, 85 and 200 times the error at the median in f64, double-double and quad-double (smallest 1.7) |
| real against pair | right in every cell: 18 pair rows `Pair`, 53 others `Real`, one `Unresolved` (D6 f64); a pair 5.7 degrees from the real axis is resolved to 1.5e-14 in f64 |
| the fit window | a window starting at order 2 fails on the cases with a farther singularity, one of the last 8 orders fails on the noisy ones |
| #183 | at its stored last steps the locator reproduces eta* to at most 0.5 f64 ulp in 31 records of the O8 regime and alpha to 2e-16 of 2/(eta*^2 + 4); two of six kappa-curve records are complex pairs |

## 1. The locator

Near its nearest singularity s a function behaves as A (1 - u/s)^alpha, u the distance from the expansion point, so
c_k = A binom(alpha, k) (-1/s)^k ~ s^-k k^(-alpha-1). `locate(c, opts)` takes c_0 ... c_K (K at least 12; in a variable scaled so
that the radius is of order one) and fits two models over a grid of cells (first order of the window 0.2, 0.3, 0.4, 0.5, 0.6 K;
degree in 1/k 1 to 9 for the ratio model, 2 to 9 for the three-term model), by least squares in the precision of the coefficients
(Householder QR; a rank-deficient fit gives no cell).

- **Ratio** (Domb and Sykes, Proc. R. Soc. Lond. A 240 (1957) 214, already cited in `singularity`): r_k = c_k/c_(k-1) = b_0 + b_1/k + ...;
  s = 1/b_0, alpha = -1 - b_1/b_0. For (1 - u/s)^alpha it is exact at degree 1, and a conjugate pair makes the ratios oscillate.
- **Three-term**: (1 + 1/k) c_(k+1) = (pi_0 + pi_1/k + ...) c_k - (theta_0 + theta_1/k + ...) c_(k-1). For two singularities
  of one exponent (a pair s, conj(s), or two real ones) it is exact with pi_0 = 1/s_1 + 1/s_2, theta_0 = 1/(s_1 s_2), pi_1 = -alpha pi_0,
  theta_1 = -(1 + 2 alpha) theta_0 (derived by hand, and confirmed by the planted pure cases C1 to C3 and T3: the error of s is at most 1.5e-14
  in f64, 1.8e-30 in double-double, 1.4e-63 in quad-double, in the output). Its roots z of z^2 - P_0 z + Q_0 give 1/s: complex, a pair; real, the nearest of two
  real singularities and the farther one. A single real singularity makes the fit rank deficient, and the ratio model is chosen.
  For two singularities of different exponents the first-order terms still give both (P_1 - Q_1 z_j = g_j (z_j - z_other)): case T1 and T2.
- **Choice and bar.** The cell whose eight neighbors differ least from it in s is the estimate; the bar of s is twice their
  largest difference, the bar of alpha twice theirs, floored at 8 u |s| and 64 u max(1, |alpha|) (u = 2^-52, 2^-104, 2^-208). Of the
  two models the smaller difference wins. A bar above 5% of |s| is `Kind::Unresolved`. The bar is a range of the estimator over fit
  choices. It would miss where all neighbors share an error that the windows do not reveal: that is what the planted suite measures.
- The code is generic over `Real` (f64, `DD`, `QD`, `OD`), uses IEEE operations only (no libm), and gives bitwise the same results on the
  default build and with `--features portable` (fingerprint 0xa5113e9ccb6776c8 on aarch64 NEON and on portable; the PR's CI passes it on x86 SSE2, AVX2 and portable, Windows and macOS).

Compared with `choptuik-large-d/hp/branch_radius.py` (read, not run): that script fits ln |c_k| over the orders above a noise floor
for Chebyshev coefficients and reads the Bernstein ellipse (a radius, no exponent, no kind, no bar); the locator reads the algebraic
singularity from the ratios and carries the fit's disagreement. They answer different questions and share no code. The walkers of #167 and #183
fit the last 10 (double-double) or 33 (quad-double) ratios in f64 with one fixed window and degree 2; the locator chooses the window and the degree
and works in the precision of the series.

## 2. Planted results

Lanes f64, double-double, quad-double with K = 30, 60, 120. 24 cases (6 pure single, 6 dressed by poles at +-4i|s|, 3 pure pairs and 2
dressed, 3 sums of two singularities, 3 noisy, one dressed and noisy, reported), 69 gated cells and 3 reported. Criteria fixed in the raw README
before the first library run: (1) the kind, (2) the error of s and of alpha within the tolerance of the cell's tier, (3) the bar covers the error
(factor 1), (4) the bar is within 100 tolerances. Tolerances: pure 1e4 u (s) and 1e6 u (alpha), pure pair and two poles 1e5 u and 1e7 u, dressed 1e-2, 1e-6, 1e-9
(s) and 0.1, 1e-4, 1e-7 (alpha) per lane, noise 2 delta and 70 delta (single), 30 delta and 1300 delta (pair), the last four ten times the amplification of the noise measured in
dev_03b.

| lane | cells meeting 1 to 4 | largest error / tolerance (s, alpha) | bar / error of s: median (min, max) |
|---|---|---|---|
| f64 | 22 of 23 | 0.15, 0.09 | 20 (1.7, 576) |
| double-double | 23 of 23 | 0.07, 0.056 | 85 (3.2, 3.4e4) |
| quad-double | 23 of 23 | 0.028, 0.042 | 202 (5.6, 1.5e8) |

Examples of the cells (error of s relative, error of alpha, in the output file; "pure" are exact at the model's degree, "dressed" are limited by the algebraic
truncation in 1/k of the poles at 4|s|):

| case | f64 | double-double | quad-double |
|---|---|---|---|
| P3 pure alpha = 1/2 | 0, 8.9e-16 | 4.6e-33, 1.8e-31 | 7.6e-65, 3.4e-63 |
| P6 pure alpha = 1e-8 | 2.2e-16, 2.4e-15 | 0, 4.1e-31 | 9.4e-65, 4.8e-63 |
| C3 pair 5.7 degrees from the real axis | 1.5e-14, 4.0e-13 | 1.8e-30, 1.0e-28 | 1.4e-63, 2.4e-62 |
| D3 dressed alpha = 1/2 | 1.5e-5, 1.4e-3 | 1.8e-12, 7.8e-10 | 7.2e-17, 6.4e-14 |
| C5 dressed pair, alpha = 1/2 | 1.6e-4, 1.8e-3 | 1.7e-11, 5.1e-9 | 7.4e-15, 2.4e-12 |
| T2 two singularities on one side (ratio 1.4) | 1.8e-5, 4.6e-3 | 6.0e-15, 2.7e-12 | 1.8e-22, 2.3e-19 |
| N1 pure with noise delta = 1e-10, 1e-19, 1e-38 | 3.1e-11, 6.3e-10 | 1.3e-20, 3.7e-19 | 3.5e-40, 1.1e-38 |
| O8 alpha = 1e-8, noise in dd and qd | 2.2e-16, 2.4e-15 | 1.4e-20, 3.9e-19 | 3.5e-40, 1.1e-38 |

**D6 in f64 is NOT MET.** alpha = 1e-8, s = 2, dressed by poles at +-8i, K = 30: kind `Unresolved`, error of s 0.41, error of alpha 12, bar of s 1.1 |s|, so
criteria 1, 2 and 4 fail and criterion 3 holds. The dressing's own coefficients (1/64)^(k/2) are larger than the singular ones (1e-8 2^-k/k)
up to order about 16 and the contamination decays like 4^-k, so a clean window would start near order 24 of 30 (0.8 K, 7 rows, degree at most 3); the default
starts (0.2 to 0.6 K) are all contaminated. A window starting near 0.8 K would resolve this regime; it was not implemented, the coordinator having decided to leave the cell as it is.
D6 in double-double and quad-double meet all four criteria. The test `d6_f64_is_not_met` records what the locator does there; it does not pass the cell.

**The noise question** (choptuik-large-d #634: a singularity or noise in low orders). With coefficient noise delta in every order, the error of s was at most
0.14 delta and of alpha 6.6 delta for a single real singularity, 2.5 delta and 127 delta for a pair (dev_03b, other cases); the suite's noisy cells (N1, N2, O8) have errors well within their tolerances (no cell of the table has an error above 0.15 of its tolerance). With a dressing the error is the truncation of the dressing, not the noise (dev_03: 5e-21 in quad-double against delta = 1e-38, the same for three seeds), and cell N3 shows it (reported, not gated).
Not tested: a series with noise and **no** singularity; the entire function exp(u) of the exploration (`dev_00_probe.log`) gives a bar of 1 (`Unresolved`).

## 3. The controls

Each check has a control that must fail it; their targets were declared after their development run (`dev_11_controls.log`, amendment A3) and are asserted by `controls_fire`:
`wrong-sigma` and `wrong-alpha` (the nearest singularity planted off by 10 tolerances) fail the tolerance in 72 of 72 cells; `pair-by-ratio` fails the kind and the tolerance in 18 of 18 pair cells (the locator
says `Unresolved`, s off by 0.7 to 2.6); `bar-shrunk` (bars times 1e-3) fails coverage in 62 of 72 cells and **is silent in 10** (P1 in f64 and double-double, where the error is 0, and eight cells whose
error is at least a thousand times below the bar); the fingerprint control (one ulp of one coefficient) was silent as first coded (coefficient 10 of P1) and fires at coefficient 20 of D3 (A5).
**The fit window:** the ratio model with the first order 2 (degree 2) fails the tolerance on T1, T2 and T3 in all three lanes (a window too early: contamination by the farther singularity), and with
the last 8 orders (degree 2) on N1, N2 and T3 in all three lanes (a window too short: the noise). Both are silent on the pure cases, where the model is exact; that is
declared, not repaired.

## 4. Comparison with #183 (reported, not gated)

#183's raw folder keeps no Taylor coefficients. Each stored point (`x_last` and its state in double-double) was expanded again by #183's own `hp::taylor_p` in double-double to order 30 in
the step variable s = (x - x_last)/l, l = 2^round(log2 est.r), and the locator run on F, F', F'' (the harness is in the raw folder, `sys/compare183`; nothing in #183's crate is changed). 51 records:
the 20 `F->0` points of largest eta* (O8's regime, eta* from 8.7e3 to 2.0e4) and 11 of their first passes (R = 1e-4 instead of 1e-8), 10 `F->0` points at eta* from 1 to 100, and 6 + 3
`a3->0` records. In the output, per record: #183's eta*, bar, fitted alpha and alpha_pred = 2/(eta*^2 + 4) beside the locator's eta* with its bar and alpha(F) (= alpha of the field + the field number) with its bar.

| group | records | kinds | eta*: locator - #183 | alpha(F) against alpha_pred | #183's fitted alpha / alpha_pred |
|---|---|---|---|---|---|
| O8 regime, last step | 20 | 20 real | median 0.28, max 0.50 f64 ulp | 1.8e-16 relative at most, bar 8e-24 | 0.942 to 1.011 |
| O8 regime, first pass | 11 | 11 real | median 0.28, max 0.50 ulp | 1.8e-16, bar 5e-24 | 0.987 to 1.024 |
| `F->0`, eta* 1 to 100 | 10 | 10 real | median 0.14 ulp, max 183 ulp (eta* = 1.012, #183's own bar 1.4e-14: 2.8 bars) | 1.6e-14 relative (largest bar 2.5e-15 absolute) | 0.9992 to 1.0000 |
| `a3->0` (kappa curve) | 6 + 3 | 4 real, 2 pair (+ 3 real first passes) | median 37 ulp, up to 2e4 ulp (at most 17 bars of #183's, 12 bars of the locator's) | no alpha_pred | |

Read: (1) in the O8 regime the exponent is no longer noisy: the fit chooses its own window and degree, in double-double, and gives alpha(F) equal to 2/(eta*^2 + 4) to the f64 rounding of that number,
where #183's 10-order f64 fit was off by up to 6%. The difference of eta* from #183's stored value is below its f64 resolution (the stored eta* is an f64: 1.8e-12 at 1e4). The bar of the locator is 1e-34 there;
it is not a statement that eta* is that accurate: the state `s_last` is rounded to double-double (1e-32), so that is a floor. #183's O8 failure was in quad-double, **not reproduced here**: only double-double
states are stored, so a regeneration in quad-double would need #183's walk to be rerun. (2) Two `a3->0` records (eta* = 1.1 and 128) are **complex pairs** with Im eta of 3e-12 and 1.2e-13 and alpha(F) = 2.09 and 2.24 (bars 0.02, 0.09): the pair near the real axis that #183's A1 inferred from a change of the stop, here found directly; their real parts differ from #183's
eta* by 12 and 9 bars of the locator (17 and 7 of #183's). The others at the kappa curve have alpha(F) = 2.0 (a logarithmic end of F'', eta* = 4.899, bar 2.4e-7; and 1.945 +- 0.022 at 31.96) or 0.5 (eta* = 2.828 and 11.958, bar 1.8e-5 and
1.9e-2): not all `a3->0` points are logarithmic ends; #183's own fit at the first pass of those two points was 0.5006. This is for #183's review, not changed here.

## 5. What was not done

- The `+avx2,+fma` build was not run on the Pi (aarch64); the pull request's CI jobs (x86 SSE2, AVX2, portable, NEON on Linux and macOS, Windows MSVC) all pass, and the fingerprint test is asserted in each.
- A window starting near 0.8 K, the cure of D6 in f64.
- A series with noise and no singularity, and series in which c_(k-1) = 0 for half the orders (an even function: use u^2).
- A fit that takes the state's own error into account (the O8 comparison is limited by the double-double state).
- octo-double: the type is supported (`Real` for `OD`) and has no planted case.
- `singularity::domb_sykes` and `recurrence` stay; the locator supersedes them for series in the lanes' precisions, and they are not removed.

## 6. Reproduce

```
cd negative-mass-i191-libs/multisolver   # branch feat/locate-singularity, 0db3966
cargo test --release --test locate -- --nocapture > planted.log      # also: --features portable
cargo test --release --lib locate; cargo test --release --doc locate
# harness (negative-mass-raw/I191-rpi5/sys/compare183), then from notes/warp:
python3 singularity_locator.py RAW/I183-rpi5/dump_run_03a HARNESS_BIN planted.log > singularity_locator_output.txt
```

The criteria, the design exploration (dev_00 to dev_03b), the first run (dev_10), the controls' development run (dev_11), the amendments and the runs of record (run_01, run_02) are in
the raw folder's README.
