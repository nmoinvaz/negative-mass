# The singularity locator: later window starts, and series with noise, an entire function or a weak singularity

November 2026, issue [#206](https://github.com/nmoinvaz/negative-mass/issues/206), a follow-up to
[#191](https://github.com/nmoinvaz/negative-mass/issues/191) ([`singularity_locator.md`](singularity_locator.md)). Pull request
[nmoinvaz/multisolver#28](https://github.com/nmoinvaz/multisolver/pull/28) (branch `feat/locate-late-windows`, commits c5b5ea7 and 21b9829 on multisolver main 63d6e4d). Raw folder: negative-mass-raw `I206-mac` at 035fc85
(README with the readings of the acceptance, the criteria and the control inputs written before the first run of record, amendment A1 with the
development runs and the controls' targets, every log). Status words as in AGENTS.md: **numerical evidence**; nothing is proved or certified. Every number quoted
below is in a log of the raw folder named in section 7.

## Answer

The first orders 0.7 K and 0.8 K are in the default grid. **D6 in f64 is still not met**, but it is no longer far: the locator now chooses the ratio model at 0.8 K, locates s to 2.1e-4 (it was 0.41)
and alpha to 1.0e-2, and its bars cover and are useful, but the bar of s is 5.3% of |s| and the resolution limit is 5%, so the kind is `Unresolved`. The other 68 gated cells still meet #191's criteria 1 to 4.
**Noise-only (300) and entire (9) series are all `Unresolved`. Weak singularities at 1, 10 and 100 times the noise are not never confidently wrong: one of 2700 is**, a cell that is as wrong on the grid of #191, so the change did not cause it.
Of the eight criteria, six are met and two (1, the D6 cell, and 6, the weak singularities) are not.

| criterion (README section 4) | verdict | number |
|---|---|---|
| 1. D6 in f64 meets #191's criteria 1 to 4 | **NOT MET** | `Unresolved` (ratio, first order 24, degree 1); error of s 2.1e-4 (tolerance 1e-2), of alpha 1.0e-2 (0.1); bars cover, bar of alpha 2.7 (at most 10); bar of s 0.053 \|s\| against the cutoff 0.05. Only criterion 1 (kind) fails; on #191's grid 1, 2 and 4 failed (error of s 0.41, bar 1.1 \|s\|) |
| 2. the other 68 gated cells still meet 1 to 4 | met | 68 of 68 (`planted_suite`); the six controls of #191 fire as declared, `bar-shrunk` silent in eight cells (ten before) |
| 3. exploration coverage at factor 2: covered at least 254 and at most 2 uncovered | met | Python instrument: 259 resolved, 257 covered, 2 uncovered (the two of #191, whose bar is exactly 0); control (bars shrunk a thousand times) covers 44. The shipped locator on the same series: 256 of 256 resolved cells covered (254 of 254 on #191's grid), control 31 |
| 4. noise-only series `Unresolved` | met | 300 of 300 (100 seeds in each lane); smallest bar of s over \|s\| 0.38 (f64) |
| 5. entire functions `Unresolved` | met | 9 of 9 (exp(u), exp(-3u), sin(sqrt v) / sqrt v; bars of s 1.5 to 2.8 \|s\|) |
| 6. weak singularities never confidently wrong | **NOT MET** | 1 of 2700: W3 (alpha 1/2, s = 1) in f64 at 100 times the noise, seed 14: `Real` by the ratio model at first order 6, degree 1, s = 0.874 (error 0.126) with a bar of 0.027; the same cell on #191's grid. The other 2699 are `Unresolved` |
| 7. the fingerprint is the same on every backend | met | 0x0f1c7236bd56966f on aarch64 NEON and with `--features portable` here; the PR's test jobs at 21b9829 all pass the same assertion: Ubuntu x86_64 SSE2, AVX2 and portable, Ubuntu ARM64 NEON, macOS ARM64 NEON, Windows MSVC |
| 8. fmt, clippy (both builds), docs, the PR's CI | met | here: `run_08_*` exit 0; CI at 21b9829: all 15 checks pass (CodeRabbit's own check reads "skipping") (Rustfmt, Clippy SSE2, NEON and portable, Docs, MSRV, the six test jobs, Coverage, codecov); CodeRabbit reviewed 21b9829 (one full review requested, after its rate limit): no actionable comments, no thread |

## 1. What changed

`Opts::default().starts` is 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 (it was to 0.6). Nothing else of the locator changes: the degrees, `bar_factor` 2, `max_rel_bar` 0.05, the floors and the rule that
chooses the cell are #191's. At K = 30 the window at 0.8 K starts at order 24 and has 7 ratio rows (degree at most 3) and no three-term row; at 0.7 K, 10 ratio rows and 9 three-term rows (degree 2). The rest of the pull request is tests
in `tests/locate.rs`: the new series, their controls, and a fingerprint that now covers them.

## 2. D6 in f64

D6 is alpha = 1e-8 dressed by poles at +-8i, K = 30. #191 found every window from 0.2 to 0.6 K contaminated by the dressing, and a clean one near order 24. With the late windows:

| grid of first orders | chosen cell | error of s | bar of s over \|s\| | kind | criteria failed |
|---|---|---|---|---|---|
| 0.2 ... 0.6 K (#191) | three-term, first order 6, degree 2 | 0.41 | 1.1 | `Unresolved` | 1, 2, 4 |
| 0.2 ... 0.7 K (without 0.8) | ratio, 21, degree 1 | 2.9e-3 | 0.22 | `Unresolved` | 1, 4 |
| **0.2 ... 0.8 K (the default now)** | ratio, 24, degree 1 | 2.1e-4 | 0.053 | `Unresolved` | **1** |
| 0.7 and 0.8 K only | ratio, 24, degree 1 | 2.1e-4 | 0.053 | `Unresolved` | 1 |
| 0.8 K only (a data point; not the default and not adopted) | ratio, 24, degree 2 | 3.9e-4 | 0.021 | `Real` | none |

So 0.8 K is what brings the error from 2.9e-3 to 2.1e-4 and the bar from 0.22 to 0.053; the bar is 5.3% against the 5% cutoff. The cell chosen on the default grid is the one whose eight neighbors agree best. **Not done:** raising `max_rel_bar`, dropping the earlier windows, or changing the choice rule, to pass the cell; the coordinator's instruction (as for #191) is to leave a cell not met and report it.
The test `d6_f64_is_not_met` records the failures `["1 kind"]`, the kind and the bar above 5%, and that the location is within the tolerances.

## 3. The suite, and what the late windows did to the other cells

`planted_suite`: 68 of the 69 gated cells meet all four criteria, with the default grid as with #191's. 28 of the 72 rows (24 cases in three lanes, N3 included) now choose a late window, 7 at 0.7 K and 21 at 0.8 K, where a late window
agrees better with its neighbors: D1 to D6 and N3 in every lane, C4 in every lane, C5 in f64 and quad-double, T1 and T2 in quad-double, that is the dressed cases and the farther singularities in the long lane. The pure, pair, two-pole and noisy-pure cases (P1 to P6, C1 to C3, T3, N1, N2, O8) keep their early windows
(37 rows at 0.2 K); the errors stay within the tolerances. The control `bar-shrunk` (bars times 1e-3) is silent in **eight** cells now, not ten: P1 in f64 and double-double,
D1 and D2 in double-double and quad-double, C4 and T1 in quad-double; it fires on C4 in double-double and T2 in quad-double, whose chosen cells changed. The control `grid-191` (D6 in f64 on #191's grid) fails exactly criteria 1, 2 and 4 as recorded there;
`grid-no-0.8` fails it again (criteria 1 and 4).

## 4. Exploration coverage (criterion 3)

The design exploration of #191 (`sys/explore` of `I191-rpi5`: 90 random cases of seeds 1 and 2, three lanes, fit degrees up to 9, the eight-neighbor spread, Python and mpmath) was rerun on this Mac with first orders 0.2 ... 0.8 K. The instrument is
checked first: on #191's grid it reproduces the 540 rows of `I191-rpi5/dev_02_rows_d9_*` exactly (every field; 256 resolved, 254 covered). It has no floor on its bar and needs one row fewer per ratio fit than the port, so the same 270 series were also dumped and run
through `multisolver::locate` as shipped.

| | grid of #191 | default grid |
|---|---|---|
| Python instrument: resolved (bar of s below 5%), covered at factor 2, uncovered | 256, 254, 2 | 259, 257, 2 |
| Python instrument: control, bars shrunk a thousand times, covered | 52 | 44 |
| shipped locator: resolved, covered, kind right | 254, 254, 254 | 256, 256, 256 |
| shipped locator: control covered | 38 | 31 |

The two uncovered cells of the instrument have a bar of exactly 0 (a rounding-level error that all neighbors reproduce); the floor of the shipped bar covers them.

## 5. Noise-only and entire series (criteria 4 and 5)

Noise-only: c_k = delta xi_k, k = 0 ... K, xi_k uniform of unit variance, delta = 1e-10, 1e-19, 1e-38, seeds 1 to 100 in each of f64 (K = 30), double-double (60) and quad-double (120): 300 series. All 300 are `Unresolved` under both grids; The smallest bar of s over |s| among them is 0.38 in f64, 0.47 in double-double and 0.62 in quad-double on the default grid (0.70 in f64 on #191's grid): 7 to 12 times the 5% cutoff, so a cutoff anywhere between D6's 0.053 and 0.38 would leave all 300
`Unresolved`. (What a cutoff in that range does to the weak series was not measured.) The entire functions exp(u), exp(-3u) and sin(sqrt v) / sqrt v (series in v) in the three lanes: 9 of 9 `Unresolved`, bars of s 1.57 to 2.22 |s| on the default grid
(the ratios of these series fall like 1/k to 0: there is no constant term to read s from).
Controls: a pole fed to the same check is `Real` in the three lanes; a locator that never says `Unresolved` (`max_rel_bar` 1e300) fails the check on 300 of 300 noise-only and 9 of 9 entire series.

## 6. Weak singularities (criterion 6)

Three real singularities at radius about 1 (the variable scaled so that the radius is of order one): W1, alpha = -1 at s = 1.25; W2, alpha = -1/2 at s = -0.8; W3, alpha = 1/2 at s = 1; each with additive noise, scaled so that the root mean square of
its coefficients c_1 ... c_K is m times the noise's, m = 1, 10, 100, 100 seeds in each lane: 2700 series. Confidently wrong means a `Pair`, or a `Real` whose error of s is above its bar.

| result at m = 1, 10, 100 (2700 series) | default grid | grid of #191 |
|---|---|---|
| `Unresolved` | 2699 | 2698 |
| `Real`, within its bar | 0 | 1 (W2, f64, m = 100) |
| **confidently wrong** | **1** (W3, f64, m = 100, seed 14) | **1** (the same) |
| no location | 0 | 0 |

At 1, 10 and 100 times the noise the locator says it cannot resolve the singularity, and that is the honest answer. The one resolved result is wrong: the ratio model at first order 6, degree 1 reports s = 0.874 (truth 1) with a bar of 0.027, and alpha 0.82 off with a bar of 0.80.

The amplitude at which a weak singularity is resolved (reported, `weak_sweep`, 100 seeds per cell; resolved right + confidently wrong, default grid):

| | m = 100 | 300 | 1e3 | 1e4 | 1e6 | 1e9 |
|---|---|---|---|---|---|---|
| W1 f64 | 0 + 0 | 1 + 1 | 10 + 2 | 99 + 1 | 99 + 1 | 99 + 1 |
| W1 double-double | 0 | 0 | 0 | 0 | 45 + 0 | 100 + 0 |
| W2 f64 | 0 | 7 + 0 | 70 + 0 | 100 + 0 | 100 + 0 | 100 + 0 |
| W2 double-double | 0 | 0 | 0 | 95 + 0 | 100 + 0 | 100 + 0 |
| W2 quad-double | 0 | 0 | 0 | 1 + 0 | 100 + 0 | 100 + 0 |
| W3 f64 | 0 + 1 | 9 + 1 | 55 + 1 | 100 + 0 | 100 + 0 | 100 + 0 |
| W3 double-double | 0 | 0 + 4 | 37 + 2 | 100 + 0 | 99 + 1 | 99 + 1 |
| W3 quad-double | 0 | 0 | 16 + 12 | 99 + 1 | 100 + 0 | 100 + 0 |

(W1 in quad-double is `Unresolved` throughout.) The first resolved results appear at m = 100 to 300 in f64 (W3 at 100, one, wrong; W1 and W2 at 300) and, in double- and quad-double, from m = 300 to beyond 1e9 by singularity (W3 at 300 in double-double and 1e3 in quad-double, W2 at 1e4, W1 at 1e6 in double-double and not by 1e9 in quad-double). **Where it is first resolved the locator is often wrong**: W3 in quad-double at m = 1e3, 12 wrong of 28
resolved; in double-double at m = 300, 4 of 4. Over the 9000 series of the sweep the default grid has 2340 right and 30 confidently wrong, the grid of #191 2202 and 30, the wrong ones in the same cells: the late windows add right results and no wrong ones.
The bar is an estimator range over fit choices, not a bound (#191); this is the size of the miss near the noise floor. #207 (the locator in the walkers) should treat a result whose bar is of the order of the noise's effect as `Unresolved`; no such rule is in the locator.

The control for the check (`wrong-location`: a singularity planted 10% off) is confidently wrong in 9 of 9 cells; `bar-shrunk-weak` fires on the noise-free series in 8 of 9 (silent on W3 in f64, where the error is 0); on the noisy series it is **vacuous**, since 2699 of 2700 are `Unresolved`.

## 7. Proposed text for `singularity_locator.md`

Replace, in "Answer" (lines 15 to 18): "**68 of the 69 gated cells meet all four criteria; one does not: D6 in f64.** There the locator says `Unresolved` with a bar of 1.1 |s| ... A window starting near 0.8 K would resolve that regime; it is not implemented." with:

> **68 of the 69 gated cells meet all four criteria; one does not: D6 in f64.** With the first orders 0.7 K and 0.8 K added (#206) the locator locates it (error of s 2.1e-4, of alpha 1.0e-2) but its bar of s, 5.3% of |s|, is above the resolution limit of 5%, so it says `Unresolved`: it gives no confident wrong answer. A series of pure noise (300) or an entire function (9) is `Unresolved` too. A weak singularity in noise (1, 10, 100 times the noise's size) is `Unresolved` but for one result of 2700, which is wrong with a bar that does not cover it: the bar is an estimator range, not a bound, and near the amplitude where a weak singularity is first resolved it is not safe to read it as one.

In the method table: "first order 0.2 K to 0.6 K" becomes "first order 0.2 K to 0.8 K". In section 1: "first order of the window 0.2, 0.3, 0.4, 0.5, 0.6 K" becomes "0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 K".
In section 2, replace the paragraph "**D6 in f64 is NOT MET.**" from "A window starting near 0.8 K would resolve this regime" on with: "With the windows at 0.7 K and 0.8 K (#206) the chosen cell is the ratio model at 0.8 K, degree 1: error of s 2.1e-4, bar 5.3% of |s|, `Unresolved` (criterion 1 fails; 2, 3 and 4 hold). Without 0.8 K it is 2.9e-3 and 22%. The test `d6_f64_is_not_met` records it."
In section 2, "The noise question", replace "Not tested: a series with noise and **no** singularity ..." with: "Tested (#206): 300 series of additive noise only and the entire functions exp(u), exp(-3u) and sin(sqrt v) / sqrt v are all `Unresolved`; real singularities whose coefficients are 1, 10 and 100 times the noise are `Unresolved` in 2699 of 2700 series, and the other is confidently wrong (W3 in f64, m = 100). The amplitude at which they are first resolved is 100 to 300 times the noise in f64 and from 300 to beyond 1e9 in double- and quad-double, by singularity; there the locator is wrong in up to 43% of its resolved results (section 6 of `singularity_locator_206.md`)."
In section 3, "`bar-shrunk` ... **is silent in 10**" becomes "**is silent in 8**" (P1 in f64 and double-double, where the error is 0, and six cells whose error is at least a thousand times below the bar: D1 and D2 in double-double and quad-double, C4 and T1 in quad-double), and the fingerprint is 0x0f1c7236bd56966f (see section 8 on the constant in section 1).
In section 5, delete "A window starting near 0.8 K, the cure of D6 in f64." and "A series with noise and no singularity" (keep "series in which c_(k-1) = 0"), and add: "A rule for a weak singularity near the noise floor, where the bar is not safe (#206, section 6)."

## 8. Where this disagrees with the repository's text

`singularity_locator.md` (line 52 and section 3) and the description of multisolver#26 quote the fingerprint 0xa5113e9ccb6776c8. The merged test (`tests/locate.rs` at eb52618, which squashed #26's review rounds) asserts 0xf3a13cba157d2e3c, and passes on this Mac
at main (63d6e4d): the note's constant is that of an earlier state of the suite. It was not changed here; the constants of this pull request are 0xc534cbdcc58395a4 (c5b5ea7, the grid alone) and 0x0f1c7236bd56966f (21b9829, with the series of #206 and the kind and the model folded in).

## 9. What was not done

- The `+avx2,+fma` build was not run on this Mac (aarch64); the PR's AVX2 job is its check, and the fingerprint test is asserted in every test job.
- Any change to make D6 in f64 or the weak cell pass: `max_rel_bar`, the grid and the choice rule are #191's.
- A rule against the confident-wrong cells near the noise floor (for example a floor on the bar from the noise estimated from the highest orders); recommended below.
- Octo-double has no planted case, as in #191.

## 10. Reproduce

```
cd negative-mass-i206-libs/multisolver   # branch feat/locate-late-windows, 21b9829
cargo test --release --test locate -- --nocapture      # also: --features portable
cargo test --release --lib locate:: ; cargo test --release --doc locate
cargo test --release --test locate late_windows_report -- --ignored --nocapture --test-threads 1
cargo test --release --test locate d6_grid_report -- --ignored --nocapture
cargo test --release --test locate weak_sweep -- --ignored --nocapture --test-threads 1
# the exploration: negative-mass-raw/I206-mac/sys/run_explore.sh <tag> 0.2,0.3,0.4,0.5,0.6,0.7,0.8 (python 3.14, mpmath 1.3.0), then
# sys/explore_rust on the dumped series
```

Logs of record in `I206-mac`: `run_01` and `run_01b` (tests, default), `run_02` and `run_02b` (portable), `run_03` (the exploration), `run_04` (the port on its series), `run_05` (the report under both grids),
`run_06` (D6 under five grids), `run_07` (the sweep), `run_08` (fmt, clippy in both builds, docs). Development runs `dev_01` to `dev_11` and amendment A1 are in the README.
