# A Taylor recursion for implicit second-order ODEs in multishoot, and #167's crate on it (issue #190)

Status words as in AGENTS.md; every number here is **numerical evidence**. Numbers are in
[implicit_taylor_lib_output.txt](implicit_taylor_lib_output.txt), which `sys/make_output.sh` of the raw folder writes from the logs
(negative-mass-raw `I190-rpi5`, README with the criteria, the development runs, the amendments A1 and A2 and the erratum E1).
Machine: rpi5-ubuntu (Raspberry Pi 5, aarch64), 2 threads, nice 10, load about 1 from other processes at each start.

## Answer

The pull request is [multishoot#73](https://github.com/nmoinvaz/multishoot/pull/73): a module `implicit` that forms the Taylor
series of N(x, y, y', y'') = 0, a polynomial in y'', order by order through one division by dN/dy'', carries the root branch, and
reports where dN/dy'' vanishes. Every planted check that is gated is met in f64, double- and quad-double, and every control fails
exactly its declared checks; three results are **not met** and recorded (the relative error of each Taylor coefficient, the
`Conserve` mode through a double root, one f64 fold estimate). With the cargo feature `implicit_lib` #167's crate runs its
recursion and the walk of #184 and #183 on the module, and #167's, #184's and #183's checks pass with the same eta\*, classes and
kinds to the printed digits.

## The facility

The caller gives the monomials (a, b, c) of N = sum_g A_g(x) y^a y'^b y''^c and, at each step's start, the Taylor coefficients of
every A_g (a jet); the module knows nothing of any one problem.

| piece | what it does |
|---|---|
| `Recursion` | order k of N set to zero and solved for y''_k, which enters linearly with the factor A0 = dN/dy'' at the order-0 values: one division per step; y_k, y'_k from y'_(k-1), y''_(k-1); order 0 of N not imposed, so N is constant along the flow. The step variable is s = (x - x0)/l, l a power of two. Cauchy sums in multidouble's dot accumulators |
| `Walker`, `WalkOpts`, `WalkOut` | steps by `taylor::step_from_norms` on the largest coefficient of orders kt - 1 and kt (each field relative to its \|c_0\| or not), checkpoints, a landing rule, Horner's sum; `form`, `planned`, `advance` for a caller that watches each step |
| `Branch::Project`, `Conserve` | the branch is the data's and is carried, never picked by its order among the roots; `Project` makes one Newton correction of y'' at each start, `Conserve` none |
| `Stop::ZeroPivot`, `Pivot`, `WalkOut::zero`, `crossings` | where A0 vanishes: no division by zero; a stop with an estimate of the zero; the double roots a branch passes, located by the zero of A0's series in the step |
| `Scalar`, `Real` | f64, double-, quad- and octo-double, and `Fl<T, 1>` of any lane (duals, for #167's matching) |

The estimate of the zero of A0, x\* = x - beta A0/A0' with beta = 1/(1 - A0 A0''/A0'^2), is exact for A0 ~ c \|x\* - x\|^beta in
exact arithmetic and an estimate otherwise; the planted folds are pure powers (exponents 1/2 and 1), and its error there is the
conditioning of A0' near the fold.

## Planted truth

Three polynomial equations with closed forms, nine branches that each meet other roots at double roots (x = 1, sqrt 3, a root of
u tanh u = 1, y = 0): N1 = y''^2 - (y + x y')^2 (quadratic), N3 = N1 (y'' - 2y) (three real roots), N4 = N1 (y''^2 - 4y^2) (four);
two folds, 16 y''^2 + 4 y'^2 + 225 x + 25 x^3 = 0 (exponent 1/2) and 3 y' + 2 y''^3 = 0 (exponent 1), ending at x = 0; a start at an
exact double root. Lanes f64, dd, qd at order `Plan::optimal_kt(tol)` = 22, 42, 80 (tol 1e-14, 1e-30, 1e-62), bar 1e6 eps
(1.1e-10, 1.2e-26, 1.5e-58), chosen before the first run; an error bar here is that bar on the max-norm relative error at three
landings, not a fit.

- **Met, gated** (per lane): `closed_form` 15 cells (worst 3.5e-14, 4.1e-30, 3.4e-62), `crossings` 12 (worst 7.8e-14, 4.5e-30,
  2.8e-62), `fold` 3, 4, 4 (the estimate within 1e-3 of the distance; worst 4.1e-5, 7.7e-19, 6.7e-51), `step_sum` 11 (below),
  `zero_pivot`, `planner` (bitwise), `lanes` (`Fl<DV,1>` and `Fl<Du,1>` equal to dd bit for bit); an octo-double smoke (6.4e-127
  against 2.3e-122).
- **`step_sum`**, declared in A1 after the first development run: at each planned step h, sum_k \|dc_k\| h^k <= 1e6 eps
  sum_k \|c_k\| h^k, with the recursion started from the closed form's state at each step's start. Worst ratio 2.9e4 f64, 7.1e3 dd,
  1.8e4 qd (bar 1e6).
- **NOT MET (1), `coeffs`**: the relative error of every Taylor coefficient against 1e6 eps: 4.9e-8 f64 (k <= 22), 2.0e-11 dd
  (k <= 42), 2.9e-15 qd (k <= 80). It grows about threefold per order (f64: 1.5e-16 at k = 3, 4.9e-8 at k = 22): order k of N is
  s^2_k minus a convolution of the y'' coefficients, and for coefficients that fall off super-exponentially the convolution
  terms exceed c_k, so c_k has an absolute error of about eps \|c_0\| and a relative one that grows. What a step uses is within
  roundoff (`step_sum`). The criterion stays as worded (coordinator).
- **NOT MET (2), `Conserve` through a double root**: K-B, R-B and R-D in every lane (9 cells) and the crossings of K-B and R-B
  (6): N is conserved at its start value, a rounding eps and not 0, and for N = eps a double root is two roots; the crossing
  becomes an avoided crossing, the Taylor series about a point at distance d has radius about d, the steps shrink toward it and the
  walk ends `Tiny` (K-B: 33, 33, 26 steps; errors 7.4e-8 f64, 1.6e-18 dd, 3.2e-50 qd). `Project` crosses. Six branches are
  met in `Conserve` (Q-A, Q-B, K-A, K-C, R-A, R-C), although Q-B, K-A, K-C, R-A and R-C cross a double root too: in these runs the
  crossing survived the rounding of N, which I did not explain further. Stated in the module's docs: near a double root use
  `Project`.
- **NOT MET (3), f64 F1 `Conserve`**: the fold estimate is 5.7e-15 from the fold against a bar of 6.3e-16: N drifts by about eps
  (sum of \|terms\| = 250) per step over 84 steps, and a drift dN moves the fold of N = const by dN/(dN/dx) = dN/225.
- Erratum E1 to the first development table: the estimate is exact for a pure power in exact arithmetic, not to eps; its error
  here is the conditioning near the fold (raw README).

## Controls

Six library controls (`SkipOrder(3)`, `RecipDD`, `NoPivotCheck`, `NoCrossing`, `ExponentHalf`, `NoPivotStop`) and the three that
need no hook (the branch choice: K-A started on K-B's root, K-A switched to K-C's root at x = 0.6, K-A against K-B's closed form;
`closed_form` fails by 1.5, 0.94, 1e4). Targets were declared in A2 after a development run with `step_sum` and are asserted by
`controls_f64`, `controls_dd`, `controls_qd`: `SkipOrder(3)` fails `step_sum`, `closed_form`, `crossings` and both folds;
`RecipDD` only in quad-double (`step_sum`, `closed_form`, `crossings`), a no-op in f64 and dd; `NoPivotCheck` only `zero_pivot`;
`NoCrossing` only `crossings`; `ExponentHalf` only the exponent-1 fold; `NoPivotStop` both folds. A prediction of the README that
was too narrow: `SkipOrder(3)` fails more than `coeffs` and `closed_form`. In the first development table, cells that were
themselves NOT MET were counted as controls firing; A2 uses gated cells only.

## #167's crate on the facility

`notes/warp/conical_global/`, Cargo feature `implicit_lib` (default off, enables multishoot's `mutate`): `poly_lib.rs` runs
`poly::implicit_taylor` (the matching with duals) on `Recursion<Fl<T, 1>>`, `hp_lib.rs` runs `hp::hp_walk` on `Walker` with the
locator, the checkpoints, the trace and the statistics of N kept in the crate; `hp.rs` and `poly.rs` differ from #183's only by
`cfg` attributes, three `pub(crate)` and a `pub use`. multishoot is taken from the PR branch by a git `rev` (eaf55a8), no local
path. The baseline is the 9debda2 tree built against multishoot 4.2.0 from crates.io; "off" is the same code against the PR's rev.

| script (arm) | exit | checks | controls | wall s | max rss MB |
|---|---|---|---|---|---|
| conical_global.py base / on | 0 / 0 | 17 / 17 | all on target / all on target | 1814 / 1856 | 181 / 181 |
| conical_loglog.py base / on | 0 / 0 | 20 / 20 | all on target / all on target | 1534 / 1201 | 98 / 98 |
| conical_offgrid.py `--grid=ctl --no-controls` base (with the series pickle) / off / on / on at eaf55a8 | 0 | 16 each | - | 680 / 317 / 162 / 162 | 99 / 72 / 72 / 72 |
| conical_offgrid.py `--grid=ctl --controls-only` n-drift, recip-dd base / on | 0 / 0 | - | both fired on target only | 1370 / 626 | 72 / 72 |

The comparison (`sys/compare_outputs.py`, rules fixed from the baseline's lines: error sizes below 1e-6, timings and step counts are
masked, everything else must be equal as a string): 126 lines of conical_global.py, 124 of conical_loglog.py and 68 of the `ctl`
run have no gated difference; in the 685 verdict records of the `ctl` run, 537 eta\* values are equal to the 10 printed digits and
the largest relative difference of any is 2.5e-15, the difference of two valid step sequences (10 records have a different step
count), not an error bar; the feature-off build equals the baseline exactly (0.0, same step counts). The feature-on runs are
not slower: the `ctl` run takes 162 s with the feature against 317 s without (same machine, 2 threads, load about 1).

**What did not go as planned.** (1) The first feature-on loglog run failed `plant_scale` (19 of 20): the double-double `conserve`
walk that the check predicts to fail ended `TinyStep` at 1430 after 397 steps instead of `End` at 3000 after 132 (all six
`--hscale` values of the baseline end `End`, one of six of the library). The cause is the Cauchy sums: `hp.rs` forms them in
multidouble's dot accumulator, the library summed plain products, and the planted regime (N's scale falls by 7e-38) cancels. The
check was not changed; `Scalar::cauchy` (the accumulators) went into the library and `W<P>::cauchy` into the crate, because the
library change alone did not reach the crate's own scalar (`run_23`, the same failure to the digit). After it the walk equals the
baseline's bit for bit. Both failing runs are kept (`run_10`, `run_23`). (2) The comparator was written after the first
feature-on offgrid run had finished, not before (the queue ran ahead), from the baseline's lines; one rule (the series pickle's
build message) was added after the first comparison. (3) The library suite and the rebase: the planted values moved by up to a
factor 16 in 21 of 147 cells when the accumulators went in (none near a bar), and the digests changed (dd 92ec461781147d77, qd
75902b0bf76ced2a, f64 b124d44e9758f1a1 unchanged); after the rebase onto multishoot main (`step_from_norms` of #70 is the same
function as the copy, which is dropped; conflicts in `control.rs` and the README, both kept) they are unchanged on the default and
portable backends, and the feature-on arms were not rerun apart from the `ctl` check once.

## Not claimed

The octo-double lane is a smoke, not a planted suite. The zero estimate is exact for a power only. Relative accuracy of each
coefficient is not given (growth about 3 per order). `Conserve` is not valid through a double root. Backends: default (aarch64
NEON) and `--features portable` were run, bitwise equal digests; `+avx2,+fma` was not run (no x86 machine). `Scalar` and
`plane::Real` overlap and should be unified.

## Reproduce

```sh
cd ~/Source/negative-mass-i190-libs/multishoot   # feat/implicit-taylor-r2 (PR #73)
cargo test --release --features mutate --lib --test implicit --test taylor --test taylor_order --test plane -- --nocapture
cargo test --release --features portable,mutate --test implicit -- --nocapture
cd ~/Source/negative-mass-i190/notes/warp/conical_global && cargo build --release --features implicit_lib
# then copy the binaries to target/release and run, from the repository root, with --nobuild:
python3 notes/warp/conical_global.py --nobuild
python3 notes/warp/conical_loglog.py --nobuild
python3 notes/warp/conical_offgrid.py --grid=ctl --no-controls --nobuild --dump=DIR
```

## Coordinator's review

Under review: multishoot pull request #73 (`feat/implicit-taylor-r2` at eaf55a8, on main 70fea18, not merged; its review *Update 2026-11-28:* merged (2827f0a) and released in multishoot 5.0.0; #216 moves the crate to the release.
and release are the user's), `wip/i190-implicit-taylor` at 550800c, raw `I190-rpi5` at 5f2ac0f. Review folder
`negative-mass-raw/I190-mac-review` (at ddf8e6f), whose README holds the criteria, committed before the first run of
`ind190`. The changed files are those the brief names. The baseline files `hp.rs` and `poly.rs` only gain `cfg` gates
and wider visibility.

- **Rerun on the Mac**, also aarch64 NEON.
  - The library suite at eaf55a8, default and `portable`: everything passes, with the Pi's digests (dd
    92ec461781147d77, f64 b124d44e9758f1a1, qd 75902b0bf76ced2a).
  - The crate at 550800c with `implicit_lib`: offgrid `ctl` passes 16 of 16, with all 537 η* equal to the Pi's baseline
    in the 10 printed digits. loglog passes 20 of 20 with every control on target, and nothing differs in its 124
    compared lines. The review's first loglog run crashed on a path of the review's own setup, and is kept.
- **ind190**, run 1, 3 of 3. It uses three implicit equations of the coordinator's choosing, none of them in the PR's
  tests, through `Walker` in f64, dd and qd:
  - E1: a fold of exponent ½ at x = 1;
  - E2: two roots that cross at x = π/2 on the cos branch;
  - E3: a cubic in y″ with three real roots, walked on the middle and the top root.

  Every closed form holds at the end point to 2.6e-16, 1.3e-32 and 4.0e-65. The crossing is located at π/2 to
  roundoff. The fold walk stops `Tiny` within 5e-13 of x = 1. The mutations `wrong-root`, `cross-at` and `fold-at`
  each fail exactly their target.
- **Not met, as ruled:** `coeffs` as worded; `Conserve` on K-B, R-B and R-D; the f64 F1 `Conserve` fold. The `plant_scale`
  failure was fixed in the library and the crate, and the check is unchanged.
- *Added by the coordinator:* **the PR's CI fails**, in `tests/slaved.rs`
  (`slaved_residual_returns_the_factorization_error`), on every platform. multishoot main at 70fea18 fails the same
  test, so the defect is on main and does not come from this PR; the runner's report did not mention the red CI. The
  CI does not assert the implicit digests, so the AVX2 build's bitwise equality is still open (#214).
- *Added by the coordinator:* `Cargo.toml` takes multishoot from the PR commit eaf55a8 for every build, the feature
  off included. The runner checked that the feature-off output equals the baseline exactly. If the PR is
  squash-merged, the commit leaves every branch, as #70's a880903 did; #216 moves the crate to a release.
- **Follow-ups:** #214 (AVX2 digests), #215 (one scalar trait and `Branch::Auto`, after #73 merges) and #216
  (`implicit_lib` as the only path, after a release).
