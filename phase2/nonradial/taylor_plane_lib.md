# taylor_plane_lib: the plane stepper in multishoot (#189)

Written by the runner of #189 (the coordinator reviews it). Raw output: `negative-mass-raw/I189-mac` (README with the criteria and
amendments A1 to A3, `dev_01` to `dev_14`, the runs of record `run_01` to `run_22`, `compare_switch.py`, `x4_digits.py`, `run_switch.sh`,
`run_lib.sh`). Status of every number below: *numerical evidence*.

**Result.** multishoot (pull request https://github.com/nmoinvaz/multishoot/pull/70, branch `feat/plane-stepper`, head a880903) has a module
`plane` that carries a k-plane of solutions of y′ = J(x) y through Taylor steps, orthonormalizes it after every step and returns the plane and
its log-volume. The planted checks P1, P2, P3 pass in f64, double-double and quad-double at the bars written before the first run, and #186's
crate runs on it behind the feature `lib_plane`: its 22 default checks pass, and every D and q token printed in the rows of #187 is the same with
the feature on and off. One criterion is met only in part (X3 below), and the issue's "skip one orthonormalization fails on a planted case" is
met on one cell, for a reason other than its intent.

## The facility

| Item | What it is |
|---|---|
| `plane::Real` | the scalar trait (`f64`, `DD`, `QD`, `OD`): `zero, one, from_f64, from_dd, to_f64, recip, sqrt, ln`, `EPS`. `Precision::Scalar` has no division, square root or logarithm and `f64` is not a `Precision`, so the facility has its own |
| `plane::JetSource` | `pattern()`: J's nonzero entries; `jets(x0, kj, out)`: coefficients 0..kj of each about x0. The caller owns the problem |
| `plane::Plane::new(n, k, basis)` | orthonormalizes the basis (Gram–Schmidt, positive diagonal) and records its log-volume |
| `Plane::advance(src, xa, xb, plan, kj)` | Taylor steps of order `plan.kt`, each the planner's rule (`taylor::step_from_norms`, `Plan`, `kth_root`) with the end-of-segment halving of `plan_path`; orthonormalization and `Σ ln R_ii` after each; `Stats { steps, hmin, hmax }` |
| `Plane::columns, column, ln_vol, ln_gram` | the orthonormal basis and the log of its volume as carried from the supplied basis |
| `plane::orthonormalize` | the Gram–Schmidt used |

Where multishoot's own types did not fit: `Work` and `Fl` hold fields of 8 B lanes on collocation points, `Jet<N>` is a series in the eight
lanes of a vector, and `Precision::Scalar` lacks the operations above; the plane's coefficients are scalars of one lane. The module shares `Plan`,
`kth_root` and the planner's step rule, which `taylor::step_size` now takes from `step_from_norms` (the same bits; its tests are unchanged).

## Planted results

Taylor tol 100 ε, order `Plan::optimal_kt`, hmax 1, x from 0. Bars (declared in the raw README before any run): plane `‖(I − Π) Y‖_F ≤ 1e4 ε`,
log-volume `≤ 1e4 ε max(1, |exact|)`, orthonormality `≤ 1e2 ε`. References are closed forms sharing no code with the stepper (the Gram terms ½ ln 5 and ln 2
are mpmath constants). Run of record: `run_17` (default backend, `mutate`), 20 tests in 7.8 s; the same digests under `--features portable` (`run_18`).

| case | system | exact | lane | plane error (bar) | ln_vol error (bar) | steps |
|---|---|---|---|---|---|---|
| P1 | J = P diag(24, 8, −8, −24) P⁻¹, cond P = 96, k = 2, basis P[:, :2], x = 12 | span P[:, :2]; ln_vol = 32 x + ½ ln 5 | f64 | 8.4e-16 (1.1e-12) | 4.5e-13 (4.3e-10) | 154 |
| | | | dd | 1.5e-31 (1.2e-28) | 6.3e-30 (4.7e-26) | 96 |
| | | | qd | 5.3e-64 (1.5e-60) | 1.8e-61 (5.8e-58) | 58 |
| | | | od (smoke) | 1.5e-129 (2.3e-124) | 2.0e-127 (8.9e-122) | 34 |
| P2 | J = Q′Qᵀ + Q Λ Qᵀ, Q = H Rot₁₂(θ₁) Rot₂₃(θ₂) a rotation, Λ = diag(20, 4, −8, −20), k = 2, basis Q(0)[:, :2] M, det M = 2, x = 12 | Φ = Q e^(Λx): span Q(x)[:, :2]; ln_vol = 24 x + ln 2 | f64 | 1.1e-16 (1.1e-12) | 5.7e-14 (3.2e-10) | 126 |
| | | | dd | 2.1e-32 (1.2e-28) | 0 (3.6e-26)¹ | 79 |
| | | | qd | 6.6e-65 (1.5e-60) | 6.7e-63 (4.3e-58) | 48 |
| P3 | J = P diag(…) P⁻¹ + x A₁ + x² A₂, k = n = 4, basis U (det 1), x = 4 | Liouville: ln_vol = ∫ tr J = −56/3 | f64 | n/a | 1.4e-14 (2.1e-11) | 116 |
| | | | dd | n/a | 1.3e-29 (2.3e-27) | 74 |
| | | | qd | n/a | 6.5e-61 (2.8e-59) | 47 |

¹ The computed dd log-volume equals the reference's parts exactly: the accumulator's resolution at 288.7 (3.5e-30) is coarser than the terms' noise.

Bitwise across backends: the digest of every part of the plane and of `ln_vol` of all ten cells above is the same on the default build (aarch64 NEON) and
with `--features portable`. The `+avx2,+fma` build was not run (no x86 machine here).

## Controls

Wrong inputs (no feature needed), each failing its check in f64, dd and qd: `J_gamma` (γ₂ scaled by 1 + 1e-6; ln_vol error 9.6e-5), `J_noQdot` (P2 without
Q′Qᵀ; plane 4.3e-3, ln_vol 8.0e-2), `J_trace` (tr A₁ off by 1; 8.0), `wrong_plane` (P1 against span P[:, 1:3]; plane 0.76), `ref_no_gram` (the exact log-volume
without the starting basis' term; 0.80 and 0.69). Library controls (feature `mutate`, `control::PlaneControl`), targets declared after development runs
(raw README A1):

| control | fails | silent by design |
|---|---|---|
| `LowSum(8)` (Taylor sum cut at order 8) | P1, P2, P3 × f64, dd, qd (ln_vol 5e-3 to 16) | |
| `ScaleOnly` (columns rescaled between steps, orthonormalized once at the end) | P1, P2, P3 × f64, dd, qd (ln_vol 40 to 510); P3 in f64 the stepper refuses (dependent columns) | |
| `LogDrop` (last column's ln R left out) | P1, P2, P3 × all lanes (48 to 129) | |
| `NoGram` (the starting volume left out) | P1, P2 × all lanes (0.80, 0.69) | P3 (starting volume 1) |
| `RecipDD` (1/k formed in double-double) | P1, P2, P3 in qd (3e-31 to 7e-30 against bars of 1e-58) | f64, dd (identical bits) |

`EndOnly` (no rescaling at all) was replaced by `ScaleOnly` after `dev_09`: the unnormalized columns grew like e^(24x), the planner's absolute tolerance
shrank the steps with them, and the f64 run of P1 took 10 114 309 steps.

## Skipping one orthonormalization

The issue's control "skip one orthonormalization fails on a planted case". In exact arithmetic it changes neither the plane nor the log-volume, and the
coordinator expected it silent unless one step's growth spread exceeds `ln(1/ε)`. `SkipQr(i, 1)` after the first, the middle and the second to last step:

- silent on P1 and P2 in f64, dd, qd, and on P3 in f64 and dd (dd within 7% of the bar at the first step: 2.2e-27 against 2.3e-27; baseline 1.3e-29);
- **on P3 in qd it exceeds the bar**, 7.8e-58, 1.5e-58, 2.6e-58 against 2.8e-59 (baseline 6.5e-61). That is not the growth spread: on one planned step (largest h) (γ₁ − γ₂) h is
  1.2 to 4.0 on P1 and P2, (γ₁ − γ₃) h 2.5 to 7.0, and (γ₁ − γ₄) h 3.0 to 6.7 on P3, against ln(1/ε) = 36.7, 73.5, 147. The error does not move with the Taylor tolerance
  (1000 ε down to ε: qd 4.3e-58, 7.8e-58, 1.6e-58, 4.6e-58; dd 6.6e-28, 2.2e-27, 1.1e-27, 1.7e-27; `dev_04`), so it is roundoff, scaled by the condition number of the
  unnormalized 4-column product (cond P = 96), 30 to 1000 times the baseline's.

The control declared in its place, by the coordinator's ruling, is `ScaleOnly` (fixed, the same in every lane). Skipping a run of m consecutive steps
(`dev_05`, `dev_07`, kept as development logs, not declared: choosing m per lane after seeing the errors is tuning) loses accuracy smoothly with the accumulated
spread, not at `ln(1/ε)`. The test `skip_one_{f64,dd,qd}` in multishoot prints these outcomes and asserts nothing.

## The crate's switch (feature `lib_plane`)

`phase2/nonradial/taylor_polar/`: `Cargo.toml` has the feature and a `[patch.crates-io] multishoot = { git, rev = a880903… }` beside the multidouble patch (the
dependency `multishoot = "=4.2.0"` is satisfied: the branch is 4.2.0 plus the module). `num.rs`: `Num` gets `type Lib` and `to_lib`, `from_lib` (identities for the
four lanes) under the feature. `plane.rs`: under the feature, `gram_schmidt` calls `plane::orthonormalize`, and `integrate` builds a `Plane`, steps it through a
bridge from the crate's `JetSource` (its jets at s² copied into the library's entries; the `TP_MAX_SECONDS` guard is kept there) and copies the plane back. With
the feature off the crate's own code is unchanged. `taylor_polar.py` is not edited and builds without features, so the binaries were built by cargo and copied
over `target/release/taylor_polar` before each run (`run_switch.sh`):

```sh
cd phase2/nonradial/taylor_polar
cargo build --release -j 2                                                # the crate's own integrator
cargo build --release -j 2 --features lib_plane --target-dir target/lib_plane   # multishoot's plane stepper
cp target/lib_plane/release/taylor_polar target/release/taylor_polar      # then, from phase2/nonradial: python taylor_polar.py ...
```

| Criterion (README of the raw folder) | Verdict | Number |
|---|---|---|
| X1 builds on and off; no absolute path in `Cargo.toml` | met | `rg /Users Cargo.toml` is empty |
| X2 default mode, on and off | met | 22 passed, 0 failed, 1 skipped (TS), exit 0, on (957 s) and off (855 s) |
| X3 every token of `--rows` equal on and off | **not met as worded** | rows and header lines 746 of 746 tokens equal (l = 2, 3, 4, 5, 6, five variants, capped rows); RA2, RA4, RA5, RA6 pass on both sides (66, 87, 87, 87 of the same number of stored tokens). **14 tokens of the two capped-row CHECK lines differ** (CAP5b 6, CAP6b 8): the load average, per-variant seconds, and the printed error bars E2 (order 24 against order 30: 1.32e-24 and 2.44e-26 off, 1.42e-24 and 8.81e-25 on) and E3 (f64 against dd: 1.53e-09 and 6.90e-09 off, 3.58e-09 and 9.60e-09 on) |
| X4 full-digit D | reported | l = 2, γ = 1: f64 2.7e-15 (14.6 digits agree), dd 2.9e-31 (30.5), qd identical; γ = 40: f64 2.2e-11, dd 1.4e-32, qd 5.6e-45; l = 4, γ = 300 dd 4.0e-28. Over all 780 dd values (D, q_int, q_ext, five variants): 333 bit-identical, worst relative difference 1.8e-24 (base), 6.5e-17 (the Taylor-tol-1e-12 variant) |
| X5 the comparison can fail | met | the on-binary with `--mutate=moncrief` against the off run: 39 mismatches, exit 1 (`run_16`) |

The 52 differing tokens of the default mode's CHECK lines (not part of X3) are errors of the planted checks against mpmath at roundoff scale, all within their bars,
the L1 check's f64-lane D (11 printed digits: −1.0204252371e-09 off, −1.0204252373e-09 on; −1.5957912910e-13 and −1.5957912954e-13), V3's error bars and timings.
The step counts of the integrations are the same on and off; their points differ by the planner's rule (`kth_root` against `powf`; `plan_path`'s end-of-segment halving).
Timing, dd, base variant, one thread, each side under its own load: l = 6, γ = 1000 6.413 s off, 6.423 s on; γ = 12 000 76.25 s and 78.29 s (loads 2.4 and 3.3); f64 γ = 12 000
2.086 s and 2.353 s; l = 5, γ = 6000 38.29 s and 38.30 s. Peak resident 57 MiB per rows run, 108 MiB (off) and 93 MiB (on) the default mode.

## For the coordinator

- X3 as I worded it ("every token … the same") is not met: 14 tokens of the CAP5b and CAP6b lines differ. They are timings, the load average and two printed error
  bars that are themselves differences between integrations at the truncation and f64-roundoff scale; every D and q token of every row is equal. I left the
  criterion as written and report it so.
- #186's `plane.rs` comment says "classical Gram–Schmidt with one reorthogonalization"; its loop updates the column between projections, so it is modified
  Gram–Schmidt applied twice. The library does what the code does and says so. I did not edit the comment (yours to correct on main).
- With the feature on, the crate's mutations `lowsum` and `noqr` (in its own `integrate`) are not applied (a message says so); multishoot's `LowSum` and `ScaleOnly` are
  their counterparts. The recip mutation (`recip_dd`, in `recips`) still applies where `recips` is used (the Riccati and the star), not in the stepper.
- The step rule differs from #186's (`powf`, no sliver rule) by an ulp or more in the step points; D differs at the truncation and roundoff scale (above), not in its printed digits.

## Reproducing

From a clone of multishoot at `feat/plane-stepper`: `cargo test --release --features mutate --test plane -- --nocapture` (the planted checks and the controls; the lines
`P1 dd: …` carry the numbers and the digest) and `cargo test --release --features portable,mutate --test plane -- --nocapture`. The crate's comparison:
`python3 compare_switch.py OFF.log ON.log` from the raw folder, with the logs of two `python taylor_polar.py --rows …` runs.

## Coordinator's review

Under review: multishoot pull request #70 (`feat/plane-stepper` at a880903, **not merged**: its review and release are
the user's), `wip/i189-plane-lib` at 6c40e02, raw `I189-mac` at bb964b6. Review folder
`negative-mass-raw/I189-mac-review`, whose README holds the criteria, committed before the first run of `ind189.py`.
*Update, 2026-11-27:* the user squash-merged multishoot#70 into multishoot main as 70fea18. No release contains it yet
(the latest is 4.2.0). The merged `plane.rs` adds one change: a column whose length falls outside 2^±256 is rescaled by
a power of two before normalizing, and other columns are normalized bitwise as before. The crate still pins a880903,
which is on no branch now. #199 repins it to main and checks the rows bit for bit.
The changed files are those the brief names.

- **ind189, a planted case the runner did not use.**
  - The system is a rotated Airy pair: u″ = x u and v″ = 4x v, coupled by a product of Givens rotations with rational
    cosines. Its growth is super-exponential, (2/3)x^{3/2} and (4/3)x^{3/2}, against the constant and rotating rates of
    P1–P3.
  - The exact plane and volume come from mpmath's Ai and Bi. The pull request's `Plane::advance` is called from a small
    harness (`airy_plane/`).
  - Run 3 meets Q1 (k = 1), Q2 (k = 2) and Q3 (k = 4, Liouville) in f64, dd and qd at X = 4 and 10, each at least 1e3
    below its bar (1e4 ε). At X = 10 the plane errors are 1.5e-16, 1.3e-32 and 3.9e-65.
  - The mutations fail exactly their targets: `wrong-scale` fails Q1 and Q2; `scaleonly` (multishoot's `ScaleOnly`)
    fails Q2 and Q3 and leaves Q1, which has nothing to orthogonalize.
  - Run 1 crashed on a formatting bug in ind189. Run 2 failed Q3 at X = 10 in dd and qd: the coordinator's reference
    carried 90 digits against a cancellation of about 73 in det(YᵀY). Amendment A1 raised it to 200 digits, with no bar
    changed. Both runs are kept.
- **Rerun of the switch.** The crate was built with `--features lib_plane` in a review worktree, and the rows for l = 2,
  4 and 6 and the default mode were run. Against the coordinator's feature-off reruns of #187, every row and header line
  is identical, with the same exit codes and verdicts. The CHECK lines differ only where the runner said they would: E2
  and E3 (CAP6b's E3 9.60e-9 on against 6.90e-9 off), the f64 lane's last digits in L1, planted errors at roundoff, and
  timings.
  - *Added by the coordinator:* the first l = 2 rerun skipped RX. The review worktree sits one level deeper, so
    `taylor_polar.py`'s default `--raw` (`../../../negative-mass-raw` from its folder) did not exist. It was rerun with
    `--raw` given. That makes two harnesses (#184's and this one) caught by a path default. A skipped check now reports
    as skipped (#187's fix), which is how it was seen.
- **X3 as worded is not met** (14 tokens of the CAP5b and CAP6b lines). The issue's "skip one orthonormalization fails"
  is met on one cell only, for a reason other than its intent; the coordinator wrote that line, and ruled `ScaleOnly` in
  its place. Both stand as the runner reports them.
- *Added by the coordinator:* the "classical Gram–Schmidt" comment in `taylor_polar/src/plane.rs` is corrected on main.
  The loop recomputes each projection from the updated column, so it is modified Gram–Schmidt applied twice.
- Follow-ups:
  - #199: `taylor_polar.py --build --features`, and the feature printed in each run;
  - #188: the complex plane stepper, added to its acceptance;
  - after multishoot#70 is reviewed and released by the user, a task to drop the crate's multishoot git pin.
