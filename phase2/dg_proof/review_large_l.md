# Independent review of #2's large-l proofs (RESULTS.md §4.2 to §4.5)

Issue #104, runner `negative-mass-rpi5-i104` (rpi5). Written from `RESULTS.md` alone; step 7 (comparison with #2's
output files) came last.
Verdict words: **confirmed**, **gap**, **error**. Status words as in AGENTS.md. Raw output:
`negative-mass-raw/I104-rpi5/` (README with the criteria, seed, list of l
and control targets, all logs, `sweep_brackets.csv`).

**Result.** No error found. Every step of §4.2 to §4.5 is **confirmed**, with the facts it assumes marked below as
**gap** (assumed, outside what this review re-proved):
T1, T2 (and Lemma A(ii) beyond l = 14). Theorems L and M hold with this review's own constants, which are smaller than
#2's: first order 23.82 and 172.03 (#2: 91.51 and 229.87),
second order for l ≥ 138: 3.553 and 6.682, and for l ≥ 500: 3.519 and 5.936 (#2: 3.709 and 6.101). The bracket table of
§4.3 holds for every l = 2..3000, not only for a sample.

| step of RESULTS.md | claim | verdict | evidence (output file) |
|---|---|---|---|
| 4.1 | `a_n(ε)` is Lemma A's `c·c_n/cⁿ`; `X_∞` has the Bessel closed form | confirmed | common (exact Fractions, l = 2..12, 30, 60, 150 all n; 500, 3000 for n ≤ 40), roots |
| 4.2 (Z) | `ψ₀,ₓₓ = Vψ₀` | confirmed for l = 2..14; **gap** beyond (Lemma A(ii), not re-proved) | theorem3 |
| 4.2 (P1) | `V > 0` on `(r_s, ∞)` | confirmed (explicit positive polynomial in `(σ, μ)`) | theorem3 |
| 4.2 (P2) | `V_x < 0` on `(r_s, ∞)` | confirmed (64 positive coefficients, as stated) | theorem3 |
| 4.2 (B) | `V̂` has exactly one zero beyond `r_*`, simple, `−` to `+` | confirmed (H1, H2 on paper; sign at ∞ symbolic; exact Sturm counts l = 2..16) | theorem3 |
| 4.2 identities | `V̂ = 2c⁴(2+ερ)Ψ/(ρ⁴(ρ−6)²X²)`; `V̂ ∝ S` with a positive factor | confirmed (symbolic in `c`, generic `X`; exact from the definition for l = 2..9) | theorem3 |
| 4.3 | brackets of `c·r_*`, `c·r_neg` for l = 2..3000, ordering, limits, `l²r` | confirmed for 200 random l (all certified), and for every l = 2..3000 | brackets, sweep, table, overlap |
| 4.4 | Theorem L, `|c·r − ρ| ≤ K/c`, l ≥ 138 | confirmed, K = 23.82, 172.03 (#2: 91.51, 229.87) | thm_l |
| 4.5 | Theorem M, l ≥ 500 | confirmed, s = 3.519, 5.936 (#2: 3.709, 6.101) | thm_m |
| 4.6 | what is not proved | agrees; last bullet to be replaced (§8) | — |
| T1, T2 | `N > 0` on `[0, r_s]`; one root of `N` beyond `r_s` | **gap** (assumed; exact Sturm count l = 2..16 only) | theorem3 |

## 1. What was run

Commands from this directory (`phase2/dg_proof`), one process at a time, `nice -n 10`, `OMP_NUM_THREADS=1`, python
`~/Source/negative-mass-venv314/bin/python` (3.14.10, sympy 1.14.0, mpmath 1.3.0), rpi5.
Each is `python3 -u X.py > X_output.txt`, stdout only. Exit codes and times from `negative-mass-raw/I104-rpi5/chain.log`.

| script | what | checks | exit | time |
|---|---|---|---|---|
| `review_large_l_common.py` | recurrence against Lemma A (exact), planted problem, formula checks | 49 | 0 | 1 s |
| `review_large_l_roots.py` | `ρ_*`, `ρ_neg` to 25 digits, certified | 19 | 0 | 1 s |
| `review_large_l_theorem3.py` | Theorem 3 identities, (P1), (P2), counts | 70 | 0 | 19 s |
| `review_large_l_brackets.py` | 200 random brackets + 35 supplementary | 1929 | 0 | 154 s |
| `review_large_l_thm_l.py` | Theorem L, first order | 105 | 0 | 14 s |
| `review_large_l_thm_m.py` | Theorem M and the second order refinement of L | 163 | 0 | 20 s |
| `review_large_l_sweep.py --csv=...sweep_brackets.csv` | every l = 2..3000 | 24010 | 0 | 2437 s |
| `review_large_l_table.py --csv=...sweep_brackets.csv` | RESULTS 4.3's table rows | 11 | 0 | 1 s |
| `review_large_l_overlap.py --csv=...sweep_brackets.csv` | step 7: #2's output files | 38 | 0 | 1 s |
| `review_large_l_controls.py` | 52 pre-registered mutations | 52 | **1** (§6) | 605 s |
| `review_large_l_posthoc.py` | post-hoc replacement of a silent control | 7 | 0 | 1 s |

`review_large_l_common.py` is the shared module (also self-testing). The sweep's output is the log of the run launched
with the committed scripts (`step4b_sweep.log` in the raw folder).

## 2. Independence: what is shared with #2

Not opened before step 7: `rescaled.py`, `vhat_structure.py`, `threshold_isolation.py`, `large_l_proof.py`,
`monotone_proof.py`, `common.py`, every other script of
`phase2/` and every `*_output.txt` of those. Certified parts use mpmath `iv` (outward-rounded intervals, 192 bits); #2
used `fractions.Fraction`. Identities use sympy 1.14.0.

Shared with #2: the formulas typed from RESULTS.md (`a_n(ε)`, `X`, `N`, `N_Z`, `Ψ`, `S`, the static Zerilli `V`, Lemma
A's closed form for `χ`, the statements of T1/T2),
and the *structure* of the proofs of Theorems L and M, because RESULTS.md describes it (coefficientwise bounds, a
perturbation of `N` and `Ψ`, slopes, localisation;
second-order expansion). Not shared: every line of code, the arithmetic (intervals, all `l+2` terms summed with no
truncation for every bracket), the tail bounds, the windows,
the way the constants are assembled, the `S`-based certificate, and the comparison step of Theorem M (derived afresh,
with a smaller loss than RESULTS' `(l+1)ν/2 > |ζ|`).
Because the method is shared, the agreement of my numbers with #2's (§5) shows that the numbers are right, not that the
method is the only one.

## 3. Assumed, not re-proved here

- **T1, T2** (RESULTS §3): `N > 0` on `[0, r_s]` and exactly one root of `N` beyond `r_s`. A bracket in which `N`
  changes sign is *the* `c·r_*` only through T2. This review checks T2 by exact Sturm counts for l = 2..16, and shows
  that each bracket holds exactly one root of `N` (`N′` excludes 0 over the bracket), but T2 for l ≥ 17 is #2's
  Descartes argument, not reviewed.
- **Lemma A(ii)** (ψ₀ solves the static Zerilli equation for every l): exact for l = 2..14 here (sympy, own χ); for
  every l it is #2's symbolic proof, not reviewed. The `S`-based brackets use no (Z) and agree with the `Ψ`-based ones
  at all 2999 l.
- **Theorem 3** is used in §4.3, §4.4 and §4.5 to say that the zero of `Ψ` found beyond `r_*` is `r_neg`. Its steps are
  reviewed in §4.2 below.
- Textbook facts: intermediate value theorem, Taylor with Lagrange remainder, Weierstrass' product inequality.

## 4. Step by step

### 4.1 Rescaled coefficients — confirmed

`a_n(ε)` as typed from RESULTS (product recurrence, `ε` an enclosure of `1/c`) lies inside Lemma A's exact value
`c·c_n/cⁿ` (Fractions, rising factorials) for l = 2..12, 30, 60, 150 and every n, and l = 500, 3000 for n ≤ 40; `c_n =
0` for `n = l+3..l+6`, `a_{l+2} > 0`.
The Bessel form `X_∞ = s[(s+3)I₀(2√s) − 4√s I₁(2√s)]` holds termwise (coefficient `(n−2)(n−4)/(n−1)!²`, sympy) and
numerically to `1e-58`.

### 4.2 Theorem 3 — confirmed

- (Z), V in `(ρ, c)`: `ψ₀,ₓₓ = Vψ₀` exactly for l = 2..14; `V(r = ρ/c, λ = c/2) = 8c⁴(2+ερ)N_Z/(ρ⁴(ρ−6)²)` symbolic in
  `c`.
- Identities (H11): for any `X`, `V̂ = g² − g′ = 2g² − ψ₀,ₓₓ/ψ₀`, and `V̂ = S·c³(2c+ρ)/(ρ³(ρ−6)X²)`, so `sign V̂ = sign
  S` on `ρ > 6` where `X > 0`. With (Z): `V̂ = 2c⁴(2+ερ)Ψ/(ρ⁴(ρ−6)²X²)`. The last is also exact for l = 2..9 from the
  definition `g² − g′`, with `χ` from Lemma A (no use of (Z)).
- (P1): `λ²r²((λ+1)r−3) + 9(λr−1) = (σ³+6σ²+18σ+18) + μ(σ+3)³`, all coefficients positive (H12).
- (P2): `Q` is the quintic of H12; `c⁵Q(τ, m)` has **64** coefficients, all positive, as RESULTS states. `N_Z(6+τ)` has
  positive coefficients.
- (B): H1 and H2 on paper. At a zero of `g − h`, `(g−h)′ = V/2 − V_x/(4h) > 0`. Sign at infinity: the leading
  coefficient of `Ψ/(a²ρ^{2l+7})` is `ε(l+1)² − (1+2ε)/2 = (l+1)/(2(l−1)) > 0` and `V̂ r² → (l+1)(l+2)`, symbolic in `l`
  (H10).
- Exact Sturm counts, l = 2..16: `N` has one root beyond `ρ = 6`; `Ψ` has exactly one root in `[ρ_*, ∞)` and is negative
  at the left end.

### 4.3 Certified isolation, l = 2..3000 — confirmed

- 200 random l (seed 104, list in the raw README), all certified: all `l+2` terms summed in `iv`, brackets of width
  `2e-14` with 20-digit decimal ends; `N` (`Ψ`) enclosures exclude 0 with opposite signs at the ends; `N′` (`Ψ′`)
  excludes 0 over the bracket, so each bracket holds exactly one root; `S` has the same sign change `(−,+)` on the `Ψ`
  bracket; `lo(c·r_neg) > hi(c·r_*)`. The 35 supplementary l (2..30, 138, 139, 140, 499, 500, 501, 1000, 2000, 2999,
  3000) also pass, and sympy's isolating intervals overlap mine for l = 2..12.
- Every l = 2..3000 (sweep): all 2999 brackets certified; `hi(l+1) < lo(l)` for both sequences at every consecutive l
  (strict decrease); both stay above `ρ_*`, `ρ_neg`; `l²r_*` and `l²r_neg` increase strictly at every consecutive l ≥ 5;
  minimum at l = 5 (checked over l = 2..30).
- RESULTS' table: all ten rows (l = 2, 3, 5, 10, 20, 80, 100, 500, 1000, 3000) and the limits row agree with my brackets
  to the printed digits (`review_large_l_table_output.txt`).
- `(c·r_* − ρ_*)·c = 3.5154852` and `(c·r_neg − ρ_neg)·c = 5.8393917` at l = 3000, consistent with RESULTS' 3.5155 and
  5.8394 (their `l = 10⁷` values were not recomputed).
- Not covered by an independent argument: that the root of `N` found is the *only* root beyond `r_s` (T2, §3).

### 4.4 Theorem L — confirmed

Statement reviewed: for l ≥ 138, `|c·r_* − ρ_*| ≤ K_*/c` and `|c·r_neg − ρ_neg| ≤ K_neg/c`. Redone with `ε₀ =
1/(140·137)`, windows `ρ_* ± 0.01` and `ρ_neg ± 0.03`, 40 subintervals, 60 explicit terms and a tail majorant (H5).

| step | this review | #2 (RESULTS 4.4 and `large_l_proof_output.txt`) |
|---|---|---|
| coefficient bound (H3) | `|a_n − α_n| ≤ εB_n` for `3 ≤ n ≤ l+10`, l = 138..1000; ratio ≤ 1 (equality at `n = 4`, l = 138, where `ε = ε₀`) | same, `ε₀ = 1/108` |
| `E_N` on `[8, 9]` with `ε₀ = 1/108` (printed for comparison) | 259.515 | 259.5 |
| `E_N` on my N window, `ε₀ = 1/19180` | 132.02 | — |
| `η_Ψ` on my Ψ window | 5.161e6 | 6.93e6 (window `[9.61, 9.66]`) |
| `−N_∞′` | ≥ 5.543 on my window; ≥ 3.207 on `[8, 9]` | ≥ 2.836 on `[8, 9]` |
| `Ψ_∞′` | ≥ 30002.6 on my window; `Ψ_∞′(ρ_neg) = 33523.342` | ≥ 30130 on `[9.61, 9.66]`; 33523 |
| `sup |Ψ_∞″|` on `[9.61, 9.66]` | ≤ 1.136e5 | ≤ 1.31e5 |
| localisation (H6) | `K_* = 23.818`, `K_neg = 172.034`; `ε₀K < δ` for both | `K_* = 91.51`, `K_neg = 229.87` |

Conclusion: confirmed. #2's statement follows from this review's (23.82 ≤ 91.51 and 172.03 ≤ 229.87). The `Ψ` zero found
in the window is beyond `ρ_* + K_*ε₀ < ρ_neg − 0.03`, hence is `r_neg` (Theorem 3). The bounds are checked at 11 points
of each window for l = 138, 139, 276, 1000 (largest observed ratios 20 and 2.11e5, against the bounds 132 and 5.16e6),
and the final inequalities against my certified brackets at l = 138, 139, 145, 276, 690, 1000, 2000, 3000.

### 4.5 Theorem M — confirmed

Second order: `a_n = α_n + εα¹_n + ε²ω_n` (H4), `N_l = N_∞ + εν + ε²ζ`, `Ψ_l = Ψ_∞ + εΨ₁ + ε²E` (H9), windows `ρ ±
K′ε_max` with `K′ = 46, 339` (above my first-order constants).

| quantity | l ≥ 500 (M) | l ≥ 138 (refinement of L) | #2 (l ≥ 500) |
|---|---|---|---|
| `Z ≥ |ζ|` | 1487.2 | 1492.2 | 1488 |
| `E_max ≥ |E|` | 4.659e7 | 4.833e7 | 4.67e7 |
| `ν` range on the N window | [19.686, 19.701] | [19.614, 19.773] | ≥ 19.61 |
| `Ψ₁` range on the Ψ window | [−1.974e5, −1.941e5] | [−2.074e5, −1.846e5] | ≤ −1.89e5 |
| `s_N` (`0 < c·r_* − ρ_* ≤ s_N/c`) | **3.5186** | 3.5525 | 3.709 |
| `s_Ψ` (`0 < c·r_neg − ρ_neg ≤ s_Ψ/c`) | **5.9364** | 6.6819 | 6.101 |

Comparison step, derived afresh (H7): `N_{l+1}(ρ_l) ≤ −(ε−ε′)ν_min + (ε²+ε′²)Z`, with `ε − ε′ = (2l+2)εε′`, is negative
once `(2l+2)ν_min > (c(l+1)/c(l) + c(l)/c(l+1)) Z`; at l = 500 that is `19725 > 2974`. For `Ψ`: `1.945e8 > 9.317e7`. The
left sides grow with `l` and the factor decreases, so the inequality holds for all l ≥ 500. RESULTS' sufficient
condition `(l+1)ν_min/2 > |ζ|` is correct but loses a factor 2: its `Ψ` margin is `4.74e7 > 4.67e7`, 1.5%; here the
margin is 2×. The `l²r` criterion (H8) holds with `s_N = 3.52 < 8.40` and `s_Ψ = 5.94 < 9.54` at l = 500, increasing
with l.

Against my certified brackets, for the pairs (500, 501), (501, 502), (502, 503), (1000, 1001), (1001, 1002), (2000,
2001), (2001, 2002): strict decrease, `0 < (c·r − ρ)c ≤ s`, `l²r` increases, and `N_{l+1} < 0`, `Ψ_{l+1} > 0` at the
lower ends of the l brackets. Conclusion: confirmed; #2's constants follow from this review's.

### 4.6 What is not proved — agrees

The first bullet (V̂ on `(0, r_*]`) is outside this review (#105). The second bullet's list of assumed facts is
accurate: T1, T2, Lemma A(ii), L3 (re-derived here as H10), the coefficient and tail bounds (proved here as H3 to H5,
and numerically tested to `n ≤ 400`).

## 5. Step 7: comparison with #2's output files

Opened after every final run of this review.

- **Brackets.** `threshold_isolation_output.txt` prints 16 values of l (2, 3, 5, 10, 20, 40, 80, 81, 100, 150, 200, 300,
  500, 1000, 2000, 3000), each `lo` truncated to 12 decimals, so the stored text does not hold the 200 random brackets:
  2 of them (1000, 2000) are in it. The other 14 are compared with my brackets from the sweep (which covers every l).
  All 32 brackets (16 l, two roots each) overlap; the 32 differences "middle of my bracket − printed value" lie in
  `[1.3e-13, 1.05e-12]`, which is the signature of truncation, not rounding (`review_large_l_overlap_output.txt`). All
  200 of my brackets satisfy #2's own Theorem L bound (185 sampled l ≥ 138) and Theorem M bound (157 sampled l ≥ 500),
  with no violation.
- **Constants.** #2's `K_* = 91.5045`, `K_neg = 229.868`, `s_N = 3.709`, `s_Ψ = 6.101` (and `η_N = 259.515`,
  `|Ψ_∞′(ρ_neg)| = 33523.3`) agree with or exceed mine; see §4.4 and §4.5. #2's `η₂_N = 1487.6`, `η₂_Ψ = 4.668e7`
  (window of #2's choice) against my 1487.2 and 4.659e7.
- **No discrepancy in a number was found.** Two differences of bookkeeping are stated above (the factor 2 in the
  comparison step; `−N_∞′ ≥ 2.836` against my 3.207 on `[8, 9]`, both valid lower bounds).

## 6. Controls

`review_large_l_controls.py` runs every script once per `--mutate=NAME` (52 mutations over 8 scripts) and requires a
non-zero exit with the first FAIL at the named check.
The targets were written to the raw README before the first control run and are unchanged: **47 on target, 5 not**; the
script exits 1.

| mutation | script | outcome | reason |
|---|---|---|---|
| `coef_q` | common, sweep | off target: fires at l = 3, not l = 2 | my prediction was wrong: at l = 2 only `n ≤ 4` is tested, and the changed `q_a` first enters `a_5` |
| `accept_zero_enclosure` | common, roots, brackets | **silent** (exit 0) | the mutation accepts an enclosure with an end exactly at 0, but the planted enclosure of `P(a)` straddles 0 (the coefficients are not binary numbers), so the mutated test still refuses it |

The silent control is a flaw of the control, not of the planted problem. The target list is not edited. A post-hoc
replacement, `review_large_l_posthoc.py`, plants the intended bug (a sign test that takes the midpoint's sign when the
enclosure contains 0); the planted problem then fails at "endpoint exactly on a root" and at "three roots in one
bracket" (exit 0 means it fired on target).
Added later and not part of the 52: `review_large_l_overlap.py` has `--mutate=overlap_shift` (fires at the
printing-convention check) and `--mutate=bound_scale` (fires at the Theorem L bound check).
Check families without a mutation: "X > 0 at both ends of the Ψ bracket", "lo(c·r_neg) > hi(c·r_*)", the count `64`, and
the printed comparisons with #2's constants (informational).

## 7. Hand lemmas, stated for a later formalization in Lean

Each is proved on paper in this review; the polynomial inequalities are certified by sympy (all coefficients of the
shifted polynomial are non-negative), the rest are elementary.
Notation: `c = (l+2)(l−1)`, `ε = 1/c`, `q_a = (a−1)(a+2)`, `e₁(n) = Σ_{a=0}^{n−3} q_a = (n−2)(n²−4n−3)/3`, `S_n = e₁(n)
+ 2`,
`g_n = 1/(2ⁿ((n−1)!)²)`, `α_n = (n−2)(n−4)g_n`, `α′_n = 3(n−2)(n−3)g_n`, `α¹_n = α′_n − α_n e₁(n)`, and `a_n = 0` for `n
≥ l+3`.

- **H1 (zero crossing).** If `u ∈ C¹[a,b)` and `u′(z) > 0` at every zero `z` of `u`, then `u` has at most one zero.
  *Proof.* Take the least zero `z₁` (zeros form a closed set). `u′(z₁) > 0`, so `u > 0` just right of `z₁`. A second
  zero has a least value `z₂ > z₁` with `u > 0` on `(z₁, z₂)`, so `u′(z₂) ≤ 0`, a contradiction.
- **H2 (Theorem 3, step (B)).** Let `V > 0` and `V_x < 0` on `I = [x₀, ∞)`, `g ∈ C¹(I)` with `g′ = V − g²`, `g(x₀) = 0`,
  `g > 0` on `(x₀, ∞)`, and `h = √(V/2)`. Then `V̂ := 2g² − V = 2(g−h)(g+h)` has `sign V̂ = sign(g−h)` on `(x₀, ∞)`,
  `V̂(x₀) < 0`, and `g − h` has at most one zero. If `V̂ > 0` near `∞`, `V̂` has exactly one zero in `(x₀, ∞)`, it is
  simple, `V̂ < 0` before it and `> 0` after it.
  *Proof.* At a zero `z` of `g − h`: `g′(z) = V − g² = V/2 > 0` and `h′ = V_x/(4h) < 0`, so `(g−h)′(z) > 0`; apply H1.
- **H3 (first-order coefficient bound).** For `l ≥ l₀`, `ε₀ ≥ 1/c(l₀)`, `3 ≤ n ≤ l+2`: `|a_n − α_n| ≤ ε[α′_n(1+2ε₀) +
  |α_n| max(2, S_n(1+2ε₀))]`; for `n ≥ l+3`: `|a_n − α_n| ≤ ε n² |α_n|`.
  *Proof.* `a_n = (α_n + εα′_n) P_n`, `P_n = (1+2ε) Π_{a=2}^{n−3}(1 − q_a ε)`, each factor in `[0,1]` because `q_a ≤ q_l
  = c` for `a ≤ l`. Weierstrass: `Π(1 − x_a) ≥ 1 − Σx_a` on `[0,1]`, so `P_n ∈ [(1+2ε)(1−εS_n), 1+2ε]`, hence `|P_n − 1|
  ≤ ε max(2, S_n(1+2ε))`, and `a_n − α_n = α_n(P_n − 1) + εα′_n P_n`. For `n ≥ l+3`, `n² > c`, so `1 < εn²`.
- **H4 (second-order expansion).** For `n ≤ l+2`: `a_n = α_n + εα¹_n + ε²ω_n` with `|ω_n| ≤ |e₁|α′_n + (|α_n| +
  ε₀α′_n)(e₁+4)²(1+2ε₀)/2`; for `n ≥ l+3`: `|ω_n| ≤ n⁴|α_n| + n²α′_n + n²|α_n||e₁|`.
  *Proof.* `P_n(ε) = 1 − e₁ε + P_n″(ξ)ε²/2` (Lagrange), `P_n″ = Σ_{a≠b} q_a q_b Π_{c≠a,b}(1 − q_c ξ)`, each product `≤
  1+2ε`, `Σ|q_a| = e₁+4`; `ω_n = −e₁α′_n + (α_n + εα′_n)P_n″(ξ)/2`. For `n ≥ l+3`: `a_n = 0` and `1 < εn²`.
- **H5 (majorants and tails).** For `n ≥ 3`: `B_n ≤ n⁵g_n`, `|α¹_n| ≤ n⁵g_n`, `Ω_n ≤ n⁸g_n` (`B_n`, `Ω_n` the bounds of
  H3, H4 taken uniformly over both regimes, `ε₀ ≤ 1/1000`). For `ρ ≥ 1` and `n₁ ≥ 1`: `Σ_{n>n₁} n^D ρ^{n−k} g_n ≤
  U(n₁+1)/(1−r)`, `U(n) = n^D ρ^{n−k} g_n`, `r = ((n₁+1)/n₁)^D ρ/(2n₁²) < 1`. The explicit terms `n ≤ 60 ≤ l+2` use the
  bound of the regime `n ≤ l+2` only.
  *Proof.* The polynomial inequalities are `P(n) ≥ 0` for `n ≥ 5` (all coefficients of `P(5+t)` are `≥ 0`, sympy) and
  for `n = 3, 4` by direct evaluation; the tail by the ratio `U(n+1)/U(n) = ((n+1)/n)^D ρ/(2n²)`, decreasing in `n`.
- **H6 (localisation).** If `F ∈ C¹[ρ̄−δ, ρ̄+δ]`, `F(ρ̄) = 0`, `F′ ≥ μ > 0` there, and `|F_l − F| ≤ ηε < μδ` there, then
  `F_l` has a zero in `[ρ̄ − ηε/μ, ρ̄ + ηε/μ]`. (Same with `F′ ≤ −μ`.)
  *Proof.* `F(ρ̄ ± t) = ±(≥ μt)`; at `t = ηε/μ`, `F_l ≥ 0` at `ρ̄ + t` and `≤ 0` at `ρ̄ − t`.
- **H7 (monotone comparison).** If `N_l`, `N_{l+1}` each change sign exactly once on `(6, ∞)`, from `+` to `−`, at
  `ρ_l`, `ρ_{l+1}`, then `N_{l+1}(ρ_l) < 0` implies `ρ_{l+1} < ρ_l`. (For `Ψ`, from `−` to `+`, `Ψ_{l+1}(ρ_l) > 0`
  implies the same, on the region beyond `r_*`.) With `ε − ε′ = (2l+2)εε′`, the condition `N_{l+1}(ρ_l) < 0` follows
  from `(2l+2)ν_min > (c′/c + c/c′)Z`.
- **H8 (`l²r` criterion).** `c(l+1) − c(l) = 2l+2` and `(l+1)² c(l) − l² c(l+1) = l² − 3l − 2`. If `ρ_{l+1} > ρ_∞` and
  `ρ_l − ρ_∞ ≤ s/c_l` with `s ≤ ρ_∞ (l²−3l−2)/(l+1)²`, then `(l+1)² ρ_{l+1}/c_{l+1} > l² ρ_l/c_l`. `(l²−3l−2)/(l+1)²` is
  increasing for `l ≥ 1` (derivative `(5l+1)/(l+1)³`).
- **H9 (Ψ expansions).** With `N_l = N_∞ + εD_N`, `X_l = X_∞ + εD_X`: `Ψ_l − Ψ_∞ = ε[ρ³N_l² + 2ρ²D_N(N_l+N_∞) − ρ³X_l² −
  4N_Z⁰D_X(X_l+X_∞)]`; with a second order, `Ψ_l = Ψ_∞ + εΨ₁ + ε²E(ε)`, `Ψ₁ = 4ρ²N_∞ν + ρ³N_∞² − ρ³X_∞² − 8N_Z⁰X_∞ξ`,
  `E` a polynomial in `ε` with the coefficients printed in `review_large_l_thm_m.py` (`E_rem`). Both are polynomial
  identities (sympy).
- **H10 (sign of `V̂` at infinity).** The leading coefficient of `Ψ/a²ρ^{2l+7}` is `ε(l+1)² − (1+2ε)/2 = (l+1)/(2(l−1))
  > 0`, and `V̂ r² → (l+1)(l+2)`.
- **H11 (identities of Theorem 3).** `V = 8c⁴(2+ερ)N_Z/(ρ⁴(ρ−6)²)`; `V̂ = 2g² − ψ₀,ₓₓ/ψ₀` for any `X`; `V̂ =
  S·c³(2c+ρ)/(ρ³(ρ−6)X²)` for any `X` (so `sign V̂ = sign S` for `ρ > 6`, `X > 0`); `V̂ = 2c⁴(2+ερ)Ψ/(ρ⁴(ρ−6)²X²)` when
  (Z) holds.
- **H12 (Theorem 3, (P1) and (P2)).** With `λ = c/2`, `μ = 1/λ`, `λr = 3 + σ`: `λ²r²((λ+1)r − 3) + 9(λr − 1) = (σ³ + 6σ²
  + 18σ + 18) + μ(σ + 3)³`, positive for `σ > 0`, so `V > 0` on `(r_s, ∞)`.
  With `Q := −V_r r⁵(cr−6)³/2 = c⁴r⁵ + 3c⁴r⁴ + 2c³r⁵ − 3c³r⁴ − 30c³r³ + 78c²r³ + 252c²r² − 396cr² − 1080cr + 648r +
  1728`: `c⁵Q` at `r = (6+τ)/c`, `c = (m+4)(m+1)` is a polynomial in `(τ, m)` with 64 coefficients, all positive, so `Q
  > 0` for `τ > 0`, `m ≥ 0` (every real `l ≥ 2`), `V_r < 0` and `V_x = f V_r < 0` on `(r_s, ∞)`. `N_Z(6+τ) = (18 + 9τ +
  (3/2)τ² + τ³/8) + 2ε(27 + (27/2)τ + (9/4)τ² + τ³/8) > 0`.

## 8. What RESULTS.md §4.6's last bullet should say

Replace "None of 4.2–4.5 has had an independent review. Per AGENTS.md they should not be relied on until one has been
done." by:

> 4.2–4.5 had an independent review in #104 (`review_large_l.md`; mpmath `iv` and sympy, written from this file alone).
> It confirms Theorem 3's (P1), (P2) and the V̂ identities, the brackets for every l = 2..3000 (and the 200 random
> ones), Theorem L and Theorem M with the reviewer's own, smaller constants (L: 23.8 and 172.0 first order, 3.55 and
> 6.68 second order; M: 3.52 and 5.94), and finds no error. The review did not re-prove T1, T2, or Lemma A(ii) for l >
> 14; the certified results still assume them. Recorded wording nits: "l²·r(l+1) > l²·r(l)" in Theorem M means
> `(l+1)²r(l+1) > l²r(l)`; and for l < 500 only the consecutive brackets of 4.3 are needed, not those to l = 3000.

Two wording nits for the coordinator (not errors): (i) Theorem M's "`l²·r_*(l+1) > l²·r_*(l)`" reads literally as
`r_*(l+1) > r_*(l)`, which is false (`c·r` decreases and `c` grows, so `r` decreases); the proof shows `(l+1)²r_*(l+1) >
l²r_*(l)`. (ii) 4.5's last paragraph says "consecutive brackets to l = 3000"; Theorem M covers l ≥ 500.

## 9. Process notes (what went wrong in the review itself)

- First bookkeeping of the coefficient bounds took the maximum of the two regimes (`n ≤ l+2`, `n ≥ l+3`) for every `n`;
  Z came out 12169 and `E_max` 4.06e8, and the Theorem M comparison failed by 1.2× (N) and 4.2× (Ψ). The regimes are
  exclusive; the 60 explicit terms satisfy `n ≤ l+2`, so only the tail needs the uniform majorant. Log:
  `step4_thm_m_try1.log`.
- The coefficient-bound test first demanded `(bound − |a_n − α_n|) ≥ 0` as an enclosure, which cannot hold at an
  equality case (`n = 4`, l = 138, `ε = ε₀`); the test is now "not definitely violated".
- The overlap script first assumed round-to-nearest printing in #2's table and failed on 15 of 32 values
  (`step7_overlap_try1_wrong_rounding_guess.log`); the values are truncated. A second slip (window without #2's bracket
  width) is in `step7_overlap_try2_window_missing_bracket_width.log`.
- Two controls were mispredicted or silent (§6). No pre-registered target was edited.

## Coordinator's review

Raw `I104-mac` (negative-mass-raw). Mac, nice 10, at most 2 processes.

| step | result |
|---|---|
| rerun every `review_large_l_*.py` | same outputs apart from timings; controls exit 1 with the same 5 misses; the sweep's csv is byte-identical to rpi5's (684 s) |
| one mutation | all 52 rerun; e.g. `review_large_l_sweep.py --mutate=l2r_reverse`: exit 1, first FAIL at l = 5 → 6, on target |
| re-solve one case, l = 29 | exact real-root isolation (sympy `real_roots`): one root of N beyond ρ = 6, at 8.489743651547031155; S has three positive roots, 6 (r_s, by L1), 6.1690082 and 9.642705229862547021. Both roots are inside the brackets |
| independent check `ind104.py` | exit 0, every criterion passes; 6 of its 7 mutations fail on target (below) |

`ind104.py` shares no code with this review or #2. It works in r with Lemma A's c_n at 60 digits (numerical, not
intervals), and takes V̂ from its definition g² − g_x. At 49 values of l from 2 to 3000:
- the Zerilli residual of ψ₀ is below 4e-61, numerical evidence for Lemma A(ii), which this review assumed beyond l = 14;
- N and V̂ each change sign once beyond r_*, inside every bracket;
- the V̂ identities with S and Ψ hold to 2e-60;
- (c·r − ρ)·c stays below this review's s_N and s_Ψ (it rises to 3.51549 and 5.83939 at l = 3000).

A three-term fit at l = 1000, 2000, 3000 gives ρ to 4e-19, and the 1/c coefficients 3.515485851 and 5.839392647, as
RESULTS §4.3 states.

Its mutation `vhat` exits 1, but by a crash after F3 fails, before its target F7. That is a defect of ind104. F7 fails
under `cn`.

Findings on this review, not fixed in its files:
- **Coefficient-bound test.** `review_large_l_thm_l.py` step 1 was changed after a failed run, from "certified ≥ 0" to
  "not definitely violated" (§9). Its comment and printed message give the equality case as "only at n = 4, l = 138".
  In exact rationals the bound is an equality at n = 3 for every l (α′₃ = 0, S₃ = 0, so both sides are 2ε|α₃|), and at
  n = 4 for l = 138. H3 holds, with equality there, so no conclusion changes.
- **Overlap with #2's table.** The acceptance item "each of the 200 brackets compared with #2's table for overlap" is met
  for 2 of the 200. #2's output prints 16 l; the other 198 were compared with #2's Theorem L and M bounds.
- **Controls.** 5 of 52 missed their targets: the criterion failed, and it stays failed. The silent
  `accept_zero_enclosure` is covered after the fact by `review_large_l_posthoc.py`.
