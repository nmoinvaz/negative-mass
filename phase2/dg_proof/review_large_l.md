# Independent review of #2's large-l proofs (RESULTS.md §4.2 to §4.5)

Issue #104, runner `negative-mass-rpi5-i104` (rpi5). Written from `RESULTS.md` alone; step 7 (comparison with #2's output files) came last.
Verdict words: **confirmed**, **gap**, **error**. Status words as in AGENTS.md. Raw output: `negative-mass-raw/I104-rpi5/`.

TBD-VERDICT

## 1. What was run

Commands run from this directory (`phase2/dg_proof`), one process, `nice -n 10`, `OMP_NUM_THREADS=1`:

TBD-RUNS

## 2. Independence: what is shared with #2

Not opened before step 7: `rescaled.py`, `vhat_structure.py`, `threshold_isolation.py`, `large_l_proof.py`, `monotone_proof.py`, `common.py`, every other script of
`phase2/` and every `*_output.txt` of those. Certified parts use mpmath `iv` (outward-rounded intervals, 192 bits); #2 used `fractions.Fraction`. Identities use sympy 1.14.0.

Shared with #2: the formulas typed from RESULTS.md (`a_n(ε)`, `X`, `N`, `N_Z`, `Ψ`, `S`, the static Zerilli `V`, Lemma A's closed form for `χ`, the statements of T1/T2),
and the *structure* of the proofs of Theorems L and M, because RESULTS.md describes it (coefficientwise bounds, a perturbation of `N` and `Ψ`, slopes, localisation;
second-order expansion). Not shared: every line of code, the arithmetic, the tail bounds, the windows, the way the constants are assembled, and the comparison step of
Theorem M (derived afresh below, with a smaller loss than RESULTS' `(l+1)ν/2 > |ζ|`).

## 3. Assumed, not re-proved here

- **T1, T2** (RESULTS §3): `N > 0` on `[0, r_s]` and exactly one root of `N` beyond `r_s`. A bracket in which `N` changes sign is *the* `c·r_*` only through T2. This review checks T2 by exact Sturm counts for l = 2..16 and shows each bracket holds exactly one root, but T2 for l ≥ 17 is #2's Descartes argument, not reviewed.
- **Lemma A(ii)** (ψ₀ solves the static Zerilli equation for every l): exact for l = 2..14 here (sympy, own χ); for every l it is #2's symbolic proof, not reviewed. The `S`-based brackets use no (Z), and agree with the `Ψ`-based ones at every l run.
- **Theorem 3** is used in §4.3, §4.4 and §4.5 to say that the zero of `Ψ` found beyond `r_*` is `r_neg`. Its steps are reviewed in §4 below.
- The textbook facts of §4.2 (intermediate value theorem, Taylor with Lagrange remainder, Weierstrass' product inequality).

TBD-STEPS

## 7. Hand lemmas, stated for a later formalization in Lean

Each is proved on paper in this review; the polynomial inequalities are certified by sympy (all coefficients of the shifted polynomial are non-negative), the rest are elementary.
Notation: `c = (l+2)(l−1)`, `ε = 1/c`, `q_a = (a−1)(a+2)`, `e₁(n) = Σ_{a=0}^{n−3} q_a = (n−2)(n²−4n−3)/3`, `S_n = e₁(n) + 2`,
`g_n = 1/(2ⁿ((n−1)!)²)`, `α_n = (n−2)(n−4)g_n`, `α′_n = 3(n−2)(n−3)g_n`, `α¹_n = α′_n − α_n e₁(n)`, and `a_n = 0` for `n ≥ l+3`.

- **H1 (zero crossing).** If `u ∈ C¹[a,b)` and `u′(z) > 0` at every zero `z` of `u`, then `u` has at most one zero.
  *Proof.* Take the least zero `z₁` (zeros form a closed set). `u′(z₁) > 0`, so `u > 0` just right of `z₁`. A second zero has a least value `z₂ > z₁` with `u > 0` on `(z₁, z₂)`, so `u′(z₂) ≤ 0`, a contradiction.
- **H2 (Theorem 3, step (B)).** Let `V > 0` and `V_x < 0` on `I = [x₀, ∞)`, `g ∈ C¹(I)` with `g′ = V − g²`, `g(x₀) = 0`, `g > 0` on `(x₀, ∞)`, and `h = √(V/2)`. Then `Ṽ := 2g² − V = 2(g−h)(g+h)` has `sign Ṽ = sign(g−h)` on `(x₀, ∞)`, `Ṽ(x₀) < 0`, and `g − h` has at most one zero. If `Ṽ > 0` near `∞`, `Ṽ` has exactly one zero in `(x₀, ∞)`, it is simple, `Ṽ < 0` before it and `> 0` after it.
  *Proof.* At a zero `z` of `g − h`: `g′(z) = V − g² = V/2 > 0` and `h′ = V_x/(4h) < 0`, so `(g−h)′(z) > 0`; apply H1.
- **H3 (first-order coefficient bound).** For `l ≥ l₀`, `ε₀ ≥ 1/c(l₀)`, `3 ≤ n ≤ l+2`: `|a_n − α_n| ≤ ε[α′_n(1+2ε₀) + |α_n| max(2, S_n(1+2ε₀))]`; for `n ≥ l+3`: `|a_n − α_n| ≤ ε n² |α_n|`.
  *Proof.* `a_n = (α_n + εα′_n) P_n`, `P_n = (1+2ε) Π_{a=2}^{n−3}(1 − q_a ε)`, each factor in `[0,1]` because `q_a ≤ q_l = c` for `a ≤ l`. Weierstrass: `Π(1 − x_a) ≥ 1 − Σx_a` on `[0,1]`, so `P_n ∈ [(1+2ε)(1−εS_n), 1+2ε]`, hence `|P_n − 1| ≤ ε max(2, S_n(1+2ε))`, and `a_n − α_n = α_n(P_n − 1) + εα′_n P_n`. For `n ≥ l+3`, `n² > c`, so `1 < εn²`.
- **H4 (second-order expansion).** For `n ≤ l+2`: `a_n = α_n + εα¹_n + ε²ω_n` with `|ω_n| ≤ |e₁|α′_n + (|α_n| + ε₀α′_n)(e₁+4)²(1+2ε₀)/2`; for `n ≥ l+3`: `|ω_n| ≤ n⁴|α_n| + n²α′_n + n²|α_n||e₁|`.
  *Proof.* `P_n(ε) = 1 − e₁ε + P_n″(ξ)ε²/2` (Lagrange), `P_n″ = Σ_{a≠b} q_a q_b Π_{c≠a,b}(1 − q_c ξ)`, each product `≤ 1+2ε`, `Σ|q_a| = e₁+4`. For `n ≥ l+3`: `a_n = 0` and `1 < εn²`.
- **H5 (majorants and tails).** For `n ≥ 3`: `B_n ≤ n⁵g_n`, `|α¹_n| ≤ n⁵g_n`, `Ω_n ≤ n⁸g_n` (`B_n`, `Ω_n` the bounds of H3, H4, `ε₀ ≤ 1/1000`). For `ρ ≥ 1` and `n₁ ≥ 1`: `Σ_{n>n₁} n^D ρ^{n−k} g_n ≤ U(n₁+1)/(1−r)`, `U(n) = n^D ρ^{n−k} g_n`, `r = ((n₁+1)/n₁)^D ρ/(2n₁²) < 1`.
  *Proof.* The polynomial inequalities are `P(n) ≥ 0` for `n ≥ 5` (all coefficients of `P(5+t)` are `≥ 0`, sympy) and for `n = 3, 4` by direct evaluation; the tail by the ratio `U(n+1)/U(n) = ((n+1)/n)^D ρ/(2n²)`, decreasing in `n`.
- **H6 (localisation).** If `F ∈ C¹[ρ̄−δ, ρ̄+δ]`, `F(ρ̄) = 0`, `F′ ≥ μ > 0` there, and `|F_l − F| ≤ ηε < μδ` there, then `F_l` has a zero in `[ρ̄ − ηε/μ, ρ̄ + ηε/μ]`. (Same with `F′ ≤ −μ`.)
  *Proof.* `F(ρ̄ ± t) = ±(≥ μt)`; at `t = ηε/μ`, `F_l ≥ 0` at `ρ̄ + t` and `≤ 0` at `ρ̄ − t`.
- **H7 (monotone comparison).** If `N_l`, `N_{l+1}` each change sign exactly once on `(6, ∞)`, from `+` to `−`, at `ρ_l`, `ρ_{l+1}`, then `N_{l+1}(ρ_l) < 0` implies `ρ_{l+1} < ρ_l`. (For `Ψ`, from `−` to `+`, `Ψ_{l+1}(ρ_l) > 0` implies the same, on the region beyond `r_*`.)
- **H8 (`l²r` criterion).** `c(l+1) − c(l) = 2l+2` and `(l+1)² c(l) − l² c(l+1) = l² − 3l − 2`. If `ρ_{l+1} > ρ_∞` and `ρ_l − ρ_∞ ≤ s/c_l` with `s ≤ ρ_∞ (l²−3l−2)/(l+1)²`, then `(l+1)² ρ_{l+1}/c_{l+1} > l² ρ_l/c_l`. `(l²−3l−2)/(l+1)²` is increasing for `l ≥ 1` (derivative `(5l+1)/(l+1)³`).
- **H9 (Ψ expansions).** With `N_l = N_∞ + εD_N`, `X_l = X_∞ + εD_X`: `Ψ_l − Ψ_∞ = ε[ρ³N_l² + 2ρ²D_N(N_l+N_∞) − ρ³X_l² − 4N_Z⁰D_X(X_l+X_∞)]`; with a second order, `Ψ_l = Ψ_∞ + εΨ₁ + ε²E(ε)`, `Ψ₁ = 4ρ²N_∞ν + ρ³N_∞² − ρ³X_∞² − 8N_Z⁰X_∞ξ`, `E` a polynomial in `ε` with the coefficients printed in `review_large_l_thm_m.py` (`E_rem`). Both are polynomial identities (sympy).
- **H10 (sign of `V̂` at infinity).** The leading coefficient of `Ψ/a²ρ^{2l+7}` is `ε(l+1)² − (1+2ε)/2 = (l+1)/(2(l−1)) > 0`, and `V̂ r² → (l+1)(l+2)`.
- **H11 (identities of Theorem 3).** `V = 8c⁴(2+ερ)N_Z/(ρ⁴(ρ−6)²)`; `V̂ = 2g² − ψ₀,ₓₓ/ψ₀` for any `X`; `V̂ = S·c³(2c+ρ)/(ρ³(ρ−6)X²)` for any `X` (so `sign V̂ = sign S` for `ρ > 6`, `X > 0`); `V̂ = 2c⁴(2+ερ)Ψ/(ρ⁴(ρ−6)²X²)` when (Z) holds.

TBD-REST
