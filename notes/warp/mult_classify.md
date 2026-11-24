# Multiplicative separable shifts u = X(x) R(ρ) + v with det R₃ = 0

November 2026, issue [#160](https://github.com/nmoinvaz/negative-mass/issues/160). Script
[`mult_classify.py`](mult_classify.py), its stdout [`mult_classify_output.txt`](mult_classify_output.txt). It follows
[`det_classify.md`](det_classify.md) §5 (#80) and uses #65's lemma, as #80 did. Imported from `det_classify.py`:
`E_mult_poly` and `E_add_poly`, only as comparison targets. Shared with it: the formula for D (retyped). Status words
as in AGENTS.md: **proved** (by hand), **checked symbolically**, **numerical evidence**, **not decided**.

## Answer

**No solution.** On a box I × J in (x, ρ), ρ > 0, with X X′ R R′ ≠ 0 (X, R real and C²) there is no u = X(x) R(ρ) + v
with D = det R₃ = 0, for v = 0 (static) and for every v ≠ 0 (traveling). The multiplicative solutions of D = 0 that
remain have X′ ≡ 0 (#65's b(ρ), which fails closure) or R′ ≡ 0 (planar). So the multiplicative row of
[`det_classify.md`](det_classify.md)'s table becomes: no solution with X′R′ ≠ 0, and there is nothing to send through
closure.

| case | verdict | status |
|---|---|---|
| v = 0 | no solution with X X′ R R′ ≠ 0 | **proved** (§3), identities and finite cases **checked symbolically** |
| v ≠ 0, \|v\| ≠ 1 and \|v\| = 1 | no solution with X X′ R R′ ≠ 0 | **proved** (§4), identities and finite cases **checked symbolically** |
| planted additive u = X + R | u = c x + g(ρ), g from the Bessel equation, found for c from 3/10 to 2 (§5) | **checked symbolically**, **numerical evidence** at 130 digits |
| perturbed additive shifts | rejected: exactly (`add_perturbed_exact`), but the thresholded `add_perturbed` **failed** at one point (§6) | exact; one failed criterion |

This is **not** the unit ideal of the prolongation that §5 of #80 started, and it is not a Gröbner computation of it.
That route did not finish (§6); the answer comes from a different exact route, the separation of E into a sum of
products. Whether the saturated prolongation ideal is the unit ideal was not computed. Its consequence, no solution,
is what is proved here.

## 1. E and its separation

With S = R″ + R′/ρ, the jets x: (X₀, X₁, X₂) = (X, X′, X″) and ρ: (R₀, R₁, S₀) = (R, R′, S), and u = X₀R₀ + v, the
entries of #65's block are u_x = X₁R₀, u_ρ = X₀R₁, u_xx = X₂R₀, u_xρ = X₁R₁, Δ_ρu = X₀S₀, and E = 8D/R₁² is a
polynomial of degree 10 in (X₀, X₁, X₂, R₀, R₁, S₀, v) that equals `det_classify.E_mult_poly` (`e_from_d`, **checked
symbolically**; also against numerical derivatives of u at ten points, relative error 4e-131). With the log variables
p₁ = X′/X, m = p₁², β = (ln X)″, c = (R′/R)², s₂ = S/R, U = X R, u = U + v, W = 3U + v:

    E/X² = s² − δ² − τ²,   s = U²c + m(W² + 1)/2,   δ = 2uU(β + m) + m(4U² + 1 − W²)/2,   τ = U s₂ + m W.

**Static, v = 0.** E = R²X⁴ [a T₋T₊ + m T₋ − (s₂ + 3m)²], a = X²R², T₋ = c + 5m − 2β, T₊ = c + 4m + 2β (#80's A = T₋T₊
and B). Hence

    Ê = Σ_{i=1}^{6} f_i(x) g_i(ρ),   f = (X², X²m, X²Θ, m, Ψ, 1),   g = (R²c², 9R²c, R², c − 6s₂, 1, −s₂²),
    Θ = 20m² + 2mβ − 4β²,   Ψ = −4m² − 2mβ,

an exact identity (`decomp_static`). The coefficient matrix of E in x-monomials × ρ-monomials has rank exactly 6 (9 ×
7 monomials), so six terms are minimal.

**Traveling, v ≠ 0.** E/X² = Σ_{k=1}^{4} U^k E_k with

    E₁ = −2mv[(1 − v²)(β + m) + s₂],   E₂ = Φ₂ + m((v² + 1)c − 6s₂) − s₂²,   E₃ = 2v(Φ₃ + 3mc),   E₄ = T₋T₊,
    Φ₂ = −4v²β² + (6v² − 2)βm + (12v² − 4)m²,   Φ₃ = −4β² + 3βm + 13m²,

and E/X² = Σ_{i=1}^{10} F_i G_i with F = (−2v(1−v²) X m(β+m), −2v X m, X²Φ₂, X²m, X², 2v X³Φ₃, 6v X³m, X⁴, X⁴m, X⁴Θ)
and G = (R, R s₂, R², R²((v²+1)c − 6s₂), −R²s₂², R³, R³c, R⁴c², 9R⁴c, R⁴) (`trav_decomp`, exact for symbolic v). The
rank of E is 10 at v = 1/3, 2, 5/7, 3/2, −2, **9 at v = ±1** (F₁ = 0 there) and 6 at v = 0 (observation, printed).

## 2. Lemma L (linear algebra, **proved**)

If Σ_{i≤r} F_i(x) G_i(ρ) = 0 on a box, then the vectors F(x) and G(ρ) in ℝʳ are orthogonal for all x, ρ, so span{F(x)}
⟂ span{G(ρ)} and

    d_F + d_G ≤ r,   d_F = dim span of the functions F_i of x,  d_G likewise.

`rank_lemma` plants it: five random draws with d_F = 4, d_G = 2, r = 6 satisfy the identity, and with d_G = 3 the sum
is nonzero at a sample point in all five.

## 3. Static case v = 0

All three lemmas use only that X, R are real and X X′ R R′ ≠ 0 on the box.

**Lemma X (x-side).** If m is not constant on any subinterval, then d_f ≥ 4. *Proof.* Pass to a subinterval with m′ ≠
0 (m is C¹), so m is a coordinate and E_x := X² is a function of m. If 1, m, E_x, mE_x are independent, d_f ≥ 4.
Otherwise c₀ + c₁m + E_x(c₂ + c₃m) = 0, E_x is not constant (E_x′ = 2p₁E_x ≠ 0), so E_x = ρ(m) = (am + b)/(cm + d), ad
− bc ≠ 0, a Möbius function. From E_x′ = 2p₁E_x and m′ = 2p₁β: β = ρ/ρ′ = (am + b)(cm + d)/(ad − bc), a quadratic
polynomial in m. For c ≠ 0: span{1, m, ρ, mρ} = span{1, m, t}, t = 1/(cm + d), of dimension 3, and Ψ = −4m² − 2mβ is a
polynomial of degree ≤ 3 with m³ coefficient −2ac/(ad − bc), and −2 for the m² coefficient at a = 0, so Ψ ∉ span{1, m,
t} (polynomials in it have degree ≤ 1). For c = 0, span{1, m, ρ, mρ} = span{1, m, m²} and ρΘ is cubic with m³
coefficient 18a/d ≠ 0. In both cases d_f ≥ 4. ∎ (`lemma_x`: the identities are **checked symbolically**, and the exact
rank of the six rational functions is ≥ 4 on 1170 parameter sets, minimum 4.)

**Lemma R (ρ-side).** d_g ≥ 3. *Proof.* If d_g ≤ 2, the span is ⟨1, R²⟩, so R′² = α + γR² (g₂ ∝ R′²) and R′⁴/R² = (α +
γR²)²/R² ∈ ⟨1, R²⟩ forces α = 0 (polynomial identity in y = R²), so R′/R = κ and R = a e^{κρ}; then s₂ = κ² + κ/ρ and
s₂² contains κ²/ρ², which is not in span{1, e^{2κρ}}; κ = 0 contradicts R′ ≠ 0. ∎ (`lemma_rho`; numerical ranks of the
six g_i for 13 sample R are 4 to 6.)

**Lemma K.** If m ≡ λ² is constant (X = X₀e^{λx}, β = 0), then d_f = 2 (f spans ⟨X², 1⟩), so the coefficient of X²
must vanish: R²(c + 4λ²)(c + 5λ²) = 0, impossible for real c ≥ 0, λ ≠ 0. ∎ (`const_m`)

**Theorem S (v = 0, proved).** If m is not constant, d_f ≥ 4 and d_g ≥ 3, so d_f + d_g ≥ 7 > 6, contradicting Lemma L.
If m is constant, Lemma K. ∎

## 4. Traveling case v ≠ 0

Let m not be constant (the constant case is `trav_exp`: the coefficient of X⁴ is R⁴(c + 4λ²)(c + 5λ²) > 0, independent
of X, X², X³).

1. **d_F ≥ 4.** F₂, F₄, F₇, F₉ are nonzero multiples of mX, mX², mX³, mX⁴ (v ≠ 0), independent since X is a
   coordinate. Put M₄ = their span.
2. **d_F ≤ 5.** For v² ≠ 1 the separation has ten terms and d_G ≥ 5 (`trav_dg5`): if d_G = 4 then S, T, S², RP, P²,
   R²P (P = R′² as a function of R, T = (v²+1)P − 6RS) lie in span{R, …, R⁴}, which forces P = a₁R + a₂R², S = s₁R +
   s₂R²; the consistency identity q² = (P_R/2)q − Pq_R, q = S − P_R/2 = R′/ρ (from ρ = √P/q and dρ/dR = 1/√P) then
   forces q = 0, i.e. R′ = 0. For v² = 1 the separation has nine terms and the same argument with span{R², R³, R⁴}
   gives d_G ≥ 4 (`trav_dg4_unit`); so d_F ≤ 9 − 4 = 5. In both cases W_F := span{F_i} = M₄ ⊕ ℝw, or M₄.
3. **m is rational.** X² = F₅ and X⁴ = F₈ lie in W_F. Write X² = mXp + y₁w, X⁴ = mXq + y₂w with p, q polynomials of
   degree ≤ 3. If y₁ = y₂ = 0, then m = X/(p₀ + p₁X). Otherwise m = X(y₂ − y₁X²)/r(X), r = y₂p − y₁q ≢ 0 of degree ≤
   3.
4. **m is a Laurent polynomial in X.** In the second case W_F ⊂ ℝX^{2 or 4} + mX·ℝ[X]_{≤3}, whose elements have pole
   order ≤ k at a pole X₁ ≠ 0 of m of order k. But F₁₀ = X⁴Θ has pole order 2k + 2 > k (`trav_pole`: the leading
   coefficient is −k²X₁²μ²). So the only poles are at X = 0. In the first case: p₁ = 0 gives m = X/p₀ (a monomial), p₀
   = 0 gives m constant, and for p₀p₁ ≠ 0 the span has d_F ≥ 6 for every v ≠ 0, contradicting step 2 (`trav_y0`: 60
   six-by-six minors of the coefficient matrix generate, with the saturation t p₀p₁v − 1, the unit ideal).
5. **Exponents.** If m has top exponent b at X = ∞ (bottom exponent a at 0), the leading exponents of the F_i in W_F
   lie in {b+1, …, b+4} plus at most one more value (a space of Laurent polynomials has as many distinct leading
   exponents as its dimension). The leading monomials of F₁, F₃, F₅, F₆, F₈, F₁₀ for m = cX^b are exact
   (`trav_exponents`), and counting the exponents outside {b+1, …, b+4} gives b ∈ {0, 1} for v² ≠ 1 and b ∈ {−2, −1,
   0, 1} for v² = 1; the same for a. So m is a Laurent polynomial with exponents in [−2, 1].
6. **No Laurent m works** (`trav_laurent`, for symbolic v). For m = k₋₂X⁻² + k₋₁X⁻¹ + k₀ + k₁X, the coefficient of X⁶
   in E/X² is 20R⁴k₁², so k₁ = 0; then the coefficient of X⁴ is R⁴(c + 4k₀)(c + 5k₀), so c = c₁² is the constant −4k₀
   or −5k₀ (k₀ ≠ 0 since c > 0), R = a e^{κρ}, s₂ = κ² + κ/ρ, and the coefficients of E/X² at X⁻², …, X³, as
   polynomials in e^{κρ} and 1/ρ, together with κ a v ≠ 0, generate the unit ideal in both cases.

**Theorem T (v ≠ 0, proved).** Steps 1 to 6 give a contradiction for every real v ≠ 0 and every real C² pair with X X′
R R′ ≠ 0. ∎

By hand are: the dichotomy of step 3, the pole-order argument of step 4, the valuation argument of step 5 and the
logic joining the steps. Checked by code are the identities, the finite exact computations (ranks, Gröbner unit
ideals, exponent tables for b ∈ [−15, 15] with the inequalities that exclude all other b) and the numerical ranks. The
rank bound alone does not decide v ≠ 0: for R = ρ² and v = 1/3 the ρ-side has d_G = 5 (`trav_gap`), which is why steps
3 to 6 are needed.

## 5. Planted truth through the same separation

Additive u = X + R: E_add = (h + X′²)² − (X + R)²X″² − S²/4 (h = R′²/2) equals `det_classify.E_add_poly` and separates
into five terms, f = (1, X′², X″², XX″², X′⁴ − X²X″²), g = (A² − B, 2A, −R², −2R, 1), rank 5 (`add_decomp`). The
x-curves with all f_i constant (d_f = 1) are found by a Gröbner elimination of {f_i − k_i} with X free: X″² = 0, X′² =
k₂ = c² free, k₅ = k₂² (`add_stratum`); then f₀·g = (A + c²)² − B = 0, the Bessel equation R″ + R′/ρ = ε(R′² + 2c²) of
#80 §3. The solutions R = −ε ln(J₀(z) + 0.3 Y₀(z)), z = √2cρ, satisfy it for c ∈ {3/10, 1/2, 1, 3/2, 2}, ε = ±1 with
relative residual ≤ 1e-130 at 50 points and d_G = 4, so d_F + d_G = 5 = r (`add_planted`): the bound is saturated, the
additive family is not excluded by Lemma L. With X = cx + x²/10 the same R has d_F ≥ 2, d_F + d_G = 6 > 5: rejected
(`add_xperturbed`). The subcases #80 recovers: with X′ = X″ = 0, E = X⁴(X²R′⁴ − S²) = f₁g₁ + f₂g₂ with f constant, and
R = −(ε/X)ln(C₁ + C₂ ln ρ) solves it; the planar case R′ = S = 0 has D = 0 (`mult_subcases`).

## 6. What failed, what was silent, what did not finish

**`add_perturbed` failed (a criterion that stays failed).** The README asked that for R + δ, δ ∈ {ρ²/10, ρ³/10}, the
relative residual |f₀·g|/Σ|f_ig_i| be ≥ 1e-3 at all 100 points (5 values of z = √2cρ for each of c ∈ {3/10, 1/2, 1,
3/2, 2}, ε = ±1, two δ). It is 3.72e-4 at one point, c = 1, ε = +1, δ = ρ³/10, z = 1.6 (ρ = 1.131); the next smallest
is 5.7e-3 and 99 of 100 are above. Nothing was changed (threshold, δ, points), by the coordinator's ruling. The reason
is the choice of the perturbation and the point, not a property of E: at that point the unperturbed (A + c²)² and B
are both 2.69 (A = R′²/2, B = S²/4, γ = g′ = 1.132) and the perturbation δ = ρ³/10 moves both to about 4.62 (4.62096
and 4.62342), so E_pert = (A + c²)² − B = −2.45e-3 is nonzero but small against the sum of the term sizes, 6.60. The
other perturbation at the same point gives 5.7e-2. A criterion in a relative residual is small wherever the
perturbation moves the two sides together. The exact check `add_perturbed_exact` (written after the failure, README
A3, not tuned) evaluates E_pert = (R₁²/2 + c²)² − S₀²/4 for u = cx + g + δ symbolically and at three fixed exact
points, and it is nonzero in all cases. The script exits 1 because of `add_perturbed`.

**A silent control.** `m-unit-sign` (README A6) was predicted silent by a hand computation before it was run (README
A7: with either sign of the P q_R term the identity forces a = 0), and it was silent; it is reported and not counted.
Its replacement `m-unit-v` fires.

**Dropped control.** `m-perturb-zero` targeted `add_perturbed`, which fails at baseline; it is dropped (A3), as #80
dropped `conical-limit`.

**The route the issue named did not finish.** The prolongation of §5 of #80 in the jets gives one irreducible
polynomial Γ (345 terms, degree 15 in (U, p₁, p₂, c₁, s₂, 1/ρ, v); at v = 0 with U² eliminated, 190 terms, degree 16
in (p₁, p₂, c₁, s₂, 1/ρ), where #80 quotes 17 counting ρ), and its x- and ρ-derivatives have 526 and 720 terms; they
do not vanish on Γ = 0 at random real points. sympy's `groebner` of Γ and its x-derivative in four variables modulo
32003 was stopped by hand after 7 min (exploration, README). That the saturated prolongation ideal is the unit ideal
was therefore not computed.

README amendments A1 to A7 (all before the run they concern): A1 read the criterion "X″ = 0 in the ideal" of
`add_stratum` as X″² ∈ ideal, corrected the `m-stratum` change (dropping f₄ alone is silent) and scoped `m-add-freq`;
A2 records `add_perturbed`; A3 the coordinator's ruling; A4 the baseline-failure rule for controls; A5 to A7 added the
traveling checks.

## 7. Checks and controls

22 checks (21 pass, `add_perturbed` fails), 21 controls that fire on their named target and on no other check (the
baseline failure `add_perturbed` is excluded from "other"; one more control is declared silent, §6). The run of record
exits 1 because of `add_perturbed`. Per check, in seconds on a Raspberry Pi 5 (4 cores shared, one process at nice
10): `e_from_d` 0.3, `lemma_x` 1.8, `add_planted` 20.1, `add_perturbed` 11.5, `add_xperturbed` 13.4, `trav_y0` 24.6,
`trav_dg4_unit` 1.1, `trav_laurent` 0.8, `trav_dg5` 0.6, every other check 0.3 or less: 76 s in all. The whole run
with the 22 control runs took 29 min 31 s and 79 MB peak resident (`/usr/bin/time -v`).

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `e_from_d` | E = 8D/R′², sum-of-squares form, v = 0 form, numerical derivatives | `m-d-sign` | Y := u_xu_ρ − u u_xρ/2 |
| `decomp_static` | six-term identity, rank 6 | `m-theta` | 19m² in Θ |
| `rank_lemma` | Lemma L on planted data | `m-rank-off` | a non-orthogonal column |
| `lemma_x` | Möbius step, exact ranks ≥ 4 | `m-fdrop` | f₃, f₅ set to 0 |
| `lemma_rho` | α = 0, κ², ranks ≥ 3 | `m-gdrop` | g₄ := c, g₆ := 1 |
| `const_m` | R²(c + 4λ²)(c + 5λ²) > 0 | `m-const-sign` | c = −n/10 |
| `mult_subcases` | b(ρ) and planar | `m-sub-drop` | S in place of S² |
| `add_decomp` | E_add, rank 5 | `m-add-drop` | u = X₀ |
| `add_stratum` | constant-f elimination | `m-stratum` | f₄, f₅ dropped |
| `add_planted` | Bessel solutions | `m-add-freq` | J₀(cρ) |
| `add_perturbed` | relative residual ≥ 1e-3 (**fails**) | none (target fails at baseline) | |
| `add_perturbed_exact` | E_pert ≠ 0 exactly | `m-pexact-zero` | δ = 0 |
| `add_xperturbed` | d_F + d_G > 5 | `m-xperturb-zero` | X = cx |
| `trav_decomp` | ten-term identity, rank 10 | `m-trav-rank` | rank 9 expected |
| `trav_gap` | R = ρ² has d_G = 5 | `m-gap-d` | R = e^ρ |
| `trav_exp` | exponential X | `m-exp-sign` | c = −n/10 |
| `trav_dg5` | d_G ≥ 5 (v² ≠ 1) | `m-dg-sign` | + P q_R |
| `trav_dg4_unit` | d_G ≥ 4 (v² = 1) | `m-unit-v` (and `m-unit-sign`, silent) | v = 1/3 |
| `trav_exponents` | exponent bookkeeping | `m-exp-allow` | exponent 4 removed |
| `trav_y0` | d_F ≥ 6 for m = X/(p₀ + p₁X) | `m-y0-p1` | p₁ = 0 allowed |
| `trav_laurent` | Laurent coefficients, unit ideals | `m-laurent-coef` | 10 for 20 |
| `trav_pole` | pole order 2k + 2 | `m-pole-theta` | no −4β² |

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/mult_classify.py --no-controls     # the 22 checks, about 80 s
    python3 notes/warp/mult_classify.py --mutate=NAME      # one control
    python3 notes/warp/mult_classify.py --logdir=DIR       # the checks, then every control, logs in DIR

## 8. What the code proves, and what is assumed

Exact identities (sympy): E from D against #80's printed E; the sum-of-squares form; the six- and ten-term
separations; the leading monomials in `trav_exponents`; the pole coefficients in `trav_pole`; the polynomial
identities of `lemma_x`, `lemma_rho`, `trav_dg5`, `trav_dg4_unit`; Gröbner unit ideals in `trav_y0` and `trav_laurent`
and the elimination in `add_stratum`; exact ranks. **Numerical evidence** (mpmath, 130 digits): the numerical ranks of
`lemma_rho`, `const_m`, `trav_dg5`, `trav_dg4_unit`, `trav_gap`, `trav_exp`, `add_planted`, `add_xperturbed`, and the
numerical-derivative check in `e_from_d`. The numerical ranks only illustrate lemmas that are proved by exact steps
(d_g ≥ 3, d_G ≥ 5, d_G ≥ 4); no proof depends on them.

Assumed: #65's lemma (so D = 0 is the algebraic stage), the formula for D as typed (it is #65's, shared with #80), X
and R real and C² on a box with X X′ R R′ ≠ 0, sympy 1.14.0. Not tested: complex solutions (several steps are
algebraic and would extend, `lemma_x` and the unit ideals are over ℚ, but `const_m`, the sign c ≥ 0 and the use of a
real coordinate are real arguments). The independence of X⁴ from X, X², X³ and of E^iw^j for the functions e^{κρ}, 1/ρ
are used as facts about functions on an interval.

## 9. Lemmas for Lean (#135)

Statements only; ℝ throughout.

1. **Lemma L.** If Σ_{i<r} F_i(x)G_i(ρ) = 0 for all x ∈ I, ρ ∈ J, then dim span{F_i} + dim span{G_i} ≤ r (linear
   algebra: orthogonality of two spans in ℝʳ).
2. **Möbius step.** If 1, m, E, mE are linearly dependent functions of a coordinate m with E ≠ const, then E = (am +
   b)/(cm + d), ad − bc ≠ 0, and E′ = 2p₁E, m′ = 2p₁β give β = E/E_m.
3. **Lemma X.** The rank of {1, m, ρ, mρ, ρΘ, Ψ} for ρ = (am + b)/(cm + d), β = ρ/ρ_m is ≥ 4 whenever ad − bc ≠ 0
   (rational functions; coefficient comparison).
4. **Lemma R.** If d_g ≤ 2 then R′ = κR, and then s₂² ∉ span{1, e^{2κρ}} unless κ = 0.
5. **Consistency identity.** For R as a coordinate, R′² = P(R), q = S − P_R/2 = R′/ρ, one has q² = (P_R/2)q − Pq_R;
   for P = a R², S = s R² it forces a = 0.
6. **Valuation lemma.** A finite-dimensional space of Laurent polynomials has as many distinct leading exponents at ∞
   (and lowest at 0) as its dimension.
7. **Pole order.** m = μ(X − X₁)^{−k} with X₁ ≠ 0 gives Θ = 20m² + 2mβ − 4β² of exact pole order 2k + 2.
8. **Laurent step.** For symbolic v, the stated polynomial system in (k₋₂, k₋₁, a, κ, v) generates the unit ideal over
   ℚ (`trav_laurent`, `trav_y0`).
9. **Theorems S and T**, from 1 to 8.

## 10. What `det_classify.md` §5 and `two_scalars.md` should now say

For §5 of `det_classify.md` (coordinator's text), replace "I did not finish the prolongation. The multiplicative case
is not classified." by: "Decided in [`mult_classify.md`](mult_classify.md) (#160): for every real v there is no
solution of D = 0 with X X′ R R′ ≠ 0 on an open box (*proved*, steps *checked symbolically*); the relation Γ of degree
16 in (p₁, p₂, c₁, s₂, 1/ρ) was not needed." In its table, the multiplicative row becomes "no solution with X′R′ ≠ 0;
X′ ≡ 0 is #65's b(ρ), R′ ≡ 0 is planar". For `two_scalars.md`'s "Not covered" bullet: the multiplicative separable
shifts u = X(x)R(ρ) + v are decided (no solution with X X′ R R′ ≠ 0, *proved*); conical shifts and the non-separable
general PDE remain, as in `det_classify.md` §12. For #158's note, which decided the one-sign multiplicative case by F
= 0 alone, the mixed-sign case now agrees (no solution), by a different route.

Not covered here: shifts that are not separable (the general static PDE, time-dependent shifts beyond planar and the
traveling reduction), the conical family's global solutions (#167), closure on the Bessel family for every member
(#162).