# An independent interval review of #133's count of the Type III points of the tanh wall

November 2026, issue [#136](https://github.com/nmoinvaz/negative-mass/issues/136). Script
[`review_he_type3_count.py`](review_he_type3_count.py), its stdout
[`review_he_type3_count_output.txt`](review_he_type3_count_output.txt). Step 7, the comparison with #133's files:
[`review_he_type3_count_compare.py`](review_he_type3_count_compare.py) and
[`review_he_type3_count_compare_output.txt`](review_he_type3_count_compare_output.txt). It reviews
[he_type3_count.md](he_type3_count.md) (#133). Nothing here claims novelty.

Status words are those of AGENTS.md. **Checked symbolically** means an exact identity in sympy or in `Fraction`s.
**Certified** means a computer-assisted proof whose trusted code is named in the last section. **Numerical evidence**
means a computation with no enclosure.

## Answer

**Confirmed.** Interval code that shares nothing with `he_type3_count.py` and `he_type3.py` reproduces #133's count.
It builds E1, E2 and P from the frame table of he_type3.md section 1, uses exact rational arithmetic (`fractions.Fraction`,
no mpmath `iv`) with its own enclosures of exp, cosh and sinh, and finds the same roots, the same ends and the same
boxes. At Le's wall (σR_b = 8, v = 1/2) P has 4 roots and at Rodal's (σ = 4, R_b = 5, v = 1) 6; y = B′ − A′ has one
root at each wall, with no common root of E1 and E2 anywhere in |w| ≤ 10⁴. Exactly three roots of P per wall carry a common root
with 0 < w < 1 and T ≠ 0, at the radii of he_type3.md section 5 to within 4 × 10⁻¹². The other roots of P have their common
root at w = 58.459064 (Le) and at w = 1.030003, 4.1006549, 8.7172638 (Rodal), as in he_type3_count.md. **No difference
in the count**; the differences of method and two wording gaps are in section 6.

## 1. What is shared with #133 and what is not

Typed from text: the frame table of he_type3.md section 1 (the entries R₀₀, R₀ₓ, R₀ρ, Rₓₓ, Rₓρ, R_ρρ, M = ηR with
η = diag(−1, 1, 1), and the printed c₂ and E1), the tanh profile and (φ, A′, B′) = (vF, v r f′, v r² f″) of he_type3.md
section 5, and from he_type3_count.md the walls, the ranges [1/4, z_t] with z_t = 40 and 60, the counts, the radii and the
w and T of each root, used only as targets.

Not shared: any code. The script did not open, import or run `he_type3_count.py`, `he_type3.py` or an output file
of `notes/warp/` before step 7. Two slips are on record: my auto-memory held a summary of #133's method (the reduction to (v, C, z) with
z = 2σr, a weight argument at small z and a z⁹ dominance at large z), which is where the idea of the two end bounds came from, and I read lines 1 to 25 of the raw
`I133-rpi5/README.md` once, which restates the reduction. Every formula was re-derived and is checked below by an exact identity.

## 2. Method

**Profile** (**proved** by hand, **checked symbolically**, `profile`). With z = 2σr and C = cosh 2σR_b,

    f = [tanh σ(r + R_b) − tanh σ(r − R_b)] / (2 tanh σR_b) = (1 + C) / (cosh z + C).

So only (v, C, z) enter. With D = cosh z + C, c = cosh z, s = sinh z,

    φ = v (c − 1)/D,    A′ = −v z (1 + C) s/D²,    B′ = −v z² (1 + C)(2 + C c − c²)/D³,    T = (A′ + B′)/(v r²).

**Algebra** (**checked symbolically**, `frame`). E1 and E2 are the characteristic-polynomial combinations c₂² − 3c₁ and
27c₀ − c₂³ of M = ηR, built by sympy from the typed frame table with (X, γ, y) = (A′², A′, B′ − A′). They equal the printed c₂
and E1; E2 is a cubic form in (X, γ, y) with 10 terms; E1 and E2 have 17 and 42 terms. The resultant in w is
(27/256) φ³ (A′ − B′)³ P with P irreducible, of 167 terms, total degrees 10, 12, …, 18. The table of P equals an exact
5 × 5 Sylvester determinant at 6 rational points.

**Reduction** (**proved**, section 3 L3). A point (z, w) with 0 < w < 1 and E1 = E2 = 0 forces P(z) y(z) = 0, since
φ > 0. The count is a count of the roots of P and of y = B′ − A′ on z > 0, and the common roots at each.

**Intervals.** Endpoints are `Fraction`s, rounded outward to a 160-bit grid after every operation (320 bits for root
refinement and Krawczyk). exp is enclosed by halving the argument k times to |y| ≤ 1/2, a Taylor sum of degree N − 1
and the remainder 2|y|ᴺ/N!, and k interval squarings; cosh and sinh come from exp(±x); forward-mode derivatives are
intervals. One unit in the last place is the grid spacing at the endpoint; a computed exp enclosure is 1 grid
spacing wide (`exp`, 300 points against 120 digits of mpmath, 0 misses).

**Cover** `Cover`: bisect [1/4, z_t]; a leaf is excluded if the natural enclosure or the mean-value form excludes 0,
or the derivative enclosure excludes 0 and the end values have one strict sign; it holds exactly one root if the derivative
enclosure excludes 0 and the end values have opposite strict signs. A leaf narrower than 2⁻²⁴ with neither property is
open and fails the check. Root leaves are refined by interval Newton at 320 bits to width below 10⁻³⁰.

**Common roots.** At a root cell Z of P (or of y) the code takes the numerical common root of E1, E2 from mpmath (a
starting point only) and tests a Krawczyk box X in (z, w) of half-widths 10⁻²⁶ and 10⁻²⁰; K(X) ⊂ int X gives exactly one
zero in X. The rest of Z × [−10⁴, 10⁴] is bisected in w until every box either excludes a zero (the enclosure of E1 or of E2
excludes 0) or lies inside a verified box. An unaccounted box fails the check.

**Ends.** (0, 1/4]: φ = z²Φ, ψ = A′ + 2φ = z⁴Ψ, δ = B′ − A′ = z⁴Δ, and P = z²⁴ G (section 3 L8). [z_t, ∞): P = (τz)⁹ S
(section 3 L9). Both are proved for the whole range, not on a grid.

## 3. The hand steps as lemmas, for Lean (#135)

Nothing is formalized here. Each statement says what is proved by hand (**proved**), what is an identity checked by sympy
(**checked symbolically**) and what is left to interval computation (**computer**). The numbers of the computer parts are in section 4.

**L1 (profile).** For σ, R, r > 0, z = 2σr, C = cosh 2σR: [tanh σ(r + R) − tanh σ(r − R)]/(2 tanh σR) = (1 + C)/(cosh z + C).
*Proof.* tanh x − tanh y = sinh(x − y)/(cosh x cosh y) and 2 cosh x cosh y = cosh(x + y) + cosh(x − y) give
2 sinh 2σR/(cosh z + C). Divide by 2 tanh σR = 2 sinh σR/cosh σR and use 2 cosh²σR = 1 + C. ∎ **Proved.**

**L2 (closed forms).** With D = cosh z + C and e = e⁻ᶻ: φ = v(c − 1)/D = v(1 − e)²/(1 + 2Ce + e²); A′ = −v z(1 + C)s/D² = −τ α̂ and
B′ = −v z²(1 + C)(2 + Cc − c²)/D³ = τ z β̂, with τ = 2v(1 + C) z e, α̂ = (1 − e²)/(1 + 2Ce + e²)²,
β̂ = [(1 + e²)² − 2Ce(1 + e²) − 8e²]/(1 + 2Ce + e²)³. **Checked symbolically** (`tail:*`).

**L3 (reduction).** If E1(z, w) = E2(z, w) = 0 for a real w, then (27/256) φ³ (A′ − B′)³ P = 0 at z, so P(z) = 0 or y(z) = 0.
*Proof.* A common root makes the Sylvester determinant of formal degrees (2, 3) vanish; use L4 below and φ = v(c − 1)/D > 0. ∎ **Proved.**
**L4 (algebra).** Res_w(E1, E2) = (27/256) φ³ (A′ − B′)³ P, where P is the 167-term polynomial of the table
(`frame`). **Checked symbolically.** For Lean the table is an input: it is in the raw folder, `I136-rpi5/P_table.tsv`.

**L5 (leaf test).** Let h be C¹ on [a, b], m in [a, b], and I an interval with h′(x) in I on [a, b]. Then h([a, b]) ⊂ h(m) + I·([a, b] − m).
If 0 ∉ I, h is strictly monotone, so it has exactly one zero when h(a), h(b) have strictly opposite signs and none when they
have the same strict sign. **Proved** (mean-value theorem).

**L6 (Krawczyk).** Let F: X → ℝ² be C¹ on a box X, m in X, Y a real 2 × 2 matrix, J an interval matrix with F′(x) in J on X,
K = m − Y F(m) + (I − Y J)(X − m). If K ⊂ int X, F has exactly one zero in X. *Proof.* T(x) = x − Y F(x) satisfies
T(X) ⊂ K ⊂ X, so T has a fixed point (Brouwer) where Y F = 0. K ⊂ int X puts the box norm of I − Y J below 1, so Y J̃ is
nonsingular for every J̃ in J. Two zeros x, x′ give 0 = Y(F(x) − F(x′)) = Y J̃ (x − x′), so x = x′. ∎ **Proved**
(the source is Krawczyk, Computing 4 (1969), *(unverified)*).

**L7 (completeness at a root).** Let Z be an interval, W₀ = [−10⁴, 10⁴], and 𝒳 a finite set of boxes Z′ × W′ with Z ⊂ Z′, each with K ⊂ int. If every
box of a finite cover of Z × W₀ either has an enclosure of E1 or of E2 that excludes 0, or lies in a box of 𝒳, then the common
zeros of E1, E2 in Z × W₀ are exactly the |𝒳| zeros of L6. **Proved** (L6 and enclosure). **Computer** (the cover).

**L8 (small z).** Put P̂(φ, ψ, δ) = P(φ, ψ − 2φ, ψ − 2φ + δ) and give φ the weight 2, ψ and δ the weight 4. Then P̂ has 358 terms, every weight is
at least 24, and the part of weight 24 is the monomial 774144 φ⁸δ² (**checked symbolically**, `small_z:*`). With
K₁ = (cosh z − 1)/z², S = sinh z/z, K₂ = (S − 2K₁)/z², K₃ = (S − 1)/z², whose power series in z² have the positive coefficients 1/(2n + 2)!, 1/(2n + 1)!,
2(n + 1)/(2n + 4)!, 1/(2n + 3)! (**proved**, Taylor series of cosh and sinh), one has φ = z²Φ, ψ = z⁴Ψ, δ = z⁴Δ with

    Φ = v K₁/D,    Ψ = v [2K₁² − (1 + C) K₂]/D²,    Δ = −v (1 + C) [−(1 + C) K₃ + (C − 2) K₁ − S K₁ − z² K₁²]/D³,    D = 1 + C + z² K₁

(**checked symbolically**), so P = z²⁴ G with G = Σ cf Φⁱ Ψʲ Δᵏ z^(2i + 4j + 4k − 24), every exponent of z non-negative, and y = z⁴ Δ.
**Computer** (**certified**): G > 0 and Δ < 0 on [0, 1/4] at both walls, from enclosures of K₁, S, K₂, K₃ by the first N terms and the tail 2c_N x_hi^N,
valid when c_(N+1) x_hi ≤ c_N/2 (the ratio of consecutive coefficients decreases).

**L9 (large z).** Every term φⁱ A′ʲ B′ᵏ of P has j + k ≥ 9 and k ≤ 9 (**checked symbolically**, the table). With L2, P = (τz)⁹ S, where S is the
polynomial Σ cf vⁱ (−1)ʲ φ̂ⁱ α̂ʲ β̂ᵏ τ^(j+k−9) s^(9−k) in (s = 1/z, τ, φ̂, α̂, β̂) with non-negative exponents (**checked symbolically**), and y = τ(z β̂ + α̂).
For z ≥ z_t ≥ 2: e ∈ (0, e_t], τ ∈ (0, τ_t] (τ = 2v(1 + C) z e⁻ᶻ decreases for z > 1), s ∈ (0, 1/z_t]. **Proved.** **Computer**
(**certified**): S > 0 on the box, from the enclosures of φ̂, α̂, β̂ over e ∈ [0, e_t] and 64 cells in s; β̂, α̂ > 0.

**L10 (exp).** For |y| ≤ 1/2 and N ≥ 1, |eʸ − Σ_(n<N) yⁿ/n!| ≤ (4/3)|y|ᴺ/N!, since the ratio of consecutive terms of the tail is at most
|y|/(N + 1) ≤ 1/4. **Proved.** exp(x) = (exp(x/2ᵏ))^(2ᵏ), so squaring an enclosure encloses. cosh and sinh are sums of enclosed exponentials.

**Theorem** (the count; **conditional** on the computer parts of L7 to L9 and on he_type3.md's theorem). At each wall, the admissible
points are the common zeros (z, w) of E1, E2 with 0 < w < 1 and T ≠ 0 at the roots of P and of y in (0, ∞); P and y have no root in (0, 1/4] and
[z_t, ∞) (L8, L9) and the roots in [1/4, z_t] are those of the leaves (L5); the common zeros at each root are those of L7.

## 4. Results

Leaves of the cover (`cover:*`, script output), 160-bit grid, one core at nice 10:

| wall | function | leaves | excluded | monotone, no root | one root | roots | seconds |
|---|---|---|---|---|---|---|---|
| Le | P | 1562 | 1554 | 4 | 4 | 4 | 123 |
| Le | y | 46 | 42 | 3 | 1 | 1 | 1 |
| Rodal | P | 2735 | 2723 | 6 | 6 | 6 | 211 |
| Rodal | y | 58 | 53 | 4 | 1 | 1 | 1 |

Every leaf excludes a root or holds exactly one, none is open, and the leaves tile [1/4, 40] and [1/4, 60]. The tanh-form P and y (sympy derivatives of the tanh profile, mpmath,
80 digits) have the sign of the leaf at a point of every leaf and change sign across each root leaf; the 4 + 1 and 6 + 1 refined root enclosures (widths below 6.3 × 10⁻⁴⁶ at Le
and 3.1 × 10⁻³¹ at Rodal for P) contain the mpmath roots of the tanh form.

The ends (`small_z:*`, `tail:*`):

| wall | G on (0, 1/4], lower bound | Δ on (0, 1/4], upper bound | S on [z_t, ∞), lower bound | limit of S, 4v³ |
|---|---|---|---|---|
| Le | 1.0911 × 10⁻⁶⁷ | −3.7453 × 10⁻⁸ | 0.49998278 (z_t = 40) | 0.5 |
| Rodal | 6.5688 × 10⁻¹⁶⁹ | −2.8278 × 10⁻¹⁸ | 3.9999386 (z_t = 60) | 4 |

One cell proves each end at each wall. y > 0 on [z_t, ∞) with β̂ and α̂ above 0.9999999.

The admissible points (`points:*`; boxes of half-widths 10⁻²⁶ in z and 10⁻²⁰ in w; T from A′ + B′ with the box proven away from 0):

| wall | r | statement | |difference| | w | T |
|---|---|---|---|---|---|
| Le | 0.990567889806 | 0.990567889806 | 4.4 × 10⁻¹³ | 0.06005176468 | −8.80794 |
| Le | 1.02393326111 | 1.02393326111 | 1.5 × 10⁻¹² | 0.04098687975 | 7.90635 |
| Le | 1.09466230932 | 1.09466230932 | 3.6 × 10⁻¹² | 0.9128666766 | 22.0305 |
| Rodal | 4.71448517823 | 4.71448517823 | 1.4 × 10⁻¹² | 0.9576393947 | −4.51919 |
| Rodal | 5.07106634302 | 5.07106634302 | 1.8 × 10⁻¹² | 0.2847666782 | 3.7259 |
| Rodal | 5.10338056385 | 5.10338056385 | 3.3 × 10⁻¹² | 0.9770299367 | 4.97171 |

The differences are the rounding of the statement's 12 digits. The other roots of P: Le r = 1.006345187, w = 58.459064; Rodal r = 5.020434866, 5.033740947,
5.076353993 with w = 4.1006549, 8.7172638, 1.030003. y: Le r = 0.9921149667, Rodal r = 4.993740859, no common root with |w| ≤ 10⁴ at either (the leading
coefficients of E1 and E2 vanish there, which is the source of the factor (A′ − B′)³).

## 5. The planted problem and the controls

Planted truth first (`planted`, 1 s). **A**: g = Π (cosh z − cᵢ) with 7 known roots in [1/4, 40], two of them 1.000 × 10⁻⁵ apart, and two outside; the same cover
finds exactly 7 root leaves, one known root in each, and no open leaf. **B**: E1 = (w − u₁)(w − u₂), E2 = (w − v₁)(w − v₂)(w − v₃) with u₁ − v₁ = z − 3, u₂ − v₂ = (z − 7)/2,
u₁ − v₃ = (z − 5)(z − 5 − 10⁻⁵), T = z − 5 − 10⁻⁵: five known common roots (z = 3, 5, 5 + 10⁻⁵, 7, 659/49), exactly three admissible, one with T = 0 (the close pair's second root)
and one with w = 7/2. The same `Cover`, `refine_root`, `krawczyk` and `common_roots` return exactly that.

Controls (`--mutate=NAME`, 13 of them; the targets were frozen in the raw README before the first run, the last two added after the first batch):

| control | exit | checks run | targets fired | other checks that fired | seconds |
|---|---|---|---|---|---|
| `sigma` | 1 | 14 | 2 of 2 | `spurious:le` | 395.68 |
| `drop-P-term` | 1 | 14 | 6 of 6 | `cover:rodal`, `points:rodal` | 332.04 |
| `exp-ulp` | 1 | 14 | 1 of 1 | `planted`, `cover:le`, `cover:rodal`, `points:le`, `points:rodal`, `spurious:le`, `spurious:rodal` | 10.53 |
| `frame-typo` | 1 | 3 | 1 of 1 | `cover:le` | 5.41 |
| `tail-zt` | 1 | 2 | 2 of 2 | - | 4.51 |
| `small-z-range` | 1 | 2 | 2 of 2 | - | 10.92 |
| `drop-leaf` | 1 | 1 | 1 of 1 | - | 135.10 |
| `krawczyk-shift` | 1 | 2 | 2 of 2 | - | 140.31 |
| `w-range` | 1 | 2 | 2 of 2 | - | 367.69 |
| `spurious-on-P-root` | 1 | 2 | 2 of 2 | - | 344.95 |
| `planted-double` | 1 | 1 | 1 of 1 | - | 1.30 |
| `profile-typo` | 1 | 1 | 1 of 1 | - | 1.10 |
| `printed-typo` | 1 | 1 | 1 of 1 | - | 2.90 |

Every control exits 1 and every frozen target fired (the raw README has the frozen table; `controls_table.py` reads the logs). `--mutate=NAME` changes: `sigma` σ + 10⁻³;
`drop-P-term` the term 4φ³B′⁹ of P; `exp-ulp` each exp endpoint one ulp inward; `frame-typo` one sign of the typed frame table; `tail-zt` the tail bound
starts at z = 4; `small-z-range` the small-z range reaches 20 (Le) and 45 (Rodal); `drop-leaf` one leaf of the cover; `krawczyk-shift` the Krawczyk box is
centred 20 half-widths off; `w-range` the admissible window is (0, 100); `spurious-on-P-root` the no-common-root test runs at a root of P;
`planted-double` the planted g gets a double root; `profile-typo` the constant 2 of B′ typed as 3; `printed-typo` the printed E1 typed with (φ² − 2)γ².

Predictions that were wrong (frozen table unchanged, in the raw README):
- `sigma` also fires `spurious:le`: the root of y moves to r = 0.99211696, 7 × 10⁻⁶ from the statement, past that check's 5 × 10⁻⁶. The radii of P's roots move by 1.5 × 10⁻⁷ to 3 × 10⁻⁵.
- `exp-ulp` fires `exp` as frozen. `profile`, `frame`, `small_z` and `tail` stay silent: cosh and sinh are sums of exp enclosures rounded outward again, which absorbs a one-ulp inward move,
  and those checks have wide margins. `planted`, `cover`, `points` and `spurious` raise an empty-interval error. A one-ulp defect of exp is only caught by `exp`.
- `frame-typo` first fired `frame` through a sympy error in the construction of P, not through the printed-E1 comparison. After that, `frame` reports a non-polynomial P itself, and the
  rerun fires through E1 minus printed = 6φA²w(w − 1)(A − B). `profile-typo` and `printed-typo` were added after the first batch because `profile` had no control.
- `drop-P-term` also fires `cover:rodal` and `points:rodal`; `cover:le` stays silent (4 roots, at other radii).

The script is frozen at b2c4ac7: the full run, its output and all 13 controls above ran on it. The first batch (script 3c7df24, 11 controls, one full run) is in the raw folder `batch1/`.
The two controls of the compare script (`compare-drop-term`, `compare-shift-leaf`) each exit 1 and fire their `tables` and `roots` comparison.

## 6. Comparison with #133 (step 7)

`review_he_type3_count_compare.py` opens `he_type3.py` and #133's `cells.tsv` (raw `I133-rpi5`) and its output after the review's own runs; nothing in it feeds back into
`review_he_type3_count.py`. It passes (15 comparisons, 348 s, exit 0).

| item | #133 | this review | result |
|---|---|---|---|
| E1, E2, P tables | from `he_type3.shadow()` (the metric route of #102); `tables` compares them with the printed frame table | built by sympy from the typed frame table | equal term by term: P 167, E1 17, E2 42 terms, ratio 1 |
| arithmetic | mpmath `iv`, 50 digits | `Fraction`, outward rounding to 160 bits, 320 for Newton and Krawczyk | same roots, boxes and counts |
| profile | ε-form (0.1), against the tanh form at 4 radii per wall | cosh form, exact identity with the tanh form, 12 radii per wall | agree |
| function covered | g = P/t⁹ and Yn (y = t Yn/D³) | P and y = B′ − A′ | same roots |
| leaves, P | 961 (Le), 2356 (Rodal) | 1562, 2735 | each of #133's 4 and 6 root leaves holds exactly one of my refined roots; both tile |
| leaves, y | 26, 38 | 46, 58 | its one root leaf per wall holds my one root |
| root enclosures | root leaves 1.5 × 10⁻⁴ to 2.4 × 10⁻³ wide (Le), 3.6 × 10⁻⁶ to 9.1 × 10⁻⁴ (Rodal), shrunk to 10⁻²⁰ | widths below 6.3 × 10⁻⁴⁶ (Le P), 3.0 × 10⁻³¹ (Rodal P), 2.7 × 10⁻⁴⁰, 1.1 × 10⁻⁴⁸ (y) | contained |
| small z | 8 cells, series by interval Horner with hand majorants; G in [1.09451, 1.09693] × 10⁻⁶⁷ on [0, 1/32] (Le), [6.58941, 6.604] × 10⁻¹⁶⁹ (Rodal) | 1 cell on (0, 1/4], series of positive coefficients; G ≥ 1.09109 × 10⁻⁶⁷, 6.5688 × 10⁻¹⁶⁹ | both bracket G(0) = 774144 Φ₀⁸ Δ₀² = 1.09451 × 10⁻⁶⁷, 6.58941 × 10⁻¹⁶⁹; #133's series of Φ, Ψ, Δ lie in my enclosures |
| large z | absolute values of the coefficients, L ≥ 0.168694 (Le), 0.52369 (Rodal) | signs kept, 64 cells in s = 1/z, S ≥ 0.499983, 3.99994 | L ≤ mine ≤ 4v³ = 0.5, 4, the limit of S |
| common roots | 2D bisection with a Krawczyk acceptance on every box | numerical start (mpmath), one Krawczyk box, bisection in w for the rest | same w to the printed digits; same T; radii within 3.6 × 10⁻¹² of its table |
| planted | path (z, z², z³), 9 roots by Sturm, 6 admissible, one carried by y | cosh product with 7 roots, pair 10⁻⁵ apart; P_pl with 5 roots, 3 admissible, T = 0 and w = 7/2 cases | both pass |
| controls | 12, one run each with all checks | 13 | all fire; each side has predictions that missed |

The statements of he_type3_count.md that the code can test hold: P̂ has 358 terms of weights 24, 26, …, 66; the (A′, B′)-degrees of P run from 9 to 15; the series S₂, N₁, N₂ of Lemma 3 give my
Φ, Ψ, Δ.

## What the code proves, and what is assumed

**Certified** (Python 3.14.10 `int` and `Fraction` arithmetic, with the code of `review_he_type3_count.py` named below): L7 to L9 and the counts of section 4, given he_type3.md's theorem.
**Checked symbolically** (sympy 1.14.0): the profile identities, the frame-table algebra, the resultant, the weights, the tail identity. **Numerical evidence**: every mpmath cross-check
(`exp`, `profile`, the tanh-form signs and roots in `cover`, `small_z`, `tail`).

Trusted code: the classes `Iv` (outward-rounded + − × ÷ on `Fraction` endpoints), `Du` (dual numbers), `exp_pt`/`cosh_iv`/`sinh_iv`, the `Poly` evaluator, `Cover`, `refine_root`,
`krawczyk`, `common_roots`, `series_iv`, `small_z_forms`, `tail_bounds`, and sympy's `resultant`, `Poly` and `factor_list` for the tables.

Assumed: the theorem of he_type3.md (E1 = E2 = 0, T ≠ 0 describes the Type III points off the equator and axis), the frame table printed in he_type3.md section 1, constant v, and the
Alcubierre metric of he_type.md. Not covered: |w| > 10⁴ (not needed, only 0 < w < 1 matters), other speeds or wall parameters (#137), and the Lean statements (#135).

## What he_type3_count.md should now say

- Review, last sentence: replace "The interval cover itself was not re-derived by independent interval code." by "#136 re-derived the cover with independent code in exact rational arithmetic
  (`fractions.Fraction`, no mpmath `iv`, its own exp, cosh, sinh enclosures, E1, E2 and P from the printed frame table, cosh form of the profile): the same 4 and 6 roots of P, one root of y per wall,
  the same ends and three admissible points per wall at the same radii to 4 × 10⁻¹², **confirmed** (review_he_type3_count.md)."
- What the code proves, Trusted: add "The cover, the ends and the boxes were reproduced without mpmath `iv` (#136)" and "the table of `he_type3.shadow()` equals, term by term, the tables built from the printed frame table
  (167, 17 and 42 terms, `review_he_type3_count_compare_output.txt`)".
- Section 0: "with 17 and 42 terms in the variables (v, F, a, b, u, c) and more in (φ, A′, B′, w)": in (φ, A′, B′, w) they also have 17 and 42 terms (`frame`, compare `tables`).
- Lemma 4: state that every term of P has k ≤ 9 (the B′-degree is at most 9); the estimate for n ≥ 10 uses it and the script does not assert it.
- he_type3.md Answer and section 5 can keep **certified**; add "reviewed independently in #136" to the Review paragraph.

## Discrepancies

No discrepancy in the count or in a number. Three wording points of he_type3_count.md, none changing a result:
1. Section 0 says E1 and E2 have more than 17 and 42 terms in (φ, A′, B′, w). They have exactly 17 and 42 there as well.
2. Lemma 4 bounds t^(n−9) |Π_n(z)| by z⁹ M_n (tz)^(n−9) for n ≥ 10, which needs k ≤ 9 for every term. That holds (the maximal B′-degree of P is 9, asserted in my `tail_table`), but neither the note nor `tail_P_lemma` states it.
3. The tail constant L (0.1687, 0.5237) is 3 to 8 times below the sign-aware bound (0.49998, 3.99994) and the true limit 4v³; it is correct and not tight.

Not a discrepancy but recorded: my review shares the frame table with #133, so a defect of the printed table itself would be shared; it equals the metric route of `he_type3.shadow()` (compare, `tables`).
