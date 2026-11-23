# A rigorous count of the Type III points of the tanh wall at Le's and Rodal's speeds

November 2026, issue [#133](https://github.com/nmoinvaz/negative-mass/issues/133). Script [`he_type3_count.py`](he_type3_count.py), its stdout
[`he_type3_count_output.txt`](he_type3_count_output.txt). It builds on [he_type3.md](he_type3.md) (#102), which found the six points and enclosed
each in a Krawczyk box, and on the interval cover of [he_tanh_bound.md](he_tanh_bound.md) (#132). Nothing here claims novelty.

Status words are those of AGENTS.md. **Checked symbolically** means an exact polynomial identity in sympy or exact arithmetic in fractions.
**Certified** means a computer-assisted proof whose trusted code is named in the last section. **Numerical evidence** means a computation
with no enclosure. Each lemma proved by hand is stated so that it can be formalized later; none is formalized here.

## Answer

At Le's wall (σ = 8, R_b = 1, so σR_b = 8, v = 1/2) and at Rodal's (σ = 4, R_b = 5, so σR_b = 20, v = 1) the Alcubierre wall

    f(r) = [tanh σ(r + R_b) − tanh σ(r − R_b)] / (2 tanh σR_b)

has **exactly three** Type III points off the equator and axis, at each wall, in r > 0 (**certified**). Each is one circle of revolution about
the motion axis, and θ → π − θ gives its mirror image. The radii are those of he_type3.md section 5, and the cover leaves no other.

| wall | r | cos θ | w = c² | T = b + a u |
|---|---|---|---|---|
| Le, v = 1/2 | 0.990567889806 | 0.245054615703 | 0.0600518 | −8.80794 |
| Le, v = 1/2 | 1.02393326111 | 0.202452166569 | 0.0409869 | 7.90635 |
| Le, v = 1/2 | 1.09466230932 | 0.955440566739 | 0.912867 | 22.0305 |
| Rodal, v = 1 | 4.71448517823 | 0.978590514299 | 0.957639 | −4.51919 |
| Rodal, v = 1 | 5.07106634302 | 0.533635342003 | 0.284767 | 3.7259 |
| Rodal, v = 1 | 5.10338056385 | 0.988448246868 | 0.977030 | 4.97171 |

(The script prints each radius with an enclosure of width below 10⁻²⁵; the table shows 12 digits. T is the value b + a u of he_type3.md.)

The proof has the shape of a count of roots of two functions of r. An admissible point, meaning a common root w ∈ (0, 1) of E1 and E2 with
T ≠ 0, has P = 0 or y = B′ − A′ = 0 (section 1). At each wall P has 4 (Le) or 6 (Rodal) roots in r > 0, three of them admissible, and y has
one root, where E1 and E2 have no common root with w in [−10⁴, 10⁴]. The other roots of P have their common root at w = 58.459 (Le) and at
w = 4.1007, 8.7173 and 1.0300 (Rodal), outside (0, 1). That fills the gap of he_type3.md section 5, where the count of three was **numerical
evidence** from a scan of about 6000 radii.

## 0. Setting

Everything is as in he_type3.md. With φ = vF, A′ = v r f′ and B′ = v r² f″ (F = 1 − f), the polynomials E1 and E2 of that note, at u = 1, are

    E1 = e₀ + e₁ w + e₂ w²,     E2 = f₀ + f₁ w + f₂ w² + f₃ w³,     e_m, f_m ∈ ℚ[φ, A′, B′],   w = c²,

with 17 and 42 terms in the variables (v, F, a, b, u, c) and more in (φ, A′, B′, w). A point of the wall off the equator and axis (0 < w < 1) with
v u ≠ 0 and (a, b) ≠ (0, 0) is Type III exactly when E1 = E2 = 0 and T = b + a u ≠ 0 (he_type3.md, **certified**). On the tanh wall a = f′ < 0
for r > 0, so (a, b) ≠ (0, 0); u = 1/r is finite; v > 0. An **admissible point** is a pair (r, w) with r > 0, 0 < w < 1, E1 = E2 = 0 and T ≠ 0.
T and A′ + B′ differ by the positive factor v r²: A′ + B′ = v r² T.

**Variables.** z = 2σr, ε = e^(−z), C = cosh 2σR_b, D = 1 + 2Cε + ε². In these variables the profile, with the two identities cosh(y + p) cosh(y − p) =
(cosh 2y + cosh 2p)/2 and cosh z = (1 + ε²)/(2ε), is

    f = (1 + C)/(cosh z + C) = 2(1 + C) ε / D,     F = (1 − ε)² / D.

With f′ = 2σ f_z, f″ = 4σ² f_zz and r = z/(2σ), one has A′ = v z f_z and B′ = v z² f_zz. Using ε ∂_ε f = 2(1 + C) ε (1 − ε²)/D² and the second
application of ε ∂_ε,

    φ = v (1 − ε)² / D,     A′ = −t ã,     B′ = z t b̃,     t = 2v(1 + C) z ε > 0,
    ã = (1 − ε²)/D²,     b̃ = [(1 + ε²)² − 2Cε(1 + ε²) − 8ε²] / D³.                                  (0.1)

The script compares (0.1) with the tanh form of he_type3.py (`jets_mp`) at four radii per wall to 10⁻⁴⁰, and its forward derivatives with
`mp.diff` (`identity`). The two walls have C = cosh 16 = 4.443 × 10⁶ and C = cosh 40 = 1.177 × 10¹⁷; v and C are the only inputs, since
(φ, A′, B′) depend on (σ, R_b, r) only through (v, C, z).

## 1. Reduction to roots of two functions

**Lemma 1** (**checked symbolically**, `tables`). The resultant of E1 and E2 in w, formal degrees 2 and 3, is

    Res = (27/256) φ³ (A′ − B′)³ P,      P ∈ ℚ[φ, A′, B′] of 167 terms, total degrees 10, 12, 14, 16, 18, (A′, B′)-degrees 9 to 15.

`tables` checks the identity exactly in sympy against `resultant`, against an exact Sylvester determinant in fractions at five rational points and, as polynomials in one variable, on three random lines, and checks
E1 and E2 against c₂² − 3c₁ and 27c₀ − c₂³ of the frame table printed in he_type3.md section 1, with X = A′², γ = A′, y = B′ − A′.

**Lemma 2** (by hand, from Lemma 1). If (r, w) is admissible, then P(φ, A′, B′) = 0 or B′ − A′ = 0 at r.
*Proof.* A common root of two polynomials in w makes their formal Sylvester determinant vanish. So φ³ (A′ − B′)³ P = 0. Now φ = vF > 0 for r > 0
since f(r) < 1 there, (1 + C)/(cosh z + C) < 1 for z > 0. ∎

So every admissible radius is a root of P or of y = B′ − A′. By (0.1), P = t⁹ g with

    g = Σ cf (−1)^j φ^i ã^j b̃^k z^k t^(j + k − 9),      (the sum is over the terms cf φ^i A′^j B′^k of P, and j + k ≥ 9 for every term),

and y = t Yn / D³ with Yn = z[(1 + ε²)² − 2Cε(1 + ε²) − 8ε²] + (1 − ε²) D. Since t > 0 and D > 0, the roots of P are those of g and the roots of y are those of Yn.
The remaining work is to list the roots of g and Yn in z > 0 and to decide the common roots at each.

## 2. The ends

**Lemma 3** (small z; by hand and **certified**). On 0 < z ≤ 1/4 neither P nor y vanishes.
*Proof.* Put ψ = A′ + 2φ and δ = B′ − A′, and write P̂(φ, ψ, δ) = P(φ, ψ − 2φ, ψ − 2φ + δ). Assign the weight 2 to φ and 4 to ψ and δ. Then P̂ has
358 terms of weights 24, 26, …, 66, and the part of weight 24 is the single monomial 774144 φ⁸ δ² (`tables`, **checked symbolically**).
With the series (cosh z − 1 = z² S₂(z) and so on)

    φ = z² Φ,    ψ = z⁴ Ψ,    δ = z⁴ Δ,
    Φ = v S₂ / D,    Ψ = v N₁ / D²,    Δ = v (1 + C) N₂ / D³,    D = C + cosh z,
    S₂ = Σ_{k ≥ 1} z^(2k−2)/(2k)!,
    N₁ = Σ_{k ≥ 2} n_k z^(2k−4),     n_k = [4^k + 2C − 2 − 2k(1 + C)] / (2k)!,
    N₂ = Σ_{k ≥ 1} m_k z^(2k−2),     m_k = [4^k (2k + 3)/2 − 2kC] / (2k + 1)!,

one has P = z²⁴ G with G = Σ cf Φ^i Ψ^j Δ^k z^(2i + 4j + 4k − 24), every exponent of z non-negative. The coefficients come from
ψ = v[−(1 + C) z sinh z + 2(cosh z − 1) D]/D², whose numerator is (cosh 2z − 1) + 2(C − 1)(cosh z − 1) − (1 + C) z sinh z, and from
δ = v(1 + C) z [sinh z cosh z + C(sinh z − z cosh z) − 2z + z cosh² z]/D³. The script encloses S₂, N₁, N₂ on eight cells of [0, 1/4] by interval Horner
in z² plus a remainder, and D by interval arithmetic. For k ≥ K = 16 the terms are bounded by the majorants (4^k + 2C + 2 + 2k(1 + C))/(2k)!,
(4^k (2k + 3)/2 + 2kC)/(2k + 1)! and 1/(2k)!, whose ratios at consecutive k are at most 5.6/((2k + 2)(2k + 3)) ≤ 10⁻² for k ≥ 16, so each remainder is at most twice its
first omitted term. The enclosure of G excludes 0 on every cell (it is positive, in [1.0945, 1.0969] × 10⁻⁶⁷ on the first cell at Le), and the enclosure of Δ is negative
(its leading coefficient is v(1 + C)(5 − C)/(3D³) < 0 for C > 5). ∎

**Lemma 4** (large z; by hand and **certified**). On z ≥ z_t, with z_t = 40 at Le and 60 at Rodal, P > 0 and y > 0.
*Proof.* Let ε_t = e^(−z_t) and enclose φ, ã, b̃ of (0.1) for ε ∈ [0, ε_t] by intervals. Collect, for each n = j + k and each k ≤ n, the coefficient
c_{n,k} = Σ_i (−1)^(n−k) cf φ^i ã^(n−k) b̃^k (an interval). Then g = Σ_n t^(n−9) Π_n(z), Π_n(z) = Σ_k c_{n,k} z^k. The coefficient of z⁹ in Π₉ is
4 φ³ b̃⁹ ≥ c₉ > 0. For z ≥ z_t ≥ 1, Π₉(z) ≥ z⁹ (c₉ − Σ_{k<9} |c_{9,k}|_max z_t^(k−9)), and for n ≥ 10, t^(n−9) |Π_n| ≤ z⁹ M_n (t z)^(n−9) with
M_n = Σ_k |c_{n,k}|_max. Now t z = 2v(1 + C) z² e^(−z) decreases for z ≥ 2, so t z ≤ s_t = 2v(1 + C) z_t² e^(−z_t). Hence g ≥ z⁹ L with
L = c₉ − Σ_{k<9} |c_{9,k}|_max z_t^(k−9) − Σ_{n ≥ 10} M_n s_t^(n−9), and the script's enclosure has L ≥ 0.1687 (Le) and L ≥ 0.5237 (Rodal).
For y: Yn = (1 − ε⁴) + 2Cε(1 − ε²) + z[(1 + ε²)² − 2Cε(1 + ε²) − 8ε²]. For z ≥ z_t the bracket is at least β = 1 − 2κ_t(1 + ε_t²) − 8ε_t² > 0 with κ_t = Cε_t
(Cε decreases in z), and the other two terms are positive, so Yn ≥ zβ > 0; the script has κ_t ≤ 1.9 × 10⁻¹¹ (Le) and 1.1 × 10⁻⁹ (Rodal), so β > 1 − 3 × 10⁻⁹. ∎

## 3. The cover of [1/4, z_t]

Both g and Yn are evaluated on intervals X with a forward-mode derivative, on mpmath `iv` at 50 digits, from (0.1). A leaf [a, b] is **excluded** if the
enclosure of g on [a, b] excludes 0, or the mean-value form g(m) + g′([a, b]) ([a, b] − m) excludes 0, and it holds **exactly one root** if the enclosure
of g′ excludes 0 and g(a), g(b), each enclosed at its exact end point, exclude 0 with opposite signs. Otherwise it is bisected, to depth 64. This is the
interval Newton and Krawczyk reasoning of #132 in one variable; no unresolved leaf may remain. Leaf end points are exact binary fractions. `cover`
checks that the leaves tile [1/4, z_t] without a gap, and, independently of the interval code, that the tanh-form P (`jets_mp`, 150 digits) has the sign of
g at every leaf end, changes sign exactly across the leaves with one root, and has the same sign as g at the midpoint of every excluded leaf. The leaves
are in the raw folder (`cells.tsv`, `--dump-cells`).

**Lemma 5** (computer-assisted, the cover). g has exactly 4 roots in [1/4, 40] at Le and 6 in [1/4, 60] at Rodal, each simple and isolated in a leaf.
Yn has exactly one root at each wall. The roots are in the table of section 5 (**certified**).
*Proof.* The leaves of `cover` and `spurious`. A leaf with one root contains exactly one by the intermediate value theorem and monotonicity, and the excluded
leaves contain none. ∎

## 4. The common roots at each root of P and of y

Let Z be a leaf with one root of g or of Yn. Interval Newton on g (or Yn), which keeps the root, shrinks Z to width at most 10⁻²⁰. On Z × [−10⁴, 10⁴] the
script bisects the wider side, discards a box when the enclosure of E1 or of E2 excludes 0, and accepts a box when the Krawczyk operator
K(X) = m − Y F(m) + (I − Y J(X))(X − m) with Y ≈ J(m)⁻¹ maps X into its interior. By Krawczyk's theorem a box so accepted holds exactly one zero of
(E1, E2) in (z, w), the interior condition proving the Jacobian is invertible there. The Jacobian is enclosed from the forward derivatives of the jets and
the polynomial in w. Accepted boxes are refined by K(X) ∩ X, and checked to be pairwise disjoint. No box may be left unresolved, so every common root
of E1 and E2 with z in Z and |w| ≤ 10⁴ lies in an accepted box. An accepted box is **admissible** if its w-range is inside (0, 1) and the enclosure of
A′ + B′ over its z-range excludes 0. It is inadmissible by w if its w-range is outside [0, 1]; a box that meets neither case is reported unresolved.

The completeness of the count is then Lemma 2 plus Lemmas 3, 4 and 5 and this search: an admissible point has a radius in a cell of P or of y, and in that
cell (z, w) is in an accepted box. Nothing needs |w| > 10⁴ since only w ∈ (0, 1) matters. The result of section 5 follows.

## 4a. The two standard lemmas the interval code rests on

**Lemma 6** (leaf test; by hand). Let h be C¹ on [a, b], m ∈ [a, b], and I an interval with h′(x) ∈ I for all x ∈ [a, b]. Then h([a, b]) ⊂ h(m) + I · ([a, b] − m) (mean-value theorem).
If 0 ∉ I, h is strictly monotone, so it has exactly one zero in [a, b] when h(a) and h(b) have strictly opposite signs, and none when they have the same strict sign.
Interval arithmetic gives outward-rounded enclosures of h(m), h(a), h(b), I; every enclosure that the script uses is of this kind.

**Lemma 7** (Krawczyk; by hand, the proof of the standard theorem). Let F: X → ℝ² be C¹ on a box X, m ∈ X, Y a real 2 × 2 matrix, J an interval matrix with F′(x) ∈ J for all x ∈ X, and
K = m − Y F(m) + (I − Y J)(X − m). If K ⊂ int X, then F has exactly one zero in X.
*Proof.* T(x) = x − Y F(x) satisfies T(x) = m − Y F(m) + (I − Y F′(ξ))(x − m) for some matrix F′(ξ) with entries between those of F′ on X, so T(X) ⊂ K ⊂ X, and by Brouwer T has a fixed point,
where Y F = 0. The strict inclusion K ⊂ int X forces the interval matrix I − Y J to have norm below 1 in the box norm of X, so every Y J̃ with J̃ from J is nonsingular, hence Y and J̃ are; two zeros x, x′ would give
Y (F(x) − F(x′)) = Y J̃ (x − x′) = 0 with J̃ a matrix from J, so x = x′. ∎ (The source of the theorem is Krawczyk, Computing 4 (1969), *(unverified)*.)

## 5. The roots

| wall | root of | r | z = 2σr | common root w | verdict |
|---|---|---|---|---|---|
| Le | P | 0.990567889806 | 15.8490862369 | 0.0600518 | admissible, T = −8.80794 |
| Le | y | 0.9921149667 | 15.87383947 | none | |
| Le | P | 1.0063451866 | 16.1015229856 | 58.459064 | outside (0, 1) |
| Le | P | 1.02393326111 | 16.3829321777 | 0.0409869 | admissible, T = 7.90635 |
| Le | P | 1.09466230932 | 17.5145969491 | 0.912867 | admissible, T = 22.0305 |
| Rodal | P | 4.71448517823 | 37.7158814259 | 0.957639 | admissible, T = −4.51919 |
| Rodal | P | 5.02043486647 | 40.1634789317 | 4.1006549 | outside (0, 1) |
| Rodal | P | 5.03374094665 | 40.2699275732 | 8.7172638 | outside (0, 1) |
| Rodal | y | 4.993740859 | 39.94992687 | none | |
| Rodal | P | 5.07106634302 | 40.5685307442 | 0.284767 | admissible, T = 3.7259 |
| Rodal | P | 5.07635399285 | 40.6108319428 | 1.030003 | outside (0, 1) |
| Rodal | P | 5.10338056385 | 40.8270445108 | 0.977030 | admissible, T = 4.97171 |

The root of y has no common root of E1 and E2 anywhere in |w| ≤ 10⁴. The coordinator's resultant scan in `I102-mac` (`ind102.py`, C8) found 5 sign
changes at Le and 7 at Rodal, against 4 and 6 of he_type3.md, and wrote that every extra one is a common root with w outside (0, 1). The extra one at each wall is
the root of y, which is a zero of (A′ − B′)³, an odd power, so the resultant changes sign there; it is not a common root.

The nearest approach of an inadmissible root to the admissible range is Rodal's r = 5.07635 with w = 1.0300, a distance 0.03 from the axis condition w = 1 (c² = 1),
resolved by the enclosure.

## 6. The planted problem

A planted problem runs through the same cover, Newton shrink, bisection and Krawczyk code with a known answer. Take the path (φ, A′, B′)(z) = (z, z², z³) and
α(z) = 5/6 − z/3, κ₁ = 2(z − 1/3)(z − 1/3 − 10⁻⁵)(z − 3/2), κ₂ = (z − 5/7)(z − 11/5)(z − 13/3), κ₃ = z³ + z + 1,

    E1 = (w − α)(w − α − 1),      E2 = Π_j (w − α + κ_j),      T = (z − 5/7)(z + 7),

written as polynomials in (φ, A′, B′, w) by substituting z = φ, z² = A′, z³ = B′. The common roots of E1 and E2 are at κ_j = 0 (w = α) and κ_j = −1 (w = α + 1), so
the resultant is Π_j κ_j (1 + κ_j), of degree 18, with 9 positive roots by a Sturm sequence in sympy. The two roots 1/3 and 1/3 + 10⁻⁵ form a close pair. Two roots have w
outside (0, 1), and the root 5/7 has w in (0, 1) with T = 0. To plant the spurious factor, the root 11/5 (of κ₂, with w = 0.1 ∈ (0, 1)) is removed from P_pl, which is the resultant
divided by z − 11/5, and carried by y_pl = (z − 11/5)(z² + 1). So one admissible point is at a root of y_pl and is not a root of P_pl. The exact answer is six admissible points,
from sympy's real roots and exact evaluation of α and T. The covers of P_pl (8 roots) and of y_pl (1 root) over [0, 23], the tail (the smallest N with Σ_{i<n}|c_i| N^(n−i) < |c_n|, which
is a Cauchy-type bound) and the classification return exactly these six, to 10⁻²⁰.

## 7. Checks and controls

The covers have 961 leaves (4 roots of g) at Le and 2356 (6 roots) at Rodal over [1/4, z_t], and 26 and 38 leaves (one root of Yn each). All 10 checks pass.

| check | what it tests | frozen targets of the controls `--mutate=` |
|---|---|---|
| `identity` | (0.1) against the tanh form, derivatives against `mp.diff`, the series of Φ, Ψ, Δ | `eps-typo` |
| `tables` | Lemma 1 exactly, the Sylvester determinant, E1 and E2 against the printed frame table, the weights | `drop-P-term`, `drop-P-term-low` |
| `planted` | section 6 | `no-tail-large`, `T-off`, `w-range`, `krawczyk-nocorr`, `double-root`, `skip-y` |
| `small_z` | Lemma 3 | `drop-P-term`, `drop-P-term-low`, `no-tail-small` |
| `tail_z` | Lemma 4 | `no-tail-large` |
| `cover` | the tiling, no unresolved leaf, the tanh-form signs at the leaf ends and midpoints and across each root cell | `cell-gap`, `eps-typo` |
| `spurious` | the same for Yn, no common root at its root | `eps-typo`, `skip-y` |
| `points` | the admissible points, the radii of he_type3.md, the other roots | `sigma-1e-3`, `drop-P-term`, `drop-P-term-low`, `eps-typo`, `w-range`, `krawczyk-nocorr` |
| `cross` | the 12-digit values of he_type3_output.txt | `sigma-1e-3`, `drop-P-term`, `drop-P-term-low`, `eps-typo` |
| `assembly` | the three ranges tile (0, ∞) | `drop-P-term`, `drop-P-term-low`, `no-tail-large`, `no-tail-small`, `cell-gap`, `eps-typo`, `skip-y` |

The 12 controls (`CONTROLS` in the script) were written to the raw README before the first run, except `drop-P-term-low`, added after the first pass of step 4. With all checks run per control,
every control exits 1 and nine fail exactly their frozen targets. Three do not:

- `drop-P-term` (the first term of the table) leaves `small_z` and `assembly` silent. The first term has total degree 18, so weight at least 36 in (φ, ψ, δ), and cannot touch the part of weight 24.
  The prediction was wrong. `drop-P-term-low` (the last term, total degree 10, weight 20) fails `small_z` and `assembly` as predicted.
- `w-range` also fails `cross`, which needs three admissible points.
- `krawczyk-nocorr` leaves `points` silent: at each real root cell there is exactly one common root, and a box accepted without the correction term still contains it. The planted cells with
  two candidate roots expose the unsound box.

The first pass of step 4 (raw `step4_pass1/`) also had `spurious` silent under `eps-typo`; the checks `cover` and `spurious` were then strengthened to require the tanh-form function to change sign across every
shrunk root cell. The final run reports a control as firing when at least one of its frozen targets fails, and prints the silent targets.

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/he_type3_count.py

It exits 1 if a check fails or no target of a control fails. `--no-controls` runs the checks alone, `--mutate=NAME` one control with all checks, `--checks=a,b --walls=Le` a subset,
`--dump-cells=F` writes every leaf of both covers. The final run runs each control as a child process on its target checks only (`--only-targets`); the run with all checks per control is
in the raw folder (`mutate_*.log`). On rpi5 (one core, nice 10, shared with choptuik-large-d) the checks take 351 s, each control with all checks about 200 s, and the full run 1714 s with a
peak rss of 124 MiB.

## What the code proves, and what is assumed

**Certified** (sympy 1.14.0 exact, fractions, and mpmath 1.3.0 `iv`): Lemma 1 and the weights (exact); Lemmas 3, 4, 5 and the Krawczyk boxes (interval arithmetic with the
remainder bounds above); the counts of section 5 given he_type3.md's theorem. **Numerical evidence**: the cross-checks against the tanh form (150 digits) and the 12-digit
values of he_type3_output.txt, which agree but prove nothing.

Trusted: mpmath's `iv` (outward-rounded + − × ÷ and `exp`), sympy (`resultant`, `Poly`, `ring`, `real_roots`), Python's `Fraction`, and the table of E1, E2 and P of
`he_type3.shadow()`. Of that table, `tables` re-derives the identity of P with the resultant, an independent determinant at five points, and E1, E2 from the printed frame
table of he_type3.md; it does not re-derive the frame table itself, which is #101's and #102's (checked there against the metric by independent code, `I102-mac`).

Assumed: the theorem of he_type3.md (the description by E1, E2, T), constant v, and the Alcubierre metric of he_type.md. Not covered: the Type III set of the accelerating bubble
(#73), other speeds or wall parameters, points with |w| > 10⁴ that are not needed, and the Lean statement of the lemmas of sections 1 to 4.

## Review

The coordinator reviewed the count (raw `I133-mac`, Mac, Python 3.14, sympy 1.14.0, mpmath 1.3.0):

- `he_type3_count.py` reran, exit 0 in 480 s, with output identical except for timings. `--mutate=eps-typo` fails exactly its six targets.
- Rodal's root at r = 5.07635 was re-solved by Newton from 10⁻³ away. It has w = 1.03000296, outside (0, 1).
- `ind133b.py` shares no code with this note's script or with he_type3.py. It derives E1 and E2 from the metric by the Ricci
  route of #102's review, and evaluates the profile in the cosh form of #132.
  - D1 and D6 are exact. D1: the resultant is (27/256) φ³ (A′ − B′)³ P, with 167 terms, the weight-24 part 774144 φ⁸ δ², and the φ³B′⁹
    coefficient 4. D6: at y = 0 both leading coefficients vanish, which is where the factor (A′ − B′)³ comes from.
  - The rest is numerical evidence. A 10⁻⁴ grid in z finds P's sign changes and, without the resultant, the common roots on E1's real
    branches: 4 and 6, at this note's radii and w, with three admissible per wall. P > 0 at small and large z, where the leading terms
    match those of Lemmas 3 and 4.
  - Five mutations each make it fail.
- The first version, `ind133.py`, failed its grid criterion through two defects of its own. Its float64 quadratic formula lost E1's small
  root, and it took the root at w = ∞ for a common root. That failure is recorded in the raw folder.
- `cmp133.py` cross-loads the table of P from `he_type3.shadow()` that this note trusts. It equals, term by term, the P of the
  independent Ricci route.

The interval cover itself was not re-derived by independent interval code.

## Discrepancies

One with the wording of the review of #102: the extra sign change of the resultant scan is the root of the spurious factor (A′ − B′)³ and not a common root with w outside (0, 1)
(section 5). No discrepancy in the count: three at each wall, at the radii of he_type3.md.
