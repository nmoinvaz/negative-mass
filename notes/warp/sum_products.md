# Sums of two products u = X₁(x)R₁(ρ) + X₂(x)R₂(ρ) + v with det R₃ = 0

November 2026, issue [#169](https://github.com/nmoinvaz/negative-mass/issues/169). Script
[`sum_products.py`](sum_products.py), its stdout [`sum_products_output.txt`](sum_products_output.txt). It follows
[`mult_classify.md`](mult_classify.md) (#160) and [`det_classify.md`](det_classify.md) §3–§5 (#80). Imported from
`det_classify.py`: `E_mult_poly` and `E_add_poly`, only as comparison targets. Shared with it: the formula for D
(retyped). Nothing is imported from `mult_classify.py`. Status words as in AGENTS.md: **proved** (by hand), **checked
symbolically**, **numerical evidence**, **not decided**.

## Answer

**Not decided.** The question is whether u = X₁R₁ + X₂R₂ + v solves D = det R₃ = 0 on an open box under a nondegeneracy
that keeps it from one product (#160) or the additive case (#80). The rank route of #160 does not reach it: the exact
separation ranks of E = 8D are 59 (v = 0), 93 (v = 1/3), 87 (v = 1) and 102 (symbolic v), against 6 and 10 in #160, and
a pair of exponentials on each side satisfies d_F + d_G ≤ r with room although it is not a solution (§3). What is
decided:

| case | verdict | status |
|---|---|---|
| X₁, X₂ dependent or R₁, R₂ dependent (one product) | no solution with X X′ R R′ ≠ 0, every v | **proved** in #160 |
| W_x W_ρ ≠ 0, 1 ∈ V_X and 1 ∈ V_R (α ≡ 0 and δ ≡ 0) | u = a + bX + cR + dXR: one product plus a constant (d ≠ 0, #160) or additive (d = 0, #80 Theorem A) | **proved** (§2), algebra **checked symbolically** |
| W_x W_ρ ≠ 0 and the ρ-side has constant ODE coefficients (r″ = γr′ + δr, γ, δ constants) | no solution, every v, with or without (N3) | **proved** (§5), identities **checked symbolically** |
| W_xW_ρ ≠ 0 and the x-side is X₁ = e^{k₁x}, X₂ = e^{k₂x}, k₁ ≠ k₂ real (one may be 0) | no solution, every v | **proved** (§5b: top frequency, then eight unit ideals); identities and unit ideals **checked symbolically** |
| W_xW_ρ ≠ 0 otherwise: (b) α ≡ 0, δ ≠ 0; (c) α ≠ 0, δ ≡ 0; (d) both ≠ 0; in particular the x-side with constant coefficients and complex exponents (X = e^{ax}cos bx, e^{ax}sin bx) or a double root (e^{kx}, xe^{kx}) | **not decided** | — |
| planted additive family u = c x + g(ρ), c ≥ 3/10 | found through the same E; perturbed X rejected by the rank bound, perturbed R by an exact test | **checked symbolically**, **numerical evidence** at 130 digits (§4) |

So the issue's second criterion (a proof under the nondegeneracy, or an explicit solution that passes closure) is **not
met**. No explicit solution was found other than the additive family, which fails closure on every member
([`bessel_closure.md`](bessel_closure.md), #162). The cases left open are (b), (c) and (d) outside the strata of
Theorems K and E; §6 says what was tried. Theorem E (the real exponential x-side) was added after step 2 at the
coordinator's request, README A2 and A3.

## 1. E, and its compact form

With the jets A_k = X₁^(k), B_k = X₂^(k) (k = 0, 1, 2), a₀, a₁, a₂ = R₁, R₁′, S₁ = R₁″ + R₁′/ρ and b₀, b₁, b₂ likewise
for R₂,

    u = A₀a₀ + B₀b₀ + v,  u_x = A₁a₀ + B₁b₀,  u_xx = A₂a₀ + B₂b₀,  u_ρ = A₀a₁ + B₀b₁,  u_xρ = A₁a₁ + B₁b₁,  u_ρρ + u_ρ/ρ = A₀a₂ + B₀b₂,

and E = 8D is a polynomial in these 12 jets and v with 376 terms (`e_two`, equal to 8D from derivatives of undefined X₁,
X₂, R₁, R₂, **checked symbolically**). It is irreducible over ℚ (`factor_list`: one factor, multiplicity 1, 68 s;
optional flag `--irreducible`, raw `dev/irreducible.log`), so there is no factor to divide out, unlike #160's R′². With
X₂ = 0 it is R₁′² times `det_classify.E_mult_poly`, and with X₂ = 1, R₁ = 1 it is 4R₂′² times `E_add_poly(X + v, …)`.

**Compact form** (**proved** by expansion, **checked symbolically** with abstract P, Q, Z, Y, h): with

    A = 2hP + Y² − Z²,  B = 2hQ − 2YZ,  C = Y² + Z² + 2h²,   4hD = C² − A² − B².

In complex form, ξ = P + iQ and η = Y + iZ give A + iB = 2hξ + η̄² and C = |η|² + 2h², so on {h ≠ 0} D = 0 is |2hξ +
η̄²| = |η|² + 2h². Nothing below uses this beyond the identity; it is the two-product analogue of #160's s² − δ² − τ²
form.

## 2. The ODE normal form, the nondegeneracy, the excluded forms

Notation: a = (X₁, X₂), r = (R₁, R₂), s = r″ + r′/ρ, p∧q = p₁q₂ − p₂q₁, V_X = span{X₁, X₂}, V_R = span{R₁, R₂}, W_x =
a∧a′, W_ρ = r∧r′.

**Lemma 1 (normal form, proved).** If W_x ≠ 0 on I and W_ρ ≠ 0 on J then a″ = αa + βa′ and s = Γr′ + δr, with α =
(a″∧a′)/W_x, β = (a∧a″)/W_x, Γ = (r∧s)/W_ρ, δ = (s∧r′)/W_ρ continuous (Cramer's rule in ℝ²). With n₀ = u − v, n₁ = u_x,
n₂ = u_ρ, n₃ = u_xρ and N = [[n₀, n₂], [n₁, n₃]]: u_xx = αn₀ + βn₁, u_ρρ + u_ρ/ρ = Γn₂ + δn₀, det N = W_xW_ρ, and

    E = n₂⁶ − 4n₂²P² − n₂²Q₂² + n₂²Y₂² + n₂²Z₂² − 2P(Y₂² − Z₂²) + 2Q₂Y₂Z₂ =: D̃(N; α, β, Γ, δ, v),
    P = −(n₁² + (n₀ + v)(αn₀ + βn₁)),  Q₂ = −(Γn₂ + δn₀),  Z₂ = n₃,  Y₂ = 2n₁n₂ + (n₀ + v)n₃

(`ode_form`, **checked symbolically**). Γ = γ + 1/ρ with r″ = γr′ + δr.

**Lemma 2 (constants, proved).** Let W_x ≠ 0. Then X₁′X₂″ − X₂′X₁″ = −αW_x (identity, `ode_form`), and 1 ∈ V_X ⟺ α ≡ 0.
*Proof.* (⇒) If pX₁ + qX₂ = 1 then 0 = pX₁″ + qX₂″ = α(pX₁ + qX₂) + β(pX₁′ + qX₂′) = α. (⇐) If α ≡ 0 then X_i″ = βX_i′,
so X_i′ = c_i e^{∫β} with constants c_i, (c₁, c₂) ≠ 0 since W_x ≠ 0. Then c₁X₂ − c₂X₁ has zero derivative, equals a
constant c₀, and c₀ ≠ 0 (else a is dependent). So 1 = (c₁X₂ − c₂X₁)/c₀ ∈ V_X. ∎ The same for δ and V_R.

**Lemma 3 (shrinking, proved).** Every solution on a box restricts to a subbox of one of three kinds: (i) W_x ≡ 0 or W_ρ
≡ 0, so u − v is one product (if W_x ≡ 0 and X₁(x₀) ≠ 0 then (X₂/X₁)′ = −W_x/X₁² = 0 near x₀, so X₂ = κX₁ there; if X₁ ≡
X₂ ≡ 0 the shift is planar); (ii) W_xW_ρ ≠ 0 with α ≡ 0 and δ ≡ 0; (iii) W_xW_ρ ≠ 0 with α ≠ 0 everywhere or δ ≠ 0
everywhere (α, δ are continuous, so on a subbox each is ≡ 0 or never 0). Kind (iii) is **(N)** of the README: W_x ≠ 0,
W_ρ ≠ 0, and 1 ∉ V_X or 1 ∉ V_R.

**Lemma 4 (the excluded forms are the decided ones, proved).** In kind (ii), V_X = span{1, X} and V_R = span{1, R} for
some functions X, R, so u = a + bX + cR + dXR with constants (the v is in a). For d ≠ 0 this is d(X + c/d)(R + b/d) + (a
− bc/d), a single product with speed v′ = a − bc/d: #160 gives no solution with X X′ R R′ ≠ 0 for every real v′. For d =
0 it is additive: #80 Theorem A gives u = c x + g(ρ) + const with g from the Bessel equation, or a planar shift
(`case_a`, **checked symbolically**). Kind (i) is #160. ∎

**What (N) leaves out:** complex X, R; functions that are not C²; boxes touching ρ = 0. Everything else is kind (i) or
(ii).

## 3. Separation, Lemma L, and why the rank bound does not decide

E = Σ_{i ≤ r} F_i(x)G_i(ρ) with x-monomials in the six x-jets and ρ-monomials in the six ρ-jets (v on the ρ-side when
symbolic). Lemma L (#160 §2, **proved**): if this vanishes on a box then span{F(x)} ⟂ span{G(ρ)}, so d_F + d_G ≤ r. The
exact ranks (`sep_rank`, **checked symbolically**; the matrix is x-monomials × ρ-monomials):

| form | v = 0 | v = 1/3 | v = 1 | v symbolic |
|---|---|---|---|---|
| E, two products | 59 (107 × 72) | 93 (180 × 106) | 87 (169 × 100) | 102 (180 × 126) |
| X₂ = 0 (one product; #160's 6, 10, 9) | 6 | 10 | 9 | 11 |
| X₂ = 1, R₁ = 1 (additive; #160's 5) | 5 | 5 | 5 | 5 |
| (b) X₂ = 1, u = Xβ(ρ) + γ(ρ) | 21 | 21 | 21 | 21 |
| (c) R₂ = 1, u = A(x)R(ρ) + B(x) | 12 | 12 | 12 | 24 |

In #160 the bound d_F + d_G ≤ 6 or 10 met lower bounds 4 + 3 and 5 + 5 almost at once. Here it does not: for (X₁, X₂) =
(eˣ, e^{κx}) and (R₁, R₂) = (e^ρ, e^{μρ}) with κ, μ transcendental the spans have exact dimensions (`slack`, over ℚ(κ),
ℚ(μ), **checked symbolically**)

| v | d_F | d_G | d_F + d_G | r |
|---|---|---|---|---|
| 0 | 12 | 22 | 34 | 59 |
| 1/3 | 22 | 36 | 58 | 93 |

and E ≠ 0 on the pair (its ι² coefficient is −u_ρ⁴, §5). So the bound is not violated by a pair that is not a solution,
and every stratum has to be decided by an argument that is not a count of dimensions. (The numerical ranks of the 12
coefficient functions of case (c) on sample curves are 4 to 9 for exponential, polynomial and trigonometric pairs and 12
for a generic pair (`dev/probe6.py`), the same picture.)

## 4. Planted truth through the same E

The additive family u = cx + g(ρ), g = −ε ln(αJ₀(√2cρ) + βY₀(√2cρ)) of `det_classify.md` §3, with c ≥ 3/10, found
through E of §1.

- `add_decomp`: E at X₂ = 1, R₁ = 1, v = 0 is Σ_{k ≤ 5} f_k(x)g_k(ρ) from the pivot columns of its coefficient matrix,
  rank 5, with f = (4A₁⁴ − 4A₀²A₂², −8A₀A₂², 4A₁², −4A₂², 1).
- `add_stratum`: the x-curves with all f_k constant (d_F = 1): the elimination of {f_k − k_k} with A₀ free (conditions
  on k from the A₀-coefficients, then a second Gröbner basis) contains A₂² (X″ = 0), and for every constant c the jets
  (A₀ free, A₁ = c, A₂ = 0) have f independent of A₀. So d_F = 1 means X = cx + d.
- `add_planted_exact`: for R = −ε ln ψ with ψ″ + ψ′/ρ + 2c²ψ = 0, E_add(X = cx) = (R′²/2 + c²)² − S²/4 = 0 identically
  (**checked symbolically**, ε = ±1).
- `add_planted`: c ∈ {3/10, 1/2, 1, 3/2, 2}, ε = ±1, R = −ε ln(J₀(z) + 0.3 Y₀(z)): the largest relative residual of f·g
  is 1e-130 at 130 digits (**numerical evidence**), d_G = 4, d_F + d_G = 5 = r, so the bound is saturated.
- `add_perturbed_exact`: with the Bessel equation for R and δ ∈ {ρ²/10, ρ³/10}, E_add(X = cx, R + δ) is a nonzero
  polynomial and is nonzero at three exact points (an exact test; #160's thresholded version failed at one point and is
  not used). `add_xperturbed`: X = cx + x²/10 with the same R has exact d_F = 3 and d_G = 4, so d_F + d_G = 7 > 5 and
  Lemma L rejects it.

## 5. Theorem K: the ρ-side with constant ODE coefficients

Write r″ = γr′ + δr, so s = Γr′ + δr with Γ = γ + ι, ι = 1/ρ. In Lemma 1 the only place ρ enters D̃ explicitly is Γ, in
Q₂ = −((γ + ι)n₂ + δn₀). Hence E is a polynomial of degree 2 in ι, and

    the coefficient of ι² in E is −n₂²·n₂² = −n₂⁴ = −u_ρ⁴   (D itself: −h²/2).

**Theorem K (proved).** Let W_x ≠ 0 on I and W_ρ ≠ 0 on J, and let γ, δ be real constants. Then D ≡ 0 on I × J is
impossible, for every real v and without (N3).

*Proof.* Let μ₁, μ₂ be the roots of μ² = γμ + δ. *Distinct roots.* The solution space of r″ = γr′ + δr over ℂ is spanned
by e^{μ₁ρ}, e^{μ₂ρ}, and R₁, R₂ span it by W_ρ ≠ 0, so r = T(e^{μ₁ρ}, e^{μ₂ρ}) for some invertible T ∈ GL₂(ℂ). Put (X̃₁,
X̃₂) = aT. E is a polynomial in the pairings of the x-jets with the ρ-jets, which are unchanged by (a, r) ↦ (aT, rT⁻ᵀ),
so E = 0 holds in the new jets, where u − v = X̃₁w₁ + X̃₂w₂, w_i = e^{μ_iρ}. With S_i = (μ_i² + ιμ_i)w_i, `rho_const`
shows E = Σ_{k ≤ 2} ι^k F_k with F₂ = −(X̃₁μ₁w₁ + X̃₂μ₂w₂)⁴ and F_k polynomials in w₁, w₂ and the x-jets. For fixed x, E
is a finite sum of terms c_{ω,k}(x) ρ^{−k} e^{ωρ}, and the functions ρ^m e^{ωρ} for distinct (m, ω) ∈ ℤ × ℂ are linearly
independent on an interval, so every c_{ω,k}(x) = 0. The coefficient of ρ⁻²e^{4μ₁ρ} comes only from F₂ and only from the
monomial w₁⁴ (a monomial w₁^i w₂^{4−i} has frequency 4μ₂ + i(μ₁ − μ₂), equal to 4μ₁ only for i = 4): it is −(X̃₁μ₁)⁴.
Likewise −(X̃₂μ₂)⁴ at e^{4μ₂ρ}ρ⁻². So X̃₁μ₁ = X̃₂μ₂ = 0 on I. The roots are distinct, so one of μ₁, μ₂ is nonzero, and
the corresponding X̃_i ≡ 0; then aT is dependent and W_x det T = 0, contradicting W_x ≠ 0. *Double root.* r = T(e^{μρ},
ρe^{μρ}), and with w = e^{μρ} and ι = 1/ρ (`rho_const`) the ρ⁻² coefficient of E, as a Laurent polynomial in ρ times
powers of w, is −w⁴(μX̃₁ + X̃₂)⁴, and no power below ρ⁻² occurs. Frequencies are kμ, and for μ ≠ 0 the frequency 4μ
arises only from w⁴; for μ = 0 the whole E is a Laurent polynomial in ρ. Either way (μX̃₁ + X̃₂)⁴ = 0 on I, so X̃₂ =
−μX̃₁, aT is dependent, a contradiction. ∎

Real or complex roots, any v, and the stratum is closed under no further assumption on the x-side. It excludes, for
example, a ρ-profile of exponentials, trigonometric functions, polynomials of degree ≤ 1 (r = (1, ρ), μ = 0) and their
products with exponentials. It does not touch the Bessel pair (Γ ≡ 0, δ constant: the ODE is r″ + r′/ρ = δr, so γ = −1/ρ
is not constant), nor the additive family (γ = g″/g′ is not constant).

## 5b. Theorem E: the real exponential x-side (added after step 2)

**Theorem E (proved).** Let X₁ = e^{k₁x}, X₂ = e^{k₂x} with k₁ ≠ k₂ real, labeled so that |k₁| ≥ |k₂| (so k₁ ≠ 0), and
let W_ρ ≠ 0 on J. Then D ≡ 0 on I × J is impossible, for every real v (for k₂ = 0 with (N3)).

*Proof.* Put z_i = e^{k_ix}. E is a polynomial in z₁, z₂ with every monomial z₁^iz₂^j of i + j ≤ 6 (`x_exp`), so for
fixed ρ it is a finite sum Σ c_{ij}(ρ) e^{(ik₁ + jk₂)x}, and the coefficient of each frequency vanishes (distinct real
exponentials are independent). (1) *Top frequency.* The frequency 6k₁ comes from (i, j) = (6, 0). For k₂ ≠ 0 another
monomial of the support shares it iff (6 − i)q = j with q = k₁/k₂, |q| ≥ 1, q ≠ 1; for i < 6 and j ≤ 6 − i this gives q
= j/(6 − i) ∈ [0, 1], so q = 1, excluded (i = 6 forces j = 0); for k₂ = 0 only i = 6 shares it, and j = 0 there. So
c₆₀(ρ) = 0, where, from E at X₂ = 0 (**checked symbolically**), c₆₀ = a₁²(a₁² + 4k₁²a₀²)(a₁² + 5k₁²a₀²), a sum of
nonnegative terms, zero only if a₁ = 0. So R₁′ ≡ 0 on J, R₁ = c₁, and c₁ ≠ 0 because W_ρ ≠ 0. (2) If k₂ = 0, then 1 ∈
V_X and 1 ∈ V_R: the excluded form (Lemma 4), here u = c₁e^{k₁x} + R₂(ρ) + v with X″ ≠ 0, which #80 Theorem A excludes.
(3) If k₂ ≠ 0: W_ρ = c₁R₂′ ≠ 0, so b₁ = R₂′ ≠ 0 on all of J. With a₀ = c₁, a₁ = a₂ = 0, E = z₂²b₁²E″ for a polynomial E″
(14 monomials; E = 0 at b₁ = 0, the planar family u = c₁e^{k₁x} + c₂e^{k₂x} + v, so that planted solutions are in the
variety of E and are removed only by dividing out b₁²). At each ρ the values (b₀, b₁, b₂) satisfy all frequency
conditions. For q = k₁/k₂ outside the finite set of ratios that make two monomials of E″ equal in frequency (iq + j =
i′q + j′), the 14 coefficients vanish separately, and together with t·c₁k₁k₂(k₁ − k₂)b₁ − 1 they generate the unit ideal
in ℚ[t, b₀, b₁, b₂, c₁, k₁, k₂, v] (Gröbner basis [1], sympy grevlex, 3 s). The resonant ratios with |q| ≥ 1, q ≠ 1 are
exactly −2, −1, 4/3, 3/2, 2, 3, 4 (printed by the script); for each, the coefficients of equal frequency iq + j are
added (k₁ = qk₂) and, with t·c₁k₂b₁ − 1, generate the unit ideal in ℚ[t, b₀, b₁, b₂, c₁, k₂, v] (0.1 s each). A unit
ideal has no complex zero, so no values c₁ ≠ 0, k₂ ≠ 0, b₁ ≠ 0, b₀, b₂, v exist: contradiction. ∎

The unit ideals are Gröbner bases over ℚ computed in sympy, with no Singular and no certificate to verify separately.
What stays open on the x-side with constant coefficients: complex exponents (X = e^{ax}cos bx, e^{ax}sin bx; in the
basis e^{(a ± ib)x} the coefficient c₆₀ can vanish with R̃₁′ = ±2ik₁R̃₁ or ±i√5k₁R̃₁, an exponential in ρ, so step (1)
does not close) and the double root (e^{kx}, xe^{kx}).

## 6. What was tried, what failed, what is open

**Open:** (b), (c), (d) outside Theorem K's stratum. A count shows why a decision needs a method of its own: in the jets
(N; α, β, γ, δ) of Lemma 1, the equations ∂ₓ^i∂_ρ^j E = 0 (i ≤ p, j ≤ q) have (p + 1)(q + 1) equations and 2(p + q) + 10
unknowns (N, ι, v, and α_k, β_k, γ_k, δ_k up to the orders), so the first overdetermined square truncation is (5, 5): 36
equations against 30 unknowns, with degrees growing to about 12. This is a count, not a result.

**Failed routes, kept in raw `I169-rpi5/dev`.**
- *Singular, prolongation ideal.* `std` over ℤ/32003 of E and its derivatives to order (2, 2) in the jets of Lemma 1,
  saturated by det N·α, did not finish in 5 min (`gen_prol.py`); case (c) (u = A(x)R(ρ) + B(x), jets A_k, B_k, R_k, ι,
  saturated by A₁B₂ − A₂B₁ and R′) did not finish at (2, 2) or (3, 3) in 2 min each (`gen_prol_c.py`). No certificate
  exists and none is claimed. #160's prolongation did not finish in sympy either; Singular did not change that at these
  orders.
- *Numerical formal series* (`formal2.py`, `dev1.log`): Taylor data for α, β, γ, δ, N₀ against the coefficients of
  s^iρ^j (i ≤ p, j ≤ q), Levenberg–Marquardt, det N₀ fixed. At (3, 3) it finds nondegenerate formal solutions (cost
  1e-28; underdetermined). At (4, 4) and (5, 5) the best costs are 3e-8 to 1e-5 in case (d) **and** in the planted
  additive case (a) (and the same minima appear for the same seeds), so the planted control is silent and the search
  says nothing about (d). It is not used.
- *Generic strata.* If the 72 ρ-monomials (v = 0) of E were linearly independent functions of ρ, each coefficient
  polynomial in the x-jets would vanish. The pure monomials a₁^i b₁^{6−i} come only from h³, with coefficient binom(6,i)
  A₀^iB₀^{6−i}/8, which forces A₀ = B₀ = 0 and so W_x = 0. This is a rank statement in the style of #160 and is as weak
  as it looks: the strata with relations among the ρ-monomials are all that matter.

**Natural next steps (recommended tasks).** (i) The remaining x-sides with constant coefficients: complex exponents and
the double root (§5b); the lexicographic top frequency of Theorem E generalizes with a case analysis of the branches
R̃₁′ = ±2ik₁R̃₁, ±i√5k₁R̃₁. (ii) Case (c), where #160's d_G ≥ 5 argument applies to the 12 ρ-functions (R′⁴, R²R′², R⁴,
RR′², R³, S², R′², RS, R², S, R, 1) and the x-side has two functions: not done. (iii) A fast differential-elimination
engine on the prolongation, with certificates that sympy checks. (iv) The strata (b), (c), (d) named in the table of the
Answer.

## 7. Checks and controls

13 checks, all pass; 15 controls fire, each on its target check and on no other, and two more (`s-rank-drop`,
`s-xexp-sat`) are declared silent (README A1, A3). Per check, in seconds on a Raspberry Pi 5 (one process at nice 10):
`e_two` 1.7, `ode_form` 1.8, `case_a` 0.0, `sep_rank` 2.6, `slack` 2.4, `rho_const` 4.2, `x_exp` 4.9, `add_decomp` 0.2,
`add_stratum` 0.1, `add_planted_exact` 0.2, `add_planted` 34.7, `add_perturbed_exact` 0.2, `add_xperturbed` 23.8: about
77 s in all. The whole run with the 17 control runs took 23 min 25 s and 81 MB peak resident (`/usr/bin/time -v`), exit
0. The first run of record (12 checks, before `x_exp` was added, raw `logs/run_01.log`) took 18 min 10 s.

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `e_two` | E from derivatives, 4hD = C² − A² − B², X₂ = 0 and additive specializations | `s-d-sign`, `s-compact-b` | Y with −uu_xρ/2; B = 2hQ + 2YZ |
| `ode_form` | Cramer form, D̃(N; α, β, Γ, δ), det N, X₁′X₂″ − X₂′X₁″ = −αW_x | `s-ode-alpha` | α with a′∧a″ |
| `case_a` | u = a + bX + cR + dXR as one product plus constant, additive, α = 0 for X₂ = 1 | `s-casea-const` | constant a in place of a − bc/d |
| `sep_rank` | exact ranks of the table of §3 | `s-rank-qyz` (and `s-rank-drop`, silent) | 2QYZ removed (ranks 56, 84, 78, 93 at the full form); hZ² removed (column counts only) |
| `slack` | d_F, d_G of the exponential pair, E ≠ 0 | `s-slack-kappa` | κ = 2 (resonant) |
| `rho_const` | ι² coefficient −u_ρ⁴, w^4 coefficients, double root ρ⁻² coefficient | `s-rhoconst-iota`, `s-rhoconst-double` | S_i without R_i′/ρ; R₂′ = μρe^{μρ} |
| `x_exp` | z₁⁶ coefficient, top frequency, E = z₂²b₁²E″, planar family, unit ideals (generic and 7 resonant ratios) | `s-xexp-coef`, `s-xexp-planar` (and `s-xexp-sat`, silent) | 6k₁²a₀² in the second factor; coefficients of E and no saturation by b₁ |
| `add_decomp` | rank 5 factorization of E (additive) | none | |
| `add_stratum` | d_F = 1 ⇒ X″ = 0 | `s-stratum-drop` | three of the five f_k only |
| `add_planted_exact` | Bessel ⇒ E_add = 0 | `s-add-exact` | wrong sign of 2c²ψ |
| `add_planted` | numerical Bessel family, d_F + d_G = 5 | `s-add-freq` | J₀(cρ) in place of J₀(√2cρ) |
| `add_perturbed_exact` | E_add(R + δ) ≠ 0 exactly | `s-pexact-zero` | δ = 0 |
| `add_xperturbed` | d_F + d_G > 5 for X = cx + x²/10 | `s-xperturb-zero` | X = cx |

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/sum_products.py --no-controls     # the 13 checks, about 80 s
    python3 notes/warp/sum_products.py --mutate=NAME     # one control
    python3 notes/warp/sum_products.py                   # the checks, then every control (about 19 min)
    python3 notes/warp/sum_products.py --only=irreducible --irreducible   # factor_list of E, 68 s

## 8. What the code proves, and what is assumed

Exact identities (sympy): E from derivatives of undefined functions; 4hD = C² − A² − B²; the specializations to
`E_mult_poly` and `E_add_poly`; the Cramer identities and D̃ of Lemma 1; X₁′X₂″ − X₂′X₁″ = −αW_x; the algebra of
`case_a`; the exact ranks of `sep_rank`; the exact dimensions of `slack` over ℚ(κ), ℚ(μ); the ι² coefficient and the
extracted coefficients of `rho_const`; the factorization f·g of `add_decomp`, the Gröbner membership of `add_stratum`;
`add_planted_exact`; the exact nonzero values of `add_perturbed_exact`; the exact ranks of `add_xperturbed`; in `x_exp`
the z₁⁶ coefficient, the divisibility E = z₂²b₁²E″, E = 0 at b₁ = 0, the support bound i + j ≤ 6 and the finite
frequency check, and the Gröbner bases [1] (sympy, exact over ℚ). **Numerical evidence** (mpmath, 130 digits):
`add_planted` (residual and the numerical d_G) and the numerical d_G of `add_xperturbed`.

By hand: Lemmas 1 to 4, Lemma L (#160 §2), Theorems K and E (the independence of distinct real exponentials, the
inequality of step (1), and the reduction of the cases k₂ = 0 and R₁ constant), and the independence of the functions
ρ^m e^{ωρ}, which is used as a standard fact about functions on an interval. Assumed: #65's lemma (so D = 0 is the
algebraic stage), the formula for D as typed (shared with #80, #160), real C² functions on a box, sympy 1.14.0. Not
tested: complex solutions (Theorem K's proof is algebraic and uses complex combinations of real solutions, so it is
stated for real functions), regularity beyond C².

## 9. Lemmas for Lean (#135)

Statements only; ℝ throughout.

1. **Compact form** (`ring`). 4hD(P, Q, Z, Y, h) = (Y² + Z² + 2h²)² − (2hP + Y² − Z²)² − (2hQ − 2YZ)².
2. **Normal form** (`field_simp`). For p∧q ≠ 0 in ℝ², c = αp + βq with α = (c∧q)/(p∧q), β = (p∧c)/(p∧q).
3. **Constants** (`ring`). With a″ = αa + βa′, X₁′X₂″ − X₂′X₁″ = −α(X₁X₂′ − X₂X₁′); and α ≡ 0 with W_x ≠ 0 gives a
   constant in the span (the integration of X_i′ = c_ie^{∫β}).
4. **D̃ identity** (`ring`). E(jets with a″, s eliminated) = D̃(N; α, β, Γ, δ, v), det N = W_xW_ρ.
5. **ι² coefficient** (`ring`). With Γ = γ + ι, the coefficient of ι² in D̃ is −n₂⁴.
6. **Independence.** For distinct (m, ω) ∈ ℤ × ℂ, the functions ρ ↦ ρ^m e^{ωρ} are linearly independent on an open
   interval of (0, ∞).
7. **Lemma L.** If Σ_{i<r} F_i(x)G_i(ρ) = 0 on a box then dim span{F_i} + dim span{G_i} ≤ r (#160 §9, statement 1).
8. **Theorem K**, from 4 to 6 and the frequency count (4 − i)(μ₂ − μ₁) = 0 ⟹ i = 4.
9. **Top frequency** (`ring`, `nlinarith`). The z₁⁶ coefficient of E at X₂ = 0 is a₁²(a₁² + 4k₁²a₀²)(a₁² + 5k₁²a₀²),
   which is zero only if a₁ = 0; and for 0 ≤ j ≤ 6 − i, i < 6, (6 − i)q = j has no solution with |q| ≥ 1, q ≠ 1.
10. **Theorem E.** Gröbner certificates of the eight unit ideals of §5b (`polyrith`/`linear_combination` once the
    cofactors are exported; not exported here).

## 10. What `det_classify.md` and `two_scalars.md` should now say

For `det_classify.md`'s Answer table, a new row after the multiplicative one: "sum of two products b = X₁R₁ + X₂R₂ + v
with W_xW_ρ ≠ 0: if 1 ∈ V_X and 1 ∈ V_R, one product plus a constant (#160) or additive (#80); with constant ODE
coefficients on the ρ-side, or with the x-side X₁ = e^{k₁x}, X₂ = e^{k₂x}, k₁ ≠ k₂ real, no solution (#169, *proved*);
otherwise **not decided** (#169)." For §12's "Not covered", add: "sums of two separable products, outside the strata of
[`sum_products.md`](sum_products.md); the rank bound of #160 gives 59 to 102 products and does not decide." For
`two_scalars.md`'s "Not covered" bullet on separable shifts: the same sentence. The remark in `det_classify.md` §5 that
the multiplicative case is decided (#160) stays as it is.
