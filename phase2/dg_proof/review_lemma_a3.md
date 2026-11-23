# Independent review of dg_proof's Lemma A(iii), L2–L6 and the Zerilli V (#142)

Reviewer: the rpi5 runner of #142. Method: sympy 1.14.0 and exact rationals only, written from `RESULTS.md` alone (§1,
§2, §4.2, the notation of §3–§4), then compared with #2's `chi_lemma.py`, `lemmas_all_l.py` and `verify_small_l.py` in
step 7. The raw folder is `negative-mass-raw/I142-rpi5` (README with criteria, planted problems, sources and the control
targets, all written before the runs; `slips.md`; step logs; control logs). Every script prints `PASS id` / `FAIL id`
and exits 0, 1 (a real check failed) or 2 (a planted check failed, before any real data); `--mutate=NAME` switches on
one deliberate error.

## 1. Verdicts

| item | verdict | evidence (check ids) |
|---|---|---|
| Lemma A(iii): χ > 0 on r > 0 | **confirmed**; two wording points (W1, W2) | A3.1–A3.9 |
| L2: r⁴V̂ = 12 + 2(l²+l+3)r + O(r²); Frobenius exponents 3, −1; limit point at 0 | **confirmed**; wording points W3, W4 | E2.1–E2.5, E2.7, E2.8 |
| L3: V̂ ≃ (l+1)(l+2)/r² at ∞; limit point at ∞ | **confirmed**; wording point W5 | E3.1, E3.2 |
| L4: H_ZΨ_DG = −(k/2)²Ψ_DG, k = c(c+2)/6 | **confirmed**; sharpened to "iff" (W6) | L4.1–L4.5 |
| L6: Ψ̂ = A*Ψ_DG is O(r³) at 0 and regular at r_s | **confirmed**; the exact leading term is (2ᵏk²/192)r³ (W7) | L6.1–L6.6 |
| V of §4.2 (static Zerilli potential, M = −1) | **confirmed** against four published forms; Zerilli's own papers **(unverified)** | V1.1–V1.9 |

No **gap** and no **error** in RESULTS.md. One **disagreement with a source**, not with RESULTS: the 1/x term of V̂
printed in Dotti–Gleiser (arXiv:0809.3615) eq. (61) differs from the exact coefficient (§4); RESULTS does not use it.

RESULTS' table lists "L1–L6", but §1 labels L1, L2/L3, L4 and L6 and gives the inequality of Lemma A(iii) no L-number.
#2's `lemmas_all_l.py` calls that inequality L5 (W8). The issue's "L2–L6" is therefore L2, L3, L4, (L5 = A(iii)) and L6.

## 2. What was checked

Notation as RESULTS: c = (l+2)(l−1), K = (l−1)l(l+1)(l+2) = c(c+2), k = K/6, r_s = 6/c, f = 1 + 2/r, ∂ₓ = f∂ᵣ, ψ₀ =
χ/(cr−6), g = f ψ₀,ᵣ/ψ₀, V̂ = g² − g,ₓ. Each check runs symbolically in l (or in c, K, k, s) where the lemma is for all l, and exactly with the actual χ for
the l-range named at the check.

### Lemma A(iii) (`review_lemma_a3_a3.py`, 114 PASS in 42 s)

- Hand derivation before code: c₅ = K(K−24)/6144 because (l−2)₆ = K(L−6) and (L−6)(L+4) = K−24 (L = l(l+1)). For q = 3/2
  − a r² + b r³ + d r⁴, a = K/32, b = K/96, d = c₅: q′ = r(4d r² + 3b r − 2a), discriminant 9b² + 32ad = (Ks/192)² with
  s = √(6(K−18)), and the positive root is **r₀ = 24/(s+6) = 4(s−6)/(K−24)**, the other root −24/(s−6).
- A3.1 c₃, c₄, c₅ from RESULTS' c_n equal −K/32, K/96, K(K−24)/6144 (polynomial identities in l; c₁ = 3/2, c₂ = 0).
- A3.2 K(2) = 24, K(l) ≥ 120 for 3 ≤ l ≤ 50, K−24 = (l−2)(l+3)(l²+l+4), c₅(2) = 0, c₅(l) > 0 for 3 ≤ l ≤ 50.
- A3.3[l] for l = 3..50, χ/r − q has zero coefficients at r⁰..r⁴ and nonnegative ones above (this uses Lemma A(i),
  reviewed in #138).
- A3.4 q′ = 4d r (r − r₀)(r − r₁) with r₀ = 24/(s+6), r₁ = −24/(s−6), exactly (symbolic in s); the product of the roots
  is −96/(K−24) < 0. So q′ < 0 on (0, r₀), q′ > 0 after, and r₀ is the global minimum on [0, ∞). This is DG's "unique
  critical point".
- A3.5 q(r₀) = 4(K−18)(Ks − 18(K−16))/(K−24)³, exactly, symbolic in s, with K = s²/6 + 18; r₀·q(r₀) equals DG09 eq. (59)
  as printed (typed from the paper, reduced modulo ρ² = 6K−108); r₀ equals DG09 (58).
- A3.6 (Ks)² − (18(K−16))² = 6(K−24)³, a polynomial identity in K once s² = 6(K−18).
- A3.7 for l = 3..50: K > 24, K > 18, 18(K−16) > 0, and (Ks)² > (18(K−16))² in exact integers, so Ks > 18(K−16) and
  q(r₀) > 0. Where the hypothesis is used: the identity gives 6(K−24)³ > 0 iff K > 24. The step from "difference of
  squares > 0" to Ks > 18(K−16) needs Ks > 0 and 18(K−16) > 0 (K > 16), and the prefactor 4(K−18)/(K−24)³ > 0 needs
  K > 24. **K > 24 is all that is used; K ≥ 120 is the range l ≥ 3.** (A3.7x: q(r₀) = 0.505718664756 at l = 3,
  0.364805798301 at l = 4, 0.0841548156722 at l = 10, 0.00383278973872 at l = 50, equal to the closed form to 40
  digits.)
- A3.8 l = 2 is separate: K = 24, c₅ = 0, χ/r = 3/2 − (3/4)r² + (1/4)r³ with minimum 1/2 at r = 2. The closed form for
  q(r₀) is 0/0 at K = 24 (s = 6, (K−24)³ = 0); its limits are r₀ = 2 and q(r₀) = 1/2, the l = 2 values.
- A3.9[l] exact Sturm sequence (own sign counter, cross-checked by sympy `count_roots`): χ/r has **no root in (0, ∞) for
  l = 2..50**. Planted PLS.1–6 (six polynomials with known positive-root counts, including a double root and a
  polynomial whose only real root is negative) run first.

### L2 and L3 (`review_lemma_a3_ends.py`, 204 PASS in 122 s; `review_lemma_a3_dg61.py`, 6 PASS and 1 DIFFERS in 3 s)

- Hand derivation before code: ψ_r/ψ = 1/r + c/6 + O(r), g = 2/r² + (1+c/3)/r + …, V̂ = 12/r⁴ + 2(c+5)/r³ + …, c+5 =
  l²+l+3; indicial 4s(s−2) − 12 = 0.
- E2.1 with a generic χ = r(3/2 + h₂r² + h₃r³): r⁴V̂ = 12 + 2(c+5)r + ((c²+5c+32h₂+6)/3)r² + …; with h₂ = −K/32 the r²
  coefficient is **l(l+1)**.
- E2.2[l], l = 2..30: the same expansion exactly with the actual χ; the pole order at 0 is exactly 4. E2.3[l]: the
  reduced denominator of V̂ is const·r⁴(χ/r)² (so V̂ is continuous on (0, ∞) because of A3.9, and V̂ → +∞ at 0, → 0 at
  ∞, hence bounded below).
- E2.4: the indicial polynomial of f(f uᵣ)ᵣ = (V̂−E)u is s² − 2s − 3 for every E, symbolically (V̂ = N/r⁴, N(0) = 12)
  and exactly for l = 2..30. The routine returns the known exponents on five planted ODEs (three built from known
  solutions by the Wronskian determinant). E2.5: x = r²/4 − r³/12 + …, so the exponents in x are 3/2, −1/2; ∫(r⁻¹)²dx
  diverges at 0 and ∫(r³)²dx converges.
- E2.7[l], l = 2..50: φ₀ = (cr−6)/χ ~ −4/r at 0, ~ (c/c_{l+2})r^{−(l+1)} at ∞, φ₀²/f ~ 8/r at 0 (not in L²) and
  integrable at ∞, φ₀′(r_s) = c/χ(r_s) > 0. E2.8: φ₀,ₓₓ = V̂φ₀ for a generic ψ, and exactly for l = 2..20.
- E3.1: ψ = rⁿ(1 + a₁/r + a₂/r² + a₃/r³), n symbolic: r²V̂ → n(n+1). E3.2[l]: deg ψ₀ = l+1 and r²V̂ → (l+1)(l+2) exactly
  for l = 2..30.
- **DG09 (61) in x** (dg61): the series of V̂ in x is 3/(4x²) + (l²+l−1)/(4x^{3/2}) + **0**/x + O(x^{−1/2}), exact in l
  (E2.6b) and numeric (E2.6c, r = 10⁻⁸, l = 2, 3, 5). DG09 prints −ℓ(ℓ+1)/(4|M|x) for the 1/x term. The printed value is
  the 1/x contribution of the first two r-terms alone (12/r⁴ and 2(l²+l+3)/r³); the O(r⁻²) term l(l+1)/r² adds
  +l(l+1)/(4x). E2.6 (the printed coefficients) therefore reports `DIFFERS` at x⁻¹ and passes at x⁻² and x⁻³ᐟ².

### L4 (`review_lemma_a3_l4.py`, 43 PASS in 70 s)

- L4.1 the intertwining algebra for a generic ψ(x): AA* = −∂² + ψ″/ψ, A*A = −∂² + g² − g′, ĤA* = A*H_Z, H_ZA = AĤ
  (symbolic).
- L4.2 the DG identity with k kept symbolic: [−Ψ_DG,ₓₓ + VΨ_DG + (k/2)²Ψ_DG]/Ψ_DG has numerator **(r+2)(c² + 2c − 6k)**
  (up to the constant of the printed factorization), so it vanishes for every c exactly when 6k = c(c+2), and not at k +
  1 (symbolic in r and c, so for every l). Planted PLL.2: Ψ_k = (r+2)ᵏe^{−kr/2} = 2ᵏe^{−kx/2} with V = 0 has eigenvalue
  −(k/2)².
- L4.3 k = K/6 is an integer divisible by 4 for l = 2..200. L4.4[l] ψ₀,ₓₓ = Vψ₀ exactly, l = 2..20 (this is H_Z = AA*,
  Lemma A(ii) of #138).
- L4.5[l] ĤΨ̂ = −(k/2)²Ψ̂, Ψ̂ ≢ 0, exactly for l = 2..20, from the definitions Ψ̂ = −Ψ_DG,ₓ + gΨ_DG and V̂ = g² − g,ₓ,
  with no intertwining and no formula for the Wronskian.

### L6 (`review_lemma_a3_l6.py`, 74 PASS in 9 s)

- L6.1 Ψ̂ = A*Ψ_DG = −f(W′χ − Wχ′)/((cr−6)χ), W = r(r+2)ᵏe^{−kr/2}, for generic W and χ (symbolic).
- L6.2 W/2ᵏ = r − (k/8)r³ + (k/24)r⁴ + (k²/128 − k/64)r⁵ + …
- L6.3 with c₂ = 0, c₃ = −K/32, c₄ = K/96, c₅ = K(K−24)/6144 and K = 6k, the coefficients of r¹..r⁴ of W′χ − Wχ′ are 0
  and that of r⁵ is 3·2ᵏk²/128, so **Ψ̂ = (2ᵏk²/192) r³ + O(r⁴)**, symbolic in k. Hand: p = 3 needs c₃ = −3k/16, p = 4
  needs c₄ = k/16, p = 5 uses c₅.
- L6.4[l], l = 2..30: the same exactly for the actual χ and k. L6.5: W′/W(r_s) = −1/(r_s(r_s+2)) = χ′/χ(r_s) (L1, #138)
  exactly when k = c(c+2)/6 (symbolic in c; fails for k+1). L6.6[l]: the reduced rational part of Ψ̂ has no pole at r_s
  and none on (0, ∞) (Sturm), l = 2..30, so Ψ̂ is smooth on (0, ∞).

### V (`review_lemma_a3_v.py`, 15 PASS in 2 s)

Sources, read as data (README of the raw folder): DG06 = Dotti–Gleiser, CQG 23 (2006) 5063, arXiv:gr-qc/0604021; DG09 =
Dotti–Gleiser, CQG 26 (2009) 215002, arXiv:0809.3615; KS = Kokkotas–Schmidt, Living Rev. Relativ. 2 (1999) 2,
arXiv:gr-qc/9909058; BCS = Berti–Cardoso–Starinets, CQG 26 (2009) 163001, arXiv:0905.2975. Zerilli's original papers
(PRL 24, 737 and PRD 2, 2141, 1970): the publisher page answered HTTP 403, **not read, (unverified)**; the four sources
above quote them.

- V1.1 RESULTS' V at M = −1, 2f[λ²r²((λ+1)r−3) + 9(λr−1)]/(r³(λr−3)²), equals each of DG06 (6), DG09 (11), KS (24) and
  BCS (15) (L → ∞) at M = −1, exactly, symbolic in r and λ (and in l). V1.2 the four published forms agree with each
  other for generic M.
- V1.3 x = r − 2 ln(1 + r/2) has dx/dr = 1/f and x(0) = 0, and equals DG06 (8) and DG09 (13) at M = −1. DG09 (13) prints
  ln((r−2M)/2M), whose argument is negative for M < 0; read as a modulus.
- V1.4 Chandrasekhar's relation V^± = ±β F,ₓ + β²F² + κF, F = Δ/(r³(μ²r + 6M)), β = 6M, κ = μ²(μ²+2), μ² = c, holds
  exactly: V⁺ = V_Z and V⁻ = f(l(l+1)/r² − 6M/r³) = BCS (14), generic M. The book was not read **(unverified)**; the
  identity itself is exact. It derives V_Z from the Regge–Wheeler potential with no input from the Zerilli polynomial.
- V1.5 V̂ from the definition equals DG09 (54) once χ″ is eliminated with DG09 (53) (M = −1, symbolic in c). V1.6 the
  χ-ODE derived from V and ψ₀ equals DG09 (53). V1.8 c_n of RESULTS equal DG09 (55) for l = 2..30, every n, and DG09
  (56), (57) (c₅ = ℓ(ℓ²+ℓ+4)(ℓ+3)(ℓ²−4)(ℓ²−1)/6144 as printed).
- V1.7 V is a rational function of M, so M = −1 is an ordinary substitution; the pole λr + 3M = 0 sits at r_s = 6/c > 0
  for M = −1 and at −6/c for M = +1.
- V1.9 DG09's Fig. 1 caption puts the double pole of the l = 2, M = −2 potential at x ≃ 0.761; the formulas give r_s = 3
  and x = 3 − 4 ln(7/4) = 0.76154.

What V1 does not cover: that the static even-parity perturbation equations of the negative-mass Schwarzschild metric
reduce to this V (DG06 eqs. (2)–(6) assert it for either sign of M; I did not derive it from the linearized Einstein
equations), and that Ĥ is the operator of the physical problem (DG's intertwining equivalence, which RESULTS §6 already
excludes).

## 3. Controls

Planted problems (run before any real data, through the same code): PLA.1–3 (quartic with known critical point), PLS.1–6
(Sturm counter), PLE.1a–e (indicial routine on five ODEs with known exponents), PLE.2 (V̂ routine on ψ = rⁿ with hand
values 4n(n+2) and n(n+1)), PLE.3 (order/leading-coefficient routine), PLE.4 (Taylor routine), PLW.1–3 (Wronskian
order), PLL.1–2 (V̂ for ψ = r by hand; eigen-residual for Ψ_k), PLV.1–2 (comparison routine on an equal and an unequal
pair).

`review_lemma_a3_controls.py`: the 6 unmutated runs exit 0 (dg61 with exactly one `DIFFERS`), and **54 of 54** mutated
(mutation, script) pairs exit non-zero with the first FAIL id equal to the target written in the raw README before the
first control run (44 distinct mutations; 9 of them change a planted-checked routine, and 13 pairs end with exit 2 at a
planted id). 723 s of run time, 60 OK, 0 MISS; the same run repeated in step 6 gives identical rows. Peak rss 125 MB
(ends), one process, nice 10, systemd units for the longer runs.

| script | PASS | seconds |
|---|---|---|
| a3 | 114 | 42 |
| ends | 204 | 122 |
| dg61 | 6 (+1 DIFFERS) | 3 |
| l6 | 74 | 9 |
| l4 | 43 | 70 |
| v | 15 | 2 |
| controls | 60 OK | 723 |

My slips (all in `slips.md`, none changed a pass criterion): a sympy `series` call with a negative order; my `k_of`
helper used Python float division, which broke exactness for l ≥ 17 in the first L6 run; a `nsimplify` of an exact
rational, a typing slip (`subs` after `cancel`), and a `simplify` of `Abs` of a complex symbol that hung for 10 minutes
(stopped). In my step 5 post I gave a wrong mutation count (23 and 8), which I had not measured; I corrected that
comment from the controls script (44 distinct, 9 planted-code). Changes after the exploratory runs and before the
controls: E2.6 moved out of the ends script into `review_lemma_a3_dg61.py` (it disagrees with DG09 (61), reported as
`DIFFERS`), E2.6b and E2.6c added, and planted Sturm, order and Taylor checks added where the routines are used.

## 4. Discrepancies

- **DG09 (61), the 1/x term of V̂ in x.** Printed −ℓ(ℓ+1)/(4|M|x); exact value 0 (E2.6b, E2.6c). The printed value omits
  the O(r⁻²) term's contribution (+l(l+1)/(4x)). RESULTS uses only r⁴V̂ = 12 + 2(l²+l+3)r + O(r²), which agrees with the
  first line of DG's (61), and the x^{3/2}, x^{−1/2} exponents, which are unaffected. No action on an external site.
- DG09 (13) as printed takes a log of a negative number for M < 0 (read as a modulus). Cosmetic.

## 5. Wording points for RESULTS.md (none changes a conclusion)

- **W1 (A(iii)).** After "(Ks)² − (18(K−16))² = 6(K−24)³ > 0 for K ≥ 120" add the two steps that turn it into q(r₀) > 0:
  Ks > 0 and 18(K−16) > 0, so Ks > 18(K−16); and 4(K−18)/(K−24)³ > 0. All three need only K > 24. Say that K ≥ 120 is
  the range l ≥ 3. Add r₀ = 24/(s+6) (= DG's (58)).
- **W2 (A(iii), l = 2).** Say that the closed form for q(r₀) is 0/0 at K = 24 (limits r₀ = 2, q = 1/2), which is why l =
  2 is separate.
- **W3 (L2).** The next coefficient is l(l+1): r⁴V̂ = 12 + 2(l²+l+3)r + l(l+1)r² + O(r³). The exponents 3 and −1 differ
  by the integer 4, so the −1 branch may carry a term ∝ x^{3/2} ln x; it does not affect "not in L²".
- **W4 (L2, limit point at 0).** Two routes work: the Frobenius exponent −1 (as written), or V̂ ≥ 3/(4x²) near 0, since
  V̂ = 3/(4x²) + (l²+l−1)/(4x^{3/2}) + O(x⁻¹ᐟ²) with (l²+l−1) > 0 (E2.6b; the theorem number of the criterion in
  Reed–Simon, X.10 as I recall it, is *(unverified)*).
- **W5 (L3).** Name the criterion for ∞: V̂ ≥ 0 near ∞, hence bounded below there (Reed–Simon X.8 as I recall it,
  *(unverified)*).
- **W6 (L4).** Add that the identity holds exactly when 6k = c(c+2): the numerator is (r+2)(c² + 2c − 6k). Add why Ψ̂ ∈
  D(Ĥ): Ĥ is essentially self-adjoint, so D(Ĥ) is the maximal domain {u ∈ L², u ∈ AC_loc, Ĥu ∈ L²}, and Ψ̂ is smooth on
  (0, ∞) (L6.6), in L² at 0 and ∞, with ĤΨ̂ = −(k/2)²Ψ̂ ∈ L². The O(r³) is more than membership needs.
- **W7 (L6).** Ψ̂ = (2ᵏk²/192)r³ + O(r⁴). Say that k = c(c+2)/6 enters twice: in L4 and in the regularity at r_s (W′/W =
  χ′/χ there). Say that c₅ = K(K−24)/6144 enters the r³ coefficient.
- **W8 (labels).** Give the Lemma A(iii) inequality its table name L5, or drop "L1–L6" from the table's
  `lemmas_all_l.py` row.

## 6. Hand steps as lemmas for Lean (#135), no Lean written

- **LL1.** For l ≥ 2: c₃ = −K/32, c₄ = K/96, c₅ = K(K−24)/6144, and K − 24 = (l−2)(l+3)(l²+l+4). (A3.1, A3.2; polynomial
  identities after unfolding the Pochhammer symbols.)
- **LL2.** For real K > 24 and s = √(6(K−18)): with a = K/32, b = K/96, d = K(K−24)/6144, 4d r² + 3b r − 2a = 4d (r −
  24/(s+6))(r + 24/(s−6)); hence q′ = 4d r (r−r₀)(r−r₁) and q ≥ q(r₀) on r ≥ 0. (A3.4.)
- **LL3.** q(r₀) = 4(K−18)(Ks − 18(K−16))/(K−24)³ and (Ks)² − (18(K−16))² = 6(K−24)³; for K > 24 this gives q(r₀) > 0.
  (A3.5–A3.7.)
- **LL4.** 3/2 − (3/4)r² + (1/4)r³ = 1/2 + (r−2)²(r+1)/4 (l = 2). (A3.8; the identity is r³ − 3r² + 4 = (r−2)²(r+1).)
- **LL5.** For l ≥ 3 the coefficients of χ/r − q are ≥ 0 (c_n ≥ 0 for n ≥ 6, Lemma A(i)), so χ/r ≥ q(r₀) > 0 on r > 0.
  (A3.3; combines with LL3.)
- **LL6.** ψ = r h/(cr−6), h = 3/2 + h₂r² + h₃r³: r⁴(g² − g,ₓ) = 12 + 2(c+5)r + ((c²+5c+32h₂+6)/3)r² + O(r³), where ∂ₓ =
  (1 + 2/r)∂ᵣ. (E2.1.)
- **LL7.** For V̂ = N(r)/r⁴ with N analytic, N(0) = 12, the ODE f(f uᵣ)ᵣ = (V̂ − E)u has indicial polynomial s² − 2s − 3
  at r = 0. (E2.4.)
- **LL8.** ψ = rⁿ(1 + a₁/r + a₂/r² + a₃/r³): r²(g² − g,ₓ) → n(n+1) as r → ∞; deg ψ₀ = l+1. (E3.1, E3.2.)
- **LL9.** For ψ twice differentiable and nonzero, with g = ψ′/ψ: (∂+g)(−∂+g) = −∂² + ψ″/ψ and (−∂+g)(∂+g) = −∂² + g² −
  g′. (L4.1.)
- **LL10.** In (r, c, k) with V of §4.2 and Ψ_DG = Eρ, E = (r+2)ᵏe^{−kr/2}, ρ = r/(cr−6): the residual numerator equals
  a nonzero rational function times (c² + 2c − 6k). (L4.2.)
- **LL11.** Wronskian coefficients: with w = r − (k/8)r³ + (k/24)r⁴ + (k²/128 − k/64)r⁵ and χ = (3/2)r − (K/32)r³ +
  (K/96)r⁴ + (K(K−24)/6144)r⁵, K = 6k: the r¹..r⁴ coefficients of w′χ − wχ′ vanish and the r⁵ coefficient is 3k²/128.
  (L6.3.)
- **LL12.** W′/W(r_s) = 1/r_s + k/(r_s+2) − k/2 = −1/(r_s(r_s+2)) with r_s = 6/c iff k = c(c+2)/6. (L6.5.)
- **LL13.** V of RESULTS = DG09 (11) at M = −1 as rational functions in (r, λ): `field_simp; ring`. (V1.1.) Optional:
  Chandrasekhar's identity (V1.4).
- Not computer algebra: limit point at 0 and ∞ (Weyl alternative with the exponents or a bound on V̂), D(Ĥ) as the
  maximal domain, Glazman's lemma. These stay textbook theorems.

## 7. Comparison with #2 (step 7)

I opened `chi_lemma.py`, `lemmas_all_l.py`, `verify_small_l.py` and their outputs only after the final runs. **No number
of #2 disagrees with mine.**

| point | #2 | this review |
|---|---|---|
| r₀ and q(r₀) | r₀ from the quadratic, equal to DG (58) at six l; q(r₀) = A₀ + A₁s reduced modulo s² = 6l(l²−1)(l+2) − 108, equal to the closed form; l = 3: 0.505718664756, l = 4: 0.364805798301, l = 10: 0.0841548156722 | r₀ = 24/(s+6) in K, s; closed form in s; DG (58), (59) typed from the paper; the same three values |
| (Ks)² − (18(K−16))² = 6(K−24)³ | factored by sympy; text "> 0 for K > 24", RESULTS says K ≥ 120 | polynomial identity, plus the sign steps A3.7 |
| χ > 0 | exact root count `count_roots(0, None)` of χ/r for l = 2..10 (`verify_small_l`) | own Sturm counter and sympy for l = 2..50 |
| r⁴V̂ | 12 + 2(l²+l+3)r + O(r²), exponents 3, −1 | the same, and the r² coefficient l(l+1); DG (61) comparison |
| V̂ at ∞ | lim r²V̂ = (l+1)(l+2) with ψ = r^{l+1}(1 + d₁/r) | same limit, n symbolic with three subleading terms, and exact for l = 2..30 |
| L4 | (H_ZΨ)/Ψ with k = c(c+2)/6 substituted after the simplification, equal to −(k/2)² | numerator factored, (r+2)(c²+2c−6k): "iff" |
| L6 | (ln h)′ − (ln χ)′ = 0 at r_s; X₁′/X₁ − H₁′/H₁ = 0 + 0·r + 0·r² + O(r³) with a symbolic c₅ | exact r³ coefficient 2ᵏk²/192; l = 2..30 exact |
| H Ψ̂ = −(k/2)²Ψ̂ | l = 2..10 from a derived formula for Ψ̂ | l = 2..20 from the definition A*Ψ_DG, and the formula checked as L6.1 |

Shared machinery: #2's exact checks take χ from `dg_spectrum/common.py` (`chi_poly`) and check the closed form of c_n
against it for l = 2..12; this review typed c_n from RESULTS and shares nothing with it. What this review shares with
#2: the formulas typed from RESULTS.md (c_n, V of §4.2, f, x, ψ₀, g, V̂, Ψ_DG, k).
#2 does not check: the sign steps after the identity; χ > 0 beyond l = 10; V against any source; the r² coefficient; the
"iff" in L4.

## 8. What RESULTS.md §4.6 and §6 should now say

§4.6, replace the last sentence of the #138 bullet ("Not reviewed: Lemma A(iii), L2–L6, and the Zerilli V itself, …")
and add: "**Independent review of the remaining inputs (#142).** `review_lemma_a3.md` (sympy and exact rationals,
written from this file alone) confirms Lemma A(iii) (closed form of q(r₀) and (Ks)² − (18(K−16))² = 6(K−24)³ exactly,
χ/r without positive root by exact Sturm for l = 2..50), L2, L3, L4 (the DG identity holds exactly when 6k = c(c+2)) and
L6 (Ψ̂ = (2ᵏk²/192)r³ + O(r⁴)), and finds the V of §4.2 equal at M = −1 to Dotti–Gleiser's eqs. (6) and (11),
Kokkotas–Schmidt's (24) and Berti–Cardoso–Starinets' (15). No error; eight wording notes (W1–W8). Not checked: Zerilli's
own papers (publisher page blocked, *unverified*), the derivation of V from the linearized Einstein equations for M < 0,
and DG's intertwining equivalence. Dotti–Gleiser's eq. (61) prints −ℓ(ℓ+1)/(4|M|x) for the 1/x term of V̂; the exact
coefficient is 0 (this file does not use it)." §6, "Caveat on scope": keep the sentence on DG's intertwining
equivalence; add that V itself now agrees with four published forms. Apply W1–W8 in §1. §5 may list "L2–L6, V" as
verified symbolically in l, and "χ > 0 exact Sturm for l = 2..50" next to the l ≤ 10 root counts.

## 9. Recommended tasks

1. **Apply W1–W8 and the §4.6 text to RESULTS.md.** Acceptance: each W changed or declined with a reason; §4.6 no longer
   says L2–L6 and V are unreviewed. Depends on #142.
2. **Feed LL1–LL13 into #135.** Acceptance: each lemma has a Lean statement that compiles; LL2–LL4 are proved; LL10 and
   LL11 are proved by `ring` after the stated factorization. Depends on #135.
3. **Derive the even-parity static perturbation equation for M < 0 from the linearized Einstein equations**
   (Regge–Wheeler gauge, ψ_Z from DG06 eq. (9)) and compare V with §4.2 symbolically; read Zerilli's papers through the
   publisher. Why: V is checked against published forms only. Acceptance: V equal exactly; a changed coefficient fails.
   Depends on none.

