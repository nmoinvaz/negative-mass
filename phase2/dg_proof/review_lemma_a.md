# Review of Lemma A, L1, T1 and T2 of RESULTS.md for every l (#138)

Written from `RESULTS.md` alone (§1 Lemma A and L1, §3 Corollary T1/T2, the notation of §4 and the static Zerilli V of §4.2),
with sympy 1.14.0 and exact rationals (`fractions.Fraction`); no floats and no tolerances. Raw folder
`negative-mass-raw/I138-rpi5` (criteria, planted problems, mutation targets written before the first run; logs).
`chi_lemma.py`, `lemmas_all_l.py`, `threshold_lemmas.py` and their outputs were opened only after the final runs (step 7, §8 below).
**Shared with #2:** the formulas typed from RESULTS.md (the closed form of c_n, V of §4.2 with f = 1 + 2/r, the T1 bracket and the T2
statement). No code, no derived formula and no output of #2 was used before step 7.

## 1. Verdicts

| item | verdict | status word | evidence | wording notes (§9) |
|---|---|---|---|---|
| Lemma A(i) signs of c_n | **confirmed** | *proved* (by hand); c_1..c_4 *checked symbolically* | hand proof (§3); c_1..c_4 symbolic in l; signs of every c_n exact for l = 2..300 | none |
| Lemma A(ii) ψ₀ solves static Zerilli | **confirmed** | *checked symbolically* (identity in (l, m)), reductions by hand | P₂, P₁, P₀ derived from V; recurrence on rᵐ; polynomial identity in (l, m) exactly 0; ψ₀,ₓₓ = Vψ₀ exact, l = 2..40 | W1, W2 |
| L1 exponents 0 and 3 at r_s | **confirmed** | *checked symbolically* | local data and exponents symbolic in c; χ′/χ and χ″(r_s) = 0 from the ODE; exact for l = 2..40 | W3 |
| T1 N > 0 on [0, r_s] | **confirmed** | *proved* (by hand), identity and minimum *checked symbolically* | termwise identity for generic coefficients; cubic bracket ≥ 9/4 + (27/4)ε²/(1+ε)² ≥ 9 − 27/4, symbolic; Sturm, l = 2..40 | W4 |
| T2 one root of N beyond r_s | **confirmed** | *proved* (by hand), Taylor identity *checked symbolically*, d_m > 0 from A(i) | Taylor identity symbolic; d_m > 0 for 4 ≤ m ≤ l+2 exact for l = 2..120; one sign change; exact Sturm, l = 2..40 | W5, W6 |

No error and no gap in the mathematics. Wording notes are listed in §9. Every check can fail: 23 mutations, 40 (mutation, script) pairs, all
fire on the target written before the first control run (§7).

## 2. Notation

c = (l+2)(l−1), L = l(l+1), K = cL = (l−1)l(l+1)(l+2), r_s = 6/c, D = cr − 6, f = (r+2)/r, ε = 1/c, ρ = cr.
c_n = (n−2)((n−4)L + n−1)(l−n+3)_{2n−4}/(2ⁿ((n−1)!)²), with (x)_{−2} = 1/((x−1)(x−2)) for n = 1. χ = Σ c_n rⁿ.
N = cχ − (cr−6)χ′. V = 2f[λ²r²((λ+1)r − 3) + 9(λr − 1)]/(r³(λr − 3)²), λ = c/2.
This V is RESULTS' static Zerilli potential for M = −1. It is the form of the standard Zerilli potential
2f[λ²(λ+1)r³ + 3λ²Mr² + 9λM²r + 9M³]/(r³(λr + 3M)²) at M = −1, as I recall that formula (*unverified* against the
publisher; the review takes V from RESULTS as given, and the DG intertwining is outside its scope).

## 3. Lemma A(i)

Hand proof (`review_lemma_a_signs.py` docstring; checks A1.1–A1.4):
- c_1 = 3/2 (Pochhammer convention), c_2 = 0 (factor n−2), c_3 = (2−L)L/32 = −K/32 since K = L(L−2) > 0, c_4 = K/96.
- For n ≥ 4: n−2 > 0 and (n−4)L + n−1 > 0. The Pochhammer is the product of the 2n−4 consecutive integers l−n+3, …, l+n−2. Its top
  factor l+n−2 > 0. If n ≤ l+2 the bottom factor is ≥ 1, so the product is > 0. If n ≥ l+3 the range contains 0, so it is 0.
- So c_n > 0 for 4 ≤ n ≤ l+2, c_n = 0 for n ≥ l+3, and c_{l+2} > 0 because l+2 ≥ 4. deg χ = l+2; c_3 is the only negative coefficient.

Checks: A1.1 c_1..c_4 symbolic in l; A1.2 both factors have all-positive coefficients after n = 4+a, l = 2+b; A1.3 sign of the Pochhammer, exact,
l = 2..300, n = 3..l+6; A1.4 sign of every c_n, exact, l = 2..300, n = 1..l+6.

## 4. Lemma A(ii)

Hand derivation (`review_lemma_a_recur.py` docstring). With ψ = χ/D, ψ_xx = f(fψ′)′ and V as above,

  ψ_xx − Vψ = −(r+2)/(r⁴D²) · (P₂χ″ + P₁χ′ + P₀χ),
  P₂ = −r²(r+2)(cr−6), P₁ = 2cr³ + 6cr² − 12r, P₀ = c²r² − 6cr + 12.

P₂ is RESULTS'. P₀ uses −c³r³ + 12c²r² − 48cr + 72 = −(cr−6)(c²r² − 6cr + 12). On rⁿ the operator gives
c(c − n(n−3)) rⁿ⁺² + (n−1)(6n − 2c(n−3)) rⁿ⁺¹ + 12(n−1)² rⁿ, hence the coefficient of rᵐ is

  12(m−1)² c_m + (m−2)(6(m−1) − 2c(m−4)) c_{m−1} + c(c − (m−2)(m−5)) c_{m−2} = 0,  c − (m−2)(m−5) = (l−m+4)(l+m−3).

Division by Π = (l−m+5)_{2m−8} uses (l−m+3)_{2m−4} = (l−m+3)(l−m+4)(l+m−3)(l+m−2)Π and (l−m+4)_{2m−6} = (l−m+4)(l+m−3)Π, and by
2ᵐ((m−3)!)². The result is a rational function of (l, m) that is 0. It covers 4 ≤ m ≤ l+4 (Π > 0 there). m = 1..4 are checked directly
(c_0 = c_{−1} = 0). For m ≥ l+5 all of c_m, c_{m−1}, c_{m−2} vanish. The recurrence with c_1 = 3/2 determines every c_m uniquely
(12(m−1)² ≠ 0 for m ≥ 2), so the closed form is the unique power-series solution with c_1 = 3/2.

| check | result |
|---|---|
| A2.1 P's computed by sympy from V (χ a generic function) equal the hand forms; the prefactor identity holds | exact, symbolic in c |
| A2.2 P₂ as RESULTS states it | exact |
| A2.3a, A2.3b recurrence extracted from the ODE (symbolic m) equals mine and RESULTS' form | exact polynomial identity in (l, m) |
| A2.4 m = 1, 2, 3, 4 | residual 0 symbolic in l |
| A2.5 identity after dividing by (l−m+5)_{2m−8} | exactly 0 |
| A2.6 Pochhammer reductions on the grid l = 2..40, 4 ≤ m ≤ l+4; Π > 0; vanishing for m ≥ l+5 | exact |
| A2.7 recurrence with the actual c_n, m = 0..l+6 | exactly 0, l = 2..60 |
| A2.8 P₂χ″ + P₁χ′ + P₀χ as a polynomial | zero polynomial, l = 2..40 |
| Z1 ψ₀,ₓₓ = Vψ₀ from the definitions (no P's), rational functions over Q | exactly 0, l = 2..40 (11 s) |

This replaces the coordinator's numerical evidence in #104 (residual below 4e-61 at 49 values of l) by an exact identity in (l, m).

## 5. L1

Hand derivation (`review_lemma_a_l1.py` docstring), with c r_s = 6 and p = P₂′(r_s), q = P₁(r_s):
p = −6r_s(r_s+2) ≠ 0, q = 12r_s(r_s+2), P₀(r_s) = 12, P₂(r_s) = 0. The indicial equation is ρ(ρ−1) + (q/p)ρ = 0 (P₀/P₂ has a simple
pole, which does not enter), so the exponents are 0 and 1 − q/p = 3. The ODE at r_s gives qχ′ + 12χ = 0, i.e.
χ′(r_s) = −χ(r_s)/(r_s(r_s+2)) with no assumption on χ(r_s). The once-differentiated ODE at r_s gives
(p+q)χ″ + (P₁′+P₀)χ′ + P₀′χ = 0 with p+q = 6r_s(r_s+2), P₁′(r_s)+P₀(r_s) = 36(r_s+2), P₀′(r_s) = 6c = 36/r_s, hence χ″(r_s) = 0.
The twice-differentiated ODE has the coefficient 2p+q = 0 on χ‴ (the resonance of 0 and 3), and the rest vanishes with the two relations:
the exponent-0 solution has no logarithm, as it must since χ is a polynomial.

| check | result |
|---|---|
| L1.1a exponents {0, 3} (routine tested on four planted ODEs: {0,3}, {0,7}, {0,1/2}, {0,0}) | exact, symbolic in c |
| L1.1b p, q, P₀(r_s), P₂(r_s) | exact |
| L1.2 χ′/χ(r_s) = −1/(r_s(r_s+2)) | exact |
| L1.3 χ″(r_s) = 0 | exact |
| L1.4 resonance coefficient 2p+q = 0 and remainder 0 (no logarithm) | exact |
| L1.5 residue of ψ₀²dx/dr at r_s is 0; φ₀,ₓₓ(r_s) = 0 | exact |
| L1.6 the actual χ: χ(r_s) > 0, ratio, χ″(r_s) = 0 | exact, l = 2..40 |

## 6. T1 and T2

**T1.** N = Σ c_n[c(1−n)rⁿ + 6n rⁿ⁻¹]. Terms: n = 1 gives 9, n = 2 gives 0, n = 3 gives 18c₃r² − 2cc₃r³, n = 4 gives 6c₄r³ + 3c₄r³(6−cr),
n ≥ 5 gives (n−1)c_n rⁿ⁻¹(6−cr) + 6c_n rⁿ⁻¹. So N − [9 + 18c₃r² + (6c₄ − 2cc₃)r³] = Σ_{m≥4}(m−1)c_m rᵐ⁻¹(6−cr) + Σ_{j≥5}6c_j rʲ⁻¹ ≥ 0
on [0, r_s] by A(i). The bracket is B(ρ) = 9 − (9/16)(1+2ε)ρ² + (1/16)(1+2ε)(1+ε)ρ³ with B′ = (3/16)(1+2ε)ρ((1+ε)ρ − 6): decreasing on
(0, 6/(1+ε)), increasing after, and 6/(1+ε) < 6. So min B = B(6/(1+ε)) = 9 − (27/4)(1+2ε)/(1+ε)² = 9/4 + (27/4)ε²/(1+ε)² ≥ 9 − 27/4 = 9/4.
Hence N ≥ 9/4 on [0, r_s] (also at r_s, where N = cχ(r_s), so **χ(r_s) > 0 follows from T1**).

| check | result |
|---|---|
| T1.1 identity for generic C_1 = 3/2, C_2 = 0, C_3..C_9 (symbols), and termwise for symbolic n | exact |
| T1.2 identity with the actual c_n as a polynomial in r | exact, l = 2..40 |
| T1.3a–f bracket from closed-form c₃, c₄ equals B; B′; minimum; numerator of min − (9 − 27/4) is (27/4)ε² (non-negative coefficients); exact minimum for each l | exact, symbolic in l; l = 2..40 |
| T1.4 N has no root in [0, r_s] (Sturm), N(0) = 9, N(r_s) ≥ 9/4 | exact, l = 2..40 |

**T2.** With t = r − r_s and cr − 6 = ct: N(r_s + t) = cχ − ctχ′ = cΣ(1−m)d_m tᵐ, d_m = χ⁽ᵐ⁾(r_s)/m! = Σ_{n≥m} c_n C(n,m) r_s^{n−m}.
d_0 = χ(r_s) > 0 (T1). The m = 1 term has factor 0. d_2 = 0 (L1). The m = 3 term is −2cd_3. For 4 ≤ m ≤ l+2 every term of d_m has n ≥ 4 so
c_n ≥ 0, and the term n = l+2 is present with c_{l+2} > 0, so d_m > 0. The list (+, 0, 0, −2cd_3, −, …, −) has exactly one sign change
whatever the sign of d_3 (there is at least one negative term since l+2 ≥ 4). The leading coefficient is −c(l+1)c_{l+2} < 0. Descartes gives at most one
root in t > 0, and N(r_s) > 0 > N(∞) gives at least one. N has no root on [0, r_s] by T1.

| check | result |
|---|---|
| T2.1 Taylor identity, generic χ of degree 8 | exact |
| T2.2 d_m = Σ c_n C(n,m) r_s^{n−m}, generic and actual χ | exact, l = 2..40 |
| T2.3 d_0 > 0, d_1 = −d_0/(r_s(r_s+2)), d_2 = 0 | exact, l = 2..120 |
| T2.4 d_m > 0 for 4 ≤ m ≤ l+2 | exact, l = 2..120 |
| T2.5 exactly one sign change of the coefficients; leading coefficient −c(l+1)c_{l+2} | exact, l = 2..120 |
| T2.6 exact Sturm: one root in (r_s, ∞), none in [0, r_s], N squarefree | l = 2..40 |

Observation, not used in the proof: d_3 > 0 for every l = 2..120, so the list is (+, 0, 0, −, …, −) in all those cases. #2's output has d_3 > 0 for l ≤ 12.

## 7. Planted problems and controls

Planted (run first by every script, through the same code, exit 2 if one fails): PL1/PL2 recurrence extractor on (1+x)²y″ − (N+M−1)(1+x)y′ + NMy = 0
with solutions (1+x)ᴺ, (1+x)ᴹ, numeric and symbolic in (N, M, k), and a corrupted sequence that must fail; PL3 exponents on ODEs built from
two known solutions by the Wronskian: {0,3}, {0,7}, {0,1/2} (cosh 2√x), {0,0}; PL4 y′/y, y″/y and the resonance condition on the same ODEs;
PL5 sign-change counter on four planted polynomials (counts changes, not roots); PL6 interval minimum with end-point minima; PL7 Zerilli residual on ψ = r².
All pass.

Controls (`review_lemma_a_controls.py`, 11.5 min): 7 unmutated runs exit 0; **40 of 40** mutated pairs exit non-zero on the first-FAIL id of the
raw README (written before the first control run). Mutations: 6 on the planted-problem code, `cn_n4to3` (the (n−4) → (n−3) change of the
issue), `cn_poch`, `cn_sign5`, `V_coeff`, `P2_claim`, `P1_6c_5c`, `P0_12to13`, `rec_12to11`, `rec_m4to3`, `poch_ratio`, `ratio_claim`, three on T1, three on T2.
Run with `python3 review_lemma_a_X.py --mutate=NAME`.

Incidents, kept in the raw folder, none changed a criterion: (1) the first `recur` run failed A2.3a on zero-valued shift keys of my extractor; (2) the first
`common` run crashed in the sign counter on a sympy boolean; (3) the first controls run reported 29 of 40 because my parser dropped the `[l=N]` suffix of
the id; the 11 logs already showed the first FAIL on target, and the rerun gives 40 of 40.

## 8. Step 7: comparison with #2

Opened after the final runs: `chi_lemma.py`, `lemmas_all_l.py`, `threshold_lemmas.py` and their outputs. No numerical or symbolic disagreement.

| # | what | #2 | this review |
|---|---|---|---|
| D1 | P₂, P₁, P₀ and recurrence coefficients | same (`a_j(n)` of `chi_lemma_output.txt` equal mine in l; P₁, P₀ expanded in l) | same |
| D2 | polynomial identity in (l, m) | `sp.simplify(sp.factor(...))` = 0 for m ≥ 4, low orders m = −1..4 | same factorization; m = 1..4 plus the convention c_0 = c_{−1} = 0; Pochhammer reductions also tested on a grid; Π > 0 shown |
| D3 | extra checks | Γ-ratio residual (`gammasimp`); exact residual l = 2..40, m = −3..l+7 | ψ₀,ₓₓ = Vψ₀ exact for l = 2..40 from the definitions; residual exact l = 2..60; P's rederived from V with generic χ |
| D4 | χ against DG's polynomial | closed form equals `common.chi_poly` for l = 2..12 | not compared; the closed form is the definition (and the unique solution with c_1 = 3/2) |
| D5 | exponents 0, 3 and χ″(r_s) = 0 | in `threshold_lemmas.py` T2, symbolic in l; `lemmas_all_l.py` L1 uses only the ratio, φ₀,ₓₓ(r_s) and the residue | same results; also the resonance (no-log) check |
| D6 | T1 identity | exact l = 2..12 only | generic coefficients and termwise (all l), exact l = 2..40 |
| D7 | T1 minimum | p(u₀) = 9 − (27/4)c(c+2)/(c+1)², u₀ = c/(c+1) (so ρ* = 6u₀) | same value; derivative factorization shows the minimum on [0, 6] is at ρ* |
| D8 | T2 d_m > 0 and root count | by argument only; d_3 printed for l ≤ 12 | exact for l ≤ 120; Sturm l ≤ 40 |
| D9 | d_0 > 0 | from χ > 0 (Lemma A(iii)) | from T1 (N(r_s) ≥ 9/4), so T2 does not depend on A(iii) |
| D10 | `threshold_lemmas.py` docstring | says T2 "needs χ″(r_s), χ‴(r_s) ≥ 0, obtained in closed form" | the code and its output say χ‴(r_s) is a free parameter; the argument needs only that d_3 has either sign. The docstring is stale |

## 9. Wording notes for RESULTS.md (not edited here)

- **W1** (§1 Lemma A(ii)) The middle coefficient of the recurrence is elided as "(…)". It is c(m−4) − 3(m−1), i.e. the recurrence reads
  12(m−1)²c_m − 2(m−2)(c(m−4) − 3(m−1))c_{m−1} + c(l−m+4)(l+m−3)c_{m−2} = 0.
- **W2** (§1 Lemma A(ii)) The identity covers 4 ≤ m ≤ l+4, because Π > 0 exactly there (this includes the termination step m = l+3, l+4). The orders
  m ≤ 3 are checked directly, and for m ≥ l+5 all three coefficients vanish. The text says only "low orders m ≤ 4".
- **W3** (§1 L1) "This forces … χ″(r_s) = 0" is imprecise. The exponents give χ′/χ from the ODE at r_s. χ″(r_s) = 0 comes from the
  differentiated ODE at r_s, which uses P₁′(r_s) + P₀(r_s) = 36(r_s+2) and P₀′(r_s) = 36/r_s. The exponents alone do not force it.
- **W4** (§3 T1) The identity is elementary and termwise (it holds for any coefficients with c_1 = 3/2, c_2 = 0), so it holds for every l; §5 says it
  was run for l ≤ 12. The bound is 9/4 + (27/4)ε²/(1+ε)², which exceeds 9 − 27/4 for every finite l.
- **W5** (§3 T2) "d_m > 0" needs the formula d_m = Σ_{n≥m} c_n C(n,m) r_s^{n−m} and c_{l+2} > 0; "d_0 > 0" can cite T1.
- **W6** (§3 T2) "Descartes' rule then gives exactly one sign change" misstates the step. The coefficient list has exactly one sign change; Descartes then
  gives at most one root beyond r_s, and N(r_s) > 0 > N(∞) (leading coefficient −c(l+1)c_{l+2}) gives at least one.

Proposed §4.6 bullet text (for the coordinator): "**Independent review of the inputs (#138).** `review_lemma_a.md` (sympy and exact rationals,
written from this file alone) confirms Lemma A(i) and (ii), L1, T1 and T2 for every l. Lemma A(ii) is an exact polynomial identity in (l, m) after the division by
(l−m+5)_{2m−8}, with ψ₀,ₓₓ = Vψ₀ exact for l = 2..40; T1's cubic bracket is ≥ 9/4 + (27/4)ε²/(1+ε)²; T2's coefficient list has one sign change whatever
the sign of d_3. No error; six wording notes (W1–W6)." The #104 bullet's last sentence ("It did not re-prove T1, T2, or Lemma A(ii) for l > 14") then points to this.

## 10. Hand lemmas for Lean (#135)

Each is a statement for l ≥ 2 an integer, r > 0 real, c = (l+2)(l−1), r_s = 6/c. "Scripted" names the checks that test it.

| id | statement | scripted |
|---|---|---|
| H1 | c_1 = 3/2, c_2 = 0, c_3 = −K/32, c_4 = K/96 | A1.1 |
| H2 | for integers n ≥ 4: n ≤ l+2 ⇒ c_n > 0; n ≥ l+3 ⇒ c_n = 0. So c_{l+2} > 0, deg χ = l+2 | A1.2–A1.4 |
| H3 | for ψ = χ/(cr−6) and any C² χ: f(fψ′)′ − Vψ = −(r+2)/(r⁴(cr−6)²)(P₂χ″ + P₁χ′ + P₀χ) (identity of rational functions over ℚ(c)) | A2.1, Z1 |
| H4 | (P₂ ∂² + P₁ ∂ + P₀) rⁿ = 12(n−1)² rⁿ + (n−1)(6n − 2c(n−3)) rⁿ⁺¹ + c(c − n(n−3)) rⁿ⁺² | A2.3 |
| H5 | for integers 4 ≤ m ≤ l+4, 12(m−1)²c_m + (m−2)(6(m−1) − 2c(m−4))c_{m−1} + c(c−(m−2)(m−5))c_{m−2} = 0 (Pochhammer reductions, Π > 0, rational identity in (l, m)) | A2.5, A2.6 |
| H6 | the same for m ≤ 4 (c_0 = c_{−1} = 0) and, trivially, for m ≥ l+5. Hence P₂χ″ + P₁χ′ + P₀χ = 0 | A2.4, A2.7, A2.8 |
| H7 | P₂(r_s) = 0, P₂′(r_s) = −6r_s(r_s+2), P₁(r_s) = 12r_s(r_s+2), P₀(r_s) = 12; exponents {0, 3}; and P₁′(r_s) + P₀(r_s) = 36(r_s+2), P₀′(r_s) = 36/r_s | L1.1 for the first four values and the exponents. The last two were checked in a scratch session (not in the outputs); their consequence χ″(r_s) = 0 is scripted (L1.3, L1.4) |
| H8 | every C³ solution of the ODE has y′(r_s) = −y(r_s)/(r_s(r_s+2)) and y″(r_s) = 0 | L1.2–L1.4 |
| H9 | for any real C_3..C_m: N − [9 + 18C_3r² + (6C_4 − 2cC_3)r³] = Σ_{k≥4}(k−1)C_k rᵏ⁻¹(6−cr) + Σ_{j≥5}6C_j rʲ⁻¹ | T1.1, T1.2 |
| H10 | with ρ = cr, ε = 1/c: B(ρ) ≥ 9/4 + (27/4)ε²/(1+ε)² for ρ ∈ [0, 6] (B′ factorization, minimum at 6/(1+ε)); so N ≥ 9/4 > 0 on [0, r_s] | T1.3 |
| H11 | N(r_s + t) = cΣ(1−m)d_m tᵐ, d_m = Σ_{n≥m} c_n C(n,m) r_s^{n−m} | T2.1, T2.2 |
| H12 | d_0 = N(r_s)/c > 0, d_2 = 0, d_m > 0 for 4 ≤ m ≤ l+2 | T2.3, T2.4 |
| H13 | a real polynomial Q(t) = d_0 − 2d_3t³ − Σ_{m=4}^{M}(m−1)d_m tᵐ with d_0 > 0, d_m > 0 (m ≥ 4), M ≥ 4 and d_3 of either sign has exactly one root on (0, ∞), and it is simple | T2.5, T2.6 test the conclusion |

H13 avoids Descartes' rule. If d_3 ≥ 0, Q′ < 0 on (0, ∞), so Q decreases strictly from d_0 > 0 to −∞. If d_3 < 0, Q′ = t²(6|d_3| − Σ m(m−1)d_m tᵐ⁻³) and the bracket
decreases strictly from 6|d_3| > 0 to −∞, so Q increases up to some t₁ and then decreases strictly to −∞, with Q(t₁) ≥ d_0 > 0. In both cases there is exactly one root and Q′ ≠ 0
there. This argument is by hand; the scripts test its conclusion, not the argument.

## 11. Reproduce

From `phase2/dg_proof` (sympy 1.14.0): `OMP_NUM_THREADS=1 python3 -u review_lemma_a_X.py > review_lemma_a_X_output.txt` for X in `common`, `signs`, `recur`,
`zerilli` (11 s), `l1`, `t1` (32 s), `t2` (63 s), `controls` (11.5 min; `--logdir=DIR` keeps each run's log). Each of the first seven exits 0.

## Coordinator's review

Raw `I138-mac` (negative-mass-raw). Mac, nice 10, at most 2 processes.

| step | result |
|---|---|
| rerun every `review_lemma_a_*.py` | the same outputs apart from timings, all exit 0 |
| rerun the controls (one mutation and all others) | exit 0, 40 of 40 on target, identical to the committed output; e.g. `t2_window` fires at `T2.6[l=2]` |
| the first control run's 29 of 40 | every miss was a bare id read against a target with `[l=N]`. In each log the first FAIL carries the targeted l, and the matching rule written before that run accepts it. No target changed |
| re-solve beyond this review's range | the χ-ODE exactly at l = 41, 45, 60, 100; T2's Taylor coefficients exactly at l = 150 |
| independent check `ind138.py` | run 1 crashed at K2 on a sympy call in my code (kept). Run 2, with that call fixed and the criteria unchanged, exits 0; its 5 mutations fire on target |

`ind138.py` derives the χ-ODE itself from ψ₀,ₓₓ = Vψ₀ and gets P₂ = r²(r+2)(cr−6), that is −1 times RESULTS' normalisation,
P₁ = −2r(cr² + 3cr − 6) and P₀ = −c²r² + 6cr − 12. Writing c_m and c_{m−1} as c_{m−2} times ratios of the rising
factorials, the recurrence is identically 0 as a rational function of (l, m) for m ≥ 5. Its middle coefficient is
−2(m−2)(c(m−4) − 3(m−1)), as W1 says. L1's three statements hold symbolically in c. The T1 cubic's minimum is at
u = r/r_s = 1/(1+ε), with value 9/4 + (27/4)ε²/(1+ε)². d₃ was positive at every l tested (2..30, 50, 150), so the case
d₃ < 0 of H13 never arose, though the argument covers it.

W1–W6 are applied to RESULTS.md §1 and §3, with a #138 bullet in §4.6.
