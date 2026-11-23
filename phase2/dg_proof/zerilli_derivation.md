# The static even-parity potential for M < 0, derived from the linearized Einstein equations (#143)

Author of the run: the rpi5 runner of #143. Method: sympy 1.14.0, exact rational-function algebra, written from the issue, `RESULTS.md` §4.2 and
the TeX of Dotti–Gleiser gr-qc/0604021 (DG06) eqs. (1)–(9). Of `phase2/` I opened only #142's `review_lemma_a3_common.py`, `review_lemma_a3_controls.py` and the head of
`review_lemma_a3_v.py` (harness and control format; the Zerilli form typed there is the target of E.5d, not an input). The raw folder is
`negative-mass-raw/I143-rpi5` (README with criteria, planted case and control targets, written before the runs; `slips.md`; run logs;
control logs). Scripts: `zerilli_derivation.py` (driver), `zerilli_derivation_lib.py` (dual-number Ricci code, covariant cross-check,
Legendre reduction, field arithmetic), `zerilli_derivation_controls.py`. Output: `zerilli_derivation_output.txt` (stdout, exit 0, 32 PASS, 208 s, peak rss 104 MB), controls in
`zerilli_derivation_controls_output.txt`. Every check prints `PASS id` or `FAIL id`; exit 0, 1 (a real check failed) or 2 (a planted check
failed, before any real data); `--mutate=NAME` switches on one deliberate error.

## 1. Result

From the metric g = g_Schw(M) + ε h, h in Regge–Wheeler gauge, axisymmetric, e^{st}, P_l(cos θ), **M symbolic**, the linearized vacuum
equations of even parity for l ≥ 2 and s ≠ 0 are equivalent, on r ∈ (0, ∞) ∖ {r_s}, to a two-dimensional first-order system in (K, H1).
DG06's master variable ψ = r(r−2M)/((λ+1)(λr+3M)) (H2 − rK′) + rK/(λ+1) turns it into

  −ψ,ₓₓ + V ψ = −s² ψ,  dx = dr/f, f = 1 − 2M/r, λ = (l−1)(l+2)/2,

with the coefficient of ψ,ₓ exactly 0, V independent of s, and

  **V = 2f [λ²(λ+1)r³ + 3λ²M r² + 9λM² r + 9M³] / (r³(λr + 3M)²)** exactly, symbolic in (r, l, M).

At M = −1 this is the V of RESULTS.md §4.2, 2f[λ²r²((λ+1)r − 3) + 9(λr − 1)]/(r³(λr − 3)²) with f = 1 + 2/r and λ = c/2,
exactly. Status: **checked symbolically**. The one outside input of the reduction is the definition of ψ (DG06 (9)); the potential is derived, not typed.
The planted case through the same code gives the Regge–Wheeler potential f(l(l+1)/r² − 6M/r³) exactly.

## 2. Setup

Coordinates (t, r, x = cos θ, φ), so the unit sphere is dx²/(1−x²) + (1−x²)dφ² and every expression is rational in (r, x, M, s, L), with
L = l(l+1) carried as one symbol (an identity in L holds for every l; the comparisons are repeated with L = l(l+1)). P(x) = P_l(x) is an abstract function
with (1−x²)P″ − 2xP′ + LP = 0; derivatives of order ≥ 2 are removed with it.

Even ansatz (DG06 (2), their signs): g_tt = −f(1 − εH0 P e^{st}), g_tr = εH1 P e^{st}, g_rr = f⁻¹(1 + εH2 P e^{st}), angular part r²(1 + εK P e^{st}) on the
unit sphere. Odd ansatz: g_tφ = ε h0 e^{st}(1−x²)P′, g_rφ = ε h1 e^{st}(1−x²)P′.

The first order in ε is carried by dual numbers (a₀ + εa₁, ε² = 0): the inverse metric g⁻¹ − εg⁻¹hg⁻¹, the Christoffel symbols and
R_μν = ∂_aΓ^a_μν − ∂_νΓ^a_μa + Γ^a_abΓ^b_μν − Γ^a_νbΓ^b_μa are built from their definitions (`zerilli_derivation_lib.py`).

## 3. The derivation

Odd parity (planted target). δR_xφ gives s h0 = f (f h1)′ (a division by s and by f); δR_rφ then gives a second-order ODE for h1, and δR_tφ is implied
(PL.O.2). With Ψ = f h1/r: Ψ,ₓₓ = (V_RW + s²)Ψ, V_RW = f(L/r² − 6M/r³) exactly (PL.O.3).

Even parity. The 10 components of δR_μν reduce to seven radial equations (tt, tr, rr, tx, rx, angular trace, traceless; A.2), the coefficients of P and xP′ in
each angular component being separate equations because P and xP′ are independent for L ∉ {0, 2} (A.1).

1. **E.1** The traceless angular equation is (H0 − H2)/2 = 0. So H0 = H2 =: H.
2. **E.2** δR_tr, δR_tθ, δR_rθ (the tr, tx, rx lines of the output; x = cos θ) contain only K′, H1′, H′ (no second derivatives), with coefficient determinant
   −s(r−2M)/(4r). They read (H0 = H2 = H)

   (tx) (f H1)′ = s (H + K)
   (tr) K′ = −(r−3M)/(r(r−2M)) K + L/(2sr²) H1 + H/r
   (rx) H′ = K′ + sr/(r−2M) H1 − 2M/(r(r−2M)) H

   and solve to the flow printed in the output (E.2 lines). The Einstein component δG_rr = δR_rr − ½g_rr g^{ab}δR_ab has no second derivative either, and on this flow
   it is an algebraic relation A = a_K K + a₁ H1 + a_H H = 0, with
   a_K = −[Lr(r−2M) + 2r⁴s² − 2r² + 6Mr − 6M²]/(2r²(r−2M)²), a₁ = (LM + 2r³s²)/(2sr³(r−2M)), **a_H = (λr + 3M)/(r²(r−2M))**.
3. **E.3** dA/dr = −(2r−M)/(r(r−2M)) A along the flow: the relation is preserved.
4. **E.4** On the flow, tr, tx and rx vanish by construction and the other three radial equations are multiples of A: tt = (f²/2)A, rr = A/2, angular trace = −r(r−2M)A
   (the multipliers are printed in the output). So the full linearized vacuum system ⇔ {flow, A = 0, H0 = H2}: a two-dimensional solution space. Nothing else is imposed.
5. **E.5** Eliminating H = −(a_K K + a₁H1)/a_H (division by λr + 3M), DG06 (9) becomes
   ψ = [r²K − (r−2M)H1/s]/(λr + 3M), the map (K, H1) → (ψ, ψ,ₓ) has determinant (r−2M)/(s(λr+3M)), and ψ,ₓₓ = aψ + bψ,ₓ with **b = 0** and
   a − s² = V independent of s (E.5b, E.5c). V equals Zerilli's form in (r, l, M) (E.5d) and RESULTS' V at M = −1 (E.5e).
6. **E.6** Converse. DG06 (3)–(4) with ∂_t → s, imposing ψ,ₓₓ = (V + s²)ψ, solves all seven radial equations (E.6a), and (9) inverts it (E.6b).

## 4. Divisions and assumptions

| step | divides by | where it vanishes | consequence |
|---|---|---|---|
| E.2 solve of (K′, H1′, H′) | s, r, r−2M | s = 0; r = 0; r = 2M | s = 0 (static) is not covered; r−2M > 0 for r > 0 when M < 0 |
| E.2/E.5 elimination of H from A | λr + 3M (through a_H), r−2M | r = −3M/λ | the reduction holds on r ≠ r_s |
| E.5 definition (9) | λ + 1 = L/2, λr + 3M, r−2M | L = 0 (l = 0, −1) | none for l ≥ 2 |
| E.5 V | (λr + 3M)² | r = −3M/λ | genuine double pole: (λr+3M)²V there is (L−2)²(L+1)²/18 ≠ 0 |
| master variable ψ | s, λr + 3M | | det (r−2M)/(s(λr+3M)) ≠ 0 off s = 0, r = r_s |

- **λr + 3M and r_s.** λr + 3M vanishes at r = −3M/λ, positive iff M < 0 (λ > 0 for l ≥ 2). At M = −1 it is 3/λ = 6/c, the **r_s of dg_proof**
  (c = (l+2)(l−1) = 2λ). For M > 0 the zero is at negative r, outside the domain. Two things happen at r_s: ψ is not defined there by (9) and V has a double pole.
  The flow of (K, H1, H) and A themselves are regular there; H is determined by (K, H1) off r_s only. This derivation does not treat r_s itself (DG06's Lemma 1, the regularity
  condition on ζ = (2λr − 3μ)ψ_Z at r = 3μ/(2λ), is not re-derived here).
- **s.** The flow (K′ and H′ carry 1/s through the term L H1/(2sr²), from δR_tr), the elimination of H from A, the coefficient of H1 in ψ and, in the odd case, h0 all divide by s.
  The static case s = 0 is not covered; V does not depend on s, so −ψ,ₓₓ + Vψ = −s²ψ is the equation for every s ≠ 0.
- **f.** r−2M = rf appears in every flow coefficient; for M < 0 it is positive on r > 0. For M > 0 it vanishes at the horizon, outside this domain.
- **l ≥ 2.** It enters through A.1 (P and xP′ independent iff L ∉ {0, 2}, shown symbolically and by rank for l = 2..30, and at l = 1, L = 2, P = x, xP′ = x are dependent).
  No step divides by L−2 = 2λ. For l = 1 the traceless equation is empty and H0 = H2 is not obtained.
- **Assumed, not derived here.** (i) Regge–Wheeler gauge is attainable for l ≥ 2 and keeps full generality (standard; DG06 say so). The checks here show the equations are gauge invariant
  (PL.G) but do not construct the gauge transformation. (ii) m = 0 loses no generality (spherical symmetry). (iii) The sign of M is arbitrary: nothing in the derivation
  assumes M > 0 or M < 0; the M < 0 content is where the singularities lie. (iv) Boundary conditions at r = 0 and the choice of ψ as the
  master variable (other master variables have other potentials) are not touched.

## 5. Checks

| ids | what | result |
|---|---|---|
| PL.0a–d | dual inverse; Ricci of static de Sitter = (3/a²)g; vacuum Schwarzschild M symbolic; Legendre reduction against explicit P_l, l = 2..6 | PASS |
| PL.G.even/odd | h = L_ξ g gives δR = 0, all 10 components | PASS |
| PL.X.even/odd/gauge_even | independent route δR = ∇δΓ − ∇δΓ equals the dual-number route | PASS |
| PL.O.1–4 | odd parity: angular form, reduction, V_RW = f(L/r² − 6M/r³) exactly, comparison routine | PASS |
| A.1, A.2 | independence of P and xP′; angular form of the seven components, two covariance identities | PASS |
| E.1–E.4, E.3b | H0 = H2; first-order system; algebraic relation; consistency; completeness | PASS |
| E.5a–e | master variable; b = 0; V independent of s; V = V_Z; M = −1 against RESULTS §4.2 | PASS |
| E.6a, E.6b | DG06 (3)–(4) solves everything; (9) inverts it | PASS |
| D.1–D.5 | divisors in the stated set; zero and double pole of λr+3M; s, r−2M, L−2; determinants nonzero | PASS |

Controls (`zerilli_derivation_controls.py`, 742 s, peak rss 120 MB): the unmutated baseline exits 0 with 32 PASS and no FAIL, and **17 of 18** mutated runs exit with the expected
code and the first FAIL id equal to the target written in the raw README before the first control run. Planted-code mutations (exit 2): `dual_inv` PL.0a, `leg_sign` PL.0d,
`drop_trace` PL.G.even (the dropped first-order part of −∂_νΓ^a_μa is invisible to the odd checks: they pass under it, `ctrl_drop_trace_odd.log`), `route2_drop` PL.X.even, `fsign_ans` PL.O.3,
`rw_typo` PL.O.3. Real-check mutations (exit 1): `l1` A.1, `ang_norm` A.2, `ans_rr` E.1, `skip_s` E.3, `mv_coef` and `fsign_mv` E.5b, `v_typo` E.5d, `res_typo` E.5e, `dg3_typo` E.6a,
`d2_zero` D.2, `div_set` D.1. **One miss:** `fsign_bg` (the background metric with f = 1 + 2M/r) was predicted to fire at PL.0c and fired at PL.O.3 (exit 2). The prediction was wrong:
that background is the Schwarzschild vacuum solution of mass −M, so its Ricci tensor is 0 and PL.0c cannot see the sign. The target is left as written and the line is a MISS in the output; the
control that does isolate a wrong sign of M in f against a correct background is `fsign_ans`, added before it was run. A first controls run (before `master_reduce` moved to field arithmetic and before `fsign_ans` existed; see the raw README) had 16 of its 17 pairs on target and the same miss.

What the checks share: the dual-number route and the covariant route (PL.X) share the background Christoffel symbols and h and nothing else; they agree in every component of three different h
(even ansatz, odd ansatz, a pure-gauge h). The odd planted case and the even derivation share the Ricci routine, the Legendre reduction and `master_reduce`.

## 6. Sources

- DG06, gr-qc/0604021, TeX source read: eqs. (1)–(9). Quoted: (6) V = 2(1 − 2M/r)[λ²r²((λ+1)r + 3M) + 9M²(λr + M)]/(r³(λr + 3M)²), which is V_Z above
  (expand to compare); (9) ψ_Z = r(r−2M)/((λ+1)(λr+3M)) (H₂ − r∂K/∂r) + rK/(λ+1); (3)–(4) the map ψ_Z → (K, H1, H2).
- Zerilli, PRL 24 (1970) 737 and PRD 2 (1970) 2141: the publisher pages (doi.org, link.aps.org, journals.aps.org) answered HTTP 403 on 2026-11-23. **(unverified)**:
  not read. The potential is taken from DG06 (6) and checked against the derivation; no claim about what Zerilli's papers contain.
- #142 compared V against DG06 (6), DG09 (11), Kokkotas–Schmidt (24), Berti–Cardoso–Starinets (15). This note derives V rather than comparing it.
- No novelty claim is made. The reduction itself is the standard Regge–Wheeler–Zerilli one; what this adds to the project is a derivation inside the project for symbolic M.

## 7. Lemmas stated for Lean (#135)

Not Lean; statements only. Everything below is an identity in the rational function field Q(r, M, s, L) (or an equivalence of ODE systems) and discharges with `field_simp; ring`
or a short ODE argument, except H0, which is the genuinely geometric input.

- **H0** (input, not liftable cheaply) The seven radial equations of §3 are the coefficients of P and xP′ in δR_μν[g_Schw(M) + ε h] for the RW-gauge even ansatz. Needs the Ricci
  tensor of a metric in Lean; the output file prints them.
- **LZ1** For L ∉ {0, 2}, a, b ∈ ℝ, (∀x, aP(x) + bxP′(x) = 0) with P ≠ 0 on an interval, Legendre's equation ⇒ a = b = 0.
- **LZ2** det of the coefficient matrix of (K′, H1′, H′) in (tr, tx, rx) is −s(r−2M)/(4r).
- **LZ3** On the flow, dA/dr = −(2r−M)/(r(r−2M)) A, tt = (f²/2)A, rr = A/2, trace = −r(r−2M)A (explicit equations in the output).
- **LZ4** With ψ = [r²K − (r−2M)H1/s]/(λr+3M) and the (K, H1) flow, f(fψ′)′ = (V + s²)ψ with V = V_Z. (the identity E.5b, E.5c, E.5d)
- **LZ5** V_Z(M = −1) = 2f[λ²r²((λ+1)r − 3) + 9(λr − 1)]/(r³(λr − 3)²) with f = 1 + 2/r.
- **LZ6** −3M/λ = 6/c at M = −1, c = 2λ; −3M/λ > 0 iff M < 0 for λ > 0; (λr+3M)²V at r = −3M/λ is (L−2)²(L+1)²/18.
- **LZ7** Equivalence: for r ∈ (0,∞)∖{−3M/λ}, s ≠ 0, r ≠ 2M, (K, H1, H0, H2) satisfy the seven radial equations iff H0 = H2 = H = −(a_K K + a₁H1)/a_H and (K, H1) satisfy the 2D flow;
  and then ψ of (9) solves −ψ,ₓₓ + Vψ = −s²ψ and (K, H1, H) are given by DG06 (3)–(4) from ψ (conversely, those formulas applied to a solution ψ give a solution).

## 8. What RESULTS.md and review_lemma_a3.md should say

See the closing comment of #143.
