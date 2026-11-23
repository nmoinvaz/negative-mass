# The static (s = 0) even-parity equation and its master variable, derived from the linearized Einstein equations (#145)

Author of the run: the rpi5 runner of #145. Method: sympy 1.14.0, exact rational-function algebra. Written from #145, #143's note `zerilli_derivation.md`,
RESULTS.md §4.2, the coordinator's observation in `negative-mass-raw/I143-mac/README.md` (the coordinator's script `ind143.py` was not opened), and
the TeX of Dotti–Gleiser gr-qc/0604021 (DG06) eqs. (2)–(4), (6), (9), fetched again on 2026-11-23. The library
`zerilli_derivation_lib.py` is #143's, **imported unchanged** (dual-number Ricci code, covariant cross-check, Legendre reduction, harness); the driver
`zerilli_static.py` is adapted from #143's `zerilli_derivation.py`, which I read. The raw folder is `negative-mass-raw/I145-rpi5` (README with
criteria, planted case and control targets, written before the runs; run logs; control logs; launchers; the DG06 TeX). Output:
`zerilli_static_output.txt` (stdout, exit 0, 35 PASS, 173 s, peak rss 100 MB); controls in `zerilli_static_controls_output.txt` (767 s, peak rss 116 MB).
Every check prints `PASS id` or `FAIL id`; exit 0, 1 (a real check failed) or 2 (a planted check failed, before any real data); `--mutate=NAME` switches
on one deliberate error.

## 1. Result

From the metric g = g_Schw(M) + ε h, h in Regge–Wheeler gauge, axisymmetric, **no time dependence**, P_l(cos θ), l ≥ 2, **M symbolic**, the linearized
vacuum equations of even parity are equivalent, for M ≠ 0 and on r ∈ (0, ∞) ∖ {r_s}, to:

- **H1 = 0.** H1 occurs only in δR_tr = (L/(2r²)) H1 and δR_tx = (f H1)′/2 and in no other equation. δR_tr alone forces H1 = 0 (L ≠ 0, r ≠ 0); δR_tx is then
  implied. At s = 0 H1 is not pure gauge and not a free function: it is zero.
- **H0 = H2 =: H** (the traceless angular equation, l ≥ 2).
- **(K, H) obey a two-dimensional first-order flow**, from δR_rx and δG_rr (λ = (l−1)(l+2)/2 = (L−2)/2, f = 1 − 2M/r):

  K′ = (λ/M) K − ((λr + 2M)/(Mr)) H,  H′ = (λ/M) K − ((λr² − 2(λ−1)Mr − 2M²)/(M r (r−2M))) H,

  with coefficient determinant −M/(2r(r−2M)). δR_tt, δR_rr and the angular trace vanish **identically** on it. Nothing else is imposed, so the static solution
  space is exactly two-dimensional (for M ≠ 0).
- **DG06's master variable (9)**, ψ = r(r−2M)/((λ+1)(λr+3M)) (H2 − rK′) + rK/(λ+1), is on the flow

  ψ = −r(λr² − 3λMr − 3M²)/((λ+1) M (λr + 3M)) K + r(r−2M)/((λ+1) M) H,

  the map (K, H) → (ψ, ψ,ₓ) has determinant −2r(r−2M)/(LM) ≠ 0, and **ψ,ₓₓ = Vψ with the ψ,ₓ coefficient exactly 0 and V equal to #143's V exactly**,
  V = 2f[λ²(λ+1)r³ + 3λ²Mr² + 9λM²r + 9M³]/(r³(λr + 3M)²), symbolic in (r, l, M); at M = −1 it equals RESULTS §4.2's V exactly.
- **Converse.** DG06 (3)–(4) at s = 0 (H1 = 0, K = f_D ψ + fψ′, H0 = H2 = (f(h_D ψ + rψ′))′ − K, with f_D and h_D the f(r) and h(r) of DG06 (4)) solves all seven radial equations whenever ψ,ₓₓ = Vψ,
  and (9)∘(3) = id on ψ and (3)∘(9) = id on (K, H): the two maps are inverse, so every static solution comes from exactly one ψ.
  The converse is an identity in M (it divides by neither M nor s).
- **Relation to #143.** (3)–(4) at s = 0 *is* #143's map evaluated at s = 0 (the formulas are polynomial in s; H1 = s(hψ + rψ′) → 0). (9) does not contain s.
  What does not survive is #143's *reduced* inverse ψ = [r²K − (r−2M)H1/s]/(λr+3M), which divides by s: at s = 0 its H1/s must be read as
  h̃1 = hψ + rψ′, the O(s) part of H1. This h̃1 is a function of (K, H) on static solutions (the s → 0 tr equation defines it), and with it (9) *equals* #143's formula
  (S.R3), the tx equation at O(s) holds on every static solution (S.R2), and h̃1 equals DG06 (3)'s hψ + rψ′ (S.R4). So **every static solution is the s → 0 limit of an
  e^{st} solution, and the static ψ is the s → 0 limit of #143's ψ** (for M ≠ 0: h̃1 divides by M).

Status: **checked symbolically** (exact simplification to 0 in Q(r, M, L), also with L = l(l+1)). The one outside input of the real reduction is the definition of ψ
(DG06 (9)); V is derived, not typed. The planted case through the same code (static odd parity) gives the Regge–Wheeler potential f(L/r² − 6M/r³) exactly.

## 2. Setup

Coordinates (t, r, x = cos θ, φ): everything rational in (r, x, M, L); L = l(l+1) is one symbol (an identity in L holds for every l, and the comparisons are repeated
with L = l(l+1)). P(x) = P_l(x) abstract with (1−x²)P″ − 2xP′ + LP = 0. The ansatz is #143's with e^{st} removed (DG06 (2), their signs): g_tt = −f(1 − εH0 P),
g_tr = εH1 P, g_rr = f⁻¹(1 + εH2 P), angular r²(1 + εK P) on the unit sphere; odd: g_tφ = εh0(1−x²)P′, g_rφ = εh1(1−x²)P′. The symbol s occurs in no static equation (S.1).
The e^{st} ansatz of #143 is rebuilt only for the relation to s → 0 (S.0, S.R*).

## 3. The derivation

**Planted case: static odd parity.** δR_rφ = ((L−2)/(2r²)) h1, so h1 = 0 for L ≠ 2; δR_xφ = (f h1)′-type (h1 jets only) is implied; δR_tφ is the only equation left,
h0″ = (Lr − 4M) h0/(r²(r−2M)). With Ψ = r h0′ − 2 h0 (the combination usually attributed to Cunningham, Price and Moncrief; its form is the one typed input of the planted case, taken from memory:
the attribution is **(unverified)**, the property that it obeys the Regge–Wheeler equation is what PL.O.3 checks), Ψ,ₓₓ = VΨ with the Ψ,ₓ coefficient exactly 0 and V = f(L/r² − 6M/r³) exactly (PL.O.3). The map (h0, h0′) → (Ψ, Ψ,ₓ) has determinant 2 − L.

**Even parity.**

1. **H.1–H.3, the H1 sector.** δR_tr = (L/(2r²)) H1 (no other jet), δR_tx = ((r−2M) r H1′ + 2M H1)/(2r²) = (f H1)′/2 (H1 jets only), and tt, rr, rx, xx, pp contain no H1 jet
   (time-reversal: H1 is the one t-odd field). So H1 decouples, and H1 = 0. The tr equation at s = 0 therefore loses the K′ relation that #143's flow takes from it
   (there K′ carries L H1/(2sr²)); (K, H) must come from other equations.
2. **E.1** The traceless angular equation is (H0 − H2)/2: H0 = H2 = H.
3. **E.2** δR_rx is first order, and δG_rr = δR_rr/2 + δR_tt/(2f²) − (angular trace)/(2fr²) has no second derivative (K″ cancels between δR_rr and the trace). Solved for (K′, H′)
   they give the flow of §1; the determinant of the coefficient matrix is −M/(2r(r−2M)). *The flow divides by M*: at M = 0 the two equations are proportional in (K′, H′).
4. **E.4** On the flow δR_tt, δR_rr and the angular trace are 0 identically (no residual algebraic relation, unlike #143's A = 0 at s ≠ 0). So the full static system is
   {H1 = 0, H0 = H2, flow}.
5. **E.5** (9) with K′ from the flow is ψ = pK + qH (§1); ψ,ₓₓ = aψ + bψ,ₓ with b = 0 and a = V_Z exactly (E.5b, E.5d); at M = −1, a is RESULTS' V (E.5e).
6. **E.6** Converse and inversion as in §1.
7. **S.0, S.R1–S.R4** The static coefficient equations equal the s = 0 values of #143's e^{st} equations (all components); with H1 = s h̃1, (δR_tr/s)|_{s=0} is −[K′ + (r−3M)K/(r(r−2M)) − H/r
   − L h̃1/(2r²)], the tr line of `zerilli_derivation.md` §3. Solving it for h̃1 and differentiating along the flow, (δR_tx/s)|_{s=0} = ((f h̃1)′ − H − K)/2 vanishes (S.R2); (9) = [r²K − (r−2M)h̃1]/(λr+3M) (S.R3);
   h̃1 = hψ + rψ′ (S.R4).

## 4. Divisions and assumptions

| step | divides by | where it vanishes | consequence |
|---|---|---|---|
| H.1: H1 = 0 | L (through c_tr = L/(2r²)), r | L = 0 (l = 0, −1); r = 0 | none for l ≥ 1, r > 0 |
| E.1: H0 = H2 | none (c₀ = 1/2); l ≥ 2 enters through A.1 | | |
| E.2 flow | **M**, r, r − 2M | M = 0; r = 0; r = 2M | M = 0 (flat) is **not covered** by this reduction; r − 2M > 0 for M < 0 |
| E.5 definition (9) | λ + 1 = L/2, λr + 3M | L = 0; r = −3M/λ | r_s, below |
| E.6 DG06 (4), f(r) and h(r) | r, r − 2M, λr + 3M | | the converse divides by no M and no s |
| E.5 p, q, det | M, L, λr + 3M | | det −2r(r−2M)/(LM) ≠ 0 off r = 0, 2M, M = 0 |
| E.5 V | r, (λr + 3M)² | r = 0; r = −3M/λ | a genuine double pole at r_s: (λr+3M)²V there is (L−2)²(L+1)²/18 ≠ 0 |
| S.R h̃1 | M, L, r − 2M | | the s → 0 relation needs M ≠ 0 |

No step divides by **s** (S.1, D.3) or by **L − 2 = 2λ** (D.4). Compared with #143's reduction: the division by s is gone, the divisions by r − 2M and λr + 3M stay, and
**M is new**, because the flow in (K, H) replaces #143's flow in (K, H1) (there the algebraic relation A has the H coefficient (λr+3M)/(r²(r−2M)) and nothing divides by M).

- **r_s.** λr + 3M = 0 at r = −3M/λ, positive iff M < 0 (λ > 0); at M = −1 it is 6/c = r_s of dg_proof. ψ is not defined there by (9) and V has a double pole. The flow of (K, H) is regular there
  (its denominators are M, r, r − 2M). This note does not treat r_s itself (DG06's Lemma 1 is not re-derived).
- **M = 0.** The converse holds for every M (an identity). The forward direction at M = 0 is **not** shown here: the flow in (K, H) is singular, and a different elimination (using the second-order
  equations) is needed. The coordinator's certificate `ind143.py` Y3 for l = 2 (M symbolic) has the denominator 3r²(2r+3M)² only, so the forward direction holds at M = 0 for l = 2 (the l = 3 certificate is at M = −1 only);
  for every l it is a recommended task of the closing comment of #145. The project's case is M < 0, where M ≠ 0.
- **Assumed, not derived here** (as in #143): (i) Regge–Wheeler gauge is attainable and keeps full generality for l ≥ 2 (the checks show gauge invariance of the equations, PL.G; they do not
  construct ξ; this is #146's task, and the static case needs the same algebraic solve); (ii) m = 0 loses no generality; (iii) the sign of M is arbitrary in the derivation; (iv) boundary conditions
  at r = 0 and the choice of ψ as the master variable are not touched; (v) "static" means a time-independent perturbation (no secular t-growth); l ≥ 2.

## 5. Checks

| ids | what | result |
|---|---|---|
| PL.0a–d | dual inverse; Ricci of static de Sitter = (3/a²)g; vacuum Schwarzschild, M symbolic; Legendre reduction against P_l, l = 2..6 | PASS |
| PL.G, PL.X | static pure gauge (even and odd ξ) gives δR = 0, all 10 components; the covariant route equals the dual-number route in every component for the even, odd and gauge h | PASS |
| PL.O.1–4 | static odd: angular form; h1 = 0, h0 ODE, implied equations; V_RW = f(L/r² − 6M/r³) exact for Ψ = r h0′ − 2h0; comparison routine | PASS |
| A.1, A.2, S.1 | independence of P and xP′; angular form of the seven components and the two covariance identities; no s in the nine static coefficient equations | PASS |
| H.1–H.3 | δR_tr = (L/(2r²)) H1; δR_tx has H1 jets only; no H1 jet elsewhere | PASS |
| E.1, E.2, E.4 | H0 = H2; the static flow (determinant −M/(2r(r−2M))); tt, rr, trace vanish on it | PASS |
| E.5a–e | (9) is ψ = pK + qH with nonzero determinant; b = 0; a = V_Z; M = −1 against RESULTS §4.2 | PASS |
| E.6a, E.6b | DG06 (3)–(4) at s = 0 solves all seven equations; (9)∘(3) = id and (3)∘(9) = id | PASS |
| S.0, S.R1–S.R4 | static = s = 0 of the e^{st} equations; the s → 0 tr and tx equations; (9) = #143's ψ with h̃1; h̃1 = hψ + rψ′ of (3) | PASS |
| D.1–D.5 | divisors in the vocabulary; zero and double pole of λr + 3M; no s; no L − 2; determinants nonzero | PASS |

Controls (`zerilli_static_controls.py`, 767 s, peak rss 116 MB): the unmutated baseline exits 0 with 35 PASS and no FAIL, and **all 39 mutated runs exit with the expected code (2 at a planted id, 1 at
a real id) and the first FAIL id equal to the target written in the raw README before the first run** (`--failfast`; the README was edited before the script's first run, for `tt_drop_K` → `tt_drop_H`, which would have
been inert, and for the vocabulary of D.1). One mutation, S.R4 with `h1p_typo`, was added after the exploratory runs and before the first control run. Every check of the table has its own mutation:
planted `dual_inv` PL.0a, `ds_typo` PL.0b, `bg_nonvac` PL.0c, `leg_sign` PL.0d, `drop_trace` PL.G (a dropped term of the linearized tensor: the first-order part of −∂_νΓ^a_μa), `route2_drop` PL.X,
`ang_norm_odd` PL.O.1, `inj_rphi` PL.O.2, `rw_typo`, `fsign_ans` and `cpm_typo` PL.O.3, `cmp_blind` PL.O.4; real `l1` A.1, `ang_norm` A.2, `s_leak` S.1, `inj_tr` H.1, `inj_tx` H.2, `inj_h1` H.3,
`ans_rr` E.1, `grr_drop` E.2, `tt_drop_H` E.4 (a dropped term of δR_tt), `psi_zero` E.5a, `mv_coef`, `fsign_mv` and `mv_lam` E.5b (a changed coefficient of (9)), `v_typo` E.5d (a changed coefficient of V),
`res_typo` E.5e, `dg3_typo` E.6a, `inv_typo` E.6b, `dep_typo` S.0, `tr_typo` S.R1, `h1t_typo` S.R2, `sr_typo` S.R3, `h1p_typo` S.R4, `div_set` D.1, `d2_zero` D.2, `div_s` D.3, `div_l2` D.4, `d5_zero` D.5.
No miss. Two predictions rested on structure (`grr_drop` at E.2: K″ does not cancel without the trace term; `tt_drop_H` at E.4: the flow changes while its determinant does not); both held.

What the checks share: the odd planted case and the even derivation share #143's Ricci routine, the Legendre reduction and the field arithmetic; the covariant route (PL.X) shares the background Christoffel symbols and h and
nothing else. A second operator (the coordinator's Lichnerowicz form, `ind143.py`) has not been run on the static equations; S.0 shows that the static equations are #143's at s = 0, which the coordinator's
check of #143 confirmed for the e^{st} equations.

## 6. Sources

- DG06, gr-qc/0604021, TeX source read (`I145-rpi5/dg06_src/`): eqs. (2) the ansatz; (3)–(4) the map ψ_Z → (K, H1, H2) with f(r) and h(r); (6) V; (9) ψ_Z. The formulas typed in the driver agree with the TeX.
- #143 `zerilli_derivation.md`: V_Z (comparison target), the s ≠ 0 reduction (typed as targets of S.R1 and S.R3 only), the `zerilli_derivation_lib.py` library (imported).
- Regge–Wheeler potential V_RW = f(L/r² − 6M/r³) (given in the task) and the Cunningham–Price–Moncrief combination r h0′ − 2h0 as the static odd master variable: **not read**, recalled **(unverified)**.
- Zerilli, PRL 24 (1970) 737 and PRD 2 (1970) 2141: not read (publisher HTTP 403 in #143) **(unverified)**; nothing here depends on them.
- No novelty claim is made. The static even-parity reduction to the Zerilli equation is standard; what this adds is a derivation inside the project for symbolic M at s = 0.

## 7. Lemmas stated for Lean (#135)

Not Lean; statements only. Each is an identity in the rational function field Q(r, M, L) (or a short equivalence of ODE systems), discharged by `field_simp; ring`, except H0.

- **H0** (input, as in #143) The static radial equations printed in the output are the coefficients of P and xP′ (P′) in δR_μν[g_Schw(M) + εh] for the static RW-gauge even ansatz.
- **LS1** δR_tr = (L/(2r²)) H1, so H1 = 0 for L ≠ 0, r ≠ 0; then δR_tx = (f H1)′/2 = 0.
- **LS2** The coefficient matrix of (K′, H′) in (δR_rx, δG_rr) has determinant −M/(2r(r−2M)).
- **LS3** On the flow of §1, δR_tt = δR_rr = (angular trace) = 0 (explicit equations in the output).
- **LS4** With ψ = pK + qH (§1) and the flow, f(fψ′)′ = Vψ with V = V_Z (the identity E.5b, E.5d).
- **LS5** DG06 (3)–(4) at s = 0 with f(fψ′)′ = Vψ satisfy the seven equations (E.6a); (9)∘(3) = id and (3)∘(9) = id (E.6b).
- **LS6** V_Z(M = −1) is RESULTS' V (as LZ5 of #143); (λr+3M)²V at r = −3M/λ is (L−2)²(L+1)²/18.
- **LS7** Equivalence: for r ∈ (0,∞) ∖ {−3M/λ, 2M}, M ≠ 0, L ≥ 6: (K, H1, H0, H2) satisfy the static radial equations iff H1 = 0, H0 = H2 = H and (K, H) solve the flow; and then ψ of (9) solves
  f(fψ′)′ = Vψ and (K, H) are given by (3)–(4) at s = 0 from ψ.
- **LS8** (planted) V_RW for Ψ = r h0′ − 2h0 and h0″ = (Lr − 4M) h0/(r²(r−2M)).

## 8. What RESULTS.md should now say

The closing comment of #145 proposed four edits. All four are applied: RESULTS.md §4.2 (after V's definition), §4.6 (the
#143 bullet and a #145 bullet), §6 (Caveat on scope), and `zerilli_derivation.md` §4. §4.6 records the coordinator's
finding on M = 0 below.

## Coordinator's review

Raw `I145-mac` (negative-mass-raw). Mac, nice 10, at most 2 processes.

| step | result |
|---|---|
| scope | the five files the issue names, nothing else |
| rerun `zerilli_static.py` | exit 0, 35 PASS, 0 FAIL; output equal to rpi5's apart from the timing lines |
| rerun the controls | 40 OK, 0 MISS, exit 0; rows equal to rpi5's apart from timings |
| re-solve of one case | l = 2, M symbolic, s = 0, forward: ψ,ₓₓ − Vψ with ψ from (9) is a combination of the static radial equations and their first two r-derivatives, with coefficients whose only denominator is r²(2r + 3M)² |
| independent check `ind145.py` | run 1 exits 0; all 6 mutations fire on target |

`ind145.py` is #145's recommended task 4, the second operator. It reuses the coordinator's `I143-mac/ind143.py` (the
Lichnerowicz form of δR, DG06 typed from the TeX, the span certificates) and shares no code with this note's scripts or
#143's library. The flow it tests is typed from the closing comment.

| id | what | run 1 |
|---|---|---|
| X0 | static planted: l = 1 odd, h0 = 1/r (Kerr, linear in the spin) and r² (rigid rotation) give δR = 0, 1/r² does not; static even pure gauge gives δR = 0 | pass |
| X1 | δR_tr = L H1 P/(2r²), δR_tx = (f H1)′P′/2, H1 in no other component, l = 2, 3, 4 | pass |
| X2 | converse: DG06 (3)–(4) at s = 0 solve all 10 components, l = 2, 3, 4 | pass |
| X3 | forward: H1, H0 − H2 and ψ,ₓₓ − Vψ lie in the span of the static equations, l = 2, 3, 4, with M symbolic and at M = 0 | pass |
| X4 | the flow of §1 solves all 10 components (l = 2, 3, 4), and every solution obeys it (certificates, l = 2, 3) | pass |
| X5 | (3)∘(9) is the identity on solutions, l = 2, 3 | pass |

So the static reduction is confirmed by a second operator, and no number of #145 disagrees. One finding answers part of
§4's open point on M = 0. X3's forward certificates have denominators r²(λr + 3M)² only, with no M, and X3 passes at
M = 0 directly. So for l = 2, 3, 4 the forward direction holds at M = 0, though the flow (whose certificates in X4 do
carry M) does not. For every l this is #148.

Wording: §1's converse wrote DG06 (4)'s f(r) as F, which the note does not define. It now reads f_D and h_D. Recommended
task 1 is #148, task 2 is added to #146 (the static case of the gauge solve), task 3 is row 18 of #135, and task 4 is
this review's `ind145.py`.
