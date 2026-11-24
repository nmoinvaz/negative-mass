# A certificate for the forward direction of the Zerilli reduction, every l, without dividing by M or s (#148)

Author of the run: the rpi5 runner of #148. Method: sympy 1.14.0, exact rational-function algebra. Written from the issue, #143's note `zerilli_derivation.md`, #145's `zerilli_static.md`,
#146's `zerilli_gauge_rs.md` (each with its Coordinator's review), and `negative-mass-raw` `I143-mac/README.md` (criterion Y3 and the observation after it) and `I145-mac/README.md`
(criterion X3 and the finding). The coordinator's scripts `ind143.py`, `ind145.py` and `ind146.py` were **not opened**. The library `zerilli_derivation_lib.py` is #143's, **imported
unchanged** (dual-number Ricci code, Legendre reduction `leg_red`, coefficient extraction `coef_eqs`, harness). The driver `zerilli_certificate.py` is new; it takes no code from the
other drivers. The raw folder is `negative-mass-raw/I148-rpi5` (README with criteria, the planted case and control targets, written before the first run, with a section for what was
added after it; run logs; launchers; the certificates as JSON). Output: `zerilli_certificate_output.txt` (stdout, exit 0, 16 PASS, 81 s, peak rss 155 MB); controls in `zerilli_certificate_controls_output.txt` (one run of the runner, 1146 s, exit 1: one pre-registered target missed, see §5).
Every check prints `PASS id` or `FAIL id`; exit 0, 1 (a real check failed) or 2 (a planted check failed, before any real data); `--mutate=NAME` switches on one deliberate error.

## 1. Result

Notation: f = 1 − 2M/r, λ = (l−1)(l+2)/2 = (L−2)/2, L = l(l+1) **one symbol**, q = (L−2)r + 6M = 2(λr + 3M), jets F_n = dⁿF/drⁿ of F ∈ {H0, H1, H2, K}, D = d/dr on jets.
The radial equations are the coefficients of the angular structures in δR_ab for g = g_Schw(M) + εh, h in Regge–Wheeler gauge (DG06 (2)), e^{st}, P(x) with (1−x²)P″ − 2xP′ + LP = 0
(`zerilli_derivation.md` §2–3; printed in the output). Five of the seven are used:

  E_tl = (H0_0 − H2_0)/2 (the traceless angular equation), E_tr = δR_tr, E_tt = δR_tt, E_rr = δR_rr, E_ang = δR_xx(P) + δR_φφ(P) (the angular trace; all as coefficients of P),

with the normalizations of #143's SPEC. With ψ = r(r−2M)/((λ+1)(λr+3M)) (H2 − rK′) + rK/(λ+1) (DG06 (9)), d/dx = f d/dr and V = 2f[λ²(λ+1)r³ + 3λ²Mr² + 9λM²r + 9M³]/(r³(λr + 3M)²):

  **ψ,ₓₓ − (V + s²)ψ = c_tl E_tl + c_tr E_tr + c_tt0 E_tt + c_tt1 D E_tt + c_rr0 E_rr + c_rr1 D E_rr + c_ang E_ang**

an exact identity in the jets, with coefficients in Q(r, M, s, L):

| coefficient | value |
|---|---|
| c_tl | −2(r−2M)/r² |
| c_tr | −8 s r²(r−2M)/(L q) |
| c_tt0 | 4r(r−2M)((L−2)r + 12M)/(L q²) |
| c_tt1 | 4r²(r−2M)/(L q) |
| c_rr0 | 4(r−2M)²(3(L−2)r² − 2(L−14)Mr − 24M²)/(L r q²) |
| c_rr1 | 4(r−2M)³/(L q) |
| c_ang | −2(r−2M)(L(L−2)r² + 4(L−2)Mr + 12M²)/(L r² q²) |

**The denominators are L, r² and q² and nothing else**: no s, no M outside q = 2(λr + 3M), no L − 2 outside q, no r − 2M. The factor q is the pole at r_s that the coordinator's certificates at l = 2, 3, 4 also have; M, s and L − 2 as factors of their own do not occur. (The tx and rx equations are not used.) Status: **checked symbolically**
(E.id: exact in Q(r, M, s, L) and at three rational points; E.L, E.M0, E.S0, S.explicit below). The same certificate gives:

- **every l ≥ 2.** L = l(l+1) ≥ 6 makes L ≠ 0 and keeps q a nonzero polynomial. The identity is exact for L symbolic and then for each integer l ≥ 2 by specialization, and was retyped from λ = (l−1)(l+2)/2 and rechecked at l = 2, 3, 4
  (L = 6, 12, 20), where the common denominators are 24r²(2r+3M)², 48r²(5r+3M)², 720r²(3r+M)² (the coordinator's have r²(λr+3M)² with other constants).
- **s = 0.** s occurs only as the factor of c_tr and in no denominator, so at s = 0 the identity is ψ,ₓₓ − Vψ ∈ span of the static equations (E.S0). #143's reduction divides by s; this one does not.
- **M = 0.** q = (L−2)r ≠ 0 for r > 0; the identity holds at M = 0, where V = L/r² (so ψ,ₓₓ = (L/r² + s²)ψ is the flat radial equation), checked with Z retyped with f = 1 (E.M0). #145's flow divides by M; this one does not.
- **Pointwise.** For r > 0, q ≠ 0 and any solution (H0, H1, H2, K) of the five equations on an interval, Z = 0 there. So ψ built from the solution by (9) obeys ψ,ₓₓ = (V + s²)ψ off r_s = −3M/λ, for every l ≥ 2, every s (including 0) and every M (including 0).

Why s and M do not enter (a reading of the certificate, from a jet-by-jet listing of each term's contribution made in a scratch script that is not in the repository; the identity itself is E.id): every
contribution of H1 to the seven terms carries a factor s, and H1_2, H1_1, H1_0 cancel among them: H1_2 between D E_tt and D E_rr, H1_1 among E_tt, E_rr and their derivatives, H1_0 with E_tr, whose H1_0 coefficient L/(2r²)
has no s, so c_tr is proportional to s and never has 1/s. The H0 jets H0_3 to H0_1 cancel among the same second-order terms, and H0_0 is closed by E_tl. M enters only through f and q.

## 2. Setup and method

Coordinates (t, r, x = cos θ, φ), as `zerilli_derivation.md` §2. The first order in ε is carried by dual numbers; P is abstract, derivatives of order ≥ 2 are removed with Legendre's equation; each component's coefficients of P and xP′ are separate equations because P and xP′
are independent for L ∉ {0, 2} (A.1, as #143's A.1; for l ≥ 2). The seven even equations are `tt`, `tr`, `rr`, `tx`, `rx`, `trace` = xx(P) + pp(P), `traceless` = xx(xP′); their jets are printed in the output.

The certificate is found by linear algebra over F = Q(r, M, s, L) (`sympy.polys.fields`). Unknowns c_ik (i over the seven equations, k = 0, 1, 2: 21 unknowns); one row per jet that occurs (18 at L = 6 and symbolic); the matrix of D^k E_i against the jets has
rank 15 (the coordinator's Y3 also reports rank 15 for l = 2 and 3). Gauss–Jordan elimination with the columns in the order O1 (k ascending; within k: traceless, tx, tr, rx, tt, rr, trace) gives a solution with the free unknowns 0: it has 7 nonzero coefficients (the table above)
and its denominators already pass E.den, so none of the other column orders in the pre-registered search (raw README) was run. The certificate is not unique (rank 15 of 21); the denominators claim is for this one. It is stored as `cert_Lsym.json`.

The check E.id does not use the solver: it differentiates the equations with sympy's `diff` on expressions in Function(r), types Z from the formulas of (9) and V, and compares in the rational function field with canonical forms, jet by jet; it also evaluates the identity at three exact rational points of (r, M, s, L) by substitution.
The one input shared by the solve and E.id is the typing of Z; its typing is tested by the mutations `v_coef`, `psi_coef` and `v_s2` (no certificate exists for a mistyped Z).

## 3. Denominators

| factor | in | vanishes at | consequence |
|---|---|---|---|
| r (exponent 2 in c_tl and c_ang, 1 in c_rr0) | c_tl, c_rr0, c_ang | r = 0 | none for r > 0 |
| q = (L−2)r + 6M (exponent 2) | c_tr, c_tt0, c_tt1, c_rr0, c_rr1, c_ang | r = −3M/λ = r_s | the identity holds on r ≠ r_s; at M = 0, q = (L−2)r; positive r_s iff M < 0 |
| L (exponent 1) | all but c_tl | L = 0 (l = 0, −1) | none for l ≥ 1 |
| s, M (alone), r − 2M, L − 2 (alone) | none | | the point of the check |

The rule applied by E.den (and PL.O.4 with r − 2M in place of q): every irreducible denominator factor is r, q, or a polynomial in L alone with no root l(l+1), l ≥ 2 integer, and not L − 2. The classifier is itself tested (PL.D.1, 11 cases,
and the mutation `den_allow_M`), and E.den fires on a target that needs 1/s (`flow_target`). A comparison: the coordinator's `I145-mac` certificates at l = 2, 3, 4 have denominators r²(λr+3M)² up to the constants 3, 1386, 225; ours have the same r- and q-structure
(r², q²) and different constants (the certificates differ; both are valid).

## 4. Checks

| id | what | result |
|---|---|---|
| PL.0c, PL.0d | Ricci of g_Schw(M) is 0 (M symbolic); `leg_red` against the explicit derivative of P_l, l = 2..6 | PASS |
| PL.D.1 | the denominator classifier on 11 constructed denominators (r−2M, s, M, L−2, L−6, L−12 refused; r, q, L, L−5 accepted) | PASS |
| PL.O.1–PL.O.4 | odd planted case: angular form, certificate exists, identity (exact and at 3 points), denominators r³ and r⁴ only (r and r−2M were allowed) | PASS |
| A.1, E.1 | P and xP′ independent; the seven even equations have the expected angular form, traceless = (H0 − H2)/2, two covariance identities | PASS |
| E.cert, E.id, E.den | the certificate exists (k ≤ 2); the identity is exact; denominators L, r, q only | PASS |
| S.explicit | δR_ab with the explicit P_2, P_3, P_4 (no Legendre reduction, no extraction) equals the reconstruction from the coefficient equations at L = 6, 12, 20: 10 components each | PASS |
| E.L | l = 2, 3, 4: Z retyped from λ, certificate and equations at L = l(l+1), exact in Q(r, M, s) | PASS |
| E.M0, E.S0 | M = 0: V(L, M = 0) = L/r², identity exact; s = 0: static identity exact | PASS |

**Planted case.** The odd sector through the same code (`certify`, `check_identity`, `den_report`): with Ψ = f h1/r and V_RW = f(L/r² − 6M/r³),

  Ψ,ₓₓ − (V_RW + s²)Ψ = −2(r−2M)²/r³ · E_rφ − 2(r−2M)(r−3M)/r⁴ · E_xφ + (r−2M)²/r³ · D E_xφ,

with the three odd equations printed in the output (E_tφ is not used). Denominators r³ and r⁴. The odd certificate is the same kind of object: exact identity, no s, no L, no r − 2M.

## 5. Controls

`zerilli_certificate_controls.py`, one complete run on the final code (raw `controls_final.log`; 1146 s; the largest process of the run reported 820 MB). The unmutated baseline exits 0 with 16 PASS. **17 of 18** required pairs exit with the code and the first FAIL id written
in the raw README before the run: 16 of the 17 pre-registered ones, and `drop_second`, which was added after control part 2 with its target written before it ran. The one miss is a pre-registered pair; it stays in the output and makes the run exit 1.

| mutation | what it changes | first FAIL | exit |
|---|---|---|---|
| `leg_sign` | Legendre reduction (library) | PL.0d | 2 |
| `den_allow_M` | the classifier also allows M | PL.D.1 | 2 |
| `odd_v`, `odd_psi` | V_RW with −5M/r³; Ψ = h1/r | PL.O.2 (no odd certificate) | 2 |
| `odd_cert` | one odd coefficient + r | PL.O.3 | 2 |
| `odd_den` | the xφ equation scaled to leading coefficient 1 in h0 (needs 1/s) | PL.O.4 | 2 |
| `l1` | A.1 at L = 2 | A.1 | 1 |
| `expl_norm` | S.explicit without the 1/(1−x²) | S.explicit | 1 |
| **`v_coef`** | V with 8λM²r in place of 9λM²r | E.cert | 1 |
| **`psi_coef`** | (9) with rK/λ in place of rK/(λ+1) | E.cert | 1 |
| `v_s2` | Z without s² | E.cert | 1 |
| **`drop_second`** | tt, rr and trace removed from the span (run at L = 6; the 3-dimensional first-order system leaves Z outside) | E.cert | 1 |
| `flow_target` | Z replaced by the K′ form of #143's flow (needs 1/s) | E.den | 1 |
| `cert_coef` | one even coefficient + r | E.id | 1 |
| `l_spec` | λ = (l−1)(l+1)/2 in the retyped Z of E.L | E.L | 1 |
| `m0_v` | V at M = 0 typed as (L+2)/r² | E.M0 | 1 |
| `s0_s1` | static target with V + 1 | E.S0 | 1 |

The three controls the task names are in bold: a changed coefficient of V and of (9) leave the span test failing (no certificate over Q(r, M, s, L)), and a dropped set of equations does (`drop_second`).

**One miss: `drop_traceless`.** Its target (written before the run, "with confidence") was E.cert; it fires at E.den. Without the traceless equation the span test still passes: the two column orders tried (O1, reversed) give rank-14 certificates, both with the
forbidden factor L − 2 in the denominators (c_rr0, c_rx, c_tt0, c_tx; the output in `controls_final_logs/drop_traceless.log`). The expectation written after control part 1 (E.den, factor L − 2) is met. The prediction was wrong. A reading, not checked separately: the traceless
equation H0 = H2 follows from the other equations through a Bianchi identity, but only with a division by L − 2. So E_tl is not needed for existence; it is what keeps L − 2 out of the denominators. (A different column order might find a certificate without
the traceless equation and without L − 2; 2 of the 23 pre-registered orders were tried, `--attempts=2`.)

**Exploratory family, at L = 6** (a single equation removed; accepted outcomes were fixed in the README: exit 0, or exit 1 at E.cert or E.den). `drop_tt`, `drop_tr`, `drop_rr`, `drop_tx`, `drop_rx` exit 0: a certificate with denominators r and 2r+3M exists without
each of them (for `drop_tt` the first column order has a forbidden factor and the second passes). `dropeq_trace` exits 1 at E.den: certificates exist (rank 15), and all 23 column orders have a factor s in a denominator. At L = 6 a factor L − 2 is the number 4,
so this family cannot see it, and nothing here says what happens at symbolic L for these five. Reading: the trace equation is what keeps 1/s out, the traceless equation what keeps 1/(L−2) out, and tt, tr, rr, tx, rx are replaceable at L = 6 (the 21 unknowns have rank 15).
`drop_trace` (the library's mutation of the Ricci code) is a different thing; the equation drop for trace is named `dropeq_trace` because of that collision, found in control part 2 (`controls_part2.log`).

The runs before the final one are kept in the raw folder: part 1 (`controls_part1.log`; `drop_traceless` stopped by me during the search over column orders), part 2 (`controls_part2.log`, the family at L = 6).

## 6. Assumed, not shown here

- The seven equations are those of #143 (its library, unchanged); their derivation from the metric and the Ricci code were checked there (PL.G, PL.X) and by the coordinator's Lichnerowicz route. S.explicit ties the Legendre-reduced coefficients to the explicit δR at l = 2, 3, 4 only.
- A.1: each coefficient of P and xP′ is its own equation for l ≥ 2. Regge–Wheeler gauge is attained for l ≥ 2 (#146). m = 0, e^{st}, the sign of M arbitrary.
- This is the forward direction: for ψ built from a solution by (9), ψ,ₓₓ = (V + s²)ψ. The converse (DG06 (3)–(4) solve the equations) is #143 E.6 for s ≠ 0 and #145 for s = 0, an identity in M. That every solution is (3) applied to a ψ, i.e. (3)∘(9) = id on solutions, is
  #143 (E.5, s ≠ 0) and #145 (s = 0); the second divides by M. **#148 does not cover (3)∘(9) = id at M = 0.**
- r = r_s itself (#146), r = 0 and r = 2M. The interval of a pointwise statement avoids q = 0 and r = 0.
- Specialization to each integer l ≥ 2 is by the identity in Q(r, M, s, L); the explicit rechecks are at l = 2, 3, 4.
- Zerilli's papers were not read (#143: publisher HTTP 403, *(unverified)*). No novelty claim is made: the certificate is a linear-algebra witness for a standard reduction.

## 7. Lemmas stated for Lean (#135)

Not Lean; statements only. LC1 to LC3 are identities in the rational function field Q(r, M, s, L) and discharge with `field_simp; ring`; H0 is the geometric input shared with #143.

- **H0** (as #143) The radial equations are the coefficients of P and xP′ in δR_μν[g_Schw(M) + εh]; the five used are printed in the output.
- **LC1** For the seven coefficients above, Σ c_ik (D^k E_i)[j] − (ψ,ₓₓ − (V + s²)ψ)[j] = 0 for every jet j of H0, H1, H2, K that occurs in the seven terms or in ψ,ₓₓ: one rational-function identity per jet.
- **LC2** (D) D acts on a linear form Σ a_n F_n as Σ (a_n′ F_n + a_n F_(n+1)); only D E_tt and D E_rr are used (k = 1).
- **LC3** (specialization) The denominators L, r, q do not vanish for l ≥ 2, r > 0, q ≠ 0; at M = 0, q = (L−2)r; at s = 0 no denominator contains s. So LC1 holds at L = l(l+1), M = 0 and s = 0 without change.
- **LC4** V(M = 0) = L/r²; V(L, M, r) is the form of the issue; (λr + 3M)² V at r = −3M/λ is (L−2)²(L+1)²/18 (as #143 LZ6).
- **LC5** (plant) Ψ,ₓₓ − (V_RW + s²)Ψ = the three-term combination above, an identity in Q(r, M, s, L) in the jets h0_0, h0_1, h1_0, h1_1, h1_2.
- **LC6** (A.1 as #143's LZ1) P and xP′ independent for L ∉ {0, 2}.

## 8. What RESULTS.md and the three notes should now say

- RESULTS.md §4.6, #145 bullet, last sentence ("For every l, a certificate free of M and s is #148."): replace by "For every l ≥ 2 it is `zerilli_certificate.md` (#148): ψ,ₓₓ − (V + s²)ψ with ψ from (9) is an exact combination of five radial equations and their first r-derivatives, with coefficients in Q(r, M, s, L) whose denominators are L, r and (L−2)r + 6M only; so the forward direction holds at s = 0 and at M = 0."
- RESULTS.md §4.6, #143 bullet ("The reduction divides by s, r − 2M and λr + 3M"): add "the forward direction itself needs no division by s or r − 2M (#148)". Caveat on scope: the forward direction holds at s = 0 and at M = 0 (#148); the converse is an identity in M (#145); (3)∘(9) = id at M = 0 is not shown.
- `zerilli_derivation.md` (#143) §4: the row "E.2 solve" and the bullet on s ("The static case s = 0 is not covered") keep their meaning for the reduction; add that the forward direction at s = 0 for every l is #148.
- `zerilli_static.md` (#145) §4, bullet "M = 0": replace "The forward direction at M = 0 is not shown here" with "It is shown for every l ≥ 2 in `zerilli_certificate.md` (#148); the flow in (K, H) still divides by M".
- `zerilli_gauge_rs.md` (#146): no change. Its assumption (v) is about the e^{st} flow at M = 0 and s = 0, which this note does not use.

## 9. Sources

- DG06 gr-qc/0604021 (TeX read in #143, #146): (2) ansatz, (6) V, (9) ψ_Z. The V of the issue is the V of #143's derived form (`zerilli_derivation.md` §1).
- #143, #145, #146 notes and `negative-mass-raw` `I143-mac`, `I145-mac` READMEs (read); `ind143.py`, `ind145.py`, `ind146.py` not opened. #143's library imported unchanged.
