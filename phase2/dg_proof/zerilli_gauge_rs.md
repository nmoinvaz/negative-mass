# The Regge–Wheeler gauge is attainable, and #143's flow at r_s against DG06's Lemma 1 (#146)

Author of the run: the rpi5 runner of #146. Method: sympy 1.14.0, exact rational-function algebra. Written from #146, the notes `zerilli_derivation.md` (#143) and
`zerilli_static.md` (#145) with their coordinator reviews, `review_lemma_a.md` §5 (#138, L1) and `review_lemma_a3.md` §6 (#142, L6.5, LL11, LL12), and the TeX of
Dotti–Gleiser gr-qc/0604021 (DG06). #143's library `zerilli_derivation_lib.py` is **imported unchanged** (harness, dual-number Ricci code, covariant cross-check,
Legendre reduction); the driver `zerilli_gauge_rs.py` is new, and its r_s part derives the even-parity radial equations and the flow again with that library (the code path of #143's E.1–E.2),
so #143's printed flow enters only as a comparison target (R.1). Not opened: the coordinator's scripts `ind143.py`, `ind145.py`. The raw folder is `negative-mass-raw/I146-rpi5`
(README with criteria, planted case and control targets, written before the first run; development, step and final logs; control logs; the DG06 TeX). Output:
`zerilli_gauge_rs_output.txt` (stdout, exit 0, 34 PASS, 124 s of which PL.G takes about 88 s, peak rss 104 MB); controls in `zerilli_gauge_rs_controls_output.txt` (762 s, peak rss 119 MB).
Every check prints `PASS id` or `FAIL id`; exit 0, 1 (a real check failed) or 2 (a planted check failed, before any real data); `--mutate=NAME` switches on one deliberate error.

## 1. Result

**Gauge, e^{st}.** A general even (l, m = 0) perturbation has seven functions of r (definitions in §2): H0, H1, H2, K and the three that Regge–Wheeler (RW) gauge sets to zero,
J_t (h_tθ), J_r (h_rθ) and G (the traceless angular part). For ξ^t = aPe^{st}, ξ^r = bPe^{st}, ξ_x = cP′e^{st} and h′ = h + L_ξ g, with f = 1 − 2M/r and ′ = d/dr,
the code reads off (from L_ξ g, not typed)

  ΔJ_t = −f a + s c, ΔJ_r = b/f + c′ − 2c/r, ΔG = 2c/r², ΔK = 2b/r − L c/r², ΔH0 = −(f′/f) b − 2s a, ΔH1 = −f a′ + s b/f, ΔH2 = −(f′/f) b + 2b′,

and solves J_t′ = J_r′ = G′ = 0 **exactly, symbolic in L, M and s**:

  **c = −r² G/2, b = −f (J_r − r² G′/2), a = (J_t − s r² G/2)/f.**

The four RW functions after the transformation are K′ = K + (L/2)G − (2f/r)(J_r − r²G′/2), H0′ = H0 + f′(J_r − r²G′/2) − 2s a, H1′ = H1 − f a′ − s(J_r − r²G′/2),
H2′ = H2 + f′(J_r − r²G′/2) + 2b′ (printed in full in the output, GS.1). The solve is algebraic: G fixes c, then J_r fixes b (c′ is known), then J_t fixes a. For l ≥ 2 it has **no
homogeneous solution** (GU.1), so RW gauge is attained and ξ is unique.

**What the solve divides by.** Only by r (through ΔG = 2c/r²) and by f = (r−2M)/r (through ΔJ_t). **Not by s, not by L − 2, not by M** (GS.3). L − 2 enters one level below the algebra:
the traceless structure T = xP′ − (L/2)P of G vanishes identically iff l ≤ 1 (GA.1), and P, xP′ are independent iff L ∉ {0, 2}. **At L = 2 (l = 1)** the traceless part of every h_ab is zero, G is not a
function of h, the G equation is empty, and Jt′ = Jr′ = 0 leaves a one-function family c(r) arbitrary, a = s c/f, b = −f(c′ − 2c/r) (GU.2). The mutation `l2` (GU.1 run at l = 1) fires.

**Gauge, s = 0.** The same solve closes, with a = J_t/f and the same b, c (GST.1, equal to the s = 0 value of the e^{st} solve): no ODE in r, no integration constant, and **no residual gauge freedom**:
the homogeneous problem has only ξ = 0 for l = 2..8 (GST.2), and a time-linear ξ = tκ would need L_κ g = 0, which forces κ = 0 (GST.3: the algebraic equations for (α, β, γ) have determinant 4f/r³ ≠ 0).
So the assumption (i) of #143 §4 and of #145 §4 is derived for l ≥ 2, with and without time dependence.

**r_s, e^{st}, s ≠ 0, M ≠ 0, l ≥ 2.** With μ = −2M, λ = (L−2)/2, r_s = −3M/λ = 3μ/(2λ):

1. r_s is a regular point of #143's three-dimensional flow in (K, H1, H): the denominators of its entries and of A's coefficients are s, r and r−2M (R.2), the Taylor coefficients at r_s are finite to order 4 (R.3). The coefficients are rational,
   so every solution is analytic at r_s. (ψ and V are not regular there; λr + 3M divides in them.)
2. ζ := (im) = (2λr − 3μ)ψ_Z, evaluated on the flow, is **ζ = 2r²K − 2(r−2M)H1/s**, with no H, and ζ′ has no H either. Write (ζ, ζ′)ᵀ = T (K, H1)ᵀ. Then det T = 4r(λr+3M)/s has a simple zero at r_s, T(r_s) has rank 1, and its second row is **c\*** times its first with

   **c\* = λ²/(3M(2λ+3)) = −2λ²/(3μ(2λ+3)) = M/(r_s(r_s − 2M))**, and in terms of L: (L−2)²/(12M(L+1)) (R.4).

   So **ζ′(r_s) = c\* ζ(r_s) for every analytic solution of the flow, with or without the relation A = 0** (R.4). At M = −1: c* = −λ²/(3(2λ+3)) = −1/(r_s(r_s+2)), #138's L1.2 (R.11, l = 2, 3, 4).
3. Conversely (Lemma 1, with c\*): on the corrected (nv1) of §3.3, the numerators of K and H1 and the double-pole numerator N2 of H2 are all multiples of ζ′ − c\*ζ at r_s, and N2′(r_s) = −(L−2)/(2M) N2(r_s)
   with no ζ″ in it (R.9). So **(K, H1, H2) is smooth at r_s iff ζ is smooth and ζ′(r_s) = c\* ζ(r_s)**. A concrete case (M = −1, l = 2) is R.12. This is the exact content of "the flow's regular solutions are those whose ζ satisfies (cz)": every
   solution of the flow is regular, every ζ from it satisfies the condition, and the condition is exactly what makes the inverse map regular. It is not a restriction on the solution space.
4. ζ on S = {A = 0} obeys DG06's (RW17) exactly (R.6), and the singular-point relation of that ODE at r_s is ζ′ = c\*ζ again; the ζ″ relation on S is ζ″(r_s) = −9s²/(L+1)² ζ(r_s) at M = −1, l = 2, 3, 4 (R.11; it is 0 at s = 0, #138's L1.3).
   The two-dimensional flow in (K, H1) (H eliminated by A) has a regular singular point at r_s with exponents {0, 1} and no logarithm (R.10), so all its solutions are analytic: an apparent singularity, as for (RW17) (DG06's Lemma 2).
5. **Static, s = 0.** The static flow in (K, H) (#145) is regular at r_s (denominators M, r, r−2M), det T_s = −8r²(λr+3M)²/(LM) (a double zero), T_s(r_s) has rank 1, and its rows are related by **the same c\***
   (R.S.2). So every static solution satisfies ζ′(r_s) = c\* ζ(r_s).

**Regular** means C^∞ in a neighborhood of r_s. For solutions of this flow that is the same as real-analytic (rational coefficients, no pole at r_s), hence "finite with finite derivatives of every order". **The flow alone decides it**: a linear
system with analytic coefficients at r_s has only analytic solutions (textbook theorem; its hypothesis is R.2, R.3). What the flow does not decide, and Lemma 1 does, is whether the *parametrization* by ζ is regular: ψ = ζ/(2(λr+3M)), V, and the map
ζ → (K, H1, H2) divide by λr + 3M (twice for H2), and ζ′ = c\* ζ is the condition that those divisions are exact.

**Discrepancies with DG06 (the arXiv TeX; the published text was not read, *(unverified)*), reported and not fixed.** (cz) is printed ζ′ = **+**2λ²/(3μ(2λ+3)) ζ at r = 3μ/(2λ). The derivation gives c\* = **−**2λ²/(3μ(2λ+3)) (R.5).
On solutions the printed condition would force ζ(r_s) = 0 (R.5), and a smooth ζ with the printed constant gives poles in K, H1 and H2 (R.12). DG06's own (RW17) already forces c\* at r_s (R.6), so the printed (cz) disagrees with the paper's own evolution equation.
Written with M instead of μ the same condition is c\* = λ²/(3M(2λ+3)); a plausible origin of the printed form is the substitution M = +μ/2 instead of −μ/2 (hand check of the algebra, not an inference about the authors). In (nv1) (R.7): K is right; H1 is (2/r) times the derived H1 (the printed one has
dimension 1/length); H2 has the opposite sign on its ζ term; its printed double-pole numerator vanishes at r_s exactly under the printed (cz) (a hand check, not scripted; the control `n2_sign` shows R.9 fires when the printed sign is used). (im) = (9)·(2λr − 3μ) is correct (R.8), and so are (RW17) (R.6) and (3)–(4), (9) (R.7, R.8). V of (6) is not tested here (#143).
`phase0/RESULTS.md` quotes the condition as DG08's Lemma 1 in the form c₀ = λ²c₋₁/(3M(3+2λ)) for ψ = c₋₁/(r−r_s) + c₀ + …; since ζ = 2λ(r−r_s)ψ this is ζ′(r_s) = c\*ζ(r_s) with the derived sign (hand check; DG08 itself not read here, *(unverified)*).

Status: **checked symbolically** (exact simplification to 0 in Q(r, M, s, L), also at M = −1 and l = 2, 3, 4), except the textbook steps named in §4.

## 2. Setup

Coordinates (t, r, x = cos θ, φ); unit sphere γ = dx²/(1−x²) + (1−x²)dφ²; L = l(l+1) is one symbol, λ = (L−2)/2; P(x) = P_l(x) abstract with Legendre's equation; E = e^{st} or E = 1. Perturbation of one (l, m = 0) mode,
sign conventions of DG06 (2):

| function | component | definition |
|---|---|---|
| H0 | h_tt | f H0 P E |
| H1 | h_tr | H1 P E |
| H2 | h_rr | H2 P E/f |
| J_t | h_tx (= h_tθ with dx = −sinθ dθ) | J_t P′(x) E |
| J_r | h_rx (= h_rθ) | J_r P′(x) E |
| K, G | h_ab (a, b angular) | r² E [K P γ_ab + G T_ab], T_ab = D_aD_bP + (L/2)Pγ_ab, so (1−x²)T_xx = xP′ − (L/2)P = −T_φφ/(1−x²) |

RW gauge: J_t = J_r = G = 0. My K is K_RW − (L/2)G with K_RW the Regge–Wheeler paper's; they agree in RW gauge (not read, *(unverified)*). Gauge vector: ξ^t = aPE, ξ^r = bPE, ξ_x = cP′E (lower index), ξ^φ = 0,
exactly the family of #143's `xi_even` and #145's. For l ≥ 2 and m = 0 any even vector field has this form (the odd part does not touch the seven functions); this is the decomposition into vector harmonics, a textbook fact assumed here.

## 3. The derivation

### 3.1 The gauge solve (e^{st})

- **PL.L** the partial-derivative form of L_ξ g equals ∇_μξ_ν + ∇_νξ_μ in all ten components; **PL.G** h = L_ξ g gives δR_μν = 0 (#143's PL.G.even; this also validates the Lie derivative against the Einstein operator).
- **PL.E** `extract` (Legendre reduction, coefficients of P, xP′, P′ read separately, the covariance identities β_xx + β_φφ = 0 and α_φφ − α_xx = L β_xx checked inside) returns the seven functions from an h built from them, with E = e^{st} and E = 1.
- **PL.P** (planted) η arbitrary, h = ∇η + ∇η (covariant route, independent of the partial form used by the solver): the solver returns ξ = −η and all seven functions of h + L_ξ g vanish; with an RW-gauge part added it returns that part exactly. Static and e^{st}.
- **GS.1–GS.3** the effect of ξ on the seven functions (printed), the solve, the four transformed RW functions, the closed forms of §1 (typed there as targets), and the divisors {r, r−2M}.
- **GA.1, GU.1, GU.2** T_l = xP_l′ − (L/2)P_l is a nonzero polynomial with leading coefficient −k_l l(l−1)/2 ≠ 0 for l = 2..30 (**proved** by hand for every l ≥ 2: P_l has degree l and leading coefficient k_l ≠ 0), and zero for l = 0, 1.
  GU.1: for l = 2..8 the homogeneous problem (read from the components of L_ξ g with explicit P_l) forces c = 0, then b = a = 0. GU.2: at l = 1 it leaves c(r) free.

### 3.2 The static gauge solve

The same code with E = 1 (GST.1): the effect of ξ has no s in it; the solve equals the s = 0 value of §3.1; H0′ = H0 + f′(J_r − r²G′/2), H1′ = H1 − f a′, H2′ and K′ as in §1 with s = 0 (printed). GST.2, GST.3 as in §1.

### 3.3 r_s

- **R.1** The flow of the e^{st} even ansatz is derived from δR_tr, δR_tx, δR_rx (H0 = H2 = H from the traceless equation) and A from δG_rr, and equals #143's printed one exactly; det of the coefficient matrix −s(r−2M)/(4r).
- **R.2, R.3** regular point and finite jets (§1, 1).
- **R.4** ζ and T (§1, 2): T = [[2r², −2(r−2M)/s], [2Mr/(r−2M), ((L−2)r+4M)/(rs)]].
- **R.5** the printed (cz) against c\*; the functional ζ(r_s) is not zero on S(r_s) = {aK K + a1 H1 = 0}, so the printed condition cuts S to a line.
- **R.6** From the flow with H eliminated by A, ζ″ = αζ′ + βζ with α = 2(Lr² − 3LMr − 2r² + 6Mr − 6M²)/(r(r−2M)((L−2)r+6M)), β printed in the output; it equals (RW17) with ∂_t² → s². The residues at r_s are α_res = 2, β_res = −(L−2)²/(6M(L+1)),
  and α_res c\* + β_res = 0.
- **R.7, R.8** (nv1) and (im). With ψ = ζ/(2λr−3μ) and (3)–(4): K = [(μ+r) r ζ′ + (λr − μ)ζ]/(r²(2λr−3μ)) (as printed); **H1 = s[2r(r+μ)ζ′ + μζ]/(2(r+μ)(2λr−3μ))** (printed with r in place of the 2 in the denominator);
  H2 = (r+μ)ζ″/(2λr−3μ) − (4λr²+6λμr−3μ²)ζ′/(2r(2λr−3μ)²) **−** (3μ²−6λμr+4λ²r²)ζ/(2r²(2λr−3μ)²) (printed with +). The corrected forms equal the flow's own map T⁻¹ and, on S, H = −(aK K + a1H1)/aH in terms of (ζ, ζ′).
  (9)∘(3) = id on a free ψ, and (im) = (9)·(2λr−3μ).
- **R.9–R.12** Lemma 1, the apparent singularity, M = −1 and the concrete ζ (§1, 3–4).
- **R.S.1, R.S.2** the static flow derived at s = 0 (δR_tr = L H1/(2r²), so H1 = 0; δR_rx and δG_rr give (K′, H′)) equals #145's; T_s and c\* (§1, 5).

## 4. Divisions and assumptions

| step | divides by | where it vanishes | consequence |
|---|---|---|---|
| gauge: c from ΔG | r² (ΔG = 2c/r²), and the structure T(x) as a function | r = 0; T ≡ 0 iff l ≤ 1 | none for l ≥ 2, r > 0; at l = 1 the solve is not unique (GU.2) |
| gauge: a from ΔJ_t, H0′, H1′, H2′ | f = (r−2M)/r | r = 2M | r − 2M > 0 on r > 0 when M < 0 |
| gauge: b from ΔJ_r | none (multiplies by f) | | |
| gauge, static | same as above | | no s, no L−2, no M |
| flow (R.1) | s, r, r − 2M (#143) | s = 0; r = 0, 2M | s ≠ 0 here; the static case has its own flow |
| ζ, T (R.4) | (λ+1) in (im); s | L = 0; s = 0 | none for l ≥ 2, s ≠ 0 |
| c\* and r_s | M, 2λ+3 = L+1 | M = 0 | M ≠ 0 (at M = 0, r_s = 0 is outside the domain) |
| ODE (R.6), 2D flow (R.10) | (L−2)r + 6M = 2(λr+3M) | r = r_s | the apparent singularity, exponents {0, 1} for the flow |
| corrected (nv1) (R.9) | r, r+μ, (2λr−3μ) once (K, H1) and twice (H2) | r = r_s | exactly the condition ζ′ = c\*ζ |
| static flow (#145) | M, r, r − 2M | M = 0 | as in #145 |

- **Textbook steps, not recomputed.** (a) A linear ODE system with coefficients analytic at r_s has only analytic solutions. (b) A smooth function vanishing at r_s is (r−r_s) times a smooth function (Taylor/Hadamard), and one vanishing to second order is (r−r_s)² times one.
  (c) The decomposition of a vector field into even and odd vector harmonics. (d) The linearized Einstein operator is gauge invariant (PL.G checks it on the family used).
- **Assumed, not derived here.** (i) l ≥ 2 and m = 0 (spherical symmetry). (ii) For the s = 0 solve: time-independent perturbations; ξ static, and ξ = tκ examined (GST.3). (iii) The typed DG06 formulas are those of the arXiv TeX
  (`I146-rpi5/dg06_src`, sha256 4da32b56…), read on 2026-11-23. (iv) Boundary conditions at r = 0 and the choice of ζ, ψ as master variables are not touched. (v) Nothing about M = 0, s = 0 in the e^{st} flow, or r ≤ 0.
- **Not done.** Whether (RW17)'s ζ-equation or (cz) are as printed in the journal version, and in DG08. A second operator for the gauge solve and for c\* (the coordinator's `ind14x.py` style). The general-l certificate of the forward direction (#148).

## 5. Checks

| ids | what | result |
|---|---|---|
| PL.0a–d | #143's: dual inverse; de Sitter; vacuum Schwarzschild, M symbolic; Legendre reduction | PASS |
| PL.L, PL.G, PL.X | two forms of L_ξ g; pure gauge gives δR = 0; dual and covariant Ricci routes agree for the even ansatz | PASS |
| PL.E, PL.P (e^{st} and static), PL.R | extraction round trip; planted pure gauge through extraction and solver; planted flows for the residue/no-log routine | PASS |
| GA.1, GS.1–GS.3, GU.1, GU.2 | structure T_l; the solve; closed forms; divisors {r, r−2M}; uniqueness l = 2..8; the l = 1 family | PASS |
| GST.1–GST.3 | static solve equals s = 0; static uniqueness; no t-linear residual | PASS |
| R.1–R.4 | flow equals #143's; regular point; finite jets; ζ, T, c\* | PASS |
| R.5–R.8 | printed (cz) = −c\*; (RW17); (nv1) relations; (im), (9)∘(3) | PASS |
| R.9–R.12 | Lemma 1 numerators; exponents {0, 1}, no log; M = −1, l = 2, 3, 4; concrete ζ | PASS |
| R.S.1, R.S.2 | static flow equals #145's; same c\* | PASS |

Controls (`zerilli_gauge_rs_controls.py`, 762 s, peak rss 119 MB): the unmutated baseline exits 0 with 34 PASS and no FAIL, and **all 37 mutated runs exit with the expected code (2 at a planted id, 1 at a real one) and the first FAIL id equal to
the target written in the raw README** (`--failfast`). 35 of the 37 targets were written before the first run; `nolog_blind` (PL.R) and `ex_typed` (R.12) were added after the first controls run, whose 35 pairs were all on target (raw README, "Added after the first run").
The three controls the task names: a changed coefficient of (cz) is `cz_coef` (R.5; `cz_sign`, the derived sign, also fires R.5); a dropped term of L_ξ g is `lie_drop` (PL.L); the L = 2 failure is `l2` (GU.1). The others: planted-code `dual_inv` PL.0a, `ds_typo` PL.0b,
`bg_nonvac` PL.0c, `leg_sign` PL.0d, `drop_trace` PL.G, `route2_drop` PL.X, `ext_G` PL.E, `solve_sign` and `solve_drop_s` PL.P, `nolog_blind` PL.R; real `l2_struct` GA.1, `gs_a` GS.1, `hand_typo` GS.2, `gdiv_set` and `gdiv_s` GS.3,
`fam_typo` GU.2, `st_Edep` GST.1, `st_l2` GST.2, `kill_dropG` GST.3, `flow_typo` R.1, `rs_typo` R.2, `pole_inject` R.3, `flow_h` R.4, `rw17_typo` R.6, `nv1_h1` and `nv1_h2` R.7, `im_coef` R.8, `n2_sign` R.9, `log_inject` R.10, `rs6c_typo` R.11,
`ex_typed` R.12, `sflow_typo` R.S.1, `sflow_h` R.S.2. **No miss.** The predictions that rested on structure all held: `drop_trace` (the even gauge h has a nonzero trace), `gs_a` and `pole_inject` (applied after the checks that would otherwise catch them),
`flow_h` (ζ′ stops being H-free), and `im_coef`, which fires at R.8 and not at R.4 because R.4 is invariant under a constant rescaling of ζ (the target was written with that in mind).

What the checks share. The solver and the planted pure-gauge case share the Legendre reduction and `extract`, but the planted h is built by the covariant route and the Lie derivative is cross-checked by PL.L and against the Einstein operator by PL.G.
The flow comes from the library's Ricci code (PL.X checks it against the covariant route, sharing only the background Christoffel symbols and h). c\* is obtained three ways: from the flow (T, R.4), from the singular-point relation of (RW17) (R.6), and from the numerators of the corrected
(nv1), which is DG06's own (3)–(4) (R.7, R.9); the last uses DG06's printed (3)–(4), (9), ζ-ODE, and shares the conventions of (K, H1, H2) and ζ = (2λr−3μ)ψ with the first two. A second operator (a Lichnerowicz form, as the coordinator's checks of #143 and #145) has not been run on the gauge transformation or on c\*.

## 6. Sources

- DG06, gr-qc/0604021, TeX read (`I146-rpi5/dg06_src/`): the section "A modified Regge-Wheeler-Zerilli formalism" (AppB1, RW17, nv1, im, Lemma 1 with cz, AppB2) and (RW2), (RW2a), (poten1), (RW3a) = (3), (4), (6), (9).
  The published (CQG) text was not read, *(unverified)*: the discrepancies above are against this TeX.
- #143, #145 notes and libraries (read; the library imported unchanged); #138 `review_lemma_a.md` §5 (L1.2: χ′/χ(r_s) = −1/(r_s(r_s+2)), L1.3: χ″(r_s) = 0, at M = −1, s = 0, i.e. c\* and the s = 0 value of §1, 4 at M = −1);
  #142 `review_lemma_a3.md` §6 (L6.5, LL11, LL12: Ψ̂ regular at r_s; consistent, not used).
- `phase0/RESULTS.md` (the DG08 Lemma 1 form quoted there, and its check c₀ = −0.533655 at M = −1, l = 2): read for the remark in §1; DG08 not read, *(unverified)*.
- Regge–Wheeler 1957, Zerilli 1970 (publisher HTTP 403 in #143), Martel–Poisson 2005: **not read**; I name no source for the gauge-invariant combinations the solve prints, *(unverified)*.
- No novelty claim is made. The gauge solve is the standard Regge–Wheeler gauge fixing; what this adds to the project is the exact symbolic solve with its divisors, and the r_s analysis with the sign of (cz).

## 7. Lemmas stated for Lean (#135)

Not Lean; statements only. Each is an identity in Q(r, M, s, L) (or a short ODE argument), discharged by `field_simp; ring`, except LG1, LR1 and LR8, which need a geometric or analytic input.

- **LG1** (input) With E = e^{st}, the seven functions of L_ξ g_Schw for ξ = (aPE, bPE, cP′E) are the ΔJ_t, ΔJ_r, ΔG, ΔK, ΔH0, ΔH1, ΔH2 of §1 (the output prints them; needs the Lie derivative of the Schwarzschild metric).
- **LG2** For f ≠ 0, r ≠ 0, the system ΔJ_t = −J_t, ΔJ_r = −J_r, ΔG = −G has the unique solution (a, b, c) of §1 (a short elimination, `field_simp; ring`).
- **LG3** The transformed H0′, H1′, H2′, K′ are the expressions of §1.
- **LG4** For l ≥ 2 the leading coefficient of xP_l′ − (L/2)P_l is −k_l l(l−1)/2 ≠ 0, so P and xP′ are independent and G is a function of h; at l = 1 it is 0 (Legendre polynomials of degrees 0, 1).
- **LG5** At l = 1 the family a = sc/f, b = −f(c′ − 2c/r), c arbitrary solves ΔJ_t = ΔJ_r = 0 and leaves ΔG empty. At s = 0: a = J_t/f; GST.3's determinant 4f/r³.
- **LR1** (input) The flow (K′, H1′, H′) of §3.3, R.1, is the system of the δR_tr, δR_tx, δR_rx equations (coefficients printed), as H0 of #143.
- **LR2** On the flow, ζ = 2r²K − 2(r−2M)H1/s and ζ′ = 2Mr K/(r−2M) + ((L−2)r+4M)H1/(rs) (no H); det T = 2r((L−2)r+6M)/s.
- **LR3** At r_s = −6M/(L−2): T(r_s) has rank 1 and row 2 = c\* row 1, c\* = (L−2)²/(12M(L+1)). (`field_simp; ring` after substituting r_s.)
- **LR4** α, β of §3.3 and (RW17) agree; α_res c\* + β_res = 0 at r_s.
- **LR5** The corrected (nv1) K, H1 equal T⁻¹(ζ, ζ′); their numerators at r_s are nonzero multiples of ζ′ − c\*ζ; N2(r_s) ∝ ζ′ − c\*ζ and N2′(r_s) = −(L−2)/(2M) N2(r_s).
- **LR6** The residue matrix of the (K, H1) flow at r_s has trace 1, determinant 0, and w^T F1 u0 = 0 (no logarithm).
- **LR7** Static: det T_s = −8r²(λr+3M)²/(LM); T_s(r_s) has rank 1 with the same c\*.
- **LR8** (analytic input) a linear ODE system with analytic coefficients has analytic solutions; a smooth function vanishing to order k at r_s is (r−r_s)^k times a smooth function. Hence Lemma 1 with c\*.

## 8. What RESULTS.md and the Zerilli notes should now say

The closing comment of #146 proposed these edits; all are applied. RESULTS.md §4.6 has a #146 bullet and the #143 bullet
points to it, §6's caveat names the gauge and r_s, and the file table has a row. `zerilli_derivation.md` §4 and
`zerilli_static.md` §4 carry a marked coordinator's sentence on each of the gauge and r_s. `notes/literature.md`'s DG06
entry records the three discrepancies, and #149 reads the other versions. `phase0/RESULTS.md` is unchanged: its DG08
Lemma 1 already has the derived sign.

## Coordinator's review

Raw `I146-mac` (negative-mass-raw). Mac, nice 10, at most 2 processes.

| step | result |
|---|---|
| scope | the five files the issue names, nothing else |
| rerun `zerilli_gauge_rs.py` | exit 0, 34 PASS, 0 FAIL; output equal to rpi5's apart from the timing lines |
| rerun the controls | 38 OK, 0 MISS, exit 0; rows equal to rpi5's apart from timings |
| DG06 at the source | the arXiv TeX (`I143-mac`'s download): (cz) prints +2λ²/(3μ(2λ+3)); (nv1), (im), (RW17) as §3 quotes them |
| re-solve of one case | the condition at r_s from DG06's (RW17) alone: ζ′/ζ = −2λ²/(3μ(2λ+3)), independent of s |
| independent check `ind146.py` | run 1 exits 0; all 5 mutations fire on target |

`ind146.py` is #146's recommended task 1 in another form. Its gauge part uses the coordinate Lie derivative
ξ^c∂_c g_ab + g_cb∂_aξ^c + g_ac∂_bξ^c and its own projection onto the seven functions, with P_l explicit. Its r_s part uses
DG06's own (3)–(4), (6), (nv1), (RW17), not the flow. It shares no code with this note's scripts, #143's library or the
coordinator's earlier checks.

| id | what | run 1 |
|---|---|---|
| Q0 | planted: on h = L_η g the solved ξ is −η and h′ = 0, l = 2, 3, 4, e^{st} and static | pass |
| Q1 | on a general h, (a, b, c) of §1 give RW gauge and the four RW functions of §1, l = 2, 3, 4, e^{st} and static | pass |
| Q2 | the seven Δ formulas of §1; at l = 1, T_ab ≡ 0 and the c-family stays free | pass |
| Q3 | (nv1)'s K = (3)'s; (nv1)'s H1 = (2/r) × (3)'s; (3)'s H2 = (nv1)'s with its ζ term reversed | pass |
| Q4 | K's numerator and (RW17) both give c* = −2λ²/(3μ(2λ+3)) at r_s; the printed (cz) is −c*; at μ = 2, c* = −1/(r_s(r_s+2)) | pass |
| Q5 | (RW17) is the Zerilli equation for ζ = (2λr − 3μ)ψ | pass |

So the gauge solve is confirmed, and so are the three discrepancies with DG06's arXiv TeX, from DG06's own equations.
No project result used the printed forms; the project's L1 already has the derived sign. Recommended task 1 is this
review's `ind146.py`, task 2 is #149, and task 3 is row 19 of #135.
