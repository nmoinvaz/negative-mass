# Independent review of #57 and #58's Alcubierre results

November 2026, issue [#59](https://github.com/nmoinvaz/negative-mass/issues/59). The reviewer had not seen the work.
Scripts written for the review: [`review_cartan.py`](review_cartan.py) with
[`review_cartan_output.txt`](review_cartan_output.txt), and [`review_spinor.py`](review_spinor.py) with
[`review_spinor_output.txt`](review_spinor_output.txt). Both run from the repository root, in the environment of
[`setup.md`](../../setup.md),

    python3 notes/warp/review_cartan.py
    python3 notes/warp/review_spinor.py

and exit 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, `--mutate=NAME` one
control. They take about 15 s and 40 s.

## Verdicts

No error was found. Both results stand. The gaps below are about conventions, wording and scope, and none of them
changes a verdict. With section 1 and the lemma, #57's verdict now has a proof by hand as well as its symbolic check,
and #58's item 5 has been reviewed. Any change of status word is for the coordinator.

| result | step | verdict |
|---|---|---|
| #57 | metric, scalar field equation R_ab = 8π(ε ∂_aφ ∂_bφ + V g_ab) | confirmed, by hand |
| #57 | Eulerian energy density ρ_E = −(v²/32π)(ρ²/r_s²) f′² | confirmed, by hand and `alc_wall` |
| #57 | every row of the frame Ricci table | confirmed, `alc_table`, and R_00, R_11, R_ρρ by hand |
| #57 | (I1) R(e_ϕ, ·) = 0 | confirmed, by hand and `axi_I1`; holds more generally (gap 6) |
| #57 | (I2) R_00 + R_11 + 2R_ρρ = 0 | confirmed, by hand and `axi_I2`; holds more generally (gap 6) |
| #57 | the identities on the axis | confirmed, by hand and `axi_axis` |
| #57 | wall curvature, R_ρρ off the axis and R_01 = v f′/r_s on it | confirmed, by hand and `alc_wall` |
| #57 | the lemma | confirmed, by hand |
| #57 | its application, off and on the axis, and the verdict | confirmed; the centre r_s = 0 needs a word (gap 5) |
| #57 | perfect-fluid corollary (algebra only) | confirmed; literature comparison still open (#60) |
| #57 | complex pair −R_00 ± i R_01 on the axis, 16 ± 4i | confirmed, by hand to the printed digits |
| #57 | Ellis–Bronnikov control, R_ll = −2a²/(l² + a²)², φ = −arctan(l/a)/(2√π) | confirmed, by hand and `planted_eb` |
| #57 | control `--mutate=eb-eps+1` | rerun, fired (exit 1, `eb_sign` fails) |
| #58 | item 1, μ and J_i = (v/16π)(∂_x∂_i f − δ_ix ∇²f) | confirmed, by hand; used by `weitzenbock` |
| #58 | item 2, E = P = 0 | confirmed; P = 0 holds exactly, not only numerically (gap 7) |
| #58 | item 3, the minimal split and μ₋ + \|J₋\| = max(0, \|J\| − μ) | confirmed, by hand |
| #58 | item 4, thin-wall limits −v²R²σ/36 and πvR²σ/32 | confirmed, by hand; the table was not recomputed |
| #58 | item 5, ∇̃ on a flat slice and D̃ = D + ½(tr K)γ⁰ | confirmed for the stated ∇̃, whose sign is a convention (gap 1) |
| #58 | item 5, D: Ḣ¹ → L² an isometry onto | confirmed, by hand and `isometry` |
| #58 | item 5, Hölder, Kato and the sharp Sobolev inequality | confirmed, by hand and `sobolev_const` |
| #58 | item 5, contraction constant κ_K = ‖tr K‖_{L³}/(2√S) | confirmed, by hand; attained in the one-step bound (`contraction_sharp`) |
| #58 | item 5, existence, uniqueness and no decaying kernel | confirmed, in the class ε − ε₀ ∈ Ḣ¹ (gap 2) |
| #58 | item 5, v_c = 3.745, 2.471, 1.562, 0.985 and the thin-wall formula | confirmed, `v_c` and by hand |
| #58 | item 6, the identity and its sign convention | confirmed, by hand and `weitzenbock` |
| #58 | item 6 against phase4/RESULTS.md §1 | consistent; phase4 §1 leaves its conventions open (gap 4) |
| #58 | control `--part=spinor --mutate=norm` | rerun, fired (exit 1, identity off by 0.5) |

## 1. #57: (I1) and (I2) by Cartan's structure equations

This derivation is by hand. `review_cartan.py` repeats it in sympy with its own exterior calculus and checks every
formula below.

Write the metric as ds² = −dt² + (dx + b dt)² + dy² + dz², with b = β^x = −v f(r_s). The coframe θ⁰ = dt,
θ¹ = dx + b dt, θ² = dy, θ³ = dz is orthonormal, and its dual frame is e_0 = ∂_t − b ∂_x = n and e_i = ∂_i. With
b_j = ∂_j b, the only nonzero exterior derivative is dθ¹ = db ∧ dt = −b_j θ⁰ ∧ θʲ.

Let K_ij = ½(δ_jx b_i + δ_ix b_j), which is ½(∂_iβ_j + ∂_jβ_i), and A_ij = ½(δ_ix b_j − δ_jx b_i). The connection
forms

    ω⁰_i = ωⁱ_0 = −K_ij θʲ,    ωⁱ_j = A_ij θ⁰

are antisymmetric once lowered, and they solve dθ^a = −ω^a_b ∧ θ^b, because K_ij + A_ij = δ_ix b_j. They are therefore
the Levi-Civita connection. The curvature Ω^a_b = dω^a_b + ω^a_c ∧ ω^c_b, contracted to Ricci, gives

    R_00 = n(K) − K_ij K_ij,    R_0i = ∂_i K − ∂_j K_ji,
    R_ij = −n(K_ij) + K_jx b_i + K_ix b_j + K K_ij − 2 K_ik K_kj,    K = tr K = b_x,

which are the Gauss–Codazzi–Ricci equations for unit lapse and flat slices (`general_3p1`). With A, B ∈ {y, z} and
∇_⊥ = (∂_y, ∂_z), they evaluate to

    R_00 = n(b_x) − b_x² − ½|∇_⊥b|²,    R_xx = −n(b_x) + b_x² − ½|∇_⊥b|²,    R_AB = ½ b_A b_B,
    R_0x = −½ Δ_⊥ b,    R_0A = ½ ∂_x b_A,    R_xA = −½ n(b_A) + b_x b_A.

So, for every b(t, x, y, z) and with no symmetry assumed,

    R_00 + R_xx = −|∇_⊥b|²,    R_yy + R_zz = ½|∇_⊥b|²,    so    R_00 + R_xx + 2(R_yy + R_zz) = 0,

and G_00 = ½(R_00 + R_xx + R_yy + R_zz) = −¼|∇_⊥b|², which is #57's ρ_E (`general_I2`).

Now let b depend on (y, z) only through ρ. The components of e_ϕ depend only on (y, z), so they commute with ∂_x and
n. Every component R(e_ϕ, ·) above is a multiple of e_ϕ(b) or of its ∂_x or n derivative. R_ϕϕ = ½ e_ϕ(b)²,
R_ρϕ = ½ e_ρ(b) e_ϕ(b), R_0ϕ = ½ ∂_x e_ϕ(b) and R_xϕ = −½ n(e_ϕ(b)) + b_x e_ϕ(b), and e_ϕ(b) = 0. That is (I1). It
also gives R_ϕϕ = 0, so R_ρρ = R_yy + R_zz, and the displayed identity becomes (I2).

Off the axis, R_ρρ = ½(∂_ρ b)² = v² f′² ρ²/(2r_s²). On the axis at r_s > 0, Δ_⊥ f(r_s) = 2f′/r_s, so
R_01 = v f′/r_s. Both match the note (`alc_wall`). Every row of the note's table, written in w = v(1 − f), equals the
Cartan result exactly (`alc_table`). R_00, R_11 and R_ρρ were also matched by hand, using n(g) = −w ∂_x g for any
g(x − vt, y, z).

The routine is first run on planted truth. de Sitter gives R_ab = 3H²η_ab, Schwarzschild gives R_ab = 0 and
Ellis–Bronnikov gives R_ll = −2a²/(l² + a²)². The structure equations themselves are checked (`structure_eqs`), and
so is ω_0j(e_i) = K_ij (`connection_K`), which section 3 uses.

## 2. #57: the lemma and its application

The lemma was checked line by line. If s = 0 or u = 0 then R = λη, and R_33 = λ = 0. Otherwise s u_a u_3 = 0 for
a ≠ 3 and λ + s u_3² = 0. If u_3 ≠ 0 then u = u_3 e³, λ = −s u_3² ≠ 0, and R_00 + R_11 + 2R_22 = −λ + λ + 2λ = 2λ ≠ 0,
which is a contradiction. If u_3 = 0 then λ = 0 and s(u_0² + u_1² + 2u_2²) = 0 forces u = 0. The reality of u is
used, and it holds for u_a = e_a(φ) with φ real. A real symmetric R − λη of rank ≤ 1 is s u⊗u with u real, so "no
real λ gives rank ≤ 1" is the same statement.

The application is right. On the axis, E_2 and E_3 satisfy the lemma's hypotheses, since R(E_A, ·) = 0 and
R_00 + R_11 = 0 there. The curvature is nonzero wherever f′(r_s) ≠ 0 and r_s > 0. The perfect-fluid corollary
follows, since R_ab = 8π[(ρ + p) u_a u_b + ½(ρ − p) g_ab] has the excluded form. On the axis, R^a_b restricted to
(E_0, E_1) is [[−R_00, −R_01], [R_01, −R_00]]. For σ = 8, R = 1, v = 1 at x − vt = 1, f′ ≈ −4, f ≈ ½ and f″ ≈ 0
give R_00 ≈ −16 and R_01 ≈ −4, so the pair is 16 ± 4i.

## 3. #58: item 5, the existence lemma

**The reduction.** With the note's gammas, {γ^a, γ^b} = 2η^ab, γ⁰ is anti-Hermitian and unitary, and the γ^j are
Hermitian (`clifford`). For symmetric K, K_ij γ^i γ^j = ½ K_ij {γ^i, γ^j} = tr K. So the stated
∇̃_i = ∂_i + ½ K_ij γ^j γ⁰ gives D̃ = γ^i ∇̃_i = D + ½ (tr K) γ⁰.

The sign of the K term is a convention, which the note does not say. The spin connection that keeps Clifford
multiplication parallel is Γ = ¼ ω_ab γ^a γ^b (`spin_restriction`). With ω_0j(e_i) = K_ij from section 1, it restricts
to the slice as ∂_i − ½ K_ij γ^j γ⁰ in these gammas. The note's ∇̃ is that connection after γ⁰ ↦ −γ⁰, an automorphism
of the Clifford algebra. Both signs give a valid identity (section 4). The existence argument uses only |tr K| and the
unitarity of γ⁰, so it holds for either sign.

**The isometry.** The symbol of D is iγ^j k_j, with (γ^j k_j)² = |k|² and γ^j Hermitian, so |D̂ψ| = |k| |ψ̂|
pointwise in k. Plancherel gives ‖Dψ‖₂ = ‖∇ψ‖₂. The inverse symbol −iγ^j k_j/|k|² maps L² onto Ḣ¹, which in three
dimensions is a space of functions contained in L⁶. `isometry` checks ‖Dψ‖₂ = ‖∇ψ‖₂ to 2×10⁻¹⁶ on a grid with
analytic derivatives.

**The contraction.** With ε = ε₀ + ψ and Dε₀ = 0, D̃ε = 0 is ψ = T(ψ) := −D⁻¹(½ tr K γ⁰(ε₀ + ψ)). T maps Ḣ¹ into Ḣ¹
when tr K ∈ L² (for the ε₀ term) and tr K ∈ L³ (for the ψ term, by Hölder and Ḣ¹ ⊂ L⁶). For h = ψ₁ − ψ₂,

    ‖∇(Tψ₁ − Tψ₂)‖₂ = ½‖tr K γ⁰ h‖₂ = ½‖tr K |h|‖₂ ≤ ½‖tr K‖₃ ‖h‖₆ ≤ ½‖tr K‖₃ ‖∇|h|‖₂/√S ≤ κ_K ‖∇h‖₂.

The steps are the isometry, the unitarity of γ⁰, Hölder (1/2 = 1/3 + 1/6), the Sobolev inequality for the scalar |h|
and Kato's inequality. Kato's inequality, |∇|h|| ≤ |∇h| almost everywhere, follows from ∂_j|h| = Re⟨h, ∂_j h⟩/|h|
where h ≠ 0 and Cauchy–Schwarz, with ∇|h| = 0 almost everywhere on {h = 0}. Banach's fixed-point theorem gives a unique
ψ ∈ Ḣ¹ when κ_K < 1. With ε₀ = 0 the map is linear with norm at most κ_K, so D̃ has no kernel in Ḣ¹.

**The constant.** S = 3(π/2)^{4/3} = ¾|S³|^{2/3} with |S³| = 2π² is the sharp constant in three dimensions, attained
by (1 + r²)^{−1/2} [Aubin 1976, Talenti 1976] *(unverified)*. `sobolev_const` computes the quotient
‖∇u‖₂²/‖u‖₆² for that function as 5.477904089531, equal to S to 10⁻¹², while three other trial functions give
6.06–8.16. In the one-step bound, Hölder is an equality for |tr K| ∝ |ψ|² and Sobolev for the bubble, so
½‖tr K ψ‖₂ = κ_K ‖∇ψ‖₂ for tr K ∝ u² and ψ = uχ (`contraction_sharp`, equal to 12 digits). The constant κ_K
cannot be improved within this argument. The threshold κ_K < 1 stays only sufficient, as the note says.

**v_c.** tr K = ∂_xβ_x = −v f′(r) cos θ, and ∫|cos θ|³ dΩ = π, so ‖tr K‖₃³ = π v³ ∫|f′|³ r² dr. Radial quadrature
gives v_c = 3.74490, 2.47097, 1.56243 and 0.98503 for σR = 2, 4, 8 and 16 (`v_c`), equal to the note's values to
their printed digits. For thin walls f′ ≈ −(σ/2) sech² σ(r − R) and ∫sech⁶ = 16/15, so ∫|f′|³ r² dr → 2R²σ²/15. That
gives the note's v_c ≈ 2√S (15/(2πR²σ²))^{1/3}, which is 0.98528 at σR = 16.

## 4. #58: item 6, the identity and its conventions

For any ε = ε₀ + ψ with ψ decaying, expanding |∇̃ε|² and |D̃ε|² gives three kinds of term.
- The algebraic terms are Σ_i |½ K_ij γ^j γ⁰ ε|² − |½ tr K γ⁰ ε|² = ¼(K_ij K_ij − (tr K)²)|ε|² = −4π μ |ε|², by the
  Hamiltonian constraint.
- The cross terms are Re⟨∂_iε, P_ij γ^j γ⁰ ε⟩ with P_ij = K_ij − δ_ij tr K and γ^jγ⁰ Hermitian. That is a divergence
  minus ½(∂_iP_ij)⟨ε, γ^jγ⁰ε⟩. The momentum constraint gives ∂_iP_ij = 8πJ_j, and ⟨ε, γ^jγ⁰ε⟩ = −ξ^j, so the cross
  terms give +4π J·ξ.
- What remains is ∫|∂ψ|² − ∫|Dψ|² = 0.

So ∫|∇̃ε|² − ∫|D̃ε|² + 4π∫(μ|ε|² − J·ξ) dV is the boundary term, 4π(E − P·ξ₀). With D̃ε = 0 this is item 6, sign for
sign. On these data E = 0 and the boundary integrals of K vanish, since K decays exponentially.

`weitzenbock` checks the integrated formula directly, with no Dirac solve and no FFT. The test spinor is ε₀ plus two
Gaussian wave packets, with item 1's closed forms for μ and J and analytic derivatives on cell-centred grids. The
residual is 4.8×10⁻¹³ (σR = 2, v = 1) and 2.5×10⁻¹² (σR = 4, v = 0.5) relative to ∫|∇̃ε|². The packets make
∫J·ξ = 2.6×10⁻² and 7.7×10⁻³, so the sign of the J term is tested, and the control `Jsign` fires. With −K in ∇̃ and D̃,
the same holds with +J·ξ (`weitzenbock_flip`). The sign of J·ξ thus follows the sign convention of ∇̃.

**Against phase4/RESULTS.md §1.** The future causal vector in the identity is X = |ε|² n + ξ^j e_j, with
|ξ| ≤ |ε|², because the γ⁰γ^j are Hermitian and anticommute (`clifford`). For the standard J_i = −n^a T_ai of item 1,
T(X, n) = μ|ε|² − J·ξ, which is item 6's integrand. phase4 §1 writes T(ξ_ε, n) = μ|ε|² + J·ξ⃗. That is the same
quantity with ξ⃗ → −ξ⃗, as the note says, or with J_i = T(e_i, n). The two forms are the two signs of ∇̃ above. Lemma A
uses only T₋(ξ_ε, n) ≥ 0 and |J₋·ξ| ≤ |J₋||ε|², and Lemma B has K = 0 and J = 0, so neither depends on the choice.
The note's statement is confirmed.

## 5. Reruns and controls

| command | exit | result |
|---|---|---|
| `phantom_source.py --mutate=eb-eps+1` | 1 | fired, `eb_sign` fails, 9.5 s, 75 MB |
| `phantom_source.py --no-controls` | 0 | 13 checks pass, output as stored |
| `alcubierre_bounds.py --part=spinor --mutate=norm` | 1 | fired, identity off by 0.5, 48 s, 3.43 GB |
| `alcubierre_bounds.py --part=spinor` | 0 | identity to 4.5×10⁻¹⁴, 4.4×10⁻¹⁴, 1.0×10⁻¹², as stored, 52 s, 3.56 GB |
| `review_cartan.py` | 0 | 12 checks pass, 7 controls fire |
| `review_spinor.py` | 0 | 9 checks pass, 5 controls fire, 2.6 GB |

| review control | change | checks that fail |
|---|---|---|
| `rotate` | (e_ρ, e_ϕ) rotated by α, tan(α/2) = 1/20 | `axi_I1`, `axi_I2` |
| `nonaxi` | shift a function of y² + 2z² | `axi_I1`, `axi_I2` |
| `lapse` | lapse 1 + b | `connection_K`, `general_3p1`, `general_I2`, `axi_I1`, `axi_I2`, `axi_axis` |
| `cartan-sign` | sign of ω ∧ ω flipped | planted truth, `general_*`, `axi_I2`, `alc_*` |
| `torsion` | wrong solution for ω | `structure_eqs` and every later check |
| `table-sign` | sign of the note's R_0ρ flipped | `alc_table` |
| `wall-v2` | 2v in the expected wall curvature | `alc_wall` |
| `sobolev` | 21/20 of S | `sobolev_const`, `contraction_sharp` |
| `gamma` | γ¹ scaled by 11/10 | `clifford`, `spin_restriction`, `isometry`, `weitzenbock`, `weitzenbock_flip` |
| `Jsign` | +J·ξ in the identity | `weitzenbock`, `weitzenbock_flip` |
| `musign` | −μ in the identity | `weitzenbock`, `weitzenbock_flip` |
| `vc-angle` | 2π for ∫\|cos θ\|³ dΩ | `v_c` |

## 6. Gaps

None of these changes a verdict. They are reported, not fixed, in the notes under review.

1. **#58 item 5, the sign of ∇̃.** In the note's own gammas, the spacetime spin connection restricts to the slice as
   ∂_i − ½K_ij γ^j γ⁰, not ∂_i + ½K_ij γ^j γ⁰. The note's choice is that connection after γ⁰ ↦ −γ⁰, and its
   identity is correct as written (section 4). The note should say which convention it uses.
2. **#58 item 5, the solution class.** Uniqueness and "no decaying kernel" hold in ε − ε₀ ∈ Ḣ¹, which the note does
   not state. The identity's boundary terms need more decay than Ḣ¹ gives in general. For these data
   ψ = −D⁻¹(source) with an exponentially decaying source, so ψ = O(r⁻²) and the terms vanish.
3. **#58 item 5, Natário's warp drive** is named without a reference. Its tr K = 0 is recalled, not read here.
4. **phase4/RESULTS.md §1** writes ξ_ε^a = ⟨ε, γ⁰γ^a ε⟩ with ξ⁰ = |ε|². In #58's gammas ⟨ε, γ⁰γ⁰ε⟩ = −|ε|², so the
   two files differ in more than the sign of ξ⃗. phase4 §1 fixes neither its gammas nor the sign of J. Nothing in
   Lemmas A or B depends on either.
5. **#57, the bubble centre.** "C² profile" should mean that f(r_s) is C² on ℝ³, which forces f′(0) = 0. Otherwise the
   metric is not C² at r_s = 0, and the table's s and c are undefined there. Both profiles used are even in r.
6. **#57, scope (a strengthening).** (I1) and (I2) hold for every shift b(t, x, ρ) that is axisymmetric about the
   x-axis, with unit lapse and flat slices (section 1). That includes Alcubierre's metric with a time-dependent
   v_s(t), which the note lists as not covered. There R_ρρ = ½(∂_ρ b)² and, on the axis, R_01 = v_s(t) f′/r_s, so the
   obstruction carries over.
7. **#58 item 2, status.** E = 0 on a flat slice, and P = 0 exactly, because K decays exponentially and the surface
   integral's limit is 0. The note calls both numerical evidence.
8. **Memory.** `alcubierre_bounds.py --part=spinor` peaked at 3.43–3.56 GB resident on this Mac, not under 2 GB as
   #58's progress post says.
9. **Not checked here.** These were not checked: the claim that the arXiv version of Alcubierre's Eq. (19) prints r_c,
   item 4's table, the 46 exact-point tests and the remaining controls of both scripts.

## What the checks share with the work under review

`review_cartan.py` shares only sympy 1.14.0 with `phantom_source.py` and `alcubierre_bounds.py`. Its jets, exterior
calculus, connection solve and curvature are its own. `review_spinor.py` shares numpy, scipy and sympy. It uses item 1's
closed forms of μ and J and the gammas of the note's "Conventions", both of which are under review. A wrong μ or J
would break `weitzenbock`, which needs both constraints. It solves no Dirac equation.

What the code proves: `review_cartan.py` proves the displayed Ricci formulas, (I1), (I2), the table and the wall
curvature as exact identities in the jets of b and the point. Each expression is zero after `sp.expand`, or after
`sp.simplify` where r_s appears as a square root. `review_spinor.py` gives numerical evidence for the Sobolev constant,
v_c and the Weitzenböck identity, and exact checks of the Clifford relations. The lemma of #57 and the existence
argument of #58 are proved by hand here, and the code proves neither.

## References

- T. Aubin, "Problèmes isopérimétriques et espaces de Sobolev", J. Differential Geom. 11, 573 (1976) *(unverified)*.
- G. Talenti, "Best constant in Sobolev inequality", Ann. Mat. Pura Appl. 110, 353 (1976) *(unverified)*.
