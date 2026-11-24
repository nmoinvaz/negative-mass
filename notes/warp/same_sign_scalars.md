# Can two scalars of one sign source the Alcubierre wall?

November 2026, issue [#78](https://github.com/nmoinvaz/negative-mass/issues/78). Script
[`same_sign_scalars.py`](same_sign_scalars.py), its stdout [`same_sign_scalars_output.txt`](same_sign_scalars_output.txt). It
builds on [`two_scalars.md`](two_scalars.md) (#65), whose "Not covered" names this case, on
[`three_scalars.md`](three_scalars.md) (#77), and on [`phantom_source.md`](phantom_source.md) (#57, #61), whose frame Ricci
tensors and identities (I1), (I2) it imports.

Status words: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails), **numerical
evidence**. #65's lemma has not been independently reviewed (#79). This result does not use it, nor #65's determinant
ideal. It uses #57's identities (I1), (I2) and frame Ricci tensor, reviewed in #59.

## Answer

No (**proved** by hand, with the identities it rests on **checked symbolically**). For every C² profile and any potential
V(φ₁, φ₂), no two minimally coupled scalars of one sign, two normal or two phantom, source Alcubierre's metric on an open
set that meets the wall, the set where f′(r_s) ≠ 0. The same holds for #61's accelerating bubble on every open set that
meets v_s(t) f′(r_s) ≠ 0.

The test fails at its algebraic stage, at λ = 0 and at the λ ≠ 0 branches that #65's weighted trace allows. With (I1) and
(I2), R_ab = λ g_ab ± (a a + b b) holds exactly where R = 0, or λ ≠ 0 and R₃ − λη₃ has rank 1, with R₃ the frame block on
(n, ∂_x, e_ρ) (Lemma S). That forces F = R₀₁(R₀ρ² + R₁ρ²) + 2R_ρρ R₀ρ R₁ρ = 0, a condition with no λ in it. On the wall
F = s² c² B/(8r³), with B an affine function of s² = sin²θ whose constant term, −2w′(u² + v²), does not vanish when
w′ ≠ 0. So at each radius the angles that pass lie in {equator, θ₀, π − θ₀} and no interval of θ does. The accelerating
bubble gives the same, by the c⁰ and c¹ coefficients of B.

At 25 of 25 off-axis exact points the set of real λ with rank(R − λg) ≤ 2 is empty, so no λ of either sign comes close.
The closure stage is never reached on the wall. Its code is tested on planted metrics, two scalars of one sign that pass
all four stages, and on three deformed ones that fail closure at their own sub-test. One of them, the shift
b = C ln ρ + b₀, passes the algebraic stage on an open set with λ = h > 0 and two phantoms, and fails closure because
its kernel plane is not integrable.

## Setup

Two minimally coupled scalars of one sign, with a potential V(φ₁, φ₂), have T_ab = s Σ_I (∂_aφ_I ∂_bφ_I − ½ g_ab (∇φ_I)²)
− g_ab V, with s = +1 (normal) or s = −1 (phantom). G_ab = 8π T_ab is equivalent to

    R_ab = λ g_ab + s (a_a a_b + b_a b_b),    λ = 8πV,    a = √(8π) dφ₁,    b = √(8π) dφ₂.

The test has four stages, as in two_scalars.md.

1. **Algebraic.** Some real λ makes M = R − λg semidefinite (positive for s = +1, negative for s = −1) of rank ≤ 2. Then
   (a, b) is fixed up to Λ(x) ∈ O(2).
2. **Closure.** Some Λ(x) makes a and b closed. For a rank-2 M with kernel plane field 𝒦 this holds iff 𝒦 is integrable,
   L_X M = 0 for X ∈ 𝒦, and the quotient 2-metric of M by 𝒦 is flat (Lemma K, below).
3. **Potential.** K(λ) = 0 for K ∈ 𝒦, so that λ is a function of (φ₁, φ₂). With a, b closed and independent, the
   contracted Bianchi identity gives dλ = (∇·a) a + (∇·b) b (**proved**, the computation of two_scalars.md with the sign
   of b flipped), so stage 3 cannot fail once stage 2 passes. Its control breaks λ directly.
4. **Residual.** G_ab − 8π T_ab = 0 with the fields substituted back.

## Lemma S: the algebraic stage

**Lemma S (proved).** Let R_ab be real and symmetric on an orthonormal frame (e₀, e₁, e₂, e₃) with R_a3 = 0 for every a and
R₀₀ + R₁₁ + 2R₂₂ = 0. Let R₃ be its block on (e₀, e₁, e₂), η₃ = diag(−1, 1, 1) and N = R₃ − λη₃. Then R = λη + s(a a + b b)
for real λ, real covectors a, b and a sign s if and only if R = 0, or λ ≠ 0 and rank N = 1. In the second case s = +1
exactly when λ < 0.

*Proof.* M = R − λη = N ⊕ (−λ). By the identity, N₀₀ + N₁₁ + 2N₂₂ = −2λ. Let s = +1, so that M ⪰ 0 of rank ≤ 2.

- Then −λ ≥ 0 and N ⪰ 0. If λ = 0, N ⪰ 0 has weighted trace 0 with positive weights (1, 1, 2), so N = 0 and R = 0.
- If λ < 0, the e₃ direction is one positive square, so rank N ≤ 1. The trace −2λ > 0 gives N ≠ 0, so rank N = 1.
- Conversely, rank N = 1 gives N = ε n nᵀ with ε wn(n) = −2λ, wn(n) = n₀² + n₁² + 2n₂² > 0. So ε = −sign λ, and M has the
  squares n nᵀ and (−λ) e³e³ᵀ, both of the sign of −λ, so M is semidefinite of rank 2. For λ < 0 that sign is +, for
  λ > 0 it is −. The case s = −1 is the mirror. ∎

So the stage passes at a point exactly where some real λ ≠ 0 makes the 3 × 3 matrix R₃ − λη₃ rank one. The sign of
the pair is then fixed by λ. The weighted trace of two_scalars.md's lemma allowed this branch; the lemma there had
mixed signs and so forced λ = 0.

`lemmaS_samples` tests the statement on 24 random tensors that satisfy the hypotheses, each in a random Lorentz frame: 8 of
the form λη₃ + ε n nᵀ (found, with λ and the pattern the lemma predicts), 8 of #65's mixed rank-2 form a aᵀ − t b bᵀ (not
found: they pass #65's test and fail here) and 8 generic (not found). The prediction comes from exact 2 × 2 minors of the
3 × 3 block, the test from the gcd of the 3 × 3 minors of the 4 × 4 tensor and 60-digit inertia.

## The wall at constant speed

With w = v(1 − f), u = r w″ − w′, v = w u + 2r w′² (this v is not the speed), z = 2w′ + u s², and s = sin θ, c = cos θ,
the frame block is, from phantom_source.py's comoving spherical route (`wall_entries`, exact):

    R₀₀ = −[2w w′ + 2w u c² + r w′² (1 + c²)]/(2r),    2r R₀₁ = −z,    2r R₀ρ = u s c,
    2r R₁ρ = v s c,    R_ρρ = w′² s²/2,    R₁₁ = −R₀₀ − 2R_ρρ.

**Proposition W (proved).** Let w′(r) ≠ 0. The angles θ at which the algebraic stage passes lie in {π/2, θ₀, π − θ₀} for
one θ₀ that depends on r. In particular no interval of θ passes.

*Proof.* Let N = R₃ − λη₃ have rank 1 for some λ. Its 2 × 2 minors vanish, among them m_d = N₀₀N₁ρ − N₀₁N₀ρ and
m_e = N₁₁N₀ρ − N₀₁N₁ρ. Since N₀₀ + N₁₁ = R₀₀ + R₁₁ = −2R_ρρ and the off-diagonal entries of N are those of R₃,

    R₀ρ m_d + R₁ρ m_e = (N₀₀ + N₁₁) R₀ρ R₁ρ − R₀₁ (R₀ρ² + R₁ρ²) = −F.

So F = 0, whatever λ is (`wall_F` checks this certificate and the factorization below as exact polynomial identities).
From the entries,

    F = s² c² B/(8r³),    B = −z (u² + v²) + 2r u v w′² s² = −2w′ (u² + v²) + s² [−u (u² + v²) + 2r u v w′²].

If u² + v² = 0 for real u, v then u = v = 0, and then v = 2r w′² ≠ 0, a contradiction. So the constant term of B does not
vanish and B is a nonzero affine function of s²: it has at most one root s₀². (`wall_ideal`: the two coefficients and
1 − zz w′ r generate the unit ideal, so they have no common zero over ℂ either.) F = 0 therefore needs s = 0, c = 0 or
s² = s₀². On the axis (s = 0), R₀₁ = −w′/r ≠ 0, so the block of N on (e₀, e₁) is [[p, R₀₁], [R₀₁, −p]] with determinant
−p² − R₀₁² < 0. It is indefinite, so M is not semidefinite and the axis does not pass. ∎

An open subset of spacetime that meets the wall contains a box in (r, θ), by axial symmetry and stationarity of the
comoving frame, where f′ ≠ 0 and the form must hold at every point. At each r of the box it would hold on an interval
of θ, which Proposition W excludes. So no (φ₁, φ₂, V) of either sign sources the metric on an open set that meets the
wall. At λ = 0 the form needs R = 0 (Lemma S), and R_ρρ = w′² s²/2 ≠ 0 off the axis.

Exact points, through `same_sign_form`, the function that recovers the planted forms (`wall_points`, `wall_lambda0`): at the
32 points of #77 (Alcubierre's tanh profile with σ = 8, R_b = 1 and v = 1, 1/2, 3; the bump (1 − r²/4)³ with v = 7/4;
random jets; 25 off the axis, 7 on it, one on the equator) no real λ gives R − λg semidefinite of rank ≤ 2. Off the axis
the gcd of the 3 × 3 minors of R − λg has no real root at all. On the axis its only real root is λ = 0, with inertia
(1,1). R has both signs at all 25 off-axis points, 60-digit inertia. For the tanh profile the jets are rounded to
55 digits, so those points are **numerical evidence**, as in #57. Proposition W does not depend on them.

## An accelerating bubble

For #61's b = −v_s(t) f(r_s), with v_s ≠ 0, route D's frame block depends on the jets only through w̃ = x_s′ − v_s f, its
r-derivatives and A = v_s′/v_s (`acc_jets`, as in #65). With w = w̃ the entries are the constant-speed ones with

    R₀₀ → R₀₀ + A w′ c,    R₁₁ → R₁₁ − A w′ c,    R₁ρ = s (v c − A r w′)/(2r)

(`acc_entries`, exact). At A = 0 route D's tensor equals route C's (`acc_const_limit`). With B₀ = A r w′, the same
certificate gives F = s² B_A(c)/(8r³),

    B_A = −z (u² c² + (v c − B₀)²) + 2r w′² u s² c (v c − B₀),

a polynomial of degree 4 in c (`acc_F`). Its extreme low coefficients are c⁰: −(r w″ + w′) B₀² and
c¹: 2B₀ [v (2w′ + u) − r w′² u]. At an instant and radius with A w′ r ≠ 0:

- if c⁰ vanishes then r w″ = −w′, so u = −2w′ and 2w′ + u = 0, and then c¹ = 4 A r² w′⁴ ≠ 0. So B_A is not identically zero
  and has at most four roots in c, hence no interval of θ passes. `acc_ideal` finds the unit ideal for the five
  coefficients with A w′ r ≠ 0.
- at A = 0 (v_s′ = 0) B_A is B and Proposition W applies to w̃, with w̃′ = −v_s f′ ≠ 0 on the wall.
- at v_s = 0 the block is [[−k c, 0, 0], [0, k c, k s/2], [0, k s/2, 0]] with k = v_s′ f′ (`acc_vs0`, exact). The 2 × 2 minors
  of R₃ − λη₃ with k s ≠ 0 generate the unit ideal, so N has rank ≥ 2 off the axis. On the axis λ = k c gives
  rank N = 1, so the stage passes there. The axis has empty interior.

Where v_s f′ ≠ 0 an open set of spacetime contains a slice of fixed t, with fixed A and fixed jets, and the argument above
applies to it. If v_s ≡ 0 on the open set the metric is flat. So the answer carries over (`acc_points`: at 6 exact points
with random jets and v_s′ ≠ 0 there is no admissible root).

## Where the stage passes on an open set: a shift b(ρ)

For a shift b(ρ) alone, R₃ = [[−h, Q, 0], [Q, −h, 0], [0, 0, h]] with h = ½b′² and Q = −½(b″ + b′/ρ) (`brho_rank1`, exact).
Then N = R₃ − λη₃ = [[λ − h, Q, 0], [Q, −λ − h, 0], [0, 0, h − λ]]. Where h ≠ 0 its 2 × 2 minors vanish exactly on
λ = h, Q = 0: the minor −(λ − h)² gives λ = h, and then −Q² gives Q = 0. Q = 0 means (ρ b′)′ = 0, so
b = C ln ρ + b₀ (**proved**). The pair is then two phantoms, since N = diag(0, −2h, 0) and λ = h > 0.

This metric, D3, passes the algebraic stage on an open set. It fails closure. The kernel is spanned by E₀ = ∂_t − b ∂_x and
∂_ρ, and [E₀, ∂_ρ] = b′ ∂_x is not in it, so the plane is not integrable. Since a closed pair a = dφ₁, b = dφ₂ would have
the kernel as its common null plane, there is none (`d3_algebraic`, `d3_closure_fails`). So no pair of phantoms sources
D3. It is the one-sign analogue of #65's b(ρ) example, where closure decides.

## Lemma K: closure of a rank-2 form

**Lemma K (proved).** Let M = s (a ⊗ a + b ⊗ b) be of rank 2 on an open set of a 4-manifold, with kernel plane field 𝒦. A map
Λ(x) ∈ O(2) with (Λ·(a, b)) closed exists if and only if (i) 𝒦 is integrable, (ii) L_X M = 0 for X ∈ 𝒦, and (iii) the
quotient metric q of M by 𝒦 on the local leaf space Q is flat.

*Proof.* If a′ = Λ·(a, b) = dφ with dφ¹ ∧ dφ² ≠ 0, then dφ(X) = 0 for X ∈ 𝒦, so 𝒦 = ker dφ¹ ∩ ker dφ² is integrable and
its leaves are the level sets of φ. M = s (dφ¹dφ¹ + dφ²dφ²) is invariant along them, and it descends to the flat q = s φ*(δ).
Conversely, (i) and (ii) make M descend to a definite 2-metric q on Q. If q is flat, its Cartesian coordinates φ¹, φ² give
M = s Σ dφ^I dφ^I. Two pairs (a, b) and dφ with the same M differ pointwise by a matrix in O(2). ∎

The code normalizes a basis (K₁, K₂) of the kernel to the identity on two coordinate rows (j₁, j₂). Then 𝒦 is integrable iff
[K₁, K₂] = 0, because the (j₁, j₂) components of the commutator are constants. L_{K_i} M = 0 is computed as a tensor. The
block of M on the other two coordinates, at fixed x^{j₁}, x^{j₂}, is a slice of M, and when M is invariant it is the pullback of
q. Its 2D Ricci tensor must vanish (`closure2`).

## Planted and deformed metrics

All go through `same_sign_form` for the algebraic stage and `closure2` for closure. In a product of a 2D Lorentzian space of
constant curvature K₁ with a flat plane, R = K₁ g₁, so λ = K₁ is a root and M = −K₁ g₂: positive for K₁ < 0 (AdS₂, two
normal scalars), negative for K₁ > 0 (dS₂, two phantoms).

| metric | source | λ | signs | what the test finds |
|---|---|---|---|---|
| P1: AdS₂ × R², (−dτ² + dz²)/z² + dR² + R² dΘ² | two normal scalars φ_I = X_I/√(8π), X = R (cos Θ, sin Θ), V = −1/(8π) | −1 | (+,+) | passes all four stages |
| P2: dS₂ × R², (−dτ² + dz²)/τ² + dR² + R² dΘ² | two phantoms, the same fields, V = 1/(8π) | 1 | (−,−) | passes all four stages |
| D1: P1 with R² → R²(1 + R²/10) | none | −1 | (+,+) | integrable, L_X M = 0, quotient not flat |
| D2: −dτ² + dz² + E²(dX² + dY²), E = 1 + (z² − τ²)/10 | none | −4γ/E, γ = 1/10 | (+,+) | integrable, L_X M ≠ 0 and K(λ) ≠ 0 |
| D3: −dt² + (dx + b dt)² + dρ² + ρ² dφ², b = 2 ln ρ + b₀ | none | C²/(2ρ²) | (−,−) | kernel not integrable |

P1 and P2 are written in coordinates (T, Z, R_n, Θ_n) sheared as τ = T + Z/3, z = Z, R = R_n, Θ = Θ_n + T/7 + Z/5, so that
the kernel vectors K₁ = ∂_T − ∂_Θ/7 and K₂ = ∂_Z − ∂_Θ/5 are not coordinate vectors. Both metrics solve G = 8π T exactly
(`p1_residual`, `p2_residual`, with the fields substituted back), and the potential stage gives K_i(λ) = 0. A numerical
factorization of M at 3 exact points, at 60 digits, gives Λ = J Aᵀ(A Aᵀ)⁻¹ with Λ A = J, Λᵀ Λ = 1 and Jᵀ J = s M to better than
10⁻⁴⁵, and Λ differs between the points (**numerical evidence** for the construction; the existence of Λ is Lemma K). For D2,
R − λg = (4γ²)(τ² − z²)(dX² + dY²) exactly, where it is positive, and E is the concircular function of the flat base
with Hessian 2γη. Whether any of these metrics is in print was not searched.

## Checks and controls

Each check has a control that breaks its input or its expectation once. The script runs every control after the checks
and confirms that its target checks fail. All 31 checks pass and all 22 controls fire (`same_sign_scalars_output.txt`).

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `planted_pass` | 12 random λη + s(a a + b b), s = ±1, transformed, are recovered | `planted-sign` | λη + a a − b b |
| `planted_fail` | 12 forms with mixed signs or three or four squares: no admissible root | `planted-fail-sign` | two squares of one sign |
| `lemmaS_samples` | Lemma S on 24 tensors with (I1), (I2) | `lemmaS-weight` | R₀₀ − R₁₁ + 2R₂₂ = 0 in place of (I2) |
| `p1_algebraic`, `p2_algebraic` | λ a root with inertia (2,0), (0,2) at 3 points | `p1-mixed`, `p2-mixed` | the plane's signature (−,+) |
| `p1_closed` | kernel integrable, L_X M = 0, quotient flat, Λ(x) | `p1-fiber` | the plane curved |
| `p1_potential`, `p2_potential` | K_i(λ) = 0 | `p1-lambda`, `p2-lambda` | λ + T/10 |
| `p1_residual`, `p2_residual` | G − 8π T = 0 | `p1-sign`, `p2-sign` | one kinetic sign flipped in T |
| `d1_closure_fails` | D1: quotient not flat | `d1-flat` | curvature 0 |
| `d2_algebraic`, `d2_closure_fails` | D2: root, L_X M ≠ 0, K(λ) ≠ 0 | `d2-gamma0` | γ = 0, so M = 0 |
| `d3_closure_fails` | D3: kernel not integrable | `d3-bracket` | the commutator without K₂(K₁^μ) |
| `wall_identities`, `wall_entries`, `wall_F` | R(e_ϕ, ·) = 0, (I2), the entries, F and the certificate | `alc-drop-w2` | the shift's square dropped from g_tt |
| `wall_F` | the F formula | `wall-F-sign` | +z in place of −z |
| `wall_ideal` | the coefficients of B and 1 − zz w′ r: unit ideal | `ideal-nosat` | without the saturation |
| `wall_points` | 32 exact points: no admissible root | `wall-planted` | R replaced by a planted λη + a a + b b |
| `wall_lambda0` | R has both signs at 25 off-axis points | `wall-psd` | R replaced by a positive definite tensor |
| `acc_entries` | the bubble's entries | `acc-drop-A`, `acc-static` | the A terms dropped; no t-derivatives |
| `acc_jets`, `acc_const_limit`, `acc_F`, `acc_vs0`, `acc_points` | the bubble's reduction, limit, F, v_s = 0 and 6 points | | |
| `acc_ideal` | the five coefficients with A w′ r ≠ 0: unit ideal | `acc-drop-odd` | only the even powers of c |
| `brho_rank1` | b(ρ): the zero set is λ = h, Q = 0 | `brho-q` | the claimed set with Q = h |
| `d3_algebraic` | D3: root λ = h, inertia (0,2) | | |

`acc-static` also fails `acc_jets`, `acc_const_limit` and `acc_vs0`, and `wall-planted` also fails `wall_lambda0`: the
README predicted one target each and nothing else needed to fail. Both extra failures are consequences of the mutation
(the jets stop depending on the t-derivatives; the planted tensor has both signs), so they are reported, not adjusted.
`lemmaS-weight` is strong: 8 of the 24 samples disagree with the lemma in that control.

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/same_sign_scalars.py

It exits 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, `--mutate=NAME` one control
and `--stop=N` the sections up to N. On rpi5 (Raspberry Pi 5, one process, nice 10) the checks took 26 s at 72 MB peak resident and the whole run with its 22 controls 351 s. Two hash seeds (`PYTHONHASHSEED` 1 and 42) give the same output apart from the timings.

### Discrepancies

- **C8 as coded: FAIL.** The first `brho_rank1` tested equality of ideals, ⟨minors, 1 − zz h⟩ = ⟨λ − h, Q, 1 − zz h⟩
  (`pilot_05_full_nocontrols.log` in `negative-mass-raw/I78-rpi5`). The inclusion of the minors in ⟨λ − h, Q⟩ holds, the
  converse does not: the Groebner basis with h ≠ 0 is {λ − h, Q²}, not radical. The criterion C8 in the README speaks of
  the solutions, so the code tested more than the README says. The coordinator had it recoded as worded: (a) every minor
  reduces to 0 modulo ⟨λ − h, Q⟩, and (b) λ − h and Q are in the radical of the saturated ideal by Rabinowitsch, with no
  exponent chosen. The basis is printed and does not gate. The recoded check passes, and `brho-q` makes (a) and (b) fail.
- The first planted geometry sheared R as well, and `closure2` on D1 had run 2 min 52 s when it was stopped by hand
  (`pilot_02_stop3.log`). The final shear leaves R alone and takes 6 s. No criterion changed.
- Two controls fire more checks than predicted, as above.

## What the code proves, and what is assumed

The code proves, as exact identities or exact ideal computations: (I1), (I2) and the entries of R₃ for the wall in terms of
u, v, z; the certificate and the factorization of F; that the coefficients of B and w′ r generate the unit ideal; the
accelerating bubble's entries, the polynomial B_A, its c⁰ and c¹ coefficients, the unit ideal of its five coefficients,
and the block at v_s = 0; that R₃ of a shift b(ρ) has the stated form and that the zero set of its rank-one minors is
λ = h, Q = 0; that the planted metrics solve the field equations and pass closure and potential, and that the deformed
ones fail at their sub-tests. The inertia of R − λg comes from 60-digit eigenvalues, so the semidefiniteness tests are
numerical evidence. The gcd of the 3 × 3 minors, and so the absence of a real root at the 25 off-axis points, is exact.

By hand: Lemma S, Proposition W and the step from it to open sets, the accelerating argument, the solution of Q = 0, the
integrability of D3's kernel and Lemma K.

Trusted: sympy 1.14.0 and mpmath 1.3.0; `ricci`, `sph_tetrad_ricci`, `lab_tetrad_ricci` and route D of phantom_source.py
(tested there and reviewed in #59); `inertia`, `lie_K`, `wall_frame` and the point lists of three_scalars.py. `same_sign_form`
and `closure2` are new, tested on the planted forms and metrics. Not independent: the wall's entries and the exact points
come from the same frame Ricci tensors, and the point lists are #77's.

Assumed: two real, minimally coupled scalars of one sign with a constant diagonal kinetic form, any potential V,
G_ab = 8π T_ab, a profile with f(r_s) C² on ℝ³ (then w ∈ C² and R₃ exists), unit lapse and flat slices for the shift, and
the local leaf space of 𝒦 regular where Lemma K is used.

## Lemmas stated for Lean (#135)

1. **(S)** Lemma S above, for real symmetric R with R_a3 = 0 and R₀₀ + R₁₁ + 2R₂₂ = 0. Status: proved.
2. **(F)** For a real symmetric 3 × 3 N = R₃ − λη₃ with R₀₀ + R₁₁ + 2R_ρρ = 0: R₀ρ m_d + R₁ρ m_e = −F, and rank N ≤ 1 implies
   F = 0. Status: proved (a polynomial identity).
3. **(W)** With w′ ≠ 0, the constant term −2w′(u² + v²) of B is nonzero. Status: proved.
4. **(A)** With A w′ r ≠ 0 the c⁰ and c¹ coefficients of B_A do not vanish together. Status: proved.
5. **(K)** Lemma K above, with the local leaf space regular. Status: proved.
6. **(E)** The entries of R₃ for the wall and for the bubble. Status: checked symbolically (#57, #59, this note).

## Not covered

- A non-axisymmetric shift and a general shift b(t, x, ρ): Lemma S and the certificate (F) apply to them, but #80's
  classification of the shifts that pass on an open set is not done here.
- The isolated angles that pass F = 0, the equator and θ₀. The stage may pass on a curve of the wall for special profiles;
  closure is a statement on a neighborhood, so it is not asked there. No such profile was found, and 40 random trials of
  a numerical search found no passing point at all (`explore_iso0` in the raw folder, not a claim).
- Three or more scalars of mixed signs: #77. Three or more of one sign, and four or more scalars, are not covered.
- Non-minimal coupling, kinetic terms that are not quadratic, and scalars with other matter.
- Bubbles whose path is not a straight line.

## What `two_scalars.md`'s "Not covered" should now say

Replace the bullet "Two scalars of one sign (#78). The lemma's argument allows λ ≠ 0 there." by:

> Two scalars of one sign: [`same_sign_scalars.md`](same_sign_scalars.md) (#78). With (I1), (I2) the algebraic stage holds iff
> R = 0 or λ ≠ 0 and rank(R₃ − λη₃) = 1, with pair (+,+) iff λ < 0. It fails on every open set that meets the wall, for
> constant speed and for the accelerating bubble (**proved** by hand, with the identities **checked symbolically**), so
> closure is never reached. The shift b = C ln ρ + b₀ passes the algebraic stage on an open set with two phantoms and
> fails closure (non-integrable kernel). Shifts b(t, x, ρ) other than the bubble and b(ρ) are not classified.

## Prior work

The record of #60 ([prior_work.md](prior_work.md) §1) found no paper that asks whether a scalar field sources the
Alcubierre metric. No further search was made for two scalars of one sign. [BD20] fit two- and three-fluid mixtures to the
Alcubierre stress tensor numerically ([prior_work.md](prior_work.md) §2, not re-read here); a minimally coupled scalar is a
perfect fluid only where its gradient is timelike, so their mixtures and these pairs overlap without being the same class.

## References

Entries as in [`../literature.md`](../literature.md), "Warp drives".

- [BD20] Béatrix-Drouhet, arXiv:2012.09941 (2020, preprint).
