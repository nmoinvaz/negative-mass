# Independent review of #64's Hawking–Ellis type map

November 2026, issue [#71](https://github.com/nmoinvaz/negative-mass/issues/71). Reviews [he_type.md](he_type.md).
Scripts [`review_he_derive.py`](review_he_derive.py) with [`review_he_derive_output.txt`](review_he_derive_output.txt),
[`review_he_compare.py`](review_he_compare.py) with [`review_he_compare_output.txt`](review_he_compare_output.txt),
[`review_he_d0.py`](review_he_d0.py) with [`review_he_d0_output.txt`](review_he_d0_output.txt), and the rerun of
`he_type.py` in [`review_he_rerun_output.txt`](review_he_rerun_output.txt). Each runs from the repository root, in the
environment of [`setup.md`](../../setup.md), and exits 1 if a check fails or a control does not fire,

    python3 notes/warp/review_he_derive.py     # 13 checks in 8 s, then 13 controls, 1.5 min in all
    python3 notes/warp/review_he_compare.py    # 13 comparisons, then 3 controls, 13 s
    python3 notes/warp/review_he_d0.py         # 3 checks, then 2 controls, 11 s

`--no-controls` runs the checks alone and `--mutate=NAME` one control.

Status words are those of AGENTS.md. Nothing here claims novelty.

## Verdicts

No error was found. All four identities reproduce exactly, by two routes that share no code with `he_type.py` or
`phantom_source.py`, and `he_type.py` reproduces its stored output. Six gaps remain, all of wording or scope, and none
changes a verdict. Any change of status word is for the coordinator.

| step | verdict | how |
|---|---|---|
| T1. D = v⁶(p₆ + p₈v² + p₁₀v⁴ + p₁₂v⁶), no other powers of v | confirmed | **proved** by hand (section 3), `T1_structure` |
| T1. p₆v⁶ = −4\|j\|⁶ and p₆ = −J⁶/(16r⁶) | confirmed, "< 0" needs j ≠ 0 (gap 1) | **proved** by hand, `T1_structure` |
| T2. p₁₂ ≥ 0 | confirmed | **proved** by hand, `T2_p12` |
| T2. p₁₂ = s⁴F²Q₁²Q₂/(16r⁶), Q₁ and Q₂ as printed | confirmed | **checked symbolically**, `T2_p12` |
| T3. axis D = −4a²v⁶Q_ax²/r⁶, so D < 0 where f′ ≠ 0 | confirmed | **proved** by hand, `T3_axis` |
| T4. equator D = v⁶[r²v²a⁴ − (a + rb)²]Q_eq²/(16r⁶) | confirmed | **proved** by hand, `T4_equator` |
| T4. Type I iff v r f′² > \|(r f′)′\| | confirmed, and the proviso Q_eq ≠ 0 can be dropped (gap 2) | **proved** by hand, `T4_equator`, `Qeq0_and_edges` |
| p₈ and p₁₀ as printed in he_type_output.txt | confirmed, term by term | **checked symbolically**, `review_he_compare.py` |
| frame Ricci table as printed | confirmed, ten entries | **checked symbolically**, `review_he_compare.py` |
| (I1), R(e_ϕ, ·) = 0 | confirmed in coordinates and in the frame | **checked symbolically**, `routes_agree` |
| the lemma | confirmed by reading and by planted tensors | `planted_types`, `type_direct` |
| complex pair −R₀₀ ± iR₀₁ on the axis | confirmed | `T3_axis` |
| the equatorial edges are Type II | confirmed, exactly | `Qeq0_and_edges` |
| `he_type.py --no-controls` | reproduces its stored output | rerun, section 4 |

Gaps and an extension.

1. **"p₆v⁶ = −4|j|⁶ < 0"** in the Answer holds where j ≠ 0. It is 0 where J = 0, which is f′ = 0 with s = 0 or f″ = 0,
   and the equatorial circle (r f′)′ = 0. The body already says so for the circle. Suggested wording, "p₆ ≤ 0, and < 0
   wherever j ≠ 0".
2. **The proviso "where Q_eq ≠ 0"** is not needed. On the equator the point is Type I iff v r f′² > |(r f′)′|, Type IV
   iff <, and Type II on equality, at every point with f′ ≠ 0 (section 3.4). Where Q_eq = 0 the discriminant vanishes
   but the point is still Type I, because the e_ρ eigenvalue then meets an eigenvalue of a block that is already
   diagonalizable.
3. **The surfaces where p₁₂ = 0** also include F = 0, that is f = 1. For the tanh profile that is the centre only, which
   he_type.md lists, but the text says "for every profile".
4. **The axis identity** is a polynomial identity in cos θ that he_type.py evaluates at c = 1, where the frame (e_ρ,
   e_ϕ) is undefined. The value is the continuous limit. Here route K evaluates the axis exactly at the regular point
   (x, 0, 0), and route C agrees with it as a polynomial identity.
5. **Shared code.** Route C of he_type.py uses `ricci` and `sph_tetrad_ricci` of `phantom_source.py`, and its two routes
   share `tanh_jets`. Route N is independent but numerical. This review replaces route C by two exact routes, so that
   gap is closed for the four identities.
6. **D = 0 off the equator** was listed as not covered. 200 sampled points are all Type II (section 5), numerical
   evidence, and the Type I subset of D = 0 is not characterized (Recommended tasks).

An identity not in the note falls out of the proof of T2. Q₂ is r²/v⁴ times the discriminant of S, so Q₂ ≥ 0 for every
profile, with equality exactly where S ∝ 1.

## What was read, and when

The task was to take the four targets from he_type.md and the issue and derive them first. Times are UTC, 2026-11-16.

| time | read | before or after the derivation ran |
|---|---|---|
| 08:43–08:49 | AGENTS.md, issue #71, `he_type.md` in full, `review.md` lines 1–149 (§1, the Cartan setup, and the start of §2–§4) | before |
| 08:49–08:57 | nothing from the work under review. `review_he_derive.py` written and run, every check passes and every control fires | |
| 08:57:57 | commit 6a70f61, the derivation and its output | |
| 08:58 | `he_type_output.txt` (outline, and §2 lines 15–45) | after |
| 08:58–09:00 | `he_type.py` (import and function outline, lines 51–135, lines 343–470), `phantom_source.md` lines 1–60 | after |
| 09:00–09:03 | `he_type.py --no-controls` rerun in the background, its output kept in the scratchpad and then copied | after |
| 09:03:28 | commit 387ba04, the comparison, the extension and the rerun output | |

Never read in this review: `phantom_source.py`, `review_cartan.py` and its output, `he_type.py` lines 136–342 and
470–633 (route N, the closed-form class and the Le and Rodal comparisons), and `phantom_source.md` after line 60. After
08:58 `review_he_derive.py` changed only in the labels of three counts, one more control (`import-work`) and a timing
row that was not a check, no identity or engine code (`git diff 6a70f61`). Nothing from the work was imported by any
script, except by that control. `review_he_derive.py` asserts that `he_type` and `phantom_source` are not in
`sys.modules`, and `review_he_compare.py` parses the printed output of `he_type.py` as text.

## 1. The derivation

The spacetime is Alcubierre's, in the comoving coordinates x → x − vt,

    ds² = −dt² + (dx + w dt)² + dy² + dz²,    w = v (1 − f(r)),    r² = x² + y² + z².

The frame is E₀ = n = ∂_t − w∂_x, E₁ = ∂_x, e_ρ and e_ϕ, as in phantom_source.md. With a = f′(r), b = f″(r), F = 1 −
f(r), c = cos θ, s = sin θ and u = 1/r, the second derivatives of w at (r c, 0, r s) are

    w_x = −v a c,   w_ρ = −v a s,   w_ρρ = −v (b s² + a u c²),   w_xx = −v (b c² + a u s²),   w_xρ = −v (b − a u) c s.

D is the discriminant of the cubic of the mixed tensor R^a_b on (E₀, E₁, e_ρ). It is the same in every frame, since it
comes from the characteristic polynomial of the 4×4 R^μ_ν divided by λ. It does not change with the sign convention of
Ricci, since D has degree 6 in the entries, nor with the sign of v, since D is even in v.

Two routes, each written for this review.

    metric ──► Ricci tensor ──► M = ηR on (n, e_x, e_ρ) ──► cubic (c₂, c₁, c₀) ──► D
                 ▲
                 ├── route K   coordinate Ricci of a general metric, own Christoffel code, Cartesian
                 │             (t, x, y, z), at the point (r c, 0, r s). The frame is read off by contraction.
                 └── route C   Cartan's structure equations, own exterior calculus, the coframe
                               (dt, dx + w dt, dρ, ρ dϕ). Solves dθ^a = −ω^a_b ∧ θ^b and checks it,
                               then Ω = dω + ω ∧ ω.

Both are exact. Each was first run on planted truth, the spatially flat FRW metric, where R_tt = −3Ä/A and R_ij = A Ä +
2Ȧ², and Painlevé–Gullstrand Schwarzschild, which is in the family of unit-lapse metrics with flat slices and has R_ab =
0 (`engines`). Both give Alcubierre's Eulerian energy density G(n, n) = −¼v²s²a², that is ρ_E = −(v²/32π)(ρ²/r²)f′²
(`alc_energy`, the formula as quoted in he_type.md, not re-read from the paper *(unverified)*).

The frame Ricci tensor, equal in the two routes as a polynomial identity (`routes_agree`),

    R₀₀ = v² [F (b c² + a s²/r) − a² (c² + ½ s²)]        R₀ₓ = (v/2) [(a/r)(1 + c²) + b s²]
    Rₓₓ = v² [−F (b c² + a s²/r) + a² (c² − ½ s²)]       R₀ρ = −(v/2) (b − a/r) c s
    R_ρρ = ½ v² a² s²                                      Rₓρ = v² c s [a² − ½ F (b − a/r)]

and R(e_ϕ, ·) = 0 (I1). These are the formulas of review.md §1 (R₀₀ = n(w_x) − w_x² − ½|∇⊥w|², and the rest) at the
radial jets, which I checked by hand against the code's output. The routes also agree on the 4×4 coordinate tensor,
whose determinant is 0 and whose cubic is that of the 3-block (`routes_agree`, 10 exact points).

The jets are checked separately against sympy's derivatives of v(1 − f(r)) for f = 1/(1 + r²) and f = r⁻³ + 2r², at five
exact points each, on and off the axis, 140 comparisons (`jets`).

## 2. The cubic and its discriminant

M = ηR on (n, e_x, e_ρ) is

    M = [ −R₀₀   −R₀ₓ   −R₀ρ ]          c₂ = −tr M,   c₁ = the sum of the principal 2×2 minors,   c₀ = −det M,
        [  R₀ₓ    Rₓₓ    Rₓρ ]          D = 18c₂c₁c₀ − 4c₂³c₀ + c₂²c₁² − 4c₁³ − 27c₀².
        [  R₀ρ    Rₓρ    R_ρρ ]

The formula for D agrees with sympy's `discriminant` as a polynomial identity, with the product of squared root
differences for five planted cubics with real and complex roots, and in sign with Sturm's real-root count at 60 random
rational cubics (`disc_tools`). In full, D has 270 terms in (v, F, a, b, u, c, s), and exactly the powers v⁶, v⁸, v¹⁰
and v¹², with 28, 71, 88 and 83 terms.

## 3. The four identities

### 3.1 T1, the powers of v and p₆

The momentum entries R₀ₓ and R₀ρ are exactly linear in v and the rest are exactly quadratic (read off the table, and
`T1_structure` checks the weights of every entry). Conjugating M by diag(−1, 1, 1) reverses the sign of the momentum
entries and keeps the characteristic polynomial, so D is even in the momentum. D is homogeneous of degree 6 in the
entries, so a term with k momentum factors has v-degree k + 2(6 − k) = 12 − k, with k even. That is v⁶, v⁸, v¹⁰ and v¹².

At k = 6 only the momentum remains, M has λ(λ² + |j|²) as its characteristic polynomial, and the discriminant of λ³ + qλ
is −4q³. So p₆v⁶ = −4|j|⁶ with |j|² = R₀ₓ² + R₀ρ², which is the pure-momentum discriminant of the note. From the table,

    4|j|²/v² = [(a/r)(1 + c²) + b s²]² + (b − a/r)² c² s² = s² (b + a/r)² + 4 c² a²/r²,

the coefficients of (a/r)², (a/r)b and b² on the two sides are 1 + 3c², 2s² and s², and the identity is checked
symbolically (`T1_structure`). So |j|² = v² J²/(4r²) with J² = s²(a + rb)² + 4c²a², and p₆ = −J⁶/(16r⁶), as printed.

### 3.2 T2, p₁₂ ≥ 0

At k = 0 the momentum is gone, M = diag(−R₀₀) ⊕ S with S the spatial block of R, and the discriminant of (λ − λ₀)χ_S(λ)
is disc(χ_S) χ_S(λ₀)², with λ₀ = −R₀₀. For a real symmetric 2×2 block disc(χ_S) = (S₁₁ − S_ρρ)² + 4S₁ρ², so

    p₁₂ v¹² = χ_S(λ₀)² [(S₁₁ − S_ρρ)² + 4 S₁ρ²] ≥ 0,

**proved**, with equality exactly when χ_S(λ₀) = 0 or S ∝ 1. The explicit factors are exact identities in both routes,

    χ_S(λ₀) = (1/4) v⁴ s² F Q₁/r²,      (S₁₁ − S_ρρ)² + 4S₁ρ² = v⁴ Q₂/r²,

so p₁₂ = s⁴F²Q₁²Q₂/(16r⁶) with Q₁ and Q₂ as in he_type.md (**checked symbolically**, `T2_p12`). The second identity says
Q₂ ≥ 0.

### 3.3 T3, the axis

On the axis s = 0 and c² = 1, so R_ρρ = Rₓρ = R₀ρ = 0 and e_ρ has eigenvalue 0. The rest is the 2×2 block

    [ A  −B ]      A = −R₀₀ = Rₓₓ = v²(a² − F b),    B = R₀ₓ = v a/r,
    [ B   A ]

with eigenvalues A ± iB, the pair of the note. Then the cubic is λ[(λ − A)² + B²], D = −4B²(A² + B²)², and A² + B² =
v²[a²/r² + v²(F b − a²)²] = v² Q_ax/r². So D = −4a²v⁶Q_ax²/r⁶. Q_ax > 0 wherever a ≠ 0, so D < 0 and the point is Type
IV wherever f′ ≠ 0, **proved**. `T3_axis` checks the identity as a polynomial identity at s = 0 in both routes, at 50
exact axis points per route, and checks that the cubic factors as stated at 20 more on route K.

### 3.4 T4, the equator, and the criterion

On the equator c = 0 and s = 1, and R₀ρ = Rₓρ = 0, so M splits into the (n, e_x) block N₂ and the e_ρ eigenvalue μ =
R_ρρ = ½v²a². The block has tr N₂ = −2v²F a/r, R₀₀ + Rₓₓ = −v²a², R₀ₓ = (v/2)(a + rb)/r, and its discriminant is

    disc(N₂) = (R₀₀ + Rₓₓ)² − 4R₀ₓ² = v² [v²a⁴ − (a + rb)²/r²].

The value of its characteristic polynomial at μ is χ₂(μ) = (v²/4)[4v²a²F(F + ra)/r² + (a + rb)²/r²] = v² Q_eq/(4r²). So

    D = disc(N₂) χ₂(μ)² = v⁶ [r²v²a⁴ − (a + rb)²] Q_eq² / (16 r⁶),

**proved**. The sign of D is the sign of v²a⁴r² − (a + rb)² wherever Q_eq ≠ 0, which is the note's criterion.

The type is settled without D. If disc(N₂) > 0, N₂ has two distinct real eigenvalues, with η-orthogonal, non-null
eigenvectors (the lemma in 1 + 1 dimensions), so M is diagonalizable and has a timelike eigenvector, whatever μ is. That
is Type I, including the points where μ equals one of them, which are the points with Q_eq = 0 and D = 0. If disc(N₂) <
0 the pair is complex, Type IV. If disc(N₂) = 0 and j ≠ 0, N₂ is a non-scalar Jordan block with a null eigenvector, Type
II, or the m = 3, rank 1 case of the lemma if μ meets it. On the equator,

    Type I  ⇔  v r f′² > |(r f′)′|,    Type IV  ⇔  <,    Type II  ⇔  =  (j ≠ 0),

at every point with f′ ≠ 0, **proved**. That is gap 2. `T4_equator` compares D and the criterion with the exact type at
300 random equatorial points on both routes. `Qeq0_and_edges` builds 76 points with Q_eq = 0 per route by solving for f′
and f″ (so F + ra = −Fq² and a + rb = ±2vaFq, with q ≠ 1), and finds each one Type I with the criterion strictly true.
It builds 80 edge points per route with v r a² = |a + rb|, and finds each one Type II. The type comes from the lemma
and, separately, from the definition (a timelike eigenvector, or null eigenvectors with the Jordan sizes of the
nullspaces), and the two agree at every point.

## 4. The comparison with the work under review

`review_he_compare.py` reads `he_type_output.txt` as text, converts its (f, f₁, f₂, r) to (1 − F, a, b, 1/u), and
compares exactly.

| quantity as printed by he_type.py | equals the quantity derived here |
|---|---|
| the ten frame Ricci entries R_(0,0) … R_(ϕ,ϕ) | yes |
| c₂ | yes |
| p₈ (71 terms) and p₁₀ (88 terms) | yes, as polynomials |

The structure check of he_type.py (lines 376–443) is the same argument as 3.1–3.4, done symbolically after a screen at
nine exact points, on the Ricci tensor of `phantom_source.py`. Its identities are those of T1–T4, and all of them hold.
It is the check that `--mutate=alc-drop-w2` breaks. The rerun of `he_type.py --no-controls` takes 1 min 53 s on this
machine at a peak of 1.5 GB ([`review_he_rerun_time.txt`](review_he_rerun_time.txt)), exits 0 with 8 checks passed, and
matches the stored output line for line, except that route N's relative error reads 8.5×10⁻¹⁶ against the stored
8.1×10⁻¹⁶ (floating-point noise) and the Python patch version reads 3.14.10 against 3.14.7. The controls of he_type.py
were not rerun.

## 5. D = 0 off the equator, an extension

he_type.md lists this as not checked. On the equator it is settled (3.4). Off it, `review_he_d0.py` takes random
rational (v, F, a, u, c, s) with c ≠ 0, solves for the real roots b of D at 50 digits, and takes the rank of M − λ at
the double root by singular values. The routine is first run on planted truth, a perfect fluid and dust, with rank 1,
and two Jordan blocks, with rank 2.

| quantity | value |
|---|---|
| sampled points of D = 0, c ≠ 0, j ≠ 0 | 200 |
| rank of M − λ = 2, Type II | 200 |
| rank 1, Type I, or rank 3 | 0 |
| worst \|D\|/scale at the roots | 6.6×10⁻⁴⁸ |
| smallest σ₂/σ₁ of M − λ | 0.0167 |

**Numerical evidence** only, for a sample. A Type I point on D = 0 needs rank(M − λ) = 1, which is a further condition,
and the equator shows such points exist there (Q_eq = 0).

## 6. Checks and controls

Every check has a control. `--mutate=NAME` makes the check fail and exits 1, and the default run confirms that each
control fires.

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `engines` | both engines on FRW and Painlevé–Gullstrand Schwarzschild | `eng-quad` | the ΓΓ terms (K) and the ω∧ω term (C) dropped |
| `jets` | the radial jets against sympy's derivatives | `jets-b` | f″ × 1.001 |
| `alc_energy` | G(n, n) = −¼v²s²a² | `k-drop-w2` | g_tt = −1, the shift's square dropped (route K) |
| `routes_agree` | routes K and C equal, (I1), the 4×4 cubic | `k-drop-w2`, `c-no-nshift` | as above, and E₀ = ∂_t in route C |
| `disc_tools` | the formula for D, three ways | `disc-formula`, `disc-c1` | 17c₂c₁c₀ for 18, and 3c₁³ for 4 |
| `planted_types` | ten tensors of known type, plain and boosted twice and rotated | `lemma-diag` | every D = 0 tensor labelled II |
| `T1_structure` | powers of v, weights, p₆, |j|² | `tgt-p6` | the target p₆ halved |
| `T2_p12` | p₁₂, χ_S(λ₀) and the discriminant of S, p₁₂ ≥ 0 | `tgt-p12` | Q₂ loses its −2ra² term |
| `T3_axis` | the axis form and the complex pair | `tgt-axis` | Q_ax with (Fb + a²) |
| `T4_equator` | the equator form, the split, the criterion | `tgt-eq`, `crit-factor` | Q_eq with 3 for 4, and the criterion with 2\|a + rb\| |
| `type_direct` | the sign of D against the eigenvectors, mpmath | `disc-c1` | as above |
| `Qeq0_and_edges` | Type I at Q_eq = 0, Type II on the edges | `lemma-diag`, `crit-factor` | as above |
| `independence` | the work under review is not imported | `import-work` | `he_type` is imported |

All 13 checks pass and all 13 controls fire. Three things the controls showed.

1. **The dropped shift's square** (`k-drop-w2`) breaks `alc_energy`, `routes_agree`, `T1_structure`, `T2_p12`,
   `T3_axis`, `T4_equator`, `Qeq0_and_edges`, as the issue expected. Route K's entries stop being polynomial in (v, F,
   a, b, u, c, s), since g^{-1} then has 1 + v²F² in its denominator, so under this control the identities are checked
   at exact points only.
2. **`c-no-nshift` does not break `alc_energy`.** The terms n(w_x) cancel in G(n, n) = ½(R₀₀ + Rₓₓ + R_ρρ + R_ϕϕ), so
   Alcubierre's energy density cannot see the shift in n. Route C's other checks catch it.
3. **`disc-formula` does not break `T1_structure` or `T3_axis`.** c₀ = 0 on the axis, and at lowest order in v the term
   c₂c₁c₀ is of higher order. So a wrong coefficient of c₂c₁c₀ cannot show on those two checks, and `disc-c1` was added
   for them.

## What the code proves, and what is assumed

**Checked symbolically**, exact, in the ring ℚ[v, F, a, b, u, c, s] with s² = 1 − c², by both routes. The frame Ricci
tensor, (I1), the cubic and D, the identities T1 to T4, χ_S(λ₀) and the discriminant of S, and the equatorial split. A
polynomial identity in that ring holds on the whole circle c² + s² = 1, the axis and equator included.

**Exact at rational points.** The weights of v, the cubic in v² at 12 points, T2 at 40 points, T3 at 50 and T4 at 300
points on both routes, the types of 156 D = 0 equatorial points per route, the 4×4 coordinate cubic at 10 points.

**Numerical evidence.** `type_direct` (mpmath at 60 digits, 300 points) and section 5.

**Trusted.** sympy 1.14.0 (`diff`, `Matrix.inv`, `cancel`, `sparse polynomial rings`, `discriminant`,
`Poly.count_roots`), mpmath 1.3.0, and python's `Fraction`.

**Assumed.** The metric is Alcubierre's with constant v and the radial profile, and the closed forms are typed from
he_type.md, not derived from it. The lemma is applied as stated, and checked against the definition of the types at the
points above and not proved again here beyond reading it. Section 3 takes the frame Ricci table from the two routes, and
matches it by hand to review.md §1.

## What the checks share with the work under review

Nothing in code. They share the statement of the problem (the metric, the frame, F = 1 − f), the definitions of the
Hawking–Ellis types, and the closed forms being tested. The scripts import sympy, mpmath, the standard library and each
other. Routes K and C share, with each other, the radial jets and the code for the cubic, its discriminant and the ring
arithmetic. Each is checked on its own, the jets against sympy's derivatives (`jets`), the discriminant three ways
(`disc_tools`) and its sign against the eigenvectors (`type_direct`). Both routes' Ricci code is checked on planted
truth (`engines`).

## Recommended tasks

1. **The type on D = 0 off the equator.** What. Characterize the set where rank(M − λ) = 1 on D = 0, which is where the
   point is Type I and not II. Why. Section 5 samples 200 points and finds only Type II, and the equator has an
   exceptional set, Q_eq = 0. Acceptance. An exact description of the Type I subset of D = 0 as polynomial conditions in
   (v, F, a, b, u, c, s), checked at exact rational points of that subset, and a sampled point of D = 0 outside it that
   is Type II, so the description fails if its conditions are loosened. It reduces to the vanishing of the 2×2 minors of
   M − λ at the double root. Depends on #64, #71.
2. **The accelerating bubble.** What. Compute D for b = −v_s(t) f(r_s) with v̇_s ≠ 0, using route C with a
   time-dependent shift, since (I1) holds there (#61) and the same cubic decides. Why. he_type.md lists it as not
   covered. Acceptance. The powers of v_s and v̇_s in D as exact identities, the reduction to the constant-speed forms
   at v̇_s = 0, and a control that drops the shift's square. Depends on #61, #64, #71.
3. **Wording of he_type.md**, gaps 1 to 4, which the coordinator edits. No task.
