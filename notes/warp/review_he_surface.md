# Independent review of #72's Type I set on D = 0 off the equator

November 2026, issue [#101](https://github.com/nmoinvaz/negative-mass/issues/101). Reviews [he_surface.md](he_surface.md) (#72).
Scripts [`review_he_surface_lib.py`](review_he_surface_lib.py), [`_ricci.py`](review_he_surface_ricci.py),
[`_typeI.py`](review_he_surface_typeI.py) and [`_points.py`](review_he_surface_points.py), with stdout in
[`_ricci_output.txt`](review_he_surface_ricci_output.txt), [`_typeI_output.txt`](review_he_surface_typeI_output.txt) and
[`_points_output.txt`](review_he_surface_points_output.txt). Nothing here claims novelty.

Status words are those of AGENTS.md. **Checked symbolically** means an exact polynomial identity in sympy, reduced modulo
c² + s² = 1, or exact arithmetic at a point, each with a control that fails. **Certified** means a computer-assisted proof
whose trusted code is named below.

## Answer

T1, T2 and T3 are **confirmed**. For c ≠ 0, s ≠ 0, v ≠ 0, a finite r and (a, b) ≠ (0, 0), a point of D = 0 is Type I
exactly when (r f′)′ = 0, v²(r f′ + F)(r f′ + c²F) + c² = 0 and 2r f′ + (1 + c²)F ≠ 0 (**certified**, section 3). Necessity
is re-derived in the unrotated frame (n, d_x, d_ρ) with an explicit cofactor certificate. Both identities of he_surface.md
hold for symbolic λ. No claim is an **error**. Two are **gaps** in the wording, neither changing T1, T2 or T3:

1. "T2 has real solutions only where −F < r a < −c²F" and "v F s² ≥ 2|c|" assume F > 0. Six of my 24 exact Type I points
   have F < 0, and there both statements fail. The general statements are that r a lies strictly between −F and −c²F in
   either order, and (v F s²)² ≥ 4c².
2. Section 2 uses m_a, m_b, m_c and m_f without defining them. With m_a = AB′ − j_r² and m_b = AC′ − j_θ² the displayed
   identity at a = 0 holds, so the step is correct.

The Gröbner elimination of section 2 was not reproduced. The cofactor certificate below replaces it.

## Independence, and what the code shares

Read: he_surface.md, he_type.md, the issue texts and comments of #72 and #64 as statements. Not opened:
`he_surface.py`, `he_type.py`, `review_he_derive.py`, `phantom_source.py`, their outputs, and raw `I72-*`, `I71-*`, `I64-*`.

One slip. #72's closing comment embeds the coordinator's `rev72.py` in a `<details>` block, and my first print of that
comment showed it before I could strip it. I did not use it. It builds a coordinate Ricci tensor with a fitted polynomial
profile, reads types from sympy's `rank`, and takes exact points from a conic. My route builds the Ricci tensor from the
jet of the shift, takes ranks by Fraction elimination, and places points by my own search, so the methods differ, but I
cannot show that what I read did not shape my choices.

Shared with #72: sympy 1.14.0, mpmath 1.3.0, the metric, the claims under test, and the rule of the lemma in he_type.md
(a repeated root with rank 3 − m is Type I, otherwise Jordan). The lemma is not re-proved here. My own type test
implements it and gets the right type on six planted tensors. The typed table of he_surface.md section 1 is compared with
my Ricci tensor and not used by it.

## Method

Comoving coordinates turn ds² = −dt² + (dx − v f dt)² + dy² + dz² into −dt² + (dx + w dt)² + dy² + dz², w = v(1 − f(r)).
At t = 0 the frame is (n, d_x, d_y, d_z) with n = d_t − w d_x, and at z = 0 these are (n, e_x, e_ρ, e_ϕ). `ricci_jet`
computes the Ricci tensor of this metric for a generic second-order jet (w, w_i, w_ij) with independent symbols, using my
own Christoffel and Ricci code, then sets the jet of w = v(1 − f) at r(c, s, 0). The 3 × 3 block is M = ηR on (n, e_x, e_ρ)
with Σ(λ) = R − λη, so that M − λ = ηΣ and det(x − M) = det Σ(x). The six distinct 2 × 2 minors of Σ are m01,01, m01,02,
m01,12, m02,02, m02,12 and m12,12, where m01,02 has rows (n, e_x) and columns (n, e_ρ).

The rotated frame (n, r̂, θ̂) enters only in `rotated_quantities`, `identity_necessity` and the A, B′ values of
`T3_multiplicity`, which compare my Ricci tensor with he_surface.md. The certificates and D are in the unrotated frame.

## Claims of he_surface.md

| # | claim | verdict | check |
|---|---|---|---|
| 1 | D = 0 off the equator is not Type II throughout | confirmed | `exact_points`: 24 exact points, D = 0, double root, rank 1 |
| 2 | Type I exactly on T1, T2, T3 | confirmed, certified | `cert_T1`, `cert_T2`, `T3_multiplicity`, `exact_points` |
| 3 | proved by hand, identities checked symbolically | confirmed | the rows below |
| 3b | a Gröbner elimination gives the same set by a second route | not reviewed | replaced by the certificate of section 3 |
| 4 | every other point of D = 0 is Type II except a Type III set | confirmed, given the lemma | `T3_multiplicity` (a triple root on T1 has rank ≤ 1), `exact_points`, `type3_exact` |
| 5 | T1, T2, T3 are polynomials; T1 and T2 are codimension 2; Type I is codimension 1 in D = 0 | confirmed | Jacobian of (b + a u, Pu) has rank 2 at all 24 points; D = 0 on T1 ∩ T2 by `D_on_T1` |
| 6 | a random sample of D = 0 never meets them | confirmed (measure zero) | #71's 200 points not rerun |
| 7 | tanh wall: Type I set empty, F + r f′ ≤ −0.054 at the root | confirmed, **numerical evidence**; **proved** in #132 ([he_tanh_bound.md](he_tanh_bound.md)), with the bound −0.0537, not −0.054 | `tanh_scan` |
| 8 | Type III points on the tanh wall at σR_b = 1, 4, 8 | confirmed, **numerical evidence** | `type3_tanh` |
| 9 | a Type III point is exact in Q(√113) | confirmed | `type3_exact` |
| 10 | axis: D = 0 only where v a = 0, M diagonal, Type I | confirmed | `equator_axis`, `axis_points` |
| 11 | the four rows of the table of types | confirmed | `exact_points` (rank 1 double, rank 2 double, rank 1 triple), `type3_exact` (rank 2 triple) |
| 12 | T1 says the momentum j is radial | confirmed | `rotated_quantities`: j_θ ∝ b + a u |
| 13 | T2 real only where −F < r a < −c²F, and then v² = −c²/[(ra + F)(ra + c²F)] | **gap**: v² confirmed, interval needs F > 0 | `exact_points`: 6 points with F < 0 |
| 14 | a triple root with rank 2 cannot lie on T1 | confirmed | `T3_multiplicity`, `cert_T2` |
| 15 | frame Ricci table, R₀₀ + R_xx + 2R_ρρ = 0, R(e_ϕ, ·) = 0 | confirmed | `frame_table`, `ricci_numeric` (4 × 10⁻⁴⁵) |
| 16 | the rotated entries j_r, j_θ, R′_rθ, R′_00, R′_rr, R′_θθ | confirmed, R′_θθ without restricting to T1 | `rotated_quantities` |
| 17 | on T1 the block splits; D = disc(N₂) χ_N₂(μ)² | confirmed | `T3_multiplicity`, `D_on_T1` |
| 18 | D(b = −a u) = v⁶a⁶u⁶(v²F²s⁴ − 4c²)P₂² | confirmed | `D_on_T1` |
| 19 | Type I with a repeated root means rank 1 and a double root | confirmed; rank 0 needs R₀ₓ = R₀ρ = 0, so a = b = 0 (by hand from claim 15) | |
| 20 | T1 is necessary, by j_θ v c a u (1 + v²F²) + m_e + c v F m_d = 0 | confirmed | `identity_necessity`, `cert_T1` |
| 21 | m_d = −j_θ(B′ + c v F j_r), m_e = −j_θ(c v F A + j_r), A + B′ = v²F a u s² | confirmed | `identity_necessity` |
| 22 | the a = 0 step | **gap**: m_a … m_f undefined; correct under my reading | `identity_necessity`, `cert_T1` |
| 23 | A, B′, AB′ − j_r², third root μ = λ + (B′ − A), B′ − A = v²au(2ra + (1 + c²)F) | confirmed | `T3_multiplicity` |
| 24 | disc N₂ = (A − B′)² on T1 ∩ T2, so v F s² ≥ 2\|c\|, equality where T3 fails | confirmed as (v F s²)² ≥ 4c²; **gap** for F < 0 | `D_on_T1` |
| 25 | the six exact-point sets of section 3, with 60, 250 and 3 points | partly reproduced | my own points, below; the counts were not |
| 26 | tanh jets, the scan of section 4, F + r f′ ≈ −σR_b/2 + 1/2 at large σR_b | scan confirmed; the asymptote not reviewed (#132: **numerical evidence** to σR_b = 10⁴) | `tanh_scan` |
| 27 | section 5, the exclusion of #57, section 6, the checks and controls of he_surface.py | not reviewed | the code was not read, by design |

## The certificate

Write Pu = u²P₂ = v²(a + Fu)(a + c²Fu) + c²u², λ* = −½ v² a s² (a + 2Fu), Θ = v² a (2a + (1 + c²)Fu) and
K₁ = ½ v² a u c s (1 + v²F²). Each line is a polynomial identity modulo c² + s² = 1, with λ symbolic (`cert_T1`, `cert_T2`,
`T3_multiplicity`).

    (C1)  −cs m01,01 + (c² − s²) m01,02 + vF c² m01,12 + cs m02,02 + vF cs m02,12 = K₁ (b + a u)

    (C2)  −s m01,12 + c m02,12 = v c a u (λ* − λ) + ¼ F c s² v³ (b + a u)²

    (C3)  m01,01 + m02,02 = −v² a² Pu + (b + a u) A₃ + (λ − λ*) B₃          A₃, B₃ polynomials by exact division

    (C4)  m_k = (b + a u) A_k + (λ − λ*) B_k + Pu C_k,   C = −v²a² (c², cs, 0, s², 0, 0)   for (m01,01, m01,02, m01,12, m02,02, m02,12, m12,12)

    (C5)  χ(x) = det Σ(x) = (x − λ*) [ (x − λ*)(x − λ* − Θ) + v² a² Pu ]  at b = −a u,   Σ(λ*) (0, −s, c)ᵀ = 0 there

    (C6)  at a = 0:  m01,01 = −Σ₀₀² − Σ₀₁²

C1 was found in the unrotated frame by linear algebra over Q at random exact points: 144 unknowns, polynomials in v, F, c, s
of weight 0 under the scaling a → ℓa, b → ℓ²b, u → ℓu, λ → ℓ²λ, and a nullspace of dimension 4 (`--derive`, 43 s, in raw
`derive_T1.log`). C2 to C5 were written from the structure and then verified as identities, not searched. The other three
nullspace vectors give (a²u² − b²) c s² v² (1 + v²F²)/4 ∈ I, which is weaker.

**Necessity.** Let all six minors vanish at a real λ, with v, c, s ≠ 0, u > 0 and (a, b) ≠ (0, 0).
- If a = 0, C6 gives Σ₀₀ = Σ₀₁ = 0, so (v/2) b s² = 0 and b = 0, which is excluded. So a ≠ 0 and K₁ ≠ 0.
- C1 gives b + a u = 0, which is T1.
- C2 at b = −a u gives λ = λ*, since v c a u ≠ 0.
- C3 at b = −a u and λ = λ* gives Pu = 0, since v²a² ≠ 0, which is T2.

**Sufficiency.** On T1, T2 and λ = λ*, C4 puts every minor in zero, so rank Σ ≤ 1, and rank 0 is excluded (claim 19). By C5
the roots are λ* twice and λ* + Θ. The double root is exactly double when Θ ≠ 0, and Θ = v² a u (2 r a + (1 + c²)F) with
v a u ≠ 0, which is T3. Then rank = 1 = 3 − 2 and the point is Type I. If Θ = 0 the root is triple with rank 1: Type II.

**Triple root with rank 2.** On T1, C5 shows a triple root λ forces λ = λ*, χ₂(λ*) = v²a²Pu = 0 and χ₂′(λ*) = −Θ = 0.
Then C4 gives rank Σ ≤ 1. So a triple root of rank 2 has b + a u ≠ 0 (claim 14).

**Discriminant.** At b = −a u, D = disc(χ₂) χ₂(λ*)² with disc(χ₂) = Θ² − 4v²a²Pu = v²a²u²(v²F²s⁴ − 4c²) and
χ₂(λ*) = v²a²Pu, which is claim 18. On T1 ∩ T2 this gives Θ² = v²a²u²(v²F²s⁴ − 4c²).

The certificate is over Q and gives the real statement through the factors v, c, s, u, a, K₁ ≠ 0 and 1 + v²F² > 0.

## Exact points and the tanh wall

All types come from my own `classify` (square-free factors of the characteristic polynomial, rank of M − λ by Fraction
Gaussian elimination, size of the largest Jordan block).

| set | how placed | points | type | rank | root |
|---|---|---|---|---|---|
| T1, T2, T3 | a search over Pythagorean (c, s), rational F, and r a between −F and −c²F with v² a rational square | 24, 6 with F < 0 | I | 1 | double |
| T1, T2, not T3 | v F s² = 2c, r a = −(1 + c²)F/2 | 1 | II | 1 | triple |
| T1, T3, not T2 | on T1 with v F s² = 2c | 1 | II | 2 | double |
| T2, T3, not T1 | b a real root of D(b) = 0 at the first Type I point's other parameters, 60 digits | 1 | II | 2 | double |
| D = 0 off T1 | the real roots of D(b) = 0 at six generic parameter sets, 60 digits | 6 | II | 2 | double |
| triple root off T1 | the recipe of he_surface.md section 3 in Q(√113), my Ricci tensor | 1 | III | 2 | triple |

The Type III point has φ = v F = (−9 + √113)/16, c = 3/5, v = 2 and the null kernel vector e = (1, −7/25, 24/25). I get
F = 0.0509 as he_surface.md does. My a, b and r differ from its a = −0.2394, b = −0.1307 and r = 7.25 by the length scale
ℓ = 5.40 of the symmetry a → ℓa, b → ℓ²b, r → r/ℓ. The invariants agree: a r = −1.735 and b r² = −6.87 in both.

The scan of `tanh_scan` takes 97 log-spaced values of σR_b from 0.01 to 100 with R_b = 1, finds the roots of (r f′)′ = 0 on
[10⁻⁴, 10⁵] by bisection at 60 digits on a 3000-point grid, and finds exactly one at every value. The largest F + r f′ there is −0.05368, at
σR_b = 0.01, r = 100.96, with r f′ = −0.6398 and F = 0.5861, which is he_surface.md's −0.0537. This is **numerical
evidence** at 97 values of one parameter, not a bound. The planted profile f = sech²(r)/10 has a root with F + r f′ = 0.8946
and r f′ < 0, and the same code finds it non-empty. The Newton solves of `type3_tanh` start from the rounded values of
he_surface.md section 4 and converge at 60 digits to v = 1.12823, 2.20885, 2.01943 and cos θ = 0.77676, 0.21157, 0.64368,
with singular values of M − λ of (0.2345, 0.1505, 2 × 10⁻⁶²), (15.12, 3.665, 3 × 10⁻⁶⁰), (40.43, 15.02, 3 × 10⁻⁶¹).

## Checks and controls

Planted truth runs first. Every check was broken once on purpose, and the target lists were written in the raw README
before the first control run and not changed. In the full runs all 25 controls fire, and each exits 1 alone.

| check | what it tests | controls `--mutate=` |
|---|---|---|
| `planted_types` | six planted tensors, 3 × 3 and 4 × 4, boosted and rotated, and deformed copies | `no-rank`, `disc-sign`, `mult-off` |
| `ricci_numeric` | my frame Ricci against 50-digit second differences of the lab metric | `drop-w2`, `no-gamma2` |
| `frame_table` | he_surface.md's typed table minus my Ricci tensor | `table-typo` |
| `equator_axis` | D of my block against he_type.md's equator and axis forms | `eq-typo`, `axis-typo` |
| `rotated_quantities` | section 1's rotated entries | `rot-sign` |
| `identity_necessity` | the necessity identity, A + B′, the a = 0 identity | `tgt-ident` |
| `D_on_T1` | D(b = −a u), disc(χ₂), the split | `tgt-DT1` |
| `cert_T1` | C1, C6 | `cof-perturb`, `cof-drop` |
| `cert_T2` | C2, C3, C4 | `cof2-perturb`, `cof2-drop` |
| `T3_multiplicity` | C5, the eigenvector, Θ, section 2's A and B′ | `tgt-T3` |
| `exact_points` | the Type I points and the points failing one condition | `no-rank`, `loosen-T1`, `loosen-T2`, `loosen-T3`, `locus-shift` |
| `axis_points` | the axis | `axis-typo` (points) |
| `type3_exact` | the Q(√113) point | `nf-modulus` |
| `tanh_scan` | the scan and the planted profile | `tanh-sign` |
| `type3_tanh` | the Newton solves | `type3-off` |

From the repository root, in the environment of [`setup.md`](../../setup.md), with `notes/warp` as the script directory:

    python3 notes/warp/review_he_surface_ricci.py
    python3 notes/warp/review_he_surface_typeI.py
    python3 notes/warp/review_he_surface_points.py

Each exits 0 only if every check passes and every control fires. `--no-controls` runs the checks alone, `--mutate=NAME`
one control (exit 1 when it fires) and, for the second script, `--derive` the search for C1. On rpi5 (one core, nice 10)
the full runs take 5 s, 10 s and 3 min 41 s, with a peak rss of 62 MiB (raw folder `I101-rpi5`, `final_*.time`).

## What the code proves, and what is assumed

**Certified** (exact identities in sympy 1.14.0, modulo c² + s² = 1): C1 to C6, the factorization of D on T1, and the
frame Ricci tensor as a polynomial in (v, F, a, b, u, c, s). The theorem they serve is the Answer above. **Checked
symbolically**: the rotated entries and identities. **Numerical evidence**: the tanh scan, the Newton solves, and the
60-digit points.

Trusted: sympy's `expand`, `rem`, `cancel`, `discriminant`, `factor_list` and `Poly`, mpmath 1.3.0, and my own code
`ricci`, `classify` and `rank_frac`. The Ricci tensor is tested against 50-digit finite differences of the lab-frame metric
at three points (agreement 4 × 10⁻⁴⁵). `classify` is tested on six planted tensors.

Assumed: the lemma of he_type.md, (I1) with e_ϕ = d_z (checked: R(e_ϕ, ·) = 0), constant v, and the hypotheses of the
theorem. Not covered: the Gröbner elimination, the 60, 250 and 3 point counts of section 3, section 4's large-σR_b
asymptote, and the accelerating bubble (#73). The bound for the tanh scan is #132's proof.

## For he_surface.md

The coordinator owns that file. Proposed changes:

1. Answer, paragraph after the table: "T2 has real solutions only where r a lies strictly between −F and −c²F (−F < r a < −c²F
   for F > 0)". Section 2, last paragraph: "(v F s²)² ≥ 4c² there, that is |v F| s² ≥ 2|c|". The table row for T1, T2, not T3:
   "|v F| s² = 2|c|".
2. Section 2: define m_a = AB′ − j_r², m_b = AC′ − j_θ², m_c = B′C′ − τ², m_f = j_r C′ − j_θ τ.
3. Section 2, "A second route": add that the first route has an explicit cofactor certificate, C1 to C6, in this note.

## Recommended tasks

1. **A bound for the tanh scan.** Done in #132 ([he_tanh_bound.md](he_tanh_bound.md)). What. Prove F + r f′ < 0 at the root of (r f′)′ = 0 for every σR_b > 0, or enclose it
   by interval arithmetic over σR_b ∈ [0.01, 100] with an asymptotic argument outside. Why. The empty Type I set of the
   tanh wall rests on a scan at 97 values. Acceptance. A rigorous enclosure, and the planted profile sech²(r)/10 found
   non-empty by the same code. Depends on #72.
2. **The Type III set (#102) by the same method.** Done in #102 ([he_type3.md](he_type3.md)). What. Triple-root conditions c₂² = 3c₁, 27c₀ = c₂³ as polynomials in
   the unrotated frame, with a cofactor certificate that rank 2 forces b + a u ≠ 0, checked at the Q(√113) point.
   Acceptance. An exact description, and a Type II point failing exactly one condition. Depends on #102.
