# The Hawking–Ellis type on D = 0 off the equator of the Alcubierre wall

November 2026, issue [#72](https://github.com/nmoinvaz/negative-mass/issues/72). Script [`he_surface.py`](he_surface.py),
its stdout [`he_surface_output.txt`](he_surface_output.txt). It builds on [he_type.md](he_type.md) (#64) and
[review_he.md](review_he.md) (#71). Nothing here claims novelty.

Status words are those of AGENTS.md. **Checked symbolically** means exact in sympy, or in exact rational and
number-field arithmetic at a point, with a control that fails.

## Answer

D = 0 off the equator is not Type II throughout. With a = f′, b = f″, F = 1 − f, u = 1/r, c = cos θ and s = sin θ, take a
point off the equator and the axis (c ≠ 0, s ≠ 0) with v ≠ 0 and (a, b) ≠ (0, 0). It lies on D = 0 and is Type I exactly
when

    T1    b + a u = 0                                 (r f′)′ = 0
    T2    v² (r a + F)(r a + c² F) + c² = 0
    T3    2 r a + (1 + c²) F ≠ 0

(**proved** by hand, section 2, with every identity **checked symbolically**; a Gröbner elimination of λ gives the same
set by a second route). Every other point of D = 0 there is Type II, except a further set of Type III points. The
conditions are three polynomials in (v, F, a, b, u, c), and T1 and T2 together are codimension 2, so the Type I points
are codimension 1 in the surface D = 0. A random sample of the surface never meets them, which is why #71's 200 points
were all Type II ([review_he.md](review_he.md) §5).

For the tanh wall the Type I set is empty. T2 needs F + r f′ > 0 at the root of (r f′)′ = 0, and a scan of
σR_b from 0.01 to 100 finds F + r f′ ≤ −0.054 there (**numerical evidence**). So on the tanh wall D = 0 off the equator is
Type II, or Type III on a curve of points that Newton at 50 digits finds at σR_b = 1, 4 and 8 (**numerical evidence**).
A Type III point of the general family is exact, in Q(√113) (**checked symbolically**).

On the axis, D = 0 only where v a = 0, and there M is diagonal, Type I.

| set, off the equator and axis, on D = 0 | type | rank(M − λ) | root |
|---|---|---|---|
| T1, T2, T3 | I | 1 | double |
| generic | II | 2 | double |
| T1, T2, not T3 (v F s² = 2\|c\|, r a = −(1 + c²)F/2) | II | 1 | triple |
| triple root with b + a u ≠ 0 | III | 2 | triple |

T1 says the momentum density j is radial. T2 has real solutions only where −F < r a < −c²F, and then it fixes
v² = −c² / [(r a + F)(r a + c² F)]. A triple root with rank 2 cannot lie on T1, because there M is block diagonal with a
2×2 block, so M − λ has rank at most 1 at a triple root.

## 1. The rotated frame

In the frame (n, e_x, e_ρ) of [review_he.md](review_he.md) §1 the Ricci block is, with (R_00 + R_xx + 2R_ρρ = 0),

    R_00 = v² [F (b c² + a u s²) − a² (c² + s²/2)]      R_0x = (v/2) [a u (1 + c²) + b s²]
    R_xx = v² [−F (b c² + a u s²) + a² (c² − s²/2)]     R_0ρ = −(v/2)(b − a u) c s
    R_ρρ = ½ v² a² s²                                      R_xρ = v² c s [a² − ½ F (b − a u)]

`frame_ricci` checks these nine entries equal route C of `review_he_derive.py` as polynomials. The radial direction
r̂ = (c, s) and the tangential θ̂ = (−s, c) in the (x, ρ) plane give the frame (n, r̂, θ̂), where

    j_r = R′_0r = v c a u                          j_θ = R′_0θ = −(v s/2)(b + a u)
    R′_rθ = (v² F c s/2)(b + a u) = −c v F j_θ    R′_00 = v² [F (a u s² + b c²) − a²(1 + c²)/2]
    R′_rr = v² [a²(1 + c²)/2 − F b c²]            R′_θθ = −½ v² a s² (a + 2F u)

(`rotated_frame`, **checked symbolically**). So the tangential momentum j_θ and the coupling R′_rθ both carry the factor
(r f′)′, and both vanish together on T1. On T1 the point splits,

      ┌─────────────────┬──────┐
      │ N₂ on (n, r̂)    │  0   │       D = disc(N₂) · χ_N₂(μ)²,   μ = R′_θθ
      ├─────────────────┼──────┤       the same split as on the equator, with e_ρ replaced by θ̂
      │ 0               │  μ   │
      └─────────────────┴──────┘

and D(b = −a u) = v⁶ a⁶ u⁶ (v²F²s⁴ − 4c²) P₂², with P₂ = v²(r a + F)(r a + c² F) + c² (`D_on_T1`, **checked symbolically**,
D taken from route C). The first factor is disc(N₂) and the second is χ_N₂(μ), the analogue of Q_eq on the equator.

## 2. The proof

Let Σ = R − λη on (n, e_x, e_ρ), so M − λ = ηΣ, and let m_a … m_f be the six distinct 2×2 minors of Σ′ in the rotated
frame, with A = R′_00 + λ, B′ = R′_rr − λ, C′ = R′_θθ − λ, τ = R′_rθ,

    m_d = j_r τ − j_θ B′ ,    m_e = A τ − j_θ j_r .

Type I with a repeated eigenvalue means rank(M − λ) = 1 with a double root. Rank 1 needs every minor to vanish.

**T1 is necessary.** Since τ = −c v F j_θ, the minors factor as m_d = −j_θ (B′ + c v F j_r) and m_e = −j_θ (c v F A + j_r),
and with A + B′ = R′_00 + R′_rr = v² F a u s² one finds, for every λ,

    j_θ · v c a u (1 + v² F²) + m_e + c v F m_d = 0 .

So rank 1 gives j_θ · v c a u (1 + v² F²) = 0. With c, u, v ≠ 0 and 1 + v²F² > 0 this is j_θ = 0 or a = 0. At a = 0,
j_r = 0 and B′ = −A, so m_a = −A², and j_θ⁴ + m_a C′² + 2 A C′ m_b − m_b² = 0 gives j_θ = 0. Hence b + a u = 0
(`necessity_T1`, the identities **checked symbolically**).

**T2 and T3.** On T1, j_θ = τ = 0, so Σ′ is the (n, r̂) block over a single entry C′, and a = 0 would give b = 0, so
a ≠ 0 and j_r ≠ 0. Rank 1 needs C′ = 0, which is λ = R′_θθ, and the block's determinant to vanish. At that λ,

    A = −v² a (a + F u c²),    B′ = v² a (a + F u),    A B′ − j_r² = −v² a² u² P₂ ,

so T2 is the vanishing of the block. The block then has rank 1, with eigenvalues 0 and B′ − A of M₂ − λ, so the third
root of M is μ = λ + (B′ − A) with B′ − A = v² a u (2 r a + (1 + c²) F). The double root λ is simple of multiplicity 2
exactly when T3 holds, and rank 1 = 3 − 2 is Type I by the lemma of he_type.md. If T3 fails the root is triple with rank
1, which the lemma calls Type II (`sufficiency_T1`, **checked symbolically**). ∎

On T1 ∩ T2 the discriminant of N₂ is (A − B′)², so v F s² ≥ 2|c| there, with equality exactly where T3 fails.

**A second route.** In the effective variables φ = vF, X = v²a², γ = v a u, y = v(b − a u), w = c², the six minors are
polynomials in (λ, γ, X, φ, w) after y = 1. A Gröbner basis of their ideal gives, for h = (φ² + 1)(w − 1),

    h · (2γ + 1),   h · (4X² − 2Xφ(1 + w) + w(φ² + 1)),   h · (λ − (w − 1)(X − φ)/2)   all in the ideal,

and the minors lie in the ideal of the three, with y = 0 giving (λ, X, φ² + 1), no real point. Here 2γ + y = v(b + a u),
and on T1 the quadratic is 4 v² a² u² P₂ (`elimination`, **checked symbolically**, sympy 1.14.0's `groebner` trusted).
This route shares only the typed table with the proof above.

## 3. Exact points, and how they were placed on D = 0

Every point below goes through the same path, route C's frame Ricci at the point, then `he_type.he_type_exact`, with the
rank of M − λ taken separately by minors in `Fraction` arithmetic and compared with the lemma. Planted truth runs
first.

**Planted** (`planted`). A Type I tensor λη + σ w wᵀ with w not null, a Type II tensor with a null eigenvector and rank 2,
a Type II tensor with a triple root and rank 1, and a Type III tensor with a 3×3 Jordan block, each boosted, rotated and
boosted again with a decoupled e_ϕ. All four get their type from `he_type_exact` and their rank by minors. Adding
R_00 + 1/7 and subtracting 1/11 from R_22 takes each off D = 0, to Type IV or Type I, so the claim fails on the
deformed copy. Over Q(√5), with a boost of γ = 3/2 and β = √5/3, the number-field routine gives the same four types.

| set | how placed on D = 0 | points | type | rank | root |
|---|---|---|---|---|---|
| T1, T2, T3 | v F s² = c (k + 1/k) and r a = [−(1 + c²)F ± (c/v)(k − 1/k)]/2, c and s a Pythagorean pair | 60 | I | 1 | double |
| generic | a null eigenvector e = (1, m₁, m₂) of R − λη | 250 | II | 2 | double |
| T1, T3, not T2 | on T1 with v F s² = 2c | 3 | II | 2 | double |
| T1, T2, not T3 | the Type I family at k = 1 | 3 | II | 1 | triple |
| T2, T3, not T1 | a real root b of the irreducible quartic factor of D(b), in its quartic field | 1 | II | 2 | double |
| triple root off T1 | the null-eigenvector ray at φ = (−9 + √113)/16, in Q(√113) | 1 | III | 2 | triple |

The Type I points solve T2 by a conic. With x = r a, T2 is x² + (1 + c²)F x + c²F² + c²/v² = 0, whose discriminant
F²s⁴ − 4c²/v² is a square d² when v F s² = c (k + 1/k) and d = (c/v)(k − 1/k). Then 2x + (1 + c²)F = ±d, so T3 holds
exactly when k ≠ 1. Each point is exact, with D = 0 checked by the cubic and discriminant of `he_type.py`.

The generic points use the Jordan lemma. If (R − λη) e = 0 with e null, then e is in the kernel of M − λ and the range lies in e⊥, so either the kernel
has dimension 2, or the range is e⊥, which contains e, and e is in a Jordan block. Either way λ is a repeated root and
D = 0. At fixed (c, e, φ) the equation is linear in
(X, γ, y, λ), so its solutions are a ray. Scaling the ray so X is a square gives rational (v, F, a, b, u). Every one of
the 250 points has D = 0 exactly and rank 2. T1 fails at all 250.

The Type III point is c = 3/5, s = 4/5, v = 2, φ = (−9 + √113)/16 and e = (1, −7/25, 24/25). In the embedding √113 = 10.630
it has F = 0.0509, a = −0.2394, b = −0.1307 and r = 7.25, with b + a u = −0.164 ≠ 0. The quadratic 8φ² + 9φ − 1 is what
tr M − 3λ leaves on the ray.

Dropping T1, T2 or T3 from the description lets one of the Type II points in, which is what the controls `loosen-T1`,
`loosen-T2` and `loosen-T3` show.

## 4. The tanh wall

The profile is f = [tanh σ(r + R_b) − tanh σ(r − R_b)] / (2 tanh σR_b), with the jets of `he_type.tanh_jets`, which agree
with mpmath's derivatives to 9 × 10⁻³¹. T1 is (r f′)′ = 0 and T2 needs (r f′ + F)(r f′ + c²F) < 0, so with F > 0 it needs
−F < r f′ < 0 at the root, that is F + r f′ > 0.

| quantity | value |
|---|---|
| σR_b scanned | 49 log-spaced values, 0.01 to 100 |
| roots of (r f′)′ = 0 on r > 0 | one at every σR_b |
| F + r f′ at the root | negative, at most −0.0537 (σR_b = 0.01, r = 101, r f′ = −0.640, F = 0.586) |
| planted profile f = (1/10) sech² r | root with F + r f′ = 0.895 > 0, found non-empty by the same scan |

So no tanh wall has a point of T1 ∩ T2, and the Type I subset of D = 0 is empty there. The scan is **numerical
evidence** at 49 values of the one parameter σR_b. It is not a bound. At the ends F + r f′ is
close to −σR_b/2 + 1/2 for large σR_b and to −0.054 for small.

Type III on the wall, by Newton at 50 digits on the two equations c₂² = 3c₁ and 27c₀ = c₂³ at fixed r, with the
singular values of M − λ and the powers of N = M − λ (`tanh_type3`, **numerical evidence**):

| σR_b | r (f) | v | cos θ | singular values of M − λ | \|N³\| |
|---|---|---|---|---|---|
| 1 | 2.0656 (0.136) | 1.1282 | 0.7768 | 0.2345, 0.1505, 1 × 10⁻⁵² | 1 × 10⁻⁵³ |
| 4 | 0.898 (0.694) | 2.2089 | 0.2116 | 15.12, 3.665, 6 × 10⁻⁵¹ | 1 × 10⁻⁴⁸ |
| 8 | 0.916 (0.793) | 2.0194 | 0.6437 | 40.43, 15.02, 6 × 10⁻⁵⁰ | 5 × 10⁻⁴⁷ |

All three radii lie inside Le's wall mask f ∈ [0.1, 0.9]. At fixed r these are isolated points, so over r they form a
curve in (v, r, θ). The extent of the curve was not mapped.

## 5. What this says about the exclusion of #57

he_type.md leaves D = 0 to the rank test, with Type I allowing a scalar and a perfect fluid, Type II allowing a scalar and
a fluid with a null flow, and Type III allowing neither. Off the equator on the tanh wall no point of D = 0 is Type I, so
there the points are Type II or Type III. The Type I points of the wall are then the band D > 0 and, on the equator,
the points with Q_eq = 0. For a general profile T1 ∩ T2 ∩ T3 can be non-empty, as the planted profile shows, and at such
points the type alone does not exclude a scalar or a fluid.

## 6. Checks and controls

All 14 checks pass and all 18 controls fire (`he_surface_output.txt`). Planted truth runs first.

| check | what it tests | control `--mutate=` |
|---|---|---|
| `planted` | four planted tensors of known type and their deformed copies | `planted-deform`, `no-rank` |
| `frame_ricci` | the typed table equals route C, R(e_ϕ, ·) = 0 | `table-typo` |
| `rotated_frame` | the entries of section 1 and the identities of section 2 | `tgt-tau`, `table-typo` |
| `necessity_T1` | the certificate of T1 | `tgt-cert` |
| `sufficiency_T1` | the split, A, B′, A B′ − j_r², the third root | `tgt-suff` |
| `D_on_T1` | D(b = −a u) from route C | `tgt-DT1` |
| `axis` | D on the axis, 20 exact points each of a = 0 and a ≠ 0 | `tgt-axis`, `no-rank` |
| `elimination` | the Gröbner elimination, y = 0, and the translation | `tgt-ideal`, `loosen-T1`, `loosen-T2`, `table-typo` |
| `locus_points` | 60 exact Type I points | `locus-T2-shift`, `no-rank` |
| `null_points` | 250 exact points of D = 0 | `null-nonnull` |
| `outside_description` | seven Type II points that fail exactly one condition | `loosen-T1`, `loosen-T2`, `loosen-T3`, `nf-modulus` |
| `number_field_points` | the field arithmetic on planted tensors, the Type III point | `nf-modulus`, `planted-deform` |
| `tanh_scan` | the scan, with a planted profile that must be found non-empty | `tanh-sign` |
| `tanh_type3` | the Type III solves on the tanh wall | `tanh-off-root` |

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/he_surface.py

It exits 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, in 11 s, and
`--mutate=NAME` one control. The full run takes about 3.5 minutes on one core.

## What the code proves, and what is assumed

**Proved** by hand: section 2, from the closed forms of section 1. **Checked symbolically**: the closed forms equal route
C, each identity of sections 1 and 2, and the Gröbner elimination. The exact points are checked in `Fraction` arithmetic
and in Q(√5), Q(√113) and the quartic field, whose elements are polynomials modulo a monic modulus, with exact zero tests.
**Numerical evidence**: the tanh scan and the Type III solves of section 4.

Trusted: sympy 1.14.0 (`groebner`, `Poly`, `factor_list`, `real_roots`, `rem`), mpmath 1.3.0, `he_type.he_type_exact` and
`he_type.cubic`, `he_type.disc` (#64), and route C of `review_he_derive.py` (#71). The closed forms are typed here, not
imported, and the proof uses only them. Shared between the proof and the points: the typed table, which `frame_ricci`
checks against route C. The exact points take their frame Ricci from route C and their types from `he_type_exact`, and the
rank of M − λ comes from minors in `Fraction` arithmetic, which shares nothing with sympy's `rank`.

Assumed: constant v and the Alcubierre metric of he_type.md. The theorem is for c ≠ 0, s ≠ 0, v ≠ 0 and (a, b) ≠ (0, 0).
At (a, b) = (0, 0) or v = 0 the tensor is zero. The equator was settled in #71 and the axis above.

Not covered: the tanh scan is a scan, not a bound. The Type III set was found, not characterized. The accelerating
bubble (#73).

## Discrepancies

1. **`he_type_exact` applied `sp.nsimplify` to an exact rational D** (`he_type.py` line 85, before #103). `sp.nsimplify`
   can return a nearby radical for a rational of large height, and it returns 0 for `Rational(1, 10**60)`. At 1 of the
   356 points that went through `classify` (an a != 0 axis point in the check `axis`), the value differed from the exact
   D, with sign and zero-ness unchanged. A tiny nonzero D would take the rank branch. No result here depended on it, since
   every point has D = 0 or |D| of order 10⁻² or more, and the script compared signs against the exact D.
   **Resolved by #103**: `he_type_exact` decides the sign with `exact_sign` and returns D as `disc` gives it. `classify`
   now compares that D with the exact D by equality at every point and records a failure of the calling check (#131).
2. **#71's 200 points and this result do not disagree.** The Type I set and the Type III set are each codimension 1 in the
   surface D = 0, so a random sample of D = 0 almost surely misses both.

## Recommended tasks

1. **Independent review of #72.** What. Re-derive T1, T2 and T3 from the metric by a route that does not use the rotated
   frame, for example the cofactor certificate of the ideal membership in section 2's second route, and confirm the
   six exact-point sets with a different rank routine. Why. The characterization is a computer-assisted
   claim with no independent review yet, which AGENTS.md asks for before it is relied on. Acceptance. The conditions reproduced, the
   identity j_θ v c a u (1 + v²F²) + m_e + c v F m_d = 0 re-derived without reading `he_surface.py`, and a control that
   makes each fail. Depends on #72.
2. **The Type III set on D = 0 off the equator.** What. Characterize the points with a triple root and rank 2 as
   polynomial conditions in (v, F, a, b, u, c), checked at exact points, and map the curve found on the tanh wall for
   the parameters of Le and Rodal in [he_type.md](he_type.md). Why. Type III excludes both a scalar and a
   fluid, and this is the only part of D = 0 where the type alone does. Acceptance. An exact description, and a Type II
   point of D = 0 that satisfies all but one of its conditions. Depends on #72.
3. **Wording of he_type.md** for the coordinator: "Off the equator, 200 sampled points of D = 0 are all Type II … Which
   points there are Type I is #72" becomes the Answer of this note, and the Type III set is added. No task.
4. **An exact sign in `he_type_exact`.** Done in #103, which added `exact_sign` and planted tensors with D = ±10⁻⁶⁰ to
   `he_type.py`, where `--mutate=nsimplify` makes `planted` fail. The task as written:
   What. Replace the `sp.nsimplify` of D by an exact comparison, and add a planted
   tensor with D = 10⁻⁶⁰. Why. Discrepancy 1. Acceptance. The planted tensor is classified by its sign, and the control
   restoring `nsimplify` makes it fail. Depends on #64.
