# The exact Type III set on D = 0 off the equator of the Alcubierre wall

November 2026, issue [#102](https://github.com/nmoinvaz/negative-mass/issues/102). Script [`he_type3.py`](he_type3.py), its stdout
[`he_type3_output.txt`](he_type3_output.txt). It builds on [he_surface.md](he_surface.md) (#72), its independent review
[review_he_surface.md](review_he_surface.md) (#101) and [he_type.md](he_type.md) (#64). Nothing here claims novelty.

Status words are those of AGENTS.md. **Checked symbolically** means an exact polynomial identity in sympy, reduced modulo
c² + s² = 1, or exact arithmetic at a point, with a control that fails. **Certified** means a computer-assisted proof whose
trusted code is named in the last section. **Numerical evidence** means a 60-digit computation with no enclosure.

## Answer

With a = f′, b = f″, F = 1 − f, u = 1/r, c = cos θ, s = sin θ and T = b + a u, take a point off the equator and the axis
(c ≠ 0, s ≠ 0) with v ≠ 0 and (a, b) ≠ (0, 0). Let x³ + c₂x² + c₁x + c₀ be the characteristic polynomial of the block M of
he_type.md, and

    E1 = c₂² − 3c₁,      E2 = 27c₀ − c₂³.

The point is Type III exactly when

    E1 = 0,      E2 = 0,      T = b + a u ≠ 0

(**certified**, section 2). E1 and E2 are polynomials in (v, F, a, b, u, c) with 17 and 42 terms, even in c, and they depend
on v only through φ = vF, A = va and B = vb. The set E1 = E2 = 0 is the Type III set together with the triple-root points of
rank 1 (T = 0), which are Type II. On D = 0 the condition E2 = 0 follows from E1 = 0, so there the description is E1 = 0 and
T ≠ 0 (`cert_disc`). The Type III set has no point on the equator or the axis, where M is block diagonal.

Each condition is needed, on D = 0 and off it. A triple root of rank 1 (E1 = E2 = 0, T = 0) and a double root p with third
root −8p (E2 = 0, E1 ≠ 0) are Type II points of D = 0 that fail exactly one condition, and a Type IV point with E1 = 0 fails
only E2 (section 4).

On the tanh wall the Type III points form a curve in (r, v, cos θ). At Le's parameters (σR_b = 8, v = 1/2) and Rodal's
(σ = 4, R_b = 5, so σR_b = 20, v = 1) the curve meets the wall in three radii each, listed in section 5, each enclosed in a
box with a unique zero (**certified**). That these are all of them is **certified** by an interval cover of every root in
r ∈ (0, ∞) ([he_type3_count.md](he_type3_count.md), #133).

| set, off the equator and axis, on E1 = E2 = 0 | type | rank(M − λ) | root |
|---|---|---|---|
| T ≠ 0 | III | 2 | triple |
| T = 0 (T1, T2, not T3 of he_surface.md) | II | 1 | triple |

## 1. The conditions in closed form

Put φ = vF, X = v²a², γ = v a u, y = v(b − a u) and w = c². The frame entries of he_surface.md section 1, in the frame
(n, e_x, e_ρ), are linear in (X, γ, y),

    R₀₀ = −X(1 + w)/2 + φ(γ + w y)      R₀ₓ = γ + (1 − w) y/2      R₀ρ = −c s y/2
    Rₓₓ = X(3w − 1)/2 − φ(γ + w y)      Rₓρ = c s (X − φ y/2)      R_ρρ = X(1 − w)/2

and M = ηR has

    c₂ = [4φ(γ + w y) − X(3w + 1)] / 2

    E1 = X² − φ[(3w − 1)γ + 2w y] X + (φ² − 3)γ² + (2φ²w + 3w − 3) γ y + ¼(φ²w² + 3φ²w + 3w − 3) y²

and E2 is a cubic form in (X, γ, y) with 10 terms, printed in he_type3_output.txt (`polynomials`). In the original variables
E1 = ¼ v² Ê₁ with Ê₁ of degree 12, and E2 = ⅛ v⁴ Ê₂ with Ê₂ irreducible of degree 17. Under (a, b, u) → (ℓa, ℓ²b, ℓu), E1 and E2
have weights 4 and 6, so only the ratios of (X, γ, y) matter, and E1 = E2 = 0 is a surface in the four effective variables
(φ, w, X : γ : y), which is codimension 2 in (v, F, a, b, u, c). The
Jacobian of (E1, E2) in (v, F, a, b, u, c) has rank 2 at every exact point of section 3.

On a double root p with third root q, E1 = (p − q)² and E2 = (p − q)²(8p + q) (`polynomials`). So E1 = 0 on D = 0 is the triple
root, and a point of D = 0 with E2 = 0 and E1 ≠ 0 has q = −8p.

## 2. The proof

Let Σ(λ) = R − λη on (n, e_x, e_ρ), so that M − λ = ηΣ and rank(M − λ) = rank Σ. Let m_k, k = 1, …, 6, be its six distinct 2×2
minors, in the order m01,01, m01,02, m01,12, m02,02, m02,12, m12,12.

**Rank ≤ 1 forces T = 0** (`cert_rank`, **checked symbolically**; the cofactors are #101's, restated and re-verified here). For every λ,

    K₁ T = −cs m01,01 + (c² − s²) m01,02 + vF c² m01,12 + cs m02,02 + vF cs m02,12,      K₁ = ½ v² a u c s (1 + v²F²).

At a = 0 one has m01,01 = −Σ₀₀² − Σ₀₁² and Σ₀₁ = ½ v b s², so rank ≤ 1 gives b = 0, which is excluded. So a ≠ 0, K₁ ≠ 0, and
a real λ with all six minors zero has T = 0.

**T = 0 with E1 = E2 = 0 forces rank ≤ 1** (`cert_T0`, **checked symbolically**). Write

    λ* = −½ v² a s² (a + 2Fu),   Θ = v² a (2a + (1 + c²)F u),   K = v² a² Pu,   Pu = v²(a + Fu)(a + c²Fu) + c²u².

All of these are polynomial identities modulo c² + s² = 1, with explicit cofactors found by exact division by T:

    χ(x) = det(x − M) = (x − λ*)[(x − λ*)(x − λ* − Θ) + K] + T H(x)
    E1 = Θ² − 3K + T g₁                          (g₁ has 9 terms)
    E2 = Θ³ + 9λ*(Θ² − 3K) + T g₂                (g₂ has 26 terms)
    c₂ + 3λ* + Θ = T h₂
    m_k(λ) = T A_k + (λ − λ*) B_k − K Γ_k,       Γ = (c², cs, 0, s², 0, 0)

At T = 0 and E1 = E2 = 0 the second and third lines give Θ³ = E2 − 9λ*E1 + T(9λ*g₁ − g₂) = 0, so Θ = 0, and 3K = Θ² − E1 + T g₁
= 0. The fourth line gives λ₀ = −c₂/3 = λ*. The fifth gives m_k(λ₀) = 0 for every k, so rank(M − λ₀) ≤ 1. A triple root of
rank 2 therefore has T ≠ 0.

**The description** (`description`). If E1 = E2 = 0, the cubic has the triple root λ₀ = −c₂/3, which is real, so N = M − λ₀ is nilpotent
and rank N ≤ 2. If also T ≠ 0, the first certificate gives rank N ≥ 2. So rank N = 2 and the point is Type III by the lemma of
he_type.md (a 3×3 Jordan block). Conversely a Type III point has E1 = E2 = 0, and T ≠ 0 by the second certificate. ∎

**E1 = 0 and D = 0 give E2 = 0** (`cert_disc`, **checked symbolically**). In (c₂, c₁, c₀),

    27 D + E2² = 27 E1 G,      G = (c₂⁴ + 3c₁c₂² + 36c₁² − 162c₀c₂) / 27,

so on E1 = 0, D = −E2²/27. Hence a point with E1 = 0 and E2 ≠ 0 has D < 0 and is Type IV, and E1 = E2 = 0 gives D = 0.

**No point on the equator or the axis** (`equator_axis`). R₀ρ and Rₓρ carry the factor cs, so at c = 0 and at s = 0 the block M is a
2×2 block plus a 1×1 block. A nilpotent 2×2 block squares to zero (Cayley–Hamilton), so a triple root there has rank ≤ 1.

## 3. The null ray and exact points

A Type III point has a one-dimensional kernel of N spanned by a null vector, and the lemma gives e = (1, e₁, e₂) with e₁² + e₂² = 1.
For fixed (c, e, φ) the equation Σ(λ)e = 0 is three linear equations in (λ, X, γ, y), solved by the cofactor vector n(c, s, e₁, e₂, φ)
(`ray_family`). On that ray, with τ = tr M − 3λ,

    E1 = τ²,      E2 = τ²(9λ + τ),      χ(x) = (x − λ)²(x − λ − τ)    (modulo c² + s² = 1 and e₁² + e₂² = 1),

so a point of the ray is Type III exactly when τ = 0 (and T ≠ 0). τ is a quadratic in φ with leading coefficient c s³/2. At
(c, e) = (3/5, (−7/25, 24/25)) it is (12/625)(8φ² + 9φ − 1), the quadratic of he_surface.md section 3, with root
φ = (−9 + √113)/16. Scaling the ray to X = 1 and choosing v = 2 gives an exact point (`exact_points`, with D = 0, E1 = E2 = 0,
T ≠ 0, rank 2 by exact elimination in Q(√d), N² ≠ 0, N³ = 0, and the Jacobian of rank 2):

| d | c | e₁, e₂ | φ = vF | F | a | b | u = 1/r | T |
|---|---|---|---|---|---|---|---|---|
| 113 | 3/5 | −7/25, 24/25 | 0.10188411 | 0.050942057 | −0.5 | −0.57042403 | 0.28816478 | −0.71450642 |
| 463 | 4/5 | 7/25, 24/25 | 0.61304831 | 0.30652416 | −0.5 | −5.854789 | 0.45536928 | −6.0824736 |
| 949 | 3/5 | 24/25, 7/25 | 0.82524348 | 0.41262174 | −0.5 | 0.56550292 | 0.055923873 | 0.53754099 |
| 1649 | 4/5 | 24/25, −7/25 | 0.066995042 | 0.033497521 | −0.5 | 0.94215683 | 0.20030931 | 0.84200217 |
| 7793 | 7/25 | 3/5, −4/5 | 0.18332578 | 0.091662891 | 0.5 | 0.061455101 | 1.0485803 | 0.58574527 |
| 120367 | 24/25 | 3/5, 4/5 | 0.50814544 | 0.25407272 | 0.5 | 1.9918423 | 0.066114942 | 2.0248997 |

The first row is #72's point: the same c, e, φ and v, with a, b and r on the scale of a ray at X = 1. #72's printed a, b and r
differ by the scale symmetry a → ℓa, b → ℓ²b, r → r/ℓ, and a r = −1.7351, b r² = −6.869 agree with −1.7357 and −6.870 from its
printed values. `he_type_exact` (#64, #103) returns Type III on this point with radicals. The seeds are
outputs of `he_type3.py --search-seeds` (the first is #72's): of 1024 pairs (c, e) with Pythagorean entries from m ≤ 6, the
points with F > 0, u > 0 and T ≠ 0 lie in 350 different fields Q(√d) and none has φ rational (`search_seeds.log` in the raw
folder), so no rational Type III point is known.

## 4. Type II points that fail one condition

All are exact, and all have D = 0 except the last (`type2_failures`, `description`).

| point | type | E1 | E2 | T | rank | root |
|---|---|---|---|---|---|---|
| T1, T2, not T3: \|vF\| s² = 2\|c\|, r a = −(1 + c²)F/2, c = 3/5, 5/13, 8/17 | II | 0 | 0 | 0 | 1 | triple |
| null ray at τ_{−6}(φ) = 0, two points in Q and one in Q(√5658629) | II | ≠ 0 | 0 | ≠ 0 | 2 | p double, q = −8p |
| null ray, φ = 1/3 and 2 | II | ≠ 0 | ≠ 0 | ≠ 0 | 2 | double |
| E1 = 0 on a rational point of the conic at c = 3/5 | IV | 0 | ≠ 0 | ≠ 0 | | D = −E2²/27 < 0 |

Dropping T ≠ 0 from the description admits the first row, dropping E1 = 0 admits the second, and dropping E2 = 0 admits the
last (controls `loosen-T`, `loosen-E1`, `loosen-E2`). The same description rejects every Type I point of #72's family
(T1, T2, T3) and every Type I and IV point of a random sample, and accepts the six Type III points.

## 5. The tanh wall

The profile is f = [tanh σ(r + R_b) − tanh σ(r − R_b)] / (2 tanh σR_b). At fixed r, E1 and E2 are polynomials in v and w = c². With
(φ, A′, B′) = (vF, v r a, v r² b) the resultant in w of E1 and E2 is (27/256) φ³ (A′ − B′)³ P(φ, A′, B′), with P irreducible, of
167 terms, whose degrees in (φ, A′, B′) are 10, 12, …, 18. On the wall P = v¹⁰(Q₀ + Q₁v² + … + Q₄v⁸), so at each r the speeds
of Type III points are among the positive roots of a quartic in v². The common root w is then read from E1 and E2 and the pair
(v, w) is polished by Newton at 60 digits (`curve_at`). Planted truth first: the solver recovers the exact points of section 3
from their jets (`tanh_planted`).

**The curve.** Over r in [0.7, 2] at Le's wall and [4, 7] at Rodal's, the number of Type III points at each r is 1, then 2,
then 1 (`tanh_curve`). They lie on two branches.

- Branch A comes from inside the wall, where v rises steeply as r decreases (v = 9.54 at r = 0.8 for Le, 4.71 at r = 4.5 for
  Rodal), and runs through the wall to the equator. It passes through v = 0 and c = 0 at r_T, the root of (r f′)′ = 0 (where
  T = 0), so that v → 0 and c → 0 on both sides of r_T. It ends on the equator at the radius where r a = −2F, with
  v = |r²b − 2F| / (4F²). There E1 = E2 = 0 with rank 1, a Type II triple root, and he_type.md's two equatorial curves Q_eq = 0
  and v r f′² = |(r f′)′| cross.
- Branch B starts at v = ∞ on the axis side, at the radius where the leading coefficient of the quartic changes sign, with
  c → 1, and then v falls as r grows (v = 0.0946 at r = 3 for Le, 0.0502 at r = 12 for Rodal, with c = 0.9947 and 0.9985).

| wall | σR_b | r_T | A ends on the equator at r (v) | B starts at r | minimum of v on A (grid) |
|---|---|---|---|---|---|
| Le, R_b = 1 | 8 | 1.007762 | 1.130295 (7.35585) | 1.043322 | 0.00727 at r = 1.008 |
| Rodal, R_b = 5 | 20 | 5.006244 | 5.377617 (19.9373) | 5.086643 | 0.0185 at r = 5.005 |

#72's three Newton points (σR_b = 1, 4 and 8 at R_b = 1) lie on the curve: Newton at their v moves r by at most 2 × 10⁻⁵ and cos θ
by at most 5 × 10⁻⁵, which is the rounding of its table (`tanh_curve`). The one at σR_b = 8 is on branch A.

**The points at the papers' speeds.** The sign changes of P(vF, v r a, v r² b) over about 6000 radii, bisected at 60 digits, give
all Type III points at fixed v (`tanh_points`). Each point has two nonzero singular values of M − λ, a third below 10⁻⁵⁹, and
|N³| below 10⁻⁵⁸. Each is also certified by a Krawczyk operator on a box of half-width 10⁻²⁵ in (r, w) at the exact v (mpmath
interval arithmetic), with T ≠ 0 and 0 < w < 1 on the box, and a box shifted by twenty half-widths is rejected
(`tanh_enclosure`). The frame Ricci tensor at each point, from second differences of the lab metric
−dt² + (dx − v f dt)² + dy² + dz² with its own Christoffel code, equals the comoving-jet tensor of the library to 10⁻⁷⁰, and the
lab block alone has a triple root of rank 2 (`lab_ricci`).

| wall, v | r | f | cos θ | θ | T | singular values of M − λ |
|---|---|---|---|---|---|---|
| Le, 1/2 | 0.990568 | 0.5377 | 0.245055 | 75.815° | −8.808 | 4.200, 1.250, 4 × 10⁻⁶⁰ |
| Le, 1/2 | 1.023933 | 0.4054 | 0.202452 | 78.320° | 7.906 | 3.811, 1.021, 5 × 10⁻⁶⁰ |
| Le, 1/2 | 1.094662 | 0.1803 | 0.955441 | 17.169° | 22.03 | 3.551, 1.491, 2 × 10⁻⁶¹ |
| Rodal, 1 | 4.714485 | 0.9076 | 0.978591 | 11.877° | −4.519 | 0.9066, 0.3476, 7 × 10⁻⁶¹ |
| Rodal, 1 | 5.071066 | 0.3616 | 0.533635 | 57.749° | 3.726 | 2.981, 1.092, 2 × 10⁻⁶⁰ |
| Rodal, 1 | 5.103381 | 0.3043 | 0.988448 | 8.717° | 4.972 | 0.9069, 0.4191, 2 × 10⁻⁶⁰ |

All six lie inside Le's wall mask f ∈ [0.1, 0.9] or Rodal's band r ∈ [ρ − 3/σ, ρ + 3/σ] = [4.25, 5.75]. Each radius and angle is one
circle of revolution about the motion axis, and θ → π − θ gives its mirror image. At each paper speed the curve crosses v
three times, once on branch B and twice on branch A, which the quartic route and the sign-change route both count
(`tanh_curve`). The scan found no near-tangent minimum of |P|, so no double crossing is hidden between grid radii, but that is a
grid statement. The count of three is **certified** by #133's interval cover over r ∈ (0, ∞)
([he_type3_count.md](he_type3_count.md)): P has 4 roots at Le and 6 at Rodal, three of them admissible at each wall, and the
spurious factor B′ − A′ has one root per wall, where E1 and E2 have no finite common root.

## 6. What this says for the exclusion of #57

he_type.md leaves D = 0 to the rank test: Type I allows a scalar and a perfect fluid, Type II a scalar and a fluid with a null
flow, and Type III neither. Type III is the only part of D = 0 where the type alone excludes both. The Type III points are
codimension 2 in the six variables of a wall point, so on a wall with a fixed speed they lie on finitely many circles of
revolution. At Le's v = 1/2 and Rodal's v = 1 that is three radii each, a set of zero volume. The volume fractions of he_type.md
(Type IV on 98.8 % of Le's wall and 90.6 % of Rodal's, the thin Type I band) do not change, and none of the points is on the
equator or the axis. So on these two walls the type alone excludes a scalar and a fluid at the Type IV points and, besides, only
on these circles.

## 7. Checks and controls

All 16 checks pass and all 28 controls fire (`he_type3_output.txt`). Planted truth runs first.

| check | what it tests | controls `--mutate=` |
|---|---|---|
| `planted_types` | six planted tensors of known type, boosted and rotated, and deformed copies; E1 and E2 on the triple roots | `no-rank`, `disc-sign`, `E1-typo`, `E2-typo` |
| `polynomials` | the typed table, det(x − M), v-freeness and weights, E1 and E2 on (p, p, q), the form in (X, γ, y) | `table-typo`, `weight-typo`, `E1-typo`, `E2-typo` |
| `cert_rank` | the rank certificate | `cof-rank-perturb`, `cof-rank-drop` |
| `cert_T0` | the chain of section 2 | `cof-g1-perturb`, `cof-g2-drop`, `cof-B-perturb`, `theta-typo`, `E1-typo`, `E2-typo` |
| `cert_disc` | 27 D + E2² = 27 E1 G | `disc27`, `cof-disc-perturb`, `E1-typo`, `E2-typo` |
| `equator_axis` | block diagonal at c = 0 and s = 0 | `axis-typo` |
| `ray_family` | the null ray and τ | `ray-typo`, `E1-typo`, `E2-typo` |
| `exact_points` | six Type III points and #72's, exact rank | `ray-typo`, `nf-radicand`, `phi-off`, `no-rank`, `E1-typo`, `E2-typo` |
| `type2_failures` | Type II points failing one condition, a Type IV point failing E2 | `loosen-T`, `loosen-E1`, `loosen-E2`, `disc-sign` |
| `description` | the description against the type, at exact points of every type | `loosen-T`, `loosen-E1`, `loosen-E2`, `nf-radicand`, `phi-off`, `no-rank`, `disc-sign` |
| `tanh_jets` | own jets against `he_type.tanh_jets` and `mp.diff` | `jets-sign` |
| `tanh_planted` | the solver on planted points | `P-term`, `T-filter-off` |
| `tanh_points` | the six points, 60 digits, rank test | `tanh-off-root` |
| `tanh_curve` | counts, both ends, r_T, crossings, #72's points | `eq-end-typo` |
| `tanh_enclosure` | Krawczyk boxes | `enclosure-shift` |
| `lab_ricci` | the lab-metric frame Ricci tensor at the six points | `lab-sign` |

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/he_type3.py

It exits 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, `--mutate=NAME` one control,
`--checks=a,b` a subset of checks, and `--search-seeds` the search for the seeds of section 3. On rpi5 (one core, nice 10) the checks
take 121 s with a peak rss of 141 MiB and the full run 236 s with 162 MiB. A control that makes a check raise (for example `E1-typo`,
which divides by zero) counts as a failure of that check, and the traceback goes to stderr.

## What the code proves, and what is assumed

**Certified** by sympy 1.14.0 identities (modulo c² + s² = 1 and e₁² + e₂² = 1): the six certificates of section 2, `ray_family`,
and with them the theorem of the Answer, which uses in addition the lemma of he_type.md (a repeated root with rank 3 − m is
Type I, otherwise Jordan), by hand. **Certified** by mpmath 1.3.0 interval arithmetic: the existence and uniqueness of a zero of
(E1, E2) in each of the six Krawczyk boxes, with T ≠ 0 there. **Numerical evidence**: the counts of points on the wall, the
extent of the two branches, and the 60-digit rank tests.

Trusted: sympy 1.14.0 (`div`, `cancel`, `resultant`, `factor_list`, `Poly`), mpmath 1.3.0 (`polyroots`, `findroot`, `iv`), the
frame Ricci tensor of `review_he_surface_lib.py` (#101), checked there against 50-digit second differences and here at the six
points by `lab_ricci`, and `he_type.he_type_exact`. Shared between the proof and the points: the library's Ricci block, which
`polynomials` compares with the typed table of he_surface.md section 1 and `lab_ricci` with the lab metric. The proof uses
only the block, and the exact points take their type from three routes (own, `classify`, `he_type_exact`).

Assumed: constant v and the Alcubierre metric of he_type.md. The theorem is for c ≠ 0, s ≠ 0, v ≠ 0, u finite and (a, b) ≠ (0, 0).

Review: the coordinator's independent check (raw `I102-mac`, `ind102.py`, no code shared with this note's scripts or
with `review_he_surface_lib.py`) derives the frame block from the metric. It reproduces section 1's table and E1, the
certificates of section 2 (the rank certificate, and the chain at T = 0), #72's point in Q(√113) with exact rank, two
tanh points at 50 digits, and the six points by its own resultant scan. That scan found 5 and 7 sign changes of the
resultant against 4 and 6 roots of P. The closing comment on #102 said every extra one is a common root with w outside
(0, 1). That is wrong for one per wall: there the extra sign change is the root of the spurious factor (A′ − B′)³, where E1's
w² coefficient φ²y²/4 and E2's w³ coefficient both vanish, so E1 and E2 share only w = ∞ (#133; raw `I133-mac`, `ind133b.py`
D5 and D6). The w that ind102 printed there, −0.358483 and −6.30241, is E1's finite root, where E2 ≠ 0.

Not covered: the Type III set of the accelerating bubble (#73). A certified count at speeds other than the papers' two. A
rational Type III point.

## Discrepancies

One in the coordinator's closing comment on #102, about the extra sign changes of its scan (Review, above). None with
he_surface.md or he_type.md. Three statements there are now settled: the Type III set is characterized (he_surface.md
"Not covered"), its curve on the tanh wall is mapped (he_surface.md section 4, "The extent of the curve was not mapped"), and
the equator ends of the curve are the crossings of he_type.md's two equatorial curves.
