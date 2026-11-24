# Can three scalars of mixed signs source the Alcubierre wall?

November 2026, issue [#77](https://github.com/nmoinvaz/negative-mass/issues/77). Script
[`three_scalars.py`](three_scalars.py), its stdout [`three_scalars_output.txt`](three_scalars_output.txt). It builds on
[`two_scalars.md`](two_scalars.md) (#65), whose "Not covered" names this case, and on
[`phantom_source.md`](phantom_source.md) (#57, #61), whose frame Ricci tensor and identities (I1), (I2) it imports.

Status words: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails),
**numerical evidence**, **conditional**. #65's lemma has not been independently reviewed (#79). This note's result
uses #65 only through the statement that D = det R₃ vanishes on a nowhere dense part of the wall (see "Rank 3 and the
two patterns"); the three-scalar lemma below is new and is proved here by hand.

## Answer

No at λ = 0, the branch the issue names, for every profile tested and every point tested. At generic wall points
R₃ ⊕ 0 has rank 3 and is indefinite, so the algebraic stage passes with three squares of mixed signs, with the pattern
(+,+,−) where D < 0 and (+,−,−) where D > 0 (**proved**). Both patterns occur at exact wall points. The closure stage
then fails. Closure under O(2,1) or O(1,2) is equivalent, for a rank-3 form with kernel K, to L_K M = 0 and flatness of
the quotient (**proved**). For the wall at λ = 0, K = ∂_ϕ and the quotient is the orbit-space metric h on (t, r, θ).
Its Ricci tensor is nonzero at every one of 25 off-axis exact points (Alcubierre's tanh profile, a bump, random jets)
and at 18 of 18 exact points of an explicit analytic profile (**checked symbolically** at those points). So no (φ₁,
φ₂, φ₃, V) of either sign pattern sources the metric on an open set of the bump profile, since an analytic function is
nonzero on an open dense set once it is nonzero at one point (**checked symbolically**, with the continuation by hand).
*Added by the coordinator:* that sentence is about the λ = 0 branch; the λ ≠ 0 branches rest on 27 roots (**numerical
evidence**, below). The continuation should use the numerator det(h)² Ric(h), which is analytic on all of r < 2 off
the axis. Ric(h) itself is analytic only where det h ≠ 0, a set that may not be connected.
For Alcubierre's tanh profile the jets were rounded to 55 digits, so that is **numerical evidence**, as in #57.

The λ ≠ 0 branches pass the algebraic stage too: 27 admissible rank-3 roots at the 25 off-axis exact points
(exact roots, inertia from 60-digit eigenvalues: **numerical evidence**). They fail closure at all 27, by L_K M ≠ 0, computed as a field at
120 digits (**numerical evidence**). The test that K(λ) ≠ 0 alone, as first coded, failed at 4 of the 27 roots, all on
the equator c = 0 (see "The λ ≠ 0 branches"). The accelerating bubble does not change the answer at the 6 exact
points tested (**checked symbolically** at those points).

What is not decided: **an arbitrary C⁴ profile.** The exact-point results are not an ODE argument. No profile for which
the orbit-space metric is flat on an open set was found, but one is not excluded. The ideal step did not finish
(C5, "Not decided" below). The verdict for an arbitrary profile is therefore **numerical evidence**, not proved.

## Setup

Three minimally coupled scalars with kinetic signs s_I and a potential V(φ¹, φ², φ³) have

    T_ab = Σ_I s_I ∂_aφ^I ∂_bφ^I − g_ab (½ Σ_I s_I (∇φ^I)² + V),

and G_ab = 8π T_ab is equivalent to R_ab = λ g_ab + Σ_I s_I a^I_a a^I_b, with λ = 8πV and a^I = √(8π) dφ^I. The pattern
is (+,+,−) or (+,−,−). A constant target metric of signature (2,1) or (1,2) is the same up to a constant linear map.
The test has four stages, as in two_scalars.md.

1. **Algebraic.** λ is a real root of det(R − λg), and M = R − λg has at most two positive and one negative square
   (pattern A, "(+,+,−)") or at most one positive and two negative squares (pattern B, "(+,−,−)"). The squares are then
   fixed up to Λ(x) ∈ O(2,1) or O(1,2), three functions of position.
2. **Closure.** Some Λ(x) makes all three covectors closed.
3. **Potential.** dλ ∈ span(a^I), so that V is a function of the φ^I.
4. **Residual.** G_ab − 8π T_ab = 0 with the fields substituted back.

**Lemma C (proved).** Let M = Σ_I s_I a^I ⊗ a^I be of rank 3 on an open set of a 4-manifold, with kernel line field K.
Then every a^I(K) = 0, and Λ(x) exists iff L_K M = 0 and the quotient metric q of M by K is flat.

*Proof.* Rank 3 makes the three covectors independent, and M K = 0 gives η_IJ a^J(K) = 0, so a^I(K) = 0. If a^I = dφ^I,
then L_K a^I = d(ι_K a^I) + ι_K da^I = 0, so L_K M = 0, φ is constant along K, and it factors through the local leaf
space Q with M = π*q and q = φ*η flat. Conversely, if L_K M = 0, M descends to a metric q of signature (2,1) or (1,2)
on Q. If q is flat, its Cartesian coordinates φ^I give q = Σ s_I dφ^I dφ^I, and dφ^I = Λ^I_J a^J with Λ(x) ∈ O(p,q),
since two frames of one form with one signature differ by that group. A slice transverse to K maps to Q, and the block
of M on the slice is the pullback of q. In 3D, flat means R_ij = 0. ∎

In the language of #65's ω: da^I = −Ω^I_J ∧ a^J is Cartan's first structure equation, Ω is the Levi-Civita connection
of q, and Λ exists iff its curvature dΩ + Ω ∧ Ω vanishes. For λ ≠ 0 the potential stage is K(λ) = 0 (**proved**: V is
constant along the fibers of φ, which are the integral curves of K).

## Rank 3 and the two patterns

**Lemma A (proved).** Let R_ab be real and symmetric on an orthonormal frame, with R_a3 = 0 and R_00 + R_11 + 2R_22 = 0
(#57's (I1), (I2)). Write R = λη + M, M = N ⊕ (−λ), N = R₃ − λ diag(−1, 1, 1). The weighted trace of N is −2λ.

- Pattern A needs λ ≤ 0. If λ > 0 the e₃ direction is a negative square, so N must be positive semidefinite, which the
  trace −2λ < 0 forbids. For λ < 0, N has at most one positive and one negative square, so det N = 0.
- Pattern B needs λ ≥ 0 by the mirror argument, and det N = 0 for λ > 0.
- At λ = 0, M = R₃ ⊕ 0. Its weighted trace vanishes with positive weights, so R₃ is zero or indefinite. If det R₃ ≠ 0
  the signature is (2,1) or (1,2), and since det is the product of the eigenvalues, (2,1) iff det R₃ < 0. So rank 3 at
  λ = 0 is pattern A where det R₃ < 0 and pattern B where det R₃ > 0. ∎

`lemma3_samples` tests this on 24 random tensors with (I1), (I2), each in a random Lorentz frame: 8 of pattern A, 10 of
pattern B, 6 of rank ≤ 2, and 30 admissible roots λ ≠ 0, each of the sign of its pattern.

On the wall, D = det R₃ has, from the frame Ricci tensor of phantom_source.md, 16 r³ D = (1 − c²) q(c²), with q a
polynomial in c² (`wall_D`, printed in the output; c = cos θ). It vanishes on the axis. By #65 (`det_jets`,
`no_open_set`, not independently reviewed), D = 0 holds on no open subset of the wall for a constant speed, so rank 3
holds on an open dense part. The accelerating bubble was treated there in the same way.

At 32 exact points of the lab-frame Ricci tensor (`wall_points`, through `form3`, the function that runs the planted
forms), λ = 0 is a root at every point. Off the axis the form has rank 3 and the pattern predicted by the sign of D at
all 25: pattern (+,+,−) at 4 (the 3 bump points and 1 random-jet point) and (+,−,−) at 21 (the 18 tanh points at σ = 8,
R_b = 1, speeds 1, ½, 3, and 3 random-jet points). On the axis the rank is ≤ 2 at 7 points, the complex-pair locus of #57.

Rank ≤ 2 remainders with two squares of one sign belong to a pattern (+,+,−) with a³ = 0, that is two normal scalars.
They are not treated here: that is #78.

## Planted controls and deformed ones

Both go through `form3`, `closure3` and the same potential and residual code as the wall.

The planted metrics are warped products ds² = dl² + a(l)² (−dt² + dR² + R² dΘ²) of a line over a flat Minkowski fiber in
polar coordinates, with the three scalars equal to the flat Cartesian coordinates of the fiber, φ^I = κ X^I/√(8π), X =
(t, R cos Θ, R sin Θ). They are written in sheared coordinates (T, R, Th, L), t = T + L/3, Θ = Th + L/5, l = L, so that
K = ∂_L − ⅓∂_T − ⅕∂_Th is not a coordinate vector and Λ(x) is not constant.

| control | warp factor a | λ | signs of (t, X, Y) | kernel |
|---|---|---|---|---|
| P1, a cone | l (κ² = 2) | 0 | (+,−,−) | ∂_l |
| P2 | cosh(l)/√2 (κ = 1) | −3 | (−,+,+) | ∂_l |

For a = αl the Ricci tensor of the fiber part is −2α² η and R_ll = 0; for a = cosh(l)/√2, R_ll = −3 and M = η on the
fiber. Both are exact solutions (`p1_residual`, `p2_residual`: G_ab − 8πT_ab = 0 with the fields substituted back).
Whether they are in print was not searched.

Each passes all four stages. `p1_closed` and `p2_closed` find K, show L_K M = 0 and that the Ricci tensor of the 3 × 3
block on the slice T = const (K^T = 1) vanishes, split M into three squares with Lagrange's reduction, build Λ = J Aᵀ(A Aᵀ)⁻¹
from the gradients J of the known fields and the split A, and check Λ A = J, ΛᵀSΛ = S′ and JᵀSJ = M. Λ is a function
of position (it depends on Θ). `form3` lists every real root of det(R − λg), including roots of lower rank; the test is applied to the root with
rank 3.

The deformed controls pass the algebraic stage and fail closure, each at its own sub-test.

| control | change | passes | fails at |
|---|---|---|---|
| D1 | P1's plane R² dΘ² becomes R² (1 + R²/10) dΘ² | algebraic, L_K M = 0 | the quotient is not flat |
| D2 | P2's warp factor gets a term e^(2l)/5 | algebraic | L_K M ≠ 0, and K(λ) ≠ 0 |

## Closure on the wall at λ = 0

K = ∂_ϕ (the adjugate of M has only its ϕϕ entry), M does not depend on ϕ, and L_K M = 0. The quotient metric is
h = R_ab θ^a θ^b on (t, r, θ), with θ⁰ = dt, θ¹ = w dt + c dr − r sin θ dθ, θ^ρ = sin θ dr + r c dθ, w = v(1 − f).
It does not depend on t. `closure3` computes the Ricci tensor of h at an exact point from the jets of w, in exact
rationals (the same function that finds the planted block flat).

| profile | points | result |
|---|---|---|
| tanh σ = 8, R_b = 1 (jets rounded to 55 digits), v = 1, ½, 3 | 18 off the axis | Ric(h) ≠ 0 at all |
| bump (1 − r²/4)³ for r < 2, v = 7/4 | 3 | Ric(h) ≠ 0 at all |
| random rational jets | 4 | Ric(h) ≠ 0 at all |
| explicit w = v (1 − (1 − r²/4)³), v = 1, ½, 7/4; 3 angles; r = 1, 3/2 | 18 | Ric(h) ≠ 0 at all |

This is `wall_flat_points` (25 of 25) and `wall_flat_explicit` (18 of 18). A rational jet set at one point is the jets of
a polynomial profile, so a nonzero Ricci tensor at one point excludes closure on an open neighborhood for that profile.
For an analytic profile, as the bump inside r < 2, the nonzero set is then open and dense. Whether some non-analytic C⁴
profile makes h flat on an open set was not decided.

The control `wall-flat-sub` replaces h by a flat block; the check then fails.

### Not decided: the arbitrary-profile ideal step (C5)

At a fixed θ₀ the numerators of det(h)² Ric(h) are polynomials in (r, w, …, w⁗), linear in w⁗. If h is flat on a box in
(r, θ), they vanish at every θ₀ and every r of the box. A Gröbner basis of those polynomials, and the ODE step for w(r),
would decide it, as in #65. At θ₀ = π/2 they have 78, 66, 87 and 42 terms; the Gröbner computation over ℚ did not
finish in 5 min, nor did the one modulo a prime with the flow derivatives in 20 min, and both were stopped. The scratch scripts are in
`negative-mass-raw/I77-rpi5/explore/`; they are not claims. **C5 is not decided.** The control `ideal-equator-only` was
not run, since the check it targets does not exist. A smooth profile needs w ∈ C⁴ for the Ricci tensor of h to exist;
for a C² profile the closure condition is not even defined classically.

## The λ ≠ 0 branches

At the 25 off-axis exact points `form3` finds 27 real roots λ ≠ 0 of det(R − λg) with rank 3 and an admissible pattern,
one at each tanh and bump point and one to three at the random-jet points. Lemma A fixes their sign (λ < 0 for pattern
A, λ > 0 for B). At each root the null vector K of N = R₃ − λη₃ and the derivative of λ from the implicit function
theorem give K(λ). Independently, `lie_numeric` takes the Taylor polynomial of w with the point's jets W0..W3, finds λ,
K and M = R − λg at the points r ± h and θ ± h (h = 10⁻³⁵) at 120 digits, and takes central differences, so that L_K M
is a field and not a value.

Result (`wall_lam_potential`): K(λ) ≠ 0 at 23 of 27 roots, K(λ) = 0 at 4, and L_K M ≠ 0 at all 27 (relative size from
1.0 to 35 at the four with K(λ) = 0). The two derivatives of λ agree to 5 × 10⁻⁶⁷, M K = 0 and (L_K M)(K, ·) = 0 hold to
9 × 10⁻⁶⁸. So closure fails at all 27 roots, at the L_K M sub-test.

**Record of C7.** C7 as written in the README was "dλ(K) ≠ 0 or L_K M ≠ 0 at every admissible root". The first script
(b56280c, log `pilot_full_nocontrols.log`) coded only K(λ) ≠ 0. **C7 as coded: FAIL at 4 of 27 roots, all at c = 0**, the
equator, where λ is even in θ − π/2 and K lies in the (E₀, E₁) plane, so K(λ) = 0 by reflection. Restricting the test
to c ≠ 0 would have changed its domain after the data; it was not done. The coordinator had the criterion computed as
written. **C7 as written: PASS**, because L_K M ≠ 0 at those four roots too. The script now tests the OR at every root.

The λ ≠ 0 result is **numerical evidence**: 27 roots at 25 exact points, with roots that are algebraic numbers (their
inertia comes from 60-digit eigenvalues). The equator is not a special locus of the verdict: L_K M ≠ 0 there too.

## The accelerating bubble

For b = −v_s(t) f(r_s) the block of route D of phantom_source.py, at z = 0, y > 0 and unit lapse, gives the same h on
(t, x, y). `acc_const_limit` is the planted truth for this path: at v_s′ = 0, x_s′ = v_s, the lab-route Ricci tensor of
h, moved by the Jacobian of (t, r, θ) → (t, x, y), equals the comoving polar one at 3 exact points. Dropping the
t-derivatives breaks it (`--mutate=acc-static`). `acc_flat_points`: with random jets of v_s (v_s′ ≠ 0), x_s and f at 6
exact points of the half-plane, the Ricci tensor of h is nonzero at all 6.

So the accelerating bubble does not change the answer at those points (**checked symbolically** at exact points). The
algebraic stage for it was not rerun: its D-structure is #65's. The same arbitrary-profile gap applies.

## Checks and controls

Each check has a control that breaks its input once. The script runs every control after the checks and confirms that
its target checks fail, stopping each after the section that holds its targets (`STOP_AT`).

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `planted3_pass` | 12 random λη + Σ s_I a_I a_I, six of each pattern, are recovered | `planted3-sign` | all squares positive |
| `planted3_fail` | 12 forms with three same-sign squares or four squares: the planted λ is not admissible | `planted3-fail-sign` | an admissible pattern |
| `lemma3_samples` | Lemma A on 24 tensors with (I1), (I2) | `lemma3-weight` | R_00 − R_11 + 2R_22 = 0 |
| `p1_algebraic` | the cone: λ = 0 a root, rank 3, inertia (1,2) | `p1-euclid` | a Euclidean fiber |
| `p1_closed` | K, L_K M = 0, flat quotient, Λ(x) in O(p,q) | `p1-fiber` | a curved fiber |
| `p1_potential` | K(λ) = 0 | `p1-lambda` | λ + L/10 |
| `p1_residual` | G_ab − 8πT_ab = 0 | `p1-sign` | one kinetic sign flipped in T |
| `p2_algebraic` | the cosh control: λ = −3, rank 3, inertia (2,1) | `p2-euclid` | a Euclidean fiber |
| `p2_closed` | the same as `p1_closed` | `p2-fiber` | a curved fiber |
| `p2_potential` | K(λ) = 0 | `p2-lambda` | λ + L/10 |
| `p2_residual` | G_ab − 8πT_ab = 0 | `p2-sign` | one kinetic sign flipped in T |
| `d1_algebraic` | D1 passes the algebraic stage at its root | `deform-euclid` | a Euclidean fiber |
| `d1_closure_fails` | D1: L_K M = 0, quotient not flat | `d1-flat` | the curvature set to 0 |
| `d2_algebraic` | D2 passes the algebraic stage at its root | `deform-euclid` | a Euclidean fiber |
| `d2_closure_fails` | D2: L_K M ≠ 0 and K(λ) ≠ 0 | `d2-onshell` | the warp factor on shell |
| `wall_identities` | R(e_ϕ, ·) = 0 and (I2) for the general profile | `alc-drop-w2` | the shift's square dropped from g_tt |
| `wall_D` | 16 r³ det R₃ = (1 − c²) q(c²) | `alc-drop-w2` | the same |
| `wall_points` | 32 exact points: λ = 0 root, rank, pattern predicted by sign D | `alc-drop-w2` | the same |
| `wall_flat_points` | Ric(h) ≠ 0 at 25 off-axis exact points | `wall-flat-sub` | h replaced by diag(−1, 1, 1) |
| `wall_flat_explicit` | Ric(h) ≠ 0 at 18 points of the explicit bump | `wall-flat-sub` | the same |
| `wall_lam_potential` | K(λ) ≠ 0 or L_K M ≠ 0 at 27 roots λ ≠ 0 | `lam-const`, `lie-dk` | dλ = 0 and L_K M = 0; the dK terms dropped |
| `acc_const_limit` | lab-route Ricci of h equals the polar one at constant speed | `acc-static` | no t-derivatives in the jets |
| `acc_flat_points` | Ric(h) ≠ 0 at 6 exact points of the accelerating bubble | `acc-flat-sub` | h replaced by diag(−1, 1, 1) |

All 23 checks pass and all 20 controls fire (`three_scalars_output.txt`). Two things the table does not say. First, the
`wall_*` and `acc_*` checks assert a failure of closure, so a flat block makes them fail, which is what the controls do.
Second, the control `lemma3-weight` is weak: 1 of the 24 samples disagrees with the lemma in that control (`pilot_ctl_lemma3-weight.log`); it fires,
but a lemma that is wrong by a small margin would need more samples.

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/three_scalars.py

It exits 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, `--mutate=NAME` one
control and `--stop=N` the sections up to N. On rpi5 (Raspberry Pi 5, one process, nice 10) the checks took 115 s
and the whole run with its 20 controls 9 min 35 s, at 88 MB peak resident.

## What the code proves, and what is assumed

The code proves, as exact identities or exact ideal-free computations: the planted metrics solve the field equations,
their quotient blocks are flat, Λ(x) ∈ O(p,q) exists, the deformed ones fail at the named sub-test, and at each exact
point of the wall the Ricci tensor of the orbit-space block is nonzero. The λ ≠ 0 results use 60- and 120-digit
numbers.

By hand: Lemma A, Lemma C, K(λ) as the potential stage, and the passage from a point to an open set for analytic
profiles.

Trusted: sympy 1.14.0 and mpmath 1.3.0; `ricci`, `sph_tetrad_ricci`, `lab_tetrad_ricci` and route D of
phantom_source.py (tested there); `form3`, `closure3`, `ricci_pt` and `lie_numeric` are new, tested on the planted
forms and metrics. `lie_numeric` has no planted case with L_K M = 0; it is tested by the identity (L_K M)(K, ·) = 0 and
by the agreement of K(λ) with the implicit-function value (the control `lie-dk` breaks the first).

Assumed: three real, minimally coupled scalars with a constant diagonal kinetic form of signs (+,+,−) or (+,−,−), any
V, G_ab = 8π T_ab, a profile with w ∈ C⁴ where the Ricci tensor of h is used, unit lapse and flat slices for the
shift, and the identities (I1), (I2) of #57.

## Lemmas stated for Lean (#135)

1. **(A)** Lemma A above, for real symmetric R with R_a3 = 0 and R_00 + R_11 + 2R_22 = 0. Status: proved.
2. **(C)** Lemma C above. Status: proved, with the local leaf space of K assumed regular.
3. **(P)** For a rank-3 M = Σ s_I dφ^I dφ^I with constant s_I, K(λ) = 0 where λ is a function of φ. Status: proved.
4. **(D)** At a rank-3 point of the wall, det R₃ < 0 iff signature (2,1). Status: proved.
5. **(E)** (I1), (I2), and the entries of R₃. Status: checked symbolically (#57, #59).

## Not covered

- An arbitrary C⁴ profile: closure is checked at exact points, not by an ODE argument (C5, not decided).
- Two scalars of one sign (#78), reached here as a³ = 0.
- Four or more scalars, and a target metric that is not constant.
- Non-minimal coupling and matter other than scalars.
- Accelerating bubbles: only the closure at 6 exact points; the algebraic stage is not rerun.
- A non-axisymmetric shift.

## What `two_scalars.md`'s "Not covered" should now say

Replace the bullet "Three or more scalars. With λ = 0, R₃ ⊕ 0 has rank 3 at generic wall points and is indefinite, so a
form with three squares of mixed signs passes the algebraic stage there. Closure would decide." by:

> Three scalars of mixed signs: [`three_scalars.md`](three_scalars.md) (#77). The algebraic stage passes on {D ≠ 0}
> with pattern (+,+,−) where D < 0 and (+,−,−) where D > 0; closure under O(2,1) or O(1,2) is flatness of the
> orbit-space metric h. It fails at every exact point tested (λ = 0) and at 27 of 27 roots λ ≠ 0, for the tanh and bump
> profiles and random jets (*checked symbolically* at points, *numerical evidence* for λ ≠ 0), not decided for an
> arbitrary C⁴ profile. Four or more scalars are not covered.

## Prior work

The record of #60 ([prior_work.md](prior_work.md) §1) found no paper that asks whether a scalar field sources the
Alcubierre metric. No further search was made for three scalars. [BD20] fit two- and three-fluid mixtures to the
Alcubierre stress tensor numerically ([prior_work.md](prior_work.md) §2, not re-read here); a minimally coupled scalar
is a perfect fluid only where its gradient is timelike, so their mixtures and these triples overlap without being the
same class.

## References

Entries as in [`../literature.md`](../literature.md), "Warp drives".

- [BD20] Béatrix-Drouhet, arXiv:2012.09941 (2020, preprint).

## Coordinator's review

Raw `I77-mac` (negative-mass-raw). Mac, nice 10, at most 2 processes.

| step | result |
|---|---|
| scope | the three files the issue names, and raw `I77-rpi5` |
| rerun `three_scalars.py` | exit 0 in 2 min 20 s: 23 of 23 checks pass, 20 of 20 controls fire with rpi5's failed-check lists, and the output equals rpi5's apart from the timings |
| re-solve | from the coordinator's own Ricci tensor, Jacobi derivative and bordered system: the 21 roots λ ≠ 0 at the bump and tanh points, their L_K M and K(λ) (B4) |
| independent check `ind77.py` | run 3 exits 0 in 6 s, and all 10 mutations fire on target. Runs 1 and 2 failed in the coordinator's own code (below) |

`ind77.py` imports nothing from the repository. It works in comoving polar coordinates (t, r, θ, ϕ) with a generic w(r)
and its own Ricci routine. It takes h as the (t, r, θ) block of the coordinate Ricci tensor, and computes Ric(h) from
the values and derivatives of h. For λ ≠ 0 it computes dλ by Jacobi's formula and dK from the bordered system
(dN) K + N dK = 0, Kᵀ dK = 0, at 60 digits, where `lie_numeric` uses central differences. The accelerating bubble goes
through a second route: Cartesian coordinates with a polynomial shift, Taylor series truncated at degree 4. It shares
with this note the metric, the frame, the profiles, the point lists, the printed q of `wall_D`, and the claims it tests.

| id | what | run 3 |
|---|---|---|
| B0 | planted: Painlevé–Gullstrand Schwarzschild is Ricci-flat; P1 and P2 have R = λg + κ² Σ s_I dX^I dX^I exactly | pass |
| B1 | planted 3D Ricci (flat polar, S² × time, the P1 and P2 blocks flat, D1's not) and L_K M (P2: 6.6e-59 where 0 is planted; D2′ with a = 1 + L²/2: 0.1) | pass |
| B2 | R_ϕμ = 0, (I2), and 16 r³ det R₃ = (1 − c²) q(c²) with the printed q | pass |
| B3 | 43 points (bump, tanh, the explicit bump, 4 random jet sets of the coordinator's): λ = 0 has kernel ∂_ϕ, the pattern follows the sign of D, Ric(h) ≠ 0 at all 43 | pass |
| B4 | one admissible root λ ≠ 0 at each of the 21 bump and tanh points; pattern A with λ < 0, B with λ > 0; L_K M ≠ 0 at all 21 (relative 0.15 to 0.32); K(λ) = 0 at the 4 roots with c = 0 and ≠ 0 at the other 17 | pass |
| B5 | accelerating bubble: Ric(h) ≠ 0 at 3 exact points; the scalar curvature of h at one point is the same rational number by both routes | pass |

So the patterns, the failure of closure at λ = 0 and the 21 L_K M ≠ 0 verdicts at λ ≠ 0 are confirmed by code that
shares nothing with `form3`, `closure3`, `ricci_pt` or `lie_numeric`. B1 also gives `lie_numeric`'s verdicts the
planted L_K M = 0 case that the section above says it lacks. Lemmas A and C were read by the coordinator, not
reviewed independently; they are row 22 of #135.

**Discrepancies.**
- C7 as coded at b56280c: FAIL at 4 of 27 roots, all at c = 0. It stays in the record; C7 as written passes.
- C5 is not decided, and the control `ideal-equator-only` was not run.
- `deform-euclid` was partly silent on its first run. The planted checks of D1 and D2 were then tightened, which
  makes them stricter. That run's log was overwritten, so raw output was lost.
- The closing comment's proposed text for two_scalars.md said "*checked symbolically* at points" for λ = 0 at every
  profile. For the tanh points it is *numerical evidence*, as this note's "Answer" says. The text applied there says
  so.
- `ind77.py`, the coordinator's check, had three runs. Run 1 crashed in mpmath's `det`. In run 2, B1 failed in two
  planted cases and B3 crashed. All three faults were in the coordinator's code: a coordinate order in the sphere
  test, a relative size that was 0/0 when L_K M = 0, and a radius truncated to 0. They are recorded in `I77-mac`.

**Recommended tasks.**
- A is #157.
- B waits for #157's outcome.
- C is #78.
- D is declined. `ind77.py` checks the closure and L_K M independently, and Lemmas A and C go to Lean (#135, row 22).
- E is declined. At λ = 0 the accelerating bubble already fails closure at exact points (C8, B5). Its λ ≠ 0 branches
  fall under B.
