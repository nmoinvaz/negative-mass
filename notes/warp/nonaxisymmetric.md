# Can a scalar source a non-axisymmetric warp shift?

November 2026, issue [#70](https://github.com/nmoinvaz/negative-mass/issues/70), after #57 and #61. Script
[`nonaxisymmetric.py`](nonaxisymmetric.py), its stdout [`nonaxisymmetric_output.txt`](nonaxisymmetric_output.txt). It uses
the lemma, `scalar_form` and Ricci routines of [phantom_source.md](phantom_source.md).

Status words: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails), **numerical
evidence**.

## Answer

No. For b = −v f(q) with q = (x − v t)² + y² + 2z², no scalar sources the metric at any point where v f′(q) ≠ 0, for
ε = −1, for ε = +1 and for any V (**checked symbolically**). The test fails at its first, algebraic stage, since no λ
makes R_ab − λ g_ab of rank ≤ 1, so the sign, closure and potential conditions are never reached. The wall is the set
where f′(q) ≠ 0, so no scalar sources any part of it.

The exact-point test agrees. At 27 exact points with both signs, 54 tests, the 26 points with R_ab ≠ 0 fail
algebraically and the one with R_ab = 0 is trivial (`nonaxisymmetric_output.txt` §5). Resting on these points alone the
result would be **numerical evidence**, and it is also that for the tanh profile, whose jets are rounded.

The issue expected (I1) and (I2) to fail for this shift, so that the lemma of #57 would not apply. They fail in the
circular frame (e_ρ, e_ϕ), but the lemma holds in any orthonormal frame, and in the frame adapted to the shift they hold
identically in the jets of f. So the lemma applies unchanged. This is a discrepancy with the issue, reported below.

The same holds for every shift b(t, x, σ(y, z)) with unit lapse (**checked symbolically**), which contains the
axisymmetric ones and this one, σ = y² + 2z². It does not cover a shift whose transverse gradient ∇_⊥∂_x b is not along
∇_⊥b. There the premises fail and nothing here decides the question, though the test still fails at all 25 points with
f′ ≠ 0 for one such shift (`--mutate=shear`, below).

## Setup

The metric is that of [phantom_source.md](phantom_source.md) with the shift b a function of q,

    ds² = −dt² + (dx + b dt)² + dy² + dz²,    b = −v f(q),    q = X² + y² + 2z²,    X = x − v t,

with v ≠ 0 constant and f a C² function of q on [0, ∞). Then f(q) is C² on ℝ³. Its level sets are the ellipsoids
q = const, with semi-axes in the ratio 1 : 1 : 1/√2, so f is not a function of r_s. The wall is where f′(q) ≠ 0.

The frame is E_0 = n = ∂_t − b ∂_x, E_1 = ∂_x, E_2 = ∂_y, E_3 = ∂_z. Below, f, f′, f″ are taken at q, m = (y, 2z) with
W = |m|² = y² + 4z², so ∇_⊥q = 2m, and Δ_⊥q = 6.

## The Ricci tensor in the frame

From `nonaxisymmetric_output.txt` §1, with A = f f′ − f′ + 2X²(f f″ + f′² − f″) and C = f f″ + 2f′² − f″,

    R_00 = −2v²(A + W f′²)           R_01 = v (3f′ + 2W f″)
    R_11 =  2v²(A − W f′²)           R_0⊥ = −2v X f″ m
    R_⊥⊥ = 2v² f′² m mᵀ              R_1⊥ = 2v² X C m

where ⊥ is the (E_2, E_3) plane. Every transverse component is along the single direction m, so the vector
e_3 = (−2z, y)/|m| in that plane has R(e_3, ·) = 0. Routes D, E and T agree on these at 8 exact points each
(`xcheck_E`, `xcheck_T`).

## The lemma in the adapted frame

Take e_2 = m/|m| and e_3 = (−2z, y)/|m|, so that (E_0, E_1, e_2, e_3) is orthonormal. Then

- **(I1)** R(e_3, ·) = 0, from R_0⊥, R_1⊥ and R_⊥⊥ all being along m.
- **(I2)** R_00 + R_11 + 2R(e_2, e_2) = −4v² W f′² + 4v² W f′² = 0.

Both are identities in X, y, z, v and the free jets f, f′, f″ (`identity_I1`, `identity_I2`), as polynomial identities
with nothing to reduce modulo. On the axis y = z = 0, where e_2 and e_3 are undefined, R(E_2, ·) = R(E_3, ·) = 0 and
R_00 + R_11 = 0 (`identity_axis`).

The lemma of [phantom_source.md](phantom_source.md) takes an orthonormal frame (e_0, e_1, e_2, e_3) with R_a3 = 0 for
every a and R_00 + R_11 + 2R_22 = 0, and gives R_ab = 0 if R_ab = λη_ab + s u_a u_b. Nothing in its proof fixes the frame
beyond orthonormality. So the scalar form holds at a point only where R_ab = 0.

R_ab ≠ 0 wherever v f′(q) ≠ 0 (`wall_curvature`). Off the axis R(m, m) = 2v² f′² W², and on it R_01 = 3v f′(q). The
sign of s never enters, so ε = −1 and ε = +1 fail alike, and so does any perfect fluid, as in #57.

R_ab = 0 needs more than f′ = 0. At a point off the axis it needs f′ = f″ = 0, and on the axis f′ = 0 and either
f″ = 0, f = 1 or X = 0. These points are not in the wall.

### The discrepancy with the issue

In the circular frame, with e_ϕ = (−z, y)/ρ, the residues printed in §2 of the output are nonzero. They are

    R(e_ϕ, E_0) ρ = −2v X f″ y z,    R(e_ϕ, e_ϕ) ρ² = 2v² f′² y² z²,    R_00 + R_11 + 2R_ρρ = −4v² f′² y² z² / ρ².

So (I1) and (I2), stated for (e_ρ, e_ϕ), fail off the planes y = 0 and z = 0 wherever f′ ≠ 0, as the issue said. They fail
because (e_ρ, e_ϕ) is the wrong frame for this shift, not because R_ab has no such structure. The circular frame is
`--mutate=circular`, which breaks `identity_I1` and `identity_I2` as the issue's premise predicts.

## Any shift b(t, x, σ(y, z))

Section 3 of the script takes the jets B_ijk = ∂_t^i ∂_x^j ∂_σ^k B and σ_y, σ_z, σ_yy, σ_yz, σ_zz as free symbols and the
frame adapted to w = ∇σ. Then (I1), (I2) and R(w, w) = ½(∂_σB)²|w|⁴ hold identically (`sigma_I1`, `sigma_I2`,
`sigma_Rperp`). Route D for this class equals route E at 6 exact points (`sigma_xcheck`).

A general b(t, x, y, z) is not covered. `--mutate=shear` below is a shift outside the class, where (I1) fails.

## Exact points

The same function `scalar_form` that passes the Ellis–Bronnikov control in [phantom_source.md](phantom_source.md) is run
on the frame Ricci tensor at exact rational points, with ε = −1 and ε = +1 (`points_fail`). Section 5 of the output lists
each point with its profile, v, position, q, f′, whether R_ab ≠ 0, whether (I1) and (I2) hold there, and both outcomes.

| profile | points | v | jets |
|---|---|---|---|
| Alcubierre's tanh, σ = 8, R = 1, in r = q^½ | 11 on q = 1, on the axis (2), the y-axis, the planes z = 0, y = 0, X = 0, general (5) | 1, 1/2, 3 | 50-digit values rounded to 45-digit rationals |
| bump (1 − q/4)³ for q < 4, 0 beyond | 6, five on q = 1, one at q = 103/36 | 7/4 | exact |
| 1/(1 + q²) | 4 | 2, 3/2 | exact |
| random jets with f′, f″ ≠ 0 | 4 on q = 1 | random | random rationals |
| random jets with f′ = 0, f″ ≠ 0 | 1 | 5/2 | exact |
| the bump, beyond it | 1 at q = 11/2 | 7/4 | exact, f = 0 |

All 26 points with R_ab ≠ 0 return `algebraic` at both signs, 25 of them with f′ ≠ 0, and the point beyond the bump,
where R_ab = 0, returns `trivial`. At every point (I1) and (I2) hold in the adapted frame (`lemma_premises`), and the
lemma's prediction agrees with `scalar_form` (`lemma_agree`). The two share nothing. The lemma is by hand and
`scalar_form` takes the gcd of the 2×2 minors of R − λη.

For the tanh profile this is **numerical evidence**, since its jets were rounded. The rounded jets are themselves the jets
of some profile, and the general result above covers every profile.

## Control: a planted scalar form

At the same 27 points, R = λη − 8π u uᵀ with random rational λ and a random u with every component nonzero goes through
the same `scalar_form` (`planted_pass`, `planted_sign`, `planted_deform`, `planted_premise`).

- With ε = −1 it returns `pass` and recovers λ exactly and u uᵀ exactly, at 27 of 27 points.
- With ε = +1 it returns `sign`, at 27 of 27 points.
- Deformed by 1/10 on the E_3E_3 component, it returns `algebraic`, at 27 of 27 points.
- It violates (I1) or (I2) of the adapted frame at 27 of 27 points, as the lemma requires of a nonzero scalar form.

## A shift that is not of this class

`--mutate=shear` adds X y/2 to q, so ∇_⊥∂_x b is no longer along ∇_⊥b. The premises then hold at only 11 of the 27
points, and the lemma says nothing where they fail. The same `scalar_form` on the same points still returns `algebraic` at all 25 points with
f′ ≠ 0 (the controls section of the output echoes both counts). That is **numerical evidence** only, from a control, not a
result of this note. The codimension count makes it unsurprising, since a symmetric 4×4 matrix has 10 components and
λη + s u uᵀ has 5 parameters.

## Checks and controls

Each check has a control that breaks its input once. The script runs every control after the checks and confirms that
its target checks fail. All 17 checks pass and all 10 controls fire (`nonaxisymmetric_output.txt`).

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `xcheck_E` | route D equals route E at 8 exact points | `static` | no t-derivatives in D's jets |
| `xcheck_T` | route D equals route T (3+1 equations) at 8 exact points | `static`, `lapse` | the same; lapse 1 + b in D and E |
| `identity_I1` | (I1) in the adapted frame | `circular`, `shear`, `lapse` | (e_ρ, e_ϕ); q + Xy/2; lapse 1 + b |
| `identity_I2` | (I2) in the adapted frame | `circular`, `lapse` | the same |
| `identity_axis` | R(E_2, ·) = R(E_3, ·) = 0, R_00 + R_11 = 0 on the axis | `shear`, `lapse` | the same |
| `wall_curvature` | R(m, m) and R_01 in closed form | `wall-v2`, `lapse` | 2v in the expected values |
| `sigma_xcheck` | route D for b(t, x, σ) equals route E at 6 points | `static` | no t-derivatives |
| `sigma_I1` | (I1) for any b(t, x, σ(y, z)) | `circular`, `lapse` | frame (y, z), not ∇σ; lapse |
| `sigma_I2` | (I2) for any b(t, x, σ(y, z)) | `circular`, `lapse` | the same |
| `sigma_Rperp` | R(w, w) = ½(∂_σB)²\|w\|⁴ | `circular`, `lapse` | the same |
| `planted_pass` | planted λη − 8π u u passes, λ and u u recovered | `planted-sign` | λη + 8π u u |
| `planted_sign` | the same with ε = +1 gives `sign` | `planted-sign`, `broken-test` | the same; a test that always passes |
| `planted_deform` | the deformed form fails | `planted-null`, `broken-test` | deformation by zero; the same |
| `planted_premise` | a planted form violates the premises | `planted-flat` | the planted form is R = 0 |
| `points_fail` | R ≠ 0 gives `algebraic`, R = 0 gives `trivial` | `planted-in-data`, `broken-test` | planted form in the data; the same |
| `lemma_premises` | (I1), (I2) hold at the 27 points | `circular`, `shear`, `lapse`, `planted-in-data` | the same |
| `lemma_agree` | the lemma's prediction equals `scalar_form` | `broken-test` | a test that always passes |

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/nonaxisymmetric.py

It exits 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, and `--mutate=NAME` one
control.

## What the code proves, and what is assumed

The code proves (I1), (I2), the axis identities and the two wall-curvature formulas as exact polynomial identities in X, y,
z, v and the jets f, f′, f″, and the same for b(t, x, σ(y, z)) in its jets. A polynomial is shown to vanish by `sympy.expand`
returning 0. The exact-point test is exact rational arithmetic on finitely many points. The lemma is by hand, and the step
from the identities to the answer is the paragraph after it in [phantom_source.md](phantom_source.md).

Route D's jet rules are tested, not proved, by agreement with two routes that share none of its code. Route E is
`ricci_at_point` of phantom_source.py, which shares only `shift_metric` with D. Route T is the 3+1 equations (momentum
constraint, Hamiltonian constraint, evolution equation for K_ij), written in the script and sharing only the explicit b
with E. D, E and T also share sympy.

Assumed, as in [phantom_source.md](phantom_source.md), one real minimally coupled scalar with
T_ab = ε(∂_aφ ∂_bφ − ½ g_ab (∇φ)²) − g_ab V(φ), G_ab = 8π T_ab, unit lapse and flat slices, v ≠ 0 constant, and f a C²
function of q.

Not covered. A shift whose ∇_⊥∂_x b is not along ∇_⊥b, a lapse other than 1 (with the lapse 1 + b, (I1) and (I2) fail,
`--mutate=lapse`), an accelerating bubble of this shape, two or more scalars, non-minimal coupling, and matter whose
Ricci tensor is not of the form λ g_ab + s u_a u_b. A normal plus a phantom scalar give R_ab − λ g_ab of rank 2 with mixed
signs, which the lemma does not exclude.

The answer rests on computer algebra, so it should get an independent review before it is relied on.
