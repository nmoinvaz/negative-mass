# Does a sheared warp shift admit the scalar form?

November 2026, issue [#81](https://github.com/nmoinvaz/negative-mass/issues/81), after #70. Script [`sheared.py`](sheared.py),
its stdout [`sheared_output.txt`](sheared_output.txt). It imports `scalar_form`, `ricci_at_point`, `shift_tetrad_ricci` and
`JetD` from [phantom_source.py](phantom_source.py), and `q_rule`, `route_T`, `lemma_residues`, `G`, `ETA` and the profile
jets from [nonaxisymmetric.py](nonaxisymmetric.py). Everything else is in `sheared.py`.

Status words: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails), **numerical
evidence**. A *by hand* step is not independently reviewed.

## Answer

The shift is b = −v f(q) with q = X² + y² + 2z² + Xy/2, X = x − vt. It is the first shift tested in this repository for which the lemma of
[phantom_source.md](phantom_source.md) does not apply.

- **No scalar form at the 29 exact points tested.** At 28 points with R_ab ≠ 0, `scalar_form` returns `algebraic` at ε = −1
  and at ε = +1, and at the one point with R_ab = 0 it returns `trivial` (`points_outcome`, §P of the output). This is
  **numerical evidence**: it is exact rational arithmetic at points, and the tanh jets are rounded. It is not a statement on
  open sets, which is #82.
- **The lemma's premises fail at 21 of the 29 points, and the test fails there anyway.** In the frame adapted to
  w = ∇⊥q, (I2) holds identically and (I1) fails exactly where v f′(q) z ≠ 0 (**checked symbolically**). The measured
  residues are in §P. At the other 8 points the premises hold and `scalar_form` agrees with the lemma (`lemma_agree`).
- **The algebraic stage can be passed at isolated points.** On the z axis, X = y = 0, with tuned jets the sheared Ricci
  tensor is exactly λη − 8π u uᵀ at that point, with ε = −1 (`zaxis_pass`, **numerical evidence**, 60 digits, largest
  2×2 minor 1.2344e−63). The z axis is a line, so by itself this says nothing about open sets. It does show that the
  codimension count of [nonaxisymmetric.md](nonaxisymmetric.md) is no argument for this shift.
- **A reduction by hand** (below, **proved** by hand, not reviewed, offered to #82) says the scalar form can hold at a
  point only where R_ab = 0 or X = y = 0. If it is right, no open set carries a scalar form with R_ab ≠ 0 anywhere on it.
  Code checks its ingredients at points and on planted matrices. It does not prove it.

## Setup

    ds² = −dt² + (dx + b dt)² + dy² + dz²,    b = −v f(q),    q = X² + y² + 2z² + Xy/2,    X = x − vt,

with v ≠ 0 constant and f a C² function on [0, ∞). The form X² + y² + Xy/2 has determinant 1 − 1/16 > 0, so q is positive
definite and f(q) is C² on ℝ³ with ellipsoidal level sets. The shift is not of the form b(t, x, σ(y, z)), because
∇⊥∂_x b = −v(f″ q_X ∇⊥q + f′(½, 0)) has the constant second term (q_X = 2X + y/2). Then

    ∇⊥∂_x b × ∇⊥b = 2 v² f′² z    (`shear_cross`, **checked symbolically**),

so ∇⊥∂_x b is along ∇⊥b exactly where v f′ z = 0.

Frame: E_0 = n = ∂_t − b∂_x, E_1 = ∂_x, E_2 = ∂_y, E_3 = ∂_z. For a transverse vector u = (u_y, u_z), R(u, ·) means
u_y R(E_2, ·) + u_z R(E_3, ·). Write w = ∇⊥q = (q_y, q_z) = (2y + X/2, 4z) and wp = (−q_z, q_y). Then (E_0, E_1, w, wp)
is orthogonal, and e_2 = w/|w|, e_3 = wp/|w| complete the orthonormal frame of the lemma.

## The Ricci tensor

For any b(t, x, y, z) with unit lapse and flat slices, with a, c ∈ {y, z}, ∂_a b the transverse gradient and
D_t b = ∂_t b − b ∂_x b, the 3+1 equations of `route_T` give, by hand,

    R_00 = −|∇⊥b|²/2 − b_x² + b_xt − b b_xx        R_11 = b_x² − |∇⊥b|²/2 − b_xt + b b_xx        R_01 = −Δ⊥b/2
    R_0a = ∂_a b_x / 2                              R_1a = b_x ∂_a b / 2 − ∂_a(D_t b)/2             R_ac = ∂_a b ∂_c b / 2

These match every component of #70's printed Ricci tensor for the unsheared q (compared by hand). For b = −v f(q) with the
sheared q the ten components of route D (`ricci` on the jets f, f′, f″) equal them as polynomial identities in
X, y, z, v, f, f′, f″ (`closed_form`, **checked symbolically**). Route D agrees with route E (Christoffel formula at a
point, sympy's own derivatives of an explicit b) and with route T (the 3+1 equations) at 8 random exact points each
(`xcheck_E`, `xcheck_T`), and with route E in all 16 components at each of the 29 sample points (`points_xcheck`).

For the sheared b: ∂_x b = −v f′ q_X, ∇⊥b = −v f′ w, D_t b = v² f′(1 − f) q_X.

## The lemma's premises

Adapted frame, w = ∇⊥q (`shear_I1`, `shear_I2`, **checked symbolically**, identities in X, y, z, v, f, f′, f″):

    (R(wp, E_0), R(wp, E_1), R(wp, w), R(wp, wp)) = (v f′ z, v² f′ (1 − f) z, 0, 0),        |w|²(R_00 + R_11) + 2R(w, w) = 0.

So (I2) always holds and (I1) fails exactly where v f′ z ≠ 0, in the E_0 component, and in the E_1 component too where f ≠ 1.
(The (I1) components R(wp, w) and R(wp, wp) vanish for every b, since R_ac = ∂_a b ∂_c b/2 and wp ⟂ ∇⊥b.)
Where z = 0 or f′ = 0 the premises hold, so the lemma gives R_ab = 0 if a scalar form holds. The points where w = 0
(z = 0, y = −X/4) have ∇⊥b = 0, and there the frame (E_2, E_3) satisfies the premises as well.

Circular frame, c = (y, z), cp = (−z, y), cpw = cp·w = z(2y − X/2) (`circular_closed`, **checked symbolically**):

    R(cp, E_0) = −(v/2)(f″ q_X cpw − f′ z/2)                  R(cp, cp) = (v² f′²/2) cpw²
    R(cp, E_1) = (v²/2)((f′² − h′) q_X cpw + h z/2)           R(cp, c) = (v² f′²/2) cpw (y q_y + 4z²)
    |c|²(R_00 + R_11) + 2R(c, c) = −v² f′² cpw²,              h = f′(1 − f), h′ = f″(1 − f) − f′²

(I1) and (I2) as stated for (e_ρ, e_ϕ) therefore fail too, and (I2) fails by an amount that vanishes where cpw = 0, at
15 of the 29 points. (I1) fails at the same 21 points as in the adapted frame, which is the set v f′ z ≠ 0 (by hand: where cpw = 0 the
circular frame is the adapted one, and otherwise R(cp, cp) ≠ 0).

How much, from the 29-row table of §P of the output. With ‖R‖_F the Euclidean norm of the 16 frame components, and the
residues taken for unit vectors:

| frame | smallest at (# in §P) | largest at | note |
|---|---|---|---|
| adapted, √(\|I1\|² + \|I2\|²)/‖R‖_F | 0.00206 (24) | 0.0704 (20) | at the 21 points where (I1) fails |
| circular | 0.0012 (18) | 0.171 (17) | larger than the adapted one at 14 of the 21 points, equal at 6, smaller at 1 |

The residue relative to ‖R‖_F is small because R is large at the wall. Its absolute size, the unit-vector (I1) norm, is
between 0.0264 and 2.7 in the adapted frame at those 21 points.

## Exact points

`scalar_form`, the function that passes the Ellis–Bronnikov control in [phantom_source.md](phantom_source.md), is run on the
frame Ricci tensor at 29 exact rational points with ε = −1 and ε = +1 (`points_outcome`). The points are in `POINTS` of
`sheared.py` and in §P of the output.

| profile | points | v | jets |
|---|---|---|---|
| Alcubierre's tanh, σ = 8, R = 1, in r = q^½ | 13, on q ≈ 1: axis (2), z = 0 plane, the critical line w = 0, X = 0, y = 0 planes, general (5), z axis (2) | 1, 1/2, 3, 3/2 | 50-digit values rounded to 45-digit rationals |
| bump (1 − q/4)³ for q < 4, 0 beyond | 6: axis, z = 0 plane, general (3), z axis | 7/4 | exact |
| 1/(1 + q²) | 3 | 2, 3/2 | exact |
| random jets, f′ f″ ≠ 0 | 2 | 8/3, 5/3 | random rationals |
| random jets with f = 1; with f′ = 0 | 1 each | 3/2, 5/2 | exact |
| jets with f″ chosen so that R(w, E_0) = 0, or R(w, E_1) = 0 | 1 each | 2 | exact |
| the bump, beyond it | 1 at q = 13/2 | 7/4 | zero, so R_ab = 0 |

Result: `algebraic` at both signs at the 28 points with R_ab ≠ 0, `trivial` at point 29. The z-axis points 12, 13 and 19, where
R(w, E_0) = R(w, E_1) = 0 automatically, are `algebraic` too, and so are 27 and 28, which satisfy one of those two on purpose.
At the 8 points where the premises hold (1 to 4, 14, 15, 26, 29) the result agrees with the lemma (`lemma_agree`). At all 29
points the adapted-frame residues equal the closed forms above (`points_residues`).

### Planted controls

At the same 29 points a planted R = λη − 8π u uᵀ, with random rational λ and u, goes through the same `scalar_form`
(`planted_pass`, `planted_sign`, `planted_deform`, `planted_premise`).

- With ε = −1 it returns `pass` and recovers λ and u uᵀ exactly, at 29 of 29 points.
- With ε = +1 it returns `sign`, at 29 of 29.
- Plus 1/10 on the E_3E_3 component (rank 2) it returns `algebraic`, at 29 of 29.
- It violates (I1) or (I2) of the adapted frame at 29 of 29, as the lemma requires of a nonzero scalar form.

The planted cases run before the real data in the output (§L before §P).

## The z-axis point

At X = y = 0 the conditions c1 and c2 below hold automatically, and c3 to c5 are three equations on (v, f, f′, f″, z). For
v = 1 and z = 1/2 they are solved by p = v²(1 − f)² = (3√2 − 4)/8 (the positive root of 32p² + 32p − 1 = 0),
f′ = −1/(256 v² (1 − f) z²) and f″ = (−(1 − f)/(64 z²) − 6f′)/(16 z²). Section Z of the script builds those jets at 60 digits
(f = 0.8258446501…, f′ = −0.08971874829…, f″ = 0.1318569451…), takes route D's R there, sets λ = R(E_3, E_3), and finds the
largest 2×2 minor of R − λη to be 1.2344e−63 (the check needs < 1e−40), every diagonal entry ≤ 1e−40 and the trace
negative: R − λη = −8π u uᵀ with ε = −1. With f″ + 1/10 the largest minor is 0.037823 (needs > 1e−6).

Status: **numerical evidence**, since the minors are zero to 60 digits and not rigorously. It shows the algebraic stage of
`scalar_form` is passed at that point of the sheared shift by some jets, so the lemma's non-applicability is real. It
says nothing about closure, potential or any neighborhood: the point lies on a line, which has empty interior. The
formulas come from the reduction below, and the check does not use them to decide, only to choose the point.

By hand, not checked by code: c3 to c5 cannot hold along a segment of the z axis for any profile. With c3, the sum
c3 + c4 = −(r0² + r1²)/λ and (I2) reduce to 1 + p = 1/(32p) at every z, so p = v²(1 − f)² is pinned to its root. Then f is
constant on the segment, f′ = 0, and c3 needs f′ ≠ 0.

## The reduction (by hand, not reviewed, for #82)

At a point with R(e_3, e_3) = R(e_3, e_2) = 0, which holds for every b, let M = R − λη = s u uᵀ in the frame
(E_0, E_1, e_2, e_3), with s ≠ 0. Write r0 = R(e_3, E_0) and r1 = R(e_3, E_1). Then M_33 = −λ and M_a3 = R(e_a, e_3).

- (A) u_3 = 0. Then λ = 0 and r0 = r1 = 0. With (I2), the lemma gives R = 0.
- (B) u_3 ≠ 0. M_23 = 0 gives u_2 = 0, so M_02 = M_12 = M_22 = 0, that is

      c1: R(e_2, E_0) = 0,   c2: R(e_2, E_1) = 0,   λ = R(e_2, e_2),

  and then M = −(1/λ)(r0, r1, 0, −λ)(r0, r1, 0, −λ)ᵀ, which is rank one exactly when

      c3: (R_00 + λ)λ = −r0²,   c4: (R_11 − λ)λ = −r1²,   c5: R_01 λ = −r0 r1.

  λ = R(e_2, e_2) > 0 makes s < 0, so only ε = −1 can pass, and ε = +1 gives `sign`.

For the sheared b, r0 = v f′ z/|w| ≠ 0 where v f′ z ≠ 0, which forces (B). Then c1 and c2, with the identity

    R(w, E_1) − v(1 − f) R(w, E_0) = v² f′² q_X |w|²    (`reduction_identity`, **checked symbolically**),

give v² f′² q_X |w|² = 0, so q_X = 0 (v f′ ≠ 0, w ≠ 0 since q_z = 4z ≠ 0). Since R(w, E_0) = −(v/2)(f″ q_X |w|² + f′ q_y/2),
c1 then gives q_y = 0. The linear map (X, y) ↦ (q_X, q_y) has determinant 4 − 1/4 ≠ 0, so X = y = 0. Where v f′ z = 0 the premises
hold and the lemma applies, with (E_2, E_3) at the points where w = 0.

**Corollary (by hand).** At a point of the sheared shift, any R_ab = λη_ab + s u_a u_b with real λ, s and u forces R_ab = 0 or
X = y = 0. On an open set U where it holds at every point, R_ab = 0 on U minus the z axis, so on U by continuity (f is C²).
Then nothing is left for a scalar, or a perfect fluid, to source. This is the question of #82, and the argument here does not
use the level surfaces of q or an ODE for f. It needs the independent review that #82 can supply.

What the code does with it. `reduction_identity` proves the identity. `reduction_pass` and `reduction_deform` run the
decision procedure of (A), (B) on 4 planted matrices R = r(η − m mᵀ), m = (α, β, 0, 1), α² + β² = 2, rotated to the
coordinate frame, against `scalar_form`: they agree at ε = ±1, and each of c1 to c5 is needed, since adding 1/10 to the
matching entry makes both return `algebraic` (20 of 20), and dropping that condition from the procedure breaks it. At the 29
points the procedure predicts `scalar_form`'s status (`reduction_agree`). None of this proves the corollary.

## Checks and controls

Each check has a control that breaks its input once. Run 2 of the script, the one in `sheared_output.txt`, passes all 20
checks and all 20 controls fire. Run 1 of the same script, before one fix, had 19 of 20 firing.

> Run 1: `--mutate=unsheared` silent on its named target `shear_cross`. Cause: `shear_cross` took its derivatives from the fixed
> sheared q, not the mutable q the README names. Fixed in run 2. Run 2's check lines are identical to run 1's, timings
> ignored, and the other 19 controls' logs are identical. The logs are in `negative-mass-raw/I81-rpi5`.

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `closed_form` | the ten route-D components equal the closed form | `form-drop`, `unsheared`, `lapse` | no −∂_a(D_t b)/2 term; q without Xy/2; lapse 1 + b |
| `shear_cross` | ∇⊥∂_x b × ∇⊥b = 2v² f′² z | `unsheared` | q without Xy/2 on the measured side |
| `shear_I1` | adapted-frame (I1) residues | `unsheared`, `circular`, `lapse` | no shear; w = (y, z); lapse |
| `shear_I2` | adapted-frame (I2) | `circular`, `lapse` | w = (y, z); lapse |
| `circular_closed` | circular-frame closed forms | `circ-swap`, `unsheared` | circle from (z, y); no shear |
| `reduction_identity` | R(w, E_1) − v(1 − f)R(w, E_0) = v² f′² q_X \|w\|² | `red-coef`, `unsheared`, `lapse` | 2v²; no shear; lapse |
| `xcheck_E` | route D = route E at 8 random points | `static` | no t-derivatives in D |
| `xcheck_T` | route D = route T at 8 random points | `static`, `lapse` | the same; lapse 1 + b |
| `points_xcheck` | route D = route E at the 29 points | `jet-g2`, `static` | f″ + 1/7 in E; no t-derivatives |
| `points_residues` | adapted-frame residues equal the closed forms at the 29 points | `circular`, `unsheared`, `planted-in-data` | circular frame; no shear; planted form in the data |
| `points_outcome` | `algebraic` where R ≠ 0, `trivial` where R = 0, both signs | `broken-test`, `planted-in-data` | a test that always passes; planted form in the data |
| `lemma_agree` | scalar_form agrees with the lemma where the premises hold | `broken-test` | the same |
| `reduction_agree` | the procedure of (A), (B) equals scalar_form at the 29 points | `broken-test` | the same |
| `planted_pass` | planted λη − 8π u uᵀ passes, λ and u uᵀ recovered | `planted-sign` | λη + 8π u uᵀ |
| `planted_sign` | the same with ε = +1 gives `sign` | `planted-sign`, `broken-test` | the same; a test that always passes |
| `planted_deform` | the deformed form fails | `planted-null`, `broken-test` | deformation by zero; the same |
| `planted_premise` | a planted form violates the premises | `planted-flat` | the planted form is R = 0 |
| `reduction_pass` | the procedure and scalar_form agree on planted rank-one forms | `red-sign` | the procedure returns `pass` for both signs |
| `reduction_deform` | a 1/10 deformation of c1…c5 returns `algebraic` | `red-no-c1` … `red-no-c5` | the procedure ignores that condition |
| `zaxis_pass` | the sheared R at the z-axis point is a phantom scalar form | `zaxis-off` | f″ + 1/10 on the pass side |

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/sheared.py

It exits 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, `--mutate=NAME` one control
and `--logdir=DIR` also writes each control's stdout to `DIR/NAME.log`. The checks take 40 s and the whole run 247 s on a
Raspberry Pi 5.

## What the code proves, and what is assumed

The code proves, as polynomial identities (`sympy.expand` returning 0) in X, y, z, v and the jets f, f′, f″, the closed form
of the ten components, `shear_cross`, `shear_I1`, `shear_I2`, `circular_closed` and `reduction_identity`. The exact-point tests
are exact rational arithmetic at 29 points. `zaxis_pass` is 60-digit arithmetic at one point. The lemma and the reduction are
by hand. The closed form of the Ricci tensor of a general b(t, x, y, z) is by hand, matched against #70's components and
checked symbolically only for the sheared q and at points for the explicit b of routes E and T.

The checks share machinery as follows. Route D and `closed_form` both take derivatives with `JetD` and `q_rule`, and both
sit on sympy. Routes E and T take sympy's own derivatives of an explicit b, and share with D only the metric's definition
(E) or nothing but sympy (T). The reduction uses `scalar_form` only as the comparison.

Assumed, as in [phantom_source.md](phantom_source.md): one real minimally coupled scalar with
T_ab = ε(∂_aφ ∂_bφ − ½ g_ab (∇φ)²) − g_ab V(φ), G_ab = 8π T_ab, unit lapse and flat slices, v ≠ 0 constant, f a C² function
of q.

Not covered: another sheared shift (this is one fixed q), a lapse other than 1, an accelerating bubble of this shape, two or
more scalars, non-minimal coupling, and matter whose Ricci tensor is not of the form λ g_ab + s u_a u_b. The closure and potential
stages are not reached at any point, since the algebraic stage fails at every sampled point, and they are not tested at the
z-axis point.

The result rests on computer algebra and a by-hand reduction, so it should get an independent review before it is relied on.

## For Lean (#135)

Candidate statements, with their hypotheses. None is written.

1. **Algebra.** Let R be a real symmetric 4×4 matrix in an η-orthonormal frame with R_33 = R_23 = 0. If R = λη + s u uᵀ with
   s ≠ 0, then either u_3 = 0, λ = 0, R_03 = R_13 = 0, or u_3 ≠ 0, u_2 = 0, R_02 = R_12 = 0, λ = R_22 and c3, c4, c5. The
   converse holds with s = −λ/u_3², so sign(s) = −sign(λ). Pure linear algebra, no hypothesis beyond the structure.
2. **Geometry.** The ten closed forms of the frame Ricci tensor of ds² = −dt² + (dx + b dt)² + dy² + dz² for b ∈ C². This is
   the heavy statement. Hypotheses: unit lapse, flat slices, the frame E_0 = ∂_t − b∂_x, and the sign convention of R.
3. **Identity.** For b = −v f(q) with the sheared q, R(w, E_1) − v(1 − f)R(w, E_0) = v² f′² q_X |w|², from item 2.
4. **Linear.** q_X = q_y = 0 implies X = y = 0, since the determinant of (X, y) ↦ (2X + y/2, 2y + X/2) is 15/4.
5. **Corollary.** Items 1 to 4 with the lemma of [phantom_source.md](phantom_source.md) give the corollary above.

## For `nonaxisymmetric.md`

Its section "A shift that is not of this class" can say: *#81 ([sheared.md](sheared.md)) runs b = −v f(q), q = X² + y² + 2z² +
Xy/2, as a result. With w = ∇⊥q, (I2) holds identically and (I1) fails exactly where v f′ z ≠ 0, at 21 of its 29 exact points,
so the lemma does not apply there, and `scalar_form` returns `algebraic` at both signs at all 28 points with R_ab ≠ 0
(**numerical evidence**, points only). The codimension count is not an argument: R_ab has 7 parameters at a point
(v, f, f′, f″, X, y, z) against the 5 of λη + s u uᵀ, and on the z axis X = y = 0 tuned jets give a phantom scalar form
exactly at that point.*
