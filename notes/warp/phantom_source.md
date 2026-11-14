# Can a phantom scalar source the Alcubierre metric?

November 2026, issues [#57](https://github.com/nmoinvaz/negative-mass/issues/57) and
[#61](https://github.com/nmoinvaz/negative-mass/issues/61) (the accelerating bubble). Script
[`phantom_source.py`](phantom_source.py), its stdout [`phantom_source_output.txt`](phantom_source_output.txt).

Status words: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails), **numerical
evidence**, *(unverified)* for a reference detail not confirmed.

## Answer

No, for every profile and either sign of the kinetic term (**checked symbolically**). At every point where
f′(r_s) ≠ 0, Alcubierre's Ricci tensor is not of the form R_ab = 8π(ε ∂_aφ ∂_bφ + V g_ab), for ε = −1, for ε = +1
and for any value of V. The test fails at its first, algebraic stage, since no λ makes R_ab − λ g_ab of rank ≤ 1. The
sign, closure and potential conditions are never reached.

Every C² profile with f(0) = 1 that decays at infinity has f′ ≠ 0 on an open interval of r_s. So no (φ, V) sources
the metric, and the conditions fail at every point of the bubble wall, the set where f′(r_s) ≠ 0. A C² profile means
that f(r_s) is C² on ℝ³, which forces f′(0) = 0 (**proved**, in "Setup"), so the centre r_s = 0 is not in the wall.

The answer is the same for an accelerating bubble, b = −v_s(t) f(r_s) with r_s² = (x − x_s(t))² + y² + z² and any
v_s(t) and x_s(t) (**checked symbolically**, in "An accelerating bubble"). There the test fails wherever
v_s(t) f′(r_s) ≠ 0. Off the axis, the identities and so the lemma hold for any shift b(t, x, ρ) axisymmetric about the
direction of motion, and the test fails wherever its R_ab ≠ 0.

The obstruction is two identities of the Ricci tensor that hold for every profile, and a five-line lemma. Prior work
that overlaps, the Type IV maps of [Rod26] and [Le26] and the fluid results for u = n, is cited after the lemma.

## Setup

The metric is Alcubierre's Eq. (8) [Alc94]. There the bubble's centre follows x_s(t), with v_s(t) = dx_s/dt and
r_s² = (x − x_s(t))² + y² + z². Up to "An accelerating bubble", v = v_s is constant and x_s = v t,

    ds² = −dt² + (dx − v f(r_s) dt)² + dy² + dz²,    r_s² = (x − v t)² + y² + z².

A profile is C² when f(r_s) is C² on ℝ³. Then f′(0) = 0, meaning f′(r) → 0 as r → 0⁺ (**proved**). On the axis,
g(s) = f(|s|) is C¹, with g′(s) = f′(s) for s > 0 and g′(s) = −f′(−s) for s < 0. Continuity of g′ at s = 0 gives
f′(0) = −f′(0), so f′(0) = 0. Only C¹ is used. So the centre r_s = 0, where θ and the table's s and c below are
undefined, is never a point where f′(r_s) ≠ 0. Alcubierre's tanh profile and the bump of "Explicit profiles" are even
in r and C² on ℝ³.

The scalar has T_ab = ε(∂_aφ ∂_bφ − ½ g_ab (∇φ)²) − g_ab V(φ), and G_ab = 8π T_ab gives
R_ab = 8π(ε ∂_aφ ∂_bφ + V g_ab). A cosmological constant is a constant in V.

The frame is orthonormal. E_0 = ∂_t + v f ∂_x is the Eulerian observer n, E_1 = ∂_x, and e_ρ, e_ϕ are the
cylindrical radial and azimuthal unit vectors about the axis of motion (ϕ is the azimuth, not the scalar φ), so
e_ρ = (y ∂_y + z ∂_z)/ρ with ρ² = y² + z². Let θ be the angle from the direction of motion, so sin θ = ρ/r_s
and cos θ = (x − v t)/r_s.

## The energy density

G_ab n^a n^b / 8π from the sympy Einstein tensor, in the lab coordinates of the issue, equals

    ρ_E = −(v²/32π) ((y² + z²)/r_s²) f′(r_s)²

exactly (check `rho_einstein`). The 3+1 Hamiltonian constraint, 16π ρ_E = K² − K_ij K^ij with
K_ij = ½(∂_iβ_j + ∂_jβ_i) and β_x = −v f, gives the same (`rho_3p1`). Doubling v on the expected side makes both
checks fail (`--mutate=rho-v2`).

This is Alcubierre's Eq. (19), whose arXiv version prints r_c in the denominator where the context requires r_s.

## The Ricci tensor in the frame

With w(r) = v(1 − f(r)), the shift in the comoving coordinate ξ = x − v t, the frame components are, from
`phantom_source_output.txt` §2 (s = sin θ, c = cos θ),

    R_00   = −[2 w w′ s² + 2 r w w″ c² + r w′² (1 + c²)] / (2r)
    R_01   = −[w′ (1 + c²) + r w″ s²] / (2r)
    R_0ρ   = (r w″ − w′) s c / (2r)
    R_11   =  [2 w w′ s² + 2 r w w″ c² + r w′² (3c² − 1)] / (2r)
    R_1ρ   = [w (r w″ − w′) + 2 r w′²] s c / (2r)
    R_ρρ   = w′² s² / 2
    R_aϕ   = 0   for every a.

The lab-coordinate computation and the comoving spherical one agree at 24 exact rational points, with random v and
jets (`ricci_xcheck`).

Two identities follow, for every C² profile and every v. "An accelerating bubble" extends both to a time-dependent
v_s(t) and x_s(t).

- **(I1)** R(e_ϕ, ·) = 0.
- **(I2)** R_00 + R_11 + 2 R_ρρ = 0.

They are checked in comoving spherical coordinates for a general w(r) (`identity_sph`), in lab Cartesian coordinates
for a general f with the jets f, f′, f″ as free symbols (`identity_lab`), and on the axis y = z = 0, where they read
R(E_2, ·) = R(E_3, ·) = 0 and R_00 + R_11 = 0 (`identity_axis`).

Since G_11 − 3G_00 = −(R_00 + R_11 + 2R_ρρ + 2R_ϕϕ) and R_ϕϕ = 0 by (I1), (I2) is equivalent to G_11 = 3G_00
([prior_work.md](prior_work.md), "The identities"). That is the identity G_zz = 3G_nn stated for Alcubierre in [SSV22]
Sec. 7.1. It is implicit in [LV04a] Eq. 16, where the average of the null energy along ± the direction of motion is
4ρ, and it is behind p = 3μ in [SAR21a] (Eq. 3.4) and p_x = 3ρ in [ABV23] (Eqs. 30–31).

The curvature does not vanish in the wall. R_ρρ = v² f′² (y² + z²)/(2 r_s²) off the axis, and R_01 = v f′(r_s)/r_s on
it (`wall_curvature`). So R_ab ≠ 0 at every point where f′(r_s) ≠ 0.

## The lemma

**Lemma (proved).** Let R_ab be real and symmetric on an orthonormal frame (e_0, e_1, e_2, e_3), with R_a3 = 0 for every
a and R_00 + R_11 + 2 R_22 = 0. If R_ab = λ η_ab + s u_a u_b with λ and s real and u a real covector, then R_ab = 0.

*Proof.* If s = 0 or u = 0 then R = λη, and R_33 = λ = 0. Otherwise R_a3 = s u_a u_3 = 0 for a ≠ 3, and
R_33 = λ + s u_3² = 0.

- If u_3 ≠ 0, then u_0 = u_1 = u_2 = 0 and λ = −s u_3² ≠ 0. Then R_00 + R_11 + 2R_22 = −λ + λ + 2λ = 2λ ≠ 0, a
  contradiction.
- If u_3 = 0, then λ = 0 and R_00 + R_11 + 2R_22 = s(u_0² + u_1² + 2u_2²) = 0. The form is positive definite, so
  u = 0 and R = 0. ∎

Apply it with (e_2, e_3) = (e_ρ, e_ϕ) off the axis and (E_2, E_3) on it, s = 8πε and λ = 8πV. The scalar form then
holds at a point only where R_ab = 0, and R_ab ≠ 0 wherever f′(r_s) ≠ 0. The sign of s never enters, so the phantom
sign does not help and the normal sign fails the same way.

The lemma excludes any R_ab = λ g_ab + s u_a u_b. A perfect fluid of any equation of state, with a timelike or null
flow, has a Ricci tensor of that form, so the same identities exclude it at every wall point. For a flow at rest in
the Eulerian frame, u = n, and any shift, this exclusion is prior work: [SAR21a] and [SAR21b] for a perfect fluid
without and with Λ, [ABV23] for an isotropic fluid and an anisotropic one with two equal pressures, and [SAR20] and
[ABV23] for dust. Their key step is G_11 = 3G_00, as above. The comparison is in [prior_work.md](prior_work.md) §2.

On the axis, (I1), (I2) and R_01 ≠ 0 make R^a_b restricted to (E_0, E_1) equal to [[−R_00, −R_01], [R_01, −R_00]].
Its eigenvalues are −R_00 ± i R_01, a complex pair. At x − v t = ±1 on the axis, for Alcubierre's profile with σ = 8,
bubble radius R = 1 and v = 1, the script prints 16 ± 4i and 0.

G^a_b differs from R^a_b by a multiple of the identity, so it has the same complex pair, and the Einstein tensor is of
Hawking–Ellis Type IV there. [Rod26] finds Type IV along the whole axis, numerically, for one parameter set. On the
axis ρ_n = S_∥ = 0, so [Le26]'s criterion 2|j| > |ρ_n + S_∥| holds wherever j ≠ 0. Off the axis, [Le26] finds the
sampled wall Type IV except at 1.3 % of points at v = 0.5 and 24.9 % at v = 2.5, while [Rod26] finds the wall annulus
of the equatorial cut Type I ([prior_work.md](prior_work.md), discrepancy 4). The two agree: #64
([he_type.md](he_type.md)) gives the type in closed form, and the Type I set is a thin band at the equator in the
middle of the wall, whose equatorial cut is Rodal's annulus and whose volume is Le's fraction. A scalar (Type I or
II) or a perfect fluid (Type I, or II for a null flow) cannot be the source at a Type IV point, so these maps already
exclude both at their Type IV points ([prior_work.md](prior_work.md) §1–2).

## An accelerating bubble

[Alc94]'s journey starts and ends at rest, v_s = 0, so its bubble accelerates. Section 6 of the script (#61) takes

    ds² = −dt² + (dx + b dt)² + dy² + dz²,    b = −v_s(t) f(r_s),    r_s² = (x − x_s(t))² + y² + z²,

with v_s(t) and x_s(t) independent free functions. Alcubierre's v_s = dx_s/dt is a special case. It also takes a
general shift b(t, x, ρ), not assumed even in ρ, off the axis. The frame is the one above, E_0 = ∂_t − b ∂_x = n,
E_1 = ∂_x, e_ρ and e_ϕ. The Ricci tensor comes from `ricci` with the jets of v_s, x_s and f, or of b, as free symbols
(route D). Time derivatives enter it. For the bubble the frame Ricci tensor contains v_s′ and x_s′, both already in
R_00, and for the general shift it contains ∂_t∂_x b and ∂_t∂_ρ b.

The results, each an exact identity in the jets and the point (**checked symbolically**), are these.

- (I1) and (I2) hold for the general shift (`gen_I1`, `gen_I2`) and for the bubble (`accel_I1`, `accel_I2`).
- R_ρρ = ½(∂_ρb)² for the general shift (`gen_Rrr`), which for the bubble is v_s² f′² ρ²/(2r_s²).
- On the axis, R(E_2, ·) = R(E_3, ·) = 0, R_00 + R_11 = 0 (`accel_axis`) and R_01 = v_s(t) f′(r_s)/r_s
  (`accel_wall`, with R_ρρ off the axis).

#59 derived the same by hand from Cartan's structure equations, for any b(t, x, y, z) and then any axisymmetric one
([review.md](review.md) §1 and gap 6). Its `review_cartan.py` checks (I1) and (I2) for b = G(t, x, y² + z²) with
code that shares only sympy with this script (`axi_I1`, `axi_I2`, `axi_axis`).

The lemma is pointwise, so its application carries over unchanged. At every point where R_ab ≠ 0, no (φ, V) of
either sign, and no perfect fluid, gives the Ricci tensor. For the bubble, R_ab ≠ 0 wherever v_s(t) f′(r_s) ≠ 0,
off the axis through R_ρρ and on it through R_01.

Route D's time derivatives are tested on planted truth. Route E evaluates g and its first and second derivatives
with sympy's own diff of explicit functions, then the Christoffel formula in exact rationals, at a point. It shares
only the metric's definition, `shift_metric`, with route D. For f = 1/(1 + r³), v_s = (1 + t)/(1 + t²),
x_s = t³/3 + t/2 and b = (t² + x) ρ/(2 + ρ² + t x²), the two routes agree in every frame component at 6 exact points
each (`accel_xcheck`, `gen_xcheck`).

The controls fire. Dropping every t-derivative from the jets breaks both cross-checks (`--mutate=accel-static`). A
lapse 1 + b breaks (I1), (I2), the axis identities and the wall curvature (`--mutate=accel-lapse`). Rotating
(e_ρ, e_ϕ) by α, tan(α/2) = 1/20, breaks (I1), (I2) and R_ρρ (`--mutate=accel-rotate`). Doubling ∂_ρb and v_s on the
expected side breaks `gen_Rrr` and `accel_wall` (`--mutate=accel-wall-v2`).

## Explicit profiles, through the control's code path

The same function that passes the control, `scalar_form`, was run at exact points of three kinds, with ε = −1 and
ε = +1 (`alc_points_fail`). Every one of the 46 tests fails at the algebraic stage.

| profile | points | v | jets |
|---|---|---|---|
| Alcubierre's Eq. (6), σ = 8, bubble radius R = 1 | 5 on r_s = 1 (2 off axis, 2 on it, 1 at x = v t) | 1, 1/2, 3 | 50-digit values rounded to 45-digit rationals |
| bump (1 − r²/4)³ for r < 2, 0 beyond | 4 on r_s = 1 | 7/4 | exact |
| none | 4 at r_s = 13/7 and 5/7 | random | random rationals |

For the tanh profile this is **numerical evidence**, since its jets were rounded. The rounded jets are themselves the
jets of some C² profile, and the general result above covers every profile.

## Control: the Ellis–Bronnikov wormhole

Paper III's Eq. (1) at D = 4 (n = 2) is ds² = −dt² + dl² + (l² + a²) dΩ², with κφ′² = 2a²/r⁴ and κ = 8π. Through
`scalar_form` with ε = −1, the test finds R_ll = −2a²/(l² + a²)² as the only nonzero component and λ = 0, so V = 0.
It extracts u_l = −a/(2√π (l² + a²)), which is closed and integrates to φ = −arctan(l/a)/(2√π). Substituted back,
G_ab − 8π T_ab = 0 (`eb_algebraic`, `eb_sign`, `eb_closed`, `eb_potential`, `eb_residual`). This φ′² = a²/(4π r⁴)
matches Paper III's κφ′² = 2a²/r⁴.

The control fails as it should. With ε = +1 the sign check fires (`--mutate=eb-eps+1`). With r(l)² = l² + a²
deformed by a⁴/(l² + a²), the algebraic check fires (`--mutate=eb-deform`).

## Checks and controls

Each check has a control that breaks its input once. The script runs every control after the checks and confirms that
its target check fails. All 22 checks pass and all 14 controls fire (`phantom_source_output.txt`).

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `rho_einstein` | Einstein-tensor ρ_E equals the formula | `rho-v2` | 2v in the formula |
| `rho_3p1` | 3+1 ρ_E equals the formula | `rho-v2` | 2v in the formula |
| `ricci_xcheck` | lab and spherical frame Ricci agree at 24 points | `xcheck-v` | 2v in the spherical one |
| `identity_sph` | (I1), (I2) in comoving spherical coordinates | `alc-drop-w2` | the shift's square dropped from g_tt |
| `identity_lab` | (I1), (I2) in lab Cartesian coordinates | `alc-drop-w2` | the same |
| `identity_axis` | the identities on the axis | `alc-drop-w2` | the same |
| `wall_curvature` | R_ρρ off the axis and R_01 on it | `wall-v2` | 2v in the expected values |
| `alc_points_fail` | 46 exact-point tests all fail algebraically | `alc-planted` | a planted λη − 8π u u at the points |
| `eb_algebraic` | the control's triple eigenvalue | `eb-deform` | r(l)² + a⁴/(l² + a²) |
| `eb_sign` | the control's rank-1 part is −8π u u | `eb-eps+1` | ε = +1 |
| `eb_closed` | du = 0 | `eb-unclosed` | u times 1 + cos²θ |
| `eb_potential` | dφ = u and dV ∧ dφ = 0 | `eb-potential` | V + t/(8π a³) |
| `eb_residual` | G_ab − 8π T_ab = 0 | `eb-residual` | φ times 11/10 |
| `gen_xcheck` | route D for b(t, x, ρ) equals route E at 6 points | `accel-static` | no t-derivatives in the jets |
| `accel_xcheck` | route D for the bubble equals route E at 6 points | `accel-static` | the same |
| `gen_I1` | (I1) for b(t, x, ρ) | `accel-rotate`, `accel-lapse` | (e_ρ, e_ϕ) rotated; lapse 1 + b |
| `gen_I2` | (I2) for b(t, x, ρ) | `accel-rotate`, `accel-lapse` | the same |
| `gen_Rrr` | R_ρρ = ½(∂_ρb)² | `accel-wall-v2`, and both above | 2∂_ρb in the expected value |
| `accel_I1` | (I1) for the accelerating bubble | `accel-rotate`, `accel-lapse` | (e_ρ, e_ϕ) rotated; lapse 1 + b |
| `accel_I2` | (I2) for the accelerating bubble | `accel-rotate`, `accel-lapse` | the same |
| `accel_axis` | the bubble's identities on the axis | `accel-lapse` | lapse 1 + b |
| `accel_wall` | the bubble's R_ρρ off the axis and R_01 on it | `accel-wall-v2`, and both above | 2v_s in the expected values |

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/phantom_source.py

It exits 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, and `--mutate=NAME` one
control.

## What the code proves, and what is assumed

The code proves the identities (I1) and (I2) and the two wall-curvature formulas as exact identities in the point and
the jets (w, w′, w″) or (f, f′, f″). A rational expression is shown to vanish by sympy's `simplify` returning 0
(spherical) or by its numerator's remainder modulo r_s² − (x − v t)² − y² − z² being 0 (Cartesian). That relation is
irreducible, so a zero remainder is exact vanishing. The lemma is by hand. The step from the two to the answer is the
paragraph after the lemma.

For the accelerating bubble, the code proves (I1), (I2), the axis identities and the two wall-curvature formulas as
exact identities in the point and the jets v_s, v_s′, x_s′, f, f′ and f″, the only ones its Ricci tensor contains, by
the same remainder modulo r_s² − (x − x_s)² − y² − z². For b(t, x, ρ) it proves (I1), (I2) and R_ρρ = ½(∂_ρb)² at
ρ > 0, as identities in the jets of b up to second order, by the remainder modulo ρ² − y² − z², which is irreducible
too.

Trusted: sympy 1.14.0 and the routine `ricci`, which the lab, spherical and section 6 routes share. It is tested by
the 3+1 route's agreement on ρ_E, by the cross-check between the two coordinate systems, by the control's R_ll and by
route E. The 3+1 route shares only the chain-rule derivative `D` with the lab route. Route E shares only the metric's
definition `shift_metric` with route D, and neither `ricci` nor the jet derivative `JetD`.

Assumed: one real, minimally coupled scalar with T_ab as above, G_ab = 8π T_ab, and a profile with f(r_s) C² on ℝ³.
Up to "An accelerating bubble", v is constant. There the shift b is C², with unit lapse and flat slices.

Not covered: a bubble whose path is not a straight line, a shift that is not axisymmetric about the direction of
motion, a lapse other than 1 or curved slices (with the lapse 1 + b, (I1) and (I2) fail), two or more scalars,
non-minimal coupling, and matter whose Ricci tensor is not of the form λ g_ab + s u_a u_b. A normal plus a phantom
scalar give R_ab − λ g_ab of rank 2 with mixed signs, which the lemma does not exclude. #65 excludes that pair on
every open set of the wall, for constant v and the accelerating bubble ([two_scalars.md](two_scalars.md), *checked
symbolically*). #70 shows that (I1) and (I2) hold for every shift b(t, x, σ(y, z)) with unit lapse, in the frame
adapted to ∇⊥σ, so the lemma covers some shifts that are not axisymmetric ([nonaxisymmetric.md](nonaxisymmetric.md)).

The answer rests on computer algebra, so it should get an independent review before it is relied on.

## Context

[CDK24] in [`../critical_ghost.md`](../critical_ghost.md) evolved a collapsing warp bubble sourced by an NEC-violating
fluid, not a scalar.

## References

- [Alc94] M. Alcubierre, "The warp drive: hyper-fast travel within general relativity", Class. Quantum Grav. 11, L73
  (1994), arXiv:gr-qc/0009013, doi:10.1088/0264-9381/11/5/001. Checked on 2026-11-13 against the arXiv abstract page
  and PDF (Eqs. (6), (8), (10) and (19) read there) and IOP's article page. For #61, the unnumbered definitions of
  v_s(t) = dx_s/dt and r_s(t) after Eqs. (2)–(5), and the journey that starts and ends with v_s = 0, were read in the
  arXiv source on 2026-11-13.

The entries below are those of [`../literature.md`](../literature.md), "Warp drives", whose journal data were checked
against INSPIRE and arXiv identifiers against arXiv for #60.

- [ABV23] Abellán, Bolívar, Vasilev, arXiv:2302.13826, GRG 55, 60 (2023).
- [Le26] Le, arXiv:2602.18023, CQG 43, 215013 (2026).
- [LV04a] Lobo, Visser, arXiv:gr-qc/0406083, CQG 21, 5871 (2004).
- [Rod26] Rodal, arXiv:2512.18008, GRG 58, 1 (2026).
- [SAR20] Santos-Pereira, Abreu, Ribeiro, arXiv:2008.06560, EPJC 80, 786 (2020).
- [SAR21a] Santos-Pereira, Abreu, Ribeiro, arXiv:2101.11467, EPJC 81, 133 (2021).
- [SAR21b] Santos-Pereira, Abreu, Ribeiro, arXiv:2108.10960, EPJ Plus 136, 902 (2021).
- [SSV22] Santiago, Schuster, Visser, arXiv:2105.03079, PRD 105, 064038 (2022).
