# Can a phantom scalar source the Alcubierre metric?

November 2026, issue [#57](https://github.com/nmoinvaz/negative-mass/issues/57). Script
[`phantom_source.py`](phantom_source.py), its stdout [`phantom_source_output.txt`](phantom_source_output.txt).

Status words: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails), **numerical
evidence**, *(unverified)* for a reference detail not confirmed.

## Answer

No, for every profile and either sign of the kinetic term (**checked symbolically**). At every point where
f′(r_s) ≠ 0, Alcubierre's Ricci tensor is not of the form R_ab = 8π(ε ∂_aφ ∂_bφ + V g_ab), for ε = −1, for ε = +1
and for any value of V. The test fails at its first, algebraic stage, since no λ makes R_ab − λ g_ab of rank ≤ 1. The
sign, closure and potential conditions are never reached.

Every C² profile with f(0) = 1 that decays at infinity has f′ ≠ 0 on an open interval of r_s. So no (φ, V) sources
the metric, and the conditions fail at every point of the bubble wall, the set where f′(r_s) ≠ 0.

The obstruction is two identities of the Ricci tensor that hold for every profile, and a five-line lemma.

## Setup

The metric is Alcubierre's Eq. (8) [Alc94], with constant v = v_s,

    ds² = −dt² + (dx − v f(r_s) dt)² + dy² + dz²,    r_s² = (x − v t)² + y² + z².

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

Two identities follow, for every C² profile and every v.

- **(I1)** R(e_ϕ, ·) = 0.
- **(I2)** R_00 + R_11 + 2 R_ρρ = 0.

They are checked in comoving spherical coordinates for a general w(r) (`identity_sph`), in lab Cartesian coordinates
for a general f with the jets f, f′, f″ as free symbols (`identity_lab`), and on the axis y = z = 0, where they read
R(E_2, ·) = R(E_3, ·) = 0 and R_00 + R_11 = 0 (`identity_axis`).

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
flow, has a Ricci tensor of that form, so the same identities exclude it at every wall point. The closest prior work
found by #60 treats only dust and perfect fluids at rest in the Eulerian frame ([prior_work.md](prior_work.md)).

On the axis, (I1), (I2) and R_01 ≠ 0 make R^a_b restricted to (E_0, E_1) equal to [[−R_00, −R_01], [R_01, −R_00]].
Its eigenvalues are −R_00 ± i R_01, a complex pair. At x − v t = ±1 on the axis, for Alcubierre's profile with σ = 8,
bubble radius R = 1 and v = 1, the script prints 16 ± 4i and 0.

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
its target check fails. All 13 checks pass and all 10 controls fire (`phantom_source_output.txt`).

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

Trusted: sympy 1.14.0 and the routine `ricci`, which the lab and spherical routes share. It is tested by the 3+1
route's agreement on ρ_E, by the cross-check between the two coordinate systems and by the control's R_ll. The 3+1 route
shares only the chain-rule derivative `D` with the lab route.

Assumed: one real, minimally coupled scalar with T_ab as above, G_ab = 8π T_ab, a C² profile and a constant v.

Not covered: a time-dependent v_s(t), two or more scalars, non-minimal coupling, and matter whose Ricci tensor is
not of the form λ g_ab + s u_a u_b. A normal plus a
phantom scalar give R_ab − λ g_ab of rank 2 with mixed signs, which the lemma does not exclude.

The answer rests on computer algebra, so it should get an independent review before it is relied on.

## Context

[CDK24] in [`../critical_ghost.md`](../critical_ghost.md) evolved a collapsing warp bubble sourced by an NEC-violating
fluid, not a scalar.

## References

- [Alc94] M. Alcubierre, "The warp drive: hyper-fast travel within general relativity", Class. Quantum Grav. 11, L73
  (1994), arXiv:gr-qc/0009013, doi:10.1088/0264-9381/11/5/001. Checked on 2026-11-13 against the arXiv abstract page
  and PDF (Eqs. (6), (8), (10) and (19) read there) and IOP's article page.
