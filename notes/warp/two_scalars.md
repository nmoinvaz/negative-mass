# Can a normal plus a phantom scalar source the Alcubierre wall?

November 2026, issue [#65](https://github.com/nmoinvaz/negative-mass/issues/65). Script
[`two_scalars.py`](two_scalars.py), its stdout [`two_scalars_output.txt`](two_scalars_output.txt). It builds on
[`phantom_source.md`](phantom_source.md) (#57, #61), whose Ricci tensor it imports, and on its review
[`review.md`](review.md) (#59).

Status words: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails), **numerical
evidence**.

## Answer

No (**checked symbolically**, with a lemma **proved** by hand). For every C² profile and any potential V(φ₁, φ₂), no
normal scalar φ₁ and phantom scalar φ₂ source Alcubierre's metric on an open set that meets the wall, the set where
f′(r_s) ≠ 0. The same holds for #61's accelerating bubble on every open set that meets v_s(t) f′(r_s) ≠ 0.

The test fails at its algebraic stage, but not at every point. Pointwise, R_ab = λ g_ab + a_a a_b − b_a b_b has real
solutions exactly where D = det R₃ = 0, with R₃ the frame Ricci tensor on (n, ∂_x, e_ρ), and there λ = 0. That set
contains the axis of motion. There a normal and a phantom gradient do give #57's complex eigenvalue pair, as #65
anticipated.

The set has empty interior in the wall. At each r_s where f′ ≠ 0, D is a nonzero cubic in cos²θ, unless 2w = r w′ and
r w″ = −w′ there, with w = v(1 − f). Those two relations hold on no interval of r_s. So every open set that meets the
wall contains points where the field equations fail, and the closure and potential stages are never reached on an open
set.

Where the algebraic stage does pass on an open set, closure decides. The shifts b(ρ) with b″ + b′/ρ = ±b′², which are
not warp bubbles, pass it with λ = 0 and fail closure. The control, an Ellis–Bronnikov wormhole with an added normal
scalar, passes all four stages.

## Setup

Two minimally coupled scalars, φ₁ normal and φ₂ phantom, with a potential V(φ₁, φ₂), have

    T_ab = ∂_aφ₁ ∂_bφ₁ − ∂_aφ₂ ∂_bφ₂ − g_ab (½(∇φ₁)² − ½(∇φ₂)² + V).

G_ab = 8π T_ab is equivalent to R_ab = λ g_ab + a_a a_b − b_a b_b, with λ = 8πV, a = √(8π) dφ₁ and b = √(8π) dφ₂.
A scalar of one sign is the case a = 0 or b = 0, so #57's test is contained in this one. The test has four stages.

1. **Algebraic.** Some real λ makes R − λg a real form with at most one positive and one negative square. The
   candidates are the real roots of the gcd of the 3 × 3 minors of R − λg, each a cubic in λ. Then a and b are fixed
   up to a sign each and a field-space boost (a, b) → (cosh χ a + sinh χ b, sinh χ a + cosh χ b), with χ any
   function.
2. **Closure.** Some χ makes both covectors closed. With k = a + b and m = a − b, this asks for e^χ k and e^−χ m
   closed, that is ω = dχ with ω ∧ k = −dk and ω ∧ m = dm. When k ∧ m ≠ 0, contracting with vectors K and M that
   satisfy k(K) = m(M) = 1 and k(M) = m(K) = 0 leaves one candidate, ω = dm(K, M) k + ι_K dk (**proved**). The stage
   checks that it solves both equations and that dω = 0.
3. **Potential.** dλ ∧ a ∧ b = 0, so that V is locally a function of (φ₁, φ₂) where dφ₁ ∧ dφ₂ ≠ 0. The condition
   is invariant under the boost.
4. **Residual.** G_ab − 8π T_ab = 0 with the fields substituted back.

With a and b closed and independent, both Klein–Gordon equations follow from the contracted Bianchi identity, and so
does stage 3 (**proved**). ∇^a G_ab = 0 with the form above gives dλ = (∇·a) a − (∇·b) b. So for a metric, stage 3
cannot fail once stage 2 passes, and its control breaks λ directly. The same holds for one scalar, as in #57.

## The lemma

**Lemma (proved).** Let R_ab be real and symmetric on an orthonormal frame (e_0, e_1, e_2, e_3), with R_a3 = 0 for
every a and R_00 + R_11 + 2R_22 = 0. Let R₃ be its block on (e_0, e_1, e_2). Then R_ab = λ η_ab + a_a a_b − b_a b_b,
with λ real and a, b real covectors, if and only if λ = 0 and det R₃ = 0.

*Proof.* M = R − λη splits as N ⊕ (−λ), with N = R₃ − λ diag(−1, 1, 1). By the identity,
N_00 + N_11 + 2N_22 = −2λ. The form a a − b b has at most one positive and one negative square.

- If λ > 0, the e_3 direction is a negative square, so N is positive semidefinite. Then N_00 + N_11 + 2N_22 ≥ 0,
  which contradicts −2λ < 0.
- If λ < 0, N is negative semidefinite, and the same sum is ≤ 0, which contradicts −2λ > 0.
- So λ = 0 and M = R₃ ⊕ 0. The weights (1, 1, 2) are positive, so a nonzero R₃ with a vanishing weighted trace is
  indefinite. It then has at most one square of each sign exactly when its rank is ≤ 2, that is when det R₃ = 0. ∎

It applies with (e_2, e_3) = (e_ρ, e_ϕ) off the axis and (E_2, E_3) on it, where #57's identities (I1) and (I2) hold.
On the axis R₃ has a zero row, so det R₃ = 0 and the form holds with λ = 0.

The proof uses only the signature of the kinetic form. So it covers any two-field sigma model whose target metric is
Lorentzian, with any potential, since that metric can be diagonalized at each point. It does not cover two scalars of
one sign. For two normal scalars M is positive semidefinite, so λ ≤ 0, and λ < 0 with N ⪰ 0 of rank ≤ 1 is allowed by
the weighted trace. Two phantoms are the mirror case.

`lemma_samples` tests the statement on 24 random rational tensors that satisfy the hypotheses, each in a random
rational Lorentz frame. All 12 with det R₃ = 0 pass with λ = 0, and all 12 with det R₃ ≠ 0 fail.

## The determinant for a constant velocity

From the frame Ricci tensor of phantom_source.md, with w = v(1 − f) and c = cos θ, 16 r³ det R₃ is a cubic in c²
(`two_scalars_output.txt` §4). Its extreme coefficients are

    c⁰:  2 w′² r (w′⁴ r² − 4 w² w′² − (w′ + r w″)²)
    c⁶: −2 (2 r w w″ − 2 w w′ + r w′²) (r w w″ − w w′ + r w′²)²

and the other two are printed there. It vanishes at c² = 1, the axis (`det_axis`).

With w′ r ≠ 0, the four coefficients generate the same ideal as 2w − r w′ and (w′ + r w″)² (`det_jets`, inclusion
both ways with Gröbner bases). So at a radius in the wall, D(r, ·) vanishes for every θ exactly where 2w = r w′ and
r w″ = −w′ there.

These relations hold on no interval with w′ ≠ 0 (**proved**, and `no_open_set`). Differentiating 2w = r w′ gives
w′ = r w″, which with r w″ = −w′ is w′ = 0.

Hence the zero set of D in the wall has empty interior (**proved**). An open subset of the wall contains a box in
(r, θ). On it D vanishes for all θ at each r of the box, since a cubic in c² that vanishes on an interval of θ is
zero, so the relations hold on an interval of r. That forces w′ = 0, which contradicts f′ ≠ 0. The field equations
must hold on an open set, so no (φ₁, φ₂, V) sources the metric on an open set that meets the wall.

## Exact points, through the control's code path

`alc_points` runs the same function that passes the control, `two_scalar_form`, on the lab-frame Ricci tensor of
phantom_source.md at #57's sample points, with fresh random jets in the last row. Every outcome agrees with D, which
is evaluated from the spherical route.

| profile | points | v | result |
|---|---|---|---|
| Alcubierre's Eq. (6), σ = 8, bubble radius R = 1 | 5 on r_s = 1, 2 of them on the axis | 1, 1/2, 3 | 9 off the axis fail, D ≠ 0; 6 on it pass with λ = 0 |
| bump (1 − r²/4)³ for r < 2 | 4 on r_s = 1, 1 on the axis | 7/4 | 3 fail, D ≠ 0; 1 passes with λ = 0 |
| none | 4 at r_s = 13/7 and 5/7 | random | 4 fail, D ≠ 0 |

On the axis at x − vt = 1, for the tanh profile with v = 1, the split is a = (0, √17, 0, 0) and b = (−4, −1, 0, 0) in
the frame (E_0, E_1, E_2, E_3). The normal gradient is spacelike and the phantom's timelike, and η(a a − b b) has the
eigenvalues 16 ± 4i and 0. So a remainder of rank 2 with mixed signs does produce #57's complex pair there. The axis
has empty interior, so this gives no source.

## An accelerating bubble

For #61's b = −v_s(t) f(r_s), with v_s ≠ 0, route D's frame R₃ depends on the jets only through the comoving shift
w̃ = x_s′ − v_s f, its derivatives in r_s and A = v_s′/v_s (`acc_jets`). At A = 0, its 16 r³ D is section 4's
polynomial exactly, so the lab route with jets and the comoving spherical route agree (`acc_static_limit`).

A brings in the odd powers c, c³ and c⁵, each with the factor A r w̃′. With A w̃′ r ≠ 0 the seven coefficients have
no common zero, since their Gröbner basis is {1} (`acc_det_jets`). So at an instant with v_s′ ≠ 0, D(r_s, ·) vanishes
for every θ at no r_s in the wall. At an instant with v_s′ = 0, section 4's relations apply to w̃, and they hold on no
interval. Either way, at fixed t the zero set of D has empty interior in the wall. An open set of spacetime that meets
the wall contains such a slice, so the answer carries over.

## A shift b(ρ): where the algebraic stage passes, closure decides

The lemma holds for every axisymmetric shift b(t, x, ρ) with unit lapse and flat slices, since (I1) and (I2) do (#61).
It is the bubble's shape that keeps D from vanishing on open sets.

For a shift b(ρ) alone, R₃ = [[−h, Q, 0], [Q, −h, 0], [0, 0, h]] with h = ½b′² and Q = −½(b″ + b′/ρ), so
D = h(h² − Q²) (`brho_det`). Where b′ ≠ 0, D vanishes on an open set exactly where b″ + b′/ρ = ε b′² with ε = ±1.
Its solutions are b = −ε ln|C − ε ln ρ| + b₀ (**proved**, a Bernoulli equation for b′). Up to x → −x,
x → x + b₀ t and a homothety, they are b = −ln(−ln ρ) on 0 < ρ < 1 and b = ln ln ρ on ρ > 1.

Both pass the algebraic stage with λ = 0 on an open set. Up to sign, the normal gradient is a = (b′/√2) dρ and the
phantom's is b = (b′/√2)((1 + εb) dt + ε dx), which is null. Both fail closure (`brho_closed`). By hand (**proved**),
k = a + b and m = a − b are proportional to ν± = dρ ± ((1 + εb) dt + ε dx), and ν± ∧ dν± = b′ dx ∧ dρ ∧ dt ≠ 0. So no
factor e^±χ makes them closed, and no normal plus phantom pair sources a shift b(ρ) on an open set where b′ ≠ 0.

## Control: Ellis–Bronnikov with an added normal scalar

The metric is ultrastatic,

    ds² = −dt² + dl² + r(l)² dΩ²,    r(l)² = (1 − β cos 2κl) / (2κ²),    0 < β < 1, κ > 0,

with φ₁ = κt/(2√π), φ₂ = ∓arctan(√((1 + β)/(1 − β)) tan κl)/(2√π) for |κl| < π/2, and V = κ²/(4π). The potential is
a constant, a cosmological constant Λ = 8πV = 2κ². The field equations give V = φ̇₁², φ₂′ = C/r² and
r′² = 1 − κ²r² − a²/r², with a² = 4πC² = (1 − β²)/(4κ²). The throat is at l = 0, where r² = (1 − β)/(2κ²), and r is
periodic in l. As κ → 0 at fixed a, with 1 − β² = 4κ²a², r² → l² + a² and φ₂ → ∓arctan(l/a)/(2√π), the
Ellis–Bronnikov control of phantom_source.md. The metric solves G = 8πT exactly (`eb2_residual`). Whether it is in
print was not searched.

The test passes all four stages. It finds λ = 2κ², a = √2 κ dt, b ∝ dl, ω = 0 and dλ = 0 (`eb2_algebraic`,
`eb2_closed`, `eb2_potential`, `eb2_residual`).

The controls fire. With cos 2κl → cosh 2κl and β > 1, both squares of R − λg are negative, as for two phantom
scalars, and the algebraic stage fails (`--mutate=eb2-ads`). Adding
(1 − β) cos(4κl)/(20κ²) to r² makes λ depend on l, and closure fails (`--mutate=eb2-deform`).

## Checks and controls

Each check has a control that breaks its input once. The script runs every control after the checks and confirms that
its target checks fail. All 16 checks pass and all 15 controls fire (`two_scalars_output.txt`).

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `planted_pass` | 12 random λη + a a − b b, transformed, pass with λ recovered | `planted-sign` | λη + a a + b b |
| `planted_fail` | 12 random forms with signs (+,+), (−,−), (+,−,+), (+,−,−) fail | `planted-fail-sign` | one + and one − square |
| `lemma_samples` | the lemma on 24 random tensors with (I1), (I2) | `lemma-weight` | R_00 − R_11 + 2R_22 = 0 |
| `eb2_algebraic` | the control's form, one + and one − square | `eb2-ads` | the two-phantom metric |
| `eb2_closed` | ω exists and dω = 0 | `eb2-deform`, `eb2-unclosed` | r² deformed; b times 1 + cos²θ |
| `eb2_potential` | dλ ∧ a ∧ b = 0 | `eb2-potential` | λ + κ² cos θ / 10 |
| `eb2_residual` | dφ_i = u_i and G_ab − 8πT_ab = 0 | `eb2-residual` | phantom gradient times 11/10 |
| `det_axis` | D = 0 on the axis, even in c | `alc-drop-w2` | the shift's square dropped from g_tt |
| `det_jets` | the coefficients' ideal is that of 2w − r w′ and (w′ + r w″)² | `alc-drop-w2` | the same |
| `no_open_set` | the relations hold on no interval | `ode-sign` | r w″ = +w′ |
| `alc_points` | 23 exact points agree with D | `alc-planted`, `alc-drop-w2` | a planted a a − b b off the axis; route C changed |
| `acc_jets` | the bubble's R₃ depends on w̃, w̃′, w̃″, A only | `acc-static` | no t-derivatives in the jets |
| `acc_static_limit` | at A = 0, route D's D is route C's | `alc-drop-w2` | route C's metric changed |
| `acc_det_jets` | with A ≠ 0 the coefficients have no common zero | `acc-drop-odd` | only even powers of c kept |
| `brho_det` | D = h(h² − Q²) for b(ρ), and the two solutions | `brho-h2` | h doubled in the expectation |
| `brho_closed` | both b(ρ) solutions pass the algebra and fail closure | `brho-drop-b2` | b² dropped from g_tt |

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/two_scalars.py

It exits 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, and `--mutate=NAME` one
control. On win2, with every core busy with other jobs, the checks took 20 s and the whole run with its controls
773 s, at 301 MB peak resident. Three hash seeds (`PYTHONHASHSEED` 1, 7 and 42) give the same output.

## What the code proves, and what is assumed

The code proves, as exact identities or exact ideal computations: the coefficients of D for a constant velocity and
for the accelerating bubble, that the first ideal equals the one of 2w − r w′ and (w′ + r w″)², that the second ideal
is the unit ideal when A w̃′ r ≠ 0, that the bubble's R₃ depends on the reduced jets only, D = h(h² − Q²) for b(ρ),
and that the two b(ρ) solutions fail closure. It also proves that the control solves the field equations. The outcomes
at the 23 sample points are exact for their jets. The tanh profile's jets were rounded to 45-digit rationals, so for
that profile they are numerical evidence, as in #57.

By hand: the lemma, the uniqueness of ω, the Bianchi step behind stage 3, the step from the relations to w′ = 0, the
empty interior of the zero set of D, the solutions of the Bernoulli equation, and the Frobenius obstruction for b(ρ).

Numerical evidence only: the signs that the control's splitting needs where sympy's assumptions do not decide them,
which are sampled at 8 random points of 0 < β < 1/2 (1 < β < 3 for `eb2-ads`).

Trusted: sympy 1.14.0 (`groebner`, `factor`, `simplify`) and the routines imported from phantom_source.py, `ricci`,
`sph_tetrad_ricci`, `lab_tetrad_ricci` and `shift_tetrad_ricci`, tested there and reviewed in #59. `acc_static_limit`
compares two of them, routes C and D, on D. The two-scalar test, the determinant analysis and the control are new code,
tested on the planted forms and the control.

Assumed: two real, minimally coupled scalars with a Lorentzian kinetic form, any potential V(φ₁, φ₂), G_ab = 8π T_ab,
and a profile with f(r_s) C² on ℝ³. The bubble moves at a constant v ≠ 0, or as in #61 with v_s ≠ 0.

## Not covered

- Two scalars of one sign. The lemma's argument allows λ ≠ 0 there.
- Three or more scalars. With λ = 0, R₃ ⊕ 0 has rank 3 at generic wall points and is indefinite, so a form with
  three squares of mixed signs passes the algebraic stage there. Closure would decide.
- Shifts b(t, x, ρ) other than the bubble and b(ρ). The lemma applies to them, but D may vanish on open sets.
- Non-minimal coupling, kinetic terms that are not quadratic, and scalars with other matter.

## Prior work

The record of #60 ([prior_work.md](prior_work.md) §1) found no paper that asks whether a scalar field, normal or
phantom, sources the Alcubierre metric. No further search was made for two scalars. [BD20] fit two- and three-fluid
mixtures to the Alcubierre stress tensor numerically and say one fluid gives an overdetermined system
([prior_work.md](prior_work.md) §2, not re-read here). A minimally coupled scalar is a perfect fluid only where its
gradient is timelike, so their mixtures and these pairs overlap without being the same class. The Type IV maps of
[Rod26] and [Le26], which include the axis, are compared in [prior_work.md](prior_work.md) §1.

## References

Entries as in [`../literature.md`](../literature.md), "Warp drives".

- [BD20] Béatrix-Drouhet, arXiv:2012.09941 (2020, preprint).
- [Le26] Le, arXiv:2602.18023, CQG 43, 215013 (2026).
- [Rod26] Rodal, arXiv:2512.18008, GRG 58, 1 (2026).
