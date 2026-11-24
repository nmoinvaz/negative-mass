# The Hawking–Ellis type of the accelerating bubble

November 2026, issue [#73](https://github.com/nmoinvaz/negative-mass/issues/73). Script [`he_accel.py`](he_accel.py),
its stdout [`he_accel_output.txt`](he_accel_output.txt). It extends [he_type.md](he_type.md) (#64, reviewed in
[review_he.md](review_he.md), #71) to the accelerating bubble of [phantom_source.md](phantom_source.md) (#61).

Status words: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails), **numerical
evidence**.

## Answer

The metric is ds² = −dt² + (dx + b dt)² + dy² + dz², b = −v_s(t) f(r_s), r_s² = (x − x_s(t))² + y² + z², with v_s(t) and
x_s(t) independent free functions. Write v = v_s, α = v_s′, u = x_s′ (Alcubierre's journey is u = v), f₁ = f′, f₂ = f″,
and c = cos θ, s = sin θ for the angle from the direction of motion. In the frame (E_0 = n, E_1 = ∂_x, e_ρ, e_ϕ) the
frame Ricci tensor contains only v, α, u, f, f₁ and f₂. Neither v_s″ nor x_s″ enters (**checked symbolically**), where
the issue expected x_s″.

The acceleration enters through two derivatives of the shift and nothing else. With T_x = ∂_t∂_x b and T_ρ = ∂_t∂_ρ b,

    T_x = −α c f₁ + v u (s² f₁/r + c² f₂),    T_ρ = −α s f₁ + v u c s (f₂ − f₁/r),

the block of the Ricci tensor on (E_0, E_1, e_ρ) is

    R_00 = T_x − A,   R_11 = A − q − T_x,   R_ρρ = q/2,   R_1ρ = B − T_ρ/2,   R_01 = j₁,   R_0ρ = j_ρ,

    A = (v²/2r)[2f (s² f₁ + r c² f₂) + r (1 + c²) f₁²],      B = (v² c s/2r)(2r f₁² + r f f₂ − f f₁),
    q = v² s² f₁²,    j₁ = v[(1 + c²) f₁ + r s² f₂]/(2r),    j_ρ = v c s (f₁ − r f₂)/(2r),

and R(e_ϕ, ·) = 0 (**checked symbolically**, `time_lemma`, `decouple`). The momentum entries do not depend on α or u.

The type is that of the 3×3 block, so it is decided by the discriminant D of its cubic and, on D = 0, by the rank test,
as in #64. With p = R_ρρ + R_00 = q/2 + T_x − A and h = R_1ρ, the block shifted by −R_00 I is

    [[0, −j₁, −j_ρ], [j₁, −q, h], [j_ρ, h, p]],

with cubic λ³ + c₂λ² + c₁λ + c₀, c₂ = q − p, c₁ = j₁² + j_ρ² − pq − h², c₀ = 2 j₁ j_ρ h − p j₁² + q j_ρ², and

    D = 18 c₂c₁c₀ − 4 c₂³c₀ + c₂²c₁² − 4 c₁³ − 27 c₀²       (52 terms, quartic in p).

The point is Type I where D > 0 and Type IV where D < 0 (**proved**, the lemma of he_type.md). So the type at a point
depends on the motion through (v, T_x, T_ρ) only, that is, through p and h. At constant speed (α = 0, u = v) the block
is #64's entry by entry and its type map is #64's at 288 000 grid points (**checked symbolically**, **numerical
evidence**).

What follows, each from these forms.

- **At v = 0 every point is Type I** (**proved**), where #64's map at v → 0 is Type IV wherever j ≠ 0. There j₁ = j_ρ =
  0, the shifted block is [[0, 0, 0], [0, 0, h], [0, h, p]] and D = α⁶ f₁⁶ s⁴/16 ≥ 0. The block is block diagonal and
  symmetric, so diagonalizable, which is Type I also where D = 0.
- **The axis is Type IV for every α and u** wherever v f₁ ≠ 0 (**proved**, **checked symbolically**). On c = ±1, j_ρ = q
  = h = 0 and D = −4 R_01² (R_00² + R_01²)², with R_01 = v f₁/r.
- **At large |α| the wall is Type I** wherever f₁ s ≠ 0 (**checked symbolically**). The α⁶ coefficient of D is f₁⁶ s⁴/16
  and D has degree 6 in α, so for fixed v, u and point, D > 0 once |α| is large enough.
- **D has only the weights 6, 8, 10, 12** under (v, u, α) → (λv, λu, λ²α), as #64's has the powers v⁶ to v¹²
  (**proved**, **checked symbolically**). The weight 6 part is −4|j|⁶ and the weight 12 part is (pq + h²)²[(p + q)² +
  4h²] ≥ 0.
- **Parity.** D(v, −α, u; −c) = D(v, α, u; c) (**checked symbolically**). So the type at (θ, α) is the type at (π − θ,
  −α), and every wall fraction is even in α.

How the map changes (**numerical evidence**, tables below). At Le's parameters and v = 0.5 the Type I fraction of the
wall is 1.15 % at α = 0, 11.99 % at |α| = 2 and 93.9 % at |α| = 32. Along Alcubierre's journey at v = 0.10 and α = ±0.39
it is 10.6 %, against 0.045 % at the same speed without the acceleration terms. The acceleration can also lower it: at v
= 2 it falls from 17.28 % to 15.19 % at |α| = 8 before it rises. On Rodal's equator the single Type I annulus of #64
(width 0.120) is three intervals at |α| = 0.5 and 1, and one of width 0.206 at |α| = 4.

## Setup

The metric, the frame and the identities (I1), (I2) are those of [phantom_source.md](phantom_source.md), "An
accelerating bubble". At z = 0 and y > 0 the Cartesian frame vectors (∂_y, ∂_z) are (e_ρ, e_ϕ), and the axisymmetry
about the x axis makes this no restriction. G^a_b and R^a_b have the same type, and (I1) makes e_ϕ an eigenvector of
R^a_b with eigenvalue 0, spacelike and orthogonal to the rest. So the type is that of M = ηR on (E_0, E_1, e_ρ), η =
diag(−1, 1, 1), exactly as in he_type.md. The discriminant is invariant under M → M − κI, so the shift by M_00 = −R_00
changes neither D nor the type. (I2) gives R_11 = −R_00 − 2R_ρρ, which makes the shifted block's (0, 0) entry 0 and its
(1, 1) entry −q.

The entries come from `ps.shift_tetrad_ricci` (route D of phantom_source.py: Christoffel and Ricci in the lab
coordinates, the jets of v_s, x_s and f as free symbols), set to z = 0 and written in (c, s) with s² reduced to 1 − c².

## The entries and the time lemma

**Lemma E (checked symbolically).** The block of route D at z = 0 equals the six entries above, as an exact identity in
(v, α, u, f, f₁, f₂, r, c, s) modulo s² = 1 − c², and R(e_ϕ, ·) = 0. It agrees with route E (exact rationals, sympy's
own derivatives of the explicit v_s = (3 + 2t)/(2 + t²), x_s = t³/5 − t²/3 + t/2, f = (2 + r²)/(3 + r⁴)) at 8 points
(`entries_E`).

**Lemma T (checked symbolically).** R(v, u, α) − R(v, v, 0) = ΔT_x (E_00 − E_11) − ½ΔT_ρ (E_1ρ + E_ρ1), where ΔT = T −
T|_{u = v, α = 0}. So the acceleration moves p by ΔT_x and h by −ΔT_ρ/2 and changes nothing else in the shifted block.
The entries' dependence on time is only through ∂_t∂_i b because R_ρρ = ½(∂_ρ b)² has no time derivative (K_ρρ = 0) and
R_0i is the momentum constraint.

## The discriminant, in closed form

**Lemma D (checked symbolically).** D of the block, through he_type's `mixed_block`, `cubic` and `disc`, equals D(j₁,
j_ρ, q, h, p) above. This is first an identity in five free symbols (the entries R_00, R_01, R_0ρ, R_ρρ, R_1ρ with R_11
= −R_00 − 2R_ρρ), then for the bubble's entries.

**Lemma W (proved, checked symbolically).** Under j → λj and q, h, p → λ²(q, h, p), D has only the powers λ⁶, λ⁸, λ¹⁰,
λ¹². D has degree 6 in the entries of M, and conjugating by diag(−1, 1, 1) flips the signs of the momentum entries only,
so D is even in them. A monomial with m momentum entries and 6 − m stress entries has weight m + 2(6 − m) = 12 − m, with
m = 0, 2, 4, 6. The weight 6 part is the discriminant of the momentum part, with cubic λ(λ² + |j|²), so −4|j|⁶. The
weight 12 part is the discriminant of the block diagonal matrix with cubic λ(λ² + (q − p)λ − pq − h²), that is (pq +
h²)²[(p + q)² + 4h²].

For the bubble, 64 r⁶ D is a polynomial in (v, α, u) with 50 monomials v^i α^j u^k. Their weights i + k + 2j are 6, 8,
10, 12, every u comes with a v (k ≤ i), the highest powers are v¹², α⁶ and u⁶, and the α⁶ coefficient is f₁⁶ s⁴/16
(`structure`). The weight 8 and 10 parts in the five symbols are printed in he_accel_output.txt §5.

## Lemmas stated for Lean (#135)

Hypotheses throughout: f is C² on ℝ³ as a function of r_s, v_s and x_s are C¹, r_s > 0, and the point has s ≥ 0, c² + s²
= 1.

1. **(E)** The six entries above, with R(e_ϕ, ·) = 0 and (I2). Status: checked symbolically; the Ricci tensor is
   sympy's, tested by route E.
2. **(S)** For a real symmetric 3×3 block with R_11 = −R_00 − 2R_ρρ, D(ηR) = D(j₁, j_ρ, q, h, p). Status: checked
   symbolically in five free symbols; by hand it is M → M − M_00 I and the characteristic polynomial of the shifted
   matrix.
3. **(T)** The classification lemma of he_type.md: D > 0 Type I, D < 0 Type IV, D = 0 by the rank test. Status: proved,
   as there.
4. **(V0)** If v = 0 the point is Type I. Status: proved.
5. **(Ax)** If c = ±1 and v f₁ ≠ 0 the point is Type IV, for every α and u. Status: proved.
6. **(W)** The weights of D. Status: proved.
7. **(P)** D(v, −α, u; −c) = D(v, α, u; c). Status: checked symbolically; by hand, (α, c) → (−α, −c) changes the signs
   of j_ρ and h only.
8. **(L)** The α⁶ coefficient of D is f₁⁶ s⁴/16. Status: checked symbolically. It gives: if f₁ s ≠ 0, there is α₀ with D
   > 0 for |α| ≥ α₀.
9. **(C)** At α = 0, u = v the entries are #64's. Status: checked symbolically, entry by entry against `ht.route_C()`.

Assumed in 1: the Ricci tensor of the metric with sympy's `ricci`, which `phantom_source.md` tests with route E and
which this script tests again at 8 points and numerically with route N′ (below).

## The map

**Le's wall** (σ = 8, R_b = 1, f ∈ [0.1, 0.9]), Type I % of the proper volume, u = v, polar 1000²:

| v \ α | −32 | −16 | −8 | −4 | −2 | −1 | −0.5 | 0 | 0.5 | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 93.92 | 87.31 | 71.32 | 34.12 | 11.99 | 2.98 | 1.14 | 1.15 | 1.14 | 2.98 | 11.99 | 34.12 | 71.32 | 87.31 | 93.92 |
| 1.0 | 87.20 | 69.88 | 31.72 | 10.55 | 5.17 | 5.44 | 5.56 | 5.60 | 5.56 | 5.44 | 5.17 | 10.55 | 31.72 | 69.88 | 87.20 |
| 2.0 | 58.46 | 25.96 | 15.19 | 16.70 | 17.13 | 17.25 | 17.27 | 17.28 | 17.27 | 17.25 | 17.13 | 16.70 | 15.19 | 25.96 | 58.46 |

The α = 0 column is #64's 1.151, 5.60 and 17.28 %. **Rodal's band** (σ = 4, ρ = 5), v = 1, u = v:

| α | −32 | −16 | −8 | −4 | −2 | −1 | −0.5 | 0 | 0.5 | 1 | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Type I % | 97.33 | 92.88 | 49.52 | 26.24 | 17.02 | 12.70 | 12.41 | 12.48 | 12.41 | 12.70 | 17.02 | 26.24 | 49.52 | 92.88 | 97.33 |

The speed of the bubble's centre matters on its own. At v = 0.5 and α = 0, the Type I fraction of Le's wall is 2.50,
1.03, 1.15, 6.56 and 19.42 % for u = 0, 0.25, 0.5, 1 and 2, where u = v = 0.5 is #64's case. That is the (v u) term of
T_x and T_ρ.

**Along Alcubierre's journey** v_s = sech²(2t), x_s = tanh(2t)/2 (P1), Le's wall, with and without the acceleration
terms:

| t | v | α | Type I %, full | Type I %, α = 0 and u = v |
|---|---|---|---|---|
| ±0.9 | 0.1036 | ∓0.3922 | 10.641 | 0.045 |
| ±0.6 | 0.3050 | ∓1.0171 | 8.597 | 0.377 |
| ±0.3 | 0.7116 | ∓1.5286 | 3.383 | 2.704 |
| 0 | 1 | 0 | 5.599 | 5.599 |

**Rodal's equator** (θ = π/2), v = 1, u = v. The Type I set on the equator has these intervals in r_s, from the sign
changes of D.

| α | intervals of Type I |
|---|---|
| 0 | [4.9458, 5.0660] |
| ±0.5 | [4.9438, 4.9656], [4.9945, 5.0191], [5.0509, 5.0683] |
| ±1 | [4.9378, 4.9510], [5.0032, 5.0090], [5.0668, 5.0747] |
| ±2 | [4.9175, 4.9209], [4.9706, 5.0386] |
| ±4 | [4.8981, 5.1036] |

At the equator j_ρ = 0 and h = α f₁/2. For α ≠ 0 the cubic no longer splits off the root p, which is why the annulus
changes. For α = 0 and any u it does: the cubic is (λ − p)(λ² + qλ + j₁²), so D = (q² − 4j₁²)(p² + pq + j₁²)² and #64's
criterion v r f′² > |(r f′)′| holds at the equator for every u (**proved**, by hand from the cubic; the script does not
check it).

The type labels with the acceleration terms dropped differ from the full ones at 1800 to 9583 of the 57 600 points of
each wall, at every slice with v ≠ 0 (`drop_detected`; the table is in he_accel_output.txt §9). At the instant of rest,
P2 at t = 0, the full block is Type I at all points and the block without the acceleration terms is the zero tensor, for
which D = 0.

## Checks and controls

The script runs every control after the checks and confirms that its target checks fail. All 13 checks pass and all 13
controls fire, and so does `phantom_source.py --mutate=accel-static`, which fails `gen_xcheck` and `accel_xcheck`
(`he_accel_output.txt`).

| check | what it tests | control `--mutate=` | the control's change | failed besides the target |
|---|---|---|---|---|
| `planted_cls` | #64's nine planted tensors get their type from D and the rank test | `disc-sign`, `no-rank` | D → −D; every D = 0 tensor labelled II | `disc-sign`: `structure`, `const_limit`, `axis`, `rest`, `closed_vs_numeric`; `no-rank`: `rest` |
| `routeN_planted` | route N′ at constant speed equals `ht.einstein_N` at 40 points | `num-drop-w2`, `num-static` | b² dropped from g_tt; every t-derivative of b zeroed | `routeN_vs_E`, `closed_vs_numeric` |
| `routeN_vs_E` | route N′ at accelerating points has route E's characteristic polynomial | the same two | the same | the same |
| `entries_E` | route D's block at z = 0 equals route E at 8 exact points | `accel-static`, `alc-drop-w2`, `accel-lapse` | phantom_source's static jets; b² dropped from g_tt in route D; lapse 1 + b | `accel-static`: `time_lemma`, `structure`, `const_limit`, `closed_vs_numeric`, `drop_detected` |
| `decouple` | (I1) and (I2) for the block | `alc-drop-w2`, `accel-lapse` | b² dropped; lapse 1 + b | see below |
| `time_lemma` | the six entries in closed form, and Lemma T | `time-lemma-sign` | sign of the ΔT_ρ term in the claim flipped | none |
| `structure` | Lemmas S and W, the powers of v, α, u, the α⁶ coefficient | `weight-break` | α counted as weight 1 | none |
| `parity` | D(−α, −c) = D(α, c) | `parity-break` | α flipped, c not | none |
| `const_limit` | at (α, u) = (0, v), #64's entries, map and fractions | `alc-drop-w2`, `accel-static` | as above | see above |
| `axis` | D on the axis, Type IV for every α and u | `axis-flip`, `alc-drop-w2` | the claim with +4 | none for `axis-flip` |
| `rest` | v = 0: R_0i = 0, D ≥ 0, Type I at exact points | `rest-v` | the slice at t = 0.05, v = 1/10 at the exact points | none |
| `closed_vs_numeric` | sign of D against route N′'s eigenvalues at 1 158 020 points | `drop-accel`, `alc-drop-w2`, `num-drop-w2`, `num-static`, `accel-static`, `disc-sign` | D evaluated at α = 0, u = v; the others as above | `drop-accel`: `drop_detected` |
| `drop_detected` | dropping the acceleration terms changes the map | `drop-accel` | the full D is the dropped one | none |

`accel-lapse` stops after section 3, where `entries_E` and `decouple` are decided, because the later symbolic sections
do not finish on its rational entries. `alc-drop-w2` and `accel-lapse` break (I2), so `structure` and `parity` fail for
them at once, by that rule and not by content. The controls cost 46 s (`accel-lapse`) to 544 s (`alc-drop-w2`) each, and
the whole script takes about 40 min on one core of a Raspberry Pi 5, 192 s of it without the controls.

`closed_vs_numeric` has 86 float64 disagreements, all on Rodal's axis at P2, t = 0, where D = 0 exactly. Each is
rechecked at 50 digits, and none is unresolved. The grids of `closed_vs_numeric` are Le's wall and Rodal's band at 240²
points in (r_s, θ) and Rodal's axis at 602 points, on 10 slices of two motions, each point with a random azimuth for
route N′.

| profile | v_s(t) | x_s(t) | t |
|---|---|---|---|
| P1, Alcubierre's journey (u = v) | sech²(2t) | tanh(2t)/2 | −0.6, −0.3, 0, 0.3, 0.6 |
| P2, x_s free, v_s through 0 | tanh(2t) | (3/10) ln cosh 2t − (2/15) cos 3t + t/4 | the same |

P1 has α = 0 and u = v at t = 0 (the constant-speed case v = 1), α of both signs and v > 0. P2 has v < 0, v = 0 at t = 0
with α = 2 and u = 1/4, and v > 0, and u ≠ v at every t. So the 10 slices give 10 distinct (v, α, u) with every sign
pattern of v and α.

## What the code proves, and what is assumed

Checked symbolically: the entries (Lemma E) and Lemma T as identities in (v, α, u, f, f₁, f₂, r, c) modulo s² = 1 − c²,
the five-symbol identities (Lemmas S and W), the weights, the powers and the α⁶ coefficient of D for the bubble, parity,
the axis form, and the equality of the block at (α, u) = (0, v) with #64's route C. Exact points: route E at 8 points,
and `he_type_exact` at 12 axis points and 32 points at v = 0, including D = 0, where the rank test gives Type I. By
hand: Lemmas T, V0, Ax and W.

Numerical evidence: the agreement of the sign of D with route N′'s eigenvalues at 1 158 020 grid points, and the
fractions.

Route N′ is a numerical Einstein tensor for an explicit b(t, x, y, z). Sympy differentiates b, so there are no jets and
no frame, and the tensor algebra is the einsum of he_type's `einstein_N`, which `routeN_planted` checks against it at
constant speed (relative error 4.6e-15) and `routeN_vs_E` against exact rationals at accelerating points (4.7e-15).
Route N′ shares with route D only the form of the metric and the tanh profile. At points where the closed-form D and the
eigenvalues are in doubt in float64 both are recomputed at 50 digits. A point is D = 0 only by an exact factor: at v =
0, D = α⁶ f₁⁶ s⁴/16 and s⁴ is evaluated in rationals.

Run 1 of the grids failed `closed_vs_numeric`: 516 points on Rodal's axis at v = 0, where D = 0 exactly and the check
read rounding noise as a sign. Run 2 changed how D = 0 is decided, raised the recheck cap from 100 to 1000 and added
`parity`, all written in the raw README before run 2 (I73-rpi5). The failed run's log is `step3.log` there.

Trusted: sympy 1.14.0, `ps.ricci`, `ht.cubic`, `ht.disc` and `ht.he_type_exact`, numpy's `einsum` and `eigvals`, mpmath
1.3.0.

Assumed: unit lapse and flat slices, a straight path, v_s and x_s free functions, an axisymmetric C² profile f(r_s).

Not covered: the surface D = 0 off v = 0 (#72's results are for constant speed; here the rank test is run at exact
points only), the time evolution of the map along a journey (slices only), a lapse other than 1, a path that is not
straight, and Le's graded grid.

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/he_accel.py

It exits 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, `--mutate=NAME` one
control, `--stop=7` the symbolic sections only, and `--quick` the subset of the grids.
