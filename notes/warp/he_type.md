# The Hawking–Ellis type of the Alcubierre wall

November 2026, issue [#64](https://github.com/nmoinvaz/negative-mass/issues/64). Script [`he_type.py`](he_type.py), its
stdout [`he_type_output.txt`](he_type_output.txt).

Status words: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails), **numerical
evidence**.

## Answer

The type of G^a_b at a point of the Alcubierre spacetime is the sign of one polynomial. In the frame of
[phantom_source.md](phantom_source.md), e_ϕ decouples by (I1), and the discriminant D of the cubic of the
(E_0, E_1, e_ρ) block decides. The point is Type I where D > 0 and Type IV where D < 0, and on D = 0 a rank test picks
Type I, II or III (**proved**). For every profile f and speed v,

    D = v⁶ (p₆ + p₈ v² + p₁₀ v⁴ + p₁₂ v⁶),    p₆ v⁶ = −4|j|⁶ < 0,    p₁₂ ≥ 0,

with each p explicit in f, f′, f″, r_s and θ (**checked symbolically**). So the wall is Type IV at small speed and turns
Type I as v grows. On the axis D < 0 wherever f′ ≠ 0, and on the equator the point is Type I exactly where
v r_s f′² > |(r_s f′)′|.

Le and Rodal describe the same map (**numerical evidence**, quadrature of the exact indicator). At Le's parameters and
v = 0.5 the Type I fraction of the wall is 1.15 %, against Le's 1.3 %, whose grid ladder gives 1.1, 1.3 and 1.0 %. At
Rodal's parameters the axis is Type IV at every radius, and the equatorial Type I annulus is r_s ∈ [4.946, 5.066], of
width 0.120 against Rodal's "∼ 2∆" = 0.14. That annulus is the equatorial cut of a thin Type I band that holds 9.4 % of
the wall at Rodal's parameters and 1.2 % at Le's. Discrepancy 4 of [prior_work.md](prior_work.md) compared a planar
cut with a volume fraction, at different speeds and wall widths, so it is not a disagreement.

## Setup

The metric, the frame (E_0 = n, E_1 = ∂_x, e_ρ, e_ϕ), the profile f(r_s) and the frame Ricci components are those of
[phantom_source.md](phantom_source.md), with w = v(1 − f). Below, a = f′(r_s), b = f″(r_s), F = 1 − f(r_s),
c = cos θ and s = sin θ, and r stands for r_s.

G^a_b = R^a_b − ½R δ^a_b, so G^a_b and R^a_b have the same eigenvectors and the same type. By (I1), e_ϕ is an
eigenvector of R^a_b with eigenvalue 0, spacelike and orthogonal to the rest. So the type is that of the mixed block
M = ηR on (E_0, E_1, e_ρ), with η = diag(−1, 1, 1). Its characteristic polynomial is λ³ + c₂λ² + c₁λ + c₀, with c₂ = −tr M,
c₁ the sum of the principal 2×2 minors, c₀ = −det M, and

    D = 18 c₂c₁c₀ − 4 c₂³c₀ + c₂²c₁² − 4 c₁³ − 27 c₀².

**Lemma (proved).** Let M be self-adjoint with respect to η on 2+1 Minkowski space.

- If D > 0, the three eigenvalues are real and distinct. Their eigenvectors are mutually η-orthogonal and span, so none
  is null and one is timelike: Type I.
- If D < 0, M has a complex pair: Type IV.
- If D = 0, a real root λ repeats with multiplicity m. If rank(M − λ) = 3 − m, M is diagonalizable, η is nondegenerate
  on each eigenspace, and an orthonormal eigenbasis contains one timelike vector: Type I. Otherwise M has a Jordan block
  at λ. With (M − λ)u = 0 and (M − λ)w = u, η(u, u) = η(u, (M − λ)w) = η((M − λ)u, w) = 0, so the eigenvector is null:
  Type II for a 2×2 block (m = 2, rank 2, or m = 3, rank 1) and Type III for a 3×3 block (m = 3, rank 2). ∎

Adding the decoupled e_ϕ does not change the type. The same cubic discriminant is in [Le26] Corollary 1, for a block
(n, ĵ, t̂) whose remaining axis ẑ decouples, where Le notes that the decoupling is a hypothesis. For Alcubierre, (I1)
is that hypothesis, with ẑ = e_ϕ, at every point.

`he_type_exact` implements the lemma, and the planted tensors test it (`planted`). A perfect fluid with
(μ, p) = (3/2, 1/5), dust, an anisotropic fluid, two canonical Type IV tensors, null dust with and without a pressure,
[Le26]'s Type III example ℓm + mℓ and λη are each written in a frame boosted twice and rotated, with a decoupled fourth
axis. All nine get their type, the fluids through D = 0 and the rank test. numpy's eigenvalues agree on the four with
D ≠ 0.

## The closed form

The momentum entries R_01 and R_0ρ are linear in v, the others quadratic, and D is even in the momentum, so D has only
the powers v⁶ to v¹² (`structure`, he_type_output.txt §2).

- **Pure momentum.** p₆ = −(J²)³/(16r⁶) with J² = s²(a + rb)² + 4c²a². This is the discriminant of the momentum part
  alone, p₆v⁶ = −4|j|⁶ with |j|² = R_01² + R_0ρ² = v²J²/(4r²).
- **Static stress.** p₁₂ = s⁴F²Q₁²Q₂/(16r⁶), with Q₁ = −c²F(a − rb)² − 4ra³ and
  Q₂ = a²(F + ra)² + c²F(a + rb)[F(rb − a) − 2ra²]. It equals χ(λ₀)²[(S₁₁ − S_ρρ)² + 4S_1ρ²] ≥ 0, the discriminant of
  the v² part, which is block diagonal, with λ₀ = −R_00, S its spatial block and χ the characteristic polynomial of S.
- **Mixed.** p₈ and p₁₀ have 160 and 347 terms. They are printed in he_type_output.txt §2.
- **Equator**, c = 0. D = v⁶[r²v²a⁴ − (a + rb)²] Q_eq²/(16r⁶), with Q_eq = 4v²a²F(F + ra) + (a + rb)². Where
  Q_eq ≠ 0 the point is Type I exactly when v r f′² > |(r f′)′|. There j is along E_1, a principal axis of the stress,
  and the rule is [Le26]'s momentum criterion 2|j| < |ρ_n + S_∥|.
- **Axis**, c = ±1. D = −4a²v⁶Q_ax²/r⁶, with Q_ax = a² + r²v²(Fb − a²)². So the axis is Type IV wherever f′ ≠ 0, with
  the pair −R_00 ± iR_01 = v²(f′² − (1 − f)f″) ± i v f′/r_s of phantom_source.md, 16 ± 4i at its sample point.

What follows for every profile, from these forms.

- As v → 0 the wall is Type IV wherever J ≠ 0, which is everywhere f′ ≠ 0 except on the equatorial circle where
  (r f′)′ = 0. That circle is Type I at every speed, where Q_eq ≠ 0.
- As v → ∞ a point is Type I wherever p₁₂ > 0, which excludes the axis, the centre and the surfaces Q₁ = 0 and Q₂ = 0.
- In between, the cubic p₁₂V³ + p₁₀V² + p₈V + p₆ in V = v² has one or three positive roots. The point is Type IV below
  the lowest, V₁. With three roots V₁ < V₂ < V₃ it is Type I on (V₁, V₂), Type IV again on (V₂, V₃) and Type I above
  V₃. For the tanh profile, three roots occur on PCT3LE % of Le's wall and PCT3RO % of Rodal's (he_type_output.txt §7,
  **numerical evidence**), and route N confirms the four bands at the widest window.

## Le's setup

[Le26] Sec. 3 and Table 3: Alcubierre's profile (their Eq. 46) with R_b = 1 and σ = 8, the slice t = 0, the wall mask
f ∈ [0.1, 0.9] (r_s ∈ [0.8627, 1.1373]), proper-volume weighted, which on these flat slices is the coordinate volume.
Their grid is graded toward the wall on ±3R_b with clustering strength a⋆ = 2 and N = 100, and their types come from a
four-dimensional eigensolver with a relative tolerance of 10⁻¹⁰.

| v | Type I %, closed form | uniform 100³ | uniform 200³ | [Le26] |
|---|---|---|---|---|
| 0.5 | 1.151 | 1.29 | 1.24 | 1.3 (Table 3), 1.1, 1.3, 1.0 at N = 80, 100, 120 (Table 12) |
| 1.0 | 5.60 | 6.48 | 5.56 | 5.7 |
| 2.0 | 17.28 | 17.31 | 17.56 | 16.7 |
| 2.5 | 22.63 | 22.35 | 22.70 | 24.9 (Sec. 3.1) |

The closed-form column is a midpoint quadrature over (r_s, θ) at 2000², which agrees with 1000² to 0.002 percentage
points. The uniform columns count the cells of a uniform cubic grid on ±3 whose centres lie in the mask.

At v = 0.5 the closed form is 0.15 percentage points below Le's 1.3 %, inside the 2-point tolerance, and close to the
mean 1.13 % of Le's ladder. At v = 2 and 2.5 it differs from Le by +0.6 and −2.3 points. Le's ladder spread was
measured at v = 0.5 only, and the paper does not give the formula of its graded grid, so its grid estimate at the higher
speeds could not be replicated. That is the one thing in Le's setup not pinned down here.

Le's momentum criterion, Type IV where 2|j| > |ρ_n + S_∥|, agrees with the exact type on 96.9 % of the wall's volume at
v = 0.5. Le reports 97.5 % of wall points, by an unweighted count on a different grid (N = 60, clustering 1.2). The
agreement falls to 66 % at v = 2.5, so the cubic is needed off the equator.

## Rodal's setup

[Rod26] Sec. 4: ρ = 5, σ = 4 and v = 1, plotted on the band r_s ∈ [ρ − 3/σ, ρ + 3/σ]. Rodal finds the axis
"Type IV for all radii", and on the equatorial cut Type IV inside the bubble and in the far field, with Type I in "an
annulus of width ∼ 2∆ about the wall", ∆ ≈ 0.014ρ = 0.07.

| cut or region | closed form |
|---|---|
| axis, θ = 0 and π | Type IV at all 150 001 radii of the band |
| equator, θ = π/2 | Type I for r_s ∈ [4.9458, 5.0660], width 0.120 = 0.86 × 2∆; Type IV inside and outside |
| Le's mask f ∈ [0.1, 0.9] | Type I on 9.35 % of the volume |
| Rodal's band | Type I on 12.48 % of the volume |

Rodal's paper gives no volume fraction to compare with.

## Discrepancy 4 settled

The Type I set is a thin band that hugs the equator near the middle of the wall, where |(r f′)′| is small. It crosses
the equatorial plane in Rodal's annulus and fills 1.2 % of the wall's volume at Le's parameters. Both papers are
consistent with the closed form at their own parameters. The fraction grows with v, from 1.2 % to 5.6 % between
v = 0.5 and v = 1 at Le's σR_b = 8, and is 9.4 % at Rodal's σρ = 20 and v = 1 (he_type_output.txt §5–7).

## Where the exclusion of #57 applies, by type

A minimally coupled scalar has a stress tensor of Type I or II, and a perfect fluid of Type I, or II for a null flow
([prior_work.md](prior_work.md) §1–2). So the map divides the wall in two.

- At Type IV points, 98.8 % of Le's wall at v = 0.5 and 90.6 % of Rodal's at v = 1, the type alone excludes a scalar and
  a perfect fluid, with any flow.
- At Type I points, the band above, and on the surface D = 0, where the points are Type II or degenerate, the type
  allows both. There only the lemma of [phantom_source.md](phantom_source.md) excludes them. On the equator the edges of
  the annulus have 2|j| = |ρ_n + S_∥| with j ≠ 0, which is Type II by [Le26] Eq. 12. Off the equator the type on D = 0
  was not checked.

## Checks and controls

The script runs every control after the checks and confirms that its target checks fail. CHECKSUMMARY

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `planted` | nine tensors of known type, I, II, III and IV, get their type from D and the rank test | `disc-sign`, `no-rank` | D → −D; every D = 0 tensor labelled II |
| `decouple` | R(e_ϕ, ·) = 0 and (I2) on route C's components with w = v(1 − f) | `alc-drop-w2` | the shift's square dropped from g_tt in route C |
| `structure` | the powers v⁶ to v¹², p₆, p₁₂ and the equator and axis factorizations | `alc-drop-w2` | the same |
| `jets` | the tanh jets against mpmath's numerical derivatives at 8 radii | `jets-f2` | f″ × (1 + 1/1000) |
| `numN` | route N's G(n, n)/8π equals Alcubierre's Eq. (19), and its axis eigenvalues equal −R_00 ± iR_01 | `num-drop-w2` | the shift's square dropped from g_tt in route N |
| `closed_vs_numeric` | the sign of D gives route N's eigenvalue type at every grid point | `alc-drop-w2`, `num-drop-w2` | as above |
| `le_fraction` | the Type I fraction at Le's v = 0.5 is within 2 points of 1.3 % | `le-speed` | v = 1 |
| `rodal_cuts` | the axis is Type IV, and the equatorial Type I set is one annulus about ρ of width within a factor 2 of 2∆ | `rodal-speed` | v = 1/4 |

Under `alc-drop-w2` the map changes, from 1.15 % Type I at Le's v = 0.5 to DROPW2 % (he_type_output.txt, Controls).

The grids of `closed_vs_numeric` are Le's wall at 240² points in (r_s, θ) for each of v = 0.5, 1, 2 and 2.5, the cells
of a uniform 100³ grid on ±3 inside Le's mask at v = 0.5, Rodal's band at 240² and Rodal's axis at 602 points, 304 778
in all. Each point gets a random azimuth for route N. The three points nearest D = 0 on each grid, and any point where
the two routes disagree in float64, are rechecked at 50 digits. None disagreed.

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/he_type.py

It exits 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, and `--mutate=NAME` one
control. It takes about a minute without the controls.

## What the code proves, and what is assumed

Checked symbolically: the decoupling, the powers of v in D, the forms of p₆ and p₁₂ and the equator and axis
factorizations, each as an exact identity in (v, f, f′, f″, r, cos θ), with sin² θ reduced to 1 − cos² θ. By hand: the
lemma, which turns the sign of D and the rank test into the type.

Numerical evidence: the agreement with route N at 304 778 points, the fractions, which are midpoint quadratures of the
exact indicator D > 0, and the root counts of §7.

Trusted: sympy 1.14.0 and phantom_source.py's `ricci` and `sph_tetrad_ricci` (route C), numpy 2.5.4's `einsum` and
`eigvals` and mpmath 1.3.0 (route N). Route N writes the metric and its derivatives by hand and shares with route C only
the metric's definition and `tanh_jets`. `tanh_jets` is checked against mpmath, and route N against Alcubierre's
Eq. (19) and the axis eigenvalues. `planted` tests `cubic`, `disc` and `he_type_exact`, which the Alcubierre map uses.

Assumed: constant v, as in phantom_source.md up to its accelerating bubble. The closed form holds for every C²
profile, and the numbers use the tanh profile on the slice t = 0, which any other t translates.

Not covered: the accelerating bubble of phantom_source.md. (I1) holds there, so the cubic still decides, but its D was
not computed. The type on D = 0 off the equator. Le's graded grid.

## References

- [Le26] H. Le, arXiv:2602.18023, CQG 43, 215013 (2026). Read in the arXiv version on 2026-11-13: Sec. 2.2 (Eq. 12 and
  Corollary 1), Sec. 3 (Tables 1–3 and the wall mask), Sec. 3.2, the type-label comparison after Lemma 3, Table 12 and
  Eq. (46).
- [Rod26] J. Rodal, arXiv:2512.18008, GRG 58, 1 (2026). Read in the arXiv version on 2026-11-13: Eq. (14), Sec. 4.1.2
  and 4.1.3 (the classification protocol and tolerances), the parameters (ρ, σ, v/c) = (5 m, 4 m⁻¹, 1) and the
  Alcubierre cuts of Sec. 4.2.

Both entries are those of [`../literature.md`](../literature.md), "Warp drives", checked against INSPIRE for #60.
