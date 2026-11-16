# Is the linearised point-mass problem ill-posed?

## Answer

Yes, in Hadamard's sense, for the Dotti–Gleiser (DG) problem with its unique self-adjoint dynamics at the singularity.
For every t > 0 the solution operator is unbounded from 𝓗^s to 𝓗^{s′} for all real s, s′ (angular Sobolev scale, radial L²(dx)).
The exact modes of dg_proof Theorem 1 give data of norm ε_l(1+L_l)^{s/2} → 0 whose solution at a fixed time has norm 1, because the
rate ω_l = (l−1)l(l+1)(l+2)/12 grows like l⁴/12. Each l is well-posed, the problem is well-posed in the weighted spaces with weight
e^{ω_l T}, and some smooth data have no solution at any t > 0. This is *proved* for the DG variable Ψ̂ (§1–2), it lives entirely at the
singularity (data supported at x ≥ x₁ > 0 obey a uniform energy bound until t = x₁, §4), and a truncated star with an l-independent
Robin parameter κ ≥ 0 has finitely many unstable l, a finite top rate Ω(R, κ), and a solution operator bounded by cosh(Ωt) in every
𝓗^s (§5). That is an ordinary instability, not ill-posedness.

Units |M| = 1, f = 1 + 2/r, x = r − 2 ln(1 + r/2), c = (l+2)(l−1), L_l = l(l+1) = c + 2, k = c(c+2)/6, **ω_l = k/2 = L_l(L_l−2)/12**.
Status words follow AGENTS.md (*proved*, *checked symbolically*, *certified*, *conditional*, *conjecture*, *numerical evidence*), plus
*argued* for reasoning with a named gap. Scripts are in this directory, each writes `X_output.txt`, and every check has a mutation
(§8). Commands run from this directory as `OMP_NUM_THREADS=1 python3 -u X.py > X_output.txt` (Python 3.14, `requirements.txt`).

## 0. Statements

| # | statement | norms | status |
|---|---|---|---|
| A | ‖C_l(t)‖ = cosh ω_l t and ‖Σ_l(t)‖ = sinh(ω_l t)/ω_l, exactly | L²(dx) → L²(dx) | *proved* |
| B1 | no bounded solution operator, any t > 0, any s, s′ | 𝓗^s ⊕ 𝓗^s → 𝓗^{s′} | *proved* |
| B2 | data → 0 with ‖u(t₀)‖ = 1, from the exact modes | 𝓗^s (all s) → 𝓗⁰ | *proved*, checked to l = 32 (*numerical evidence* for the PDE evolution) |
| B3 | smooth data with no solution for any t > 0 | 𝓗^{−s} (all s) | *proved* |
| B4 | well-posed iff the data weight is ≳ cosh ω_l t | weights e^{ω_l T} → 𝓗⁰ | *proved* |
| C | uniform energy bound for data in {x ≥ x₁}, t < x₁ | energy norm ℰ | *proved* given Theorem 3 and certified r_neg → 0, plus finite speed (standard, not re-proved) |
| D | after t = x₁ the growth is e^{ω_l(t−x₁)}, onset → x₁ | L²(dx) | *numerical evidence*; the overlap lower bound is *argued* |
| E | star, κ ∈ [0, ∞]: ‖C_l‖ ≤ cosh Ωt for all l, Ω < ∞ | 𝓗^s → 𝓗^s, all s | *proved* given Theorem 2 and the certified bound c·r_* ≤ 9.1152; checked l ≤ 60 |
| F | the same for the metric perturbation | metric Sobolev | *argued* (§7) |

## 1. Setting and hypotheses

For each (l, m) the even-parity perturbation, after DG's intertwiner, is a function Ψ̂_{lm}(t, x) on x ∈ (0, ∞) with
∂ₜ²Ψ̂ = −Ĥ_l Ψ̂, where Ĥ_l = −∂ₓ² + V̂_l = A*A is the operator of `phase2/dg_proof/RESULTS.md`.

- H1 (DG 2008, *not re-proved*, dg_proof "Caveat on scope"). This equation is the linearised Cauchy problem, smooth metric data
  corresponding to Ψ̂ data. Everything below is a statement about Ψ̂.
- H2 (dg_proof L2/L3, *proved*). Ĥ_l is essentially self-adjoint on C_c^∞(0, ∞) in L²((0, ∞), dx), limit point at both ends. No
  boundary condition is chosen at the singularity, and the solution in C([0, T]; L²) is unique. The excluded branch is x^{−1/2}.
- H3 (dg_proof Theorem 1, *proved*). The part of the spectrum below 0 is the single simple eigenvalue −ω_l², eigenfunction
  ψ_l ∝ Ψ̂_l = 2^k e^{−ω_l x} r⁴Π_l(r)/χ_l(r), unit-normalised in L²(dx). The closed form with Π a polynomial is derived and checked in
  `check_modes.py` (no pole at r_s).

Norms. ‖·‖ is the norm of L²((0, ∞), dx), dx = dr/f. For u = Σ u_{lm}(x)Y_{lm}, ‖u‖_s² = Σ_{lm}(1+L_l)^s ‖u_{lm}‖², the space 𝓗^s,
s ∈ ℝ. Data are (u₀, v₀) = (Ψ̂, ∂ₜΨ̂)(0). Hadamard well-posedness means existence, uniqueness and ‖u(t)‖_{s′} ≤ C(t)(‖u₀‖_s + ‖v₀‖_s).

## 2. Theorem A and the Hadamard conclusion

With the spectral calculus of Ĥ_l put C_l(t) = cos(t√Ĥ_l) and Σ_l(t) = sin(t√Ĥ_l)/√Ĥ_l, so the solution is u(t) = C_l(t)u₀ + Σ_l(t)v₀.

**Theorem A** (*proved*). For l ≥ 2 and t ≥ 0, in L²(dx),
**‖C_l(t)‖ = cosh(ω_l t)** and **‖Σ_l(t)‖ = sinh(ω_l t)/ω_l**.

*Proof.* By H3, Ĥ_l ≥ 0 on the orthogonal complement of ψ_l. There |cos t√λ| ≤ 1 and |sin(t√λ)/√λ| ≤ t, and cosh ≥ 1 and
sinh(y)/y ≥ 1, so the bounds are attained on ψ_l, where the multipliers are cosh ω_l t and sinh(ω_l t)/ω_l. ∎

**Corollary B** (*proved*; full problem, t > 0, s, s′ ∈ ℝ).

- B1 (no continuous dependence). For data (ψ_l Y_{l0}, 0) the ratio ‖u(t)‖_{s′}/‖u₀‖_s equals cosh(ω_l t)(1+L_l)^{(s′−s)/2}, unbounded
  in l. So ‖S(t)‖_{𝓗^s⊕𝓗^s → 𝓗^{s′}} = ∞. The same holds for any norm in which ψ_lY_{l0} has size at most polynomial in l, for
  example D(|Ĥ|^{n/2}) in the radial variable, whose norm of ψ_l is exactly (1+ω_l²)^{n/2}. A polynomial loss in translating
  metric data to Ψ̂ data cannot repair it.
- B2 (the sequence). Fix t₀ > 0 and put ε_l = 1/cosh(ω_l t₀), u₀^{(l)} = ε_l ψ_l Y_{l0}, v₀ = 0. Then ‖u₀^{(l)}‖_s = ε_l(1+L_l)^{s/2} → 0 for
  every s, since ln ε_l ≈ −l⁴t₀/12, while u^{(l)}(t₀) = ψ_l Y_{l0} has ‖u^{(l)}(t₀)‖₀ = 1 for every l. For t > t₀ it is cosh(ω_l t)/cosh(ω_l t₀) → ∞.
- B3 (no existence). Let u₀ = Σ_{l≥2} e^{−l²} ψ_l Y_{l0} (in every 𝓗^s), v₀ = 0. If u ∈ C([0, t]; 𝓗^{−s}) were a weak solution, the
  coefficient ⟨u(τ), ψ_lY_{l0}⟩ would satisfy y″ = ω_l² y (ψ_l ∈ D(Ĥ_l)), so it equals e^{−l²}cosh(ω_l τ), which is not O((1+l)^p).
  So no solution exists at any t > 0, in any Sobolev space.
- B4 (where it is well-posed). Let ‖(u₀,v₀)‖_{W_T}² = Σ e^{2ω_lT}(‖u_{0,lm}‖² + ‖v_{0,lm}‖²). For 0 ≤ t ≤ T, ‖u(t)‖₀ ≤ ‖u₀‖_{W_T} + ¼‖v₀‖_{W_T}
  (ω_l ≥ 2). Conversely a diagonal weight w with ‖u(t)‖₀ ≤ C‖u₀‖_w needs w_l ≥ cosh(ω_l t)/C. A weight e^{a l^p} works iff p > 4, or
  p = 4 and a > t/12. Data must decay faster than e^{−l⁴T/12} in l, which is faster than analytic.
- Uniqueness holds (H2), so of Hadamard's three conditions the second holds and the first and third fail.

Hadamard's own example (the Cauchy problem for the Laplace equation) has growth e^{lt}, linear in the frequency. Here it is l⁴.

## 3. The sequence, advanced in time

`hadamard_sequence.py` advances u_tt = −Ĥu on a grid (P1 Galerkin, lumped mass, leapfrog, `evolve.py`) with the potential V̂_l from χ and
data ε_l ψ_l, 40 points per 1/ω_l, t₀ = 0.01. It reports G_l(t) = ‖u(t)‖/‖u(0)‖ against cosh(ω_l t). Thresholds were fixed from
`planted.py` before it was run (§8). They are |ln G − ln cosh ωt| ≤ 0.01 + 5·10⁻⁴ωt and a Rayleigh quotient within 1.5·10⁻³.

| l | ω_l | ω_l t₀ | log₁₀ ‖data‖ | ln G(t₀) | ln cosh ω_l t₀ | diff |
|---|---|---|---|---|---|---|
| 2 | 2 | 0.02 | −0.000 | 0.000 | 0.000 | +0.0000 |
| 5 | 70 | 0.70 | −0.099 | 0.227 | 0.227 | +0.0001 |
| 8 | 420 | 4.20 | −1.523 | 3.508 | 3.507 | +0.0009 |
| 12 | 2002 | 20.0 | −8.394 | 19.331 | 19.327 | +0.0045 |
| 16 | 6120 | 61.2 | −26.28 | 60.521 | 60.507 | +0.0137 |
| 20 | 14630 | 146.3 | −63.24 | 145.640 | 145.607 | +0.0328 |
| 24 | 29900 | 299.0 | −129.6 | 298.374 | 298.307 | +0.0671 |
| 28 | 54810 | 548.1 | −237.7 | 547.530 | 547.407 | +0.1229 |
| 32 | 92752 | 927.5 | −402.5 | 927.035 | 926.827 | +0.2080 |

The data norm falls strictly with l (e^{−927} at l = 32) and ‖u(t₀)‖ stays 1. The error in ln G is 2.2·10⁻⁴ of the exponent at every l = 2..32, the
same relative error as the planted toy (2.0·10⁻⁴). I have not explained why it is so uniform in l. Over T = 20/ω_l (growth e^{19.3}) the error is
+0.0045 for every l = 2..12, so ω_l itself is tested. The exact mode has Rayleigh quotient −ω_l² to 4.5·10⁻⁴ (planted 4.0·10⁻⁴) and the grid norm
equals the mp norm to 3·10⁻¹⁰. This is *numerical evidence* for the PDE evolution. The statement B2 is *proved* without it.

## 4. Data away from the singularity

```
  t
  ▲   growth e^{ω_l (t − x₁)} at and after the arrival
  │
  x₁ ●
  │    ╲
  │      ╲    ray x + t = x₁ from the data edge
  │        ╲
  └─────────●──────────────► x
  r = 0     x₁   data supported in x ≥ x₁, uniform energy bound while t < x₁
```

**Proposition C** (*proved*, given Theorem 3 and the certified r_neg → 0 of dg_proof §4, and finite speed of propagation for the sector wave
equation away from the singularity, which is standard and not re-proved here). Let x₁ > 0 and 0 < θ < 1, and put x₀ = (1−θ)x₁. There is an L with
x(r_neg(l)) ≤ x₀ for l ≥ L. Take l ≥ L and data supported in {x ≥ x₁}. Define the energy on x ≥ x₀ by
E_l(t) = ½∫_{x₀}^∞ (u_t² + u_x² + V̂_l u²) dx, with norm ℰ_l = (2E_l)^{1/2}. For t ≤ θx₁ the solution vanishes near x₀, so the boundary term
−u_t u_x(x₀) is 0 and E_l(t) = E_l(0). Since V̂_l ≥ 0 on x ≥ x₀ (Theorem 3), E_l(t) ≥ ½(‖u_t‖² + ‖u_x‖²). Hence there is no growth in the
energy norm ℰ_l, uniformly in l ≥ L, for t ≤ θx₁. The finitely many sectors l < L are bounded by Theorem A.

So the ill-posedness is the singularity's. Data have to be signalled to r = 0 first.

### Growth after x₁

This is *numerical evidence* (`causal.py`). For g ≥ 0 supported at x ≥ x₁ the overlap a_l = ⟨ψ_l, g⟩ obeys
−ln a_l/ω_l → x₁ from above. For the indicator of [0.02, 0.05]:

| l | 8 | 12 | 16 | 24 | 32 | 40 |
|---|---|---|---|---|---|---|
| ω_l(m_l − x₁), m_l = −ln a_l/ω_l | 2.79 | 3.95 | 4.62 | 5.48 | 6.07 | 6.49 |
| ½ ln ω_l | 3.02 | 3.80 | 4.36 | 5.15 | 5.72 | 6.16 |

so ln a_l = −ω_l x₁ − ½ ln ω_l + O(1). A smooth bump on [0.004, 0.006] at l = 12, 16, 20, 24 has G ≤ 0.78 for t ≤ 0.9x₁ (G ≤ 1.05 required),
the spectral bracket a_l cosh ω_l t ≤ G ≤ √(a_l² cosh² + 1) holds at all 20 times, the slope of ln G after 1.5x₁ equals ω_l to 1%, and the onset time
−ln a_l/ω_l is 1.16, 1.16, 1.135, 1.106 x₁, decreasing toward x₁. It is argued, not proved, that a_l ≥ e^{−ω_l(x₁+δ)} for every δ > 0 (so growth
at every t > x₁), by a Riccati comparison w = −ψ_l′/ψ_l ≤ sup √(V̂+ω²) on x ≥ x₁, which needs V̂_l ≤ poly(l) on {x ≥ x₁} (numerically V̂ ≈ (l+1)(l+2)/r² there).
It gives non-existence for t > x₁ for smooth compactly supported data as in B3.

## 5. The truncated star (the control)

The star is x ≥ x_R with Robin u_x = κu at x_R (κ = ∞ Dirichlet), dg_proof Theorem 2.

**Theorem E** (*proved*, given Theorem 2 and the certified bound c·r_*(l) ≤ 9.1152 for every l ≥ 2, dg_proof §4.3–4.5; not independently
reviewed). Let R > 0 and κ ∈ [0, ∞]. Let U(R) = {l : R < r_*(l)} ⊂ {l : c(l) < 9.1152/R}, a finite set, and Ω_l = √(−E₀(l)) (0 if
Ĥ_{R,κ} ≥ 0), Ω(R, κ) = max_{l∈U(R)} Ω_l < ∞. Then for every l, in L²((x_R, ∞), dx), ‖C_l^{R,κ}(t)‖ ≤ cosh(Ωt) and ‖Σ_l^{R,κ}(t)‖ ≤ max(t, sinh(Ωt)/Ω), and
‖u(t)‖_s ≤ cosh(Ωt)‖u₀‖_s + max(t, sinh(Ωt)/Ω)‖v₀‖_s for every s ∈ ℝ.

*Proof.* For l ∉ U(R), R ≥ r_* so κ ≥ 0 ≥ κ_c(R) and Theorem 2 gives no negative eigenvalue, hence Ĥ_{R,κ} ≥ 0, and the multipliers are
bounded by 1 and t. For l ∈ U(R) the spectral theorem gives cosh(Ω_l t) and sinh(Ω_l t)/Ω_l, finitely many. ∎

Two explicit bounds follow from the form (min–max) argument and are *proved*. First Ω(R, κ) ≤ Ω(R, 0) for κ ≥ 0, since the form increases with κ.
Second Ω^D(R) ≤ ω_{l_D(R)} with l_D = max{l : c(l) < 6/R}, since the Dirichlet form domain sits inside the full one and Dirichlet is unstable iff R < r_s.

A bounded operator maps null sequences to null sequences, so **no sequence of data with ‖u₀‖ → 0 and ‖u(t)‖ ↛ 0 exists on the star** with κ ≥ 0
(the alternative the issue allowed). `star_control.py`, nothing shared with `dg_spectrum` (no Chebyshev, no `common.py`):

| R | Dirichlet: unstable l, Ω^D | proved bound ω_{l_D} | κ = 0 (Neumann) | κ = 10 |
|---|---|---|---|---|
| 0.05 | 2..10, 657.3 | 990 | 2..12, 2620 | 2..12, 2613 |
| 0.12 | 2..6, 117.7 | 140 | 2..8, 395.5 | 2..8, 391.5 |
| 0.30 | 2..4, 19.63 | 30 | 2..5, 43.68 | 2..5, 38.26 |
| 1.00 | 2, 1.853 | 2 | 2, 2.239 | 2, 1.994 |
| 3.00 | none | | none | none |

- SC1. The number of negative eigenvalues of the discretised problem (Sylvester inertia) equals [R < r_s] + [κ < κ_c(R)] in all 1180 cases
  (R above, κ ∈ {Dirichlet, 0, 1, 10}, l = 2..60). The unstable l are exactly {2..l_U} in each row. For Dirichlet this is l ≤ l_D(R).
  Ω agrees with the dg_spectrum table where both exist (R = 0.3 Neumann 43.68 vs 43.64, Dirichlet 19.63 vs 19.63, R = 1 Neumann 2.239 vs 2.239, Dirichlet 1.853 vs 1.853).
- SC2. V̂ ≥ 0 on [R, 50] in all 268 (R, l ≤ 60) cases where 10.7/c < R, the sufficient condition of the proof.
- SC3. r_*, r_neg by float root finding equal the certified brackets of `threshold_isolation_output.txt` to 1.2·10⁻¹³ (l = 2..100).
- SC0. The limit R → 0 (c·r = 0.2, Dirichlet and Neumann) has one negative eigenvalue equal to −ω_l² within 5.8·10⁻⁴, l = 2..60.
- SC4 (same data in both problems, R = 0.12, t₀ = 0.01, data = exact mode restricted to x ≥ x_R = 0.00346, unit norm, l = 2..24):

| l | ln a_l | ln G, point mass | spectral bracket | ln G, star κ = 0 |
|---|---|---|---|---|
| 6 | −0.00 | 0.764 | [0.763, 0.861] | 1.025 |
| 8 | −0.40 | 3.104 | [3.103, 3.104] | 1.452 |
| 12 | −6.95 | 12.380 | [12.376, 12.376] | −0.348 |
| 16 | −21.57 | 38.945 | [38.935, 38.935] | −0.265 |
| 20 | −51.19 | 94.443 | [94.420, 94.420] | −0.322 |
| 24 | −104.10 | 194.250 | [194.203, 194.203] | −0.355 |

  The star obeys ln G ≤ ln cosh(Ωt₀) = 3.262 for every l and every recorded time (max 3.231, and 0.568 vs 0.574 for Dirichlet), while the point
  mass reaches e^{194}. The values near −0.347 = −½ ln 2 are plain wave splitting, no growth.

## 6. What a finite star changes

| | point mass | star, R > 0, κ ≥ 0 |
|---|---|---|
| unstable l | all l ≥ 2 | l with R < r_*(l), i.e. c(l) ≲ 8.49/R (Dirichlet: c < 6/R) |
| top rate | ω_l ≈ l⁴/12, unbounded | Ω(R, κ) < ∞, Dirichlet ≤ ≈ 3/R², Neumann 6.6/R² at R = 0.05 |
| solution operator | unbounded on every 𝓗^s | ≤ cosh(Ωt) on every 𝓗^s |
| null data, solution not null | exists | none |
| growth dependence on data | needs data decay e^{−l⁴T/12} | none beyond the finite set |

- The cap scales like R⁻². Over R ∈ [0.05, 1] the Neumann value Ω·R² rises from 2.2 to 6.6 (*numerical evidence*, not converged). As R → 0 the constants
  diverge and the point-mass problem is recovered.
- The Phase 0 heuristic "all l ≥ 2 removed once k ≳ 8(|M|/R)²" is right in scaling. The sharp thresholds are k ≈ c²/6, with c < 6/R for Dirichlet
  (k < 6/R²) and c < 8.49/R for κ = 0 (k < 12/R²), which bracket 8.
- The caveat on κ is that an l-dependent negative interface condition loses the cap. SC5 with κ_l = −l² (R = 0.3) gives √(−E₀) = 332, 842, 1546,
  3547 at l = 20, 30, 40, 60, i.e. 0.83 to 0.985 of |κ_l|. The control assumes l-independent κ ≥ 0.
- It is argued, not proved, that a regular star supplies κ ≥ 0 at large l. A regular centre selects the r^{l+1} branch, whose logarithmic derivative in x is positive,
  so in the vacuum variable the interface acts like a Robin condition with κ of order l f(R)/R > 0 (Ψ̂ differs by the intertwiner, not worked out).
  dg_spectrum notes real stars give an energy- and l-dependent condition.

## 7. What is argued and what is not covered

- Metric variables (F, argued). H1 identifies Ψ̂ data with metric data through differential maps (Zerilli–Moncrief, DG intertwiner), which I take to have l-dependence
  polynomial for data supported at fixed r. Then B1–B3 transfer, since ε_l is super-exponentially small. An exponential-in-l loss in that map would break the transfer.
  Not checked.
- Boundary condition. The result is for the DG self-adjoint dynamics (H2). dg_spectrum finds that cut-offs at R → 0 with κx_R > −1/2 give the same modes
  (*numerical evidence*), so the rates ω_l are robust within that class.
- Hypotheses beyond this note. Theorem 1, Theorem 2 and the certified §4 results of dg_proof (no independent review yet) support C, D and E. A and B use only Theorem 1.
- Even parity only, linear only, no matter fields. The sectors l = 0, 1 are finitely many and cannot affect boundedness in l. Time is the Schwarzschild t, so "t > x₁" is in the static frame.
- Novelty is not claimed. No literature check was made for this note.

## 8. Checks

Thresholds fixed in the planted script before the real-data scripts were written (§3), each script exits non-zero on a failed check,
`check_mutations.py` runs every script unmutated and once per mutation (output `check_mutations_output.txt`).

| script | what it checks | mutations | time on rpi5 |
|---|---|---|---|
| `check_modes.py` | χ vs `common.chi_poly`; ψ₀ solves static Zerilli; V̂ two ways; closed form vs (−∂ₓ+g)Ψ_DG; ĤΨ̂ = −ω²Ψ̂ at 66 points; Π, χ/r have no positive root (l ≤ 40); Ψ̂ ~ r³; norm in two quadratures | chi, omega, jac, vhat, nodes, pow | 9 s |
| `check_fastmodes.py` | float vs mp for x(r), V̂, Ψ̂, κ_c | coef, xr | 2 s |
| `planted.py` | exact toy mode V = 3/(4x²) − 3ω/x grows as cosh ωt (error 2·10⁻⁴ωt, second order); V = 3/(4x²) bounded; Robin κ = ∓ω; sequence at fixed t₀ | omega, signV, label, kappa, dt | 2 s |
| `hadamard_sequence.py` | §3: Rayleigh quotient, growth l = 2..32, T = 20/ω_l, data norm, grid norm | omega, vhat, sign, data, dt | 170 s |
| `star_control.py` | §5: SC0–SC5 | omega0, kappa, rs, sign, bound, nodes, kl | 80 s |
| `causal.py` | §4: overlap exponent, bracket, uniform bound, slope, onset | x1, omega, near, far | 11 s |

Shared with dg_proof: the formula for χ (Lemma A) and the statement of Theorem 1. Not shared: the eigensolver, the Chebyshev discretisation, `common.py`, the sympy
construction of V̂ (V̂ here is built by Horner evaluation in mpmath and float64 and cross-checked against 2g² − V_Zerilli). The evolution tests the Cauchy problem, not
the count in Theorem 1. The count is tested by SC0/SC1 through inertia, not through Glazman.

## 9. Discrepancies and changes made to my own checks

- No result disagrees with dg_proof or dg_spectrum. The Ω values of §5 agree with the dg_spectrum growth table where both exist (to 10⁻³).
- The first SC1 prediction counted Dirichlet as [κ < κ_c] = 1. Theorem 2 counts Dirichlet as κ = +∞, so the term is 0. The check caught my bug on its first run,
  and the prediction was corrected, not the threshold.
- SC0 first used c·r = 0.02. The count was 1 throughout but E₀ was off by 0.4 at l = 24 and 2.2 at l = 60, because 10⁻¹⁶·2/h² at that cut exceeds ω². It now uses c·r = 0.2 with the
  same 10⁻² tolerance (errors ≤ 5.8·10⁻⁴).
- SC4 first ran at 20 points per 1/ω and failed the bracket by 0.19 at l = 24, which is the planted error of that resolution (7·10⁻⁴ωt). It now runs at the calibrated 40.
- SC5 first required ratio ≥ 0.9 at every l, a guess. At l = 20 it is 0.83. The check now states what is claimed (≥ 0.8, ≥ 0.98 at l = 60, increasing).
- CA1 first required |m_l/x₁ − 1| ≤ 5% from l = 8 and 1% from l = 16, a guess without the (ln ω)/(2ω x₁) rate. It failed at 33%, 9.9%, 3.8%. It now states
  0 < ω(m_l − x₁) ≤ ln ω and monotone decrease.
- CA2 first required ln G ≥ 0.9ω(t − x₁) at l = 20. The offset is −8.6 (the smooth edge of the bump, O(√ω)), so it fails. It now checks the slope (the spectral theorem)
  and the onset time.
- The first causal mutations `bound` (threshold 1.0001) and `thr` (onset shifted) did not fire, because G ≤ 0.78 and the onset test were insensitive to them.
  They were replaced by `near` and `far`, which move the data, and both fire.
- The three changes to SC5, CA1, CA2 are to my own guessed constants, made after seeing the data. They weaken nothing that A–E state.
- The phase 0 text says "looks closer to ill-posedness than an ordinary instability". It is ill-posedness, in the sense above.

## 10. Recommended tasks

1. Metric transfer (F). Reconstruct the metric perturbation of the exact mode (Zerilli–Moncrief inverse and DG's intertwiner, Regge–Wheeler gauge) for l = 2..12
   and bound its norm in an H^s norm of the slice data against ‖Ψ̂‖. Acceptance: ratio ≤ C(1+l)^p with p fitted, with a mutation that breaks it, or an exponential ratio.
   Depends on this note.
2. Overlap lower bound (D). Prove a_l ≥ e^{−ω_l(x₁+δ)} for g ≥ 0 at x ≥ x₁ via the Riccati comparison, which needs sup_{x ≥ x₁}V̂_l ≤ poly(l). Acceptance: explicit
   constants checked against the CA1 table for l ≤ 40. Depends on dg_proof Theorems 3 and L.
3. Independent review of dg_proof §4. Propositions C and Theorem E use r_neg → 0 and c·r_* ≤ 9.1152 beyond l = 3000. Acceptance: a reviewer who did not write §4
   signs the claims in §6 of dg_proof, or lists what fails.
4. Interface condition of a real star. Compute κ_eff(l) of the phantom star's interface (the Dirichlet-to-Neumann data of the interior, all l up to 40) and test
   κ_eff ≥ 0 or κ_eff ≥ κ_c. Acceptance: κ_eff(l) for l = 2..40 against κ_c(R, l), with a mutation, or a counter-example. Depends on `phase2/nonradial`.
5. Update PLAN.md Phase 2 (coordinator). The question "is the point-singularity problem ill-posed" is answered here, for Ψ̂ (B1–B4), and the star answer is §5–6.
