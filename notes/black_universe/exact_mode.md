# Why Y = R e^{K(t+x)} is exact for every black universe: an exact time-dependent family (#269)

Code [exact_mode.py](exact_mode.py), output of record [exact_mode_output.txt](exact_mode_output.txt), assembled from negative-mass-raw `I269-mac/runs`; the criteria, written before the
runs of record, are in `I269-mac/README.md`. Builds on #262 ([cminus1_mode.md](cminus1_mode.md)), #212 ([README.md](README.md)) and #232 ([linear.md](linear.md)). Status words as in AGENTS.md: the
residuals of the family are *checked symbolically*, the derivation of Sec. 2 is *proved* given the reduced equations (which S6 and N0 tie to a brute-force Ricci tensor), the numbers are *numerical
evidence*. No sentence here says a result is new (#260 reads the literature; this note was not searched against it).

## Answer

**1. The solution is the linearization of an exact family.** In the ingoing null chart with the scalar as the radial coordinate (v = t + x, φ = √2 arctan ρ, g_{vρ} = 1, g_{ρρ} = 0), member c,

```
ds² = −F dv² + 2 dv dρ + r² dΩ²,    r = λ(v) R(ρ),    F = A − 2ρ λ′/λ + 1/λ² − 1,    λ² = L + β e^{Kv},    K = −2c/π,    φ = √2 arctan ρ  (unchanged)
```

solves the five field equations identically in c, L and β (S3, *checked symbolically*; 50-digit residuals 2.8e-49, N1). With L = 1 it is `ds² = −[A − (1 + Kρ) s/(1 + s)] dv² + 2 dv dρ + (1 + s) R² dΩ²`, s = β e^{Kv}:
the throat radius at fixed φ is b(v) = λ(v), and the static member is β = 0. To first order in β: δr = (β/2) R e^{K(t+x)}, which is #262's Y_K (S5, N2).

**2. Why for every member (proved by hand, Sec. 2).** Three facts of the static problem are used and nothing else: the Ellis pair (R, φ) obeys R″ = ½φ′²R and R²φ′ = √2; the potential enters only through V and V′ on the
unchanged profile φ(ρ), which the static equations eliminate; and A − 1 solves (1 + ρ²)y″ = 2y. Then

```
 E_ρρ   r″ = ½φ′² r  is homogeneous in r          → any r = λ(v) R is allowed
 scalar (R²φ′ = √2)                                  → (F − A)′ = −2λ′/λ             → F = A − 2ρλ′/λ + h(v)
 E_θθ                                                → h = 1/λ² − 1
 E_vρ                                                → holds identically
 E_vv                                                → (λλ′)′ = K λλ′,   K = ½(R²A′ − 2ρA) + ρ
```

K is constant on the static solutions (the Wronskian of A − 1 with the other solution R²: K = ½ W[R², A − 1]), and equals −2c/π on the family. It is the coefficient of ρ + R² arctan ρ in A, and it is also 3m
(b = 1). So the rate does not come from asymptotic flatness: with c and K independent the family exists with rate K (S4b). It is a consistent scaling reduction, not a symmetry of the field equations.

**3. Fate (Sec. 5).** β > 0: no singularity at any finite v, r = √(1 + βe^{Kv}) R grows like e^{Kv/2}, and for c < 0 F tends to R²(c + K arctan ρ) < 0 at every ρ. β < 0 (collapse): λ² reaches 0 at the finite
v\* = −ln|β|/K, where r → 0 at every ρ and the Ricci scalar diverges like 1/(v\* − v) at every ρ ≠ −1/K (*proved* from the closed form; the Kretschmann invariant ∝ λ⁻⁴, N4, *numerical evidence*).

**4. What it does not give.** For β ≠ 0 the family is not asymptotically flat (F − A grows like ρ at fixed v), and its linearization is Y_ref, not the solution Y_∞ = Y_ref − Y_b/I₀ that decays at infinity
(#262). It is the fate of the Y_ref direction, which a datum with incoming exponentially growing flux excites; it is not shown to be the fate of asymptotically flat data (#250 stays open).

## 1. The chart and the reduced equations

2D metric h = −F dv² + 2 dv dρ, r(v, ρ), φ(ρ), ε = −1, R_ab = εφ_aφ_b + V g_ab, ε□φ = V′; ′ = ∂_ρ, a dot is ∂_v. Christoffels: Γ^v_vv = F′/2, Γ^ρ_vv = (FF′ − Ḟ)/2, Γ^ρ_vρ = −F′/2, the
2D curvature R₂ = −F″. With H_ab = ∇_a∇_b r: H_ρρ = r″, H_vρ = ṙ′ + F′r′/2, H_vv = r̈ − F′ṙ/2 − (FF′ − Ḟ)r′/2, (∇r)² = 2ṙr′ + Fr′², □r = 2H_vρ + FH_ρρ:

```
(ρρ)  −2r″/r = −φ′²              (vρ)  −F″/2 − 2H_vρ/r = V          (vv)  F F″/2 − 2H_vv/r + V F = 0
(θθ)  1 − r□r − (∇r)² = V r²     (φ)   −r⁻²[(r²)˙ φ′ + (r²Fφ′)′] = V′(φ)
```

N0 compares these with the brute-force Ricci residuals on a generic non-solution metric: 2.4e-51 (the residuals are 4.4); S6 repeats the symbolic comparison at arbitrary λ(v).

Static member (λ = 1, F = A, r = R, RR′ = ρ): (vρ) gives V = −A″/2 − ρA′/R²; (θθ) gives V R² = 1 − (ρA)′; together R²A″ = 2(A − 1), whose solutions are R² and Z = ρ + R² arctan ρ, so
A − 1 = cR² + KZ. With these, (R²A′ − 2ρA)′ = −2, and K := ½(R²A′ − 2ρA) + ρ is constant; W[R², Z] = 2 gives K = the coefficient of Z (S4a).

## 2. The derivation

Take r = λR, G = 1, F = A + f. Write μ = λ̇/λ.

1. (ρρ): r″ = λR″ = ½φ′²r. Holds for every λ.
2. (φ): V′ is the static value, so (2μ)R²φ′ + (R²fφ′)′ = 0. With R²φ′ = √2: f′ = −2μ, f = −2ρμ + h(v).
3. (θθ): r□r + (∇r)² = 4λλ̇ρ + λ²(ρF)′. The static (θθ) removes V: 1 − λ² − λ²(ρf)′ − 4λλ̇ρ = 0. With f from 2: (ρf)′ = −4ρμ + h, so h = 1/λ² − 1.
4. (vρ): −F″/2 − 2μR′/R − F′R′/R with F′ = A′ − 2μ, F″ = A″, equals −A″/2 − A′R′/R = V.
5. (vv): H_vv/r = λ̈/λ − F′μ/2 − (FF′ − Ḟ)R′/(2R) and V = −A″/2 − A′R′/R, so (vv) reads
   −2λ̈/λ + F′μ + (ρ/R²)[F(F′ − A′) − Ḟ] = 0. Insert F = A − 2ρμ + h, Ḟ = −2ρμ̇ + ḣ, ḣ = −2μ/λ² and multiply by R². The terms in A combine to 2μ(K − ρ). The ρ and ρ² terms cancel
   (the ρ terms by h = 1/λ² − 1), and what is left is −2μ̇ − 4μ² + 2Kμ = 0, that is (λλ̇)˙ = Kλλ̇.

The ρ dependence of the leftover is the overall factor 1/(1 + ρ²) (S2). Integrating, (λ²)˙ = Kβe^{Kv}, λ² = L + βe^{Kv}.

## 3. The family

- **c = 0** (Ellis, K = 0, A = 1): λ² = L + βv, F = (1 − βρ)/λ²; the linearization is Y = (β/2) R (t + ρ), the secular mode of the massless case.
- **β = 0, L ≠ 1** (N6): r = √L R, F = A + 1/L − 1 solves the equations with the same potential, for every member. This is a static variation of the throat radius at fixed φ, for every c, not only c = 0. #212's
  operator does not see it: at W = 0 the tρ equation it was derived from is empty, and the operator on Y = R gives 1.4142 (c = −1, ρ = 1) and 0.3031 (c = −3/10).
- **c > 0** (K < 0): the same family, with s decaying toward v → +∞; not examined further.
- **Horizon at c = −1** (ρ_h = 0): the throat ρ = 0 has r = λ(v) = b(v), r′ = 0 and ṙ = λ̇ ≠ 0 there for all v, and the horizon value Y|_h = (β/2)e^{κv} is b(v) − 1 to first order: the growth at κ of #232 is the growth of the throat radius.

## 4. The linearization in the static gauge

First order in β (L = 1): δr = (β/2)R e^{Kv}, δF = −β(1 + Kρ)e^{Kv}, δG = 0, δφ = 0. The gauge δg_{tρ} = 0 is reached by v = t + x + εξ, ξ = β e^{Kt} Ξ(ρ), Ξ′ = (1 + Kρ)e^{Kx}/A²
(ξ_ρ = −δF/A²), and then

```
Y = (β/2) R e^{K(t+x)},     H₂ = β(1 + Kρ) e^{K(t+x)}/A,     H₀ = −β(1 + Kρ)e^{K(t+x)}/A + 2Kβ e^{Kt} Ξ(ρ)        (g_tt = −A(1 + H₀), g_ρρ = (1 + H₂)/A)
```

Ξ is a quadrature; its additive constant is the time reparametrization t → t + c₀e^{Kt}, which changes H₀ alone (#262's review: at W = K the equations do not involve H₀'s amplitude). N2 evaluates the
brute-force linearized E_tt, E_tρ, E_ρρ, E_θθ and the scalar equation at this perturbation (β = 2) at 50 digits: 1.2e-50 relative to the sum of the terms at c = −1, −3/10, −3/2, −10, 8 points each, Ξ = 0 and 3/7
(bar 1e-25). E_tρ is linear in H₂ and fixes it, so this H₂ is the b of the coordinator's `ind262c.py`. S5 is the symbolic statement for Y: #212's operator annihilates R e^{Wx} at W = the rate.

## 5. Fate

F − A = [1 − L − s(1 + Kρ)]/(L + s). Take L = 1.

| | β > 0 | β < 0 |
|---|---|---|
| range of v | all of ℝ | v < v\* = −ln\|β\|/K |
| areal radius | r² = (1 + s)R² → ∞, d ln r/dv → K/2 (N5: 0.31830989 = K/2 at c = −1) | r → 0 at every ρ as λ² = 1 + s → 0, λ² ≈ K(v\* − v) |
| F | F → A − 1 − Kρ = R²(c + K arctan ρ) = (2c/π)R² arccot ρ, negative at every ρ for c < 0 (−0.88104 at c = −1, ρ = ½) | F ≈ (1 + Kρ)/λ² → ±∞ |
| curvature | finite at every v; Kretschmann → a finite limit (2.58708 at c = −1, ρ = ½; N5) | Ricci scalar ≈ −2(1 + Kρ)/(R⁴λ²), diverging for ρ ≠ −1/K; Kλ⁴ → 8.5423 at c = −1, ρ = ½ (N4) |
| verdict | expansion, no singularity at finite v (late slices have F < 0: ∂_v spacelike everywhere for c < 0) | collapse to a curvature singularity on the null hypersurface v = v\* |

The Ricci scalar is the trace of the field equations, −Fφ′² + 4V with φ′² = 2/R⁴ and V bounded, and F = A − 1 + (1 − Kρs)/λ² from λλ̇ = Ks/2: *proved*. The line ρ = −1/K is exceptional: F = A there for all v, and the
Kretschmann invariant stays finite as λ² → 0 (69.979 at c = −1, 123.96 at c = −3/2; info lines of N4), where r → 0 as well; it was not examined further. The sign of β is the sign of the horizon coefficient of Y_ref
(c₂ ∝ β at c = −1), so the two signs are the two outcomes along the Y_ref direction.

## 6. Checks (derive: S0–S6; check: N0–N6; mutations: 11)

| id | what | result |
|---|---|---|
| S0 | static member solves the five equations in the null chart (planted truth) | 0 |
| S1 | arbitrary λ(v): E_ρρ, E_vρ, E_θθ, scalar vanish | 0 |
| S2 | E_vv = 2(Kλλ̇ − λλ̈ − λ̇²)/((1 + ρ²)λ²) | equal |
| S3 | λ² = L + βe^{Kv}: all five residuals | 0 |
| S4a, S4b | K = ½(R²A′ − 2ρA) + ρ = −2c/π; c, Ks independent: rate Ks | equal; equal |
| S5 | #212's operator annihilates R e^{Wx} at W = rate | 0 |
| S6 | hand-reduced equations at arbitrary λ(v) = brute force | equal |
| N0 | hand-reduced against brute-force residuals, generic metric (residual 4.4) | 2.4e-51 (bar 1e-30) |
| N1 | exact family, 5 members, 120 points, both signs of β | 2.8e-49 (bar 1e-30) |
| N2 | static-gauge brute-force linearization, 4 members, 64 point-gauge pairs | 1.2e-50 (bar 1e-25) |
| N3 | Kretschmann routine on Schwarzschild (planted) | exact to 40 digits |
| N4 | β < 0: Kλ⁴ converges to a positive limit, K > 1e12 at λ² = 1e-8, 9 cases | pass |
| N5 | β > 0: finite limit (< 1e-5), F < 0, d ln r/dv = K/2 (1e-10), 12 cases | pass |
| N6 | static family L ≠ 1, 5 members | 3.3e-49 (bar 1e-30) |

All 15 pass in the runs of record (derive 8 of 8, check 7 of 7; 7 s and 6 s). Eleven mutations, each with the checks it must fire declared before the runs of record: in the runs of record the fired sets equal the
declared ones (11 of 11; `exact_mode_output.txt`). The perturbed candidates of the issue are `kf` (2 → 2.02), `kh` (1 → 1.01), `rate` (K → 1.001K) and `shape` (R → ρ): each fires S3, N1 and N2 or more.
`naive` (b → b(v) with F = A unchanged) fires ten checks, and `eps` does not fire N2 because the first-order equations at δφ = 0 do not contain ε. The double-null solver of #232 was not used (the issue's alternative).

## 7. Proposed text (the files are not touched here)

`cminus1_mode.md` Sec. 1, replacing "Why it exists for every c is not understood here (recommended task 3)": *Update from #269 ([exact_mode.md](exact_mode.md)): Y_K = R e^{K(t+x)} is the linearization of the exact
family r = λ(v)R, F = A − 2ρλ′/λ + 1/λ² − 1, λ² = L + βe^{Kv} in the ingoing null chart (checked symbolically); K = ½(R²A′ − 2ρA) + ρ is a constant of the static equations.* In Sec. 6, the first bullet
loses its first clause; task 3 is done.

`README.md` Sec. 5, after the Update from #262: *Update from #269: the exact solution Y_ref is the linearization about the member of an exact nonlinear solution, a time-dependent throat radius b(v)² = L + βe^{Kv}
with the scalar unchanged and F = A − 2ρλ′/λ + 1/λ² − 1 (checked symbolically for every member; [exact_mode.md](exact_mode.md)). For β < 0 it collapses at finite advanced time with a diverging Ricci scalar, for
β > 0 it expands; it is not asymptotically flat and linearizes to Y_ref, not to Y_∞.*

Paper IV (`paper/4-black-universes/main.tex`), replacing "Why such a solution exists for every member is open." (paper signature):

```latex
It is the linearization about the member of an exact solution of the field equations (\#269; checked symbolically). In the ingoing null chart
$\dd s^2=F\,\dd v^2-2\,\dd v\,\dd\rho-r^2\dd\Omega^2$, $v=t+x$, with the scalar unchanged, take $r=\lambda(v)R(\rho)$ and $F=A-2\rho\dot\lambda/\lambda+\lambda^{-2}-1$.
Four of the five field equations hold for every $\lambda(v)$; the fifth is $\dd(\lambda\dot\lambda)/\dd v=K\lambda\dot\lambda$, so $\lambda^2=L+\beta e^{Kv}$, and to first
order in $\beta$ the perturbation is $\delta r=\tfrac\beta2R\,e^{K(t+x)}$. The rate is $K=\tfrac12(R^2A'-2\rho A)+\rho$, a constant of the static equations. For $\beta<0$ the areal radius
vanishes at the finite $v_*=-K^{-1}\ln|\beta|$ with a diverging Ricci scalar; for $\beta>0$ it grows as $e^{Kv/2}$. The family is not asymptotically flat for $\beta\ne0$ and linearizes to the
solution that grows toward spatial infinity, not to the one that decays there, so it does not give the fate of asymptotically flat data.
```

## 8. What this does not show

- The fate of asymptotically flat data, or of the physical mode Y_∞ = Y_ref − Y_b/I₀: Y_b is not in the family, and nothing is said about its nonlinear completion (#250).
- The global causal structure: where the event horizon lies for β ≠ 0, what the line ρ = −1/K is at β < 0 (finite Kretschmann as r → 0), and whether v = v\* ends the spacetime or can be continued.
- Whether other static Einstein–scalar solutions admit r = λ(v)R. The derivation used R²φ′ = const and R″ = ½φ′²R. No other (R, φ) pair was tried.
- That the double-null solver of #232 reproduces the first-order mode: not run.
- The values of the constants are those of the runs of record; the Kretschmann limits are 50-digit evaluations at the stated points, not a proof.

## Recommended tasks

1. **The family as planted truth for the nonlinear evolution code.** *What:* feed #217's kernel (`evolve.py`) data on two characteristic rays from the exact family (β = ±0.05, c = −1.5 and −1), in its double-null chart,
   and compare r(t, ρ) with λ(v)R. *Why:* an exact nonlinear time-dependent solution with a known singularity time v\* for β < 0 tests the kernel's growth, convergence and horizon crossing beyond linear order.
   *Acceptance:* error ratio 4.0 ± 0.3 per halving of h; v\* recovered to 1e-3 for β < 0; a sign mutation of β fails. *Depends on:* #217, #269.
2. **Causal structure of the family.** *What:* the event and apparent horizons, the null hypersurface v = v\*, the line ρ = −1/K and the ρ → ±∞ ends, for β of both signs, with a conformal diagram. *Acceptance:*
   the apparent horizon where (∇r)² = 0 (closed form in ρ and s) against a numerical root, to 1e-12; the diagram's claims each tested by a null-ray integration. *Depends on:* #269.
3. **Second order about Y_ref and Y_∞.** *What:* at c = −1 the second-order response to Y_∞ (with Y_b): does it have a closed form, and does it show collapse or expansion? *Acceptance:* the second-order
   metric perturbation in the static gauge solves the linearized equations with the quadratic source to 1e-25 at 40 digits; the sign of its throat term decides which of β ≷ 0 it resembles. *Depends on:* #249, #269.

## Reproduce

From this directory, `~/Source/negative-mass/.venv/bin/python`, one thread, about 15 s in all:

```
exact_mode.py derive [--mutate=kf|kh|rate|shape|naive|pot|eps|charge|hand_drop]   # S0 - S6, 8 s
exact_mode.py check [--only=N4,N5] [--mutate=kf|kh|rate|kappa|shape|naive|pot|eps|hand_drop|kret_drop]   # N0 - N6, 6 s
exact_mode.py mutations DIR PREFIX      # fired sets of stored logs DIR/PREFIX<name>.log against the declared ones
```

`exact_mode_output.txt` is `derive`, `check` and `mutations` of the runs of record (`I269-mac/runs`), in that order.
