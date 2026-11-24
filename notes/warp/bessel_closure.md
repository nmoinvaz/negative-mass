# Closure fails on every open set for the Bessel shift family

November 2026, issue [#162](https://github.com/nmoinvaz/negative-mass/issues/162). Script
[`bessel_closure.py`](bessel_closure.py), its stdout [`bessel_closure_output.txt`](bessel_closure_output.txt). It builds
on [`det_classify.md`](det_classify.md) §3 (#80) and [`two_scalars.md`](two_scalars.md) (#65, reviewed in #79). It
imports only `closure` and `dform` from [`two_scalars.py`](two_scalars.py), read-only; the Ricci tensor, the coframe,
the ODE reduction and the Frobenius coefficient are its own code.

Status words as in AGENTS.md: **proved** (by hand), **checked symbolically**, **numerical evidence**, **not decided**.

## Answer

Closure of [`two_scalars.md`](two_scalars.md) fails on every nonempty open set, for every member u = c x + g(ρ) of the
Bessel family with c ≠ 0, for ε = ±1 and for every α, β, d (**proved**; every identity used is **checked
symbolically**). The obstruction is the Frobenius integrability of the two null covectors k = a + b and m = a − b of R.
In the coframe θ⁰ = dt, θ¹ = dx + u dt, θ² = dρ, with W = √(g′² + 2c²), the coefficients of k̃ ∧ dk̃ and m̃ ∧ dm̃ (k̃ =
√2 k, m̃ = √2 m) on dt ∧ dx ∧ dρ are, on the ODE g″ + g′/ρ = ε (g′² + 2c²),

    O_k = g′ (g′² + 4c² + 5cεW − 2c g′/(ρW)),     O_m = g′ (g′² + 4c² − 5cεW + 2c g′/(ρW)),
    O_k + O_m = 2 g′ (g′² + 4c²).

The sum needs no ODE and has the only real zero g′ = 0 when c ≠ 0. There g″ = 2εc² ≠ 0, so the zeros of g′ are isolated,
and k or m fails Frobenius at every point of a member except on a discrete set of ρ. Closure on an open set would need
O_k = O_m = 0 on it.

**This is not "nonzero for all (c ≠ 0, g′, ρ, ε)" as the issue worded it.** The pair vanishes at g′ = 0, and each of
O_k, O_m vanishes alone on a hypersurface of the jets. What is true, and what closure needs, is that the zero set
contains no open set of any member (§4).

| claim | status |
|---|---|
| R = a a − b b on the ODE, ε = ±1 (own Ricci tensor), a ∧ b ≠ 0 | **checked symbolically** |
| closure on an open set ⟹ k ∧ dk = 0 and m ∧ dm = 0 | **proved** (one line, §1) |
| k̃ ∧ dk̃ = (βα′ − αβ′ + g′β² − cβγ) dt ∧ dx ∧ dρ for k̃ = αθ⁰ + βθ¹ + γθ² | **proved**; **checked symbolically** |
| O_k, O_m on the ODE; O_k + O_m = 2g′(g′² + 4c²) | **checked symbolically** (and **proved**, §3) |
| closure fails on every open set, c ≠ 0 (Theorem 1) | **proved** |
| O_k alone, O_m alone vanish on no interval (Theorem 2) | **proved**; the elimination polynomial **checked symbolically** |
| the Bessel evaluation agrees with the ODE reduction | **numerical evidence**, 80 digits |
| c = 0: O = g′³, the Frobenius argument of #65 | **checked symbolically** |

## 1. The test, and why failing it is enough

[`two_scalars.md`](two_scalars.md) stage 2: closure holds when some χ makes e^χ k and e^−χ m closed, with k = a + b, m =
a − b. By #65's lemma (**proved**; reviewed in #79) the only λ with R − λg = a a − b b is λ = 0, for every axisymmetric
shift. So R = a a − b b = ½(k ⊗ m + m ⊗ k) and R(ξ, ξ) = k(ξ) m(ξ). Unique factorization of a quadratic form gives: k
and m are fixed by R up to k → λk, m → m/λ and a swap. So the two null lines of R are an invariant of R, and so is every
condition below.

**The test used here is a necessary condition for closure, not an equivalent one.** If e^χ k is closed then dk = −dχ ∧
k, so k ∧ dk = 0. If e^−χ m is closed then dm = dχ ∧ m, so m ∧ dm = 0. These are the two 2-form equations ω ∧ k = −dk
and ω ∧ m = dm of #65 wedged with k and with m, with ω = dχ. So failure of k ∧ dk = 0 or m ∧ dm = 0 at one point of an
open set gives failure of closure on that set. The converse fails: k and m can each be integrable while the factors e^±χ
cannot be matched (planted P4 of §6: k = dt, m = (1 + tρ) dρ, and the planar u = c x + d of §5). The test has a second
layer, a compatibility of the two integrating factors, that is not used here.

## 2. The family and its null pair

ds² = −dt² + (dx + u dt)² + dρ² + ρ² dφ², u = c x + g(ρ). D = det R₃ = g′² (2c² + g′² − S)(2c² + g′² + S)/8 with S = g″
+ g′/ρ (#80). On the ODE S = ε (g′² + 2c²) the solutions are g = −ε ln(α J₀(√2 c ρ) + β Y₀(√2 c ρ)) + d, defined where
the argument of the logarithm does not vanish. The additive constant d is an x-translation by d/c, so it does not enter.

R = a a − b b on the ODE for ε = ±1, with a = (g′ dρ + 2c θ¹)/√2, b = W ν/√2, ν = θ⁰ + ε θ¹ (`frame`, exact, with the
script's own Ricci tensor; #80 checked the same identity with `phantom_source.ricci`). The minors of (√2 a, ν) in the
coframe are (−2c, −g′, −εg′), with squares summing to 2 (g′² + 2c²) > 0, so a ∧ b ≠ 0 and k ∧ m ≠ 0. The structure
equations are

    dθ⁰ = 0,   dθ¹ = du ∧ dt = c θ¹ ∧ θ⁰ + g′ θ² ∧ θ⁰,   dθ² = 0.

(du = c dx + g′ dρ = c θ¹ − c u θ⁰ + g′ θ², and θ⁰ ∧ θ⁰ = 0.)

## 3. The obstruction

**Lemma 1 (proved).** For k̃ = α θ⁰ + β θ¹ + γ θ² with α, β, γ functions of ρ, and ′ = d/dρ,

    k̃ ∧ dk̃ = (β α′ − α β′ + g′ β² − c β γ) θ⁰ ∧ θ¹ ∧ θ²,   θ⁰ ∧ θ¹ ∧ θ² = dt ∧ dx ∧ dρ.

*Proof.* dk̃ = (α′ + g′β) θ² ∧ θ⁰ + β′ θ² ∧ θ¹ + c β θ¹ ∧ θ⁰ by the structure equations (dα = α′ θ², and θ² ∧ θ² = 0).
Wedge with k̃ and keep the three terms with all of θ⁰, θ¹, θ²: α θ⁰ ∧ β′ θ² ∧ θ¹ = −αβ′ θ⁰θ¹θ², β θ¹ ∧ (α′ + g′β) θ² ∧
θ⁰ = β(α′ + g′β) θ⁰θ¹θ², and γ θ² ∧ cβ θ¹ ∧ θ⁰ = −cβγ θ⁰θ¹θ². The components with dφ vanish since nothing depends on φ
and k̃ has no dφ part. ∎

For k̃ = √2 (a + b): α = W, β = 2c + εW, γ = g′. For m̃ = √2 (a − b): α = −W, β = 2c − εW, γ = g′. Then

    O_k = 2c W′ + g′ (2c + εW)(c + εW),     O_m = −2c W′ + g′ (2c − εW)(c − εW),

since βα′ − αβ′ = (2c + εW) W′ − W εW′ = 2c W′, and g′β² − cβγ = g′β (β − c). With W W′ = g′ g″ and g″ = εW² − g′/ρ (the
ODE), 2c W′ = 2c g′ (εW² − g′/ρ)/W, and the two closed forms of the Answer follow (`frobenius` (b), exact for ε = ±1).

**Lemma 2 (proved).** O_k + O_m = 2 g′ (g′² + 4c²), with no use of the ODE.

*Proof.* The W′ terms cancel, and (2c + εW)(c + εW) + (2c − εW)(c − εW) = 2 (2c² + W²) = 2 (g′² + 4c²) using ε² = 1 and
W² = g′² + 2c². ∎ (`sum` (a), exact for ε = ±1 for a symbolic g, with no ODE.)

## 4. Where it vanishes, and why no open set

**Theorem 1 (proved).** Let c ≠ 0 and let g ∈ C² on an open interval J ⊂ (0, ∞) solve g″ + g′/ρ = ε (g′² + 2c²), ε = ±1.
Then closure fails on every nonempty open set U ⊂ {ρ ∈ J}, for u = c x + g(ρ).

*Proof.* Suppose closure holds on U. By §1, O_k = O_m = 0 on U. Their values depend on ρ only, and the ρ-values of U
contain an interval I. By Lemma 2, g′ (g′² + 4c²) = 0 on I. Since c ≠ 0, g′² + 4c² > 0 (its discriminant in g′ is
−16c²), so g′ = 0 on I, so g″ = 0 on I. At a point of I the ODE gives g″ = ε (0 + 2c²) − 0 = 2εc² ≠ 0. Contradiction. ∎

(`isolated` (a): g″ = 2εc² at g′ = 0, from the script's own ODE rule. `isolated` (b), planted: g = const solves the ODE
at c = 0 and at no c ≠ 0, which is where the step fails for #65's b(ρ).)

Every choice of ε, α, β, d is covered: g = −ε ln|ψ| with ψ = αJ₀ + βY₀ solves the ODE where ψ ≠ 0, and Theorem 1 uses
only the ODE.

**What the zero set is.**

| expression | vanishes where | on a member |
|---|---|---|
| O_k + O_m | g′ = 0 | at the zeros of g′, which are simple (g″ = 2εc² ≠ 0): a discrete set of ρ |
| (O_k, O_m) both | g′ = 0 (Lemma 2, c ≠ 0) | the same |
| O_k alone | g′ = 0, or B_k := g′² + 4c² + 5cεW − 2c g′/(ρW) = 0 | discrete, Theorem 2 |
| O_m alone | g′ = 0, or B_m := g′² + 4c² − 5cεW + 2c g′/(ρW) = 0 | discrete, Theorem 2 |

B_k + B_m = 2 (g′² + 4c²) > 0, so B_k and B_m never vanish together. Each of B_k = 0 and B_m = 0 is a hypersurface of
the jet space (c, g′, ρ) with a point on it: the exact jets (c, g′, ρ, ε) = (−2, 1, 4/39, +1) have O_k = 0 and O_m = 34,
and (2, 1, 4/39, +1) have O_k = 34 and O_m = 0 (for k̃, m̃; the table of §6). At both, `two_scalars.closure` still fails
('no chi'), through the other covector. So one covector alone can be integrable at a jet; the pair never is, except
where g′ = 0.

**Theorem 2 (proved; the elimination polynomial checked symbolically).** For c ≠ 0, neither O_k nor O_m vanishes
identically on any interval of a member.

*Proof.* Suppose O_k ≡ 0 on an interval I. Zeros of g′ are isolated, so shrink I to have g′ ≠ 0. With G = g′ and r =
1/ρ, H_k = (W/G) O_k = (G² + 4c²) W + 5cεW² − 2c G r ≡ 0 on I, and along the member (G, r)′ = (εW² − G r, −r²) and W′ =
G G′/W. The coefficient −2cG of r in H_k is nonzero, so r = F(G) on I with F rational in (G, W), and L_V H_k ≡ 0 on I as
well. Substituting r = F(G) turns L_V H_k into a function of G alone whose numerator is p₀(G) + p₁(G) W with polynomial
p₀, p₁. So N := p₀² − p₁² (G² + 2c²) = 0 on I. The script computes (`each`, exact, for k and m and ε = ±1, all four
giving the same N)

    N = −(G² + 2c²) P(G, c),
    P = G¹² − 21c² G¹⁰ − 78c⁴ G⁸ + 720c⁶ G⁶ + 5936c⁸ G⁴ + 13056c¹⁰ G² + 8704c¹²,

monic of degree 12 in G. So G takes values in a finite set on I, hence is constant by continuity. Then G′ = 0 and r =
F(G) is constant, so ρ is constant on I. Contradiction. The same for O_m. ∎

Observation, not a criterion: the real roots of P at c = 1 are G = ±2.9388 and ±4.7402. At such a jet the elimination
polynomial vanishes, so a member can meet {O = 0} tangentially there (a double zero); whether any root is a genuine
tangency of the flow was not tested in the script.

## 5. The two neighbours

**c = 0** is #65's b(ρ). The same code gives O_k = O_m = g′³ for both signs of g′ and of ε, and for ν± = dρ ± ((1 + εb)
dt + εdx) the coefficient of ν± ∧ dν± is g′ = b′ (`c0`, exact; #65: ν± ∧ dν± = b′ dx ∧ dρ ∧ dt). The sum identity gives
2g′³. So the argument for c ≠ 0 contains #65's Frobenius argument, and the one step that differs is Theorem 1's last
line: at c = 0 the ODE has the solution g = const (flat), which is why #65 excludes b′ = 0.

**g′ ≡ 0**, u = c x + d with c ≠ 0, is outside the ODE (D ≡ 0 as a planar shift, #80 §2). R = c² (θ¹θ¹ − θ⁰θ⁰) exactly
(`planar`, own Ricci), so a = cθ¹, b = cθ⁰. Both null covectors are integrable (O = 0), and `two_scalars.closure` fails
at the second stage: the candidate ω = (−c (c x + d), −c, 0, 0) has (dω)_tx = c² ≠ 0. So closure fails here too, for the
same reason as #80's Θ′ = −2c² for the planar shift u = c x + d. **Corollary (proved, using #80 §3 and the above).** For
c ≠ 0, let u = c x + g(ρ) be C² with D = 0 on an open set U. On the open part of U where g′ ≠ 0 the ODE holds with
locally constant ε (S = ±(g′² + 2c²) are disjoint, and g″ is continuous), and Theorem 1 applies there; where g′ ≡ 0 on
an open subset the shift is planar and fails as above. So closure fails on every open set where u = c x + g(ρ) has D = 0
and c ≠ 0 (an open set whose ρ-range has g′ ≢ 0 contains an open subset with g′ ≠ 0; closure on a set gives closure on
its open subsets).

## 6. Planted truth, and what the points show

`planted`, through the same functions (`O_of`, `two_scalars.closure`):

| pair (a, b) | O_k, O_m | `closure` | expected |
|---|---|---|---|
| (dρ, dt) | 0, 0 | ok, ω = 0 | pass |
| k = e^−χ d(t + x), m = e^χ d(t − x + ρ²/2), χ = c x + g(ρ) | 0, 0 | ok, ω = dχ = (0, c, g′, 0) | pass |
| (dt, x dρ) | 1, −1 | fails, 'no chi' | fail (Frobenius) |
| k = dt, m = (1 + tρ) dρ | 0, 0 | fails, 'omega = d chi needed' | fail (compatibility only) |

The last row is the reason §1 says "necessary": the Frobenius test is blind to it, and `closure` catches it.

`points`: exact jets, ε = ±1 each (20 points, x₀ = 1/3, t₀ = 1/10). #80's six jets (c, g′) = (2, 1), (6, 17), (0, 1/2)
at ρ = 5/7 are reproduced: `closure` fails at all six, and the coordinate O_k, O_m equal the closed forms exactly at all
20 points. (For c = 0, 1/8 = (1/2)³.)

| (c, g′, ρ) | ε | O_k | O_m | `two_scalars.closure` |
|---|---|---|---|---|
| (2, 1, 5/7) | +1 | 677/15 | −167/15 | fails |
| (2, 1, 5/7) | −1 | −223/15 | 733/15 | fails |
| (6, 17, 5/7) | +1 | 1595569/95 | −196979/95 | fails |
| (6, 17, 5/7) | −1 | −245531/95 | 1644121/95 | fails |
| (0, 1/2, 5/7) | ±1 | 1/8 | 1/8 | fails |
| (−2, 1, 4/39) | +1 | 0 | 34 | fails |
| (−2, 1, 4/39) | −1 | 60 | −26 | fails |
| (2, 1, 4/39) | +1 | 34 | 0 | fails |
| (2, 1, 4/39) | −1 | −26 | 60 | fails |
| (1, 1, 5/7) | +1 | 5 + 61√3/15 | 5 − 61√3/15 | fails |
| (1, 1, 5/7) | −1 | 5 − 89√3/15 | 5 + 89√3/15 | fails |
| (3, 7/2, 2/3) | +1 | 38507/88 | −8785/88 | fails |
| (3, 7/2, 2/3) | −1 | −12313/88 | 42035/88 | fails |
| (−6, 7, 3/11) | +1 | −763 | 3465 | fails |
| (−6, 7, 3/11) | −1 | 3857 | −1155 | fails |
| (1, 0, 5/7), (−3, 0, 1/2) | ±1 | 0 | 0 | pointwise: 'omega closed' |

At the four jets with g′ = 0 both covectors are integrable and `closure` reports 'omega closed'. That is pointwise and
not decisive: `closure` evaluates ω at the point, so dω = 0 is not tested there (as #80 noted for its own use). It does
not matter for open sets, because the zeros of g′ are isolated.

`bessel`: five actual members (c, ε) = (½, +1), (3/2, −1), (−2, +1), (1, −1), (5/2, +1) with random rational α, β and no
ODE substituted anywhere. At 80 digits and 60 points the coordinate O_k, O_m of the explicit expression equal the closed
forms to 1.77e-74 (threshold 1e-50), and O_k + O_m equals 2g′(g′² + 4c²) to 4.03e-76 (**numerical evidence**; it tests
that the reduction on the ODE and the Bessel solutions agree, not the proof). Observation, not a criterion: on the grid
0.3 < ρ < 8.95 the numbers of sign changes of (g′, O_k, O_m) are (2, 2, 5), (9, 21, 9), (12, 28, 12), (3, 11, 3), (16,
16, 36) for the five members. In each, one of O_k, O_m has more sign changes than g′, so members do cross the
hypersurface B_k = 0 or B_m = 0; this is why the sum, not either alone, is the test at a point.

## 7. Checks and controls

10 checks, all pass; 13 controls, each exits 1 with its target check failing and no other (`bessel_closure_output.txt`;
each control's own stdout is in the raw folder, `final_mut/`). Per check, in seconds: `frame` 1, `frobenius` 1, `sum` 0,
`isolated` 0, `each` 0, `points` 38, `bessel` 59, `planted` 0, `c0` 0, `planar` 1; 101 s in all.

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `frame` | R = a a − b b on the ODE, a ∧ b ≠ 0 | `frame-nu` | ν = θ⁰ − εθ¹ |
| `frobenius` | coordinate k̃ ∧ dk̃ = the coframe formula; closed forms on the ODE | `frob-struct`, `frob-ode` | drop the c β γ term; ODE with −ε |
| `sum` | O_k + O_m = 2g′(g′² + 4c²); discriminant; real roots | `sum-factor` | expects 2g′(g′² + 2c²) |
| `isolated` | g″ = 2εc² at g′ = 0; planted: g = const at c = 0 | `isolated-c0` | the step sampled at c = 0 |
| `each` | N = −(G² + 2c²) P, P monic; H = (W/G) O | `each-w2`, `each-friction` | W² = G² + c²; drop −G r in G′ |
| `points` | 20 exact points, closure fails where (O_k, O_m) ≠ 0 | `points-closed` | (a, b) = the closed (dρ, dt) |
| `bessel` | 80-digit Bessel members against the closed forms | `bessel-freq` | J₀(cρ) for J₀(√2 cρ) |
| `planted` | closed pairs pass, (dt, x dρ) fails, compatibility-only fails | `plant-pass`, `plant-fail` | non-integrable m; closed (dt, ρ dρ) |
| `c0` | O = g′³, ν± ∧ dν± = g′ | `c0-sign` | expects −g′³ |
| `planar` | R = c²(θ¹θ¹ − θ⁰θ⁰), O = 0, (dω)_tx = c² | `planar-closed` | the closed (c dx, c dt) |

`points-closed` and `planar-closed` are the controls that make the obstruction vanish. All 13 also fired on their
targets in a scratch run of the draft script before the README of the raw folder was written
(`scratch_draft_ctl_*.log`); that run is not of record.

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/bessel_closure.py --no-controls      # the 10 checks, about 100 s on a Raspberry Pi 5
    python3 notes/warp/bessel_closure.py --mutate=NAME       # one control
    python3 notes/warp/bessel_closure.py                     # the checks, then every control (about 23 min)

Times are on a Raspberry Pi 5 (aarch64, 4 cores shared with other jobs), one process at nice 10; the full run took 22
min 43 s wall with a peak resident memory of 72 MB (`/usr/bin/time -v`).

## 8. What the code proves, and what is assumed

The code proves, as exact identities (sympy): R = a a − b b on the ODE for ε = ±1 (`frame`); the coordinate coefficient
of k̃ ∧ dk̃ equals the structure-equation formula, and equals the closed forms on the ODE (`frobenius`); O_k + O_m =
2g′(g′² + 4c²) (`sum`); the elimination polynomial N (`each`); and, at 20 exact jets, that `two_scalars.closure` fails
exactly where (O_k, O_m) ≠ (0, 0) (`points`).

By hand: that a∧b ≠ 0 from the minors; the uniqueness of {k, m} up to rescaling and a swap (factorization of R(ξ, ξ));
closure ⟹ Frobenius (§1); Lemma 1 from the structure equations; Theorem 1; Theorem 2's reduction (the step from N = 0 to
ρ constant); the corollary of §5.

Imported and trusted: `two_scalars.closure` and `two_scalars.dform` (#65, reviewed in #79), used in `points`, `planted`
and `planar`. The Frobenius coefficient route and the closure route are different computations of the same failure, and
agree at the 20 points; they share the pair (a, b), which `frame` checks against an own Ricci tensor, and the exterior
derivative of the components. Not independent of #80: the family, its ODE and the identity R = a a − b b are #80's,
re-derived here; the coordinator's `ind80.py` K4 confirms the factorization of D and the reduction.

Assumed: two real minimally coupled scalars with a Lorentzian kinetic form and any potential (so λ = 0 by #65's lemma);
G_ab = 8π T_ab; g ∈ C² on the interval; sympy 1.14.0, mpmath 1.3.0.

Not tested: closure on the members at the isolated points g′ = 0 themselves (pointwise unreachable); non-axisymmetric or
time-dependent deformations of the family; whether a root of P is a genuine tangency (§4).

## 9. Lemmas for Lean (#135)

Statements only; ℝ throughout; do not write Lean here.

1. **Structure** (assumed for the differential geometry): in a 3-dimensional coframe with dθ⁰ = 0, dθ¹ = c θ¹ ∧ θ⁰ + g′
   θ² ∧ θ⁰, dθ² = 0 and dα = α′ θ² (same for β), k̃ ∧ dk̃ = (βα′ − αβ′ + g′β² − cβγ) θ⁰θ¹θ².
2. **Closed forms** (`ring`/`field_simp`). With W > 0, W² = g′² + 2c², W′ = g′ g″/W and g″ = εW² − g′/ρ: 2cW′ + g′(2c +
   εW)(c + εW) = g′ (g′² + 4c² + 5cεW − 2c g′/(ρW)), and the m-version with the signs of W flipped.
3. **Sum** (`ring`, using ε² = 1). (2c + εW)(c + εW) + (2c − εW)(c − εW) = 2 (g′² + 4c²) when W² = g′² + 2c².
4. **Positivity** (`nlinarith`). c ≠ 0 ⟹ g′² + 4c² > 0.
5. **Isolation**. If g ∈ C² solves g″ + g′/ρ = ε (g′² + 2c²), c ≠ 0, ε² = 1, then g′(ρ₀) = 0 ⟹ g″(ρ₀) = 2εc² ≠ 0, so g′
   ≡ 0 on an interval is impossible.
6. **Closure ⟹ Frobenius**. dk = −ω ∧ k ⟹ k ∧ dk = 0, since k ∧ k = 0. Same for m.
7. **Uniqueness of the null pair**. k(ξ) m(ξ) = k′(ξ) m′(ξ) as polynomials in ξ ⟹ (k′, m′) = (λk, m/λ) or (λm, k/λ)
   (`UniqueFactorizationMonoid`).
8. **Theorem 1**. 3, 4, 5 and 6 combine: the O_k = O_m = 0 hypothesis on an interval gives g′ ≡ 0 there and contradicts
   5.
9. **Theorem 2**. The four polynomial identities N = −(G² + 2c²) P(G, c) (`ring`), P monic of degree 12 so N ≠ 0
   (`Polynomial.card_roots'`), and a continuous function into a finite set on an interval is constant.

## 10. What `det_classify.md` and `two_scalars.md` should now say

**`det_classify.md` §3**, last paragraph, "A Frobenius argument for c ≠ 0, which would give failure on every member, was
not made." Replace with: "**Update from #162** ([`bessel_closure.md`](bessel_closure.md)): decided. For every member
with c ≠ 0, every ε, α, β, closure fails on every open set (*proved*, identities *checked symbolically*). The null
covectors k = a + b, m = a − b of R have k̃ ∧ dk̃ = O_k dt ∧ dx ∧ dρ and m̃ ∧ dm̃ = O_m dt ∧ dx ∧ dρ with O_k + O_m =
2g′(g′² + 4c²), so the pair is Frobenius-integrable only where g′ = 0, and there g″ = 2εc² ≠ 0 makes the zeros isolated.
For c = 0 this is #65's ν ∧ dν = b′ dx ∧ dρ ∧ dt." In the table (row "b(ρ) of #65 and its linear-in-x deformation"),
"fails at 6 exact jets (c = 0, 2, 6)" becomes "fails on every open set, every member (#162)", status **proved**,
**checked symbolically**. In its §10, add to "By hand" the Bessel closure (Theorem 1). In its §12, the suggested bullet's "fails
closure at the jets tested" becomes "fails closure on every open set (#162, *proved*)".

**`two_scalars.md`**, "Not covered", the bullet on shifts b(t, x, ρ): replace "which is not localized in x and fails
closure at the jets tested (*proved*, *checked symbolically*)" with "which is not localized in x and fails closure on
every open set ([`bessel_closure.md`](bessel_closure.md), #162, *proved*)". The Frobenius sentence for b(ρ) in "A shift
b(ρ)" is unchanged.
