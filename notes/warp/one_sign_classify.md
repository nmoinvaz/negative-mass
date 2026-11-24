# Which shifts pass the one-sign scalar stage on an open set?

November 2026, issue [#158](https://github.com/nmoinvaz/negative-mass/issues/158). Script
[`one_sign_classify.py`](one_sign_classify.py), its stdout
[`one_sign_classify_output.txt`](one_sign_classify_output.txt). It builds on
[`same_sign_scalars.md`](same_sign_scalars.md) (#78: Lemma S, Lemma K, the planted metrics), on
[`det_classify.md`](det_classify.md) (#80: the block R₃, the reductions, the mixed-sign analogue) and on
[`two_scalars.md`](two_scalars.md) (#65). The script imports `rank1_lams`, `same_sign_form`, `closure2` and `simp_r`
from `same_sign_scalars.py`, `ricci` from `phantom_source.py` and `rq` from `two_scalars.py`. The block formulas, Lemma
R, every reduction, the planted metrics and the controls are new here.

Status words as in AGENTS.md: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails),
**numerical evidence**, **conjecture**, **not classified**.

## Answer

The one-sign stage is an overdetermined system for one function. Inside each reduction it is decided; for the general
static shift it is not.

| reduction | solutions on an open set | localized in x | closure (Lemma K) | status |
|---|---|---|---|---|
| b(ρ) | b = C ln ρ + b₀, λ = h = C²/(2ρ²), two phantoms (the planted case) | no: x-independent | **fails**: kernel not integrable | **proved**, **checked symbolically** |
| planar b(t, x) | **every** shift with K = −a_x ≠ 0, λ = K: two normal scalars if K < 0, two phantoms if K > 0 | none when K is constant (Theorem P); a localized shift has K not constant and fails closure | **passes iff K is constant**: then #78's P1 (AdS₂ × ℝ²) or P2 (dS₂ × ℝ²) | **proved** (Theorem P), **checked symbolically** |
| static or traveling u(x, ρ) with u_xρ ≡ 0 | u = c x + C ln ρ + d (+ b(ρ), planar), λ = C²/(2ρ²) + c², two phantoms | no: linear in x | **fails** for C ≠ 0; passes for C = 0 (planar K = c²) | **proved**, **checked symbolically** |
| multiplicative u = X(x) R(ρ) + v | none with X′ R′ ≠ 0 (F = 0 alone) | — | — | **proved**, **checked symbolically** |
| conical u = F(x/ρ) + v | none with F′ ≠ 0 | — | — | **proved**, **checked symbolically** |
| spherical profile u = w(r) (the bubble at constant speed) | none: F ≠ 0 off finitely many angles | — | not reached | **proved** (as #78), **checked symbolically** |
| static or traveling u(x, ρ) with u_xρ ≠ 0 | resolved to a finite-type system with a compatibility condition C ≢ 0; ends are b(ρ) solutions; the tangent space at b(ρ) is w = c x + g(ρ) | **not decided** | not run | **not classified** |
| b(t, x, ρ) with t-dependence beyond the planar and the gauge x → x + D(t) | — | — | — | **not classified** |

So the answer to the issue is: **within the reductions b(ρ), planar, additive, multiplicative, conical and spherical, no
bubble-like shift (u_ρ ≢ 0, u_x ≢ 0, limits in x) passes the one-sign stage and closure.** The planar shifts pass the
stage and, when K is constant, closure. They are not bubbles (u_ρ ≡ 0), and a planar shift with constant K is not
localized. The only ρ-dependent families that pass the algebraic stage on an open set are b(ρ) and u = c x + C ln ρ + d,
both two phantoms, and closure kills both. No counterexample family is found. The general static class with u_xρ ≠ 0 and
the time-dependent shifts beyond the planar and gauge ones are open (**not classified**), so the statement "no bubble
passes" is **not** proved for them.

A shift b = c x + C ln ρ + d(t) with any d(t) is the static member u = c x + C ln ρ + d up to the gauge x → x + D(t)
(`travel`), so it passes and fails in the same way.

## Setup

ds² = −dt² + (dx + b dt)² + dρ² + ρ² dφ², n = ∂_t − b ∂_x. With (I1), (I2) of #57 the frame block R₃ on (n, ∂_x, e_ρ)
is, as in #80,

    R₃ = [[P − h, Q, Z], [Q, −P − h, Y], [Z, Y, h]],
    P = n(b_x) − b_x² = a_x,  Q = −½(b_ρρ + b_ρ/ρ),  Z = ½ b_xρ,  Y = −½ n(b_ρ) + b_x b_ρ,  h = ½ b_ρ²,

with a = b_t − b b_x. `block` re-derives it from the coordinate Ricci tensor of an undefined b(t, x, ρ) and confirms
R(e_φ, ·) = 0 and (I2) (**checked symbolically**). A static or traveling shift b = B(x − vt, ρ) has the jets of the
static u = B + v, and b(t, x − D(t)) − D′(t) has the same five jets as b at (t, x − D) (`travel`, **checked
symbolically**). So the classification is of static shifts u(x, ρ), up to the gauge x → x + D(t).

By Lemma S of #78 the stage passes at a point iff R = 0 or some real λ ≠ 0 has rank(R₃ − λη₃) = 1. R = 0 forces h = R_ρρ
= 0, so b_ρ = 0, so the shift is planar with K = 0 and flat. Where it is not flat the stage reads as follows.

**Lemma R (proved).** Let λ ≠ 0 be real.

- If h > 0, rank(R₃ − λη₃) = 1 iff there is v ∈ ℝ³ with P = −(v₀² + v₂²), Q = −v₀v₁, Z = −v₀v₂, Y = −v₁v₂, h = ½(v₀² +
  v₁²). Then λ = h + v₂² > 0, N = R₃ − λη₃ = −v vᵀ, and the pair is two phantoms.
- If h = 0, it holds iff Q = Z = Y = 0 and P ≠ 0, with λ = −P.

*Proof.* N = ε v vᵀ with ε = ±1, since N is real symmetric of rank 1. Its diagonal gives N₀₀ + N₁₁ = −2h = ε(v₀² + v₁²).
If h > 0 then ε = −1 and v₀² + v₁² = 2h. Then N₂₂ = h − λ = −v₂² gives λ = h + v₂², N₀₀ = P − h + λ = −v₀² gives P =
−(v₀² + v₂²), and the off-diagonal entries give Q, Z, Y. N₁₁ = −P − h − λ = −v₁² then holds automatically. The converse
is the same computation. If h = 0 then v₀ = v₁ = 0, so Q = Z = Y = 0, and N = diag(P + λ, −P − λ, −λ) has rank 1 iff λ =
−P ≠ 0. ∎

`lemmaR` checks N = −v vᵀ identically, F = Q(Z² + Y²) + 2hZY = 0 on it, #78's certificate R₀ρ m_d + R₁ρ m_e + F = 0 for
a generic block with (I2), and the planted forms through `rank1_lams`: 16 random v-blocks give exactly [(h + v₂², MM)],
8 planar blocks give [(−P, PP if P > 0 else MM)], P = 0 gives [], and 16 generic blocks give [] (**checked
symbolically**).

**Dichotomy.** On a box where u ∈ C², either b_ρ ≡ 0 (planar) or b_ρ ≠ 0 on a dense open subset, where h > 0. In the
second case, for a static shift, either Z = ½u_xρ ≡ 0 on a sub-box (additive, u = X(x) + R(ρ)) or Z ≠ 0 on a sub-box.
Every classification below shrinks the box in this way.

## b(ρ): the planted case, and Lemma E

With P = Z = Y = 0, Lemma R's equations give v₀² + v₂² = 0, so v₀ = v₂ = 0, so Q = 0 and λ = h. `brho_solve` finds v₀
and v₂, hence Q, in the radical of the ideal (Rabinowitsch, exact), and `dsolve` of b″ + b′/ρ = 0 gives b = C₁ + C₂ ln ρ
(**checked symbolically**). So b = C ln ρ + b₀ (**proved**), as #78 found. In 4D, at three exact points, λ = h is a root
of the minors with inertia (0, 2); the perturbed b = 2 ln ρ + ρ²/10 has no admissible root (`brho_points`). The kernel
is spanned by n and ∂_ρ, [n, ∂_ρ] = b′ ∂_x is not in it, so the plane is not integrable and closure fails
(`brho_closure`; #78's D3). Status: **proved**, **checked symbolically**.

**Lemma E (proved).** The set of (P, Q, Z, Y, h) that satisfy Lemma R's equations for some v is closed, since v₀² + v₁²
= 2h and v₂² ≤ |P| bound v. So if along a sequence of points of a solution (P, Z, Y) → 0, every limit satisfies v₀ = v₂
= 0 and Q = 0. For a shift u(x, ρ) whose second jets converge as x → ±∞, the limit profile u_±(ρ) is a solution of the
b(ρ) reduction: u_± = C_± ln ρ + b_±. The usual bubble, with u_± constant, is the case C_± = 0.

## Planar shifts: λ = K, closure iff K is constant, and none is localized

For b(t, x): h = Q = Z = Y = 0, so R₃ = diag(P, −P, 0) and N = diag(P + λ, −P − λ, −λ), which has rank 1 only at λ = −P
=: K, where the φ direction supplies the second square. So **every planar shift with K ≠ 0 passes the algebraic stage**,
with M = R − λg = −K(dy² + dz²) and the pair (+, +) iff K < 0 (`planar_stage`: the minors, and λ = K with inertia (2,
0), (0, 2), (0, 2) at three exact points each for b = 3x/2, tan(t) x, −tanh(t) x; **checked symbolically**). This is #80
§12's remark with the sign M = −K(dy² + dz²) the coordinator gave.

K = −a_x (`block`). The kernel of M is span(∂_t, ∂_x), integrable, and L_X M = −X(K)(dy² + dz²) for X in it, so Lemma K
gives: **closure holds iff K is constant** (**proved**). The quotient −K(dy² + dz²) is then flat. For K ≡ K₀ ≠ 0 the 2D
metric −dt² + (dx + b dt)² has constant curvature K₀, so the shift is #78's P1 (K₀ < 0, two normal scalars) or P2 (K₀ >
0, two phantoms) in a Painlevé–Gullstrand chart, and passes all four stages. `planar_closure`: kernel dimension 2,
integrable, L_X M = 0 and flat for the three metrics above, and L_X M ≠ 0 for b = x²/10 (K = 3x²/50) and b = 1/(1 + x²).

**Theorem P (proved).** Let b ∈ C²(I × ℝ), I an interval of positive length, with K = −a_x ≡ K₀ ≠ 0. There are no t₁ <
t* in I, with t* − t₁ < π/(2κ) if K₀ = −κ² < 0, such that b(t₁, ·) and b(t*, ·) are both bounded. So no such shift tends
to a constant, or has compact support, in x.

*Proof.* a(t, x) = b_t − b b_x has a_x = −K₀, so a = −K₀ x + α(t) with α continuous. Along a comoving curve dx/dt = −b
the chain rule gives x″ = −b_t + b b_x = −a = K₀ x − α(t) (`planar_thm`). The right side is linear in x, so the curve
through (t₁, X) exists on [t₁, t*]. Write κ² = |K₀|, Δ = t* − t₁, β(X) = −b(t₁, X) and x_p a particular solution. For K₀
= κ²,

    x_X(t*) = x_p(t*) + (X − x_p(t₁)) cosh κΔ + (β(X) − x_p′(t₁)) sinh(κΔ)/κ,
    b(t*, x_X(t*)) = −x_X′(t*) = −κ (X − x_p(t₁)) sinh κΔ − (β(X) − x_p′(t₁)) cosh κΔ − x_p′(t*).

If b(t₁, ·) is bounded then β is bounded, so |x_X(t*)| → ∞ and |b(t*, x_X(t*))| ≥ κ sinh(κΔ)|X| − const → ∞ as X → ∞. So
b(t*, ·) is unbounded. For K₀ = −κ² the same holds with cos and sin, since cos κΔ > 0 and sin κΔ > 0 for 0 < κΔ < π/2. ∎

`planar_thm` checks the chain rule, that the displayed x_X solves x″ = ±κ² x − α(t) for both signs, and that the
coefficients of X in x_X and x_X′ are cosh and κ sinh (resp. cos and −κ sin) with the ratio κ tanh(κΔ) (resp. −κ
tan(κΔ)) (**checked symbolically**, at fixed β). The bound on β(X) and the passage to b are by hand. A localized planar
shift is bounded at every t, so it has K ≢ const: it passes the algebraic stage where K ≠ 0 and fails closure.

## Additive static shifts u = X(x) + R(ρ) (+ v)

This is u_xρ ≡ 0. With S = R″ + R′/ρ the block is Z = 0, Y = X′R′, Q = −S/2, P = −u X″ − X′², h = R′²/2
(`additive_block`, **checked symbolically**). A traveling B = X(ξ) + R(ρ) is the same with v absorbed in R.

**Theorem A (proved).** Let u = X(x) + R(ρ) + v be C² on a box with R′ ≠ 0 and rank(R₃ − λη₃) = 1. Then X′ ≡ 0 (b(ρ)),
or X″ ≡ 0, S ≡ 0, and u = c x + C ln ρ + d. The shift passes the algebraic stage with λ = C²/(2ρ²) + c² and two
phantoms.

*Proof.* If X′ ≡ 0 on a sub-box, u = R(ρ) + v is b(ρ). Otherwise shrink to X′ ≠ 0. Then Y = X′R′ ≠ 0, so v₁v₂ ≠ 0, and Z
= 0 = −v₀v₂ gives v₀ = 0. Then S = 0 and u X″ = 0 (`additive_solve`: S, v₀ and u X″ lie in the radical of Lemma R's
ideal saturated by R′X′; X″ alone does not). S = 0 gives R = C ln ρ + R₀ with C ≠ 0. u = X + R + v cannot vanish on an
open set because R′ ≠ 0, so X″ = 0 on a dense open set, hence on the box. Conversely u = c x + C ln ρ + d has v = (0,
±C/ρ, ∓c) and λ = h + c². ∎

`additive_points` finds λ = C²/(2ρ²) + c² as a root with inertia (0, 2) at three exact points for (c, C, d) = (3/2, 2,
1/5) and (1, 3, 0), and no admissible root for 3x/2 + 2 ln ρ + x²/10 (**checked symbolically**).

**Closure.** The kernel of M = R − λg is span(n, K₂) with K₂ = v₂∂_x − v₁∂_ρ, v₁ = ±C/ρ, v₂ = ∓c. Since n(v_i) = 0 and
K₂(u) = v₂c − v₁C/ρ = ∓(c² + C²/ρ²),

    [n, K₂] = K₂(u) ∂_x = ∓(c² + C²/ρ²) ∂_x,

which is not in span(n, K₂) when C ≠ 0, because the ∂_t coefficient of an element of the plane forces it to be a
multiple of K₂, and K₂ is not along ∂_x. So the plane is not integrable and **closure fails for every C ≠ 0**
(**proved**; `additive_closure` finds the 2D non-integrable kernel at c = 3/2, C = 2 and c = 1, C = 3, **checked
symbolically**). For C = 0 the shift is planar with K = c² and passes (`additive_closure`: integrable, L_X M = 0, flat).
No member with c ≠ 0 is localized.

## Multiplicative shifts u = X(x) R(ρ) + v

With W = 3XR + v and S = R″ + R′/ρ,

    F = −X X′² R′² (S (1 + W²) − 2 X R′² W)/8   (`mult`, identity).

On a box with X X′ R′ ≠ 0 and R ≠ 0, X ranges over an interval, so the quadratic in X, S(1 + W²) − 2XR′²W, has all its
coefficients zero. They generate, with 1 − t R R′, the unit ideal, so there is no solution (**proved**, **checked
symbolically**). Only X′ R′ = 0 remains: b(ρ) or planar. #80 left the multiplicative case of its D = 0 problem
unclassified. Here F = 0 alone decides it.

## Conical shifts u = F(x/ρ) + v

Let U = F + v, G = F′, H = F″, η = x/ρ, w = G + ηH. At x = ηρ,

    ρ² (P, Q, Z, Y, h) = (−(G² + U H), −η w/2, −w/2, −(η G² + U w/2), η² G²/2)   (`conical_entries`),

so Q = ηZ. A traveling B = F(ξ/ρ) is the case U = F + v. With h ≠ 0 (η G ≠ 0), shrink to U ≠ 0 (U ≡ 0 gives G = 0).

**Proposition C (proved).** There is no solution with F′ ≠ 0 on an open set.

*Proof* (`conical_solve`, ideal computations exact). Lemma R's equations with η G U w ≠ 0 imply

    E_a = w (1 + U²) + 2 η U G² = 0,    E_b = (η² + 1) η H + (η² − 1) G = 0

(both reduce to 0 modulo the Gröbner basis). The branch w = 0 with η G U ≠ 0 is the unit ideal. E_b is G′ = −(η² − 1)
G/(η(η² + 1)), so G = C η/(η² + 1), C ≠ 0, and w = 2G/(η² + 1). Then E_a = (2G/(η² + 1))(1 + F² + C η² F), so 1 + F² +
Cη²F = 0 on an interval. Its derivative, with F′ = G, is linear in F and vanishes only at F = −Cη²/(2(η² + 2)). For that
F,

    F′ − G = −C η (η⁴ + 6η² + 6)/((η² + 1)(η² + 2)²) ≠ 0   for C η ≠ 0,

a contradiction. ∎

#80 left the conical form of its D = 0 problem undecided. The one-sign condition has two equations for one function, and
that is what decides it here.

## The spherical profile: the bubble at constant speed

For u = w(r), r = (x² + ρ²)^½, s = ρ/r, c = x/r, with a₁ = r w″ − w′, a₂ = w a₁ + 2 r w′², z = 2w′ + a₁ s²,

    8 r³ F = s² c² (−z (a₁² + a₂²) + 2 r a₁ a₂ w′² s²)   (`bubble_F`, identity),

which is #78's B. Its s⁰ and s² coefficients with 1 − t w′ r generate the unit ideal, and the s⁰ coefficient is −2w′(a₁²
+ a₂²) ≠ 0 over ℝ when w′ ≠ 0 (a₁ = a₂ = 0 would give a₂ = 2 r w′² = 0). So at each r only finitely many angles pass and
no open set does (**proved**, as #78's Proposition W; **checked symbolically**). `bubble_points`: at 20 exact points
`rank1_lams` is empty. This is the planted case "the bubble at constant speed must not pass". It agrees with #78, so no
discrepancy.

## General static shifts with u_xρ ≠ 0: a finite-type system, not classified

On a box with h > 0 and Z ≠ 0, Lemma R forces v₀v₂ ≠ 0. Shrink to u ≠ 0 and u_x ≠ 0 (u_x ≡ 0 on an open set would give Z
= 0). With p = Y/Z = u + 2u_xu_ρ/u_xρ = v₁/v₀ the five equations reduce to two (`static_system`: A = Q(Z² + Y²) + 2hZY
and B = PQY − Z(Q² + Y²) lie in the radical):

    u_ρρ + u_ρ/ρ = 2p u_ρ²/(1 + p²),    u u_xx + u_x² = u_ρ²/(1 + p²) + u_xρ² (1 + p²)/(4u_ρ²).

These resolve u_ρρ and u_xx in terms of (u, u_x, u_ρ, u_xρ), with u_xρ free. With a = u_x, β = u_ρ and s = tan(θ/2), tan
θ = p, set v₀ = β cos θ, v₁ = β sin θ, v₂ = a/(u cos θ − sin θ) (this solves Y = −v₁v₂ identically). Then u_xx = Φ₁,
u_xρ = Ψ, u_ρρ = Φ₃ are rational in (ρ, u, a, β, s), and the block built from them equals Lemma R's (`static_system`,
identity). Equality of the mixed derivatives of a and β gives two equations E₁, E₂, linear in (s_x, s_ρ):

    a_xρ = a_ρx,  β_xρ = β_ρx  ⟹  s_x = σ_x(ρ, u, a, β, s),  s_ρ = σ_ρ(ρ, u, a, β, s).

So a C⁴ solution is an integral surface of the rank-2 distribution spanned by X₁ = D_x, X₂ = D_ρ on the 5-dimensional
space (ρ, u, a, β, s), and [X₁, X₂] = C ∂_s with C = D_ρσ_x − D_xσ_ρ. A solution through a point exists only if C = 0
there and on the surface.

- **Planted.** On the jets of u = c x + C ln ρ + d (s = 1, v₀ = 0), E₁ = E₂ = 0 with σ = 0, and C = 0 at 6 exact points
  (`static_system`). So the machinery sees the known family, which lies in u_xρ = 0.
- **Result** (`static_compat`). At three exact points with s ∉ {0, ±1} the 2 × 2 matrix of (E₁, E₂) in (σ_x, σ_ρ) is
  nonsingular and C ≠ 0 (the exact values are in the output). So the system is not Frobenius closed: a generic
  algebraically admissible jet does not extend to a solution. Solutions with u_xρ ≠ 0 have their jets in the
  4-dimensional hypersurface C = 0 of the 5-dimensional jet space, so they form a family with at most 3 parameters
  including the translation in x (the family u_xρ = 0 attains that bound at s = 1, where C vanishes identically).
  **Whether any exists is not classified.** Further compatibilities (X₁C = X₂C = 0, and so on) were not computed: C is a
  large rational function, and a Gröbner computation of the chain was not attempted.

**Corollary of Lemma E (proved).** A localized shift, if one exists in this class, connects two b(ρ) solutions at x →
±∞.

## The linearization at b = C ln ρ + b₀

Let u = u₀ + ε w, u₀ = C ln ρ + b₀ (C u₀ ≠ 0). At ε = 0, (P, Q, Z, Y) = 0 and h = C²/(2ρ²). With v₁ ≠ 0 Lemma R gives Z
v₁² + QY = 0, P v₁² + Q² + Y² = 0 and v₁⁴ − 2h v₁² + Q² = 0 (`linearized`: in the radical, saturated by v₁). At O(ε) the
first two read (w_xρ/2) v₁² = 0 and −u₀ w_xx v₁² = 0 with v₁² = C²/ρ² (`linearized`, **checked symbolically**). So the
tangent space of the solution set at u₀ is w_xρ = 0, u₀ w_xx = 0, i.e. w = c x + g(ρ). **No first-order perturbation of
the log solution is localized in x** (**proved** for the linearized system; this says nothing about finite
perturbations).

## Closure, family by family

| family | algebraic stage | Lemma K | status |
|---|---|---|---|
| b = C ln ρ + b₀ | passes, λ = h | kernel not integrable | fails (#78's D3) |
| planar, K constant ≠ 0 | passes, λ = K | passes; #78's P1, P2; all four stages | not localized (Theorem P) |
| planar, K not constant | passes, λ = K | L_X M ≠ 0 | fails |
| u = c x + C ln ρ + d, C ≠ 0 | passes, λ = C²/(2ρ²) + c² | kernel not integrable | fails |
| multiplicative, conical, spherical | no open set | not reached | — |
| static, u_xρ ≠ 0 | not classified | not run | — |

## Checks and controls

21 checks, all pass (`one_sign_classify_output.txt`). The script exits 1 if a check fails or a control does not fire.
Per check, by section: 1 `block`, `travel`; 2 `lemmaR`; 3 `brho_solve`, `brho_points`, `brho_closure`; 4 `planar_stage`,
`planar_closure`, `planar_thm`; 5 `additive_block`, `additive_solve`, `additive_points`, `additive_closure`; 6 `mult`; 7
`conical_entries`, `conical_solve`; 8 `bubble_F`, `bubble_points`; 9 `static_system`, `static_compat`; 10 `linearized`.

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `block` | R₃ of an undefined b(t, x, ρ) against the coordinate Ricci tensor, R(e_φ, ·) = 0, (I2) | `block-Y-sign` | Y built as +n(b_ρ)/2 + b_x b_ρ |
| `travel` | traveling = static, the gauge x → x + D(t) | `travel-offset`, `gauge-drop` | u = B without + v; without − D′(t) |
| `lemmaR` | N = −v vᵀ, F = 0 on it, the certificate, 41 planted blocks | `lemmaR-sign`, `lemmaR-F` | P = +(v₀² + v₂²); F with h Z Y |
| `brho_solve` | b(ρ): v₀, v₂, Q in the radical; dsolve | `brho-Q` | Q = h claimed |
| `brho_points` | b = 2 ln ρ + 1/3 passes at 3 points; b + ρ²/10 does not | `brho-pert` | the planted solution perturbed |
| `brho_closure` | the kernel is not integrable (#78's D3) | `brho-bracket` | commutator without K₂(K₁^μ) |
| `planar_stage` | R₃ = diag(P, −P, 0), λ = −P only; 3 planar metrics at 3 points each | `planar-sign` | λ = +P expected |
| `planar_closure` | closure iff K constant, 3 pass and 2 fail | `planar-constK` | the two failing shifts replaced by b = x, 2x |
| `planar_thm` | x″ = K₀x − α, the family x_X, the growth ratio | `planar-ode` | x″ = −K₀x + α |
| `additive_block` | R₃ of X(x) + R(ρ) + v | `add-block` | P = −X″ − X′² claimed |
| `additive_solve` | S, v₀, u X″ in the radical | `add-quad` | X″ claimed in place of u X″ |
| `additive_points` | c x + C ln ρ + d passes at 3 points; + x²/10 does not | `add-pert` | the planted solution perturbed |
| `additive_closure` | non-integrable kernel for C ≠ 0, passes for C = 0 | `add-C0` | the two C ≠ 0 members replaced by u = x |
| `mult` | F identity, unit ideal | `mult-W` | W = 2XR + v |
| `conical_entries` | the five scaled entries, Q = ηZ | `conical-Q` | Q = −ηw |
| `conical_solve` | E_a, E_b reduce to 0, w = 0 unit ideal, the ODE and the contradiction | `conical-Eb` | (η² + 1) in place of (η² − 1) |
| `bubble_F` | chain-rule jets, F = s²c²B/(8r³), unit ideal | `bubble-sign` | +z |
| `bubble_points` | 20 exact points: no rank-1 λ | `bubble-planted` | planted rank-1 blocks |
| `static_system` | A, B in the radical, the parametrization identity, the planted family, C = 0 at 6 points | `compat-Phi1` | opposite sign of a² in Φ₁ |
| `static_compat` | det ≠ 0 and C ≠ 0 at 3 exact points (a result) | `compat-s1` | the points drawn on s = 1, where C = 0 |
| `linearized` | the first-order parts, the planted tangent family | `lin-claim` | +u₀ w_xx |

All controls fire on their target check and on no other: 23 of 23 (`controls: all fired`), one or two for each of the 21
checks. **One control was silent in the first batch.** `add-bracket` (the kernel commutator without K₂(K₁^μ), target
`additive_closure`) exited 0: for the additive family the nonzero part of the commutator is K₁(K₂^μ), so the mutated
commutator is still nonzero and the plane still looks non-integrable. It is reported as silent (`ctl_add-bracket.log` in
the raw folder, README amendment A6) and was replaced after seeing the data by `add-C0`, which feeds the check the
integrable C = 0 member and fires. Two controls were added after the first full run (`add-block`, `compat-s1`) because
two checks had none (README A7). The original `brho-bracket` (D3) fired. The other 20 controls of the first batch fired
on their named target alone, and so did the two added later.

`planar_thm`'s statement in the README as worded holds at a fixed initial velocity; the bound on β(X) in Theorem P is by
hand (README A1). Failed development runs are kept in the raw folder and described in its README (the float zero in
`blocks`, the positivity test in `conical_solve`, a zero draw in `static_compat`): all were bugs in the check code,
found before any run of record, and none changed a criterion.

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/one_sign_classify.py                 # the checks, then every control
    python3 notes/warp/one_sign_classify.py --no-controls   # the 21 checks
    python3 notes/warp/one_sign_classify.py --mutate=NAME   # one control

The checks take 26 s on a Raspberry Pi 5 (aarch64, one process at nice 10, shared with other jobs) and the whole run
with its 23 controls 203 s, at 82 MB peak resident for the largest process (`/usr/bin/time -v`: 99% of one core).

## What the code proves, and what is assumed

The code proves, as exact identities or exact ideal computations: the block of R₃ for an undefined b and the gauge and
traveling reductions; N = −v vᵀ, F = 0 on it and #78's certificate; for b(ρ), the radical membership of v₀, v₂, Q; the
planar block and the minors; the additive block, the radical membership of S, v₀ and u X″; the F factorization and the
unit ideal for multiplicative shifts; the conical entries, the ideal membership of E_a, E_b, the unit ideal of the w = 0
branch, and the ODE and contradiction; the spherical identity and its unit ideal; the radical membership of A, B, the
parametrization, the planted compatibility and C ≠ 0 at three exact points; the first-order parts at the log solution.
Exact outcomes at exact points: `brho_points`, `planar_stage`, `additive_points`, `bubble_points`, `lemmaR` (d).

By hand: Lemma R, the dichotomy, Lemma E, the closure criterion for planar shifts, Theorem P (the bound on β(X) and the
passage to b), Theorem A's last step (u ≠ 0 on a dense set), the kernel and commutator of the additive family, the case
analysis of Proposition C, the passage from "C ≠ 0" to "not Frobenius closed", the dimension bound, the linearization.

Trusted: sympy 1.14.0; `ricci` (reviewed in #59); `rank1_lams`, `same_sign_form` and `closure2` of #78 (reviewed in
#78's coordinator review, with `ind78.py`). The checks that call them are not independent of #78. `block` shares `ricci`
with every earlier note.

Assumed: two real, minimally coupled scalars of one sign with a constant diagonal kinetic form and any potential; G_ab =
8πT_ab; axisymmetric shifts b ∈ C² (C⁴ for the compatibility argument, C³ in Theorem P's ODE step) with unit lapse and
flat slices; the local leaf space of the kernel regular where Lemma K is used. Not independent of #78: the planted
metrics and the closure code. No independent check (`ind158.py`) was run: that is for the coordinator.

## Lemmas stated for Lean (#135)

1. **(R)** Lemma R above: for real symmetric R₃ = [[P − h, Q, Z], [Q, −P − h, Y], [Z, Y, h]], h > 0, λ ≠ 0: rank(R₃ −
   λη₃) = 1 iff ∃ v with the five relations, λ = h + v₂². Status: proved (polynomial identities, `ring`, and the sign of
   ε).
2. **(F)** F = Q(Z² + Y²) + 2hZY = 0 when rank(R₃ − λη₃) ≤ 1: #78's lemma (F). Status: proved.
3. **(B)** With h > 0, Z Y Q ≠ 0: PQY = Z(Q² + Y²). Status: proved.
4. **(E)** The set of (P, Q, Z, Y, h) of Lemma R is closed, and (P, Z, Y) = 0 on it gives Q = 0. Status: proved.
5. **(P)** Theorem P: for x″ = K₀x − α(t) with bounded initial velocity, the displayed growth of x′ in X. Status:
   proved.
6. **(A)** Additive: Z = 0, X′R′ ≠ 0 gives v₀ = 0, S = 0, u X″ = 0. Status: proved.
7. **(M)** Multiplicative: a polynomial of degree ≤ 2 in X vanishing on an interval has zero coefficients, and the three
   coefficients with 1 − t R R′ are a unit ideal. Status: checked symbolically (the Gröbner certificate); give explicit
   cofactors in Lean.
8. **(C)** Conical: E_a and E_b are in the ideal; the ODE and the contradiction. Status: checked symbolically.
9. **(K)** Lemma K, as #78. Status: proved.

## Not covered, and what `same_sign_scalars.md`'s "Not covered" should now say

Not classified here: static or traveling shifts with u_xρ ≠ 0 (§General static); time-dependent shifts other than planar
b(t, x) and the gauge x → x + D(t); the Galilean-dilation similarity forms b = μ ln ρ + F((x − μ t ln ρ)/ρ, t/ρ), which
contain the planted b(ρ) and are a PDE in two variables; non-axisymmetric shifts; three or more scalars; non-minimal
coupling.

Replace the first bullet of `same_sign_scalars.md`'s "Not covered" ("A non-axisymmetric shift and a general shift b(t,
x, ρ): Lemma S and the certificate (F) apply to them, but #80's classification of the shifts that pass on an open set is
not done here.") by:

> A non-axisymmetric shift, and a general shift b(t, x, ρ) beyond the reductions of
> [`one_sign_classify.md`](one_sign_classify.md) (#158). There, with (I1), (I2), the stage holds iff R = 0 or λ ≠ 0 and
> rank(R₃ − λη₃) = 1, and for h > 0 iff R₃ − λη₃ = −v vᵀ (Lemma R), two phantoms. Within b(ρ), planar b(t, x), additive
> u = X(x) + R(ρ), multiplicative, conical and spherical profiles it is decided (*proved*, with the identities *checked
> symbolically*): the ρ-dependent families that pass are b = C ln ρ + b₀ and u = c x + C ln ρ + d, and closure fails for
> both; every planar shift passes with λ = K = −a_x, closure holds iff K is constant (#78's P1, P2), and no such shift
> is bounded in x at two close times. Static shifts with u_xρ ≠ 0 and time-dependent shifts beyond the planar and the
> gauge x → x + D(t) are *not classified*.

## References

None beyond the notes cited. No new claim of novelty is made, and no literature search was done.

## Discrepancies

- `add-bracket` was a silent control, replaced by `add-C0` after the data (above).
- Three development runs failed on bugs in the check code (above).
- Nothing in #78, #80 or #65 disagrees with this note. #80 §12's remark for #158 (planar shifts pass with λ = K) is
  confirmed, with the sign M = −K(dy² + dz²) the coordinator gave. #80's mixed-sign multiplicative and conical cases are
  still open there; the one-sign cases are closed here by F = 0 alone.

## Recommended tasks

- **Finish the static class u_xρ ≠ 0.** Prolong the finite-type system of §General static: compute X₁C, X₂C and their
  consequences on C = 0 until the constraint set is empty or a family appears. Acceptance: either a Gröbner or resultant
  certificate that no real jet with s ∉ {0, ±1} survives, or an explicit solution that passes Lemma R at exact points;
  planted: the s = 1 family; a control that breaks Φ₁. Depends on this note. Cost: C is a large rational function, so
  work modulo a prime or with a random-point rank test first.
- **Time-dependent b(t, x, ρ) with Z ≡ 0, b = X(t, x) + R(t, ρ).** The rank-1 equations are a short ODE system in t.
  Acceptance: the solution set, with closure, or a proof that every solution is gauge-equivalent to u = c x + C ln ρ +
  d. Depends on Lemma R.
- **The Galilean-dilation similarity form.** b = μ ln ρ + F((x − μ t ln ρ)/ρ, t/ρ), a two-variable PDE that contains
  b(ρ). Acceptance: its rank-1 system, and whether F ≡ const is the only solution.
- **Independent check `ind158.py`** (coordinator): Lemma R, the additive, multiplicative and conical eliminations and
  Theorem P's growth, by code that shares nothing with this script.
- **The Lean statements** of §Lemmas for Lean are rows of #135.
