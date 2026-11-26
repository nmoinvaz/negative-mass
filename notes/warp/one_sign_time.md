# Time-dependent additive shifts b = X(t, x) + R(t, ρ) at the one-sign stage

November 2026, issue [#166](https://github.com/nmoinvaz/negative-mass/issues/166). Script
[`one_sign_time.py`](one_sign_time.py), its stdout [`one_sign_time_output.txt`](one_sign_time_output.txt). It builds on
[`one_sign_classify.md`](one_sign_classify.md) (#158: the block R₃, Lemma R, the Dichotomy, planar shifts, Theorem P,
the static additive class) and on [`same_sign_scalars.md`](same_sign_scalars.md) (#78: Lemma S, Lemma K). The script
imports `blocks`, `block3` and `v_block` from `one_sign_classify.py`; `same_sign_form`, `rank1_lams`, `closure2` and
`simp_r` from `same_sign_scalars.py`; `ricci` from `phantom_source.py`; `rq` from `two_scalars.py`. Lemma Z, the solve
of each branch, the families, the identity of the whole family and the controls are new here. The criteria were written
in the raw folder (`negative-mass-raw/I166-rpi5/README.md`) before the script existed; the code tests them as worded,
with the alignments of its amendment A1.

Status words as in AGENTS.md: **proved** (by hand), **checked symbolically** (sympy, exact, with a control that fails),
**numerical evidence**, **conjecture**, **not classified**.

## Answer

On a box with Z = ½ b_xρ ≡ 0, b = X(t, x) + R(t, ρ). Either R_ρ ≡ 0 (planar, decided in #158) or h = ½ R_ρ² > 0, and
then the one-sign stage holds on an open set iff

    b = c(t) x + C(t) ln ρ + d(t),   C ≠ 0,   c′ = 2γ c − γ²,   γ = C′/(2C),

with λ = C²/(2ρ²) + (γ − c)² > 0, two phantoms, and R_ab = λ(g_ab − ρ² dφ_a dφ_b) − α_a α_b for α = (C/ρ)(dx + b dt) +
(γ − c) dρ (**proved**; the identity **checked symbolically** for the whole family). Closure fails for every member, and
none is localized in x.

**The solution set is not the static family up to gauge.** C(t) is any C¹ function without zero, k₀ any constant, and c
= |C| (k₀ − K) with K′ = γ²/|C|. The function d(t) is pure gauge (D′ + cD = d). c(t) = b_x and C(t) = ρ b_ρ are gauge
invariants, so a member with C′ ≢ 0 is gauge-equivalent to no static member. The static members are C′ = 0, which forces
γ = 0 and c′ = 0 (#158's u = c x + C ln ρ + d).

| branch (Lemma Z) | solutions on an open set | closure (Lemma K) | localized in x | status |
|---|---|---|---|---|
| planar, R_ρ ≡ 0 | all shifts with K = −a_x ≠ 0 (#158) | passes iff K is constant (#158) | never (Theorem P of #158) | **proved** (#158) |
| A: v₀ = 0 | b = c(t) x + C(t) ln ρ + d(t), c′ = 2γc − γ², one free function C(t), one constant k₀, and d(t) | **fails**: kernel plane not integrable | no: affine in x | **proved**, **checked symbolically** |
| B: v₂ = 0, v₀ ≠ 0 somewhere | none | — | — | **proved**, **checked symbolically** |
| B ∩ A: v₀ = v₂ = 0 | the sub-family γ = c (so c′ = c²), e.g. b = −x/t + t⁻² ln ρ | **fails** | no | **proved**, **checked symbolically** |

So within Z ≡ 0, the shifts that pass the stage and closure are the planar ones with constant K, which are never
localized. No bubble-like shift (u_ρ ≢ 0, limits in x) passes the one-sign stage and closure with Z ≡ 0. The
time-dependent ρ-dependent shifts with Z ≢ 0 are still **not classified**.

## Setup

ds² = −dt² + (dx + b dt)² + dρ² + ρ² dφ², n = ∂_t − b∂_x. A C² shift on a box has Z = ½ b_xρ ≡ 0 iff b = X(t, x) + R(t,
ρ). With S = R_ρρ + R_ρ/ρ, #158's block R₃ = [[P − h, Q, Z], [Q, −P − h, Y], [Z, Y, h]] has

    P = X_tx − (X + R) X_xx − X_x²,   Q = −S/2,   Z = 0,   Y = −R_tρ/2 + X_x R_ρ,   h = R_ρ²/2.

`block` re-derives it from the coordinate Ricci tensor of the metric with undefined X(t, x) and R(t, ρ), checks R(e_φ,
·) = 0 and (I2), compares it with #158's `blocks` on the jets, and checks that R = r(t) gives diag(P, −P, 0) (**checked
symbolically**). The gauge x → x + D(t), b → b(t, x − D) − D′, keeps the class (X̃ = X(t, x − D) − D′, R̃ = R) and the
block at (t, x, ρ) equals that of b at (t, x − D, ρ) (`gauge`, undefined X, R, D, **checked symbolically**).

**Regularity.** b ∈ C² for the algebraic stage and for Lemma Z, Branch A and Branch B (b_tx and b_tρ exist and are
continuous). The closure computation differentiates K₂ = (γ − c)∂_x − (C/ρ)∂_ρ along n, so it needs C ∈ C² (b ∈ C³).

**Dichotomy.** On a box either R_ρ ≡ 0 (planar, #158) or R_ρ ≠ 0 on a sub-box, where h > 0. The set {h > 0} is open and
the interior of {h = 0} is planar, so the two cases cover an open dense subset.

## Lemma Z (proved)

Let h > 0 and Z = 0. By Lemma R, rank(R₃ − λη₃) = 1 for a real λ ≠ 0 iff R₃ − λη₃ = −v vᵀ, which gives Z = −v₀v₂ = 0. So
v₀ = 0 or v₂ = 0.

- **(A)** v₀ = 0: P = −v₂², Q = 0, Y = −v₁v₂, 2h = v₁². Then v₁² = 2h, 2hP + Y² = −v₁²v₂² + v₁²v₂² = 0, and λ = h + v₂²
  = h − P. Conversely Q = 0 and 2hP + Y² = 0 give v₁ = √(2h), v₂ = −Y/v₁.
- **(B)** v₂ = 0: P = −v₀², Q = −v₀v₁, Y = 0, 2h = v₀² + v₁². Then P ≤ 0, 2h + P = v₁² ≥ 0 and Q² = v₀²v₁² = −P(2h + P),
  and λ = h. Conversely Y = 0, P ≤ 0, 2h + P ≥ 0 and Q² + P(2h + P) = 0 give v₀ = √(−P), v₁ = ±√(2h + P) with Q = −v₀v₁.

`lemmaZ` computes the elimination ideals of Lemma R's equations at v₀ = 0 and at v₂ = 0 (lex Gröbner bases): ⟨Q, 2hP +
Y²⟩ and ⟨Y, Q² + P(2h + P)⟩, equal to the claimed ones, and the sign-flipped claims are not (**checked symbolically**).
The inequalities are by hand. `rank1_lams` of #78 (which does not use Lemma R) returns exactly [(h + v₂², MM)] on 12
A-blocks, [(h, MM)] on 12 B-blocks and [] on 12 generic Z = 0 blocks.

(A) and (B) are closed conditions on the jets, so on {h > 0} = A ∪ B the union of the interiors of A and B is dense. On
an open dense subset of {h > 0} every point has a neighborhood where (A) holds throughout or (B) holds throughout. The
classification below is on such neighborhoods.

## Branch A: b = c(t) x + C(t) ln ρ + d(t) (proved)

Let (A) hold on a box. Q = 0 is R_ρρ + R_ρ/ρ = 0, so R = C(t) ln ρ + d(t) (`rho_ode`: `dsolve` returns C₁ + C₂ ln ρ,
**checked symbolically**), with C = ρ R_ρ ≠ 0 because h > 0, and C′ = ρ R_tρ. Then v₁ = ±C/ρ and, with R_tρ = C′/ρ,

    ρ²(2hP + Y²) = C² P + (X_x C − C′/2)²,   P = X_tx − (X + C ln ρ + d) X_xx − X_x²

(`branchA` (i)). This is the whole of (A) once Q = 0. As a polynomial in ℓ = ln ρ it has the ℓ-coefficient −C³ X_xx. It
vanishes for all ρ in an interval, so X_xx = 0 (C ≠ 0), hence X = c(t) x + e(t) with c ∈ C¹. The ℓ⁰ part becomes

    C² c′ − c C C′ + C′²/4 = C² (c′ − 2γc + γ²) = 0,   γ = C′/(2C)

(`branchA` (ii), (iii)). With e absorbed in d this is the family above. On it v₂ = −Y/v₁ = γ − c and h − P = h + v₂² =
C²/(2ρ²) + (γ − c)² (`branchA` (iv)).

**The ODE.** c′ − 2γc + γ² = 0 is linear in c. With A = |C| (A′ = 2γA, since A′/A = C′/C) and K′ = γ²/A, c = A (k₀ − K)
satisfies it: c′ − 2γc + γ² = A′(k₀ − K) − A K′ − 2γA(k₀ − K) + γ² = 0 (`c_ode`, **checked symbolically**; `dsolve` of
the ODE returns only an implicit integral relation, so the closed form is verified by substitution). For C′ = 0, γ = 0
and c′ = 0.

**The whole family satisfies the stage.** `identity` computes the Ricci tensor of ds² for b = c(t) x + C(t) ln ρ + d(t)
with c, C, d undefined, substitutes c′ = 2γc − γ², and finds R_ab − [λ(g_ab − ρ² dφ_a dφ_b) − α_a α_b] = 0 identically,
with M n = 0 and M K₂ = 0 for M = −α⊗α − λρ² dφ⊗dφ (**checked symbolically**). With the sign-flipped ODE the difference
is nonzero (control).

**Gauge.** d(t) is removed by x → x + D(t): b̃ = c x + C ln ρ + (d − cD − D′), and D = e^{−∫c} ∫ d e^{∫c} gives D′ + cD
= d (`gauge` (b), **checked symbolically**). c = b_x and C = ρ b_ρ do not change under the gauge (`gauge` (a), the
entries at (t, x, ρ) equal those of b at (t, x − D, ρ)), so the family modulo gauge is parametrized by (C(t), k₀).

| member | b | C | γ | c | γ − c | λ |
|---|---|---|---|---|---|---|
| P1, static | 3x/2 + 2 ln ρ + 1/5 | 2 | 0 | 3/2 | −3/2 | 2/ρ² + 9/4 |
| P3, gauge copy, d = tanh t | 3x/2 + 2 ln ρ + tanh t | 2 | 0 | 3/2 | −3/2 | 2/ρ² + 9/4 |
| P4, gauge copy, D = tanh t | 3(x − tanh t)/2 + 2 ln ρ + 1/5 − (1 − tanh²t) | 2 | 0 | 3/2 | −3/2 | 2/ρ² + 9/4 |
| T1, new | x/2 + e^{2t} ln ρ | e^{2t} | 1 | 1/2 | 1/2 | e^{4t}/(2ρ²) + 1/4 |
| T2, new | (1/2 + e^{2t}/3) x + e^{2t} ln ρ | e^{2t} | 1 | 1/2 + e^{2t}/3 | 1/2 − e^{2t}/3 | e^{4t}/(2ρ²) + (1/2 − e^{2t}/3)² |
| T3, new | −tanh t · x + e^{2t}(1 − tanh²t) ln ρ | e^{2t}(1 − tanh²t) | 1 − tanh t | −tanh t | 1 | C²/(2ρ²) + 1 |
| T4, new, v₂ = 0 | −x/t + t⁻² ln ρ | t⁻² | −1/t | −1/t | 0 | t⁻⁴/(2ρ²) |

Each is a root of the minors of R − λg in 4D with inertia (0, 2), two phantoms, at 3 exact points (`planted`, `family`),
and has Q = 2hP + Y² = 0 identically. The ODE-violating members F1 = P1 + x²/10 and F2 = (1/2 + e^{2t}/3 + 1/5) x +
e^{2t} ln ρ have no root with an admissible pattern at 3 of 3 points, and 2hP + Y² ≠ 0. A draw with tanh t = ±1 gives C
= 0 for T3; the script redraws C = 0 points (the first probe hit one).

## Branch B: no solution with v₀ ≠ 0 (proved)

Let (B) hold on a box, and suppose P < 0 at a point; shrink to a box with P < 0, where v₀ ≠ 0.

1. Y = 0 means X_x R_ρ = R_tρ/2. Since R_ρ ≠ 0, X_x = R_tρ/(2R_ρ) is a function of (t, x) and of (t, ρ), so X_x = κ(t)
   and X = κ(t) x + e(t), and R_tρ = 2κ R_ρ. For each fixed ρ, G(t) = R_ρ(t, ρ) solves G′ = 2κG, so R_ρ = m(t) f(ρ) with
   m > 0, m′ = 2κ m (`t_ode`: `dsolve` returns C₁ exp(2∫κ), **checked symbolically**).
2. Then P = κ′ − κ² = −p with p := κ² − κ′ > 0 (p = v₀²), Q = −m(f′ + f/ρ)/2 and h = m²f²/2. With g := (f′ + f/ρ)²/4 the
   equation Q² + P(2h + P) = 0 of (B) is (`branchB` (i), **checked symbolically**)

       m² g − p m² f² + p² = 0   for all (t, ρ) in the box.

3. If f² is not constant, take ρ_a, ρ_b with F_a = f(ρ_a)² ≠ F_b. For any two times the four equations (two radii, two
   times) give p(t₁) = p(t₂) and m(t₁)² = m(t₂)² (`branchB` (ii): p₁ − p₂ and m₁² − m₂² lie in the ideal saturated by
   p₁p₂m₁m₂(F_a − F_b), **checked symbolically**). So p and m are constant, κ = m′/(2m) = 0, and p = κ² − κ′ = 0,
   contradicting p > 0. (Only continuity of p is used: no derivative of p.) A control that claims m₁ − m₂ in the ideal
   fails: only m² is determined.
4. If f² is constant, f = f₀ ≠ 0 and g = f₀²/(4ρ²) depends on ρ while p f₀² − p²/m² does not: subtracting the equations
   at ρ_a ≠ ρ_b gives f₀²(1/ρ_a² − 1/ρ_b²) = 0, a contradiction (`branchB` (iii): the ideal with the saturation p m f₀
   (ρ_a − ρ_b)(1/ρ_a + 1/ρ_b) is the unit ideal, **checked symbolically**).

So P ≥ 0, with P ≤ 0 from (B): P ≡ 0 on the interior of B, hence Q² = −P(2h + P) = 0, Q = 0, Y = 0, and (A) holds there.
Branch B adds nothing to branch A. Its overlap with A is v₀ = v₂ = 0, that is γ = c, and then c′ = 2γc − γ² = c²: c =
−1/(t − t₀) or c = 0 (`family`: T4 has γ = c and λ = h).

## Closure (proved, `closure_general`, `closure_points`)

On branch A, M = R − λg = −α⊗α − λρ² dφ⊗dφ has rank 2 (`identity`) and the kernel plane field 𝒦 = ker α ∩ ker dφ =
span(n, K₂), K₂ = w ∂_x − (C/ρ) ∂_ρ, w = γ − c. The fields are independent, α(n) = 0 and α(K₂) = 0. Since b_x = c and
b_ρ = C/ρ,

    [n, K₂] = (w′ + wc − C²/ρ²) ∂_x − (C′/ρ) ∂_ρ = (C′/C) K₂ + r ∂_x,    r = w′ + wc − 2γw − C²/ρ².

An element of 𝒦 with no ∂_t and no ∂_ρ component is zero, so [n, K₂] ∈ 𝒦 iff r = 0. But w′, wc and γw are functions of t
alone and ∂_ρ r = 2C²/ρ³ ≠ 0, so r vanishes at no more than one ρ for each t and on no open set. So 𝒦 is not integrable
on any open set, and by Lemma K there is no Λ(x) ∈ O(2) that makes (a, b) closed: **closure fails for every member**,
whatever c(t), C(t), d(t) are. The ODE is not used. `closure_general` (a) computes [n, K₂] symbolically with undefined
c, C and finds [n, K₂] − (C′/C)K₂ = (0, r, 0, 0) with r as stated and ∂_ρ r = 2C²/ρ³. (b) is a second route through
forms, sharing α with (a) and nothing else: the coefficient of dt∧dx∧dρ∧dφ in α∧dα∧dφ is F = C³/ρ³ + (C/ρ)(−w′ + (2γ −
c)w), and ∂_ρ(ρF/C) = −2C²/ρ³ ≠ 0, so α∧dα∧dφ ≠ 0 (**checked symbolically**). The step needs C ∈ C² (b ∈ C³), because w′
contains C″.

`closure_points` runs #78's `closure2` (kernel from the nullspace of M, with no formula for K₂) at P1, P3, P4, T1, T2,
T3 and T4: kernel dimension 2, [K₁, K₂] ≠ 0 for each (the nonzero bracket components are in the output). The positive
control of the same call, the planar b = x (C = 0, K = 1), is integrable with L_X M = 0 and a flat quotient (#158's A0).

## No member is localized in x (proved)

b_xx = 0 and b_x = c(t) for every member (`localized`). So b(t, ·) is affine in x with slope c(t) and bounded on a whole
line only when c(t) = 0. Setting c = 0 in c′ = 2γc − γ² gives γ² = 0, so C′ = 0: the only member with c ≡ 0 on a
t-interval is b = C ln ρ + d(t), which is x-independent and static up to the gauge D′ = d (#158's b(ρ), the planted
case). The members with C′ ≠ 0 or c ≠ 0 are not localized, and on a box that reaches x → ±∞ the limits of a member exist
only for c(t) = 0 and are then equal. A localized shift is not affine in x on a box that reaches from its core to its
limits, so it has Z ≢ 0 there (cf. #158's Lemma E).

## Checks and controls

14 checks, all pass (`one_sign_time_output.txt`); the script exits 1 if a check fails or a control does not fire.

| check | what it tests | control `--mutate=` | the control's change |
|---|---|---|---|
| `block` | R₃ of X(t, x) + R(t, ρ) from the coordinate Ricci tensor, (I2), `blocks` of #158, planar limit | `block-P` | P without −X_x² |
| `gauge` | (a) the entries at (t, x, ρ) of the gauge copy equal those of b at (t, x − D); (b) d̃ = d − cD − D′, D solves D′ + cD = d | `gauge-drop`, `gauge-ode` | X̃ without −D′; D′ − cD = d |
| `lemmaZ` | the two elimination ideals; 36 planted blocks through `rank1_lams` | `lemmaZ-A`, `lemmaZ-B`, `lemmaZ-blocks` | ⟨Q, 2hP − Y²⟩; ⟨Y, Q² − P(2h + P)⟩; B-blocks with P = +v₀² |
| `rho_ode` | `dsolve` of R″ + R′/ρ = 0 | `rho-dsolve` | R″ + 2R′/ρ = 0 solved |
| `branchA` | Q = 0, the identity ρ²(2hP + Y²) = …, the ℓ-coefficient, the ODE, v₂ and λ | `A-P`, `A-ode` | P without R; −C′²/4 |
| `c_ode` | c = A(k₀ − K) solves the ODE; static; T2 | `c-ode` | K′ = −γ²/A |
| `t_ode` | `dsolve` of G′ = 2κG | `t-dsolve` | exp(∫κ) claimed |
| `branchB` | E_B identity; separation (p₁ = p₂, m₁² = m₂²); the f-constant unit ideal | `B-claim`, `B-const` | m₁ − m₂ claimed; 1/ρ² dropped |
| `planted` | P1 to P4 through `run_family`; F1 rejected | `plant-pert` | F1 among the members that must be found |
| `family` | T1 to T4 through the same code; C′ ≠ 0; the ODE; T4; F2 rejected | `fam-ode` | F2 among the members that must be found |
| `identity` | R_ab = λ(g − ρ²dφ²) − αα for the whole family | `ident-ode`, `ident-alpha` | c′ = 2γc + γ²; (c − γ)dρ |
| `closure_general` | [n, K₂] and α∧dα∧dφ | `closure-bracket`, `closure-form` | no K₂(n^μ) term; α without b dt |
| `closure_points` | `closure2` at seven members, positive control b = x | `closure-C0` | two members replaced by b = x |
| `localized` | b_xx = 0; c = 0 forces γ = 0; d̃ = d − D′ | `loc-claim` | c′ = 2γc − γ²c |

All 21 controls fire on their named check and on no other (`controls: all fired`). **One control was silent in the first
dev batch:** `gauge-drop` exited 0, because the check code `to_jets` replaced X(t, x) inside its own derivatives, so
every entry of `gauge` (a) was 0 and the check compared 0 = 0 in the first two dev runs. This is a bug in the check
code, not in the criterion (README A1). It was fixed (derivatives are replaced first, and the script asserts that the
compared P entries contain the jet symbol X₁₁), after which `gauge-drop` fires on `gauge`; the other 20 controls fired
in that first batch. All logs are kept.

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/one_sign_time.py                 # the checks, then every control
    python3 notes/warp/one_sign_time.py --no-controls   # the 14 checks
    python3 notes/warp/one_sign_time.py --mutate=NAME   # one control

The checks take 19 s on a Raspberry Pi 5 (aarch64, one process at nice 10, shared with other jobs) and the whole run
with its 21 controls 3 min 10 s, at 77 MB peak resident (`/usr/bin/time -v`, 99% of one core). The script needs no
Singular and runs on macOS.

## What the code proves, and what is assumed

The code proves, as exact identities, exact ideal computations or `dsolve` and substitution checks: the block for X + R
and its gauge behavior; the two elimination ideals of Lemma Z and the planted blocks; the ρ-ODE and t-ODE solutions; the
identity for ρ²(2hP + Y²), its ℓ-coefficient, the ODE and v₂; the closed form of c; the E_B identity, the separation
ideal and the f-constant unit ideal; the Ricci identity of the whole family; [n, K₂] and α∧dα∧dφ; `closure2` at seven
members; b_xx = 0 and the c = 0 consequence. Exact outcomes at exact points: `planted`, `family`, `lemmaZ` (c). The
inertia of `same_sign_form` comes from 60-digit eigenvalues of rational matrices (as in #78 and #158).

By hand: Lemma Z's inequalities and the density of the interiors of A and B; the passage from "X_xx ℓ-coefficient" to
X_xx = 0 on an interval; the solution of the ODE as a linear ODE; Branch B's steps 1, 3 (the passage from two times to
all times, and m constant ⟹ κ = 0 ⟹ p = 0) and 4's sign facts; the closure argument (that [n, K₂] ∈ 𝒦 iff r = 0 and that
r is not zero on any open set); that ∂_ρ r ≠ 0 contradicts integrability on an open set; Lemma K (as #78).

Shared machinery: the 4D Ricci tensors in `block`, `identity`, `run_family` and `closure_points` all come from `ricci`
(reviewed in #59), as in every earlier note. The algebraic stage at exact points (`same_sign_form`, minors of R − λg)
uses the 4D Ricci tensor and no Lemma R or Lemma Z, so it is independent of the branch equations; `rank1_lams` likewise.
The two closure routes of `closure_general` share α; `closure2` shares no formula for K₂ with them. No independent check
(`ind166.py`) was run: that is for the coordinator.

Assumed: two real, minimally coupled scalars of one sign with a constant diagonal kinetic form and any potential; G_ab =
8πT_ab; axisymmetric shifts b ∈ C² (C³ for closure) with unit lapse and flat slices; the local leaf space of the kernel
regular where Lemma K is used. Not independent of #78 and #158: `closure2`, `same_sign_form`, `rank1_lams`, `blocks`,
`v_block`.

## Lemmas stated for Lean (#135)

1. **(Z)** For real symmetric R₃ with Z = 0 and h > 0: ∃ λ ≠ 0 with rank(R₃ − λη₃) = 1 iff (A) Q = 0 ∧ 2hP + Y² = 0, or
   (B) Y = 0 ∧ P ≤ 0 ∧ 2h + P ≥ 0 ∧ Q² + P(2h + P) = 0. Status: proved (Lemma R, polynomial identities).
2. **(Sep)** If g(ρ) = p(t) F(ρ) − B(t) for all t, ρ in a product of intervals and F is not constant, then p and B are
   constant (two radii with F_a ≠ F_b determine p and B). Status: proved.
3. **(Ell)** If a + b ln ρ = 0 for ρ in an interval, with a and b independent of ρ, then a = b = 0. Status: proved.
4. **(ODE)** c′ = 2γc − γ² with A′ = 2γA, A > 0, K′ = γ²/A has the solution c = A(k₀ − K), and every solution is of this
   form (a linear ODE). Status: proved.
5. **(K)** With r = w′ + wc − 2γw − C²/ρ² and C ≠ 0, r vanishes on no open set (∂_ρ r = 2C²/ρ³). Status: proved.

## Changes to other notes

**`one_sign_classify.md`, Answer table.** Replace the last row ("b(t, x, ρ) with t-dependence beyond the planar and the
gauge x → x + D(t)") by two rows:

| reduction | solutions on an open set | localized in x | closure (Lemma K) | status |
|---|---|---|---|---|
| time-dependent additive b = X(t, x) + R(t, ρ), Z ≡ 0, R_ρ ≢ 0 | b = c(t) x + C(t) ln ρ + d(t), C ≠ 0, c′ = 2γc − γ², γ = C′/(2C), λ = C²/(2ρ²) + (γ − c)², two phantoms; d(t) is gauge; a member with C′ ≢ 0 is gauge-equivalent to no static one; v₂ = 0, v₀ ≠ 0 has no solution ([`one_sign_time.md`](one_sign_time.md), #166) | no: affine in x | **fails** for every member (kernel not integrable) | **proved**, **checked symbolically** |
| b(t, x, ρ) with Z ≢ 0 and t-dependence beyond the planar and the gauge | — | — | — | **not classified** |

**Answer paragraph.** In "So the answer to the issue is", replace "The only ρ-dependent families that pass the algebraic
stage on an open set are b(ρ) and u = c x + C ln ρ + d, both two phantoms, and closure kills both." by "The only
ρ-dependent families that pass the algebraic stage on an open set with Z ≡ 0 are b = c(t) x + C(t) ln ρ + d(t) with c′ =
2γc − γ² (it contains b(ρ) and u = c x + C ln ρ + d as the case C′ = 0), all two phantoms, and closure kills all of
them." and replace "the time-dependent shifts beyond the planar and gauge ones are open" by "the time-dependent shifts
with Z ≢ 0 are open". The sentence after the table ("A shift b = c x + C ln ρ + d(t) with any d(t) is the static member
… up to the gauge") stays true; add: "A shift b = c(t) x + C(t) ln ρ + d(t) with C′ ≢ 0 is not gauge-equivalent to a
static member, and fails closure too (#166)."

**"Not covered" of `one_sign_classify.md`.** Replace "time-dependent shifts other than planar b(t, x) and the gauge x →
x + D(t)" by "time-dependent shifts with Z ≢ 0 other than the planar ones and the gauge x → x + D(t) (the Z ≡ 0 class is
decided in #166)".

**"Not covered" of `same_sign_scalars.md`.** In the first bullet (the #158 replacement), replace the last sentence
("Static shifts with u_xρ ≠ 0 and time-dependent shifts beyond the planar and the gauge x → x + D(t) are *not
classified*.") by: "Time-dependent additive shifts X(t, x) + R(t, ρ) (Z ≡ 0) are decided in
[`one_sign_time.md`](one_sign_time.md) (#166): the ρ-dependent ones that pass are b = c(t) x + C(t) ln ρ + d(t) with c′
= 2γc − γ², γ = C′/(2C), two phantoms, and closure fails for every member (*proved*, with the identities *checked
symbolically*). Static shifts with u_xρ ≠ 0 (conditionally excluded in #165) and time-dependent shifts with Z ≢ 0 beyond
the planar and the gauge x → x + D(t) are *not classified*."

## References

None beyond the notes cited. No new claim of novelty is made, and no literature search was done.

## Discrepancies

- The issue's second alternative ("every solution is gauge-equivalent to the static family") is false for this class:
  the members with C′ ≢ 0 are not, and they pass the algebraic stage. The first alternative (the solution set with
  closure on each family) holds.
- `gauge-drop` was silent in the first dev batch, from a vacuous comparison in the check code (above, README A1); fixed
  before the run of record, no criterion changed.
- The README's "no floating point decision" is not exact: the inertia in `same_sign_form` is a 60-digit eigenvalue
  computation on rational matrices (README A1 (1)).
- The closure argument needs b ∈ C³ (C ∈ C²); #158 assumed C² for its stages. Nothing in #158, #78 or #80 disagrees.

## Recommended tasks

- **Independent check `ind166.py`** (coordinator). Lemma Z, Branch A's ODE, Branch B's separation and the closure
  bracket by code that shares nothing with this script. Acceptance: Lemma Z's two ideals and the family identity
  re-derived from the coordinate Ricci tensor of an independent implementation, at least five mutations on target.
  Depends on this note.
- **Time-dependent shifts with Z ≢ 0.** The finite-type system of #158 and #165 with t-dependence: jets (t, x, ρ), the
  same Lemma R equations (five equations for one function of three variables). Why: it is the remaining time-dependent
  class, and the static class is excluded only conditionally (#165, #177, #178). Acceptance: either the prolongation to
  an empty set (a certificate, or modular evidence stated as such), or an explicit family that passes Lemma R at exact
  points and then closure; planted: members T1 to T4 of this note (Z = 0) must pass the same jet equations, and a
  control that flips the sign of the c ODE must reject them. Depends on #177 or #178 for method, not result.
- **Lean rows** for Lemmas Z, Sep, Ell, ODE and K of the section above (#135).

## Coordinator's review

Reviewed 2026-11-26 at 1df365a (raw `I166-rpi5` 6ee0ec3); the review's files are in `negative-mass-raw/I166-mac` (raw
63a0d43), with its criteria committed before the first run of `ind166.py` (raw 65308b0).

| item | result |
|---|---|
| scope | the three files the task named, and the raw folder |
| rerun on the Mac | exit 0, 14 of 14 checks, 21 controls all fired, 106 s; the output is identical to `one_sign_time_output.txt` apart from timings and the Python patch version |
| independent check `ind166.py` | 7 of 7 criteria pass, 5 s; all 7 mutations make exactly their targets fail |
| hand proofs | Lemma Z and its inequalities, the density of the branches, branch A (the ℓ-coefficient, the ODE and its closed form), branch B's separation and the f-constant case, the gauge and its invariants, the closure bracket and the localization argument were read and rederived; no gap found |
| amendment A1 | no criterion was weakened. The dropped "degree 1" test was extra to the README. The README's "no floating point decision" is not met by `same_sign_form`'s 60-digit inertia, as A1 says; ind166's J2 decides the inertia exactly for T1, T2, T4 and the static member |

`ind166.py` has its own Ricci tensor from Christoffel symbols, imports nothing from this repository, and shares with #166
only sympy, the metric, and #158's block formulas and Lemma R's equations (in its Lemma Z and branch B criteria).
- **Planted case:** b = x gives the 2D curvature-1 Ricci tensor.
- **Family identity:** R_ab = λ(g_ab − ρ²dφ_adφ_b) − α_aα_b holds for undefined c, C, d on the ODE, and fails with the
  ODE's sign flipped.
- **Stage at exact points:** at three exact points each, the characteristic polynomial of R − λg is X²(X² − τX + e)
  with τ < 0 < e for T1, T2, T4 (shifted to t = 1) and the static member. For F1 the 3 × 3 minors have no common root.
- **Lemma Z and the bracket:** Lemma Z's two elimination ideals, and the kernel and bracket of closure, match.
- **Re-solved:** branch B's separation by `solve`, and the c-ODE by `dsolve`.

Findings:
- **The answer stands as stated.** The issue's second alternative is false: c(t) = b_x and C(t) = ρb_ρ are gauge
  invariants, so a member with C′ ≢ 0 is gauge-equivalent to no static one. Closure fails for every member, so
  within Z ≡ 0 only the planar shifts with constant K pass both stages, and none is localized.
- **The silent control** `gauge-drop` exposed a vacuous comparison in the check code. It was fixed before the run of
  record and is reported as silent.

Follow-ups: #179 (time-dependent shifts with Z ≢ 0, the remaining time-dependent class); #135 row 32 (Z, Sep, Ell,
ODE, K).
