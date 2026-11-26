# The static one-sign class with u_xρ ≠ 0: prolongation of #158's finite-type system

November 2026, issue [#165](https://github.com/nmoinvaz/negative-mass/issues/165). Script
[`one_sign_static.py`](one_sign_static.py), its stdout [`one_sign_static_output.txt`](one_sign_static_output.txt). It
builds on [`one_sign_classify.md`](one_sign_classify.md) (#158: Lemma R, §General static shifts with u_xρ ≠ 0). The
script imports `static_blocks`, `v_block` and `rad_member` from `one_sign_classify.py`, and so `phantom_source.py`,
`same_sign_scalars.py` and `two_scalars.py` through it. The coordinates (u, b, v, s), the chain, the planted cases, the
toy system and the verification of the modular certificates are new here. Raw output, Singular inputs and outputs,
certificates and every log are in `negative-mass-raw/I165-rpi5` (its README holds the criteria, the development probes
and the amendments).

Status words as in AGENTS.md, plus **modular evidence** (the coordinator's word for this task): an identity Σ cᵢ gᵢ = κ
≠ 0 checked exactly in GF(p) for the stated primes, by code that does not call Singular, on generators that this script
computed exactly over ℚ. It is a statement about the reduction modulo p at one value of s, not about ℚ, and it is never
called certified.

## Answer

**Outcome 3: not finished.** The chain was computed exactly through level 2. At five fixed values of s (2, 3, 5, 7, 11)
the saturated level-2 ideal is the unit ideal modulo each of the primes 2147483629 and 2147483587: **modular evidence**
(10 identities checked in GF(p)) that no jet with s ∉ {0, ±1} survives for generic s, and no more than that. **Not
met:** an exact certificate over ℚ at any fixed s (D1), and one for all s (D2). Every full-ring Gröbner computation
stopped at its time limit, and the rational reconstruction of the cofactors at one slice failed with up to 50 primes.
There is no surviving set and no explicit solution, and no solution family of any dimension was found.

| item | result | status |
|---|---|---|
| system in (u, b, v, s), ρ = 1 | [X̄₁, X̄₂] + X̄₁ = C̄ (∂_s − (v w_s/w) ∂_v); C̄ = −ρ²C of #158 at exact points | **checked symbolically** |
| f₀, f₁, f₂, f₁₁, f₁₂ | 319, 1388, 1818, 3837, 4892 terms; denominators 8 s u Δ′, 8 s u (1+s²) Δ′ | **checked symbolically** |
| planted u = c x + C ln ρ + d (s = ±1) | every unsaturated F_k ∈ (s²−1); removed only by saturation | **checked symbolically**, **proved** (induction) |
| toy system with C = q²+1 | rejected, cofactors (1, 0, −q) over ℚ | **checked symbolically** (exact over ℚ) |
| level 0, 1 at s = 2, 3, 5, 7, 11 | V(f₀) a surface, V₁ saturated 0-dimensional in each slice (a curve in 4-space) | **modular evidence** (p = 2147483629) |
| level 2 at those s, primes 2147483629, 2147483587 | saturated ideal is the unit ideal | **modular evidence** (10 identities in GF(p)) |
| exact certificate over ℚ at a fixed s (D1), for all s (D2) | not obtained | **NOT MET** |
| (H) dim V₁sat ≤ 1 ⟹ no solution (Lemma Q, Prop. C) | implication proved; (H) has modular evidence only | **conditional** |
| s ≡ 0 on an open set | empty | **proved** (by hand) |
| s ≡ s* constant, s* ∉ {0, ±1} | empty (k′ = d′ = 0 forces u_x = 0) | **proved**, **checked symbolically** |
| surface inside Δ′ = 0 | not treated | **not classified** |

## Setup

The class, from #158. A static shift u(x, ρ) of ds² = −dt² + (dx + b dt)² + dρ² + ρ² dφ² with u_xρ ≠ 0, on a box with u
u_x u_ρ ≠ 0, passes the one-sign stage iff the block R₃ equals λη₃ − v vᵀ (Lemma R). With v₀ = β cos θ, v₁ = β sin θ, v₂
= a/(u cos θ − sin θ) (a = u_x, β = u_ρ, s = tan(θ/2)) the five equations reduce to

    u_xx = Φ₁ = (v₀² + v₂² − a²)/u,   u_xρ = Ψ = −2 v₀ v₂,   u_ρρ = Φ₃ = −β/ρ + 2 v₀ v₁,

and the equality of the mixed derivatives of a and β gives s_x = σ_x, s_ρ = σ_ρ (#158's E₁, E₂), rational in (ρ, u, a,
β, s). A C⁴ solution is an integral surface of X₁ = D_x, X₂ = D_ρ on (ρ, u, a, β, s).

**Scaling (proved).** (x, ρ) → k (x, ρ) maps solutions to solutions and (ρ, a, β) → (kρ, a/k, β/k) on the jets. The jets
(u, A = ρa, B = ρβ, s) are invariant, and X̄ᵢ = ρ Xᵢ act on them. For a weighted-homogeneous f of weight W the ρ = 1
restriction of X₂ f is W f + (a ∂_a + β ∂_β) f + …, and W f lies in the ideal of f, so the chain needs no weight term.
All ideals below are in ℚ[u, b, v, s] with b = B, v = V = A/w, w = u c − σ = D/(1+s²), c = (1−s²)/(1+s²), σ = 2s/(1+s²),
D = u (1−s²) − 2s. In these coordinates

    ρ²Φ₁ = (b² c² + v² − v² w²)/u,   ρ²Ψ = −2 b c v,   ρ²Φ₃ + b = 2 b² c σ,   σ̄_x = ρσ_x,  σ̄_ρ = ρσ_ρ.

The fields X̄₁ = a ∂_u + … and X̄₂ = b ∂_u + … have s components σ̄_x, σ̄_ρ, and (`fields`, **checked symbolically**)

    [X̄₁, X̄₂] + X̄₁ = C̄ (∂_s − (v w_s/w) ∂_v),    C̄ = −ρ² C,   C = D_ρσ_x − D_xσ_ρ  (#158's C).

The sign: #158 writes [X₁, X₂] = C ∂_s, but the s component of [X₁, X₂] is −C. The zero set is the same. `Cpoints`
compares C̄ with −ρ²C computed in 5D with explicit ρ at three exact points (random u, b, v, ρ, s ∈ {1/2, 2, 3}).

**Nonvanishing quantities** (the open set of the reduction): u (u ≡ 0 on an open set is flat); b = ρβ (h > 0, from Lemma
R's dichotomy); v (Z ≠ 0, so v₀ v₂ ≠ 0); s²−1 (u_xρ ≠ 0, given b v ≠ 0); s²+1; D (a ≠ 0, since A = v D/(1+s²)); s (Y, Q
≠ 0; the cases s ≡ 0 and s ≡ s* are decided below); Δ′ (below). The denominators in lowest terms (`denoms`, **checked
symbolically**): σ̄_x has 8 s u Δ′, σ̄_ρ has 8 s u (1+s²) Δ′, and the numerator of C̄ is −v b (s²−1) D f₀. Δ′ is a
polynomial of 26 terms and degree 8 in s, quadratic in (b, v) (stored in the raw folder, `chain/delta.txt`). **The
branch where the surface lies inside Δ′ = 0 is not covered**: there E₁, E₂ do not determine σ and the system is a
different one.

## The chain

f₀ is the numerator of C̄ with the nonvanishing factors removed; f₁ = X̄₁f₀, f₂ = X̄₂f₀, f₁₁ = X̄₁f₁, f₁₂ = X̄₂f₁, each
a numerator in lowest terms with the nonvanishing factors removed, primitive over ℤ. I_{k+1} = I_k + (X̄ᵢ f). All are
even in v.

| generator | terms | degrees in (u, b, v, s) |
|---|---|---|
| f₀ | 319 | (6, 4, 4, 18) |
| f₁ | 1388 | (9, 6, 6, 31) |
| f₂ | 1818 | (10, 7, 6, 33) |
| f₁₁ | 3837 | (12, 8, 8, 44) |
| f₁₂ | 4892 | (13, 9, 8, 46) |

(`chain`: recomputed in sympy, with SHA-256 equal to the stored `chain/chain.json`; the chain builds in 73 s.) In the old
jets (u, a, β, s) the numerator of C̄ has 656 terms and degree 37 (development probes 2 and 3 of the raw README, not part
of the script), so the change a = v w, with D = (1+s²) w, shortened it by about a factor 2.

## Planted cases (through the same chain code)

**The family u = c x + C ln ρ + d (s = 1; also s = −1).** At s = ±1 the numerators of σ̄_x and σ̄_ρ vanish, and u, 1+s²,
D, Δ′, s do not vanish there as polynomials. So X̄ᵢ((s²−1) g) = (s²−1) X̄ᵢ g + 2 s σ̄ᵢ g stays in (s²−1) over the local
ring, and every unsaturated chain numerator F_k lies in (s²−1) (**proved** by induction from the checked premises; F₀,
F₁, F₂ checked to vanish at s = ±1, and at five exact jets of the explicit family, where Lemma R's block is checked
directly). The saturated f₀, f₁, f₂ do not vanish at s = 1: the saturation by s²−1 is what removes the family
(`planted_s1`, **checked symbolically**).

**A toy system with a known inconsistency** (`toy`): X₁ = ∂_p, X₂ = ∂_q + p (q²+1) ∂_r on (p, q, r) has C = q²+1. The
chain gives f₀ = q²+1, f₁ = 0, f₂ = q, and Singular's `liftstd` returns the cofactors (1, 0, −q): (q²+1) − q·q = 1. The
check recomputes the generators and the sum exactly over ℚ (**checked symbolically**). The flat toy X₂ = ∂_q has an
empty chain, so it survives.

## What Singular showed, and what it did not

Modulo a prime p and at a fixed integer s = s₀ the generators are polynomials in (u, b, v), and `std` or `liftstd` takes
seconds. In ℚ[t, u, b, v] with the saturation generator 1 − t u v b D(s₀) Δ′(s₀) (all five factors are nonvanishing
quantities):

| ideal in the slice (mod 2147483629) | dimension in (u, b, v) | degree (`mult`) |
|---|---|---|
| (f₀) | 2 | 10 |
| (f₀), saturated by u v b D Δ′ | 2 | 65 |
| (f₀, f₁, f₂) | 1 | 22 |
| (f₀, f₁, f₂), saturated | 0 | 194 |
| (f₀, f₁, f₂, f₁₁), saturated | empty: unit ideal | |
| (f₀, f₁, f₂, f₁₂), saturated | empty: unit ideal | |
| (f₀, f₁, f₂, f₁₁), saturated by u only | 1 | 26 |
| (f₀, f₁, f₂, f₁₁), saturated by u v b D (not Δ′) | 1 | 32 |

All rows are the same at s₀ = 2, 3, 5, 7, 11 (`--compute=dims`, `certs/dims_p2147483629.json`, printed by the default
mode; the dimension of the unsaturated rows is Singular's dimension in (t, u, b, v) minus 1, since t is free there). So
V(f₀) is a surface in each slice, V₁ off the excluded set is a curve of degree 194 in 4-space (194 points per slice,
with multiplicity), and the leftover of level 2 without saturating by Δ′ lies in Δ′ = 0.

`liftstd` of (f₀, f₁, f₂, f₁₁, 1 − t u v b D Δ′) returns G₁ = 1 at each of the 10 pairs (p, s₀), p ∈ {2147483629,
2147483587}, s₀ ∈ {2, 3, 5, 7, 11}, with cofactors of 1568, 658, 496, 234 and 1754 terms and total degree at most 23, in
about 5 s. **Each identity Σ cᵢ gᵢ = 1 was verified in GF(p)[t, u, b, v] by this script** (`slice_modp`, **modular
evidence**), with g recomputed from the exact chain. Smaller sets already suffice at s = 5 and s = 7 (dev probe 10, mod
32003): f₀ f₁ f₂ with one of f₁₁, f₁₂, saturated by u and Δ′, and no subset without f₁ or f₂. The unsaturated ideal is
not the unit ideal: without saturating by Δ′ a curve of degree 12 to 32 survives in each slice, and it lies in Δ′ = 0
and the other excluded sets.

What this is not. (1) It is a statement modulo p at five values of s. A proper ideal over ℚ is proper modulo all but
finitely many primes, so the unit ideal modulo two primes is evidence for the unit ideal over ℚ and not a proof. (2) It
says nothing about the roots of the elimination polynomial e(s) ∈ J ∩ ℚ[s]: at finitely many s* the slice may be
nonunit, and these could be real. *Added by the coordinator:* these roots do not matter for the classification. If s is
not constant on any open subset, the set where it takes one of finitely many values has empty interior, so the jets of a
solution meet slices that are empty; if s is constant on an open subset, the `sconst` argument below excludes it. So
e(s) ≠ 0, that is the unit ideal over ℚ(s), already gives outcome 1; D2, the unit ideal in the full ring, asks for more. (3) The route to an exact certificate failed at every step that was tried (all logged
in the raw folder):

| attempt | result |
|---|---|
| `std(f₀, f₁, f₂)` mod 32003, 4 variables | no result at 25 min, 1.9 GB (probe 5) |
| the same with s a parameter, F_p(s)-coefficients, with f₁₁ and saturation | no result at 20 min (probe 7) |
| `liftstd` at u = 3 (the other variable) | no result at 10 min (probe 9) |
| `liftstd` over ℚ at s = 7 | no result at 15 min, 1.0 GB (probe 11) |
| 12 to 50 primes at s = 2, CRT with rational reconstruction | fails at every k up to 50 (467 digits of modulus); the "successes" have height near half the modulus, so the coefficients exceed about 230 digits or the lifts differ between primes (probe 13) |
| `std` of f₀ f₁ f₂ f₁₁ and 1 − t u v b s (s²−1)(s²+1) D Δ′ mod 2147483629, block order, full ring | stopped at the 60 min limit: exit 124 at 3600.3 s, 1.8 GB, no result (`certs/fullring.json`, `logs/compute_fullring_1.*`) |

The lift mod p is stable across s: the supports of the cofactors are the same at 60 consecutive values of s (cofactors
2, 3, 4) or differ in at most 3 terms (1 and 5), so the cofactors are rational functions of s with fixed support, but
their degree in s exceeds 20 (a reconstruction from 40 values did not predict held-out values), and each coefficient is
big. An interpolation of the full certificate would need (2d + 2) lifts per prime at about 5 s and many primes, with d
the s-degree.

## Why a surviving set is unlikely, and what (H) would give

At each tested s the saturated ideal of f₀, f₁, f₂ is zero-dimensional (0-dimensional, 194 points per slice, table
above), so V₁ off the excluded set is a curve in 4-space for those s. A solution surface would project to a point or a
curve of this V₁ under φ: (x, ρ) ↦ (u, ρu_x, ρu_ρ, s). A point is impossible (u_x ≠ 0). A curve is covered by the
following.

**Lemma Q (proved).** Let u ∈ C² on a connected open Ω ⊂ {ρ > 0}, u_x ≠ 0, u_ρ ≠ 0, and suppose (u, ρu_x, ρu_ρ) has rank
≤ 1 at every point of Ω. Then on a smaller open set u = F((x − x₀)/ρ) with F′ ≠ 0.

*Proof.* du ≠ 0, so the rank is 1 and A = ρu_x, B = ρu_ρ are C¹ functions Ā(u), B̄(u) locally. The level curves of u
have slope dx/dρ = −u_ρ/u_x = −B̄/Ā = −m(u), constant along each curve, so they are lines x + m(c) ρ = k(c).
Differentiating x + m(u) ρ − k(u) = 0 in ρ gives ρu_ρ = −m ρ/(m′ρ − k′), which must equal B̄(u), constant along the
line, so m k′ = 0, and m ≠ 0 (B̄ ≠ 0). So k = x₀ is constant, all the lines pass through (x₀, 0), and (x − x₀)/ρ =
−m(u), a function of u with nonzero derivative (otherwise u would be constant). ∎

With the translation x → x − x₀ this is #158's conical class u = F(x/ρ), which Proposition C there excludes. Hence:

**Corollary (conditional on H).** If (H): *the saturated ideal (f₀, f₁, f₂) : (u v b s (s²−1)(s²+1) D Δ′)^∞ has
dimension ≤ 1 over ℚ̄*, then no C⁴ solution of the static class with u_xρ ≠ 0 exists on a box where u v b s (s²−1) D Δ′
≠ 0. *Proof.* The image of φ lies in a real algebraic set of dimension ≤ 1, so rank dφ ≤ 1 on a dense open subset
(constant rank theorem), and Lemma Q with #158's Proposition C excludes it. ∎ (H) has modular evidence at five slices
(0-dimensional, 194 points per slice, table above) and is **not proved**; a component of dimension 2 inside one of
finitely many slices s = s* is not excluded by the slices (a surface with s constant is excluded by the `sconst`
argument below). The level-2 computation is stronger than (H). *Added by the coordinator:* for the same reason (H)
can be weakened to (H′): the saturated ideal is zero-dimensional over ℚ(s). Components inside finitely many slices
carry only solutions with s locally constant, which `sconst` excludes.

**s ≡ s* constant is empty (proved, checked symbolically).** If s = tan(θ/2) ≡ s* on an open set, s* ∉ {0, ±1}, put c =
cos θ*, σ = sin θ* (both ≠ 0) and κ = sin 2θ* = 2 c σ ≠ 0. Then v₀ = β c, v₁ = β σ, Q = −v₀v₁ gives u_ρρ + u_ρ/ρ = κ β²,
so B = ρu_ρ has B_ρ = κ B²/ρ and u = −ln(W)/κ + d(x) with W = k(x) − κ ln ρ. Z = −v₀v₂ and v₂ = a/(u c − σ) give u_xρ =
−2 c β a/(u c − σ). The residual times ρ W² is −(k′(c u − σ) + 2 c k′/κ − 2 c W d′). It is linear in the independent
functions ln W, W and 1 of ρ, so the coefficient of ln W, c k′/κ, vanishes, so k′ = 0, and then 2 c d′ = 0, so d′ = 0
and u_x = 0, contradicting u_xρ ≠ 0. The code (`sconst`) checks the solution of the Q equation, the residual and the two
coefficients symbolically; the last step is by hand. ∎

**s ≡ 0 is empty (proved).** If s ≡ 0 on an open set then v₁ = 0, so Q = −v₀v₁ = 0 and Y = −v₁v₂ = 0, that is u_ρρ +
u_ρ/ρ = 0 and u u_xρ/2 + u_x u_ρ = 0. The first gives u = f(x) ln ρ + g(x) with f ≠ 0. The second becomes 3 f f′ ln ρ +
g f′ + 2 f g′ = 0 for all ρ in an interval, so f f′ = 0, so f′ = 0, so f g′ = 0, so g′ = 0, so u_x = 0, contradicting
u_xρ ≠ 0 (Z ≠ 0 forces u_x ≠ 0 on a dense set). ∎

## What the code proves, and what is assumed

Exact identities and ideal computations over ℚ: the block identity and A, B ∈ rad (imported from #158, re-run), the
commutator identity and C̄ ↔ C (`fields`, `Cpoints`), the denominators and the factorization of C̄ (`denoms`), the chain
and its hashes (`chain`), the planted family (`planted_s1`), the toy certificate (`toy`), the s-constant branch
(`sconst`). Modulo p (modular evidence only): the 10 slice identities (`slice_modp`). Everything about the dimensions
and degrees of the Singular ideals (the table above) is evidence read from Singular's output and is not checked by this
script's default mode; it is stored in `certs/dims_p*.json` and printed. Trusted: sympy 1.14.0, Python integer
arithmetic, the imported #158 functions (reviewed in #158). Singular 4.4.1 is trusted for nothing. Assumed: u ∈ C⁴ and
the one-sign stage of #158; the open set where the listed quantities do not vanish; scaling.

## Checks and controls

10 checks run in 155 s on a Raspberry Pi 5 (one process at nice 10, shared with other jobs) and the default mode with
its 11 controls in about 33 min; the script exits 1 if a check fails or a control is not on target. The two `[NOT MET]`
lines (D1, D2) do not change the exit code. From the repository root, in the environment of
[`setup.md`](../../setup.md):

    python3 notes/warp/one_sign_static.py                 # the checks, then every control
    python3 notes/warp/one_sign_static.py --no-controls   # the checks
    python3 notes/warp/one_sign_static.py --mutate=NAME   # one control
    python3 notes/warp/one_sign_static.py --compute=STEP  # Singular: chain | toy | slices | dims | fullring (not run by default)

`--raw` defaults to `../negative-mass-raw/I165-rpi5` beside the repository; the stored certificates are in its `certs/`
and `chain/`.

| check | what it tests | control `--mutate=` | the control's change | failing checks |
|---|---|---|---|---|
| `blocks` | block from (Φ₁, Ψ, Φ₃) equals Lemma R's | `Phi1-sign` | Φ₁ with the opposite sign of a² | `blocks`, `chain`, `denoms`, `planted_s1`, `slice_modp` |
| | | `Psi-sign` | Ψ with the opposite sign | `blocks`, `chain`, `slice_modp` |
| `AB` | A, B in the radical of Lemma R's ideal | `AB-sign` | B built as PQY + Z(Q²+Y²) | `AB` |
| `fields` | [X̄₁, X̄₂] + X̄₁ = C̄ (∂_s − (v w_s/w) ∂_v) | `fields-sx` | the s component of X̄₁ doubled in this check only | `fields` |
| `Cpoints` | C̄ = −ρ²C at 3 exact points, C̄ ≠ 0 | `Cpoints-s1` | the three points drawn on s = 1 | `Cpoints` |
| `denoms` | denominators 8 s u Δ′, 8 s u (1+s²) Δ′; C̄ = c v b (s²−1) D f₀ | `denoms-claim` | 4 s u Δ′ claimed | `denoms` |
| `planted_s1` | the family at s = ±1 survives, saturation removes it | `plant-s2` | the test at s = 2 | `planted_s1` |
| `sconst` | s constant: u = −ln W/κ + d, residual, coefficients | `sconst-kappa` | κ = cσ in place of 2cσ | `sconst` |
| `toy` | the toy with C = q²+1 rejected, the flat toy survives | `toy-flat` | the flat toy passed to the rejection test | `toy` |
| `chain` | term counts, degrees, v even, SHA-256 of f₀ … f₁₂ | the two Φ controls | | |
| `slice_modp` | 10 identities in GF(p), modular evidence | `cert-corrupt` | one cofactor coefficient changed by 1 | `slice_modp` |
| | | `cert-satdrop` | the saturation generator verified with u dropped from P | `slice_modp` |

All controls fire on their targets and on no other check (`one_sign_static_output.txt`). Two controls had mispredicted
targets before the runs of record (the Φ controls: a mutated Φ is solved consistently into σ, so `fields` does not fire,
and `denoms`, `planted_s1`, `slice_modp` do), and `fields` had no control until the raw README's amendment A3. Both are
README amendments (A1, A3). Failed development runs are kept in the raw folder (`logs/dev_*`, `logs/ctl_dev_*`, `dev/`).

## Lemmas stated for Lean (#135)

1. **(Scale)** the jets (u, ρu_x, ρu_ρ, s) are invariant under (x, ρ) → k (x, ρ), and the chain ideals need no weight
   term. Status: proved.
2. **(Comm)** [X̄₁, X̄₂] + X̄₁ = C̄ (∂_s − (v w_s/w) ∂_v) with C̄ = −ρ²(D_ρσ_x − D_xσ_ρ). Status: checked symbolically
   (a rational-function identity).
3. **(Q)** rank ≤ 1 of (u, ρu_x, ρu_ρ) with u_x u_ρ ≠ 0 gives u = F((x − x₀)/ρ). Status: proved.
4. **(S0)** s ≡ 0 gives u_x = 0. Status: proved.
5. **(Sconst)** s ≡ s* ∉ {0, ±1} gives u_x = 0. Status: proved, with the key identity checked symbolically.
6. **(Plant)** the chain ideal of the unsaturated C̄ is contained in (s²−1) over the local ring at s = ±1. Status:
   proved from checked premises.
7. **(H)** as stated above, as a hypothesis, or the unit-ideal statement at level 2 over ℚ(s). Status: not proved.

## What the neighbouring notes should now say

`one_sign_classify.md`, §General static shifts: replace the last bullet's sentence "Further compatibilities (X₁C = X₂C =
0, and so on) were not computed: C is a large rational function, and a Gröbner computation of the chain was not
attempted." by

> Further compatibilities were computed in [`one_sign_static.md`](one_sign_static.md) (#165), in the scale-invariant jets
> (u, ρu_x, ρu_ρ, s): C̄ has 319 terms after removing the nonvanishing factors, X̄₁C̄ and X̄₂C̄ 1388 and 1818, and at five fixed
> values of s the saturated ideal of the second-level chain is the unit ideal modulo two primes (modular evidence). That
> is evidence for emptiness for generic s, not a proof; the statement for all s, the branch Δ′ ≡ 0 and the branch s ≡ s* are
> not decided. A conditional exclusion (H ⟹ none, by Lemma Q and Proposition C) is proved.

and keep the table row "static u_xρ ≠ 0" as **not classified**, adding "(evidence: none; #165)". *Added by the
coordinator:* the quoted text calls s ≡ s* not decided, but this note proves it empty; the text applied to
`one_sign_classify.md` says so. `same_sign_scalars.md`,
"Not covered": keep "static shifts with u_xρ ≠ 0" in the list, now with the pointer to `one_sign_static.md`.

## Discrepancies

- #158 writes [X₁, X₂] = C ∂_s; with its C = D_ρσ_x − D_xσ_ρ the s component is −C (sign only).
- README criteria S3 and S4 as first worded failed in dev run 1 (the missing ∂_v term and the sign), and the control
  targets of `Phi1-sign` and `Psi-sign` were mispredicted. D1 as first worded (exact over ℚ) is **not met**; the modular
  identity is the separate check `slice_modp` (modular evidence). S3, S4 and D1 were changed (README A1, A2) without
  asking the coordinator, who judges S3 and S4 in review. No threshold or point set was changed after a failing run of
  record.
- The step "decide one of three outcomes" ends in outcome 3. The unit-ideal certificate over ℚ(s) that outcome 1 needs
  was not obtained.

## Recommended tasks

- **Exact unit-ideal certificate over ℚ(s) for level 2 (D2), or a proof of (H).** What: Σ cᵢ gᵢ = κ ≠ 0 in ℚ[t, u, b, v,
  s] for g = (f₀, f₁, f₂, f₁₁, 1 − t u v b s (s²−1)(s²+1) D Δ′), generators in `chain/chain.json`, or the dimension
  statement (H). Why: it removes the *conditional* status and decides outcome 1. Not tried: msolve, Magma or F4 on a
  machine with 32 GB or more; interpolation in s of the modular lifts (5 s per lift and prime, fixed support, s-degree
  above 20) with many primes. Acceptance (can fail): the identity checked exactly in sympy, or an exact point of the
  variety (outcome 2). Depends on: nothing.
- *Added by the coordinator:* not needed. The roots of e(s) do not decide anything (see "What this is not"), and this
  task is not filed.
  **The roots of the elimination polynomial e(s).** What: the s* at which the saturated level-2 slice is not the unit
  ideal: scan s₀ over F_p for p near 2000 (about 0.2 s a slice), find e(s) by interpolation, decide which roots are
  real, not in {0, ±1}, and carry real points. Why: these are the only s at which a solution could sit, and the modular
  slices do not see them. Acceptance: the list of roots with, for each real one, the real points of its slice or "none".
  Depends on: nothing.
- **The branch Δ′ ≡ 0.** What: surfaces inside Δ′ = 0, where E₁, E₂ give one compatibility g = 0 and the tangency of X̄ᵢ
  to Δ′ = 0 fixes σ_x, σ_ρ; run the chain of that system. Why: it is the one open branch of the reduction. Acceptance: a
  unit-ideal identity (modular evidence, then exact) or a surviving set with an explicit solution. Depends on: nothing.
- **Independent jet check** (coordinator, as `ind158.py`): extend real points of V₁ at a slice to fifth order with a
  different tool and see where the extension fails. Acceptance: the same conclusion from code that shares nothing with
  this script.
- **Lean rows** for Lemma Scale, Q, S0, Sconst and Plant above (#135).

## Coordinator's review

Reviewed 2026-11-26 at f58b0b8 (raw `I165-rpi5` 1b6fed7); the review's files are in `negative-mass-raw/I165-mac` (raw
3e1e24b), with its criteria committed before the first run of `ind165.py` (raw 567c56b, 94fe654).

| item | result |
|---|---|
| scope | the three files the task named, and the raw folder |
| rerun on the Mac | `one_sign_static.py --raw=…/I165-rpi5`: exit 0, 10 of 10 checks, D1 and D2 `[NOT MET]`, 11 of 11 controls on target, 930 s; the output is identical to `one_sign_static_output.txt` apart from the raw path and timings (rerun 2; rerun 1 named a Python path that does not exist on the Mac and did not run) |
| independent check `ind165.py` | 8 of 8 criteria pass, 26 s; all 7 mutations make exactly their targets fail |
| hand proofs | Scaling, Lemma Q, the corollary from (H), s ≡ 0, s ≡ s* and the planted induction were read and rederived; no gap found |
| amendments A1 to A3 | S3 and S4 as amended are right: [ρX₁, ρX₂] = ρ²[X₁, X₂] − X̄₁, and ∂_s at fixed (u, ρu_x, ρu_ρ) is ∂_s − (v w_s/w)∂_v in (u, b, v, s). D1 is reported not met, and the modular identity is the separate check `slice_modp` |

`ind165.py` derives the system by another route. It works in the invariant jets from X̄ᵢ = ρXᵢ and the chain rule. It
solves the A and B components of [X̄₁, X̄₂] + X̄₁ = 0 for σ̄ and takes C̄ as the s component, where #165 solves E₁, E₂ in
5D and sets ρ = 1. It imports nothing from this repository. It shares with #165 sympy, #158's five equations in Lemma R's
parametrization, and the stored objects under test. Its findings:
- The denominators of σ̄_x, σ̄_ρ are 8 s u Δ′ and 8 s u (1+s²) Δ′, with the stored Δ′, which is irreducible over ℚ.
- C̄'s numerator, stripped of the nonvanishing factors, is −f₀.
- X̄₁f₀, X̄₂f₀, X̄₁f₁ and X̄₂f₁, stripped, are f₁/4, f₂/4, f₁₁/8 and f₁₂/8 exactly.
- Two stored certificates, re-expanded over GF(p), sum to 1 (modular evidence).
- `dsolve` re-solves the s ≡ s* and s ≡ 0 cases, and the planted family has σ̄ = C̄ = 0 at 12 exact jets.

Findings:
- **The answer stands:** outcome 3, with an exclusion *conditional* on (H) outside the branch Δ′ ≡ 0. The level-2 unit
  ideal and (H) are supported by *modular evidence* at five slices.
- **Less is needed than the note asked for** (the two sentences marked above). A solution with s not locally constant
  meets slices that are empty, and `sconst` excludes s locally constant. So outcome 1 needs only the unit ideal over
  ℚ(s), not over the full ring (D2), and (H) weakens to zero-dimensionality over ℚ(s) (H′). The roots of e(s) decide
  nothing, so that recommended task is not filed.
- **#158's sign slip** ([X₁, X₂] = C∂_s where the s component is −C) is reported in `one_sign_classify.md`.
- **Process.** The development probes ran for 2 h before the criteria were committed; the README records them.
  - S3, S4 and D1 were changed after development failures without asking the coordinator. S3 and S4 corrected
    identities that were stated wrong.
  - D1 at first passed a modular identity under its own name, with the status "certified over GF(p)". After the
    coordinator's message it was reported as not met (A2).
  - The suggested text for `one_sign_classify.md` called s ≡ s* undecided although this note proves it empty. The text
    applied there says it is decided.

Follow-ups: #177 (the unit ideal over ℚ(s), or (H′)); #178 (the branch Δ′ ≡ 0); #135 row 31 (Scale, Q, S0, Sconst,
Plant).
