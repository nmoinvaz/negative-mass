# F + r f′ at the roots of (r f′)′ = 0 on the tanh wall, for every σR_b

November 2026, issue [#132](https://github.com/nmoinvaz/negative-mass/issues/132). Script
[`he_tanh_bound.py`](he_tanh_bound.py), its stdout [`he_tanh_bound_output.txt`](he_tanh_bound_output.txt). It builds on
[he_surface.md](he_surface.md) section 4 (#72) and [review_he_surface.md](review_he_surface.md) (#101). Nothing here claims
novelty.

Status words are those of AGENTS.md. **Checked symbolically** means an exact identity in sympy. **Certified** means a
computer-assisted statement by interval arithmetic with the trusted code named below.

## Answer

For every σR_b > 0 and every root r > 0 of (r f′)′ = 0 of the tanh profile f = [tanh σ(r + R_b) − tanh σ(r − R_b)] /
(2 tanh σR_b), with F = 1 − f,

    r f′ < 0     and     F + r f′ < 0 .

**Proved** by hand (section 2), with every identity **checked symbolically**, and **certified** in Lean: the whole
statement, with Mathlib's derivatives, is `tanh_wall_no_typeI` in [lean/NegativeMass/TanhWall.lean](../../lean/NegativeMass/TanhWall.lean)
(#134). The proof is one variable. With w = 2σr and
C = cosh 2σR_b the condition (r f′)′ = 0 is linear in C, and at a root

    F + r f′ = − w sinh w (cosh w − 1) m(w) / [ Q(w) (cosh w + C)² ],
    m(w) = w (cosh w + 2) − 3 sinh w = Σ_{k≥2} 2(k − 1) w^(2k+1) / (2k + 1)!,     Q(w) = sinh w + w cosh w ,

and every factor is positive for w > 0. The statement holds for every real C > −1, so it does not depend on C being
cosh 2σR_b.

An independent **certified** cover by interval arithmetic on the profile's own jets (section 3) confirms it for
σR_b ∈ [0, 100] in 794 cells, with F + r f′ ≤ −0.05 on every root enclosure. The cover reaches σR_b = 0 through
the variable y = σr, where f = sech² y. The largest upper bound over the cells is −0.050001 (σR_b ∈ [0.2289, 0.2319]).
The supremum over σR_b is the limit σR_b → 0, F + r f′ = −0.0536756 at y = 1.0095565 (**checked** by enclosure to
25 digits, section 4). It is not attained for σR_b > 0.

So T2 of he_surface.md has no real solution at a root of T1 on the tanh wall. At a root F > 0 and r f′ < −F, so
r f′ + F < 0 and r f′ + c²F ≤ r f′ + F < 0, and v² (r f′ + F)(r f′ + c²F) + c² > 0 for every c and every v ≠ 0. **The
Type I set of D = 0 off the equator and off the axis is empty on the tanh wall for every σR_b, R_b and v.** It replaces
the numerical evidence at 49 and 97 values of σR_b (#72, #101).

## 1. Variables

With y = σr, p = σR_b, w = 2y and C = cosh 2p ≥ 1, the profile is

    f = cosh² p / [cosh(y + p) cosh(y − p)] = (1 + C) / (cosh w + C) ,     F = (cosh w − 1) / (cosh w + C) > 0 .

The first form has no 0/0 at p = 0 and tends to sech² y there. Both equalities are exact identities in e^y and e^p
(`identity`). Roots r > 0 of (r f′)′ = 0 are roots w > 0, and F and r f′ depend on (w, C) alone, so the problem has two
parameters (w, C) and one equation. With L = ln f and primes in y, f′/f = L′ = −(tanh(y + p) + tanh(y − p)) and

    r f′ = y f L′ = −(1 + C) w sinh w / (cosh w + C)² ,      (r f′)′ = σ f Ψ ,      Ψ = L′ + y (L′² + L″) .

## 2. Route 1: the proof

**The root condition is linear in C.** With P = sinh w cosh w + w (2 − cosh² w) and Q = sinh w + w cosh w,

    Ψ (cosh w + C)² = −2 (P + C Q) ,

so (r f′)′ = 0 for r > 0 exactly where N(w; C) = P + C Q vanishes, w > 0 (`reduction`). Q > 0 for w > 0.

**F + r f′ is linear in C.** With A = cosh² w − cosh w − w sinh w and B = cosh w − 1 − w sinh w,

    (F + r f′) (cosh w + C)² = A + C B .

**Eliminate C at a root.** There C Q = −P, so Q (A + C B) = A Q − B P =: T, and

    T = A Q − P B = − w sinh w (cosh w − 1) m(w) ,    m(w) = w (cosh w + 2) − 3 sinh w

(an exact identity in w and e^w, `reduction`). The coefficient of w^(2k+1) in m is 1/(2k)! − 3/(2k + 1)! + 2 δ_{k0}, which
is 2(k − 1)/(2k + 1)! for k ≥ 1 and 0 at k = 0. It vanishes at k = 0 and 1 and is positive for k ≥ 2, so m(w) > 0 for
w > 0. Hence T < 0 and F + r f′ = T / [Q (cosh w + C)²] < 0 at every root. ∎

The code checks each identity exactly, the coefficient formula for k ≤ 40 in `Fraction` arithmetic and sympy's series
of m to w^21. The positivity for all k is the hand argument above. There is no search and no tolerance in this route.

## 3. Route 2: the cover, independent of the elimination

The cover shares with route 1 only the profile. It evaluates f, L′, L″, L‴ (the jets of ln f, closed forms from `identity`
and `jets`) on boxes in (y, p) by interval arithmetic, and never uses N, P, Q, A, B or T. For a cell of p:

1. **r → 0.** Ψ / y = L′/y + L′² + L″ < 0 on (0, y_a] for all p in the cell, y_a = 1/16, with L′/y = −2 sinhc(2y) / (cosh(y + p)
   cosh(y − p)) enclosed by monotonicity (sinhc(z) = sinh z / z, increasing, equal to 1 at 0).
2. **r → ∞.** For y ≥ y_t = p_max + 2 every p in the cell has y − p ≥ 2, so tanh(y ± p) ≥ tanh 2, L′ ≥ −2,
   L′² ≥ 4 tanh² 2 and L″ ≥ −2 sech² 2. Then Ψ ≥ −2 + y (4 tanh² 2 − 2 sech² 2) > 0, since 4 tanh² 2 − 2 sech² 2 ≥ 3.5760951
   and y_t ≥ 2. There is no root in [y_t, ∞).
3. **In between.** [y_a, y_t] is bisected until each leaf has Ψ > 0, Ψ < 0 or width ≤ max(2⁻¹⁴, Δp / 8). The leaves must
   tile [y_a, y_t], and with the two exclusions tile (0, ∞) with no gap. Their signs must read − (small y), − … −,
   candidates, + … +, + (tail): one group of candidate leaves with a sign change from − to +. Ψ is continuous, so every p in
   the cell has a root in the group. Ψ_y > 0 on every candidate leaf (Ψ_y = 2L″ + L′² + y (2 L′ L″ + L‴)), so the root is
   unique for each p.
4. **The bound.** On every candidate leaf, F + r f′ = 1 − f (1 − y L′) ≤ −0.05 and r f′ = y f L′ < 0. A cell that fails any
   step is split at its midpoint and retried, down to a width of 10⁻⁶, where the cover reports failure.

The 794 cells tile [0, 100] exactly. The table below samples them (all cells are in `cells_full.tsv` of the raw
folder `I132-rpi5`).

| σR_b cell | root y enclosure | upper bound of F + r f′ | candidate leaves |
|---|---|---|---|
| [0, 0.00153] | [1.00740, 1.01178] | −0.05006 | 37 |
| [0.0992, 0.1007] | [1.01169, 1.01604] | −0.05071 | 35 |
| [0.2289, 0.2319] | [1.02805, 1.03679] | −0.050001, the largest | 33 |
| [0.4944, 0.5005] | [1.10531, 1.12317] | −0.05731 | 30 |
| [0.9766, 1.0010] | [1.35956, 1.42843] | −0.07270 | 24 |
| [1.9531, 2.0508] | [2.13454, 2.35265] | −0.13819 | 28 |
| [4.8828, 5.0781] | [4.92685, 5.29681] | −0.72232 | 27 |
| [9.9609, 10.156] | [9.93591, 10.31384] | −2.21109 | 16 |
| [29.883, 30.078] | [29.84266, 30.15531] | −8.48579 | 20 |
| [99.805, 100] | [99.76016, 100.05881] | −30.28088 | 24 |

The cells are 266 in [0, 1), 32 in [1, 3), 36 in [3, 10) and 460 in [10, 100), with widths 0.00153 to 0.195. A cell is split
only until its bound reaches −0.05, so the largest upper bound sits just below that criterion. The true supremum, −0.0536756
(section 4), is 0.0037 below it.

The closed form of route 1 is not used to prove anything here. It checks the cover: at the left end, the middle and the right
end of every cell, the root w of N(w; cosh 2p) = 0, enclosed by certified bisection, and F + r f′ there both lie inside the
cell's enclosures. The scan of #72 and #101 agrees as well: −0.05368 at σR_b = 0.01.

The same `cell_proof` runs on the planted profile f = sech²(r)/10 of #101. Its root is enclosed in r ∈ [1.009491, 1.009609]
with F + r f′ ∈ [0.894622, 0.894642] > 0 and r f′ ≤ −0.064 < 0, so T1 and T2 have real solutions there and the code finds the
set non-empty. It is the p = 0 profile (C = 1) scaled by 1/10, so its root equals the root of the limit of the tanh wall.

## 4. The ends

**σR_b → 0.** C → 1 and f → sech² y. The root of N(w; 1) is w₀ = 2.0191131, y₀ = 1.0095565 and F + r f′ = −0.0536756 there,
enclosed to 25 digits by `sup`. This is the supremum of F + r f′ over the tanh walls, approached as σR_b → 0 (the scan's
−0.05368 at σR_b = 0.01). The cover includes the cell p ∈ [0, 0.00153].

**σR_b → ∞.** Route 1 needs no asymptotics. For the asymptote he_surface.md quotes, `large_p` encloses the root at
σR_b = 100, 10³ and 10⁴ and finds F + r f′ in (−p/2 + 1/2, −p/2 + 1/2 + 1/(4p)):

| σR_b | y at the root | F + r f′ | F + r f′ + p/2 − 1/2 | 1/(8p) |
|---|---|---|---|---|
| 100 | 100.00499979168 | −49.4987500833 | 0.00125 | 0.00125 |
| 1000 | 1000.0004999998 | −499.499875 | 0.000125 | 0.000125 |
| 10⁴ | 10000.00005 | −4999.4999875 | 1.25 × 10⁻⁵ | 1.25 × 10⁻⁵ |

At leading order F + r f′ = −(w − 3)(w + 1) / (4 w) with w = 2p + 1/p + …, which is −p/2 + 1/2 + 1/(8p). This is
**numerical evidence** at three values of p. It is a remark on the size, and no conclusion above depends on it.

## 5. Checks and controls

All 7 checks pass and all 10 controls fire exactly on their targets (`he_tanh_bound_output.txt`). Planted truth
runs before the cover. Each target list was written to the raw README (`I132-rpi5`) before the first control run and is
unchanged.

| check | what it tests | control `--mutate=` |
|---|---|---|
| `identity` | 13 exact identities: the three forms of f, L′, L″, L‴, Ψ, Ψ_y, E, the planted jets | `identity-break` |
| `jets` | the interval evaluators against mpmath's `diff` of the definition at 7 points, p = 0 included | `identity-break`, `jets-typo` |
| `reduction` | the exact identities of section 2 and the Taylor coefficients of m | `proof-typo`, `factor-typo` |
| `planted` | f = sech²(r)/10 found non-empty by `cell_proof` | `planted-eps`, `no-tail`, `no-small` |
| `cover` | the cells of section 3 and the cross-check with the closed form | `shift`, `identity-break`, `no-tail`, `no-small`, `cell-gap` |
| `large_p` | the closed form at σR_b = 100, 10³, 10⁴ | `asym-typo` |
| `sup` | the limit σR_b → 0 | none |

`shift` adds 0.06 to F + r f′ in the cover, so the bound fails at small p, where the true value is −0.0537 + 0.06 > 0.
`no-tail` and `no-small` drop one of the two exclusions, so every cell leaves a gap in (0, ∞) and `cell_proof` refuses it.
`identity-break` puts cosh p for cosh² p in the cosh form of f, which breaks the identity. `sup` has no control of its own.
Its value is cross-checked by the cover cell that contains p = 0.

From the repository root, in the environment of [`setup.md`](../../setup.md):

    python3 notes/warp/he_tanh_bound.py

It exits 1 if a check fails or a control does not fire. `--no-controls` runs the checks alone, `--pilot` the cover of
p ∈ [0, 0.01] only, in 13 s, and `--mutate=NAME` one control. The full run, the 7 checks and then the 10 controls,
took 8150 s (2 h 16 min) on one core of a machine shared with another runner. The checks alone took 950 s there, and
372 s on a quieter core, where the 10 controls took about 44 min.

## What the code proves, and what is assumed

**Proved** by hand: section 2, from the identities of `reduction`. **Checked symbolically**: every identity in
`identity` and `reduction`, in sympy 1.14.0 as rational functions of a = e^y and q = e^p, and the Taylor coefficients of
m. **Certified**: section 3, by mpmath 1.3.0 `iv` for p ∈ [0, 100], and the enclosures of section 4. The statement of the
Answer holds for every σR_b > 0 on the strength of route 1 alone. Route 2 is a second proof on [0, 100]. **Certified** in
Lean (#134): the Answer as `tanh_wall_no_typeI`, from `m_pos` (m > 0 by m(0) = m′(0) = m″(0) = 0 and m‴ = w sinh w),
`sign_at_root` and `profile_cosh`, each depending only on the axioms propext, Classical.choice and Quot.sound; trusted
code: the Lean 4 kernel and Mathlib v4.34.1 (`lake build` in lean/, raw `I134-mac`).

Trusted: mpmath 1.3.0 `iv` (outward-rounded + − × ÷ and `exp` on intervals, from which cosh, sinh and tanh are built here),
mpmath's `diff` in `jets`, sympy 1.14.0 `cancel`, `diff` and `series`, and Python's `Fraction`. The cover shares with route 1
the closed forms of f, L′, L″, L‴, which `identity` checks against the definition and `jets` against numerical derivatives.
The cross-check shares the closed form with route 1.

Assumed: the profile is the tanh wall with R_b, σ > 0, and r > 0. For C = cosh 2σR_b > −1 the denominators are positive. T2,
T3 and the Type I set are those of he_surface.md and are not re-proved here.

Not covered: the accelerating bubble (#73), other wall profiles, and a bound on F + r f′ uniform in σR_b beyond [0, 100]
(route 1 gives the sign, not −0.05). The count of Type III points on the wall is #133.

## For he_surface.md and review_he_surface.md

The coordinator owns those files. Proposed changes:

1. he_surface.md, Answer, third paragraph: "a scan of σR_b from 0.01 to 100 finds F + r f′ ≤ −0.054 there (**numerical
   evidence**)" becomes "F + r f′ < 0 at every root for every σR_b > 0 (**proved**, [he_tanh_bound.md](he_tanh_bound.md))".
   The old "≤ −0.054" is a rounding slip: its own section 4 table says "at most −0.0537" (σR_b = 0.01), and −0.0537 >
   −0.054. The supremum is −0.0536756. The same "≤ −0.054" is in claim 7 of review_he_surface.md.
2. he_surface.md section 4, the paragraph after the table: "The scan is **numerical evidence** … It is not a bound" becomes a
   pointer to he_tanh_bound.md, with the supremum −0.0536756 as σR_b → 0. "At the ends F + r f′ is close to −σR_b/2 + 1/2
   for large σR_b" is **numerical evidence** at σR_b = 100, 10³, 10⁴ (section 4 above), with the next term 1/(8σR_b).
3. he_surface.md "Not covered": delete "the tanh scan is a scan, not a bound (#132)". Section 5, "Off the equator on the
   tanh wall no point of D = 0 is Type I", then rests on a proof.
4. review_he_surface.md, claim 7 and claim 26, and its "Not covered": the tanh scan as a bound and the large-σR_b asymptote
   are answered by this note. Recommended task 1 is done.

## Recommended tasks

1. **Independent review of #132.** What. Re-derive P, Q, A, B and the factorization T = −w sinh w (cosh w − 1) m(w) from the
   profile, and the positivity of m, without reading `he_tanh_bound.py`, and evaluate F + r f′ at the root by an independent
   root finder at 50 digits for σR_b = 10⁻³ to 10³. Why. AGENTS.md asks for an independent review of a computer-assisted proof
   before it is relied on, and this one replaces the scan that the empty Type I set rested on. Acceptance. The identities
   reproduced, the root values agreeing with the table of section 3, and a control that breaks each. Depends on #132.
2. **A uniform bound F + r f′ ≤ −0.05 for every σR_b.** What. Show that the value of F + r f′ at the root, a function of w alone
   on the curve C = −P(w)/Q(w), is at most −0.0536 for w ≥ w₀, by an interval cover in w and a hand bound for w → ∞. Why. The
   cover stops at σR_b = 100 and route 1 gives only the sign. Acceptance. An enclosure below −0.05 on every w-cell and an
   explicit bound beyond, with a control that raises the value by 0.06. Depends on #132.
