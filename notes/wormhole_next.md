# The wormhole questions: what is tractable next

November 2026. Which open questions of [`critical_ghost.md`](critical_ghost.md) (CG), [`bh_to_wormhole.md`](bh_to_wormhole.md)
(BW), [`wormhole_largeD/`](wormhole_largeD/README.md) and Paper III can be moved now, with this repository's tools and those
of the sister project [choptuik-large-d](https://github.com/nmoinvaz/choptuik-large-d) (CLD), and what each costs. Pilot
scripts, logs and records are in [`wormhole_next/`](wormhole_next/).

Status words: **literature** (with a reference), **numerical evidence**, **checked symbolically** (sympy, with a control
that fails), **argued** (a derivation, not a proof), **conjecture**, *(unverified)* for a reference detail not confirmed.

## Summary

- **The "untested escape" of BW §3 and §6.2 was tested in CLD after BW was written.** BW (2026-10-24) says the
  ghost-sector rates μ of Choptuik's solution are known only at D = 4 and 5. CLD's
  [`notes/GHOST_ADMIXTURE.md`](https://github.com/nmoinvaz/choptuik-large-d/blob/main/notes/GHOST_ADMIXTURE.md)
  (2026-10-25 to 27) counts the zeros of the test field's determinant by the argument principle. It finds only the boost
  μ = 0 at D = 8.25, 16.76, 30.07 and 45.05 for −0.1 ≲ Re μ ≲ 6.2, and near the imaginary axis at D = 4, each in two node
  layouts (**numerical evidence**). BW and CG are out of date here. §5 drafts the replacement text.
- **Pilot (§1.3): the least-damped real ghost-sector rate at 12 values of D from 4 to 67**, on CLD's stored
  double-double states with CLD's test-field determinant. Every value is negative, and μ falls steadily with D (table in
  §1.3). Resolution and layout changes move it by less than 10⁻¹¹. D = 4 reproduces CLD's −0.3651202372 and lies within
  Gundlach and Martín-García's −0.362 ± 0.012 [GMG96]. A planted zero is recovered and a wrong operator fails as it
  should.
- **One pre-stated check fails. At D = 5 the rate is −0.44304, outside the window −0.450 to −0.446 of Mittal, Gundlach and
  Martel [MGM26]** (a periodicity fit to evolutions, with no error bar). This is reported, not resolved (§1.4).
- **At D = ∞ no symmetry fixes μ (argued, §1.5), but the data approach −1.** Fits give μ∞ between −1.06 and −0.995 (a
  spread over fit choices), and μ + 1 ≈ 45/n² over n = 28–65. That μ∞ = −1, the rate of the flat-space mode ψ = t, is a
  **conjecture**. Then BW's fade-out exponent 2|μ|γ would tend to 1 (0.27 at D = 4, 0.95 at D = 45).
- **Recommended first step** (§4): finish the ghost-sector question on CLD's existing code and states. That means
  winding counts at D = 5–7 and 50–78, a box between the real zero and the origin at a few D, the line Im μ = ω at every
  pilot D, and a planted-truth extrapolation of μ∞. It needs about 10 h of CPU at 2–4 threads and no new code. It turns
  BW's argument into numerical evidence for every D in [4, 78]. Second (§2.1): a proof that the symmetric wormhole has
  exactly one unstable mode in every D, which Paper III's type I reading needs.

## 1. The ghost sector of Choptuik's solution

### 1.1 The problem (argued in CG §2.3, restated)

A second massless scalar ψ vanishes in the critical solution and enters the stress tensor quadratically, ghost or
normal. At linear order it does not act on the metric or on φ. It obeys □ψ = 0 on the fixed background, the same
equation for either sign. Its perturbations are ψ ∝ e^{μτ} × (periodic in τ), with τ = −ln(−t), regular at the centre
and at the past light cone. They are the imaginary-part perturbations of a complex scalar about the real solution, so the
ghost spectrum is the charge-mode spectrum of [GMG96, MGM26]. ψ = φ_Ch is an exact zero: the boost here, the U(1) phase
for the complex field.

In CLD's variables the test field is the p, m pair of the scalar's equations, with a and c held at the background
(GHOST_ADMIXTURE §1):

```math
x(x-c)\,p_x = -x\,(\partial_\tau + \mu)\,p + G p - c p + \tfrac{n}{2}\,c\,(p+m),\qquad
x(x+c)\,m_x = -x\,(\partial_\tau + \mu)\,m - G m + c m - \tfrac{n}{2}\,c\,(p+m),
```

with n = D − 2 and G = (n − 1)(c a² − c) taken from the background. CLD's crate `hp/modecount_cx` builds the
multiple-shooting matrix of exactly this system (`CSolver.test_det`, mode 1) and removes a's and c's tangents during the
integration. Its mode 2 (`TEST_MUT`) keeps a's response to p and m. That is a wrong operator, used here as the control.

So, to the task's question: the test field decouples, CLD's code already computes its modes, and double-double states at
M = 32–64 suffice. Real-axis zeros converge to 10⁻¹¹ (§1.3). One determinant evaluation takes 5–90 s at 2 threads on
this machine, depending on D and M.

### 1.2 What CLD has (literature: GHOST_ADMIXTURE §§2–4, CLD README, "A second scalar field adds no unstable mode")

| D | M | region counted | zeros of det M_ψ | layouts |
|---|---|---|---|---|
| 4 | 48, 64 | −0.114 < Re μ < 0.193, \|Im μ\| < 6.38 | 0 only (the 9-node layout's growing zeros are artifacts) | 2 |
| 4 | 64 | the lines Im μ = 0 and ω, −0.7 < Re μ < −0.02 | −0.3651202 on the real axis, none on Im μ = ω | 2 |
| 8.25 | 32 | −0.08 < Re μ < 6.70, \|Im μ\| < 1.5ω | 0 only | 2 |
| 16.76 | 32 | −0.06 < Re μ < 6.26, \|Im μ\| < 1.5ω | 0 only | 2 |
| 30.07 | 32 | −0.15 < Re μ < 6.18, \|Im μ\| < 1.5ω | 0 only | 2 |
| 45.05 | 32 | −0.13 < Re μ < 6.28, \|Im μ\| < 1.5ω | 0 only | 2 |

The determinant's zeros repeat at μ + 2ikω, with ω = 2π/Δ, in CLD's folded layout. So the strip |Im μ| ≤ ω holds every
Floquet exponent once, and the line Im μ = ω is the other half-period parity sector. Each count has a planted control.
Each also passes CLD's port and identity checks, each with a control that fires (GHOST_ADMIXTURE §2).

For BW, this settles that no ghost-sector mode grows inside those boxes, at five values of D from 4 to 45. It leaves
three things open:
- D between 4 and 8.25 and above 45, where nothing was counted;
- the region left of Re μ ≈ −0.1, off the two lines, which decides whether the real zero is the least-damped mode;
- the value of μ, which BW's exponent needs.

### 1.3 Pilot: the real-axis zero at D = 4–67

**Method.** `wormhole_next/ghost_mu.py` imports CLD's `hp/modecount_common.py` read-only and evaluates `test_det` on a
grid of real μ. It refines each sign change by the Illinois method, on the determinant with its exponential trend divided
out (det M_ψ grows like e^{~300μ}). The real axis carries no candidate poles: the horizon poles sit at odd multiples of
iω (CLD `hp/modecount_common.py`), so every sign change is a zero.
- **States.** From CLD's `hp/data` and from the raw repository's `T79-win/` (choptuik-large-d-raw at 3137102).
- **Environment.** CLD's `.venv` (hpkernel of 2026-10-28) with `hp/modecount_cx/build/modecount_cx.so` of 2026-10-25.
  That `.so` is byte-identical to the build CLD's T66 and T70 ghost runs used.
- **Load.** Two threads per job at nice 10, at most two jobs at once. About an hour of wall time for everything.
- **Output.** Logs in `wormhole_next/logs/`, records in `wormhole_next/ghost_mu.jsonl`. The table and checks come from
  `wormhole_next/ghost_mu_table.py`, with its output in `ghost_mu_table.log`. The D = 4 log has no closing "exit" line
  because its queue's wrapper was stopped early to move the controls to D = 5. Its record is complete.

**Pass criteria, written before the runs:**

1. D = 4 reproduces CLD's −0.3651202372 to 10⁻⁸. This checks this environment against CLD's own runs.
2. D = 4 lies within [GMG96]'s −0.362 ± 0.012 (perturbation theory, an independent code).
3. D = 5 lies within [MGM26]'s window −0.450 to −0.446 (evolutions, an independent method).
4. M = 48 and 64 agree to 10⁻⁶ at D = 5 and 6, as do M = 32 and 48 at D = 8.25 and a second node layout at D = 5.
5. A zero planted at −0.17 is returned to 10⁻⁹. It is planted by multiplying det by μ + 0.17, with the same evaluations
   and the same finder.
6. Control. With mode 2, the smallest pivot at μ = 0 must be more than 10³ times mode 1's, since the boost zero is lost.
   No real zero may lie within 0.012 of mode 1's, so that criterion 2 or 3 would fail with the wrong operator.

**Results.** Criteria 1, 2, 4, 5 and 6 pass and criterion 3 fails. The script exits 1 for that reason.
- Criterion 4: the zero moves by at most 2×10⁻¹² between M = 48 and 64. It does not move at all, at the printed 12
  digits, between the 10-node and 11-node layouts at D = 5 (the 11-node layout uses kc 80 and ks 200 and is the one that
  removed D = 4's artifacts). The second layout was made by CLD's `hp/relayout.py`, with Δ unchanged to the last digit.
- Criterion 5: the planted zero comes back at −0.1699999999998.
- Criterion 6: mode 2's pivot at μ = 0 is 2.2×10⁻³ against mode 1's 2.7×10⁻²⁶. Its real zeros are −0.4024 and −0.8707,
  and the nearer one is 0.041 from mode 1's.

| D | μ, least-damped real zero | other real zeros in the window | states (M) | 2\|μ\|γ | δ_Q = (D − 2 − μ)γ |
|---|---|---|---|---|---|
| 4 | −0.3651202372 | (window −0.40 to −0.30) | 64, 11 nodes | 0.2731 | 0.8845 |
| 5 | −0.4430361582 | none above −1.0 | 48, 64; 2 layouts | 0.3661 | 1.4224 |
| 6 | −0.5052240827 | none above −1.0 | 48, 64 | 0.4346 | 1.9376 |
| 7 | −0.5573230924 | none above −1.2 | 48 | 0.4904 | 2.4452 |
| 8.25 | −0.6124708660 | −1.5110 | 32, 48 | 0.5489 | 3.0750 |
| 10 | −0.6761242034 | none above −1.4 | 48 | 0.6161 | 3.9527 |
| 12 | −0.7345316442 | −1.3994 | 48 | 0.6778 | 4.9530 |
| 16.76 | −0.8312825960 | −1.3743 | 32 | 0.7814 | 7.3274 |
| 20.04 | −0.8743130935 | −1.3813 | 32 | 0.8284 | 8.9625 |
| 30.07 | −0.9423919782 | −1.4401 | 32 | 0.9063 | 13.9494 |
| 45.05 | −0.9755511996 | −1.5354 | 32 | 0.9486 | 21.4059 |
| 67.01 | −0.9895615865 | none above −1.4 | 32 | (no γ stored) | |

- Every window runs to −0.04, except D = 4's (CLD's scan already covers −0.7 to −0.02 there) and the M = 64 and second
  layout checks.
- γ is CLD's: `hp/data/local/gamma_dd.jsonl` for D ≥ 8.25 and `T79-win/gamma.jsonl` below.
- The exponent 2|μ|γ is BW §3's fade-out exponent. δ_Q is the charge exponent of a charged scalar in D dimensions
  ([MGM26], Eq. (66)). At D = 4, (2 + 0.36512)(0.37396) = 0.8845, against [GMG96]'s 0.883 ± 0.007.
- The other parity sector at D = 5 (CLD's `hp/modecount_test_axis.py`, line Im μ = ω, Re μ from −1.0 to −0.02) has
  one zero, at −0.94997 + iω, and nothing in (−0.92, −0.02). On both lines the least-damped rate at D = 5 is therefore
  the real one.

What the pilot does not cover:
- complex zeros between the real zero and Re μ ≈ −0.1 off the two lines;
- the line Im μ = ω at D ≠ 4, 5;
- winding counts at D = 5–7, 10, 12, 20 and 67.

So "least damped" in the table means least damped on the real axis.

### 1.4 The D = 5 discrepancy

Criterion 3 fails. The rate −0.44304 lies 0.003 beyond the edge of [MGM26]'s window −0.450 to −0.446, a difference of
1.1% of μ.

On this side, four things are excluded:
- the τ resolution (M = 48 and 64 agree to 2×10⁻¹²);
- the node layout and the centre and horizon series orders (second layout, same 12 digits);
- the zero finder (the planted zero comes back to 2×10⁻¹³);
- the operator (the wrong operator moves the zero by 0.04, so the comparison discriminates).

The same code at D = 4 lands within [GMG96]'s perturbative error and reproduces CLD's independent scan.

On MGM26's side, the number is μ chosen to make the scaled imaginary field "most strictly periodic (with almost equal
peaks)" over the self-similar part of an evolution. They give no error bar. Their Δ = 3.224 and γ = 0.41325 differ from
CLD's 3.2217019 and 0.4131269 by 7×10⁻⁴ and 3×10⁻⁴ relative.

A candidate cause (**conjecture**): faster-decaying ghost-sector modes in the early echoes. The other parity sector's
mode at −0.950 + iω flips sign every half period and still looks periodic over a full one, and any further real-axis
mode lies below −1.0. Either kind inflates early peaks, and equalizing the peaks then pulls the fitted μ downward, which
is the direction of the difference.

The consequence for the charge exponent is small: δ_Q(5) = 1.4224 here, against their 1.425 (0.2%). The coordinator
decides whether to raise this with the authors. No public action has been taken.

### 1.5 Large D: is there a Proposition 5 for μ?

CLD's Proposition 5 derives λ∞ = 2 (γ∞ = 1/2) from a symmetry of the D = ∞ inner equations: an s-translation combined
with a time reparametrization. The ghost sector has no such mode (**argued**):

- The inner equations are Y_t = cS, S_t = c(Y_s + Y³/4) and c_s = cY²/4, with s = r²/(2n). A test field obeys
  ψ_s = c⁻¹(c⁻¹ψ_t)_t with c frozen. With Ψ = 2ψ_t/c and Σ = Ψ_t/c this becomes Ψ_t = cΣ and
  Σ_t = c(Ψ_s + Y²Ψ/4) = (cΨ)_s. That is the inner system linearized with δc = 0, with the cubic's coefficient 3/4
  replaced by 1/4 (**checked symbolically**, `wormhole_next/inner_test_field.py`; its control, the coefficient 3/4,
  leaves a residual).
- The two-field inner system has the target's symmetries (translations, and the rotation or boost) and the spacetime ones
  (s- and t-translation, scaling, reparametrization). The spacetime symmetries map ψ = 0 to ψ = 0, so their generators
  have δψ = 0. The target translation gives δψ = const, which is invisible in Ψ. Only the rotation gives a ghost mode,
  δψ = φ with μ = 0.
- Proposition 5's argument that the 1/n expansion fixes no growth rate should carry over to the test field: it is the
  same parity induction with δc = 0. That has not been checked symbolically here. If it holds, μ∞ is selected beyond all
  orders at the light cone, like Δ∞.

The data nonetheless point to a simple value:
- μ(D) falls from −0.365 to −0.990 between D = 4 and 67.
- n²(μ + 1) is 45.4, 45.3 and 44.1 at n = 28.1, 43.1 and 65.0.
- Least-squares fits in powers of 1/n on D ≥ 8.25 or D ≥ 16.76, with 2 to 4 terms, give μ∞ from −1.057 to −0.995
  (`ghost_mu_table.log`). That is a spread over fit choices, not checked against planted truth.

So μ∞ = −1 (**conjecture**), apparently with a vanishing 1/n term. Without gravity the inner test equation is the heat
equation ψ_s = ψ_tt, whose self-similar polynomial solutions 1, t, t² + 2s, … have μ = 0, −1, −2, …. The conjecture says
gravity's O(1) inner coupling leaves the rate of ψ = t unchanged at D = ∞. No mechanism for that was found. If it holds,
with γ∞ = 1/2, the fade-out exponent tends to 1 and δ_Q to (n + 1)/2.

There are two ways to test it:
- extend μ(D) to D = 78 with CLD's states and apply `hp/planted_dinf.py`'s planted-truth protocol to the extrapolation;
- solve the linear Stokes problem for the test field on the inner background (new code on CLD's `solver/inner_stokes`).

CLD's `matter/LITERATURE.md` §8 already lists μ(D) as a "cheap side test".

## 2. Other avenues

### 2.1 Exactly one unstable mode in every D

*Done since: proved for the symmetric wormhole in every D ≥ 4, §7 and `wormhole_largeD/one_mode.md`.*

- **Known.** [GGS09a]'s abstract claims one growing mode for every member of the D = 4 static family. Their proof, read
  through a summarizing fetch, uses the zero mode 1/Ψ₀ and covers γ₁ ≠ 0. For the symmetric member (γ₁ = 0, Paper III's
  case) that zero mode is not normalizable, and uniqueness is asserted "by continuity". CG §1.1 says the whole family is
  proved. [CPS20] treat the symmetric case in D = 4 by a Sturm argument; that is known here from a search summary only,
  not read. [TS13] find one negative mode numerically in D = 4–11, 20, 50 and 100.
- **What it would add.** A proof for every D, carrying CPS20's Sturm count, or a nodal count of a static zero mode such as the
  family's dilation, over to Paper III's master equation. Two searches found no D-dimensional proof.
- **Cost.** About a day: the D-dimensional Schrödinger form, the node count, and sympy checks with controls.
- **Value.** Paper III's type I reading needs exactly one mode. That would become **proved** instead of numerical
  evidence, and the γ₁ = 0 gap would close in D = 4 as well.

### 2.2 Paper III's next order, and a derivation of the third term

- **Known (this repository).** W²a² = 2(D − 3) − 4 ln 2 + Ω₁/(2(D − 3)) + O(D⁻²). Ω₁ = 2π²/3 + 8 ln²2 is identified by
  PSLQ at 10⁻⁹. The remainder is about 73/(2(D − 3))².
- **What it would add.**
  - Ω₁ derived. The order-3 matching integral involves K(u) = ∫₀^u ln(1 + t²)/(1 + t²) dt, which reduces to
    polylogarithms.
  - Ω₂ computed from the order-4 integral to about 20 digits and identified against weight-3 constants (ζ(3), π² ln 2,
    ln³ 2, …).
- **Cost.** A day for each. `matched.py` and `order3.py` generate the near-zone orders mechanically.
- **Value.** Modest: the expansion already matches the numerics to 4×10⁻² at D = 20 (Paper III's table).

### 2.3 The type I lifetime law across D

- **Known.** Paper III's Table 3 compares the predicted slope 1/(√2 W) with [SH02] at D = 4 and with two-run slopes from
  [ST17] at D = 5 and 6, all consistent. No bisection exists in any D. The threshold mass is unsettled: 0.30 a [SH02]
  against 0.22 a [GGS09b].
- **What it would add.** A bisection with lifetime fits in D = 5 and 6, and the D-dependence of the threshold mass.
- **Cost.** High. Neither repository can evolve through a throat (CG §5.2: in an areal-radius gauge a minimal sphere
  shows up as a → ∞). It needs a two-ended 1+1 code in double-null or proper-radius form, days of work before the first
  run.

### 2.4 Charged and asymmetric wormholes (branch `wip/wormhole-general`)

- **State of that branch** (not on main):
  - asymmetric growth rates at D = 4–40, with a large-D fit whose free slope matches the leading √ω to 4×10⁻⁸ (γ₁ = 1)
    and 4×10⁻⁷ (γ₁ = 0.5);
  - GGS09a's Table I reproduced at D = 4;
  - for the symmetric charged wormhole, the pair of real rates merges into a complex, oscillating growing pair at ε₂*,
    with ε₂*(D − 2)² extrapolating to 1.00.
- **Next there, cheap.** Derive the merge point's leading-order limit, ε₂*(D − 2)² → 1, from `leading.py`'s leading
  equation (hours), then write it up. This belongs to that branch's session and is not duplicated here.

### 2.5 Hayward's statements in D dimensions

- **Known.** BW §1's slope formula holds in every D. For any warped product M² × Sⁿ the Einstein tensor's ++ component
  is G₊₊ = −(n/r)∇₊∇₊r, which with g₊₋ = −e^{−f} gives BW's expression for every n (BW checked n = 2, 3 and 4). Nozawa and
  Maeda prove the signature law in n dimensions for Einstein–Gauss–Bonnet gravity [NM08, Proposition 10]. In the GR
  branch, under the null energy condition, an outer trapping horizon is non-timelike.
- **What it would add.** Nothing to compute. BW §1 should cite [NM08] for D > 4.

### 2.6 Two-field and evolution tests (CG §4.4–4.5)

- **Boost-rotating self-similar solutions** (CG §5.3 item 3). These are the O(1,1) real form of Hirschmann–Eardley's
  equations, with CLD's `css_lateral.py`. Their existence is unknown and the search takes days. It is the only route
  found to a genuinely two-field critical solution.
- **A two-field 1+1 bisection.** This repository's `phase5/spherical` code with a second, ghost field. It is the direct
  nonlinear test that a ghost admixture leaves γ and Δ unchanged. It takes days and cannot follow throats.
- **Non-spherical ghost modes** (ψ ∝ Y_ℓ, ℓ ≥ 1). These are still test fields, with the ℓ(ℓ + n − 1)/r² term added to
  §1.1's equations.
  - [MGG99] show that all non-spherical perturbations of the coupled problem decay at D = 4. The ghost's ℓ ≥ 1 test
    field is a different operator and has not been computed.
  - It is expected to decay faster than ℓ = 0 (**conjecture**).
  - It needs the ℓ term in CLD's crate or an independent f64 shooting code, about a day.

## 3. Ranked list

| rank | step | what it settles | cost | what exists |
|---|---|---|---|---|
| 1 | Correct BW §3 and §6.2, CG §2.3, §4.3 and §5.3, and PLAN with CLD's counts and the pilot's μ (draft in §5) | the stale "untested" | an hour of editing, no CPU | §1.2, §1.3 |
| 2 | Finish the ghost sector on CLD's code (§4) | the escape, as numerical evidence on all of 4 ≤ D ≤ 78; BW's exponent and δ_Q as functions of D; a test of μ∞ = −1 | about 10 h of CPU at 2–4 threads; no new code | `ghost_mu.py`, CLD's `modecount_test.py` and `modecount_test_axis.py`, states at D = 4–67 in CLD's `hp/data` |
| 3 | One unstable mode in every D for the symmetric wormhole (§2.1) | Paper III's type I reading on a proof | about a day, analytic with sympy checks | Paper III's master equation, [CPS20]'s D = 4 argument |
| 4 | Non-spherical ghost modes, ℓ ≥ 1 (§2.6) | whether a non-spherical admixture also fades | about a day | §1.1's equations |
| 5 | μ∞ from the inner problem's linear Stokes condition (§1.5) | an independent μ∞, and whether it is exactly −1 | days, new code on CLD's `solver/inner_stokes` | the inner background, §1.5's equations |
| 6 | Paper III: Ω₁ derived, Ω₂ computed and identified (§2.2) | the next order | a day each | `matched.py`, `order3.py` |
| 7 | The charged pair's merge, ε₂*(D − 2)² → 1, at leading order (§2.4) | a closed form for the merge | hours | `wip/wormhole-general` (its session) |
| 8 | Two-field 1+1 bisection near Choptuik's threshold (§2.6) | γ and Δ unchanged by a ghost admixture, nonlinearly | days | `phase5/spherical` |
| 9 | Boost-rotating self-similar solutions (§2.6, CG §5.3) | a genuinely two-field critical solution, if one exists | days, existence unknown | CLD's `css_lateral.py` |
| 10 | Type I lifetime bisection in D = 5 and 6 (§2.3) | Paper III's lifetime law, properly | a week, with a two-ended evolution code first | nothing |
| 11 | Hayward in D dimensions (§2.5) | nothing new | a citation | [NM08] |

The D = 5 discrepancy (§1.4) is not a step of its own. Our side is checked, and step 2's other-sector scans and boxes at
D = 5 are the remaining diagnostics.

## 4. Recommendation

Start with step 2, after the one-hour step 1.

It is the question the escape argument actually rests on. It runs on code and states that exist, and the pilot has
validated the method against [GMG96], CLD's own runs, a planted zero and a wrong operator. It also gives the two numbers
BW and a charged-collapse reader want, 2|μ|γ and δ_Q, at every D. Concretely, from CLD's `hp/` with its `.venv`, with
2 threads per job, at most two jobs, nice 10:

1. Winding counts in the D gaps (about 100 to 350 evaluations each, 0.5–2 h at M 48):
   `python3 modecount_test.py <state> --threads=2 --data=<folder> --allow_unraw`, on D = 5, 6, 7 (`T79-win/down_D*`)
   and D = 52, 60 and 67 (`hp/data/states/up9_D52.0414.npz`, `up13_D60.0020.npz`, `hp/data/large_d/up70_D67.0115.npz`).
   Each run gets a `--plant` control.
2. The leading-mode boxes, from just right of the real zero to the origin, at D = 5, 16.76 and 45.05. Run the same script
   with `--boxes=B0,U0 --aLL=<μ + 0.03>` and a `--plant` inside.
3. The other parity sector: `modecount_test_axis.py <state> --s=1 --lo=-1.6 --hi=-0.02` at the pilot's D.
4. μ(D) at about 10 more D between 50 and 78 with `wormhole_next/ghost_mu.py`. Then extrapolate μ∞ with the
   planted-truth protocol (synthetic μ(n) with a known limit through the same fit code, criteria written first) and test
   μ∞ = −1.

Output belongs either here (logs beside the scripts, as in the pilot) or in a CLD raw-data folder if the coordinator
makes it a CLD task. CLD's writers need `--allow_unraw` for the first. Step 3 of the ranked list, the one-mode proof, is
the best next analytic task and can run in parallel, since it needs no CPU.

## 5. Proposed text for BW, CG and PLAN (for the coordinator)

These are drafts. The numbers come from §1.3.

- **BW, Summary, fourth bullet.** Replace the last two sentences by: "Its energy at horizon formation, relative to the
  black hole's, falls like (p − p*)^{2|μ|γ}, with μ the least-damped ghost-sector rate: exponent 0.273 in D = 4, 0.366 in
  D = 5, rising to 0.95 at D = 45 (wormhole_next.md §1.3). So the first horizon is always a black-hole horizon. The escape,
  a ghost-sector mode that grows in some D, has been looked for and not found at D = 4, 8.25, 16.76, 30.07 and 45.05
  (choptuik-large-d, notes/GHOST_ADMIXTURE.md; numerical evidence)."
- **BW §3, the list of known values.** "μ = −0.3651 at D = 4 (choptuik-large-d; [GMG96]: −0.362 ± 0.012), −0.4430 at
  D = 5 ([MGM26]'s evolutions give −0.448; see wormhole_next.md §1.4), and the values at D = 6–67 in wormhole_next.md
  §1.3." In item 5, use the exponents from that table. The closing sentence becomes: "It rests on μ < 0 for every
  ghost-sector mode except the boost: numerical evidence from counts at five values of D and real-axis rates at 12 from
  D = 4 to 67."
- **BW §6.2.** "Tested since (numerical evidence). choptuik-large-d's argument-principle counts of the test field's
  determinant find only the boost for −0.1 ≲ Re μ ≲ 6.2 at D = 8.25, 16.76, 30.07 and 45.05, and near the imaginary axis
  at D = 4 (notes/GHOST_ADMIXTURE.md there). The least-damped real rates at D = 4–67 are in wormhole_next.md §1.3. Still
  open: counts for 4 < D < 8.25 and D > 45, and the region between the real zero and Re μ = −0.1, which decides whether
  the real zero is the least-damped mode." Also note that the computation used CLD's `hp/modecount_cx`
  (test = TEST_ON), not a frozen-metric `dss_ell.py`.
- **CG §2.3, item 3.** Replace "In D = 5 the leading imaginary-perturbation rate is μ = −0.448 [MGM26] … not re-read" by
  "In D = 5, [MGM26] fit μ = −0.448 to evolutions (window −0.450 to −0.446). The test-field determinant gives −0.4430
  (wormhole_next.md §1.4)."
- **CG §4, item 3.** Replace "Unknown in D > 5: whether the decay … persists" by "Tested at linear order in D = 4–67
  (choptuik-large-d GHOST_ADMIXTURE.md; wormhole_next.md §1.3): no growing ghost-sector mode in the regions counted. The
  nonlinear test (an evolution) and the counts listed in wormhole_next.md §1.2 remain."
- **CG §5.3, item 1.** Mark it done, with the pointer to CLD's GHOST_ADMIXTURE.md.
- **CG §1.1, first row** ("proven for the whole static family"). Say instead that GGS09a's proof covers γ₁ ≠ 0, that the
  symmetric member is argued by continuity there, and that it is treated by [CPS20]. Read CPS20 before fixing the wording.
- **PLAN.md, "Later", "Untested prediction".** "Tested at linear order (numerical evidence, D = 4–67): no ghost-sector
  mode other than the boost grows (choptuik-large-d notes/GHOST_ADMIXTURE.md; notes/wormhole_next.md). Next:
  wormhole_next.md §4."

## 7. Step 3: one unstable radial mode in every D (proved)

The proof is in [`wormhole_largeD/one_mode.md`](wormhole_largeD/one_mode.md), and its checks are in `one_mode.py`
(log `one_mode.log`: all pass, all controls fire).

- **What is proved.** The massless, reflection-symmetric Ellis wormhole has exactly one unstable spherically symmetric
  mode in every D ≥ 4. The mode is even across the throat, and there is no odd one.
- **What is assumed.**
  - Paper III's master equation in the gauge δφ = 0, checked symbolically there for symbolic n.
  - Standard spectral facts for −d²/dl² + U with U bounded and vanishing at infinity (refs. in one_mode.md): self-adjointness,
    σ_ess = [0, ∞), and exponential decay of bound states.
- **How the proof goes.**
  - With Y = rZ, q = Z′/w and χ = s q, where s = r′ r^{(n−2)/2} and w = s², Paper III's equation becomes
    −χ″ + Uχ = −W²χ with U = (n − 2)(n − 4)/(4r²) − 3n²/(4r^{2n}), regular at the throat. At D = 4 this is [CPS20]'s
    −3/r⁴.
  - The odd zero-energy solution s is positive for l > 0, so the odd sector has no unstable mode.
  - The even zero-energy solution s·H has exactly one positive zero. For D ≥ 6 this rests on the sign of a Beta function,
    H(∞) = −B(1/2 − 1/m, −1/2)/(2m) < 0 with m = D − 3. For D = 4 and 5, H diverges to −∞.
  - A Wronskian comparison lemma and a variational argument then give exactly one even negative eigenvalue.
- **Checks.**
  - Each identity is checked in sympy with a control: S1–S5.
  - Numerically (N1–N4): the Beta formula against the integrated ODE at m = 3–25; the zero counts at D = 4–50; and the
    bound state's −E against Paper III's W² to 3×10⁻⁹ (D = 4, 10, 20, 50). The control, a wrong potential, moves it by
    0.14–11%.
- **Relation to the literature** (read for this):
  - [GGS09a] prove D = 4 for γ₁ ≠ 0 and argue γ₁ = 0 by continuity.
  - [CPS20] sketch γ₁ = 0 with a non-normalizable zero-energy solution.
  - Two searches found no proof for D > 4.
  - critical_ghost §1.1 overstated GGS09a and is corrected.

## Files

All in `wormhole_next/`. `ghost_mu.py` runs under CLD's `.venv` (it needs hpkernel and `modecount_cx`), as the
queue scripts do. The other two need only this repository's setup (numpy, sympy).
- `ghost_mu.py`: real-axis zeros of CLD's test-field determinant (mode 1) or its control (mode 2), with a planted zero
  (`--plant`).
- `ghost_mu_table.py`: the table of §1.3, its six checks and the μ∞ fits. Run it as
  `python3 ghost_mu_table.py [--cld=… --raw=…]`; it exits 1 because check 3 fails. Its output is in
  `ghost_mu_table.log`.
- `inner_test_field.py`: the D = ∞ test-field equations of §1.5 (sympy, with a control).
- `run_pilot.sh` (queues a and b), `run_pilot2.sh` (c and d), `run_pilot3.sh` (e, the second layout), `run_pilot4.sh`
  (f, the line Im μ = ω), `run_pilot5.sh` (g): the queues as run. Also `logs/`, `ghost_mu.jsonl` and
  `modecount_test_axis.jsonl`.
- The second-layout state is not kept, since `*.npz` is ignored here. `run_pilot3.sh` remakes it.

## References

- [GMG96] C. Gundlach, J. M. Martín-García, "Charge scaling and universality in critical collapse", PRD 54, 7353 (1996),
  arXiv:gr-qc/9606072 (as in BW). μ₁ = 0.362 ± 0.012 with τ = ln(t/r₀), that is μ = −0.362 here.
- [MGM26] S. Mittal, C. Gundlach, L. Martel, "Critical and extremal gravitational collapse of a spherical charged scalar
  field in 4+1 spacetime dimensions", arXiv:2609.17207 (2026; title from arXiv, journal none). Read for this note through
  the arXiv HTML page via a summarizing fetch, so the quotes are close but not guaranteed byte-exact.
  - μ is "the least negative λᵢ associated with χ", taken as the value that makes χ_rot e^{−μτ} "most strictly periodic";
    "μ ∈ [−0.450, −0.446]", "the midpoint μ = −0.448".
  - Eq. (66): δ_Q = (D − 2 − μ)γ. Also γ ≃ 0.41325, Δ ≃ 3.224, δ_Q ≃ 1.425.
  - No error bars; μ only at D = 5.
- [GGS09a], [GGS09b], [SH02], [ST17], [TS13]: see CG.
- [CPS20] F. Cremona, L. Pizzocchero, O. Sarbach, "Gauge-invariant spherical linear perturbations of wormholes in
  Einstein gravity minimally coupled to a self-interacting phantom scalar field", PRD 101, 104061 (2020),
  arXiv:1911.13103 (journal from the DOI on the arXiv page). Its body was not read.
- [NM08] M. Nozawa, H. Maeda, "Dynamical black holes with symmetry in Einstein–Gauss–Bonnet gravity", arXiv:0710.2709
  (journal *(unverified)*). Proposition 10 read through ar5iv via a summarizing fetch.
- [MGG99] J. M. Martín-García, C. Gundlach, "All nonspherical perturbations of the Choptuik spacetime decay", PRD 59,
  064031 (1999), arXiv:gr-qc/9809059 (journal from a search summary of the arXiv listing).
- CLD: choptuik-large-d at main fece0ba, in particular notes/GHOST_ADMIXTURE.md, the README (Proposition 5) and
  matter/LITERATURE.md §8.
