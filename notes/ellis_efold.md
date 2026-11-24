# The Ellis family's linear e-folding time at m/a = 1/2 (#151)

November 2026. Code `ellis_efold.py`, output `ellis_efold_output.txt`, criteria and logs in negative-mass-raw `I151-rpi5`
(README written before the first run; the run numbers below are its logs). Run it from this directory with the repository's
Python (numpy, scipy, mpmath, sympy): `python3 ellis_efold.py > ellis_efold_output.txt`, about 2.5 minutes on one core.

Status words as in AGENTS.md: the unit map, the background, the linearization and the printed potential are **checked
symbolically** (30 digits where numeric); every eigenvalue below is **numerical evidence** (three methods, one spread).

## Result

**[SN26]'s linear time T = 0.7584 is [GGS09a]'s T at γ₁ = m/a = 1/2, and it is right.** Solving the linear problem here gives
β = 1.875127511869 in [GGS09a]'s units (b = 1, γ₀ = 0) and **T = 0.75835715**, which is 0.006 % below 0.7584 and
0.0004 above [GGS09a] Table I's 0.758. In [SN26]'s time (evolution time t, units of m) the e-folding time is τ_lin = 5.1308 M
against their 5.13 M. The error bar is 5.6e-12 (a spread over discretizations and methods, below), far inside the 0.2 %
the issue asks for.

| quantity | here | source |
|---|---|---|
| γ₁ = 0, β², [GGS09a] units | 1.39705243 | [TS13] −1.397 (n = 4); repo 1.39705243 (`one_mode.py`, N4) |
| γ₁ = 0, T = τ/r_throat | 0.846045 | [GGS09a] Table I 0.846 |
| γ₁ = 1/2, β | 1.875127511869 | none printed |
| γ₁ = 1/2, T = τ α_th / R⋆ | 0.75835715 | [GGS09a] Table I 0.758; [SN26] linear 0.7584 (T here is 0.006 % below) |
| γ₁ = 1/2, τ in [SN26]'s t | 5.1308 M | [SN26] τ_lin = 5.13 M |
| [SN26] measured τ = 5.12, 5.23 → T | 0.75676, 0.77302 | [SN26] 0.76, 0.77; −0.211 % and +1.933 % from T here (they print 0.2 % and 1.9 %) |
| 16 values γ₁ = 0 … 2.0 | all within 0.0005 of Table I (worst 0.00048) | [GGS09a] Table I (3 digits) |

## Normalizations

All from the papers' TeX (negative-mass-raw `I151-rpi5/src_*`); pages and equation numbers are the arXiv PDFs'.

**[GGS09a]** (arXiv:0806.0608 v2). Metric Eq. (13) p. 3 with Eqs. (11)–(12):
−e^{2γ₁ arctan(x/b)+2γ₀} dt² + e^{−2γ₁ arctan(x/b)−2γ₀}[dx² + (x²+b²) dΩ²]. Throat at x = γ₁ b, with
r_throat = B √(1+γ₁²) e^{−γ₁ arctan γ₁}, B = b e^{−γ₀} (p. 4). Masses m_{+∞} = B γ₁ e^{−γ₁π/2}, m_{−∞} = −B γ₁ e^{γ₁π/2}.
Rescaled to b = 1, γ₀ = 0 (p. 6) the mode is χ e^{βt} with −χ_ρρ + Wχ = −β²χ, ∂_ρ = e^{2γ₁ arctan x} ∂_x, W of Eq. (32) p. 7.
The time scale is (p. 10, the unnumbered formula before Table I):

```math
T = \frac{\tau_{\rm unstable}}{r_{\rm throat}} = \frac{e^{2\gamma_1\arctan\gamma_1}}{\sqrt{1+\gamma_1^2}\,\beta},\qquad
\tau = e^{\gamma_1 \arctan\gamma_1}\,t \;\text{ the proper time at the throat}.
```

At γ₁ = 0 this is T = 1/β with r_throat = b = 1, and β² = 1.397 is [TS13]'s ω²a² (the n = 4 entry of `critical_ghost.md` §1.1).
(p. 10 writes r_throat = b√(1+γ₁²)e^{−γ₁ arctan γ₁} with b where p. 4 has B. They agree in the units b = B = 1 of §III B.)

**[SN26]** (arXiv:2610.13206 v2). Eq. (2) p. 2: X = (r − a²/4r)/a, u = (m/a)(arctan X − π/2), and Eq. (3): α = e^u,
γ_ij = e^{−2u} Ω² δ_ij. Throat radius R⋆ = √(a²+m²) exp[(m/a) arctan(a/m)] (p. 2; 3.8895 at a = 2, m = 1). Sec. IV B p. 7:
T ≡ τ α_th/R⋆, α_th = e^{u(r_t)} = 0.575, "the γ₁ = m/a = 0.5 member of [GGS09a]'s family (their metric with b = a,
γ₀ = −πγ₁/2)", "T = 0.758 [GGS09a]", "T = 0.7584, i.e. τ_lin = 5.13 M" from their own solvers (their `linear_mode.py` is not in
the arXiv source and was not read).

**The map** (`ellis_efold.py` U1, sympy, a and m symbolic; mutation `units` flips γ₀ and makes it fail). With ρ_E = r − a²/4r
and x = ρ_E, [SN26]'s metric is [GGS09a] Eq. (13) with b = a, γ₁ = m/a, γ₀ = −πγ₁/2, hence B = a e^{πγ₁/2}:

| | [GGS09a] | [SN26] |
|---|---|---|
| ADM mass of the + end | B γ₁ e^{−γ₁π/2} | m |
| throat areal radius | B√(1+γ₁²) e^{−γ₁ arctan γ₁} | R⋆ |
| lapse at the throat | e^{γ₁ arctan γ₁ + γ₀} | α_th = e^{u(r_t)} |
| time | T̃ = e^{2γ₀} t/b, mode e^{βT̃} | t, proper time of the + end (α → 1 at X → ∞); τ = a e^{πγ₁}/β |

So τ α_th/R⋆ = e^{2γ₁ arctan γ₁}/(√(1+γ₁²) β): [SN26]'s T and [GGS09a]'s T are the same quantity, the proper time at the throat
over the areal throat radius. [SN26]'s τ is an evolution time in units of m; T = τ α_th/R⋆ takes the factor α_th/R⋆ = 0.1478 per m.
With the printed numbers, 5.13 × 0.575 / 3.8895 = 0.7584 (U1), so [SN26] applied the same definition.

## What ran

| method | what it solves | shares with the others |
|---|---|---|
| A, shooting | −χ_ρρ + Wχ = −β²χ, W = printed Eq. (32); DOP853 rtol 1e-13, inward from ρ̃ = ±30 | the background, ρ(x) |
| B, matrix | the same operator, sixth-order differences on h = L/1500, W = −V + 2(∂Ψ₀/Ψ₀)² derived with sympy from Eq. (24) and (28) at 40 digits | the background, ρ(x) |
| C, gauge δΦ = 0 | no W: the linearized Einstein–ghost equations derived from the Ricci tensor with sympy, reduced to one ODE for C = δc, decaying solutions continued across the throat on semicircles in the complex x plane | the background functions only |

Method C's ODE, with c' = (x−γ₁)/(1+x²), d' = γ₁/(1+x²):

```math
C'' = -\frac{2x}{1+x^2}\,C' + \beta^2 e^{-4\gamma_1\arctan x}\,C + \frac{2}{x-\gamma_1}\Big(C' - \frac{\gamma_1}{1+x^2}\,C\Big),
```

with the throat x = γ₁ a regular singular point (exponents 0 and 3). The ODE and the reconstruction δa = [C' + (c'−d')C]/c', δd =
δa − 2C satisfy the tt, tx, xx, θθ and scalar equations to 1e-26 (L1). Both semicircles give the same roots and a real
determinant to 1e-12 (no monodromy). Methods A and B work in ρ̃ = λρ, λ = e^{2γ₁ arctan γ₁} ([GGS09a]'s β̄ scaling, p. 9), only for
conditioning; the values returned are in ρ units.

**What the error bar is.** At γ₁ = 1/2, β varies by at most 3.5e-13 (A: rtol 1e-13 → 1e-11, domain +10), 7.4e-12 (B: N 3000 →
6000, domain +10) and 2.1e-13 (C: semicircle radius, domain) over those solver choices, and A, B, C differ by 5.8e-13 and
1.0e-13. The reported ±5.6e-12 on T is the largest of these. It is a spread over discretization and solver choices at fixed
physics, not a rigorous enclosure. The root-finder tolerance is 1e-14. The potential itself is *not* in the spread: W is
checked separately (W1: the printed Eq. (32) equals −V + 2(∂Ψ₀/Ψ₀)² to 1e-30, and method C, which has no W, gives the same β).

## Checks and controls

The final run (`notes/ellis_efold_output.txt`) passes 46 checks (rc 0, 2 min 16 s, peak rss 128 MB). Every check has a
mutation that must make it fail. The 12 mutated runs are summarized in `negative-mass-raw I151-rpi5/mutations_summary.md`
(generated by `summarize_mutations.py` from the logs; the note's numbers are checked by `check_note_numbers.py`).
The 12 controls:

| mutation | what it breaks | target checks | result |
|---|---|---|---|
| `potential` | W + 0.1/(1+x²)² in A and B: **the γ₁ = 0 control fails** (T = 0.8657, β² = 1.3343) | G0 G2 G3 G4 G5 G6 | all six fire; nothing else fails |
| `printedW` | one coefficient of Eq. (32), A only | W1 G1 G2 G3 G4 G5 | fire (G6 also fails, not predicted) |
| `matrix` | B's centre stencil weight | P-B G1 G0 | fire (G7, D1 also fail) |
| `planted` | the exact answer β = ℓ | P-A P-B | fire |
| `xmap` | the sign in A's x(ρ) | P-A G1 G2 G3 G4 G5 | **first run silent on G5**: it crashed in the table loop (no root) before G5 ran. Rerun after the code reports an exception as a failed check: fires |
| `units` | γ₀ = +πγ₁/2 in the [SN26] map | U1 | fires |
| `tfactor` | T without e^{2γ₁ arctan γ₁} | U2 G3 G4 G5 G6 | fire; G0 passes, as predicted (γ₁ = 0) |
| `count` | a second well in the G7 count | G7 | fires |
| `domain` | L = 4 | D1 G0 G1 | **first run silent on D1**: crashed in the table loop. Rerun fires (P-A, P-B, G2, G5, G7 also fail) |
| `deltaphi` | the 2 in method C's ODE | L1 G2 | fire (D1 also fails) |
| `cpath` | semicircle radius 0 | P-C G2 | **first two runs hung** (solve_ivp on a NaN right-hand side), stopped by me. After the code rejected a non-finite right-hand side: fires |
| `background` | κΦ₁² = −2(1+2γ₁²) | S1 | fires |

Two things to know. A control that crashes or hangs is not a check that failed: those runs are kept as the silent first
runs, and the code now reports an exception inside a check as that check's failure (commits fcbad54, f9e3cee; no
criterion, tolerance or target list changed, and the final output equals run 4's apart from the timing line). The
README's predictions that a check would *not* fail were all right; the failures outside the target lists are in the summary.

## What `critical_ghost.md` should now say

These are proposals; the file is the coordinator's. Line numbers are those of `critical_ghost.md` at main cde03d9.

```diff
@@ §1.1, row 1, last sentence (line 92)
-Time scale τ/r_throat = 0.846 (massless symmetric) → 0.590 (large asymmetry γ₁)
+Time scale T = τ/r_throat, τ the proper time at the throat, T = e^{2γ₁ arctan γ₁}/(√(1+γ₁²) β) in the units b = 1, γ₀ = 0:
+0.846 (γ₁ = 0), 0.758 (1/2), 0.675 (1), 0.618 (2) → 0.590 (large γ₁). The first row of Table I is reproduced at 16 values to
+its 3 digits (`ellis_efold.md`, #151); the limit 0.590 is quoted from [GGS09a].

@@ §1.2, [SN26] bullet (line 139)
-    - one exponential mode, proper-time e-folding time 0.76–0.77 R⋆ measured against its linear 0.7584;
+    - one exponential mode, proper-time e-folding time T = τ α_th/R⋆ = 0.76–0.77 (α_th the throat lapse, R⋆ the throat radius:
+      the same T as [GGS09a]'s, their γ₁ = m/a = 1/2) measured against its linear 0.7584; recomputed here as 0.758357 (#151);

@@ §1.3, [SN26] bullet (lines 177–178)
-... Its linear e-folding time is not reconciled
-  here with [GGS09a]'s normalization.
+... Its linear e-folding time 0.7584 is [GGS09a]'s T at γ₁ = m/a = 1/2, recomputed as 0.758357 with τ_lin = 5.1308 m in
+  its evolution time (their 5.13); its measured 0.76 and 0.77 are 0.21 % below and 1.93 % above it (#151, `ellis_efold.md`).
```

The n = 4 remark of §1.1 ("1/√1.397 = 0.846, checked here") now rests on three methods: β² = 1.39705243, T = 0.846045.
Nothing else in `critical_ghost.md` needs to change.

## Not claimed

- Nothing about the second row of [GGS09a] Table I (γ₁ ≥ 2.5, the limit 0.590): the code stops at γ₁ = 2. The value 0.590 is quoted from [GGS09a].
- [SN26]'s nonlinear runs and its own linear solver were not re-run or read; only its stated definitions and numbers were checked.
- [TS13] was not re-read; its 1.397 is checked against this repository's own D = 4 value and against T(0) = 0.846.
- Uniqueness of the unstable mode (one bound state) is checked numerically at γ₁ = 0, 1/2, 1, 2 (G7: one negative eigenvalue in B,
  one root in A's scan). The proofs are [GGS09a] (γ₁ ≠ 0), `wormhole_largeD/one_mode.md` (γ₁ = 0).
