# The Ellis family's linear e-folding time at m/a = 1/2 (#151)

November 2026. Code `ellis_efold.py`, output `ellis_efold_output.txt`, criteria and logs in negative-mass-raw `I151-rpi5`
(README written before the first run; the run numbers below are its logs). Run it from this directory with the repository's
Python (numpy, scipy, mpmath, sympy): `python3 ellis_efold.py > ellis_efold_output.txt`, about 6 minutes on one core with #163's
additions (`--section=old` is #151's 2.5 minutes). #163 (second row of Table I and the large-γ₁ limit) is the last part of this note.

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

- (#151, superseded by #163 below) The second row of Table I and the limit were not computed in #151; #163 does.
- [SN26]'s nonlinear runs and its own linear solver were not re-run or read; only its stated definitions and numbers were checked.
- [TS13] was not re-read; its 1.397 is checked against this repository's own D = 4 value and against T(0) = 0.846.
- Uniqueness of the unstable mode (one bound state) is checked numerically at γ₁ = 0, 1/2, 1, 2 (G7: one negative eigenvalue in B,
  one root in A's scan). The proofs are [GGS09a] (γ₁ ≠ 0), `wormhole_largeD/one_mode.md` (γ₁ = 0).

## Coordinator's review

Merged as 9de9a04. Raw negative-mass-raw `I151-mac` at bf3b3e5 (criteria at 75d19a6, before the first run).

| step | result |
|---|---|
| scope | the three files the brief names, and raw `I151-rpi5` |
| rerun on the Mac | exit 0, 25 s, 46 of 46 checks. The output equals rpi5's apart from the timings and B's last digits (LAPACK). β_A is the same to 12 digits. B's spread at γ₁ = 1/2 is 8.7e-12 here, so the error bar of T is 6.6e-12 on the Mac against 5.6e-12 on rpi5 |
| re-solve | β at γ₁ = 1/2 by two own routes: 1.8751275118668 (printed W, Numerov) and 1.8751275118690 (no W) |
| independent check `ind151.py` | E8 fails as written; the other nine criteria pass; 10 of 10 mutations fire on target in run 2 |

ind151 shares no code with `ellis_efold.py` and uses no scipy. Route 1 solves the printed W (retyped from the TeX) by
Numerov with Richardson extrapolation. Route 2 derives the linearized equations in the gauge δΦ = 0 from the
coordinator's own Ricci tensor and reduces them to a first-order system for C = δc and G = δa − δc:
C′ = d′C + c′G and (x − γ₁)G′ = PC + QG. It crosses the throat on the real axis by Frobenius series (exponents 0 and 1)
and needs no complex path.

| id | what | result |
|---|---|---|
| E0 | planted Ricci | pass |
| E1 | background, ghost sign κΦ₁² = −2(1 + γ₁²) | pass |
| E2 | the [SN26] → [GGS09a] map, M₊ = m, R⋆, α_th, the T identity | pass |
| E3 | the δΦ = 0 reduction satisfies all five linearized equations; #151's printed ODE gives the same G′ | pass |
| E4 | typed W = −V + 2(∂Ψ₀/Ψ₀)², ∂²Ψ₀ = VΨ₀ | pass |
| E5, E6 | planted Pöschl–Teller (route 1) and a planted analog of the throat structure with exact β = 9/2, 3/2 (route 2) | pass |
| E7 | route 1 against #151's β_A at γ₁ = 0, 1/2, 1 | pass, rel ≤ 3.1e-12 |
| E8 | route 2 against route 1 | **fail**, see below; β agrees to ≤ 2.9e-12 |
| E9 | T(0), T(1/2), T(1) round to Table I; T(1/2) = 0.758357150 rounds to [SN26]'s 0.7584 | pass |

**E8.** It fails on one clause: \|Q₀ − 1\| = 2e-23 at γ₁ = 1/2, where the bound was 1e-25. Q₀ = 1 exactly by hand, and
the derived Q is 1 to 5e-31 at the throat. The 2e-23 is aliasing of the check's own 64-point Cauchy sum,
(1/2 ÷ √(1 + γ₁²))⁶⁴ ≈ 4e-23, so the bound was tighter than the method allows. The other clauses of E8 pass: route 2's
β agrees with route 1 to 1.2e-12 at γ₁ = 1/2, there is one sign change at each γ₁, and the throat has no logarithm. The
no-log residual of the exponent-1 resonance is ≤ 1e-22 at three γ₁ and two β. That last point is an independent
confirmation of method C's trivial monodromy.

**Findings.**
- The note's proposed §1.3 sentence said "its measured 0.76 and 0.77 are 0.21 % below and 1.93 % above it". The
  percentages belong to the unrounded T = 0.7568 and 0.7730; the rounded 0.76 is above 0.7584. `critical_ghost.md` §1.3
  says this.
- [SN26]'s `numbers.tex` gives the level-3 τ unrounded as 5.119 m. Then T = 0.7566, 0.23 % below the linear value; the
  0.2 % [SN26] prints still holds.
- The three controls that crashed or hung in their first runs (`xmap`, `domain`, `cpath`) are recorded as silent there,
  with logs kept. The fix reports an exception as the check's failure, which is stricter, and no criterion moved.

`critical_ghost.md` §1.1, §1.2 and §1.3 now carry the proposals, with the §1.3 correction marked as the coordinator's.


# #163: the second row of Table I and the large-γ₁ limit

November 2026. Same code and output files; criteria, amendments and logs in negative-mass-raw `I163-rpi5` (README with the criteria at 91e4990, before
the first run of the new code; amendments 1–8 there). The final run (`ellis_efold_output.txt`) runs 100 checks in 6 min 15 s on one core (rpi5, peak rss
134 MB, script commit d3d2030, exit 1: **two checks fail, P-C2 and LE5, both as registered and reported below**). Status words as in AGENTS.md: the identities N2, LE0, N3a are
*checked symbolically* or to 16 digits; every eigenvalue is *numerical evidence* (three methods, one spread); the limit of T(γ₁) is *numerical evidence*.

## Result

**All 16 entries of the second row of [GGS09a] Table I are reproduced** (γ₁ = 2.5 … 10): worst |T − printed| = 0.000489, and every value rounds to the printed
digits. **The limit equation, solved here from the TeX's Eq. (43), gives 1/β̂ = 0.588351371**, which rounds to the paper's 0.588 and is 0.00165 below
the table's 0.590; the issue's criterion "within 0.001 of 0.590" fails as registered. The finite-γ₁ time keeps falling after γ₁ = 8, where the table's
digits stop changing: T(10) = 0.589697, T(160) = 0.588357, with the gap of β̂ to its limit falling like γ₁⁻² (measured order 1.990, 1.997, 1.999, 2.000). So 0.590 is T
at γ₁ ≈ 8–10 rounded, not the limit.

| quantity | here | source |
|---|---|---|
| 16 values γ₁ = 2.5 … 10 | all within 0.00049 of Table I (worst 0.000489 at γ₁ = 4.0 and 6.5) | [GGS09a] Table I, second row |
| error bar of T (second row) | ≤ 4.5e-12 (spread over discretizations and solvers at γ₁ = 2.5, 5, 10; not an enclosure) | D1b |
| 1/β̂ of Eq. (43), three solvers | 0.588351371 (shooting and matrix 3e-13 apart; finite volume with Richardson 5e-11 from them) | [GGS09a] Sec. IV: "yields 1/β̂ = 0.588" |
| 1/β̂ against 0.590 | −0.001649 (criterion ≤ 0.001: **fails**) | [GGS09a] Sec. IV, after Table I: "converges to a value of 0.590" |
| T(10), T(20), T(40), T(80), T(160) | 0.589697, 0.588689, 0.588436, 0.588373, 0.588357 | this run (LE6) |
| first row, unchanged | legacy path equals #151's β_A at 16 values with rel 0.0 (R1); the new path equals them to 2.6e-14 (converged, L = 60) and, at γ₁ ≤ 1.8, to 2.3e-13 (as run, L = 30). At γ₁ = 2.0 it differs from #151's as-run value by 6.29e-6 | R1, R2 |

## The far field: what the problem was

**The printed W does not lose its digits.** The issue and #151 assumed that the bracket of Eq. (32) loses its digits where x ~ ρe^{πγ₁} is large. Measured: at
γ₁ = 2.5 and 10, x = 10 … 1e12, the bracket in double precision agrees with the derived W (60 digits) to 8.7e-16 (DIAG N1, `prefeas_conditioning.log`). The cancellation
of −3/q against 3(x−γ₁)²/q² is against a term 2(F_x/F)² of the same size.

**What stopped #151's code at γ₁ ≥ 2.5 was the position of the throat.** #151 works in ρ̃ = λρ, λ = e^{2γ₁ arctan γ₁}, with the origin at x = 0 and the domain ρ̃ ∈ [−30, 30]. The mode
lives around the throat, at ρ̃ = ∫₀^{γ₁} e^{2γ₁(arctan γ₁ − arctan y)} dy = 2.2 (γ₁ = 1), 23.8 (2), **84.8 (2.5)**, 321 (3), 5.2e3 (4), 1.8e6 (6), 7.1e8 (8), 3.0e11 (10).
The legacy path finds no root at γ₁ = 2.5 and 3.0 (DIAG L0). At γ₁ = 2 the throat is 6 units from the domain's end, which is why #151's γ₁ = 2.0 value is truncated (below).

**Handling.** The new path (`gshoot`, `gmatrix`, `C2_det`) puts the origin of u = (ρ̃ − ρ̃_throat)/sc, sc = max(1, γ₁), at the throat (x = γ₁), domain |u| ≤ 20:

```math
-\chi_{uu} + sc^2\,a(x)^2\,[\text{bracket of Eq. (32)}]\,\chi = -b^2\chi,\quad \frac{dx}{du} = sc\,a(x),\ x(0) = \gamma_1,\quad
a(x) = e^{2\gamma_1[\arctan x - \arctan\gamma_1]},\quad T = \frac{sc}{\sqrt{1+\gamma_1^2}\,b}.
```

Then b = sc β̄ = γ₁β̄ = β̂ for γ₁ ≥ 1, the mode is O(1) wide in u at every γ₁, and x(+L′) = 1350 at γ₁ = 10 (the x(−L), x(+L) of the output), never near 1e12. Beyond x_sw = 2 max(1, γ₁) the code
evaluates the bracket with its first two terms combined, 3(γ₁² − 1 − 2γ₁x)/q², exact algebra (N2, symbolic) and equal to the printed form to 8e-16 (N1a). It is kept as the far-field
form so that the large-x part of W is code that a control can break (`largex`); it is **not** needed for accuracy (R3: with and without it, β differs by 4.6e-15). The left start χ′ = +bχ and the
right start χ′ = −bχ ignore the 2/r² tails; their effect is below the spread of D1b (L′ → L′ + 8: 3.6e-13 at most).

## Second row

| γ₁ | b = β̂ | T | printed | T − printed |
|---|---|---|---|---|
| 2.5 | 1.526710344 | 0.608155 | 0.608 | +0.000155 |
| 3.0 | 1.574655248 | 0.602470 | 0.602 | +0.000470 |
| 3.5 | 1.605497568 | 0.598895 | 0.599 | -0.000105 |
| 4.0 | 1.626362434 | 0.596511 | 0.597 | -0.000489 |
| 4.5 | 1.641074421 | 0.594846 | 0.595 | -0.000154 |
| 5.0 | 1.651808826 | 0.593641 | 0.594 | -0.000359 |
| 5.5 | 1.659867460 | 0.592740 | 0.593 | -0.000260 |
| 6.0 | 1.666064333 | 0.592050 | 0.592 | +0.000050 |
| 6.5 | 1.670928007 | 0.591511 | 0.592 | -0.000489 |
| 7.0 | 1.674813055 | 0.591081 | 0.591 | +0.000081 |
| 7.5 | 1.677964155 | 0.590732 | 0.591 | -0.000268 |
| 8.0 | 1.680554378 | 0.590447 | 0.590 | +0.000447 |
| 8.5 | 1.682708832 | 0.590209 | 0.590 | +0.000209 |
| 9.0 | 1.684519717 | 0.590010 | 0.590 | +0.000010 |
| 9.5 | 1.686056155 | 0.589841 | 0.590 | -0.000159 |
| 10.0 | 1.687370771 | 0.589697 | 0.590 | -0.000303 |

Three methods: A (shooting, printed W), B (matrix, W derived with sympy at 40 digits) and C (δΦ = 0, no W). A against B at the 16 values and at γ₁ = 0.5, 1, 2: worst 3.6e-12 (G1b). A against C at γ₁ = 2.5, 4, 6, 8, 10:
worst 1.5e-12, upper against lower semicircle equal to 12 digits (G2b). One bound state: one root in A's scan and one negative eigenvalue in B at γ₁ = 2.5, 5, 10 (G7b).
Spread of b over solver choices (D1b; A: rtol and L′, B: N and L′, C: radius and L′): at most 7.5e-12 relative (B at γ₁ = 5), so the error bar of T is at most 4.5e-12; it is measured at three values and assumed for the rest of the row.

## The limit equation

[GGS09a] Sec. IV (TeX label `Eq:DiffEqLargeGamma`, the PDF's Eq. (43)):
−(e^{2(1−1/z)}∂_z)²χ + e^{4(1−1/z)}[2z⁻² − 10z⁻³ + 3z⁻⁴]χ = −β̂²χ, z > 0, and T → 1/β̂. Checked here: the printed Ŵ(z) of the TeX (retyped) equals γ₁² W̄(γ₁z) of the code to 6e-16 (N3a),
and tends to the printed limit potential, the gap falling by a factor of 100 per decade of γ₁, 1.7e-3 at γ₁ = 1e2 to 1.7e-11 at 1e6 (N3b); the zero mode (2/π)(z−1)/z² e^{1−1/z}
solves the printed equation exactly (LE0, symbolic). Solved by shooting and by a sixth-order matrix in u = ρ̂ (sharing z(u) and Ŵ), and by a finite-volume scheme on a uniform grid in ln z (sharing only Ŵ):

| solver | β̂ | 1/β̂ |
|---|---|---|
| A, shooting | 1.699664604104 | 0.588351371 |
| B, matrix | 1.699664604104 | 0.588351371 |
| finite volume, N = 8000 and 16000 (Richardson in h²) | 1.6996646040 (1.6996646791, 1.6996646228) | 0.588351371 |

Planted truth through the same code (sech²(ρ̂(z)), ℓ = 0.9, 1.6): A and B to 1.5e-12, finite volume to 1.1e-9 (P-LA, P-LB, P-LC). One bound state in all three (LE4).
The finite-γ₁ problem tends to the limit equation: β̂ = b at γ₁ = 10, 20, 40, 80, 160 is 1.687370771, 1.696569783, 1.698889552, 1.699470757, 1.699616137, against 1.699664604; the gap falls by a factor 4 per doubling (LE6).
Extrapolating T(80), T(160) with the measured order 2 gives 0.588351, equal to 1/β̂.

**The gap between 0.588 and 0.590, as it is.** [GGS09a] write both: 1/β̂ = 0.588 from Eq. (43) and 0.590 from the table. The limit equation gives 0.588351, so their 0.588 is reproduced and the
table's 0.590 is the finite-γ₁ value T = 0.5904 … 0.5897 at γ₁ = 8 … 10 rounded to three digits. T(γ₁) is still falling there (0.58970 at 10, 0.58869 at 20).
Their sentence "the time scale converges to a value of 0.590 r_throat" reads the three digits of the table; their own limit equation does not support 0.590.

## Checks and controls (#163)

New checks: N1a, N1b, N2, N3a, N3b, LE0, U3, P-A2, P-B2, P-C2, P-C2b, P-LA, P-LB, P-LC, R2, R3, G5b, G1b, G2b, G7b, D1b, LE3, LE4, LE5, LE6, plus R1 in #151's section.
Every one passes in the final run except P-C2 and LE5. Controls (`--mutate=`, 17 mutations; the new ones are `largex`, `centre`, `scale`, `lelimit`, `lefv`), from
`I163-rpi5/mutations_summary.md`, generated by `summarize_mutations.py` from the logs:

| mutation | target | fired | SILENT | failed outside the target list | exit, wall |
|---|---|---|---|---|---|
| `potential` | G0, G2, G3, G4, G5, G6, R1 | G0, G2, G3, G4, G5, G6, R1 | none | P-C2, LE5 | rc 1, 6:04.96 |
| `printedW` | W1, G1, G2, G3, G4, G5, R1, N1a, N1b, N3a, R2, G5b, G1b, G2b | W1, G1, G2, G3, G4, G5, R1, N1a, N1b, N3a, R2, G5b, G1b, G2b | none | G6, P-C2, LE5, LE6 | rc 1, 5:54.54 |
| `matrix` | P-B, G1, G0, P-B2, G1b, P-LB, LE3 | P-B, G1, G0, P-B2, G1b, P-LB, LE3 | none | G7, D1, P-C2, G7b, D1b, LE4, LE5 | rc 1, 6:07.36 |
| `planted` | P-A, P-B, P-A2, P-B2, P-LA, P-LB, P-LC | P-A, P-B, P-A2, P-B2, P-LA, P-LB, P-LC | none | P-C2, LE5 | rc 1, 6:11.15 |
| `xmap` (run 1) | P-A, G1, G2, G3, G4, G5, R1 | P-A, G1, G2, G3, G4, G5 | R1 | G6, G7, D1, P-C2, LE5 | rc 1, 9:01.97 |
| `units` | U1 | U1 | none | P-C2, LE5 | rc 1, 6:12.66 |
| `tfactor` | U2, G3, G4, G5, G6, U3 | U2, G3, G4, G5, G6, U3 | none | P-C2, LE5 | rc 1, 6:07.85 |
| `count` | G7, G7b, LE4 | G7, G7b, LE4 | none | P-C2, LE5 | rc 1, 13:15.35 |
| `domain` | D1, G0, G1, D1b, G5b, R2, P-A2, P-B2, G1b, LE3 | D1, G0, G1, D1b, R2, P-A2, P-B2, G1b, LE3 | G5b | P-A, P-B, G2, G5, G7, P-C2, P-LA, P-LB, G2b, LE5 | rc 1, 2:19.89 |
| `deltaphi` | L1, G2, G2b, P-C2b | L1, G2, G2b, P-C2b | none | D1, P-C2, D1b, LE5 | rc 1, 5:50.13 |
| `cpath` | P-C, G2, P-C2b, G2b | P-C, G2, P-C2b, G2b | none | P-C2, D1b, LE5 | rc 1, 4:18.68 |
| `background` | S1 | S1 | none | P-C2, LE5 | rc 1, 6:10.36 |
| `largex` | N1a, N2, R2, R3, G5b, G1b, G2b | N1a, N2, R2, R3, G5b, G1b, G2b | none | P-C2, LE5, LE6 | rc 1, 4:35.11 |
| `centre` (run 1) | G5b, G1b, G7b, D1b | none (stopped at the 1 h cap, rc 124, after N1a … P-B2) | **G5b, G1b, G7b, D1b: never ran** | | rc 124, 1:00:00 |
| `scale` | U3, G5b | U3, G5b | none | P-C2, LE5 | rc 1, 3:52.56 |
| `lelimit` | LE0, N3b, LE6 | LE0, N3b, LE6 | none | P-C2, LE5 | rc 1, 3:56.06 |
| `lefv` | LE3, P-LC | LE3, P-LC | none | P-C2, LE5 | rc 1, 3:52.07 |
| `centre` (run 2) | G5b, G1b, G7b, D1b | G5b, G1b, G7b, D1b | none | P-C2, R3, G2b, P-C2b, LE5, LE6 | rc 1, 9:28.78 |
| `xmap` (run 2) | P-A, G1, G2, G3, G4, G5, R1 | P-A, G1, G2, G3, G4, G5, R1 | none | G6, G7, D1, P-C2, LE5 | rc 1, 8:51.77 |

Every row fails P-C2 and LE5 as well, which fail without any mutation; those two are not in any target list. **16 of the 17 controls fire on every target; `domain` does not**: it is
silent on G5b (the second row stays within 0.001 of the printed values at L′ = 2, worst 0.000552, while the truncation is 5.3e-4 in β and R2, D1b, G1b, P-A2, P-B2, LE3 fire on it).
Two controls were silent in their first run and fire after a code change that touches only mutated paths (README amendment 8): `xmap` was silent on R1 (the legacy table loop raised before R1;
R1 now reports FAIL when it does) and `centre` hit the 1 h cap because its planted problems at x₀ = 0 integrate extremely slowly (the mutation now leaves planted problems alone). Both first logs are kept
(`mut_xmap.log`, `mut_centre.log`). The unmutated rerun (`final_run2`, the file `ellis_efold_output.txt`) equals the first final run apart from the timing line. Unpredicted extra failures are in the table's fifth column; none of the checks
that the README predicted would pass failed on a control, except the ones listed there (for example `matrix` also fails G7, D1, G7b, D1b, LE4, which use the same mutated B).

**The two failed checks.**
- *LE5(a)*: 1/β̂ = 0.588351 against 0.590, −0.001649, criterion 0.001. Reported as it is above; LE5(b) (rounds to the paper's 0.588) passes.
- *P-C2* (registered; stays failed): a manufactured solution C* = e^{γ₁ arctan x}e^{−(x−γ₁)²/2} carried on a semicircle of radius 0.5√(1+γ₁²) with b/sc = 1.3. Error 2.3e-11, 9.3e-8, 4.8e6
  at γ₁ = 2.5, 5, 10. Cause (`diag_pc2_growth.py`): the unwanted homogeneous solution grows along that arc by 1.7e2, 3.4e4, 1.4e9, so the test amplifies its own 1e-12 integration error. It
  is also not the real path (the real radius is 3.83 at γ₁ = 10, not 5.02). It is not a failure of method C, which agrees with A to 1.5e-12 (G2b).
- *P-C2b* (added after the coordinator's ruling, criterion written before its first run): a manufactured solution decaying to the right like the real mode, carried through exactly
  the path, radius and direction of the real roots (right segment, both semicircles, left segment), with the source from the ODE retyped independently. Passes at γ₁ = 2.5, 4, 6, 8, 10 with worst error 1.4e-11
  (criterion 1e-8). `cpath` (r = 0) and `deltaphi` (coefficient 1 in the integrated ODE) both fire it (computed before they were declared).

## What `critical_ghost.md` should now say

Proposals; the file is the coordinator's. §1.1, row 1, last sentences (line 92 at main 2455a0d):

```diff
-... 0.618 (2) → 0.590 (large γ₁). The first row of Table I is reproduced at 16 values to its 3 digits (`ellis_efold.md`, #151); the limit 0.590 is quoted from [GGS09a]
+... 0.618 (2), 0.608 (2.5), 0.594 (5), 0.5897 (10) → 0.58835 (γ₁ → ∞). Both rows of Table I are reproduced at 32 values to their 3 digits (`ellis_efold.md`, #151, #163).
+The limit is 1/β̂ = 0.588351 from [GGS09a]'s own limit equation, solved here from the TeX ([GGS09a] print 0.588); the table's 0.590 is T at γ₁ ≈ 8–10 rounded,
+and T falls like γ₁⁻² to 0.58835 (numerical evidence)
```

## Failed, and what disagrees

- The premise that the printed W loses its digits at large x is false for the bracket (8.7e-16 at x = 1e12). #151's code stopped at γ₁ = 2 because its domain misses the throat.
- #151's γ₁ = 2.0 value is truncated: β = 60.662389104967 as run (L = 30) against 60.662007295832 converged (L = 45 and 60 agree), rel 6.29e-6; T = 0.6178881 against 0.6178919. No printed digit changes.
  The criterion "the first row is unchanged at every γ₁ ≤ 2 to rel 1e-9" holds exactly for the legacy path (R1, rel 0.0), for the new path against the three β that #151 printed (γ₁ = 0, 0.5, 1, 12 digits)
  and against its other values at γ₁ ≤ 1.8. **It fails at γ₁ = 2.0 for the new path against #151's as-run value, by 6.29e-6**, because that value was truncated. The new path equals the converged value to 1.3e-14.
- 1/β̂ = 0.588351 against the table's 0.590: −0.001649 (above).
- P-C2 fails at γ₁ = 5 and 10 (above). Run 1 of the new code crashed in its first check, on a sample point that is exactly the throat of γ₁ = 10 (README amendment 1).

## Not claimed

- The convergence T(γ₁) → 1/β̂ is *numerical evidence* (LE6, five values of γ₁, order 2.000); [GGS09a]'s argument is pointwise convergence of the equation, not of its eigenvalue.
- Uniqueness of the bound state for γ₁ ≥ 2.5 is checked at three values (G7b) and for the limit equation (LE4, with the paper's nodal argument for the zero mode); the proof is [GGS09a].
- [SN26]'s nonlinear results are not touched; nothing here changes the §1.2–§1.3 statements about T(1/2).
