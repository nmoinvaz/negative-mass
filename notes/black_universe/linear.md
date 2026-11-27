# The linear double-null response of the black universes, and c = −1 (#232)

Code [linear.py](linear.py) (solver, comparators, checks), derivation [linear_derive.py](linear_derive.py) with
[linear_derive_output.txt](linear_derive_output.txt). Output of record [linear_output.txt](linear_output.txt), assembled from
negative-mass-raw `I232-mac/runs`; the criteria, written before the runs of record, are in `I232-mac/README.md`. Status words as in AGENTS.md: rates and
ratios are *numerical evidence*, the algebra of `linear_derive.py` is *checked symbolically*. No sentence here says a result is new.

## Answer

**At c = −1 the linear double-null solver grows at κ**, for data that perturb the horizon (rate 0.636627 against κ = 0.636620, 1.2e-5), and does not for waves
incident on an unperturbed horizon. The growth is carried by the horizon values of Y, which obey `Y_h(p) = c₁ + c₂ e^{bp}`, `b = √2 κ` exactly; it comes
from the term `c₂`. #212's count of zero at c = −1 is the count under the condition `u = Y/ρ → 0` at the horizon, which at ρ_h = 0 sets `Y_h = 0`, so
`c₂ = 0`. With `Y` regular and nonzero on the horizon the Y equation has the solution of rate W = κ (the exponent `1 − W/(2κ)` is smooth there,
without a logarithm, only at W = κ). **`README.md`'s statement that c = −1 has no growing mode is corrected, not confirmed**: it holds for `Y_h = 0`.

At c = −1.5 and −0.3 the linear rate equals #212's W to 9.2e-5 and 2.0e-7; the planted c = 0 case gives `W²a²` to 8.0e-6 (linear solver), 2.1e-10 (frequency-domain
master equation) and 1.4e-5 (#39's kernel, differential). The disagreement is not in the nonlinear kernel (it reproduces the linear response to 2.3e-3
while its own drift is below 5e-4), and not in the equation or in `Y` for pulse data (the 1D master equation in characteristic form reproduces the double-null `Y`
to 9e-4 of its maximum, and grows at κ at c = −1 when fed the same data). It lies in the horizon boundary condition of the mode problem, which at ρ_h = 0 is a condition on
`u = Y/ρ`. The nonlinear code's growth at κ is consistent with its truncation error seeding that sector: its drift and the linear solver's seed `Y_h` both fall by 4 per halving of h.

## Method

**The linear system.** #217's three evolution equations with the potential, linearized about the static member in sympy (`linear_derive.py` A0, A1, 40 digits):
`U_pm = B0 U + Bp U_p + Bm U_m`, `U = (δr, δf, δφ)`; the matrices are printed in `linear_derive_output.txt`. Gauge perturbations `L_ξ(member)` with arbitrary `ξ(p)`, `η(m)`
solve the system and both constraints on the exact member (B1, B2, six members) and `Y = δr − (ρR/√2) δφ` vanishes on them (B3).

**The chart.** The member in null coordinates `(p, σ)`: `p = x⁺` is the advanced tortoise time of #217's code and σ a null coordinate of the other family with
`Q(σ) = sinh(b(σ* − σ))/sinh(bσ*)`, σ* = 16. For σ well below σ* it is the tortoise `m` (`Q = e^{−bm}`); `Q` vanishes linearly at σ* (the future horizon) and is
negative beyond it (the interior). The Kruskal product `Π = e^{bp} Q = e^{2κx}` determines ρ analytically through the horizon: `ρ − ρ_h = Π E(ρ)` with `E` the solution
of the regular equation `(ln E)′ = 1/δ − 2κ/A`, `δ = ρ − ρ_h`, and `e^{−f} = −(Q′/b) e^{bp} E A/δ`. At c = 0 the chart is the tortoise one. The grid covers `p ∈ [0, P]` and σ up to
the interior point where `Π(P) = −0.05`; the horizon is a node. The chart solves the nonlinear equations to 5.6e-10 and carries `r_p = φ_p = 0`, `f_p = −√2 κ` on the horizon (C2, C2h).

**The scheme.** One cell, `N = E + W − S + h_p h_σ (B0 Ū + Bp Ū_p + Bm Ū_σ)` with the cell-centre averages of the four corners, the background at the exact cell centre, linear in `N`: a 3 × 3
solve per cell, second order. Pure-gauge solutions are reproduced with errors that fall by 4.00 per halving of h (C3: c = −1.5, −1, −0.3, 0), through the horizon and into the
interior; both constraints of an evolved solution fall by 4.0 per halving (C5).

**Data.** Characteristic data on the two rays through the corner ρ_c = 1: `δf = 0`, `δφ` a Gaussian pulse (amplitude 10⁻³, centre 3, width 0.5, zero at the corner), `δr` from the
linearized constraint along the ray with corner constants. Classes: R (pulse on the right, ingoing ray only), L (left only), G (both), A0, A1R, A1L (corner constants
`δr`, `δr_p`, `δr_σ` = 10⁻³, no pulse). A pulse leaves `Y` nonzero on the ray beyond it: the constraint makes `δr` respond.

**Comparators.** (i) The 1D master equation of #212 in characteristic form, `u_pm = −V̂ u/2`, with `u = Y/ρ` of the same data on the two rays, marched on the tortoise grid (c ≤ −1).
(ii) #217's `evolve_master1d.py` (compact data on a `t = 0` slice). (iii) #217's nonlinear kernel (`evolve.py`) on member ± the same data (ε = 10⁻⁶): the odd part of `Y` on the diagonal; the even part is the untuned member's own drift.
(iv) At c = 0 the frequency-domain master equation, shot from the throat.

## Checks (18 in `linear.py check`, 12 in `linear_derive.py`; each fires under at least one mutation)

| id | what | result |
|---|---|---|
| derive A1 | linearization against the central difference of the nonlinear right-hand sides, 40 digits | 1.8e-25 |
| derive B1–B3 | `L_ξ(member)` solves the evolution equations and constraints, six members; `Y` is invariant | ≤ 5.8e-41 |
| derive C1, C2 | `M = u_pm + V̂u/2` contains `δr, δφ, δr_p − δr_m, δf` only, and `(∂_t) M = 0` on solutions: the master equation holds up to a static term | 4.6e-41, 9.1e-40 |
| derive D1–D3 | horizon: `δr_pp + f_p δr_p = 0`; exponents `±W/A′` (c = −1.5); exponents `1 ∓ W/(2κ)` and logarithm coefficient `π w(1 − w)/8` at c = −1 | exact |
| C1 | `coef()` against the sympy matrices at 40 chart points | 6.4e-17 |
| C2, C2h | chart against the nonlinear equations across the horizon; horizon node | 5.6e-10, exact |
| C3 | exact pure-gauge solutions, ratios per halving | 4.00 at all four members |
| C5 | constraints along the evolution, ratio per halving | 3.99 to 4.00 |
| P1a, P1b, P1c | planted c = 0: `W²a²` against 1.39705243433350684 | 8.0e-6, 2.1e-10, 1.4e-5 |
| R1 | rate against #212's W at c = −1.5, −0.3 (and −1.01, −0.99) | 9.2e-5, 2.0e-7 (1.3e-5, 1.1e-5) |
| R3 | c = −1 rate against κ | 1.2e-5 |
| L1 | nonlinear kernel against the linear solver on the same data | 2.3e-3 (bar 3e-3) |
| L2 | 1D master equation (characteristic form) against the double-null `Y`, pulse data | 8.8e-4 (bar 2e-3) |
| L3a–L3c | horizon sector at c = −1; projection; control at c = −1.5 | see below |
| L4 | profile exponent at the horizon | 0.4995 (c = −1), 0.4475 (c = −1.5) |
| L5 | seeds ∝ h² | ratios 4.00, 4.00; 3.97, 4.01 |

Mutations: 18 for `linear.py` and 11 for `linear_derive.py`, each with the checks it must fire declared before its run of record (the six for A0, B3, C1, D1, D2, D3a were added after the first run of record,
README amendment A2); in the runs of record the fired sets equal the declared ones (18 of 18 and 11 of 11; `linear_output.txt`).

## Rates (Table 1)

Rate of `Y` on the diagonal x = 0 (time t = (p + m)/√2, units b = 1). "linear" is the linear double-null solver for data G, rate over t ∈ [16, 20], Richardson of h = 0.01 and 0.005.
"1D char" is #212's equation from the same data, "1D time" `evolve_master1d.py` over [45, 60]. "nl slope" is #217's lifetime law, 1/slope at h = 0.0025 (`evolve_output.txt`); "nl seed" the untuned
member's own growth over [8, 14] at h = 0.01; "nl same − lin" the kernel's odd response to the same data minus the linear solver's over [6, 10] (not asymptotic).

| c | κ | #212 W | 1D time | 1D char | linear | linear / W − 1 | nl slope | nl seed | nl same − lin |
|---|---|---|---|---|---|---|---|---|---|
| −1.5 | 0.551682 | 0.494171 | 0.49417 | 0.494217 | 0.494217 | 9.2e-5 | 0.49422 | 0.49053 | −5.8e-5 |
| −1.01 | 0.635123 | 0.633148 | 0.63315 | 0.633156 | 0.633156 | 1.3e-5 | | 0.62996 | −3.1e-4 |
| **−1** | 0.636620 | none | **−0.03702** | 0.636627 | **0.636627** | 1.2e-5 (of κ) | 0.63677 | 0.63339 | −3.2e-4 |
| −0.99 | 0.638098 | 0.640122 | | | 0.640129 | 1.1e-5 | | 0.63685 | −3.3e-4 |
| −0.3 | 0.601300 | 0.970080 | | | 0.970080 | −2.0e-7 | 0.96992 | 0.93403 | −2.8e-3 |

The 1D equation has no time-domain value at c = −0.99 and −0.3: the throat lies in the static region and `u = Y/ρ` is singular there.
At c = −1 two columns disagree, 1D time (−0.037) and linear / 1D char (κ); they differ in the data, not in the equation.

## c = −1

**Data classes** (rate of `Y(x = 0)` over [16, 20], Richardson; `Y_h(0)` the horizon value on the left ray; 1D/DN the 1D master equation's `Y` over the double-null `Y` at t = 14):

| class | `Y_h(0)` | rate of `Y(x = 0)` | rate, 1D char | 1D/DN at t = 14 |
|---|---|---|---|---|
| R: pulse on the right ray, left ray zero | 0 | 0.0115 (the seed) | −0.0086 | 0.986 |
| L: pulse on the left ray | −1.5e-4 | 0.63661 | 0.63661 | 1.0000 |
| G: both | −1.5e-4 | 0.63663 | 0.63663 | 1.0000 |
| A0: corner `δr` | 1.3e-3 | 0.63662 | 0.63662 | 1.316 |
| A1R: corner `δr_p` | 0 | 0.63662 | 0.0047 | 0.0004 |
| A1L: corner `δr_σ` | 2.5e-3 | 0.63660 | 0.63662 | 2.75 |

For class R, `Y_h` is zero in the continuum; the values printed in the output are the truncation seed (L5). The pulse classes have `M = 0` (derive C2) and the 1D equation reproduces them;
the corner constants carry a static `M ≠ 0` and it does not (A1R: the horizon sector reaches the 1D equation only through `δr_p` along the left ray, which `u` on the rays does not hold).

**The horizon sector.** On the horizon `r_p = φ_p = 0`, so the constraint is `δr_pp + f_p δr_p = 0` (derive D1) with `f_p = −b`: `δr_h = c₁ + c₂ e^{bp}`, and at c = −1 `Y_h = δr_h`
(the factor `ρR/√2` vanishes). For three left pulses `Y_h` follows this form with a residual that falls by 3.99 per halving of h (3.0e-4 at h = 0.005), and the late exponent
equals `b = 0.900316` to 7.4e-8 (L3a). `c₁ ≈ 10⁻¹⁰` against `c₂ ≈ 10⁻⁴`. Combining the three pulses so that `c₁ = c₂ = 0` (cancelling `Y_h` at p = 0 and 2) leaves `|Y(x = 0, t = 14)|`
smaller by 2.4e4 than one pulse, with a rate of −0.177 over [10, 14]: it decays, as #212's equation does (L3b). The same operation at c = −1.5 reduces by 40 and leaves the rate at 0.478 (W = 0.494; L3c): the
horizon sector is specific to c = −1.

**The exponents.** At a horizon with ρ_h ≠ 0 the Y equation has the exponents `±W/A′(ρ_h)` (derive D2), so `Y ~ e^{Wv}` is regular and nonzero on the horizon: the ingoing branch is the regular one.
At c = −1 they are `1 ∓ W/(2κ)` (D3a); `Y ~ ρ^{1 + W/(2κ)} e^{Wt}` (the branch `u → 0`, #212's) vanishes on the horizon, while `ρ^{1 − W/(2κ)} e^{Wt}` is smooth in Kruskal coordinates only when
the logarithm of the Frobenius recursion vanishes, with the coefficient `π w (1 − w)/8` (w = W/κ): at W = κ exactly (D3). There the branch is `e^{κv}`, nonzero on the horizon. The late-time profile of the solution near
the horizon has the exponent `0.4995` in δ = ρ − ρ_h at c = −1 (the branch `1/2`; `3/2` for the branch `u → 0`) and `0.4475` at c = −1.5, where `W/A′(ρ_h) = 0.4479` (L4).

**Where the disagreement lies.**

| candidate | check | verdict |
|---|---|---|
| the nonlinear kernel | L1: its odd response to the same data equals the linear solver's to 2.3e-3 until its own drift reaches 5e-4 (c = −1.5 and −1, data G and A0) | not the cause |
| the data | L5: the untuned member's drift (nonlinear) and `Y_h` of class R (linear) both fall by 4 per halving of h: the truncation seed carries `c₂` at order h². A right-ray pulse (class R) has `c₂ = 0` in the continuum; the truncation error does not | the nonlinear code's growth is its seed in the horizon sector |
| the master equation's variable and gauge | `Y` is gauge invariant (B3); the 1D equation reproduces the double-null `Y` to 9e-4 for pulse data at c = −1.5, −1.01, −1, growing at κ at c = −1 (L2); the equation holds up to a static term (C1, C2) | not the cause for pulse data; `u = Y/ρ` is singular on the horizon for `Y_h ≠ 0` |
| the boundary at the horizon | L3, L4, derive D3: growth is `c₂ e^{bp}` of the horizon values; the regular branch is `ρ^{1/2}`; the condition `u → 0` takes `ρ^{3/2}` | the cause |

The cause is one thing seen two ways: the condition `u → 0` is the condition `Y = O(ρ^{3/2})`, which excludes `Y_h ≠ 0`; the variable `u = Y/ρ` is what turns a regular `Y_h` into `u ∼ 1/ρ`.

**What #217's c = −1 control shows.** The control "no divergent lifetime at A = 0 along #212's family" was not met: lifetime slope 1/κ to 2.4e-4, A\*ₕ at second order. The family `member + A × δ̂`, with `δ̂` taken
from the code's own truncation-seeded departure, takes its growing part from the horizon sector (`c₂ ≠ 0`), where the linear solver puts the growth at κ; that fits the slope 1/κ.
This is read from the checks above; it is not a separate run of that family.

## Proposed text (the files are not touched here)

`README.md`, Result, first bullet: *…the member c = −1, whose horizon sits at the area minimum, has **none under the condition used here, u = Y/ρ → 0 at the horizon, which sets
Y = 0 on it**. With Y regular and nonzero on the horizon it has one, W = κ (#232; `Y_h ∝ e^{κv}`, numerical evidence, exponents checked symbolically).* Third bullet: *…for c → −1 from either side the growth rate tends to
κ, and at c = −1 it is κ for perturbations with `Y_h ≠ 0`; the table's 0 counts modes with `Y_h = 0`.* Sec. 4, the c = −1 row: add "(`Y_h = 0`); W = κ with `Y_h ≠ 0` (#232)". Sec. 5, after "the computation finds no bound state": *The solution
`u ∝ e^{−κx}`, growing toward the horizon, is the limit of the bound states; `Y = ρu ∝ e^{κv}` is regular there. #212's condition excludes it, regularity of `Y` admits it (W = κ is the one value at which the branch
`ρ^{1 − W/(2κ)}` of Y is smooth, with no logarithm).* Sec. 8, "knife edge": *The rate is κ at c = −1 as well; the member is stable only against perturbations with `Y_h = 0`, such as waves incident on an unperturbed horizon
(class R: no growth).*

`critical_ghost.md` §1.1, the [BKZ12] row: replace "The single member c = −1 … has none, so the exception is one member, not an open set. As c → −1 from either side the rate tends to the surface gravity" with
*At c = −1 there is none for perturbations with zero Y on the horizon (#212's condition); the linear double-null response grows at κ for perturbations that change the horizon's area radius (`Y_h ∝ e^{κv}`, #232, numerical evidence), so
the exception holds only against waves incident on an unperturbed horizon. The rate tends to κ as c → −1 from either side and is κ at c = −1.*

## What this does not show

- A linear statement about the horizon sector. That the sector is a datum of physical initial data (a spacelike slice through the horizon) is not tested; the characteristic data here carry it on the ray `p = 0` and through `δr_p`
  at the corner. Nothing about formation from collapse or the nonlinear fate beyond #217's runs.
- The 1D equation was not evolved at c = −0.99 and −0.3 (`u` singular at the throat); #212's W is quoted there. Radial perturbations of D = 4, one potential; nothing about polar modes (#218).
- The static term `M` of the master equation is derived (C2) and its effect shown (A0, A1R, A1L), but it is not added to the 1D equation as a source.
- Class R at c = −1 is shown to decay only to t ≈ 20 above a truncation seed that falls as h²; the 1D time-domain decay (−0.037 over [45, 60]) is #217's.

## Recommended tasks

1. **The count with a regular horizon condition.** *What:* `stability.py`'s count with both branches admitted where `Y` is smooth in Kruskal coordinates, at the 20 members. *Why:* #212's counts are for `Y_h = 0`.
   *Acceptance:* the 19 members with c ≠ −1 unchanged (rates to 1e-9); at c = −1 the regular decaying solution at W = κ with `Y ∼ ρ^{1/2}` and no root elsewhere on [0.01, 8]; mutations that fire on a wrong exponent. *Depends on:* #212, #232.
2. **#223's c = −1 member at zero charge.** *What:* the same condition for the sourced black bounce at its c = −1 member. *Why:* #223 reports no mode there by the same master equation. *Acceptance:* the count with the regular condition; its rate against κ.
   *Depends on:* #223, #232.
3. **A family that leaves the horizon sector.** *What:* #217's lifetime law at c = −1 along `member + A × (class R data)` and along a family with `c₂` tuned. *Why:* the control as worded ("no divergent lifetime at A = 0") should hold for `c₂ = 0`
   above the seed. *Acceptance:* no divergent lifetime to t = 30 for class R at h = 0.005 and 0.0025 with the seed falling by 4 per halving; the slope 1/κ along the tuned family. *Depends on:* #217, #232.
4. **The horizon sector from Cauchy data.** *What:* count the free data on a spacelike slice through the horizon at c = −1 and show whether `c₂` is free there. *Why:* decides whether the growth belongs to physical initial data.
   *Acceptance:* `c₂` as a functional of slice data, checked on the linear solver's solutions to 1e-6. *Depends on:* #232.
5. **The static term as a source.** *What:* add `M(x)` to the 1D equation and reproduce classes A0, A1R, A1L. *Acceptance:* `Y` to 2e-3 of its maximum; a wrong `M` fails. *Depends on:* #232.

## Reproduce

From this directory, `~/Source/negative-mass/.venv/bin/python`, one thread, about 4 minutes in all:

```
linear_derive.py                       # 12 checks, 3 s; --mutate=lin_Vpp|lin_cross|gauge_sign|frob_log|V2_coef|lin_const|Y_coef|M_ep|hz_rp|exp_W|exp_one
linear.py check [--only=C1,L3] [--mutate=NAME]    # 18 checks, 2 min
linear.py table                        # Table 1, 1 min
linear.py c1                           # the c = -1 tables, 1.5 min
linear.py mutations DIR PREFIX         # fired sets of stored logs DIR/PREFIX_<name>.log against the declared ones
```

`linear_output.txt` is `check`, `table`, `c1` and `mutations` of the runs of record (`I232-mac/runs`), in that order.

## Coordinator's review

Review of branch `wip/i232-c1-linear` at c4552dd and raw `I232-mac` at f5db91a, recorded in negative-mass-raw `I232-mac-review`
(criteria 889f52d; amendment A1; runs 25b0c3a). On the Mac, nice 10.

- **Scope.** The branch adds `linear_derive.py`, `linear.py`, their outputs and this note, nothing else.
- **Rerun.** `linear_derive.py` gives 12 of 12 and `linear.py check` 18 of 18. Every line of `table` and `c1` appears
  verbatim in `linear_output.txt`.
- **Independent check** (`ind232.py`). The coordinator worked it by hand, then in sympy, from #212's printed master
  equation, sharing no code.
  - At c = −1 (ρ_h = 0, A₁ = 4/π) the horizon exponents of Y are 1 ± W/A₁, so {1/2, 3/2} at W = κ = 2/π.
  - ρ^{1/2}e^{κt} equals a constant times e^{κv}(1 + O(ρ)), so the s = 1/2 branch at W = κ is V: smooth and nonzero on
    the future horizon.
  - At c = −1.01 and −0.99 the exponents at W = κ are ±1/2, so there that branch diverges.
  - So the regular solution growing at κ exists at c = −1 only. This is why #212's condition u = Y/ρ → 0, which keeps
    1 + W/A₁, removes it there and nowhere else.
  - One check missed its bar by the review's numerics (K3, amendment A1). K3b passes.
- **What it changes.**
  - The c = −1 statement of #212 and #223 is corrected, not refuted: it holds for Y_h = 0.
  - Whether the c = −1 black universe is unstable against physical perturbations turns on #248.
  - #217's control is explained: its family lies in the horizon sector.
- **Applied.** `README.md` (Result and the knife-edge item) and `critical_ghost.md` §1.1 carry the correction, marked as
  from #232.
- **Follow-ups.**
  - #248: is the horizon sector free in Cauchy data? (task 4)
  - #249: recount with a regular horizon condition, at #212's and #223's members and with the source term (tasks 1, 2, 5).
  - #250: #217's lifetime law along class R and a tuned family (task 3).
