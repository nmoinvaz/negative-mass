# The linear spectrum of the symmetric Ellis wormhole beyond its unstable mode, D = 5 and D = 4 (#272)

Driver: [ellis_spectrum.py](ellis_spectrum.py) (self-contained; uses `notes/wormhole_largeD/w2_crate` only for the planted W²). Output of record: [ellis_spectrum_output.txt](ellis_spectrum_output.txt)
(`python3 ellis_spectrum.py run`, exit 1: 8 of 57 checks fail, listed under "Checks"). Raw output, with the criteria written before the runs of record:
negative-mass-raw `I272-mac` (README, `dev/`, `runs/`, `runs/mut/`). Status words as in AGENTS.md: every number is *numerical evidence* unless marked *checked symbolically*; no sentence says a result is new.
Units: throat radius a = 1, proper time t = τ (g_tt = −1 on the background), rates per unit t.

## Answer

**No pole of E2's resolvent lies in the region of κ, at either D, in either parity (on the sheet reached from Re W > 0, |W| > 0.2): the rates 0.495 (D = 5) and 0.244 (D = 4) are not quasinormal modes.** Beyond the unstable mode the resolvent of the master equation has
(i) a pole at W = 0 (simple at D = 4, double with a logarithm at D = 5), (ii) a cut along negative real W (a power-law tail), and (iii) a few zeros hugging the cut at |Re W| ≥ 2.6, far faster than κ.
The fitted κ is read here as the local exponent of a power law in the fit window, not a rate (*numerical evidence* on the linear response; the nonlinear data were not re-fitted, and the second sheet and |W| < 0.2 were not scanned).

| acceptance (issue #272) | verdict | numbers |
|---|---|---|
| the transient's rates against κ = 0.495 ± 0.01 (D = 5) and 0.244 (D = 4) | **no mode matches** (S2); the response is a power law (T7), its local rate passes through 0.245 at t = 10 at D = 4 and is 0.03 at D = 5 | no zero with Re W ∈ [−1.4, −0.2], Im W ∈ [0, 1] (D = 5, 4; even, odd); the #255 estimator on the linear response gives k = 0.338 (D = 4), 0.185 (D = 5) |
| the scale mode's secular term against ṙ∞/E_net at B = 0.02, or the reason no bounded Y exists | **no bounded Y exists** (Z2, Z3); the secular speed is v ≈ 0.28 P at t = 10 (D = 5, slowly falling), none at D = 4 (the shift saturates at Y∞ ≈ P); ṙ∞/E_net is not B independent (T6) | `Y^(1)` grows as l³ ln l / 2 (D = 5), l³/2 (D = 4); ṙ∞/E_net = 0.0136, 0.0231, 0.0476, 0.068 at B = 0.02, 0.04, 0.07, 0.1 |
| planted: #195's W² at D = 5 to 1e-8 | **met** | 8.0e-14 (D = 5), 1.4e-13 (D = 4) by shooting; 6.4e-5, 7.6e-5 by the collocation solver |
| planted: Pöschl–Teller, first three overtones, to 1e-8 | **met** | 7.9e-14 (n = 0, 1, 2: even, odd, even); the second solver 2.4e-7 |
| a second method for the least-damped quasinormal mode at D = 5 to 1e-4 | **not met as worded** | the zero nearest the axis, W = −4.2347 + 0.3743 i (even), can be exposed only by a rotation |θ| > 1.5, where the collocation solver is not usable; its zeros agree over four ray angles and two far points to 1.6e-6 (one code); the collocation solver confirms the unstable mode and finds no other isolated eigenvalue (S4) |
| #267's decay of the mass excess at D = 4 | compared, **not reproduced as a rate** (T2m) | #267: 0.205 to 0.222 per proper time; the linear response's local rate is 0.216 at t = 12, 0.182 at t = 15, falling as 3/t |

## The linear system, and what its zero modes are

Gauge and notation of `notes/wormhole_largeD/derive.py`: n = D − 2, throat radius 1, r₀′² = 1 − r₀^{−2(n−1)}, r₀″ = (n−1) r₀^{1−2n}, h = diag(−e^{2ν}, e^{2λ}), ν = εN, λ = εL, r = r₀ + εY, δφ = 0.
Keeping the time dependence (no e^{Wt}), the first-order equations (ellis_spectrum.py `sym_system`) are:

| equation | content |
|---|---|
| E1tl | L_t = Y_{tl}/r₀′ |
| S | (N − L + nY/r₀)′ = 0 |
| E1tt | n(r₀′N′ − Y_tt) + r₀(N″ − L_tt) = 0 |
| E1ll | n r₀(r₀′L′ − Y″) + n r₀″ Y + r₀²(L_tt − N″) = 0 |
| E2 | 2(n−1) r₀′(r₀′L − Y′) + r₀(r₀′L′ − r₀′N′ + 2 r₀″L − Y″ + Y_tt) − r₀″ Y = 0 |

Integrating E1tl and S: **L = Y′/r₀′ + g(l), N = L − nY/r₀ + c(t)**, with g a function of l alone and c(t) a time reparametrization. E2 then reads **Y_tt = D Y − 2(n−1) g/r₀**,
with D Y = Y″ − pY′ + qY, p = 2r₀″/r₀′ + n r₀′/r₀, q = n r₀′²/r₀² + r₀″/r₀: the operator of `master.py`. E1tt and E1ll contain only g (*checked symbolically*, n = 2 and 3: the difference is 0, and
E1tt = E1ll = 0 at g = 0; Z1). The solutions with g = 0 are exactly the solutions of Y_tt = DY. The static solutions:

| mode | Y | g | what it is | secular partner |
|---|---|---|---|---|
| scale | Y₁ = r₀ (L = 1, N = 1 − n) | 0 | r → λ r(l/λ) | t Y₁ solves all five equations |
| asymmetry | Y₂ = r₀ ∫₀ˡ r₀′² r₀^{n−2} dl (∼ l³ at the throat, ∼ l^n/(n−1) at infinity) | 0 | the second parameter of the static family: the Misner–Sharp mass perturbation δm = (n−1)Y/(2r₀^n) + r₀^{n−1}r₀′²g tends to +½ at +∞ and −½ at −∞ (by hand, not in the script) | t Y₂ solves all five |
| translation | Y_φ = r₀ⁿ r₀′ | −(n−1) r₀^{1−n}/r₀′ | l → l + c with φ → φ + c: a gauge mode, absent from E2's W → 0 limit because g ≠ 0 | none: t Y_φ is not a solution (residual ∝ t − 1, Z2) |
| time | N = const | 0 | t → (1 + ε)t | |

Each of the first three is a solution of all five equations with its own N and L (sympy residual exactly 0 at n = 2, 3; Z2). `t² Y₁/2 + Y^(1)`, with D Y^(1) = Y₁ (the derivative of the even throat solution with respect to W² at 0), is one too.

## The generalized kernel, and why no bounded Y exists

An exact solution `δ = tZ + Y` needs D Y = 0, so Y ∈ span{Y₁, Y₂}; these grow as l and lⁿ (table). So there is no bounded Y in the gauge's variable, and the chain extends upward (Y^(1) is unbounded too).

| function | Y | u = Y/(r₀′ r₀^{n/2}) (canonical variable) | where |
|---|---|---|---|
| Y₁ (scale) | ∼ l | ∼ l^{1−n/2}: l^{−1/2} (D = 5), const (D = 4) | l → ∞ |
| Y₂ (asymmetry) | ∼ l^n/(n−1) | ∼ l^{n/2} | l → ∞ |
| Y^(1) (D Y = Y₁, Y(0) = 0) | l³/2 (D = 4); ½ l³ ln l + O(l³) (D = 5) | ∼ l^{3−n/2} ln l | l → ∞ |

Z3: Z/l³ = 0.4991, 0.5000, 0.5000 at l = 10², 10³, 10⁴ (D = 4, predicted ½); at D = 5 it is 2.015, 3.133, 4.279 and the slope per e-fold over the last decade is 0.4978 (predicted ½).
In the canonical variable Y₁ is bounded and decays at D = 5 but is not square integrable (∫u² ∼ ln l: a zero-energy resonance, marginal), and at D = 4 it tends to a constant.

What a localized perturbation does instead is set by the pole of the resolvent at W = 0. With Δ(W) = Y_out′/Y_out − Y_th′/Y_th at l_m = 2 (the mismatch of the throat and the outgoing solutions), d ln|Δ|/d ln W
at real W = 0.01 … 0.08 is 1.03, 1.06, 1.10 at D = 4 (a simple pole of the Green function: a static shift) and 1.77, 1.73, 1.69 at D = 5 (Δ/W² = 11.7, 9.9, 8.2, 6.6 against ln(1/W) = 4.6, 3.9, 3.2, 2.5: a double pole with a logarithm: t/ln t growth; Z4).

## The poles of the resolvent

Even (Y(0) = 1) and odd (Y ∼ l³) throat solutions, outgoing at infinity along a ray from l_m = 2 (the ray angle chosen so that e^{−Wl} decays along it), Wronskian M(W) analytic in W; zeros counted by the argument principle on boxes
(Re W ∈ [−8.2, −0.2], Im W ∈ [0, 8]), refined to 0.2 and polished by secant (S2, S3); on Re W > 0 the boxes around and beside the unstable root have windings 1, 0, 0 (S1).

| D | parity | unstable W | zeros with Re W < 0 | spread over four ray angles and two far points |
|---|---|---|---|---|
| 5 | even | 1.727703 (W² = 2.9849589) | −4.2347 + 0.3743 i | 1.6e-6 |
| 5 | odd | none | −6.781 + 0.467 i (candidate) | 1.8e-1: **not a zero** |
| 4 | even | 1.181970 (W² = 1.3970524) | −2.6187 + 0.1527 i; −5.6668 + 0.2475 i | 2.2e-10; 3.5e-5 |
| 4 | odd | none | −4.1310 + 0.2149 i; −7.218 + 0.266 i (candidate) | 6.9e-8; 1.8e-2: **not a zero** |

- **None lies in the region of κ** (Re W ∈ [−1.4, −0.2], Im W ∈ [0, 1]; S2 met in four of four). The confirmed zeros decay in less than 0.4 time units (rates 2.6 to 5.7).
- The two candidates move with the ray by 1.8e-2 and 0.18: S3 fails for the odd sectors as the README predicted; they are not reported as modes.
- The region |W| < 0.2 and the second sheet (Im W < 0 continued through the cut) are not covered. Zeros hugging the cut at small |W| would be invisible to a contour rotation (it needs |θ| above π/2 − arctan(Ω/γ), and the continuation beyond the axis is the other sheet).
- **Complex scaling** (Chebyshev collocation on l = z(s), rotation −0.7, −0.85, −1.0, N = 481, S = 26, Dirichlet ends, parity reduced) exposes the unstable mode (relative 6.4e-5 at D = 5, 7.6e-5 at D = 4; W2) and no other isolated eigenvalue with |W²| < 80 in either parity (S4);
  the rest is the rotated continuum. This cannot reach the zeros near the cut (their Ω/γ ≤ 0.09 needs |θ| > 1.48; at −1.53 and −1.555 the eigenvalues nearest the targets are continuum), so the shooting zeros are confirmed by one code only.

## The time domain

E2 as Y_tt = DY on l > 0 (staggered grid, Y even about 0, RK4, dt = 0.4 h, the unstable mode's discrete left and right eigenvectors used to remove its component each step), data Y(0, l) = 0, Y_t(0, l) = Gaussian
(centers 3, 6, 1.5, widths 0.5, 1). Y(t, 0) is the throat radius perturbation. Second order in h (observed order 2.00, 1.99 at D = 4, 1.99, 1.97 at D = 5 over h = 0.02, 0.01, 0.005; T1); the series are Richardson combinations of 0.01 and 0.005 (correction 4e-4 of Y at D = 4, 3.6e-3 at D = 5).
P = ∫ w Y₁ σ₀ dl, w = 1/(r₀′² r₀ⁿ), the momentum conjugate to the scale mode, for data far from the throat; Y(t, 0)/P is the same function for centers 3 and 6 to 1.8e-3 (D = 4) and 7.7e-3 (D = 5) at t = 20 (T4).

| t | D = 4: Y/P·P, v, rate = −d ln v/dt | D = 5: Y, v, rate |
|---|---|---|
| 8 | 0.4323, 5.13e-3, 0.283 | 0.4088, 3.08e-2, 0.032 |
| 10 | 0.4403, 3.03e-3, 0.245 | 0.4687, 2.92e-2, 0.024 |
| 12 | 0.4451, 1.92e-3, 0.216 | 0.5258, 2.80e-2, 0.019 |
| 15 | 0.4494, 1.06e-3, 0.182 | 0.6077, 2.66e-2, 0.014 |
| 20 | 0.4530, 4.7e-4, 0.144 | 0.7367, 2.51e-2, 0.010 |
| 25 | 0.4547, 2.4e-4, 0.119 | 0.8592, 2.40e-2, 0.008 |

(data: center 3, width 0.5, P = 0.4577 at D = 4, 0.1036 at D = 5; v = dY/dt.)

- **D = 4: the shift saturates and the approach is a power law.** Y(30)/P = 0.9954: the throat radius gains P (the simple pole). The speed falls as t^{−p} with p = rate × t = 2.7, 2.9, 3.0 at t = 15, 20, 25 (T7: rate × t = 2.74, 2.89, 2.97, spread 8 % against the 50 % of an exponential).
  Its local rate is 0.245 at t = 10 (κ = 0.244), 0.216 at t = 12, 0.182 at t = 15 (#267's 0.205 to 0.222 lies between t = 12 and 14): it is not a constant, so an exponential fitted in [6.5, 12] returns whatever the window gives (k = 0.338 here, T2-D4 fails against 0.244 ± 0.024),
  and the late rate depends on when the data arrived (0.144, 0.166, 0.144, 0.140 at t = 20 for the four data; spread 0.025 against the bar 0.02: T3-D4 fails).
  *Tested in the nonlinear data (#291, [throat_tail.md](throat_tail.md); numerical evidence):* family G at D = 4, B = 0.02, the mean of eight pairs per h. The throat's speed falls as a power law with rate × t = 2.11 at t = 10 and 2.32 at t = 16, not as an exponential, with an exponent about 0.4 below the linear response's (2.74 at t = 15). The data end at t = 18.2, where the roundoff offset of p*ₕ takes over, so the linear t⁻³ is neither confirmed nor excluded at t = 20 to 25.
- **D = 5: the speed is nearly constant and decays as a slow power.** v/P = 0.28 at t = 10, 0.24 at t = 20 (T5-D5 fails: a constant speed would have rate 0); rate × t = 0.215, 0.205, 0.198 at t = 15, 20, 25 (T7), so v ∼ t^{−0.2}. The exponential fitted in [6.5, 12] gives k = 0.185, against 0.495 (T2-D5 fails).
  This is consistent with the double pole with a logarithm (a t/ln t growth has a speed that falls slowly), not with a Jordan block of constant speed: the speed falls by 14 % between t = 10 and 20. The match to t/ln t is not tested beyond the sign and slowness.
- **Against #255's ṙ∞.** E_net is second order in the amplitude, so no linear map takes it to P. The one consequence the linear theory has is that a speed linear in the injected energy has ṙ∞/E_net independent of B; #255's table gives 0.0136, 0.0231, 0.0476, 0.068 at B = 0.02, 0.04, 0.07, 0.1 (a factor 5.0; T6 fails). The speed measured at D = 5 (ṙ∞ ∝ E_net^{3/2}) is not the linear secular term;
  at B = 0.02 the linear term would need P = ṙ∞/(v/P) = 7.5e-4, which is 0.049 E_net (0.24 E_net at B = 0.1): P/E_net would have to grow as B.

## Checks

57 checks, 49 pass, 8 fail (exit 1). The 8 are the ones the README predicted before the runs of record: S3 for the two odd sectors, T2-D5, T2-D4, T2m-D4, T3-D4, T5-D5, T6. They are results; no bar was moved.

| group | checks | verdict |
|---|---|---|
| planted: PT1 to PT3, W1 to W3 at D = 5, 4 | 9 | all pass |
| zero modes: Z1 (2), Z2 (seven candidates at n = 2, 3: six solutions and the one that must fail), Z3 (2), Z4 (2) | 20 | all pass |
| spectrum: S1 (2), S2 (4), S3 (4), S4 (4) | 14 | 12 pass; S3 fails for D = 5 odd (spread 0.18) and D = 4 odd (1.8e-2): a candidate zero in each moves with the ray |
| time domain: T1 (2), T2 (2), T2m, T3 (2), T4 (2), T5 (2), T6, T7 (2) | 14 | 8 pass; T2-D4, T2-D5, T2m-D4, T3-D4, T5-D5, T6 fail |

Controls (`--mutate=NAME`), each computed in development before it was declared (README), and on the runs of record (`runs/mut/`) each fires exactly its declared set:

| mutation | fires (record) |
|---|---|
| `wrong_p` (n → n + 1 in the Y′ coefficient) | W1, W2, W3 at D = 5 and 4 |
| `pt_V` (V₀ = 2.6) | PT1 |
| `ray_flip` | PT1, PT2 |
| `ecs_noscale` (θ = 0) | PT3 |
| `wrong_secular` (L = 0 in t Y₁) | Z2 `t Y₁` at n = 2, 3 |
| `pole_off` (odd throat solution in Z4) | Z4 at n = 2, 3 |
| `no_projection` (the unstable mode kept) | T4, T5-D4, T7 at both D (T2, T2m, T5-D5, T6 fail without it) |

Discrepancies and things to know:
- T7's message prints "67 %" as what an exponential's rate × t does over t = 15 to 25; the right figure is 50 % ((25 − 15)/20). Text only; the bar is 12 %.
- The README says "eight candidates each" for Z2; there are seven (six solutions, one that must fail).
- The development series were first read at h = 0.02 on a time grid whose output step was 12 steps of 0.008 (0.096, labelled 0.1): the D = 4 rate then appeared to plateau at 0.22 to 0.24. On the exact grid and Richardson series (second order, measured) the rate falls (0.283 at t = 8 to 0.119 at t = 25). The plateau was an artifact; the dev logs `dev_02` and `dev_03` hold both.
- PT3 first failed in development (N = 201, rotation −1.2: 6.1e-2 for the third overtone, which is barely exposed); N = 481, S = 24 gave 2.4e-7 and is the setting of record. Z3 at D = 5 first used the slope over two decades (0.4918, 1.6 % off); the last decade (0.4978) is the declared window.
- Two of my own development processes were stopped by name (`pkill -f`) before I switched to PIDs. Both were mine.


## What `throat_drift.md` and `critical_ghost.md` §1.3 should now say (proposals for the coordinator, not applied)

**`throat_drift.md`, the row "a generalized zero mode (a Jordan block) of the static family; a slow tail of the ghost field; a drift of the lapse relative to the far region"**, replace the verdict "not established" by:
*the Jordan-block reading is ruled out as the cause of a constant speed; a slow tail is present in the linear problem.* #272 (*numerical evidence* and *checked symbolically*): the scale and asymmetry zero modes and their secular partners t·Y are exact linear
solutions, `Y^(1)` (D Y^(1) = Y₁) is unbounded (l³ ln l / 2 at D = 5, l³/2 at D = 4), so no bounded Y exists; the resolvent has a double pole with a logarithm at W = 0 at D = 5 and a simple pole at D = 4, which in the time domain is a throat speed
that falls slowly (v ∼ t^{−0.2} at D = 5, 28 % of P at t = 10 and 24 % at t = 20) at D = 5 and a shift that saturates with a t^{−3} speed at D = 4; and ṙ∞/E_net is not B independent (0.0136 at B = 0.02, 0.068 at B = 0.1), so the speed is not linear in the injected energy.

**`throat_drift.md`, recommended task 1** is done (#272). **The transient's rate κ** (0.495 at D = 5, 0.244 at D = 4): no zero of the linear resolvent lies in Re W ∈ [−1.4, −0.2], Im W ∈ [0, 1] (D = 5 and 4, even and odd); the zeros that exist have rates 2.6 to 5.7. At D = 4 the linear throat speed
falls as t^{−3} (rate × t = 2.7, 2.9, 3.0 at t = 15, 20, 25) and its local rate passes through 0.245 at t = 10, 0.216 at t = 12 and 0.182 at t = 15: an exponential fitted to such a power law returns the window's rate. κ is not a mode (not shown for the nonlinear data, whose arrival time of the pulses differs).

**`critical_ghost.md` §1.3**, the bullet on the generic one-sided family, add after "the cause is not established": *Linear analysis (#272, `wormhole_evolution/ellis_spectrum.md`; numerical evidence): the symmetric wormhole's linear
spectrum beyond the unstable mode has no discrete mode with a rate near the transient's κ (0.495 at D = 5, 0.244 at D = 4); the late linear response of the throat radius is a power law (speed ∼ t^{−3} at D = 4, with the shift saturating; ∼ t^{−0.2}, nearly constant, at D = 5, from a double pole with a logarithm at W = 0), so the
transient's "rate" is the local exponent of a tail and the near-constant D = 5 speed is the threshold response of the scale mode. It is not linear in the injected energy (ṙ∞/E_net varies by 5 over B = 0.02 to 0.1).*


## Limits

- The linear problem is the even sector for the time domain (the throat radius); the odd sector has zeros and the candidates above but no time-domain run.
- The response depends on the data: the linear speed and rates are for Gaussian σ₀ and are not the response to #255's pulses (a matter source of second order). The statement is about the structure of the response, not its amplitude.
- The cut region |W| < 0.2, the second sheet and the lower half plane (the mirror images) were not scanned; zeros there would not change the structure above.
- The near-cut zeros rest on the shooting code alone (checked over ray angle and far point, not by a second code).
- #255's fits are of the nonlinear runs over t ∈ [6.5, 12]; the linear response here has a different arrival time of the data, and the local rate is time dependent: the agreement of the D = 4 local rate with 0.244 at t = 10 is a coincidence of the window, not shown to be more.

## Recommended tasks

1. **The linear response to #255's pulses at B = 0.02.** What: put the pulses' stress as a source into `Y_tt = DY + S` (the ghost pulse A* at first order, the normal pulse at second order in B) and evolve with this note's time-domain code. Why: P cannot be got from E_net; at B = 0.02 the amplitudes are small enough (r_th − 1 = 0.002) for the linear theory to hold, and #255 measures ṙ∞ = 2.1e-4 there. Accept (can fail): the linear speed over t ∈ [6.5, 12] within 30 % of 2.1e-4 and its B scaling between B = 0.02 and 0.04 within 20 % of the measured exponent 2.8; a flipped source sign must fail it. Depends on: this note, #255.
2. **The second sheet and |W| < 0.2.** What: continue M below the negative real axis (rays beyond −π/2 at l_m = 2, θ down to −2.4) and scan Re W ∈ [−1.4, −0.1], Im W ∈ [−0.4, 0.05] at D = 4, 5, with a planted pole just across a cut (a potential with a c/l⁴ tail and a known second-sheet pole) through the same code. Why: a pole just across the cut near −0.2 to −0.5 would give a quasi-exponential for a window and is the one candidate for κ left open here. Accept: the planted pole is found to 1e-6; zero or one zero per box reported at each D. Depends on: this note.
3. **The tail in the nonlinear data.** What: G's pair at D = 4, B = 0.02, to t = 25 (the pair at 1e-12), and the log-log slope of the throat speed over t ∈ [15, 25]. Why: this note predicts rate × t = 2.7 to 3.0 (a t^{−3} speed) and not a constant rate. Accept: the slope within 20 % of −2.9 against an exponential's growing rate × t. Depends on: #255, this note.
4. **An independent code for the near-cut zeros.** What: confirm W = −4.2347 + 0.3743 i (D = 5, even) and −2.6187 + 0.1527 i (D = 4, even) with a method that shares nothing with the shooting code (a Leaver-type series about the throat, or hyperboloidal collocation). Why: complex scaling cannot reach them and the issue's second-method criterion is open. Accept: 1e-4, and the same planted Pöschl–Teller overtones through it. Depends on: this note.


## Files

`ellis_spectrum.py`, `ellis_spectrum_output.txt`, this note; raw `negative-mass-raw/I272-mac` (README: method, criteria, mutations; `dev/`; `runs/`, `runs/mut/`).

## Coordinator's review

Merged at df53984. Review folder `negative-mass-raw/I272-mac-review` at 52f6e7bf, criteria committed (e7b556e5) before the first run of
`ind272.py`. The changed files are those the brief names.

- **Rerun on the Mac** from 6373c3c, with `w2_crate` built in the review tree: `ellis_spectrum.py run` equals the stored output apart
  from timings and the header's commit (exit 1, the same 8 failures). A first attempt without the crate built stopped at two stages
  and is kept.
- **ind272**, 4 of 4. It shares only the formula of D, the data and the definition of P with this note's code. It has its own
  integrator (velocity Verlet), outer boundary (Dirichlet at l = 75), eigen-solve and rate estimate. Its W² matches #195's to 2.5e-8,
  which anchors D. It reproduces the time-domain table: Y(t, 0) to 1e-4, the D = 4 local rates to 3.5e-3, rate × t = 2.735, 2.886,
  2.965 (a power law), and v/P = 0.2815, 0.2421 at D = 5. The mutations `wrong-n`, `no-proj` and `width` fire.
- **Not checked independently:** the zero search (S1 to S4) and the symbolic zero-mode checks (Z1 to Z4) rest on this note's code.
  The second-method criterion stays not met, as worded.
- On T7's message: with the spread taken over the mean, as T7 computes it, an exponential gives 50 %, so the note's correction is right.
  The 67 % is the rise from t = 15 to 25.
- Applied: the proposed text for `throat_drift.md` and `critical_ghost.md` §1.3, and a paragraph in Paper V's §on the drift.
