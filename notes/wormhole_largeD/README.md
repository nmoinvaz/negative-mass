# The ghost-field wormhole's instability at large D

October 2026. A follow-up to [`../critical_ghost.md`](../critical_ghost.md) (CG) §1.1, which noted from Torii and
Shinkai's table [TS13] that the wormhole's single unstable mode grows faster in higher dimensions, roughly like
√(2D)/a. This note derives that dependence.

**The system.** The static, massless, symmetric wormhole held open by a ghost scalar (Ellis–Bronnikov) in D spacetime
dimensions, with throat radius a:

```math
ds^2 = -dt^2 + dl^2 + r(l)^2\,d\Omega_{D-2}^2, \qquad r'^2 = 1 - (a/r)^{2(D-3)} .
```

Its radial perturbations have exactly one growing mode, e^{Wt} [GGS09a for D = 4; TS13 numerically for D = 4–11, 20,
50, 100]. Since g_tt = −1, t is proper time everywhere, and 1/W is the e-folding time of the throat's collapse or
expansion.

Status words: **derived** (closed-form matched asymptotics, each step checked numerically), **computed** (numerics,
converged), **identified** (a closed form found by an integer-relation search on a computed number).

## Result

```math
W^2a^2 = 2(D-3) - 4\ln 2 + \frac{\pi^2/3 + 4\ln^2 2}{D-3} + O(D^{-2}) .
```

- The first two terms are **derived**.
- The third is **computed** to ten digits from the third-order matching integral, Ω₁ = 10.4233603788. It is
  **identified** as 2π²/3 + 8 ln²2 = 10.42336037874, with 3Ω₁ − 2π² − 24 ln²2 = 0 found by PSLQ at tolerance 10⁻⁹.
- **Update from #200:** against #195's values (growth.csv, regenerated in #200), the remainder of the three terms, times
  s² = (2(D − 3))², is 72.63, 73.70, 74.04 and 74.21 at D = 200, 500, 1000 and 2000. It tends to Ω₂ + Ω₃/s + …, as an
  O(1/D²) remainder should: at D = 1000, 74.04 = 74.37 − 0.33 − 0.003 (`logs_paper_numbers.txt`). The earlier "about
  73" and "66 ± 6" came from growth.csv's low `modes_r.py` rows.
- **Update from #195** ([w2_highprec.md](w2_highprec.md), W² in quad-double at 64 D from 50 to 10⁴; quad- against
  octo-double to at least 60 digits). Ω₁ = 10.42336037873851714322648088 by a fit, equal to 2π²/3 + 8 ln²2 to 8e-50.
  The fourth term is Ω₂ = 74.3693867279333175953847, not "about 73". It is **identified** as
  80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2 (PSLQ at 2.7e-44; the coordinator's independent extrapolation agrees to 1e-40).
  Both remain identifications, and the derivation is #25's and #26's.
  The `modes_r.py` values of the table below and of growth.csv (D ≥ 20) are low by 2.7e-9 (D = 30) to 6e-6
  (D = 2000); `w2_output.txt` has the 60-digit values.
  *Added by the coordinator:* the "66 ± 6" at D = 1000 above came from those low values.
- **Update from #26** ([omega2.md](omega2.md)). Ω₂ = 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2 is **derived**. G₄ = 0, written
  as one integral of closed forms, reduces by hand to the Clausen integral ∫₀^{π/2} tan θ Cl₂(π − 2θ) dθ = π³/12 and the
  moments ∫₀^{π/2} lnⁿcos θ dθ. Each step is checked in sympy, or in mpmath at 43 to 45 digits. An evaluation from the
  hierarchy's definitions, which uses none of the reduction, gives the same Ω₂ to 64 digits, and #195's fit to its 38
  printed digits. The same steps at order 3 give Ω₁ = 2π²/3 + 8 ln²2 with no Clausen term. Both reductions use the
  exchange identities E1 and E2, derived in omega2.md and checked numerically; their independent re-derivation is #209.
  In the form above, the next term is (20ζ(3) − (2/3)π² ln 2 − (8/3) ln³2)/(D − 3)². The remainder after it is
  Ω₃/(2(D − 3))³, with Ω₃ = −652.84010757006 a fit (#195), not the "about 73/(2(D − 3))²" above.
  *Added by the coordinator:* the review derived the order equations again from #195's ℓ-system and integrated them by
  its own Taylor stepper (raw `I26-mac-review`, ind26). Its Ω₀, Ω₁ and Ω₂ equal the closed forms to 7e-45, 2e-44 and
  4e-40.
- *Added by the coordinator:* each coefficient found so far is of uniform weight, counting ln 2 as weight 1, π² as 2
  and ζ(3) as 3, with products adding. Ω₀ = −4 ln 2 has weight 1, Ω₁ = 2π²/3 + 8 ln²2 weight 2, and Ω₂ weight 3. On the
  wider basis with all weights up to 3, #26's PSLQ gives zero coefficients to the lower weights at 1e-30. This is an
  observation about three terms, not a claim about the series. A pure weight-4 search for Ω₃ at 35 digits found
  nothing (#195), and #208 repeats it with more digits.
  **Update from #208:** Ω₃ breaks the pattern. Its closed form below has four terms of weight 4 and one, −384ζ(3), of
  weight 3. The pure weight-4 search fails at 112 digits too, and the search with all weights up to 4 finds the relation.
- **Update from #208** ([omega3.md](omega3.md)). The fifth term is **derived**:
  Ω₃ = 22π⁴/45 + (16/3)π² ln²2 + (32/3) ln⁴2 − 320ζ(3) ln 2 − 384ζ(3) = −652.8401075700561254352393037996473610804617….
  G₅ = 0 reduces by hand to Clausen integrals and log-cosine moments, using #26's E1 and E2 and one new identity, E3,
  which removes Y₄. Two evaluations from the hierarchy's definitions agree on Ω₃ to 112 digits, and #195's fit agrees to
  4e-38. In the form above, the next term is (Ω₃/8)/(D − 3)³ = −81.60501344625701568…/(D − 3)³. The remainder after it is
  Ω₄/(2(D − 3))⁴. **Update from #220** ([omega4.md](omega4.md)): the sixth term is **derived** too:
  Ω₄ = −10432ζ(5) − 160π²ζ(3) + 640ζ(3) ln²2 − (64/9)π² ln³2 − (88/45)π⁴ ln 2 − (128/15) ln⁵2 − (2432/45)π⁴
  + 1536ζ(3) ln 2 + 3072ζ(3) = −12794.49986002319463…, so (Ω₄/16)/(D − 3)⁴ = −799.6562412514…/(D − 3)⁴.
  It comes from G₆ = 0 by E1, E2, E3 and a new step E4. Two method-B grids agree to 209 digits, and #195's fit agrees
  to its 14. *Added by the coordinator:* the review's own Taylor hierarchy, sharing no code with method B, gives it
  to 1.1e-39. Ω₄ is not of uniform weight either: its terms have weights 5, 4 and 3. The remainder after six terms is
  Ω₅/(2(D − 3))⁵. **Update from #263** ([omega5.md](omega5.md)): the seventh term is **derived** as well:
  Ω₅ = (256/45)ln⁶2 + (64/9)π²ln⁴2 − (2560/3)ζ(3)ln³2 + (176/45)π⁴ln²2 − 3072ζ(3)ln²2 + 640π²ζ(3)ln 2 + (9728/45)π⁴ln 2
  − 12288ζ(3)ln 2 + 41728ζ(5)ln 2 − (38812/315)π⁶ + 3200ζ(3)² − (256/3)π²ζ(3) + 74240ζ(5) + (39424/45)π⁴ − 36864ζ(3)
  = 40860.4553545256873…, so (Ω₅/32)/(D − 3)⁵ = 1276.8892298289277…/(D − 3)⁵.
  - Method: G₇ = 0 by E1–E3, E4 one level up, and two new steps, E2W for L₄ with a weight and E5 for Y₄. Then finite-part
    integration by parts to Beta moments, Parseval and seven classical Euler sums, which are checked numerically at 34 digits and not
    proved here.
  - Checks: two method-B grids agree to 382.7 digits, and #195's W² values give it to 31 digits with nothing shared. The
    coordinator's evaluation of the closed form matches to 1e-106 (`I263-mac-review`).
  - **Update from #264** ([identity_review.md](identity_review.md)): E4, E2W and E5 are re-derived independently with every boundary term (*proved*). They hold on an independent
    hierarchy to 62–63 digits, and nothing differs from `omega4.md` or `omega5.md`. At Ω₄ and Ω₅, |G₆| = 1.7e-59 and |G₇| = 6.9e-57 on that hierarchy, the level of its cutoff. Hypotheses:
    G₂ = G₃ = 0 for E4 and E5's constant, and G₄ = 0 for E2W and E4 on Y₅. The finite-part chain and the seven Euler sums are checked numerically, not re-derived.
  - Ω₅'s terms have weights 6, 5, 4 and 3. No Li₄(½), Li₅(½), Li₆(½) and no weight-7 constant enters.
  - **All the ln 2 of Ω₀ … Ω₅ comes from one factor:** W²a² = s 2^{−4/s} F(1/s), with F(ε) = 1 + (2π²/3)ε² + 80ζ(3)ε³ + (22π⁴/45 − 384ζ(3))ε⁴
    + … free of ln 2 through ε⁶. That is checked symbolically for these terms, both here and in the coordinator's review; nothing is
    claimed beyond them.
  - The eighth coefficient is computed (#282, [omega6.md](omega6.md)): Ω₆ = 3672259.16259439…, from the order-8 matching integral on two grids that agree to 430.8 digits. Its closed form (weights 7 to 3; Ω₆/64 = 57379.0494155373514…) is identified by PSLQ at 390 digits, not derived. The remainder after eight terms is Ω₇/(2(D − 3))⁷ with Ω₇ ≈ 2.55 × 10⁷, an estimate. F is free of ln 2 through ε⁷: ε² … ε⁶ checked symbolically, ε⁷ identified in the full basis with ln 2 to ln⁷2. The "pure ζ(3) tower" −384, 3072, −36864 of f₄ … f₆ does not continue (f₇'s is 565248, not 589824).

**What it means.**
- **Faster in higher D.** The e-folding time is a/√(2(D − 3)) at large D. The exact values are 0.846 a (D = 4),
  0.579 a (D = 5) and 0.462 a (D = 6), falling to 0.072 a at D = 100.
- **Where the mode lives (derived).** Near the throat, in the shell of width about a/D where the geometry differs from
  flat space, the throat moves rigidly: δr is constant there, with δr′/δr = −1/a at its outer edge. Outside the shell
  space is nearly flat, and the mode is a Gaussian, δr ∝ (r/a) e^{−r²/a²}, to leading order. So at large D the
  instability is a flat-space mode with the throat acting as a boundary.
- **A prediction for threshold evolutions.** CG §1.3 reads the knife edge between collapse and expansion as type I
  critical behaviour, with the static wormhole as critical solution. Then near the threshold the time before the
  wormhole collapses or expands grows like T ≈ (1/W) ln(1/|p − p*|). The slope 1/W should be 0.579 a at D = 5 and 0.462 a
  at D = 6. Shinkai and Torii evolved D = 5 and 6 [ST17] but did not measure this slope. At D = 4 it is the 0.846 a that
  Shinkai and Hayward's 2002 runs already show (CG §1.3).

## Method

### 1. The master equation

Gauge δφ = 0, which is regular everywhere because φ′ never vanishes. The perturbations are ν = N e^{Wt}, λ = L e^{Wt}
and r = r₀ + Y e^{Wt} in ds² = −e^{2ν}dt² + e^{2λ}dl² + r²dΩ². With n = D − 2 and a = 1 (`derive.py`):
- the tl equation gives L = Y′/r₀′;
- the scalar equation gives N = L − nY/r₀ (an integration function of t is removed by reparametrizing t);
- the angular equation gives

```math
Y'' - \left(\frac{2r_0''}{r_0'} + \frac{n\,r_0'}{r_0}\right)Y' + \left(\frac{n\,r_0'^2}{r_0^2} + \frac{r_0''}{r_0}\right)Y = W^2Y .
```

The tt and ll equations then hold identically (`master.py check`, symbolic). Its control, a wrong coefficient of Y′,
leaves nonzero residuals. The throat (r₀′ = 0) is a regular singular point with exponents 0 and 3. The mode that moves
the throat is the even one, with Y(0) = 1.

### 2. Numerics

Two independent integrations:
- in proper distance l (`master.py modes`, for D ≤ 10);
- in the areal radius r (`modes_r.py`, for D ≥ 20). It uses a Frobenius start at the throat and matches to the exact
  far-field solution r^{(D−1)/2} K_{(D−3)/2}(W r).

Both find exactly one growing mode, and they agree with [TS13]'s Table I. TS13 used a different master variable, a
gauge-invariant combination regularized at the throat.

| D | W²a² here (D ≤ 10: the l-integration; D ≥ 20: #195's crate, by #200) | [TS13] | 2(D − 3) − 4 ln 2 + Ω₁/(2(D − 3)) |
|---|---|---|---|
| 4 | 1.39708 | 1.39705243 | (not asymptotic) |
| 10 | 12.0455 | 12.0442650 | 11.972 |
| 20 | 31.5751101284 | 31.5751101285 | 31.53398 |
| 50 | 91.3457759136 | 91.3457759137 | 91.33830 |
| 100 | 191.2830177304 | 191.2830177297 | 191.28114 |
| 200 | 391.2543343808 | | 391.25387 |
| 500 | 991.2379721479 | | 991.23790 |
| 1000 | 1991.2326572613 | | 1991.23264 |

- The l-integration is accurate to about 10⁻³ at D ≤ 10. At large D its start near the throat drifts, which is why the
  r-integration exists.
- The r-integration's values move by 4×10⁻¹⁰ to 1.5×10⁻⁶ from D = 20 to 1000 when its numerical parameters are changed
  (`logs_modes_r.txt`).
  **Update from #200:** that spread understates the r-integration's error. Its values were low by 8×10⁻¹⁰ (D = 20) to
  1.9×10⁻⁶ (D = 1000) against #195's crate (`w2_highprec.md`). The table above now holds the crate's values.
- A fit of the r-integration values in 1/D gives W²a² − 2(D − 2) = −4.77259 ± 10⁻⁵, against −2 − 4 ln 2 = −4.7725887.

### 3. Large D: matched asymptotics

**Variables.** Near the throat, R = (r/a)^s with s = 2(D − 3), ε = 1/s, ℓ = ln R and r = e^{εℓ}. Write
W² = s + Ω with Ω = Ω₀ + εΩ₁ + …. The master equation becomes exactly

```math
Y_R = \varepsilon\,\frac{L\,r}{R}, \qquad (R-1)L_R = (1+\varepsilon\Omega)\,rY + \frac{L}{2R} + \left(\tfrac12+\varepsilon\right)\frac{L(R-1)}{R} - \left(\tfrac12+\varepsilon\right)\frac{Y(R-1)}{Rr} - \frac{Y}{2rR},
```

where L = dY/dr. At each order in ε this gives (R − 1)L_k′ − L_k/2 = h_k(R) (`matched.py`).

**The condition at each order.** The solution analytic at the throat is

  L_k = −2h_k(1) + (R − 1)^{1/2} ∫₁^R (h_k − h_k(1))(R′ − 1)^{−3/2} dR′.

At large R it contains G_k R^{1/2}, with G_k = ∫₁^∞ (h_k − h_k(1))(R − 1)^{−3/2} dR. But R^{1/2} = (r/a)^{D−3} is the
far zone's *growing* solution. So the mode requires G_k = 0 at every order.

**The orders** (closed forms; u = √(R − 1)):
- **Order 0:** Y₀ = 1, L₀ = −1.
- **Order 1:**
  - h₁ = Ω₀ + ln R − 2(R − 1)/R.
  - G₁ = 0 identically (∫ ln R (R − 1)^{−3/2} = 2π = ∫ 2(R − 1)/R · (R − 1)^{−3/2}).
  - L₁ = −2Ω₀ − 2 ln R exactly, and Y₁ = −ln R.
- **Order 2:**
  - Y₂ = −2Ω₀ ln R − (3/2)ln²R.
  - h₂ = Ω₁ − 2Ω₀(R − 1)/R − Ω₀ ln R − 2 ln²R.
  - G₂ = −4πΩ₀ − 2J, with J = ∫₁^∞ ln²R (R − 1)^{−3/2} dR = 8∫₀^∞ ln(1 + u²)/(1 + u²) du = 8π ln 2.
  - So **Ω₀ = −4 ln 2.**
- **Order 3:**
  - L₂ = −2Ω₁ + 2Ω₀ ln R + 4 ln²R − 8Ω₀ u arctan u − 16u K(u), with K(u) = ∫₀^u ln(1 + t²)/(1 + t²) dt. Its R^{1/2}
    term vanishes exactly when Ω₀ = −4 ln 2, as it must.
  - G₃(Ω₁) is linear, with slope −4π (computed: −12.5663706143, against −4π = −12.5663706144). G₃ = 0 gives
    **Ω₁ = 10.4233603788** (`order3.py`; the two finer of three grids agree to 10⁻⁹).
- **Update from #26, orders 3 and 4** ([omega2.md](omega2.md), `order4.py`). Y₃ is closed in Cl₂ and Cl₃, and L₃ also
  needs Li₂ and Li₃, through ∫₀^θ ln²cos φ dφ. G₄ = 0 gives Ω₂ = 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2. The integral is
  evaluated two ways: from the closed forms after the exchange identities E1 and E2, and from the definitions on
  Chebyshev panels in √ℓ. With E2, the order-3 integral needs no K(u): it reduces to ∫ lnⁿcos θ dθ, n ≤ 2.

**Matching.** The decaying far solution r^{(n+1)/2}K_ν(Wr), expanded in the overlap 1 ≪ R ≪ e^s, gives Y′/Y =
−1 − (2Ω₀ + 3 ln R)/s + …. The near zone gives the same, so the O(1/s) matching holds automatically. The constraint
comes only from the absence of the growing term.

## Limits

- Radial (spherically symmetric) perturbations only.
- The massless symmetric member only. The asymmetric (two masses), charged and massive-potential families are not done.
- The third term's closed form is an identification, not a derivation.
  **Update from #195:** this now holds for the fourth term too (Ω₂, identified above). What changed is the precision.
  **Update from #26:** both Ω₁ and Ω₂ are now derived (Result). The derivation rests on the exchange identities E1
  and E2, which are derived in omega2.md and checked numerically, but not yet re-derived independently (#209). Ω₃ is a
  fit, not identified (#195; #208 computes it from the order-5 integral).
  **Update from #208:** Ω₃ is now derived too, resting on E1, E2 and a third identity, E3, which #209 should also
  re-derive. Ω₄ is derived (#220); Ω₅ is derived (#263; E2W and E5 re-derived independently in #264); Ω₆ is computed and its closed form identified, not derived (#282); Ω₇ is an estimate (≈ 2.55 × 10⁷), not computed.
  **Update from #209:** E1, E2, E3 and the two Clausen integrals are re-derived independently, without the notes
  ([exchange_review.md](exchange_review.md)), with every boundary term stated. No sign, constant or result differs.
  So Ω₁, Ω₂ and Ω₃ are **derived**, and #25 is closed. *Added by the coordinator:* the review checked E1 to E3 on a
  hierarchy that shares no code with method B, and showed that they need G_m = 0 (raw `I209-mac-review`).
- **Prior work** (searched only briefly, three web searches). [TS13] give the numbers and say the time scale "becomes
  shorter as n becomes large", but give no formula. No large-D analysis of these wormholes was found. This is not a
  full literature check, and nothing here claims novelty.

## Files

- `derive.py`: the first-order equations.
- `master.py`: the symbolic check and its control (`logs_check.txt`), and the l-integration (`logs_modes_lowD.txt`).
- `modes_r.py`: the r-integration with its convergence check (`logs_modes_r.txt`).
- `matched.py`: the h_k from the series (`logs_matched.txt`).
- `order3.py`: Ω₁ (`logs_order3.txt`).
- `w2_crate/`, `w2_highprec.py`, `w2_highprec.md`: W² to 60 digits at 64 D, and the fit of Ω₁ … Ω₄ (#195, `w2_output.txt`).
- `order4.py`, `omega2.md`: orders 3 and 4, Ω₂ from G₄ = 0, and its reduction by hand (#26, `order4_output.txt`).
- `order5.py`, `omega3.md`: order 5, Ω₃ from G₅ = 0, its identification and its reduction by hand (#208, `order5_output.txt`).
- `exchange_review.py`, `exchange_review.md`: E1, E2, E3 and the Clausen integrals re-derived independently (#209, `exchange_review_output.txt`).
- `omega_graph.py`: `omega_coeffs.svg`, the coefficients Ω₀ … Ω₅, the estimate of Ω₆, and the error of each truncation against #195's W².
- `profiles.py`: the mode profiles and the table of growth rates behind Paper III's figure (`profiles.csv`,
  `growth.csv`).
- `paper_numbers.py`: every number quoted in Paper III, with checks (`logs_paper_numbers.txt`).

Written up as Paper III, [`paper/3-wormhole-largeD/main.pdf`](../../paper/3-wormhole-largeD/main.pdf).

## References

- [TS13] T. Torii, H. Shinkai, "Wormholes in higher dimensional space-time: exact solutions and their linear stability
  analysis", PRD 88, 064027 (2013), arXiv:1309.2058. Table I read for this note (full-precision values quoted above).
- [GGS09a], [SH02], [ST17]: see CG's reference list.
