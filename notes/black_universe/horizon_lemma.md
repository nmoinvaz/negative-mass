# A Killing horizon at an extremum of the areal radius (#261)

Written by the coordinator. The algebra is checked in [horizon_lemma.py](horizon_lemma.py), output
[horizon_lemma_output.txt](horizon_lemma_output.txt). Status words as in AGENTS.md: the lemma is *proved* (by hand), and its
algebra is *checked symbolically*. No sentence here says it is new; #260 reads the literature.

## Setting

- A static spherically symmetric solution in D = n + 2 dimensions, n ≥ 2, of Einstein gravity with a scalar field φ (kinetic sign ε,
  any potential). It has a bifurcate Killing horizon of surface gravity κ > 0.
- In Kruskal coordinates the metric is −2e^{−F} dV dU + r² dΩₙ². The Killing field is κ(V∂_V − U∂_U), so r, F and φ depend on (V, U)
  only through Π = −VU. The future horizon H⁺ is U = 0, and in the exterior V = e^{κv}, with v the advanced time.
- Write r = R(ρ(Π)) with ρ the radial coordinate of the static form −A dt² + dρ²/A + R² dΩₙ². The horizon is a simple zero ρ_h of A,
  κ = A′(ρ_h)/2, and ρ′(0) = dρ/dΠ at Π = 0 is finite and nonzero.

## Lemma

For every linear perturbation (δr, δF, δφ) that solves the linearized field equations:

1. On H⁺, δr_VV = 0, so δr = c₁ + c₂V along the whole line U = 0.
2. Under a linear diffeomorphism ξ, ∂_V(L_ξ r) on U = 0 equals −R′(ρ_h) ρ′(0) ∂_V(V ξᵁ). So c₂ (and c₁) are the same in every gauge
   if R′(ρ_h) = 0. If R′(ρ_h) ≠ 0, the shift U → U + η(U) changes c₂ by −η(0) R′(ρ_h) ρ′(0), and c₂ is pure gauge.
3. If R′(ρ_h) = 0, the horizon's cross-sections change in area by Ωₙ n r_hⁿ⁻¹ (c₁ + c₂ e^{κv}) to first order. With c₂ ≠ 0 the
   change grows at the surface-gravity rate.
4. If A = f(R(ρ)) with f differentiable, then R′(ρ_h) = 0 forces A′(ρ_h) = 0, that is κ = 0. A non-extremal horizon at an extremum of
   the areal radius needs a metric function A that is not a function of the areal radius alone.

*Added from #301* ([static_c1_members.md](static_c1_members.md)): by item 2, L_ξ r = −V ξᵁ R′(ρ_h) ρ′(0) on U = 0 vanishes at V = 0, so
c₁ = δr(S) is the same in every gauge at every member, also where R′(ρ_h) ≠ 0 and c₂ is pure gauge; δA_S = 8π r_S c₁ with
r_S = R(ρ_h). #301 checks it on exact pure gauge at c = −1, −1.5, −0.3 and on Schwarzschild.

## Proof

1. For this metric R_VV = −(n/r)(r_VV + F_V r_V) (L0, by a brute-force Ricci tensor at n = 2 and 3). Since g_VV = 0, the VV equation is
   the constraint C_VV ≡ r_VV + F_V r_V + (r/n) ε φ_V² = 0.

   A background function g(Π) has ∂_V g = −U g′ and ∂_V² g = U² g″, both zero on U = 0. Linearize:

   δC_VV = δr_VV + F_V δr_V + r_V δF_V + (ε/n) φ_V² δr + (2r/n) ε φ_V δφ_V.

   On U = 0 every term but the first has a factor r_V, F_V or φ_V, so δC_VV = δr_VV there (L1). Then δr_VV = 0 on U = 0, which
   integrates to δr = c₁ + c₂V.
2. L_ξ r = ξⱽ r_V + ξᵁ r_U. On U = 0, r_V = r_VV = 0 and r_U = −V R′(ρ_h) ρ′(0). The term ξⱽ_V r_V vanishes too. So
   ∂_V(L_ξ r) = ∂_V(ξᵁ r_U) = −R′(ρ_h) ρ′(0) ∂_V(V ξᵁ) (L2), and similarly L_ξ r = −V ξᵁ R′(ρ_h) ρ′(0).

   If R′(ρ_h) = 0, both vanish for every ξ, so c₁ and c₂ are gauge invariant. If R′(ρ_h) ≠ 0, ξ = η(U)∂_U keeps the double-null form
   and shifts c₂ by −η(0) R′(ρ_h) ρ′(0), which can cancel it.
3. In the exterior V = e^{κv}, so δr = c₁ + c₂e^{κv} on H⁺. The perturbed horizon lies at U = O(δ). Since r_U = 0 on U = 0 when
   R′(ρ_h) = 0, r there equals r on U = 0 to first order. The area Ωₙ rⁿ then changes by Ωₙ n r_hⁿ⁻¹ δr.
4. A′(ρ_h) = f′(R(ρ_h)) R′(ρ_h) = 0 (L3), and κ = A′(ρ_h)/2. □

The matter enters only through the linearized T_VV on U = 0. That is εφ_V(2δφ_V) for the scalar, and it vanishes there. The proof holds
unchanged for any matter whose linearized T_VV vanishes on the background horizon, for example a radial Maxwell field (T_VV = 0 in
spherical symmetry).

## The members

| spacetime | R′(ρ_h) | κ | the lemma's verdict |
|---|---|---|---|
| c = −1 black universe (b = 1) | 0 | 2/π | c₂ gauge invariant, growing as e^{κv} on H⁺ (#248 finds it free in regular data) |
| c = −1.5 black universe | 0.321578 | 0.551682 | c₂ pure gauge; #248 uses this as its planted case (0.3216) |
| c = −0.3 black universe | −0.679735 | 0.601300 | c₂ pure gauge; c₁ invariant (#301) |
| Schwarzschild (2M = 1) | 1 | 1/2 | c₂ pure gauge |
| Simpson–Visser, a = 2m | 0 | 0 | extremal: outside the lemma (item 4) |

All 20 of #212's black universes have A not a function of R alone. Only c = −1 puts the horizon at the throat.

## Checks (`horizon_lemma.py`, 5 of 5)

| check | what | mutation that fails it |
|---|---|---|
| L0 | R_VV of the warped double-null metric at n = 2, 3 | `wrong-ricci` (n + 1 for n) |
| L1 | the linearized C_VV on U = 0 is δr_VV, for any n, ε and background functions of Π | `nonstatic` (r gains V² w(U)): a term 2V w(0) δF_V remains |
| L2 | ∂_V(L_ξ r) on U = 0 for an arbitrary ξ | `l2-short`: the form without ξᵁ_V r_U. This was the check's first draft, which failed on its first run because it omits that term. The lemma's conclusion was unaffected. |
| L3 | A = f(R(ρ)) at an extremum of R has A′(ρ_h) = 0 | `a-not-fr` (A + (ρ − ρ_h)) |
| L4 | the members' R′(ρ_h) and κ at 30 digits | `rp-nonzero` |

## What this does not show

- That c₂ is free in Cauchy data. For the c = −1 member that is #248's result, numerical evidence.
- Growth off the horizon at fixed ρ. The lemma is about H⁺. #248 shows that the c = −1 mode Y = c₂VΨ(Π) is smooth across S and grows at
  fixed ρ.
- Anything nonlinear. Nor whether c₂ ≠ 0 ends in collapse or expansion (#250).
- Which other known solutions have a non-extremal horizon at an area extremum. A search of the literature is part of #260.
