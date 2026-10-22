# Phase 4: mass bounds with a signed source T = T₊ − T₋ (gap 6)

**Summary.** Witten's identity gives a bound on each slice: the mass is at least minus the phantom
sector's energy on that slice, weighted by |ε|². The bound is exact up to the gravitational gradient
term. On the Shirokov phantom star it is tight to 0.7 %. But it is **not a bound in terms of E₋ alone**,
and it is **not uniform in time**:
- The weight |ε|² is unbounded on compact phantom configurations. There are DEC-satisfying data with
  E₋ → 0 and M → −∞.
- E₋ is not conserved. The Bondi runaway pair is a configuration in which E₋ (and E₊) grow without
  bound at fixed total mass.

So as a *dynamical* statement ("a radiating system cannot lose more than X") the Witten argument
reduces to something close to trivial. It bounds the loss at time u by the phantom content at time u,
and nothing bounds that content. Unbounded loss is neither excluded nor shown.

Units G = c = 1. Scripts: `witten_bound.py`, `sobolev_bound.py` (outputs `*_output.txt`).

## 1. Witten identity with a signed source

Setup. (Σ, g, K) is a complete asymptotically flat initial data set (one end, g − δ = O(r^−q), q > ½,
K = O(r^−q−1), μ, J ∈ L¹). It is spin automatically in 3D. The source is μ = μ₊ − μ₋ and J = J₊ − J₋,
where (μ±, J±) are the normal projections of T±. **DEC sector by sector:** μ± ≥ |J±|. ε₀ is a constant
spinor at infinity with |ε₀| = 1, aligned with P so that the boundary term is E − |P|. ε solves the
hypersurface Dirac equation D̃ε = 0 with ε → ε₀.

The Witten / Parker–Taubes identity (checked numerically below to 1.6×10⁻⁷, including normalisation) is

  E − |P| = (1/4π) ∫_Σ |∇̃ε|² dV + ∫_Σ T(ξ_ε, n) dV,  ξ_ε^a = ⟨ε, γ⁰γ^a ε⟩ (future causal, ξ⁰ = |ε|²).

Each T± satisfies the DEC, so T₊(ξ_ε, n) ≥ 0 and 0 ≤ T₋(ξ_ε, n) = μ₋|ε|² + J₋·ξ⃗ ≤ (μ₋ + |J₋|)|ε|².

**Lemma A (per slice).** If the Witten spinor exists, then
  E − |P| ≥ −∫_Σ T₋(ξ_ε, n) dV ≥ −∫_Σ (μ₋ + |J₋|) |ε|² dV ≥ −2 ∫_Σ μ₋ |ε|² dV.
The first inequality is an identity up to the two non-negative terms (1/4π)∫|∇̃ε|² and ∫T₊(ξ_ε, n).

*Existence caveat.* Parker–Taubes prove existence from the DEC of the *total* T (coercivity). With a
signed source, existence is an extra hypothesis: D̃ must have no decaying kernel. In the time-symmetric
case, existence follows from the κ < 1 condition of Lemma B.

*Normalisation caveat.* ξ is not normalised. |ε|² → 1 at infinity, but in the interior it can be
anywhere in (0, ∞). Nothing gives |ε|² ≤ 1 or |ε|² ≥ 1 in general. The examples below have |ε|² ≫ 1
inside compact phantom matter.

**Lemma A′ (conformally flat time-symmetric data, exact).** For g = ψ⁴δ and K = 0, ε = ψ⁻²ε₀ is an
explicit Witten spinor, by conformal covariance of the Dirac operator in n = 3. Then
  M = (3/2π) ∫ |∂ψ|² ψ⁻⁴ d³x + ∫ μ₊ ψ⁻⁴ dV − ∫ μ₋ ψ⁻⁴ dV,
so M ≥ −W₋, with W₋ := ∫ μ₋ ψ⁻⁴ dV. No existence assumption is needed, and all spherical data are
covered. The weight ψ⁻⁴ = (r_iso / r_areal)² is the local "redshift" factor.

**Lemma B (a priori bound, time-symmetric).** Take K = 0, so J = 0 and ε is a harmonic spinor. Let
S_g be the Sobolev constant, ∫|∇h|² ≥ S_g ‖h‖²_{L⁶(g)}. Let a = ‖μ₋‖_{L^{3/2}(g)} and E₋ = ∫ μ₋ dV_g. If
  κ := 8πa / (3 S_g) < 1,
then the Witten spinor exists and
  **M ≥ −E₋ / (1 − κ).**

*Proof sketch.* Put u = |ε| and v = u − 1 ∈ L⁶. The refined Kato inequality for harmonic spinors in 3D
gives |∇ε|² ≥ (3/2)|∇u|². The identity gives 4πM ≥ (3/2)∫|∇v|² − 4π∫μ₋u². For any t > 0,
u² ≤ (1 + 1/t) + (1 + t)v². Hölder and Sobolev then give ∫μ₋v² ≤ (a/S_g)∫|∇v|². Choose 1 + t = 1/κ: the
gradient terms cancel and 1 + 1/t = 1/(1 − κ). The same estimate with ε₀ = 0 shows that D has no
decaying kernel, which gives existence (Fredholm, index 0).

For g = ψ⁴δ, S_g ≥ 3(π/2)^{4/3} (inf ψ / sup ψ)². With the plain Kato inequality the constant is
4πa/S_g instead. This is the spinor analogue of Shi–Tam's L^{3/2} bound (§5), and it carries the same
weakness: it needs the phantom sector to be dilute in L^{3/2}. Once κ ≥ 1 it says nothing, and §3(b)
shows that this is unavoidable.

## 2. Bondi version

On an asymptotically null slice Σ_u ending on the cut S_u of 𝒥⁺, the Horowitz–Perry /
Ludvigsen–Vickers / Reula–Tod identity has the same structure, with T(ξ_ε, n) in the bulk. So, given a
Witten spinor on Σ_u,

  **M_B(u) ≥ −∫_{Σ_u} T₋(ξ_ε, n) dV =: −W₋(u)**, and in the time-symmetric-like regime of Lemma B,
  M_B(u) ≥ −E₋(u)/(1 − κ(u)).

Mass loss. Shirokov's fields are massive (m = 1), so no matter flux reaches 𝒥⁺ and
dM_B/du = −(1/4π)∮|N|² ≤ 0. The decrease is monotone and purely gravitational, whatever the sign of the
matter. (For a *massless* ghost the phantom's own flux would *raise* M_B: then
dM_B/du = −(1/4π)∮(|N|² + |∂_u(rφ₊)|² − |∂_u(rφ₋)|²), which is no longer monotone.)

Combining the two: the energy radiated between u₀ and u satisfies M_B(u₀) − M_B(u) ≤ M_B(u₀) + W₋(u).
That is a genuine bound **only if W₋(u) is bounded**. It is not:

- **E₋ is not conserved.** ∇_a T₋^{ab} = 0 holds sector by sector, but without a Killing field
  ∫T₋(∂_t, n) changes through the gravitational coupling.
- **The runaway pair is the explicit "pair production at fixed M".** Bondi's pair, Shirokov's
  simulation, and our phase-3 1PN integrations all show m₊ and m₋ = −m₊ accelerating together at fixed
  total mass ≈ 0. The sector energies grow like E± ≈ γ m, and the Witten weight on the phantom grows like
  ∫(μ₋ + |J₋|) ≈ γ(1 + v)E₋⁽⁰⁾. So W₋(u) → ∞ along the runaway (rough numbers for the Shirokov star are
  in `witten_bound_output.txt`: 0.016 at v = 0.1, 0.205 at v = 0.99).
- **The constraints allow unbounded E± at M = 0.** Take T₊ = T₋ pointwise, i.e. two fields with
  identical data. Then T = 0, the metric is flat for all time, M = 0, and E₊ = E₋ is arbitrary. This
  case is stable: each field evolves freely, and E± is constant.

**Conclusion for task 2.** The proposed claim, "M_B(u) ≥ −E₋(u) and E₋ is conserved, so extraction is
bounded", is **refuted** on both halves. The best correct statement is
M_B(u) ≥ −∫_{Σ_u} T₋(ξ_ε, n). It bounds the loss at time u by the phantom sector's weighted energy *at
time u*, and in the runaway that energy grows without bound. Whether M_B can actually go to −∞ depends on
the gravitational-wave luminosity of the runaway pair. That luminosity has not been computed (relativistic
regime, Bičák-type boost-symmetric radiation), and it is the open question this phase leaves.

## 3. Numerical checks (`witten_bound.py`, `sobolev_bound.py`)

**(a) Shirokov stars.** Spherical data use Lemma A′ (areal → isotropic, ψ⁻⁴ = (r_iso/r)²).

| star | M_ADM | E = ∫ρ dV | W = ∫ρψ⁻⁴ dV | G = (1/4π)∫\|∇ε\|² | M − (G ± W) |
|---|---|---|---|---|---|
| canonical ω = 0.750 | +0.0146260 | 0.0146623 | 0.0145169 | 1.09e−4 | +1.6e−7 |
| phantom ω = 0.7603 | −0.0145805 | 0.0145445 | 0.0146883 | 1.08e−4 | −1.7e−7 |

- The naive bound M ≥ −E₋ is **violated**: |M|/E₋ = 1.0025. The phantom's negative "binding energy"
  makes |M| > E₋. The canonical star is the mirror case, M/E₊ = 0.9975.
- The Witten bound −W₋ = −0.014688 holds with slack G = 0.74 % of W₋.
- The bound sup ψ⁻⁴ · E₋ = 1.0133 E₋ gives −0.014738.
- The Lemma B bound −E₋/(1 − κ) gives −0.014597, with κ = 3.6×10⁻³. It is tight to **0.11 %**. At
  ω = 0.905: −0.008392 vs M = −0.008387.

**Mixed pair.** M_tot = +4.5×10⁻⁵. While the phantom star stays slow and dilute, the bound is
M_B ≥ −W₋ ≈ −0.0147. That is about 320 × |M_tot|, so the bound allows the total mass to fall far below
zero even before any growth of W₋. During the runaway the floor itself moves down like γ(1 + v).

**(b) Phantom shells: no bound M ≥ −C·E₋.** These are conformally flat, time-symmetric shells with
μ₋ ≥ 0 and J = 0, so the DEC holds sector-wise. Analytically, for a thin shell at isotropic radius R:
ψ_in = 1 − |M|/2R, E₋ = |M|ψ_in, W₋ = |M|ψ_in⁻³. As |M| → 2R the proper energy E₋ → 0 while M = −2R.
Taking R large then gives **M → −∞ at E₋ → 0** (the areal radius of the shell → 0). Numerically, with
Gaussian shells of width 0.05 at R = 1:

| \|M\| | ψ_min | E₋ | \|M\|/E₋ | W₋ | identity residual |
|---|---|---|---|---|---|
| 0.1 | 0.950 | 0.0952 | 1.05 | 0.116 | 1e−12 |
| 1.0 | 0.501 | 0.515 | 1.94 | 7.33 | 1e−8 |
| 1.8 | 0.102 | 0.229 | 7.85 | 1006 | 2e−5 |
| 1.95 | 0.027 | 0.107 | 18.3 | 27582 | 1e−3 (resolution) |

The Witten bound −W₋ holds throughout but becomes useless (W₋ ≫ |M|). Lemma B holds where κ < 1
(|M| = 0.1: −0.110 vs −0.100; |M| = 0.5, width 0.3: −1.10 vs −0.50) and is void beyond that. Caveat: these
shells are generic DEC matter. I did not build a scalar-field realisation, although a static profile with
V ≥ 0 should be able to supply any smooth μ₋.

**(c) Mixed shells.** The bound holds in all 5 cases tested, including M_tot = −1.5 with E₋ = 0.32.

## 4. What is proven / what is open

Proven (modulo standard Witten-spinor analysis, which is only sketched here):
- Lemma A: E − |P| ≥ −∫T₋(ξ_ε, n), under DEC per sector and existence of the spinor.
- Lemma A′: an exact identity for conformally flat time-symmetric data.
- Lemma B: M ≥ −E₋/(1 − κ) for time-symmetric data with κ < 1, with existence included.
- The Bondi analogue of Lemma A on each cut. For massive fields, M_B is monotone.

Refuted:
- Any bound M ≥ −C E₋ with universal C.
- Conservation of E₋.
- Hence any uniform-in-time floor for M_B from this method.

Open:
1. Does the runaway pair radiate enough for M_B → −∞? This needs the GW luminosity of a relativistically
   accelerating ± pair.
2. A non-time-symmetric Lemma B. Kato fails for ∇̃ when K ≠ 0, so |K| control is needed.
3. Whether a sharper weight exists. One candidate is a different spinor or a sector-adapted ε₀ that
   captures T₊'s positivity against T₋.
4. Existence of the Witten spinor without κ < 1.

## 5. Literature status (searched Oct 2026)

No paper found that states a positive-energy / Witten bound for a signed sum of two DEC tensors, or for
ghost / phantom matter. The inequality of Lemma A is a one-line corollary of Witten's identity and would
be regarded as folklore. **Do not claim novelty.** Related work:

- **Shi & Tam 2004**, math/0406559, "Some lower estimates of ADM mass and Brown–York mass", Thm 3.1. It
  gives a lower bound on ADM mass without R ≥ 0, in terms of ∫R₊, ∫R₋ and ‖R₋‖_{L^{n/2}} vs the Sobolev
  constant. It needs ‖R₋/4‖_{n/2} < Λ/2. This is the Riemannian (time-symmetric) analogue of Lemma B and
  has the same structure. Our Lemma B is essentially a spinor re-derivation.
- **Cecchini, Lesourd & Zeidler 2023**, arXiv:2307.05277. Positive mass for spin initial data with the
  DEC violated inside, provided a "dominant energy shield" (μ − |J| ≥ σ > 0 on a collar) is present.
  Witten argument with a modified connection.
- **de Aguiar Alves, Landulfo & Arderucio Costa 2024**, arXiv:2408.00154, "Positive mass in GR without
  energy conditions". Uses Penrose–Sorkin–Woolgar, ANEC and the dynamical instability of
  barotropic negative-mass stars. Qualitative only, with no bound in terms of the negative energy.
- **Deffayet, Held, Mukohyama & Vikman 2023**, arXiv:2305.09631. Classical ghost + positive-energy
  systems with bounded motion for all data (integrable classes); runaway in other classes. No gravity.
- Tafel 2015, arXiv:1505.07775: a Trautman–Bondi lower bound that assumes the DEC.
- Cline, Jeon & Moore 2004 (quantum ghost decay via gravitons), already in the notes.

Not found: any analysis of whether gravitationally coupled classical ghost matter can drive the Bondi mass
to −∞.
