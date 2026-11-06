# Exactly one unstable radial mode, in every dimension

November 2026. A proof that the massless, reflection-symmetric Ellis wormhole has exactly one unstable radial mode in
every D ≥ 4. It is Paper III's premise for the type I reading (the wormhole is a critical solution with one unstable
mode), which until now rested on Torii and Shinkai's numerics in D = 4–11, 20, 50, 100 [TS13] and, at D = 4, on the
arguments of [GGS09a] and [CPS20]. Checks: `one_mode.py` (log `one_mode.log`).

Status words: **proved** (by hand, here), **checked symbolically** (sympy, with a control that fails), **checked
numerically** (mpmath, with a control), **assumed** (a standard fact or an earlier result, named).

## Statement

**Theorem.** Let n = D − 2 ≥ 2 and consider the static wormhole of Paper III with throat radius 1,

```math
ds^2 = -dt^2 + dl^2 + r(l)^2\,d\Omega_n^2, \qquad r'^2 = 1 - r^{-2(n-1)}, \quad r(0) = 1 ,
```

held open by a massless ghost scalar. Its spherically symmetric linear perturbations of the form e^{Wt} × (a profile
that is smooth and decays exponentially as |l| → ∞) exist for exactly one W > 0. The profile of that mode is even in l.
There is none that is odd in l (antisymmetric between the two sides).

The proof works for real n > 1 as well; only n ≥ 2 is claimed.

**What is assumed.**
1. Paper III's master equation for these perturbations in the gauge δφ = 0 (regular because φ′ never vanishes), with
   Y = δr:

   ```math
   Y'' - \Big(\frac{2r''}{r'} + \frac{n\,r'}{r}\Big)Y' + \Big(\frac{n\,r'^2}{r^2} + \frac{r''}{r}\Big)Y = W^2\,Y ,
   ```

   and that the remaining Einstein equations then hold identically (Paper III, `master.py check`, checked symbolically
   with n symbolic, with a control).
2. Standard facts on H = −d²/dl² + U for U bounded, continuous and tending to 0 at ±∞ [Wei87; BS91]: H is self-adjoint
   on H²(ℝ) and σ_ess(H) = [0, ∞) (Weyl); below 0 the spectrum consists of eigenvalues, whose eigenfunctions are C² and
   decay exponentially together with their derivative; the min-max principle.

**What is not claimed.** Non-radial perturbations, the asymmetric (two-mass), charged and massive-potential families,
and nonlinear stability are all outside this theorem.

## Proof

**Step 1: a regular master equation (checked symbolically, S1–S4).** Put s = r′ r^{(n−2)/2} and w = s².
- With Y = rZ the Z term cancels for any r(l) (S1):

  ```math
  Z'' - \frac{w'}{w}\,Z' = W^2 Z .
  ```

- With q = Z′/w the equation reads q′ = W² Z/w, so Z = w q′/W² and, differentiating, (w q′)′ = W² w q. Conversely a
  solution q of that equation gives Z = w q′/W² with Z′ = w q, which solves the Z equation (W ≠ 0; S2).
- With q = χ/s, it becomes a Schrödinger equation (S2):

  ```math
  -\chi'' + U\chi = -W^2\chi, \qquad U = \frac{s''}{s} = \frac{(n-2)(n-4)}{4r^2} - \frac{3n^2}{4r^{2n}} \quad\text{(S3).}
  ```

U depends on l only through r(l), so it is even, smooth and bounded, with U(0) = −(n + 4)(n − 1)/2 (S4). It tends to 0
like (n − 2)(n − 4)/(4l²). At D = 4 it is −3/r⁴, the potential of [CPS20] (their Eq. (71) with b = 1), and χ = (Y/r)′/r′
is their master variable (A − C)/r in the gauge δφ = 0.

**Step 2: modes are eigenvalues (proved).** Exponential decay passes both ways through the chain Y = rZ, Z = w q′/W²,
χ = s q:
- Given a mode, χ = Z′/s is smooth and decays exponentially, so −W² is an eigenvalue of H. Near the throat
  s ≈ (n − 1) l. For an even profile Z′ is odd and O(l). For an odd profile Y = O(l³), since Paper III's throat
  exponents are 0 and 3, so Z′ = O(l²). In both cases χ is smooth at the throat.
- Given an eigenvalue −W² < 0 with eigenfunction χ, the profile Z = (sχ′ − s′χ)/W² and Y = rZ are smooth and decay
  exponentially (s grows only like a power of l), and they solve Paper III's equation.

The map is one-to-one for W ≠ 0: its kernel is Y ∝ r, which solves the Y equation only at W = 0. Parity is kept, since
s is odd and r is even: even Y gives even χ. So the theorem says that H has exactly one negative eigenvalue, and that
it is even.

**Step 3: the two zero-energy solutions (proved; S5, N1–N3).**
- **Odd.** χ₀ = s solves Hχ = 0 (U = s″/s). It vanishes only at l = 0, since r′ = 0 only at the throat, and it is
  positive for l > 0.
- **Even.** χ_N = sH, with H′ = −1/s² (S5). H is fixed by H(l) = 1/(m²l) − ∫₀^l (1/s² − 1/(m²y²)) dy, with m = n − 1.
  The integrand is even and bounded, so H is odd and χ_N is even, with χ_N(0) = 1/m > 0.
  - On (0, ∞), χ_N vanishes exactly where H does. H is strictly decreasing from +∞, so χ_N has at most one positive
    zero, and it has one exactly when H(∞) < 0.
  - For n = 2 or 3, 1/s² ~ l^{2−n} is not integrable at ∞, so H → −∞.
  - For n ≥ 4, the substitution t = r^{−2m} turns ∫dl/s² into (1/2m)∫t^{α−1}(1 − t)^{−3/2}dt with α = 1/2 − 1/m. Near
    the throat 1 − t = m²l² + O(l⁴), and the divergence 2(1 − t)^{−1/2}/(2m) is exactly the subtracted 1/(m²l). So
    H(∞) = −B(α, −1/2)/(2m), with B the analytically continued Beta function.
  - B(α, −1/2) = Γ(α)Γ(−1/2)/Γ(−1/m). For m ≥ 3, Γ(α) > 0 (since 0 < α < 1/2), Γ(−1/2) = −2√π < 0, and Γ(−1/m) < 0
    (Γ is negative on (−1, 0)). So B > 0 and H(∞) < 0.
- **Numerical checks.** The formula for H(∞) is checked against the integrated ODE at m = 3, 4, 6, 10, 25 (N1, with
  the sign flipped as control). B > 0 is checked for m = 3–200 (N2). The zero counts, one even and none odd, are checked
  by integration at D = 4, 5, 6, 8, 12, 20, 50 (N3).

So χ_N has exactly one zero on (0, ∞) for every n ≥ 2, and χ₀ none.

**Step 4: a comparison lemma (proved).** Let −u″ + Uu = Eu with E < 0, and −v″ + Uv = 0. Let u > 0 on an interval
(α, β), whose ends are zeros of u, the point 0 with u′(0) = v′(0) = 0 (Neumann), the point 0 with u(0) = v(0) = 0, or
±∞. Then v changes sign in (α, β).
- Proof. Suppose v ≥ 0 on (α, β) (otherwise replace v by −v), and let 𝒲 = u′v − uv′.
- Then 𝒲′ = −E uv, so 𝒲(β) − 𝒲(α) = −E∫uv > 0.
- At a zero β of u, 𝒲(β) = u′(β)v(β) ≤ 0; at a zero α, 𝒲(α) = u′(α)v(α) ≥ 0.
- At a Neumann end, or an end where u = v = 0, 𝒲 = 0.
- At ±∞, 𝒲 → 0: u and u′ decay exponentially (assumption 2), while χ₀, χ_N and their derivatives grow at most like
  a power of l (s ~ l^{(n−2)/2}, and |H| grows at most like l).
- Every case gives 𝒲(β) − 𝒲(α) ≤ 0, a contradiction. ∎

**Step 5: at most one, and none odd (proved).**
- *Odd sector.* An odd eigenfunction u with E < 0 has u(0) = 0. On its first nodal interval (0, β) the lemma with
  v = χ₀ (χ₀(0) = 0) forces a zero of χ₀ in (0, β), but χ₀ > 0 there. So there is no odd negative eigenvalue.
- *Even sector.* Suppose there were two even eigenvalues E₁ < E₂ < 0. Their eigenfunctions are orthogonal, so at least
  one of them, u, changes sign on (0, ∞).
  - Its nodal intervals on [0, ∞) include [0, x_first), with u′(0) = 0, and (x_last, ∞).
  - The lemma with v = χ_N gives a zero of χ_N in each of them, so two positive zeros of χ_N. That contradicts Step 3.
  - So there is at most one even negative eigenvalue.

**Step 6: at least one (proved).** Let x₁ be the positive zero of χ_N, and let ψ = χ_N on [−x₁, x₁] and 0 outside.
- ψ is even and in H¹(ℝ). Integrating by parts, with ψ(±x₁) = 0, its quadratic form is ∫(ψ′² + Uψ²) = 0.
- If H had no negative spectrum in the even sector, ψ would minimize the form there. The form vanishes at ψ, so the
  first variation q(ψ, φ) vanishes for every even φ, and for every odd φ by parity.
- Then ψ is a weak, hence classical, solution of Hψ = 0 on ℝ. That is impossible, since ψ = χ_N on (−x₁, x₁) but
  ψ = 0 beyond x₁, which uniqueness forbids.
- So the even sector has spectrum below 0, which is an eigenvalue because σ_ess = [0, ∞).

Steps 5 and 6 give exactly one negative eigenvalue, and it is even. By Step 2 that is exactly one unstable radial mode,
even in l. ∎

**Cross-check (N4).** Shooting the even bound state of −χ″ + Uχ = Eχ gives −E = W² of Paper III's table at
D = 4, 10, 20 and 50 to 10⁻⁶. The control, U with 3n² → 3n² + 1, moves it by more than 10⁻³.

This shares Paper III's master equation but not its integration codes, master variable or end conditions.

## Relation to the literature (read for this note)

- **[GGS09a]** (D = 4). The abstract states that all static ghost-scalar wormholes have "precisely one unstable,
  exponentially in time growing mode". The body proves it for γ₁ ≠ 0, by a nodal count of the zero mode 1/Ψ₀ of a
  factorized operator. For the symmetric member (γ₁ = 0) that zero mode, x/√(1 + x²), is not normalizable. A
  Rayleigh–Ritz bound shows a negative eigenvalue exists, and uniqueness is asserted "by continuity". (Read through an
  ar5iv page via a summarizing fetch.)
- **[CPS20]** (D = 4). The gauge-invariant master equation used in Step 1, with potential −3b²/r⁴ in the symmetric case.
  - Sec. VI.1 applies "the Sturm oscillation theorem" to the zero-energy solution (x − bγ₁)/(rF), which has one zero,
    "for each γ₁ including γ₁ = 0".
  - For γ₁ = 0 that solution is not normalizable, a threshold case, and Sec. VII.1 states the unique negative eigenvalue
    citing general facts [BS91] and [GGS09a]'s analysis, without deriving it there.
  - Read through ar5iv via a summarizing fetch. Their cited oscillation theorems are [Sim05, Theorem 3.4] and
    [Wei87, Theorem 14.8], used for their AdS and dS cases.
- **[TS13]** find one negative mode numerically in D = 4–11, 20, 50, 100.
- **Found by two searches:** nothing proving the count in D > 4.
- **What differs here:** the D-dimensional potential and the treatment of the two parity sectors separately. The
  threshold case is handled by the even zero-energy solution χ_N, whose positive zero is counted exactly (Step 3), so no
  oscillation theorem at the edge of the continuum is needed. This also covers D = 4, γ₁ = 0.

## References

- [GGS09a] J. A. González, F. S. Guzmán, O. Sarbach, "Instability of wormholes supported by a ghost scalar field. I.
  Linear stability analysis", CQG 26, 015010 (2009), arXiv:0806.0608 (journal as cited in [CPS20]).
- [CPS20] F. Cremona, L. Pizzocchero, O. Sarbach, "Gauge-invariant spherical linear perturbations of wormholes in
  Einstein gravity minimally coupled to a self-interacting phantom scalar field", PRD 101, 104061 (2020),
  arXiv:1911.13103 (journal from the DOI on the arXiv page).
- [TS13] T. Torii, H. Shinkai, PRD 88, 064027 (2013), arXiv:1309.2058 (see `README.md`).
- [Wei87] J. Weidmann, *Spectral Theory of Ordinary Differential Operators*, Lecture Notes in Math. 1258, Springer,
  Berlin, 1987 (as cited in [CPS20]; *(unverified)* against the publisher).
- [BS91] F. A. Berezin, M. A. Shubin, *The Schrödinger Equation*, Springer (Kluwer), Dordrecht, 1991 (as cited in
  [CPS20]; *(unverified)* against the publisher).
- [Sim05] B. Simon, "Sturm oscillation and comparison theorems", in W. O. Amrein, A. M. Hinz, D. B. Pearson (eds.),
  *Sturm–Liouville Theory: Past and Present*, Birkhäuser, Basel, 2005 (as cited in [CPS20]; *(unverified)* against the
  publisher).
