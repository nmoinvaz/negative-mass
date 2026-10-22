# The Dotti–Gleiser operator has exactly one negative eigenvalue for every l ≥ 2

Units |M| = 1 (M = −1), f = 1 + 2/r, x = r − 2 ln(1 + r/2), c = (l+2)(l−1) = 2λ, r_s = 6/c,
k = (l−1)l(l+1)(l+2)/6 = c(c+2)/6. χ is the polynomial of DG eq. (55), ψ₀ = χ/(cr − 6), g = ψ₀′/ψ₀ (′ = d/dx),
V̂ = g² − g′, Ĥ = −∂ₓ² + V̂ = A*A with A = ∂ₓ + g. φ₀ = 1/ψ₀ = (cr − 6)/χ, so Aφ₀ = 0 and Ĥφ₀ = 0.

Scripts are in this directory. Each `X.py` writes `X_output.txt`.

| script | what it establishes |
|---|---|
| `chi_lemma.py` | Lemma A: coefficients of χ in closed form with their signs; ψ₀ solves static Zerilli for **all l** (recurrence identity, symbolic in l); χ > 0 (DG's bound, closed form) |
| `lemmas_all_l.py` | L1–L6, symbolic in l: regularity at r_s, Frobenius exponents at r = 0, behaviour at ∞, DG mode (all l), DG positivity inequality |
| `threshold_lemmas.py` | T1–T4, symbolic in l: N > 0 on [0, r_s], unique root r_*, κ_c asymptotics |
| `ground_state_identity.py` | the Jacobi/ground-state identity |
| `verify_small_l.py` | exact sympy checks, l = 2..10, of every ingredient; closed-form admissible E = 0 solution and its zero |
| `kappa_c_threshold.py` | κ_c closed form against dg_spectrum; Sturm/Descartes root isolation of N for l = 2..50 |
| `closed_form_F.py` | partial-fraction structure of the E = 0 solution, l = 2..15 |
| `large_l_limit.py` | l → ∞ limit (Bessel closed form), ρ_* and ρ_neg |
| `numeric_check.py` | numerical spectra at l = 9, 10, plus 108 truncated cases at l = 6, 8, 10 |

## 1. Lemmas (all l ≥ 2; algebra checked symbolically in l)

**Lemma A (χ).** χ = Σ_{n=1}^{l+2} c_n rⁿ with c_n = (n−2)((n−4)l(l+1)+n−1)·(l−n+3)_{2n−4}/(2ⁿ(n−1)!²).
(i) c₁ = 3/2, c₂ = 0, c₃ = −K/32 with K = (l−1)l(l+1)(l+2) = 6k, c₄ = K/96. Every c_n with n ≥ 4 is ≥ 0, and c_{l+2} > 0.
(ii) The static Zerilli equation for ψ₀ = χ/(cr−6) is equivalent to P₂χ″ + P₁χ′ + P₀χ = 0 with P₂ = −r²(r+2)(cr−6).
On rⁿ this ODE gives a three-term recurrence 12(m−1)²c_m − 2(m−2)(…)c_{m−1} + c(l−m+4)(l+m−3)c_{m−2} = 0. Dividing by
(l−m+5)_{2m−8} turns it into a polynomial identity in (l, m), which is verified to be 0; the low orders m ≤ 4 are
checked symbolically in l. The series terminates because (l−n+3)_{2n−4} = 0 for n ≥ l+3. So **ψ₀ solves the static
Zerilli equation for every l**. DG checked this only for small l; phase 0 checked l ≤ 5.
(iii) **χ > 0 on r > 0.** DG's argument (eqs. 57–59), checked and put in closed form:
- For l = 2, χ/r = 3/2 − (3/4)r² + (1/4)r³, whose minimum is 1/2 at r = 2.
- For l ≥ 3, χ/r ≥ q(r) := 3/2 − (K/32)r² + (K/96)r³ + c₅r⁴, since the dropped terms are non-negative by (i).
  q′ = r(−2a + 3br + 4c₅r²) has exactly one positive root r₀, which is DG eq. (58). We find
  q(r₀) = 4(K−18)(K s − 18(K−16)) / ((l−2)(l+3)(l²+l+4))³ with s = √(6(K−18)), and
  (Ks)² − (18(K−16))² = 6(K−24)³ > 0 for K ≥ 120. So q(r₀) > 0.
- DG's "unique critical point" claim is correct (it is the positive root of a quadratic with roots of opposite sign).
  So is their non-negativity claim for the dropped terms.

**L1 (r_s).** The χ-ODE is singular at r_s (P₂(r_s) = 0) and has exponents 0 and 3 there. This forces
χ′/χ(r_s) = −1/(r_s(r_s+2)) and χ″(r_s) = 0. Consequences:
- φ₀,ₓₓ(r_s) = 0, so V̂ = φ₀,ₓₓ/φ₀ is analytic at r_s.
- The residue of ψ₀² dx/dr at r_s is 0.
- The DG mode mapped by A*, Ψ̂ = (−∂ₓ + g)Ψ_DG, is regular at r_s (L6).

Combined with Lemma A, V̂ is real-analytic on (0, ∞).

**L2/L3 (ends).**
- At r = 0: r⁴V̂ = 12 + 2(l²+l+3)r + O(r²). The Frobenius exponents in r are 3 and −1 for every l, i.e. x^{3/2} and
  x^{−1/2}. The x^{−1/2} branch is not in L²(dx), so Ĥ is **limit point at 0**.
- At ∞: V̂ ≃ (l+1)(l+2)/r², so Ĥ is limit point at ∞.
- So Ĥ is essentially self-adjoint on C_c^∞(0,∞): Weyl alternative, Reed–Simon II Thm X.7, or Weidmann,
  *Spectral Theory of ODOs*, LNM 1258, §5. This is DG's "unique self-adjoint extension", with domain functions
  ~ x^{3/2}. V̂ is continuous, → +∞ at 0 and → 0 at ∞, so Ĥ ≥ inf V̂ > −∞.
- φ₀ ∼ −4/r is the excluded branch at 0. At ∞, φ₀ ∝ r^{−(l+1)} lies in L². Since χ > 0, φ₀ has exactly one zero, at
  r_s, and it is simple (φ₀,ᵣ(r_s) = c/χ(r_s) ≠ 0).

**L4 (DG mode, all l).** Ψ_DG = r(r+2)^k e^{−kr/2}/(cr−6) satisfies H_Z Ψ_DG = −(k/2)²Ψ_DG identically in l, with k kept
symbolic and then k = c(c+2)/6. Since Ĥ A* = A* H_Z (because Ĥ = A*A and H_Z = AA*), Ψ̂ = A*Ψ_DG satisfies
ĤΨ̂ = −(k/2)²Ψ̂. Ψ̂ = O(r³) at 0 (L6, all l), it is regular at r_s, and it decays like e^{−kr/2}r^{k+…}. So Ψ̂ ∈ D(Ĥ).
Direct check: ĤΨ̂/Ψ̂ = −(k/2)² for l = 2..10 (`verify_small_l`).

**Ground-state identity.** If V = φ″/φ and u = φh, then u′² + Vu² = φ²h′² + (φφ′h²)′ (`ground_state_identity.py`).

## 2. Theorem 1: full problem

**Theorem 1.** For every l ≥ 2, Ĥ (the closure on C_c^∞(0,∞) = the DG self-adjoint operator) has exactly one negative
eigenvalue. It is simple and equals −(k/2)², with eigenfunction A*Ψ_DG. 0 is not an eigenvalue.

*Proof.*
- **At most one.** By Glazman's lemma (e.g. Birman–Solomyak, *Spectral Theory of Self-Adjoint Operators in Hilbert
  Space*, §10.2), dim Ran P_{(−∞,0)}(Ĥ) is the largest dimension of a subspace F of a form core on which q[u] = ⟨u, Ĥu⟩
  is negative definite. Take C_c^∞(0,∞) as the core (an operator core, hence a form core). Suppose dim F ≥ 2. Pick
  u ∈ F\{0} with u(x_s) = 0, where x_s = x(r_s). Then h = u/φ₀ is smooth on (0,∞): u vanishes at the simple zero of the
  analytic φ₀. Integrating the ground-state identity over supp u (compact), the boundary terms vanish, so
  q[u] = ∫φ₀²h′² ≥ 0. This is a contradiction.
- **At least one.** The explicit eigenfunction of L4. Alternatively, without it: w = φ₀·1_{x > x_s} is in the form
  domain and q[w] = [φ₀φ₀′]_{x_s}^∞ = 0. If there were no negative spectrum, q ≥ 0 and w would minimize q, so Ĥw = 0
  weakly, hence w ∈ D(Ĥ) ⊂ C¹. But w has a kink at x_s (φ₀′(x_s) ≠ 0), a contradiction.
- **Identification.** The eigenvalue is therefore −(k/2)².
- **No zero eigenvalue.** The E = 0 solution admissible at x = 0 is v below. It grows like r^{l+2}, so it is not in L². ∎

Only standard facts are used: limit-point classification, the form/operator correspondence (Kato VI §2) and Glazman's
lemma. No oscillation theorem at the continuum edge is needed.

**Closed-form admissible E = 0 solution (consistency with Sturm oscillation, all l).** v = φ₀F, with
F(r) = f.p.∫₀^r ψ₀² (dx/dr) dr. Since dF/dx = ψ₀², Ĥv = 0.
- Partial fractions: ψ₀² dx/dr = polynomial + A/(r−r_s)² + 0/(r−r_s) − (1/2)/(r+2), with A = r_sχ(r_s)²/(c²(r_s+2)) > 0.
  The B = 0 term is L1. The −1/2 coefficient was checked for l ≤ 15, and A for all l from the formula. So F is rational
  plus −½ ln(1 + r/2), and v is analytic at r_s.
- v ≃ −r³/32 at r = 0, the admissible branch: F = r⁴/128 + …
- F > 0 on (0, r_s). On (r_s, ∞), F increases from −∞ to +∞. Hence **v has exactly one zero, r_z ∈ (r_s, ∞), for every l**.
  Values: l = 2: 2.78312; l = 3: 1.05257; l = 10: 0.09276. These match the dg_spectrum ODE counts.
- This is what renormalized oscillation theory predicts: one zero ↔ one eigenvalue below 0 (Gesztesy–Simon–Teschl,
  Amer. J. Math. 118 (1996) 571; Teschl, *Mathematical Methods in QM*, 2nd ed., Ch. 9). We use it only as a check;
  the proof above does not depend on it.

## 3. Theorem 2: truncated problem

Ĥ_{R,κ} acts on L²(x_R, ∞) with the Robin condition Ψ̂′(x_R) = κΨ̂(x_R); κ = ∞ means Dirichlet. Its form is
q[u] = ∫_{x_R}^∞ (u′² + V̂u²) dx + κu(x_R)², with core C_c^∞([x_R,∞)). Define
**κ_c(R) = (d ln φ₀/dx)(R) = (1 + 2/R)·N(R)/((cR − 6)χ(R)),  N := cχ − (cr − 6)χ′.**

**Theorem 2.** For every l ≥ 2 and R ≠ r_s: #{negative eigenvalues of Ĥ_{R,κ}} = [R < r_s] + [κ < κ_c(R)],
where Dirichlet counts as κ = +∞. At R = r_s: 1 for every finite κ and 0 for Dirichlet. When κ = κ_c(R), φ₀ is an L²
zero-energy eigenfunction.

*Proof.*
- **Upper bounds (Glazman).** Every u in the core satisfies q[u] = ∫φ₀²h′² + (κ − κ_c)u(x_R)² as long as φ₀ ≠ 0 on
  supp u, or u(x_s) = 0 when R < r_s. Imposing u(x_s) = 0 (only if R < r_s) and u(x_R) = 0 (only if κ < κ_c) costs at
  most [R<r_s] + [κ<κ_c] linear conditions and leaves q ≥ 0. So no negative-definite subspace has larger dimension.
- **Lower bounds.** Use the test functions w₁ = φ₀·1_{[x_R, x_s]} and w₂ = φ₀·1_{[x_s,∞)} (for R < r_s), or φ₀ itself
  (for R > r_s). They give q[aw₁ + bw₂] = a²φ₀(x_R)²(κ − κ_c) and q[w₂] = 0.
- Suppose dim P₋ were smaller than claimed. Then some w in the span with w ⊥ Ran P₋ has q[w] = 0. Ran P₋ is finite
  dimensional and consists of eigenvectors, and q ≥ 0 on its orthogonal complement, so q(w, ·) = 0 and w ∈ D(Ĥ) with
  Ĥw = 0. Then w would have to be C¹ (the kink of w₂ excludes it) and satisfy the boundary condition (which forces
  κ = κ_c). This is a contradiction.
- The R = r_s statements follow from the same argument. ∎

**Corollary (thresholds, all l).** T1 and T2 give:
- N > 0 on [0, r_s]. Proof: N − [9 + 18c₃r² + (6c₄ − 2cc₃)r³] = Σ_{m≥4}(m−1)c_m r^{m−1}(6 − cr) + Σ_{j≥5}6c_j r^{j−1} ≥ 0,
  and the cubic bracket is ≥ 9 − 27/4 on [0, r_s].
- N has exactly one root r_* on (r_s, ∞). Proof: N(r_s + t) = cΣ(1−m)d_m t^m with d_m = χ^{(m)}(r_s)/m!. Here d₀ > 0,
  d₂ = 0 by L1, and d_m > 0 for 4 ≤ m ≤ l+2. Descartes' rule then gives exactly one sign change.

Hence κ_c < 0 on (0, r_s), κ_c → ∓∞ at r_s^∓, κ_c > 0 on (r_s, r_*) and κ_c < 0 on (r_*, ∞). Therefore:
- **Dirichlet is unstable iff R < r_s.**
- **Neumann is unstable iff R < r_*.**
- **Every κ ≥ 0 is stable when R ≥ r_*.**
- For R < r_s and κ ≥ 0 there is exactly one negative eigenvalue.

These hold for all l. They were previously only numerical, for l ≤ 20.

**κ_c asymptotics (symbolic in l).**
- κ_c = −(l+1)/R − (l⁴+3l³+4l²+4)/(l(l−1)(l+2)) R⁻² + O(R⁻³) as R → ∞. The second coefficient comes from
  c_{l+1}/c_{l+2} = (l+1)(l²−2l−2)/(l(l−1)).
- κ_c = −2/R² − (l²+l+1)/(3R) + O(1) as R → 0, i.e. κ_c·x_R = −1/2 − (l²+l−1)R/12 + O(R²).

## 4. Large-l thresholds (derived; limit not proven rigorously)

With ρ = c·r, the expansion of c·χ(ρ/c) converges coefficient by coefficient to X(ρ) = Σ(n−2)(n−4)(ρ/2)ⁿ/(n−1)!². In
closed form, X = s[(s+3)I₀(2√s) − 4√s I₁(2√s)] with s = ρ/2 (checked numerically against the series). So
φ₀ → (ρ − 6)/X(ρ) on the scale r ∼ 1/l². The limiting values are:
- **c·r_* → ρ_* = 8.48570097**
- **c·r_neg → ρ_neg = 9.63598944**
- the lower edge of the V̂ < 0 region tends to ρ = 6.16881.

Finite-l values of c·r_neg, from exact polynomial roots: 10.696 (l=2), 9.834 (5), 9.689 (10), 9.650 (20), 9.637 (80).
So **r_neg ≈ 9.636/((l+2)(l−1))**. The observed "9/l²" is this law seen at moderate l: l²r_neg = 8.78 at l = 5,
9.23 at l = 20 and 9.52 at l = 80, approaching 9.636 slowly from below. Likewise r_* ≈ 8.486/((l+2)(l−1)).
Not proven: that the root locations converge (this would need uniform convergence plus simple limit roots, which is
plausible but not written out), and that r_* < r_neg for all l (it holds for l ≤ 80 and in the limit).

## 5. Verification summary

- **Symbolic in l (all l):** Lemma A(i–iii), L1–L6, T1, T2 (except the explicit identity check, which was also run
  for l ≤ 12), T3, T4.
- **Exact for l = 2..10 (`verify_small_l_output.txt`):**
  - χ/r has no positive root, and den V̂ = r⁴(χ/r)²·const.
  - Ĥφ₀ = 0, and Res_{r_s} = 0. F′ = ψ₀² dx/dr, F > 0 on (0, r_s), and v has exactly one zero beyond r_s.
  - Ψ̂ has no pole at r_s and ĤΨ̂/Ψ̂ = −(k/2)² exactly. Ψ̂ ∼ r³.
  - N has exactly one root on (r_s, ∞).
- **Exact root isolation of N for l = 2..50:** no root in (0, r_s), one root in (r_s, ∞), one Descartes sign change.
- **Numerical:** at l = 9 and 10 the full problem has exactly one negative eigenvalue, E₀ = −(k/2)² to displayed
  precision, at two (N, a) settings. The truncated count formula holds in 108/108 cases at l = 6, 8, 10, which adds to
  the earlier 1000-case scan for l = 2..5.

## 6. Status

- **Rigorous for all l ≥ 2:** Theorem 1 (exactly one negative eigenvalue, simple, −(k/2)², no zero eigenvalue) and
  Theorem 2 (the truncated count formula), together with the threshold corollary (r_s, r_*, signs of κ_c, uniqueness of
  r_*) and the κ_c asymptotics. The proofs need the computer algebra listed above to be correct (symbolic identities in
  l) plus textbook theorems: Weyl limit-point/limit-circle; Kato's representation theorem; Glazman's lemma.
- **Corrected relative to dg_spectrum:** "χ > 0 for l ≤ 20 by root counting" is now a proof for all l (DG's argument,
  verified). "ψ₀ solves static Zerilli", previously checked for l ≤ 5, now holds for all l.
- **Not proven:** the l → ∞ limits ρ_* and ρ_neg as rigorous limits; r_* < r_neg for all l; any statement about the
  shape of V̂ (number of sign changes) for general l.
- **Caveat on scope:** these are statements about Ĥ. Their meaning for Zerilli/metric perturbations relies on DG's
  intertwining equivalence, which is not re-proved here.
