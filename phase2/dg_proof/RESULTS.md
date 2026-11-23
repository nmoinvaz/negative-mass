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
| `rescaled.py` | module: the rescaled coefficients a_n(ε) of §4.1 and exact rational interval arithmetic (series plus tail enclosure) |
| `vhat_structure.py` | Theorem 3 (§4.2): V_Z > 0 and V_Z′ < 0 beyond r_s, the quadratic form Ψ, planted toy, zero counts |
| `threshold_isolation.py` | certified isolation of c·r_* and c·r_neg for every l = 2..3000, sympy Sturm cross-check, monotonicity tables; `--csv=PATH` writes every bracket as exact fractions |
| `compare_brackets.py` | #139: compares that csv with #104's `sweep_brackets.csv` for overlap at every l (`compare_brackets_output.txt`) |
| `zerilli_derivation.py` | #143: V from the linearized Einstein equations, M symbolic, s ≠ 0 (§4.2, §4.6; note `zerilli_derivation.md`; controls in `zerilli_derivation_controls.py`) |
| `large_l_proof.py` | Theorem L (§4.4): explicit constants for the limits ρ_*, ρ_neg |
| `monotone_proof.py` | Theorem M (§4.5): monotone approach for every l ≥ 500 |
| `check_mutations.py` | runs every new script unmutated and once per `--mutate=...`, which must exit non-zero |

The new scripts run from this directory as `OMP_NUM_THREADS=1 python3 -u X.py > X_output.txt` (`threshold_isolation.py` takes
12 min on one core of a Raspberry Pi 5, the others under a minute each). `check_mutations.py` writes `check_mutations_output.txt`.

## 1. Lemmas (all l ≥ 2; algebra checked symbolically in l)

**Lemma A (χ).** χ = Σ_{n=1}^{l+2} c_n rⁿ with c_n = (n−2)((n−4)l(l+1)+n−1)·(l−n+3)_{2n−4}/(2ⁿ(n−1)!²).
(i) c₁ = 3/2, c₂ = 0, c₃ = −K/32 with K = (l−1)l(l+1)(l+2) = 6k, c₄ = K/96. Every c_n with n ≥ 4 is ≥ 0, and c_{l+2} > 0.
(ii) The static Zerilli equation for ψ₀ = χ/(cr−6) is equivalent to P₂χ″ + P₁χ′ + P₀χ = 0 with P₂ = −r²(r+2)(cr−6).
On rⁿ this ODE gives a three-term recurrence 12(m−1)²c_m − 2(m−2)(c(m−4) − 3(m−1))c_{m−1} + c(l−m+4)(l+m−3)c_{m−2} = 0. Dividing by
(l−m+5)_{2m−8} turns it into a polynomial identity in (l, m), which is verified to be 0. It covers 4 ≤ m ≤ l+4, which
includes the termination orders m = l+3, l+4. The orders m ≤ 3 are checked symbolically in l, and for m ≥ l+5 every
term vanishes. The series terminates because (l−n+3)_{2n−4} = 0 for n ≥ l+3. So **ψ₀ solves the static
Zerilli equation for every l**. DG checked this only for small l; phase 0 checked l ≤ 5.
(iii) **χ > 0 on r > 0** (`lemmas_all_l.py` calls this inequality L5). DG's argument (eqs. 57–59), checked and put in
closed form:
- For l = 2, χ/r = 3/2 − (3/4)r² + (1/4)r³ = 1/2 + (r−2)²(r+1)/4, whose minimum is 1/2 at r = 2. The closed form below
  is 0/0 at K = 24 (l = 2), with limits r₀ = 2 and q(r₀) = 1/2; that is why l = 2 is separate.
- For l ≥ 3, χ/r ≥ q(r) := 3/2 − (K/32)r² + (K/96)r³ + c₅r⁴ with c₅ = K(K−24)/6144, since the dropped terms are
  non-negative by (i). With a = K/32, b = K/96 and s = √(6(K−18)), q′ = r(−2a + 3br + 4c₅r²) has exactly one positive
  root r₀ = 24/(s+6) = 4(s−6)/(K−24), which is DG eq. (58); the other root is −24/(s−6). We find
  q(r₀) = 4(K−18)(K s − 18(K−16)) / (K−24)³, K − 24 = (l−2)(l+3)(l²+l+4), and
  (Ks)² − (18(K−16))² = 6(K−24)³. For K > 24 the right side is > 0. Since Ks > 0 and 18(K−16) > 0, this gives
  Ks > 18(K−16), and 4(K−18)/(K−24)³ > 0. So q(r₀) > 0. Only K > 24 is used; K ≥ 120 is the range l ≥ 3.
- DG's "unique critical point" claim is correct (it is the positive root of a quadratic with roots of opposite sign).
  So is their non-negativity claim for the dropped terms.

**L1 (r_s).** The χ-ODE is singular at r_s (P₂(r_s) = 0) and has exponents 0 and 3 there. The ODE at r_s
gives χ′/χ(r_s) = −1/(r_s(r_s+2)), and the differentiated ODE at r_s gives χ″(r_s) = 0 (#138; the exponents alone do not). Consequences:
- φ₀,ₓₓ(r_s) = 0, so V̂ = φ₀,ₓₓ/φ₀ is analytic at r_s.
- The residue of ψ₀² dx/dr at r_s is 0.
- The DG mode mapped by A*, Ψ̂ = (−∂ₓ + g)Ψ_DG, is regular at r_s (L6).

Combined with Lemma A, V̂ is real-analytic on (0, ∞).

**L2/L3 (ends).**
- At r = 0: r⁴V̂ = 12 + 2(l²+l+3)r + l(l+1)r² + O(r³). The Frobenius exponents in r are 3 and −1 for every l, i.e.
  x^{3/2} and x^{−1/2}. They differ by an integer, so the −1 branch may carry a term ∝ x^{3/2} ln x; this does not
  change what follows. The x^{−1/2} branch is not in L²(dx), so Ĥ is **limit point at 0**. In x,
  V̂ = 3/(4x²) + (l²+l−1)/(4x^{3/2}) + 0·x⁻¹ + O(x^{−1/2}), so V̂ ≥ 3/(4x²) near 0, which gives the same conclusion.
  DG's eq. (61) (arXiv:0809.3615) prints −ℓ(ℓ+1)/(4|M|x) for the x⁻¹ term. The exact coefficient is 0 (#142, and an
  exact series reversion at l = 2 in raw `I142-mac`); nothing here uses that term.
- At ∞: V̂ ≃ (l+1)(l+2)/r² ≥ 0, so V̂ is bounded below near ∞ and Ĥ is limit point at ∞.
- So Ĥ is essentially self-adjoint on C_c^∞(0,∞): Weyl alternative, Reed–Simon II Thm X.7, or Weidmann,
  *Spectral Theory of ODOs*, LNM 1258, §5. This is DG's "unique self-adjoint extension", with domain functions
  ~ x^{3/2}. V̂ is continuous, → +∞ at 0 and → 0 at ∞, so Ĥ ≥ inf V̂ > −∞.
- φ₀ ∼ −4/r is the excluded branch at 0. At ∞, φ₀ ∝ r^{−(l+1)} lies in L². Since χ > 0, φ₀ has exactly one zero, at
  r_s, and it is simple (φ₀,ᵣ(r_s) = c/χ(r_s) ≠ 0).

**L4 (DG mode, all l).** Ψ_DG = r(r+2)^k e^{−kr/2}/(cr−6) satisfies H_Z Ψ_DG = −(k/2)²Ψ_DG identically in l. With k kept
symbolic, (H_Z + (k/2)²)Ψ_DG/Ψ_DG = (r+2)(c² + 2c − 6k)/(r²(cr−6)), so the identity holds exactly when k = c(c+2)/6
(#142). Since Ĥ A* = A* H_Z (because Ĥ = A*A and H_Z = AA*), Ψ̂ = A*Ψ_DG satisfies
ĤΨ̂ = −(k/2)²Ψ̂. Ψ̂ = O(r³) at 0 (L6, all l), it is regular at r_s, and it decays like e^{−kr/2}r^{k+…}. Ĥ is essentially
self-adjoint (L2/L3), so D(Ĥ) is the maximal domain {u ∈ L²: u, u′ ∈ AC_loc, Ĥu ∈ L²}. Ψ̂ is smooth on (0, ∞) (L6), in
L² at both ends, and ĤΨ̂ = −(k/2)²Ψ̂ ∈ L². So Ψ̂ ∈ D(Ĥ); the O(r³) is more than membership needs.
Direct check: ĤΨ̂/Ψ̂ = −(k/2)² for l = 2..10 (`verify_small_l`), and for l = 2..20 from the definition of Ψ̂ (#142).

**L6 (Ψ̂ at 0 and r_s).** With W = r(r+2)^k e^{−kr/2}, Ψ̂ = f(Wχ′ − W′χ)/((cr−6)χ).
- At 0, Ψ̂ = (2ᵏk²/192) r³ + O(r⁴). The lower orders cancel because c₃ = −3k/16 and c₄ = k/16, and c₅ = K(K−24)/6144
  enters the r³ coefficient.
- At r_s, W′/W = 1/r_s + k/(r_s+2) − k/2 equals χ′/χ = −1/(r_s(r_s+2)) (L1) exactly when k = c(c+2)/6, so Ψ̂ has no pole
  there. k = c(c+2)/6 therefore enters twice, in L4 and here. The other candidate poles are the zeros of χ, and there
  are none (Lemma A(iii)), so Ψ̂ is smooth on (0, ∞).

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
  and the cubic bracket has its minimum on [0, r_s] at r = r_s/(1+ε), equal to 9/4 + (27/4)ε²/(1+ε)² ≥ 9 − 27/4 (#138).
  The identity holds term by term, for any coefficients with c₁ = 3/2 and c₂ = 0.
- N has exactly one root r_* on (r_s, ∞). Proof: N(r_s + t) = cΣ(1−m)d_m t^m with d_m = χ^{(m)}(r_s)/m!
  = Σ_{n≥m} c_n C(n,m) r_s^{n−m}. Here d₀ = N(r_s)/c > 0 by T1, d₂ = 0 by L1, and d_m > 0 for 4 ≤ m ≤ l+2, since every
  c_n with n ≥ 4 is ≥ 0 and c_{l+2} > 0. The coefficient list then has exactly one sign change, whatever the sign of d₃.
  Descartes' rule gives at most one positive root, and N(r_s) > 0 > N(∞) gives one.

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

## 4. Large-l thresholds: r_* < r_neg for every l, and the limits

Notation: ρ = c·r, ε = 1/c, X(ρ) = cχ(ρ/c), N(ρ) = X − (ρ−6)X′ (the N of §3 in the variable ρ), φ₀ = c(ρ−6)/X. r_* is
the root of N beyond r_s (§3). **r_neg is the largest zero of V̂ on (0,∞)**, the upper edge of {V̂ < 0}. Theorem 3 shows it is the only
zero of V̂ beyond r_*.

### 4.1 Rescaled coefficients

With L = l(l+1) = c + 2, pairing the factors l−a and l+1+a gives (l−n+3)₍₂ₙ₋₄₎ = Π_{a=0}^{n−3}(L − a(a+1)), and
L − a(a+1) = c − (a−1)(a+2). Hence X = Σ_{n=1}^{l+2} aₙ(ε)ρⁿ with a₁ = 3/2, a₂ = 0 and, for n ≥ 3,

  aₙ(ε) = (n−2)[(n−4) + 3(n−3)ε] / (2ⁿ(n−1)!²) · Π_{a=0}^{n−3} (1 − (a−1)(a+2)ε).

The product vanishes at a = l, so the sum stops at n = l+2. At ε = 0 it is X_∞ = Σ(n−2)(n−4)(ρ/2)ⁿ/(n−1)!², which equals
s[(s+3)I₀(2√s) − 4√s I₁(2√s)] with s = ρ/2 (`large_l_limit.py`, and the roots agree to 18 digits in `threshold_isolation.py` B). Checked
exactly against `common.chi_poly` for l = 2..30, every n (`threshold_isolation.py` A).

### 4.2 Theorem 3: the shape of V̂ beyond r_*

**Theorem 3.** For every l ≥ 2, V̂ has exactly one zero r_neg in (r_*, ∞). It is simple, V̂ < 0 on [r_*, r_neg) and V̂ > 0 on
(r_neg, ∞). In particular **r_* < r_neg for every l ≥ 2**.

*Proof.* Let V = 2f[λ²r²((λ+1)r−3) + 9(λr−1)]/(r³(λr−3)²), λ = c/2, be the static Zerilli potential for M = −1,
g = ψ₀′/ψ₀ (′ = d/dx) and h = √(V/2). V is derived from the linearized vacuum Einstein equations in Regge–Wheeler
gauge, for time dependence e^{st} with s ≠ 0 and the master variable of Dotti–Gleiser gr-qc/0604021 eq. (9), in
`zerilli_derivation.md` (#143, checked symbolically in (r, l, M)). Its double pole at λr = 3 is r_s = 6/c. The static
case s = 0, which (Z) uses, is #145.
- (Z) ψ₀,ₓₓ = Vψ₀ (Lemma A(ii), and checked exactly with this V for l = 2..14). So g′ = V − g² and
  **V̂ = g² − g′ = 2g² − V**.
- (P1) V > 0 on (r_s, ∞). With λr = 3 + σ and μ = 1/λ, the numerator of V is a polynomial in (σ, μ) with positive
  coefficients. Also f > 0.
- (P2) V_x < 0 on (r_s, ∞). V_r = −2Q/(r⁵(cr−6)³) with Q a quintic in r. With r = (6+τ)/c and l = 2+m, c⁵Q is a polynomial in (τ, m)
  whose 64 coefficients are all positive, so Q > 0 for τ ≥ 0, m ≥ 0, that is for every real l ≥ 2. And dx = dr/f.
- (B) On (r_*, ∞), κ_c < 0 (§3), so g = −κ_c > 0, and g(r_*) = 0. Since V̂ = 2(g−h)(g+h), sign V̂ = sign(g−h). At a zero z of
  g − h we have g′ = V − g² = V/2 > 0 and h′ = V_x/(4h) < 0, so (g−h)′(z) > 0. A C¹ function whose every zero is a strict upward
  crossing has at most one zero. At r_*, g − h = −h < 0. At infinity V̂ ≃ (l+1)(l+2)/r² > 0 (L3; the leading coefficient of Ψ below gives
  the same sign). So there is exactly one zero. ∎

The quadratic form used from here on: with g = −κ_c = −c f N/((ρ−6)X), f = c(2+ερ)/ρ and V = 8c⁴(2+ερ)N_Z/(ρ⁴(ρ−6)²),

  **V̂ = 2c⁴(2+ερ) Ψ / (ρ⁴(ρ−6)²X²),  Ψ = ρ²(2+ερ)N² − 4N_Z X²,  N_Z = (1+2ε)ρ³/8 − 3ρ²/4 + 9ρ/2 − 9.**

N_Z > 0 on ρ ≥ 6 (positive coefficients in ρ = 6+τ). So sign V̂ = sign Ψ on ρ > 6, and at r_*, Ψ = −4N_ZX² < 0.
(`vhat_structure.py`: the identities symbolic in c for a generic X and exact against `chi_poly` for l = 2..9; P1, P2 as coefficient
positivity; a planted toy ψ = x²+3 whose V̂ has the known single zero x = 1; exact Sturm counts: for l = 2..16, N has one root beyond ρ = 6
and S, an equivalent polynomial of V̂, has exactly one root in [ρ_*, ∞) and three positive roots in all.) Status: *proved* (by hand), with the
polynomial positivity and the identities *checked symbolically*.

### 4.3 Certified isolation, l = 2..3000

For each l the sign of N (and of V̂ through the polynomial S = ρX[εN+(2+ερ)N′] − (2+ερ)N(X+2ρX′), which has the sign of V̂ by an identity
checked in `vhat_structure.py`) is evaluated in exact rational arithmetic on a bracket of width 2·10⁻¹³ (series truncated at 90
terms, omitted terms enclosed by a ratio test). Uniqueness of the root is T2 and Theorem 3, so each bracket contains r_* c, r_neg c.
Cross-checks that share no machinery: sympy Sturm isolation of the polynomials built from `chi_poly` (l = 2..24, brackets inside
sympy's intervals), the stored table of `large_l_limit_output.txt` (l = 2, 3, 5, 10, 20, 40, 80), and the Bessel closed form for the
limits. Output: `threshold_isolation_output.txt`. Status: *certified* (trusted code: `rescaled.py`, `threshold_isolation.py`).

#139 wrote every bracket (raw `I139-rpi5/brackets_i2.csv`) and compared it with #104's independent sweep: all 2999 pairs overlap,
for c·r_* and for c·r_neg (`compare_brackets_output.txt`; widest bracket 2.0·10⁻¹³ here, 2.0·10⁻¹⁴ in #104, largest distance
between midpoints 4.6·10⁻¹⁵). The ends of a bracket depend on the machine, since `locate` seeds each one with a float from
scipy's `brentq`. On the coordinator's Mac, 35 of the 2999 rows differ from rpi5's, and each pair overlaps. The printed table
does not change (raw `I139-mac`).

| l | c·r_* | c·r_neg | l²·r_* | l²·r_neg |
|---|---|---|---|---|
| 2 | 9.115122 | 10.696286 | 9.115122 | 10.696286 |
| 3 | 8.789041 | 10.143519 | 7.910137 | 9.129167 |
| 5 | 8.604503 | 9.833910 | 7.682592 | 8.780277 |
| 10 | 8.517779 | 9.689316 | 7.886832 | 8.971589 |
| 20 | 8.494079 | 9.649909 | 8.128306 | 9.234363 |
| 80 | 8.486244 | 9.636891 | 8.384063 | 9.520855 |
| 100 | 8.486049 | 9.636568 | 8.403693 | 9.543046 |
| 500 | 8.485715 | 9.636013 | 8.468845 | 9.616856 |
| 1000 | 8.485704 | 9.635995 | 8.477244 | 9.626388 |
| 3000 | 8.485701 | 9.635990 | 8.482876 | 9.632781 |
| ∞ | 8.485700967 | 9.635989443 | 8.485700967 | 9.635989443 |

For every l = 2..3000 (interval separation of consecutive brackets):
- r_* < r_neg;
- c·r_* and c·r_neg decrease strictly with l and stay above their limits;
- l²·r_* and l²·r_neg both have their minimum at l = 5 (7.6826 and 8.7803), increase strictly from l = 5 on, and stay below the limits for l ≥ 3.

So the old "9/l²" is l²·r_neg seen at moderate l, approaching 9.636 slowly from below, and it is non-monotone only below l = 5. Sparse larger l
(refined to 10⁻²²): (c·r_* − ρ_*)·c = 3.515486 and (c·r_neg − ρ_neg)·c = 5.839393 at l = 10⁷, so c·r_* ≈ ρ_* + 3.5155/c and
c·r_neg ≈ ρ_neg + 5.8394/c (*numerical evidence* for the coefficients, *proved* as bounds in 4.4–4.5).

### 4.4 Theorem L: the limits, with a rate

**Theorem L.** For every integer l ≥ 138,
  |c·r_* − ρ_*| ≤ 91.51/c  and  |c·r_neg − ρ_neg| ≤ 229.87/c,
with ρ_* = 8.485700967373530… and ρ_neg = 9.635989442586144…, the simple roots of N_∞ = X_∞ − (ρ−6)X_∞′ and Ψ_∞ = 2ρ²N_∞² − 4N_Z⁰X_∞².
Hence c·r_* → ρ_* and c·r_neg → ρ_neg, i.e. **r_* ≈ 8.486/((l+2)(l−1))** and **r_neg ≈ 9.636/((l+2)(l−1))** as limits, not only as fits.

*Proof.* Four steps, the constants from `large_l_proof_output.txt`.
- Coefficients. For n ≥ 3, aₙ = (αₙ + εα′ₙ)Pₙ(ε), αₙ = (n−2)(n−4)/(2ⁿ(n−1)!²), α′ₙ = 3(n−2)(n−3)/(2ⁿ(n−1)!²) and
  Pₙ = Π(1 − qₐε), q_a = (a−1)(a+2). For n ≤ l+2 every factor with a ≥ 1 lies in [0,1] (q_a ≤ c) and the a = 0 factor is 1+2ε, so
  0 < Pₙ ≤ 1+2ε, and by Weierstrass' inequality Pₙ ≥ (1+2ε)(1 − εSₙ), Sₙ = (n−2)(n²−4n−3)/3 + 2. Hence |aₙ − αₙ| ≤ ε βₙ with
  βₙ = α′ₙ(1+2ε₀) + |αₙ|max(2, Sₙ(1+2ε₀)), ε₀ = 1/108 (l ≥ 10). For n ≥ l+3, aₙ = 0 and c < n², so |aₙ − αₙ| ≤ ε n²|αₙ|. Checked against the exact aₙ for
  l = 10..240 (`large_l_proof.py` 1).
- Perturbation. N_l − N_∞ = Σ(aₙ − αₙ)ρ^{n−1}(6n − (n−1)ρ), bounded coefficientwise: |N_l − N_∞| ≤ 259.5 ε on [8, 9], and
  (Ψ_l − Ψ_∞)/ε is an explicit polynomial in X, N and their ε-differences (checked exactly at 20 random rational points),
  enclosed by rational interval arithmetic: |Ψ_l − Ψ_∞| ≤ 6.93·10⁶ ε on [9.61, 9.66].
- Slopes. −N_∞′ = (ρ−6)X_∞″ ≥ 2.836 on [8, 9]. Ψ_∞′(ρ_neg) = 33523 and |Ψ_∞″| ≤ 1.31·10⁵ on [9.61, 9.66], so Ψ_∞′ ≥ 30130 there
  (the closed forms of Ψ_∞′, Ψ_∞″ are checked against difference quotients).
- Localisation. If F is C¹ with F′ ≥ μ on [ρ̄−δ, ρ̄+δ], F(ρ̄) = 0, and |F_l − F| ≤ ηε < μδ there, then F_l changes sign inside
  [ρ̄ − ηε/μ, ρ̄ + ηε/μ]. For N the root in that interval is r_*c (T2: unique root). For Ψ it lies beyond r_*c once ε ≤ 5.2·10⁻⁵, so it is
  r_neg c (Theorem 3: unique zero beyond r_*). K_* = 259.5/2.836, K_neg = 6.93·10⁶/30130, valid for ε < 5.2·10⁻⁵, which is l ≥ 138. ∎

The constants are rigorous but 25 to 40 times the true coefficients 3.5155 and 5.8394. Status: *certified*. What the code proves: the
inequalities listed (coefficient bound at the l tested, interval enclosures with tail bounds, the final constants). It assumes T2,
Theorem 3, the identity Ψ ↔ V̂, and the hand-proved coefficient bound beyond l = 240. It checks the final inequalities against the certified
brackets of 4.3 for l = 138, 139, 145, 276, 690, 1000, 2000, 3000.

### 4.5 Theorem M: monotone approach

**Theorem M.** For every integer l ≥ 500: c·r_*(l+1) < c·r_*(l) and c·r_neg(l+1) < c·r_neg(l), with
0 < c·r_*(l) − ρ_* ≤ 3.709/c and 0 < c·r_neg(l) − ρ_neg ≤ 6.101/c, and (l+1)²·r_*(l+1) > l²·r_*(l), (l+1)²·r_neg(l+1) > l²·r_neg(l)
(this read "l²·r_*(l+1) > l²·r_*(l)" before #104's review, which is false as written).

*Proof.* Expand to second order: aₙ = αₙ + εα¹ₙ + ε²ωₙ with α¹ₙ = α′ₙ − αₙe₁(n), e₁(n) = Σqₐ = (n−2)(n²−4n−3)/3 and
|ωₙ| ≤ (|αₙ| + εα′ₙ)·½(e₁+4)²(1+2ε) + α′ₙ|e₁| (Taylor remainder of Pₙ, using |Pₙ″| ≤ (Σ|qₐ|)²(1+2ε) and Σ|qₐ| = e₁+4).
So N_l = N_∞ + εν + ε²ζ and Ψ_l = Ψ_∞ + εΨ₁ + ε²E, with explicit ν, Ψ₁ (entire functions) and |ζ| ≤ 1488, |E| ≤ 4.67·10⁷ on the windows
ρ_* ± 92/c₅₀₀, ρ_neg ± 230/c₅₀₀ (c₅₀₀ = c(500)).
- ν ≥ 19.61 on the N window, and Ψ₁ ≤ −1.89·10⁵ on the Ψ window.
- By Theorem L, ρ_l = c·r_*(l) lies in the N window (and c·r_neg(l) in the Ψ window) for l ≥ 138. Since N_l(ρ_l) = 0, N_{l+1}(ρ_l) = (ε′−ε)ν + ε′²ζ′ − ε²ζ, with ε′ = ε_{l+1} and ε − ε′ = (2l+2)εε′. This is < 0 as soon as
  (l+1)ν_min/2 > |ζ|max, which is 501·19.61/2 = 4913 > 1488. So ρ_*(l+1) < ρ_*(l). The same computation for Ψ: Ψ_{l+1}(ρ_neg(l)) > 0 since
  501·1.89·10⁵/2 = 4.74·10⁷ > 4.67·10⁷. So ρ_neg(l+1) < ρ_neg(l).
- Localisation: N_l(ρ_*) = εν + ε²ζ > 0 gives ρ_*(l) > ρ_*, and N_l(ρ_* + sε) < 0 for s = 3.709, since −N_∞′ ≥ 5.599 on the narrow
  window and 1.05(ν_max + εζ_max)/5.599 = 3.709. The same for Ψ with Ψ_∞′ ≥ 33523 − 20ε|Ψ_∞″|/2 and s = 6.101.
- The decrease of l²r: l²r_{l+1} > l²r_l iff ρ_l − ρ_{l+1} < ρ_l(l²−3l−2)/((l+1)²c_l). Since ρ_{l+1} > ρ_∞, the left side is below
  s/c_l, and s < ρ_∞(l²−3l−2)/(l+1)² (3.709 < 8.40, 6.101 < 9.54, increasing in l). ∎

With 4.3 (consecutive brackets, needed for l < 500) this gives: **c·r_* and c·r_neg decrease strictly for every l ≥ 2, so they approach ρ_* and ρ_neg
from above; l²·r_* and l²·r_neg increase strictly for every l ≥ 5.** Status: *certified*, with the same reading as 4.4. `monotone_proof.py`
also checks the expansions at l = 500, 501, 502, 1000, 1001 against exact evaluation, and N_{l+1} < 0 at the lower end of r_*(l)'s
bracket (Ψ_{l+1} > 0 at the lower end of r_neg(l)'s) for l = 500..502.

### 4.6 What is not proved here

- The shape of V̂ on (0, r_*]: V̂ → +∞ at 0 and V̂(r_*) < 0, so there is a zero below r_*, and exact Sturm counts give exactly one zero in
  (r_s, r_*) for l = 2..16 (three positive roots of S in all), but this is not proved for all l. The lower edge ρ → 6.16881 of the limit
  is *numerical evidence* (`large_l_limit.py`).
- What the code proves, apart from the theorems: the coefficient identity of 4.1 (exact, l ≤ 30), the polynomial identities and coefficient
  positivity (sympy 1.14.0, exact), and the numerical inequalities of 4.3–4.5 (exact rational arithmetic, `fractions.Fraction`; floating point
  only guesses brackets, every bracket is certified in rationals). The facts assumed beyond that: T1, T2 and Lemma A(ii) of this file, L3 (the
  sign of V̂ at infinity, also seen from the leading coefficient of Ψ), the hand-proved coefficient and tail bounds (stated in the proofs,
  tested numerically for n ≤ 400), and elementary analysis (intermediate value, Taylor, Weierstrass' product inequality).
- **Independent review (#104).** `review_large_l.md` (mpmath `iv` and sympy, written from this file alone) confirms 4.2–4.5 and
  finds no error: Theorem 3's (P1), (P2) and V̂ identities, the brackets for every l = 2..3000, and Theorems L and M with its own,
  smaller constants (L: 23.82 and 172.03; M: 3.519 and 5.936). It assumed T1, T2 and Lemma A(ii) for l > 14, which #138 then
  reviewed (next bullet). Five of its 52 controls missed their targets (`review_large_l.md` §6). The coordinator's check of the review is in
  negative-mass-raw `I104-mac`; its findings are in `review_large_l.md`, Coordinator's review.
- **Independent review of the inputs (#138).** `review_lemma_a.md` (sympy and exact rationals, written from this file alone)
  confirms Lemma A(i) and (ii), L1, T1 and T2 for every l, and finds no error. Lemma A(ii) is an exact identity in (l, m).
  ψ₀,ₓₓ = Vψ₀ is exact for l = 2..40, and at l = 41, 45, 60 and 100 in the coordinator's check (raw `I138-mac`). Its six
  wording notes are applied in §1 and §3. Lemma A(iii), L2–L6 and the Zerilli V were reviewed next (#142, next bullet).
- **Independent review of the remaining inputs (#142).** `review_lemma_a3.md` (sympy and exact rationals, written from
  this file alone) confirms Lemma A(iii), L2, L3, L4 and L6 and finds no error. By exact Sturm counts, χ/r has no
  positive root for l = 2..50. The V of §4.2 at M = −1 equals four published forms: Dotti–Gleiser gr-qc/0604021 eq. (6)
  and arXiv:0809.3615 eq. (11), Kokkotas–Schmidt gr-qc/9909058 eq. (24) and Berti–Cardoso–Starinets arXiv:0905.2975
  eq. (15). Its eight wording notes are applied in §1. The coordinator's check (raw `I142-mac`) passes, and no number
  of #142 disagrees with #2. Not checked: Zerilli's own papers (the publisher page was blocked, *(unverified)*) and DG's
  intertwining equivalence (§6). The derivation of V is the next bullet (#143).
- **Derivation of V (#143).** `zerilli_derivation.md` (sympy, exact) derives V from the linearized vacuum Einstein
  equations for symbolic M: Regge–Wheeler gauge, e^{st} with s ≠ 0, l ≥ 2. The even-parity system is a two-dimensional
  first-order flow, and DG06's master variable (9) turns it into −ψ,ₓₓ + Vψ = −s²ψ; at M = −1 this V equals §4.2's
  exactly. The reduction divides by s, r − 2M and λr + 3M, so it holds on r ≠ r_s. It assumes the gauge is attainable.
  The coordinator's check (raw `I143-mac`, its own Lichnerowicz operator) confirms the converse for l = 2, 3, 4 and the
  forward direction by an exact certificate for l = 2 (M symbolic) and l = 3 (M = −1). Still not checked: Zerilli's own
  papers (publisher HTTP 403, *(unverified)*), DG's intertwining equivalence, the gauge's attainability and the behavior
  at r_s (#146), and the case s = 0 (#145).

## 5. Verification summary

- **Symbolic in l (all l):** Lemma A(i–iii), L1–L6, T1, T2 (T1's identity holds term by term, hence for
  every l; #138), T3, T4.
- **Exact for l = 2..10 (`verify_small_l_output.txt`):**
  - χ/r has no positive root, and den V̂ = r⁴(χ/r)²·const.
  - Ĥφ₀ = 0, and Res_{r_s} = 0. F′ = ψ₀² dx/dr, F > 0 on (0, r_s), and v has exactly one zero beyond r_s.
  - Ψ̂ has no pole at r_s and ĤΨ̂/Ψ̂ = −(k/2)² exactly. Ψ̂ ∼ r³.
  - N has exactly one root on (r_s, ∞).
- **Exact Sturm counts for l = 2..50 (#142, `review_lemma_a3_a3_output.txt`):** χ/r has no positive root.
- **Exact root isolation of N for l = 2..50:** no root in (0, r_s), one root in (r_s, ∞), one Descartes sign change.
- **Certified, exact rational arithmetic (§4):**
  - c·r_* and c·r_neg bracketed for every l = 2..3000 (`threshold_isolation_output.txt`), with the cross-checks of 4.3.
  - Theorem 3's polynomial positivity (V_Z > 0 and V_Z′ < 0 beyond r_s, 64 coefficients), the identities V̂ ↔ Ψ and V̂ ↔ S, exact against
    `chi_poly` for l = 2..9 (`vhat_structure_output.txt`).
  - Theorem L (l ≥ 138) and Theorem M (l ≥ 500) (`large_l_proof_output.txt`, `monotone_proof_output.txt`).
  - Every check can fail: `check_mutations_output.txt` lists each script run unmutated (exit 0) and once per `--mutate=...` (exit 1).
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
- **Proved in §4 (status words as there):**
  - *proved*: r_* < r_neg for every l ≥ 2, from Theorem 3 (V̂ has exactly one zero beyond r_*, simple, from − to +).
  - *certified*: c·r_* → ρ_* = 8.485700967… and c·r_neg → ρ_neg = 9.635989443… (Theorem L, rates 91.5/c and 229.9/c for l ≥ 138);
    both sequences decrease strictly in l for every l ≥ 2 and stay above their limits, and l²·r_* and l²·r_neg increase strictly for
    every l ≥ 5 (certified brackets to l = 3000, Theorem M beyond l = 500). So r_neg ≈ 9.636/((l+2)(l−1)) and r_* ≈ 8.486/((l+2)(l−1)) are limits.
  - *numerical evidence*: the next coefficients, c·r_* ≈ ρ_* + 3.5155/c and c·r_neg ≈ ρ_neg + 5.8394/c (certified only as the bounds 3.709/c and
    6.101/c for l ≥ 500).
- **Not proven:** the number of zeros of V̂ on (0, r_*) for general l (exactly one in (r_s, r_*) is exact for l = 2..16 only), and the
  convergence of its lower edge to 6.16881. The certified parts had an independent review in #104 (§4.6).
- **Caveat on scope:** these are statements about Ĥ. Their meaning for Zerilli/metric perturbations relies on DG's
  intertwining equivalence, which is not re-proved here. V itself agrees with four published forms of the Zerilli
  potential at M = −1 (#142), and it is derived from the linearized Einstein equations for s ≠ 0 on r ≠ r_s in
  `zerilli_derivation.md` (#143). The static case s = 0 is open (#145).
