# The exchange identities and the Clausen integrals behind Ω₂ and Ω₃, derived again

November 2026. Issue #209 (milestone G5, Paper III), the independent re-derivation asked for by #26's task 4 and #208's E3.
The derivations of §2 to §5 were written from the order equations and the definitions alone, before `omega2.md` and
`omega3.md` were opened; §8 is the comparison made afterwards. The script is
[`exchange_review.py`](exchange_review.py), its output [`exchange_review_output.txt`](exchange_review_output.txt), the raw
logs `I209-rpi5` of negative-mass-raw.

Status words as in AGENTS.md: *derived* (by hand, each step checked), *checked symbolically* (sympy residual 0),
*numerical evidence* (agreement of integrals). A digit count is −log₁₀ of a relative difference between two evaluations.
It is not an error bar.

## Result

E1, E2, E3 and the two Clausen integrals are **derived** as stated, each integration by parts with its boundary terms (§2 to §5),
and checked.

| step | by parts on | boundary terms (both ends vanish) | where the hypothesis enters |
|---|---|---|---|
| E1 | ∫Y_k/sin²θ | [−cot θ Y_k]: O(θ) at 0, O(x l^p) at π/2 | Y_k = O(θ²), Y_k polynomial in l |
| E2 | ∫tan θ C_m; then ∫ℓ/sin²θ = −π/2 | [−ℓ C_m]: O(θ³) at 0, O(x l^{q+1}) at π/2; [−cot θ ℓ − θ]: 0 and −π/2 | G_m = 0 at π/2 |
| E3 | ∫Y_m P′, P = cot θ (ℓ − Ω₀) | [Y_m P]: 2Ω₀Ω_{m−2}θ at 0, O(x l^{m+2}) at π/2; m = 4 at k = 5 | the same two |

- **Hypotheses.** G_j = 0 for j < k and Y_k = 0 at l = 0. Ω_{k−2} stays free, and all three identities hold for any value of it.
  The slope of G_k in Ω_{k−2} is −4π, from E1 alone.
- **E3 at every order.** It removes Y_{k−1} for any k ≥ 2, not only Y₄ at k = 5. It collects two Y_{k−1} terms: ½Y_{k−1} inside E2's
  h_{k−1}, and the j = k − 1 term Y_{k−1}c₁ of h_k.
- **Checks.** Sympy residual 0 for every algebraic step (S1 … S11). The integrals agree to 40 to 45 digits on a panel hierarchy
  (E1 at k = 3, 4, 5; E2 at m = 2, 3, 4; E3 at k = 3, 4, 5; G₃, G₄, G₅ without Y_{k−1}), and to 50 digits from closed forms at
  k = 2 (E1) and m = 1, 2 (E2). The boundary terms decay or scale as predicted (table in §6).
- **Clausen.** ∫tan θ Cl₂(π − 2θ) = π³/12 from the moments of ln cos, and again from the sine series (two routes);
  ∫Cl₂(π − 2θ)² = π⁵/360 from orthogonality and uniform convergence, without completeness. Quadrature agrees to 60 digits (bar 40).
- **Third grid** (N = 64, w = 1/3, V = 24): Ω₂ against #26's A70 to 50.3 digits (A70 prints 48), Ω₃ against #208's B4 to 87.8 digits
  (bar 40 each).
- **Controls.** 27 of 27 on target (§9); the silences are listed there.
- **Comparison with omega2.md and omega3.md** (§8): no difference in a sign, constant or result. Differences in route and in what is
  stated: E2's route; boundary terms of E1, E2 and E3 at j = 4 stated and checked here; one incomplete boundary argument in
  omega2.md step 1 (the θ = 0 end of [D₃P]); E3 valid at every k; slope −4π derived.

## 1. Setting

Variables on the throat zone. R = r^s = sec²θ, θ ∈ [0, π/2), so that

```math
u = \sqrt{R-1} = \tan\theta,\qquad l = \ln R = -2\ell,\qquad \ell = \ln\cos\theta,\qquad S = \sin^2\theta = 1 - 1/R .
```

The order equations (README.md "Method", `order4.py`'s docstring), for k ≥ 0:

```math
\frac{dY_k}{dl} = F_{k-1} = \sum_{j<k} L_{k-1-j}\,\frac{l^j}{j!},\qquad (R-1)L_k' - \tfrac12 L_k = h_k,\qquad
h_k = (1 - 1/R)L_{k-1} + \sum_{j<k} Y_j c_{k-j} + \tfrac12 Y_k ,
```

with Y₀ = 1, h₀ = Y₀/2, and Y_k = 0 at l = 0 for k ≥ 1 (the mode is normalized by Y = 1 at the throat). The functions c_m are
the table of `cc` (order4.py, order5.py). Used here: c₁ = Ω₀ + 3l/2 − 1 + 1/R, and c_m(l = 0, R = 1) = Ω_{m−1}.
The solution analytic at the throat, and the condition at order k:

```math
L_k = -2\Omega_{k-1} + u\,C_k,\qquad C_k(\theta) = 2\int_0^\theta \frac{h_k - \Omega_{k-1}}{\sin^2\varphi}\,d\varphi,\qquad G_k = C_k(\pi/2) = 0 .
```

**Hypothesis H_k.** G_j = 0 for 2 ≤ j ≤ k − 1 (G₁ = 0 identically). That is, Ω₀ … Ω_{k−3} have their values. Ω_{k−2}
is free, and Ω_{k−1} cancels in G_k. Every identity below is an identity in Ω_{k−2}, and §7 tests it at three values.

**Four facts, used at every boundary.**

- **F1, the substitution.** R − 1 = tan²θ, dR = 2 sec²θ tan θ dθ, so (R − 1)^{−3/2}dR = 2 dθ/sin²θ, which is the second
  form of C_k above. Also dl/dθ = 2 tan θ, hence cot θ · dl/dθ = 2, and tan θ = −dℓ/dθ. *Checked symbolically* (S1, S2).
- **F2, the throat.** h_k(1) = Ω_{k−1}: at R = 1, l = 0, 1 − 1/R = 0, Y_j = δ_{j0}, and c_k(0, 1) = Ω_{k−1}. So
  h_k − Ω_{k−1} = (R − 1)a_k(R) with a_k analytic at R = 1 (induction: l = ln R, 1/R, Y_j, L_{k−1} are analytic there). Then
  C_k = 2∫₀^θ a_k R dφ = O(θ), L_k = −2Ω_{k−1} + O(θ²) is analytic in R − 1 = tan²θ, and Y_k = ∫₀^l F_{k−1}dl′ = O(l) = O(θ²)
  for k ≥ 1. At the throat ℓ = −θ²/2 + O(θ⁴) (S7).
- **F3, growth at θ → π/2** (x = π/2 − θ, cot θ ≈ x, tan θ ≈ 1/x, l ≈ 2 ln(1/x)). Under H_k, every Y_j, L_j, h_j, F_j with
  j ≤ k − 1 (and Y_k) is O(l^p) for some p. Induction: L₀ = −1. If the L_i, i < j, are O(l^p), then Y_i for i ≤ j, F_{j−1}
  and h_j are polynomially bounded in l. Because G_j = 0, C_j(θ) = −2∫_θ^{π/2}(h_j − Ω_{j−1})/sin²φ dφ, and with 1/sin²φ ≤ 2 on
  [π/4, π/2] and l(φ) increasing, |C_j| ≤ K x l^q, so |u C_j| ≤ K l^q. Hence L_j is O(l^q) and Y_{j+1} = ∫F_j dl is O(l^{q+1}).
  Without G_j = 0 the term u G_j ~ G_j/x grows like e^{l/2}, and every integral below diverges. *Numerical evidence*
  (GROW): Y_k ∝ l^{k+1}, L_k ∝ l^k at l = 256.
- **F4, convergence.** By F2 and F3 each integral below converges at both ends, so each split of an integral is legal.
  At θ → 0 the integrands are bounded (O(θ²)/sin²θ), at θ → π/2 they are polynomials in l times at most a bounded function,
  and ∫₀^{π/2} lⁿ dθ < ∞.

## 2. E1: Y_k out of G_k

By F1, (1 − 1/R)/sin²θ = 1, so with S_k = Σ_{j<k} Y_j c_{k−j}:

```math
G_k = 2\int_0^{\pi/2}\frac{h_k - \Omega_{k-1}}{S}\,d\theta
    = 2\int L_{k-1}\,d\theta + 2\int\frac{S_k - \Omega_{k-1}}{S}\,d\theta + \int\frac{Y_k}{S}\,d\theta .
```

Integrate the last term by parts with d(−cot θ) = dθ/sin²θ:

```math
\int_0^{\pi/2}\frac{Y_k}{\sin^2\theta}\,d\theta = \Big[-\cot\theta\,Y_k\Big]_0^{\pi/2} + \int_0^{\pi/2}\cot\theta\,\frac{dY_k}{d\theta}\,d\theta,
\qquad \frac{dY_k}{d\theta} = F_{k-1}\frac{dl}{d\theta} = 2\tan\theta\,F_{k-1} .
```

- **Boundary at θ = 0:** cot θ ≈ 1/θ and Y_k = O(θ²) (F2), so cot θ Y_k = O(θ). Its leading term is F_{k−1}(0)θ with
  F_{k−1}(0) = L_{k−1}(0) = −2Ω_{k−2}.
- **Boundary at θ = π/2:** cot θ ≈ x and Y_k = O(l^p) (F3), so cot θ Y_k = O(x l^p) → 0.
- The integrand is cot θ · 2 tan θ F_{k−1} = 2F_{k−1}.

So ∫Y_k/S = 2∫F_{k−1}dθ, and F_{k−1} = L_{k−1} + Σ_{j=1}^{k−1} L_{k−1−j} l^j/j!:

```math
\boxed{G_k = 4\int L_{k-1}\,d\theta + 2\int\sum_{j=1}^{k-1}L_{k-1-j}\frac{l^j}{j!}\,d\theta + 2\int\frac{S_k - \Omega_{k-1}}{\sin^2\theta}\,d\theta .}
```

This is E1. It holds under H_k, with Ω_{k−2} free. Consequence: Ω_{k−2} enters the right side only through
L_{k−1} = −2Ω_{k−2} + u C_{k−1} (C_{k−1} does not contain it) and through Y₀c_k + Y₁c_{k−1}, where it appears as
Ω_{k−2}(l + Y₁) = 0. So the slope of G_k in Ω_{k−2} is 4·(−2)(π/2) = −4π, the value stated in `order4.py`.

## 3. E2: ∫L_m by a second exchange, m = k − 1

Write L_m = −2Ω_{m−1} + tan θ C_m and use tan θ = −dℓ/dθ, C_m′ = 2(h_m − Ω_{m−1})/sin²θ:

```math
\int_0^{\pi/2}L_m\,d\theta = -\pi\Omega_{m-1} + \Big[-\ell\,C_m\Big]_0^{\pi/2} + 2\int_0^{\pi/2}\frac{\ell\,(h_m - \Omega_{m-1})}{\sin^2\theta}\,d\theta .
```

**Boundary terms of this step.**

- **θ = 0:** ℓ = O(θ²) and C_m = O(θ) (F2), so ℓ C_m = O(θ³).
- **θ = π/2:** this is where G_m = 0 is used. C_m(θ) = G_m − 2∫_θ^{π/2}(…) = −2∫_θ^{π/2}(h_m − Ω_{m−1})/sin²φ dφ = O(x l^q)
  (F3), and ℓ = −l/2, so ℓ C_m = O(x l^{q+1}) → 0. With G_m ≠ 0 the term would be −ℓ G_m, which diverges like ln(1/x), and
  so would ∫L_m (L_m ≈ G_m tan θ).

The remaining integral splits (F4) into 2∫h_m ℓ/sin² − 2Ω_{m−1}∫ℓ/sin². The antiderivative of ℓ/sin²θ is
A(θ) = −cot θ ℓ − θ (d/dθ: ℓ/sin²θ + 1 − 1, since d(−cot θ ℓ) = ℓ/sin²θ + cot θ tan θ). Its boundary values:

- **θ = 0:** cot θ ℓ ≈ −θ/2 → 0, so A(0) = 0.
- **θ = π/2:** cot θ ℓ ≈ x ln x → 0, so A(π/2) = −π/2.

So ∫₀^{π/2} ℓ/sin²θ dθ = −π/2 (S6), and −2Ω_{m−1}(−π/2) = +πΩ_{m−1} cancels −πΩ_{m−1}:

```math
\boxed{\int_0^{\pi/2}L_{m}\,d\theta = 2\int_0^{\pi/2}\frac{h_{m}\,\ln\cos\theta}{\sin^2\theta}\,d\theta .}
```

This is E2 at m = k − 1. It needs H_{m+1} through G_m = 0, nothing about Ω_{m−1}, and it holds for every constant Ω_{m−1} in
h_m. The boundary terms are: [−ℓ C_m] and [−cot θ ℓ], each zero at both ends. #26 checked them numerically only.

## 4. E3: Y_{k−1} out of G_k, with m = k − 1

Apply E1, then E2 at m to 4∫L_m = 8∫h_m ℓ/S. Split h_m = H_m + ½Y_m, with H_m = (1 − 1/R)L_{m−1} + Σ_{j<m}Y_j c_{m−j} free
of Y_m, and S_k = S_k′ + Y_m c₁, with S_k′ = Σ_{j<m}Y_j c_{k−j}. Y_m then enters in two places:

- from E2, through the ½Y_m of h_m: 8∫(½Y_m)ℓ/S = 4∫Y_m ℓ/S;
- from the j = m term of S_k: 2∫Y_m c₁/S.

By F1, c₁ = Ω₀ + 3l/2 − 1 + 1/R = Ω₀ − 3ℓ − S, so c₁/S = (Ω₀ − 3ℓ)/S − 1 (S3). The two terms add to

```math
\int Y_m\Big[\frac{4\ell + 2\Omega_0 - 6\ell}{S} - 2\Big]d\theta = -2\int Y_m\Big[\frac{\ell - \Omega_0}{\sin^2\theta} + 1\Big]d\theta = 2\int Y_m P'\,d\theta,
\qquad P = \cot\theta\,(\ell - \Omega_0) .
```

Here P′ = −[(ℓ − Ω₀)/sin²θ + 1] (S4, S5). Integrate by parts, with dY_m/dθ = 2 tan θ F_{m−1} and
P · 2 tan θ = 2(ℓ − Ω₀) (S4):

```math
2\int_0^{\pi/2} Y_m P'\,d\theta = 2\Big[Y_m P\Big]_0^{\pi/2} - 2\int P\,\frac{dY_m}{d\theta}\,d\theta = 2\Big[Y_m P\Big]_0^{\pi/2} - 4\int_0^{\pi/2}(\ell - \Omega_0)\,F_{m-1}\,d\theta .
```

**The boundary term [Y_m P]** (for k = 5, m = 4, the j = 4 boundary):

- **θ = 0:** ℓ = O(θ²), so P = −Ω₀/θ + O(θ), and Y_m = O(θ²) (F2): Y_m P = O(θ). Its leading term is
  −Ω₀F_{m−1}(0)θ = 2Ω₀Ω_{m−2}θ.
- **θ = π/2:** P = tan x (ℓ − Ω₀) = O(x l), and Y_m = O(l^{m+1}) (F3): Y_m P = O(x l^{m+2}), which is O(x l⁶) → 0 at m = 4.

Both ends vanish, and no Y_m is left:

```math
\boxed{G_k = 8\int\frac{H_m\,\ell}{\sin^2\theta}\,d\theta + 2\int\big(F_m - L_m\big)\,d\theta + 2\int\frac{S_k' - \Omega_{k-1}}{\sin^2\theta}\,d\theta - 4\int(\ell - \Omega_0)\,F_{m-1}\,d\theta,\qquad m = k - 1 .}
```

Also the by-parts identity by itself: ∫Y_m[(ℓ − Ω₀)/sin²θ + 1]dθ = 2∫(ℓ − Ω₀)F_{m−1}dθ. This is E3. It is valid for every
m ≥ 1 (the argument above does not use m = 4), so at k = 4 it removes Y₃ and at k = 5 it removes Y₄.

**At k = 5** (m = 4, F₃ = L₃ + L₂l + ½L₁l² + ⅙L₀l³, F₄ − L₄ = L₃l + ½L₂l² + ⅙L₁l³ + (1/24)L₀l⁴):

```math
G_5 = 8\int\frac{H_4\,\ell}{\sin^2\theta} + 2\int\Big(L_3 l + \tfrac12 L_2 l^2 + \tfrac16 L_1 l^3 + \tfrac1{24}L_0 l^4\Big)
 + 2\int\frac{\sum_{j\le3}Y_jc_{5-j} - \Omega_4}{\sin^2\theta} - 4\int(\ell - \Omega_0)\big(L_3 + L_2 l + \tfrac12 L_1 l^2 + \tfrac16 L_0 l^3\big).
```

Only Y₀ … Y₃ and L₀ … L₃ remain. Y₄, Y₅ and L₄ are gone. The slope of G₅ in Ω₃ is again visible: it comes only from
8∫(Y₀Ω₃)ℓ/S = −4πΩ₃, because Y₀Ω₃l + Y₁Ω₃ = Ω₃(l − l) = 0 in the third term.

```text
G_5 = 0
├─ E1   Y₅ → F₄                       parts on 1/sin²; boundary [−cot θ Y₅]
├─ E2   ∫L₄ → ∫h₄ ℓ/sin²             parts on tan θ; boundaries [−ℓ C₄], [−cot θ ℓ]; needs G₄ = 0
│   └─ h₄ contains ½Y₄
└─ E3   Y₄ → F₃                       parts against P; boundary [Y₄ P]; Y₄ from E2's h₄ and from h₅'s j = 4 term
remains: Y₀ … Y₃, L₀ … L₃
```

## 5. The Clausen integrals

With Cl₂(x) = Σ_{n≥1} sin(nx)/n² (Clausen), C(θ) = Cl₂(π − 2θ) = Σ (−1)ⁿ⁺¹ sin(2nθ)/n², using sin(nπ − 2nθ) =
(−1)ⁿ⁺¹ sin 2nθ. The series converges uniformly, since Σ1/n² < ∞. Cl₂′(x) = −ln(2 sin(x/2)), so
C′(θ) = −2Cl₂′(π − 2θ) = 2 ln(2 cos θ) (S9). C(0) = Cl₂(π) = 0 and C(π/2) = Cl₂(0) = 0.

**∫₀^{π/2} tan θ C dθ = π³/12, route 1 (by parts).** With tan θ = −dℓ/dθ:

```math
\int_0^{\pi/2}\tan\theta\,C\,d\theta = \Big[-\ell\,C\Big]_0^{\pi/2} + \int_0^{\pi/2}\ell\,C'\,d\theta = 2\int_0^{\pi/2}\ell\,(\ell + \ln 2)\,d\theta .
```

Boundary: at 0, ℓ = 0 and C = 0. At π/2, C = Cl₂(2x) = 2x(1 − ln 2x) + O(x³) and ℓ = ln sin x ≈ ln x, so ℓC = O(x ln²x) → 0 (S7).
The moments come from B(a) = ∫₀^{π/2}cos^aθ dθ = (√π/2)Γ((a + 1)/2)/Γ(a/2 + 1) by differentiating in a at a = 0:
∫ℓ = B′(0) = −(π/2) ln 2 and ∫ℓ² = B″(0) = (π/2)(ln²2 + π²/12) (S8). So

```math
2\Big[\tfrac\pi2\big(\ln^2 2 + \tfrac{\pi^2}{12}\big) - \tfrac\pi2\ln^2 2\Big] = \frac{\pi^3}{12} .
```

**Route 2 (the series).** I_n = ∫₀^{π/2}tan θ sin(2nθ)dθ. From sin 2(n + 1)θ − sin 2(n − 1)θ = 2 cos 2nθ sin 2θ and
tan θ sin 2θ = 2 sin²θ = 1 − cos 2θ:

```math
\tan\theta\,[\sin 2(n+1)\theta - \sin 2(n-1)\theta] = 2\cos 2n\theta - \cos 2(n+1)\theta - \cos 2(n-1)\theta .
```

For n ≥ 2 the right side integrates to 0 over [0, π/2], so I_{n+1} = I_{n−1}. Direct integrals give I₁ = ∫2 sin²θ = π/2 and
I₂ = ∫4 sin²θ cos 2θ = −π/2. So I_n = (−1)ⁿ⁺¹π/2 (S10). Termwise integration is legal by Tonelli: with x = π/2 − θ,
|tan θ sin 2nθ|/n² ≤ min(2/n, 1/(x n²)), whose sum over n is O(1 + ln(1/x)), integrable. Then
∫tan C = Σ (−1)ⁿ⁺¹/n² · (−1)ⁿ⁺¹π/2 = (π/2)ζ(2) = π³/12 (S11). The two routes share nothing but the definition of Cl₂.

**∫₀^{π/2} C² dθ = π⁵/360.** The functions sin(2nθ), n ≥ 1, are orthogonal on [0, π/2] with ∫sin²(2nθ) = π/4 (S10).
The series for C converges uniformly, so C² integrates termwise:

```math
\int_0^{\pi/2} C^2\,d\theta = \sum_{n\ge1}\frac{1}{n^4}\cdot\frac{\pi}{4} = \frac\pi4\,\zeta(4) = \frac{\pi^5}{360} .
```

This is Parseval's identity for the orthogonal system sin 2nθ, but only orthogonality and uniform convergence are used, not
completeness.

## 6. Checks

All from `exchange_review.py`; the declared run is raw `I209-rpi5/run_01_sym_exact_clausen_exch_grid.log` (39 of 39 pass, 203 s on
rpi5, code 886f764). "Digits" is −log₁₀ of the relative difference of two evaluations.

**Symbolic** (sympy, residual 0, *checked symbolically*): S1 the substitution and (R − 1)^{−3/2}dR = 2dθ/sin²θ; S2 1 − 1/R = sin²θ,
ln R = −2ℓ, dl/dθ = 2 tan θ; S3 c₁/sin²θ = (Ω₀ − 3ℓ)/sin²θ − 1; S4 P′ and P · dl/dθ; S5 the Y_{k−1} coefficient 4ℓ/S + 2c₁/S =
2P′; S6 the antiderivative of ℓ/sin²θ and its end values; S7 the boundary limits (cot θ ℓ at both ends, ℓθ, x ln⁶x, x ln²x);
S8 the Beta moments; S9 C′ = 2 ln(2 cos θ) and 2(∫ℓ² + ln 2 ∫ℓ) = π³/12; S10 orthogonality and I_n; S11 the two series sums.

**Closed forms, no panel code** (*numerical evidence*, mpmath at 50 digits, integrands at 150 digits):

| id | check | digits | bar |
|---|---|---|---|
| P0b | K(u) = ∫₀^u ln(1 + t²)/(1 + t²)dt against 2θ ln 2 − Cl₂(π − 2θ) at three θ | 49.5 | 40 |
| P0c | the table's h₁, h₂ against README's closed forms | 48.2 | 40 |
| P1 | E1 at k = 2, Ω₀ = 3/10: G₂ from the definition and E1's right side against −4πΩ₀ − 16π ln 2 | 50.0, 50.0 | 40 |
| P2_1 | E2 at m = 1: ∫L₁ = 2∫h₁ℓ/sin² (and against −πΩ₀ − 2π ln 2) | 50.0 (50.0) | 40 |
| P2_2 | E2 at m = 2 with the closed L₂ (Cl₂ in K) | 50.0 | 40 |

**The Clausen integrals and the moments** (mpmath at 60 digits, tanh-sinh, error estimates 1e-99):

| id | check | digits | bar |
|---|---|---|---|
| CM | ∫ℓ, ∫ℓ² against −(π/2) ln 2 and (π/2)(ln²2 + π²/12) | 60.0, 60.0 | 40 |
| CL1 | ∫₀^{π/2} tan θ Cl₂(π − 2θ)dθ = 2.583856390024985014623026255591782933519… against π³/12 | 60.0 | 40 |
| CL1a | the by-parts step: ∫tan C = 2∫ℓ ln(2 cos θ) | 60.0 | 40 |
| CL2 | ∫C² = 0.85005467995911514795205919456509890689… against π⁵/360 | 60.0 | 40 |

**E1, E2, E3 on the hierarchy** (method B's `PanelK`, E-grid N = 36, w = 1/2, V = 16, dps 80, 32 panels; its Ω₀ … Ω₂ meet the
closed forms to 44.2, 43.8, 43.7 digits; *numerical evidence*). Each identity is tested with the free constant at its root, at 0
and at 1/2, because at the root G_k = 0 and a relative digit count is empty:

| id | identity | at the root | shifted | bar |
|---|---|---|---|---|
| P0 | Y₁, L₁, Y₂, L₂ against README's closed forms, 576 nodes with l ≤ V²/4 | worst difference 8e-42 (L₁), 4e-43 (Y₂), 1e-33 (L₂) | | 1e-30 |
| E1_3 | G₃ = E1 (k = 3) | right side 6.5e-43 against terms of size 117 | 43.0 digits | 30 |
| E1_4 | G₄ = E1 (k = 4) | 4.3e-42 against 1.4e3 | 42.2 | 30 |
| E1_5 | G₅ = E1 (k = 5) | 5.7e-42 against 1.1e4 | 42.0 | 30 |
| E2_2 | ∫L₂ = 2∫h₂ℓ/sin² | 45.0 | 45.1, 45.2 (Ω₁ = 0, 1/2) | 30 |
| E2_3 | ∫L₃ = 2∫h₃ℓ/sin² | 43.8 | 43.3, 43.3 | 30 |
| E2_4 | ∫L₄ = 2∫h₄ℓ/sin² | 43.7 | 43.9, 43.9 | 30 |
| E3_3, E3_4, E3_5 | ∫Y_m[(ℓ − Ω₀)/sin² + 1] = 2∫(ℓ − Ω₀)F_{m−1}, m = 2, 3, 4 | | 43.9, 43.6, 44.9 | 30 |
| E3G_3, E3G_4, E3G_5 | G_k against the boxed form of §4 with no Y_{k−1}, no Y_k | | 42.7, 41.6, 40.4 | 30 |

The four terms of E3G_5 at Ω₃ + 1/2 are −5489.5031, 1292.1853, −60.078104 and 4251.1127; they sum to G₅ = −2π.

**Boundary terms** (BND, at the E-grid's first two nodes θ = 2.4e-4 and 2.1e-3, and at three nodes along l and the last node
x = π/2 − θ = 2.6e-56):

| term b | throat end, b/θᵖ at the two nodes | prediction | outer end, last node | b/(x l^q), l = 30, 60, 100 |
|---|---|---|---|---|
| Y₃ cot θ | −20.8467, −20.8467 (p = 1) | −2Ω₁ = −20.8467 | 2.2e-49 | 0.025, 0.011, 0.006 (q = 4) |
| Y₄ cot θ | −148.739, −148.739 | −2Ω₂ | 1.8e-46 | 0.053, 0.028, 0.017 (q = 5) |
| Y₅ cot θ | 1305.68, 1305.68 | −2Ω₃ | 1.8e-44 | 0.004, 0.007, 0.006 (q = 6) |
| ℓ C₂ | 8.31777, 8.31777 (p = 3) | O(θ³) | 3.2e-51 | 2.18, 2.09, 2.05 (q = 3) |
| ℓ C₃ | −46.6446, −46.6446 | O(θ³) | 8.1e-49 | 1.64, 1.83, 1.90 (q = 4) |
| ℓ C₄ | −107.509, −107.510 | O(θ³) | −1.1e-46 | −1.95, −1.46, −1.27 (q = 5) |
| cot θ ℓ | −0.5, −0.5 (p = 1) | −θ/2 | −3.3e-54 | −0.5, −0.5, −0.5 (q = 1) |
| Y₄ P | −412.391, −412.392 (p = 1) | 2Ω₀Ω₂ = −412.39 | −2.3e-44 | −0.022, −0.013, −0.008 (q = 6) |

The outer-end values of ℓC_m come from L_m at the last node, where the grid's L_m is not exact (§10, limit 2); the ratios at
l ≤ 100, where it is, carry the evidence. GROW (reported, no bar): Y_k/l^k at l = 256 is −1.0, −1.46, 0.51, 1.65, 0.63 for
k = 1 … 5, so Y_k grows no faster than l^k; the proof's bound l^{k+1} is not sharp.

**Remark.** The slope of G_k in Ω_{k−2} is −4π at every k (E1 alone gives it, §2). It is also the coefficient of Ω_{k−2} in
the boxed form of §4, where it comes from 8∫(Y₀Ω_{k−2})ℓ/S = −4πΩ_{k−2}. GR4 below finds it to 88.8 digits.

## 7. The third grid

Declared before the run (raw README): N = 64, w = 1/3, V = 24, dps 130, 72 panels. #26's grids are (40, 1/2, 18) and
(50, 2/5, 20), #208's (80, 3/10, 28) and (90, 1/4, 30). Run 98 s on rpi5, peak 100 MB, through the unchanged `order5.run_B`.

| check | against | digits | bar |
|---|---|---|---|
| GR1 Ω₂ = 74.369386727933317595384670295846617266734143415020337… | #26's A70 (`order4_output.txt`, 48 printed digits) | 50.3 | 40 |
| GR2 Ω₃ = −652.84010757005612543523930379964736108046174754557825736… | #208's B4 (`order5_output.txt`, 113 printed digits) | 87.8 | 40 |
| GR3 Ω₀, Ω₁ | −4 ln 2, 2π²/3 + 8 ln²2 | 87.8, 87.3 | 40 |
| GR4 slope of G_k in Ω_{k−2}, k = 2 … 5 | −4π | 88.8 (each) | 40 |
| reported: Ω₂, Ω₃ | 80ζ(3) − (8/3)π² ln 2 − (32/3) ln³2, #208's closed form | 87.2, 87.8 | |

GR1's 50.3 digits are set by A70's 48 printed digits.

## 8. Comparison with omega2.md and omega3.md

Made after §1 to §7 and the output were committed (branch commit a2e647e). I read `omega2.md` and `omega3.md` for the first
time then. The list below covers every step, boundary term, sign and constant of E1, E2, E3 and the Clausen steps in both notes.

**No difference in a sign, a constant or a result.** The two notes and this one agree on: E1 and its three integrals; E2's
statement and its hypothesis G_{k−1} = 0; c₁ = Ω₀ − 3ℓ − S; P = cot θ (ℓ − Ω₀) and P′ = (Ω₀ − ℓ − S)/S; the Y₄ terms
(4ℓ + 2c₁)Y₄/S = 2P′Y₄ (omega3.md step 2); the by-parts result −4∫(ℓ − Ω₀)F₃dθ; the factor 4 and the coefficient −32 of P′ for the
D₃ term of Y₃ (omega2.md step 1: δ = −32P′, equal to my 2P′ times Y₃'s −16D₃); C′ = 2ℓ + 2ln 2; ∫ℓ = −(π/2)ln 2 and
∫ℓ² = (π/2)(ln²2 + π²/12); ∫tan θ C = π³/12; ∫C² = π⁵/360; and the slope −4π.

**Differences in route or in what is stated** (none changes a result):

1. **E2's route.** omega2.md writes L = −2h + 2 tan θ B, B = ∫₀^θ h_φ cot φ dφ (one integration by parts inside C), then
   exchanges ∫tan θ B against ℓ. Here one integration by parts acts on C_m itself, with ∫ℓ/sin²θ = −π/2 closing the constants.
   Both need G_{k−1} = 0 at the same place: in omega2.md's route B(π/2) = G/2 must vanish for the term [−2ℓB] to go, here
   [−ℓC_m] goes only if G_m = 0.
2. **Boundary terms.** omega2.md states none for E1 or E2 and checks E2's numerically (N4a at k − 1 = 2, N4bg at 3); omega3.md
   says so as well ("not derived here: E2's boundary terms"). §2 and §3 derive them: E1's [−cot θ Y_k]; E2's [−ℓC_m] (O(θ³)
   at 0, O(x l^{q+1}) at π/2, which needs G_m = 0) and [−cot θ ℓ − θ]. omega3.md's remaining gap, "the boundary argument of E3
   at j = 4 (analytic, not checked at j = 4 numerically)", is checked here: [Y₄P] is O(θ) with coefficient 2Ω₀Ω₂ = −412.39 at
   the throat and ~x l⁶ at the outer end (BND), and E3_5 and E3G_5 agree to 44.9 and 40.4 digits.
3. **omega2.md step 1, the boundary term at θ = 0.** It reads "D₃(0) = 0, P(π/2) = 0, so the boundary terms vanish". P has a pole
   at 0 (P = −Ω₀/θ + O(θ)), so D₃(0) = 0 does not suffice: it needs D₃ = O(θ²), which holds since D₃′ = 2C and C = 2 ln 2 · θ +
   O(θ³) give D₃ = 2 ln 2 · θ² + O(θ⁴), and then D₃P = −2 ln 2 Ω₀ θ → 0. The conclusion stands (H2b, 45 digits). This is the one
   place where a written argument is incomplete; it is not an error.
4. **E3's scope.** omega3.md states E3 for Y₄ at k = 5. The argument holds for Y_{k−1} at every k ≥ 2 (§4), and E3_3, E3_4, E3_5 and
   E3G_3, E3G_4, E3G_5 check it at k = 3, 4, 5. #26's step 1 is the same by-parts step applied only to the D₃ part of Y₃; applying
   it to all of Y₃ would give another, equivalent, integrand.
5. **The slope −4π.** Both notes quote "G_k is affine in Ω_{k−2} with slope −4π" as a fact checked numerically. §2 derives it from
   E1 (and §4 from the boxed form: 8∫(Y₀Ω_{k−2})ℓ/S = −4πΩ_{k−2}, the other Ω_{k−2} terms cancelling since Y₁ = −l).
6. **∫C² = π⁵/360.** omega3.md uses Parseval's identity on [0, π], ½∫₀^π Cl₂² = (π/4)ζ(4). Here the substitution x = π − 2θ is
   not made: C = Σ(−1)^{n+1}sin(2nθ)/n² is integrated against itself on [0, π/2], which needs only orthogonality and uniform
   convergence, not completeness. The value is the same.
7. **Growth at π/2.** omega3.md has Y₄ = O(ℓ⁴). The induction of §1 proves only O(l^{k+1}). The grid shows Y_k/l^k bounded
   (Y₄/l⁴ = 1.6 at l = 256), so omega3.md's order is the right one; the boundary term is O(x l⁵) there instead of my O(x l⁶).

Nothing here disagrees with the repository's text or with another computation.

## 9. Controls

Every check has a `--mutate=` control (27 of them, `exchange_review.py --mutate=NAME`). Their targets were computed in
development runs (raw `dev_ctl_*`, `dev2_ctl_*`, `dev3_ctl_cl_moment.log`), declared in the raw README and in `TARGETS`, and
then run (`ctl_*`, serial, exit 1 each). The summary (`CT_*` in the output) reads the logs: **27 of 27 on target**, meaning every
declared target fired and no declared-silent check fired.

| control | breaks | fires | declared silent, and silent |
|---|---|---|---|
| sym_S1, sym_S2, sym_S6 … sym_S11, sym_P | one sympy claim each | the check named (S1, S2, S6, S7, S8, S9, S10, S11, S4) | S9 for sym_S8 |
| sym_c1 | c₁'s 3l/2 taken as l | S3, S5 | |
| exact_K | K = θ ln 2 − Cl₂ | P0b, P2_2 | P0c, P1, P2_1 |
| exact_L1 | L₁ = −2Ω₀ − l | P0c, P1, P2_1, P2_2 | P0b |
| exact_h2 | README's h₂ with l/10 | P0c | P1, P2_1, P2_2 |
| e1_sign, e1_factor | E1's sign, E1's 4 as 2 | P1, E1_3, E1_4, E1_5 | E2_*, E3G_* |
| e2_weight | E2's ℓ as ℓ + ln 2 | P2_1, P2_2, E2_2, E2_3, E2_4 | E1_*, E3G_* |
| bnd_G | Ω₀ + 1/100 | P0, BND | E1_*, E2_*, E3_*, E3G_* |
| bnd_Y0 | Y_k + 1e-3 | P0, BND, E1_*, E2_*, E3_*, E3G_* | |
| e3_Om0, e3_factor | P without Ω₀; the 4 as 2 | E3G_3, E3G_4, E3G_5 | E3_*, E1_*, E2_* |
| cl_moment | ℓ + ln 2 in the moments | CM | CL1, CL1a, CL2 |
| cl_arg | Cl₂(π − θ) | CL1, CL1a, CL2 | CM |
| cl_arg2 | Cl₂(2θ) | CL1, CL1a | CM, CL2 |
| grid_coarse, grid_short | N = 10; V = 6 | GR1 … GR4 | |
| grid_wrongh | the sign of −l/R in c₂ | GR1, GR2, GR3 | GR4 |
| grid_wrongOm1 | Ω₁ + 1/100 in the G₄ step | GR1, GR2 | GR3, GR4 |

Silences that are results:

- **bnd_G.** The shifted Ω₀ leaves E1, E2, E3 true: the grid's L_m has no growing term at any Ω₀ (§10, limit 2). Only the closed
  forms (P0) and the throat boundary test (BND: C₂(0) = −G₂) see it.
- **e3_Om0, e3_factor.** The by-parts identity E3_k holds with any constant in place of Ω₀ in P. Only the combination with c₁ in
  G_k fixes it to Ω₀, and E3G_k sees that.
- **cl_arg2.** ∫C² is blind to the sign pattern (−1)^{n+1}, because the orthogonality of sin 2nθ does not use it. The first
  integral sees it. `cl_arg` (the task's example, Cl₂(π − θ)) makes ∫tan θ Cl₂(π − θ) diverge at π/2, where Cl₂(π/2) = G ≠ 0
  (Catalan's constant); the quadrature returns a large number and CL1, CL1a fire.
- **grid_wrongh** leaves the slope of G_k at −4π (GR4): the slope comes from c_k's constant and Y₁ = −l, not from c₂'s l/R term.

## 10. Limits

1. The hierarchy of the E checks, and of the third grid, is method B's code, shared with #26 and #208. The grid changes the
   discretization, not the method. The closed-form checks P0, P1, P2 and the derivations of §2 to §5 do not use it.
2. `PanelK` builds L_m = −2Ω − u·(tail from the node to V). That is the solution without the growing term for any Ω, and it is
   the throat-analytic solution only when G_m = 0. E1, E2, E3 are therefore tested for that solution at every Ω₀, and the grid
   cannot test that G_m = 0 is needed. The hypothesis shows at the throat: with Ω₀ shifted by 1/100 (control bnd_G), C₂(0) =
   −G₂ ≠ 0 and ℓC₂/θ³ is no longer constant (272 and 38 at the first two nodes). L_m is also inexact at the nodes with
   e^{(l − V²)/2} not small (l > V²/4 = 64), because the tail beyond V is cut and u = e^{l/2} amplifies rounding. The θ-integrals
   carry the weight e^{−l/2} and are not affected.
3. F2 and F3 are asymptotic arguments (F3 an induction, sketched in §1), checked by the boundary table, not computer-assisted
   proofs.
4. Agreement of 40 to 45 digits is that of two evaluations on one grid, limited by the grid (the E-grid's Ω₀ … Ω₂ meet the closed
   forms to 44 digits). It is not an error bar.
5. The hypotheses of E1 to E3 are G_j = 0 for j < k (Ω₀ … Ω_{k−3} at their values) and Y_k = 0 at l = 0. The checks use the
   solved values of Ω₀ … Ω₂ from the same grid, not the closed forms, so they do not rely on #26's or #208's results.

## Files

- `exchange_review.py`, `exchange_review_output.txt`: the checks and the declared run with the controls' summary (66 checks, 198 s
  on rpi5; run `python3 exchange_review.py sym exact clausen exch grid ctlsum --ctldir=<raw folder> --ctlprefix=ctl` from this
  directory after the controls, `run_controls.sh` of the raw folder, have run).
- Raw `I209-rpi5` (negative-mass-raw): README with the declarations, development logs `dev*`, declared runs `run_*`, control runs
  `ctl_*`, `units.txt`.
