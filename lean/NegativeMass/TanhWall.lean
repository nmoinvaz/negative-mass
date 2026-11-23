import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp
import Mathlib.Analysis.Calculus.Deriv.MeanValue

/-!
# The tanh wall has no Type I point on D = 0 off the equator (#132, #134)

notes/warp/he_tanh_bound.md proves that at every root r > 0 of (r f')' = 0 of the tanh profile
f r = (tanh σ(r + R_b) - tanh σ(r - R_b)) / (2 tanh σR_b), r f' < 0 and F + r f' < 0 with F = 1 - f, for every
σR_b > 0. T2 of he_surface.md needs F + r f' > 0 there, so the tanh wall has no Type I point on D = 0 off the equator.

With w = 2σr and C = cosh 2σR_b:
* `m_pos`: m w = w (cosh w + 2) - 3 sinh w > 0 for w > 0, from m(0) = m'(0) = m''(0) = 0 and m''' = w sinh w.
* `sign_at_root`: the algebra at a root, in s = sinh w, c = cosh w (c² = s² + 1), with the root condition P + C Q = 0.
* `profile_cosh`: the tanh form equals (1 + C) / (cosh w + C).
* `tanh_wall_no_typeI`: the statement in r, with Mathlib's `deriv`.
-/

open Real

namespace NegativeMass.TanhWall

/-- If `g 0 = 0` and `g` has a derivative that is positive on `(0, ∞)`, then `g` is positive on `(0, ∞)`. -/
theorem pos_of_deriv_pos {g g' : ℝ → ℝ} (h0 : g 0 = 0) (hd : ∀ x, HasDerivAt g (g' x) x)
    (hpos : ∀ x, 0 < x → 0 < g' x) {w : ℝ} (hw : 0 < w) : 0 < g w := by
  have hmono : StrictMonoOn g (Set.Ici 0) := by
    apply strictMonoOn_of_deriv_pos (convex_Ici 0)
    · exact fun x _ => (hd x).continuousAt.continuousWithinAt
    · intro x hx
      rw [interior_Ici] at hx
      rw [(hd x).deriv]
      exact hpos x hx
  have h := hmono (Set.mem_Ici.mpr le_rfl) (Set.mem_Ici.mpr hw.le) hw
  rwa [h0] at h

/-- m'' = w cosh w - sinh w is positive for w > 0 (its derivative is w sinh w). -/
theorem m2_pos {w : ℝ} (hw : 0 < w) : 0 < w * cosh w - sinh w := by
  refine pos_of_deriv_pos (g := fun w => w * cosh w - sinh w)
    (g' := fun x => 1 * cosh x + x * sinh x - cosh x) (by simp) (fun x => ?_) (fun x hx => ?_) hw
  · exact ((hasDerivAt_id' x).mul (hasDerivAt_cosh x)).sub (hasDerivAt_sinh x)
  · have : 0 < x * sinh x := mul_pos hx (sinh_pos_iff.mpr hx)
    linarith

/-- m' = w sinh w - 2 (cosh w - 1) is positive for w > 0 (its derivative is m''). -/
theorem m1_pos {w : ℝ} (hw : 0 < w) : 0 < w * sinh w - 2 * (cosh w - 1) := by
  refine pos_of_deriv_pos (g := fun w => w * sinh w - 2 * (cosh w - 1))
    (g' := fun x => 1 * sinh x + x * cosh x - 2 * sinh x) (by simp) (fun x => ?_) (fun x hx => ?_) hw
  · exact ((hasDerivAt_id' x).mul (hasDerivAt_sinh x)).sub (((hasDerivAt_cosh x).sub_const 1).const_mul 2)
  · have := m2_pos hx
    linarith

/-- m w = w (cosh w + 2) - 3 sinh w is positive for w > 0 (its derivative is m'). -/
theorem m_pos {w : ℝ} (hw : 0 < w) : 0 < w * (cosh w + 2) - 3 * sinh w := by
  refine pos_of_deriv_pos (g := fun w => w * (cosh w + 2) - 3 * sinh w)
    (g' := fun x => 1 * (cosh x + 2) + x * sinh x - 3 * cosh x) (by simp) (fun x => ?_) (fun x hx => ?_) hw
  · exact ((hasDerivAt_id' x).mul ((hasDerivAt_cosh x).add_const 2)).sub ((hasDerivAt_sinh x).const_mul 3)
  · have := m1_pos hx
    linarith

/-- The sign at a root, in s = sinh w, c = cosh w. With P = s c + w (2 - c²), Q = s + w c and C > -1, a root of
P + C Q = 0 has F = (c - 1)/(c + C), g = r f' = -(1 + C) w s/(c + C)², g < 0 and F + g < 0. The certificate is
Q (A + C B) + w s (c - 1) m = B (P + C Q) + 2 w (c - 1) (c² - s² - 1), A = c² - c - w s, B = c - 1 - w s. -/
theorem sign_at_root {s c w C : ℝ} (hcs : c ^ 2 = s ^ 2 + 1) (hw : 0 < w) (hs : 0 < s) (hc : 1 < c)
    (hC : -1 < C) (hm : 0 < w * (c + 2) - 3 * s)
    (hroot : s * c + w * (2 - c ^ 2) + C * (s + w * c) = 0) :
    -(1 + C) * w * s / (c + C) ^ 2 < 0 ∧ (c - 1) / (c + C) + -(1 + C) * w * s / (c + C) ^ 2 < 0 := by
  have hD : 0 < c + C := by linarith
  have hD2 : 0 < (c + C) ^ 2 := by positivity
  have hQ : 0 < s + w * c := by nlinarith
  constructor
  · apply div_neg_of_neg_of_pos _ hD2
    have : 0 < (1 + C) * w * s := by
      have h1 : 0 < 1 + C := by linarith
      positivity
    linarith
  · set N := (c ^ 2 - c - w * s) + C * (c - 1 - w * s) with hN
    have key : (s + w * c) * N = -(w * s * (c - 1) * (w * (c + 2) - 3 * s)) := by
      rw [hN]
      linear_combination (c - 1 - w * s) * hroot + 2 * w * (c - 1) * hcs
    have hneg : N < 0 := by
      have hT : 0 < w * s * (c - 1) * (w * (c + 2) - 3 * s) := by
        have h1 : 0 < c - 1 := by linarith
        positivity
      by_contra h
      rw [not_lt] at h
      nlinarith [mul_nonneg hQ.le h]
    have hsplit : (c - 1) / (c + C) + -(1 + C) * w * s / (c + C) ^ 2 = N / (c + C) ^ 2 := by
      rw [hN]
      field_simp
      ring
    rw [hsplit]
    exact div_neg_of_neg_of_pos hneg hD2

/-- The tanh profile in cosh form: with w = 2y and C = cosh 2p, f = (1 + C) / (cosh w + C), which has no 0/0 at
p = 0. -/
theorem profile_cosh {y p : ℝ} (hp : 0 < p) :
    (tanh (y + p) - tanh (y - p)) / (2 * tanh p) = (1 + cosh (2 * p)) / (cosh (2 * y) + cosh (2 * p)) := by
  have hA : tanh (y + p) - tanh (y - p) = sinh (2 * p) / (cosh (y + p) * cosh (y - p)) := by
    rw [tanh_eq_sinh_div_cosh, tanh_eq_sinh_div_cosh, div_sub_div _ _ (cosh_pos _).ne' (cosh_pos _).ne', ← sinh_sub]
    congr 2
    ring
  have hB : cosh (y + p) * cosh (y - p) = (cosh (2 * y) + cosh (2 * p)) / 2 := by
    rw [cosh_add, cosh_sub, cosh_two_mul, cosh_two_mul]
    linear_combination (cosh p ^ 2 - 1 / 2) * cosh_sq y + (sinh y ^ 2 + 1 / 2) * cosh_sq p
  have hC : 1 + cosh (2 * p) = 2 * cosh p ^ 2 := by
    rw [cosh_two_mul]
    linear_combination (-1 : ℝ) * cosh_sq p
  have hsp : sinh p ≠ 0 := (sinh_pos_iff.mpr hp).ne'
  have hcp : cosh p ≠ 0 := (cosh_pos p).ne'
  have hden : cosh (2 * y) + cosh (2 * p) ≠ 0 := by
    have := cosh_pos (2 * y)
    have := cosh_pos (2 * p)
    positivity
  rw [hA, hB, hC, sinh_two_mul, tanh_eq_sinh_div_cosh]
  field_simp

/-- #132's result with Mathlib's derivatives: for the tanh profile f with σ > 0 and R_b > 0, every r > 0 where
(r f')' = 0 has r f' < 0 and F + r f' = 1 - f + r f' < 0. So T2 of he_surface.md, which needs F + r f' > 0 at a root of
T1: (r f')' = 0, never holds there, and the tanh wall has no Type I point on D = 0 off the equator. -/
theorem tanh_wall_no_typeI {σ R r : ℝ} {f : ℝ → ℝ} (hσ : 0 < σ) (hR : 0 < R) (hr : 0 < r)
    (hf : f = fun x => (tanh (σ * (x + R)) - tanh (σ * (x - R))) / (2 * tanh (σ * R)))
    (hroot : deriv (fun x => x * deriv f x) r = 0) :
    r * deriv f r < 0 ∧ 1 - f r + r * deriv f r < 0 := by
  set C := cosh (2 * (σ * R)) with hCdef
  have hC1 : 1 < C := one_lt_cosh.mpr (by positivity)
  have hDpos : ∀ x, 0 < cosh (2 * (σ * x)) + C := fun x => by
    have := cosh_pos (2 * (σ * x))
    linarith
  have hfc : f = fun x => (1 + C) / (cosh (2 * (σ * x)) + C) := by
    rw [hf]
    funext x
    rw [mul_add, mul_sub]
    exact profile_cosh (y := σ * x) (p := σ * R) (by positivity)
  have hw : ∀ x, HasDerivAt (fun x => 2 * (σ * x)) (2 * σ) x := fun x => by
    simpa using ((hasDerivAt_id x).const_mul σ).const_mul 2
  have hD : ∀ x, HasDerivAt (fun x => cosh (2 * (σ * x)) + C) (sinh (2 * (σ * x)) * (2 * σ)) x :=
    fun x => ((hw x).cosh).add_const C
  have hS : ∀ x, HasDerivAt (fun x => sinh (2 * (σ * x))) (cosh (2 * (σ * x)) * (2 * σ)) x :=
    fun x => (hw x).sinh
  have hderiv : ∀ x, deriv f x = (0 * (cosh (2 * (σ * x)) + C) - (1 + C) * (sinh (2 * (σ * x)) * (2 * σ)))
      / (cosh (2 * (σ * x)) + C) ^ 2 := fun x => by
    rw [hfc]
    exact ((hasDerivAt_const x (1 + C)).div (hD x) (hDpos x).ne').deriv
  have hKfun : (fun x => x * deriv f x) = fun x => x * ((0 * (cosh (2 * (σ * x)) + C)
      - (1 + C) * (sinh (2 * (σ * x)) * (2 * σ))) / (cosh (2 * (σ * x)) + C) ^ 2) := by
    funext x
    rw [hderiv x]
  -- (r f')' = -(1 + C) (2σ) P' / (cosh w + C)³ with P' = (s + w c)(c + C) - 2 w s², w = 2σr
  have hK : HasDerivAt (fun x => x * deriv f x)
      (-(1 + C) * (2 * σ) * ((sinh (2 * (σ * r)) + 2 * (σ * r) * cosh (2 * (σ * r))) * (cosh (2 * (σ * r)) + C)
        - 2 * (2 * (σ * r)) * sinh (2 * (σ * r)) ^ 2) / (cosh (2 * (σ * r)) + C) ^ 3) r := by
    rw [hKfun]
    have hnum := ((hasDerivAt_const r (0 : ℝ)).mul (hD r)).sub (((hS r).mul_const (2 * σ)).const_mul (1 + C))
    have hden := (hD r).pow 2
    have hE := hnum.div hden (pow_ne_zero 2 (hDpos r).ne')
    convert (hasDerivAt_id' r).mul hE using 1
    simp only [Pi.pow_apply, Pi.mul_apply, Pi.sub_apply]
    norm_num
    have := (hDpos r).ne'
    field_simp
    ring
  rw [hK.deriv] at hroot
  set w := 2 * (σ * r) with hwdef
  set s := sinh w with hsdef
  set c := cosh w with hcdef
  have hP' : (s + w * c) * (c + C) - 2 * w * s ^ 2 = 0 := by
    have h3 : (c + C) ^ 3 ≠ 0 := pow_ne_zero 3 (hDpos r).ne'
    have h1 : -(1 + C) * (2 * σ) ≠ 0 := by
      have : 0 < 1 + C := by linarith
      have : 0 < 2 * σ := by linarith
      intro h
      nlinarith
    rcases div_eq_zero_iff.mp hroot with h | h
    · rcases mul_eq_zero.mp h with h' | h'
      · exact absurd h' h1
      · exact h'
    · exact absurd h h3
  have hcs : c ^ 2 = s ^ 2 + 1 := cosh_sq w
  have hwpos : 0 < w := by positivity
  have hPQ : s * c + w * (2 - c ^ 2) + C * (s + w * c) = 0 := by
    linear_combination hP' - 2 * w * hcs
  obtain ⟨h1, h2⟩ := sign_at_root hcs hwpos (sinh_pos_iff.mpr hwpos) (one_lt_cosh.mpr hwpos.ne') (by linarith)
    (m_pos hwpos) hPQ
  have hrf : r * deriv f r = -(1 + C) * w * s / (c + C) ^ 2 := by
    rw [hderiv r]
    ring
  have hF : 1 - f r = (c - 1) / (c + C) := by
    rw [hfc]
    show 1 - (1 + C) / (c + C) = (c - 1) / (c + C)
    have hcC : c + C ≠ 0 := (hDpos r).ne'
    field_simp
    ring
  exact ⟨by rw [hrf]; exact h1, by rw [hF, hrf]; exact h2⟩

end NegativeMass.TanhWall
