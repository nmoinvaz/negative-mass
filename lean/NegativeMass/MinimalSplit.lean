import Mathlib.Analysis.InnerProductSpace.PiL2

/-!
# The minimal split of a signed source (#58, #174)

notes/warp/mass_bounds.md item 3: write T = T₊ − T₋ with μ± ≥ |J±|, where μ is the Eulerian energy density and J the
momentum density at a point. Over all such splits the least μ₋ is max(0, −μ, (|J| − μ)/2). It is attained with
J₋ = −μ₋ Ĵ when |J| ≥ −μ and with J₋ = −J otherwise, and along that split μ₋ + |J₋| = max(0, |J| − μ).

The statements hold for J in any real normed space; `minimal_split_R3` is the case ℝ³ of the note.
* `split_lower`: every split has μ₋ ≥ max(0, −μ, (‖J‖ − μ)/2), from μ₊ = μ + μ₋ ≥ ‖J₊‖ ≥ 0 and
  ‖J‖ ≤ ‖J₊‖ + ‖J₋‖ ≤ μ + 2μ₋.
* `split_attained`: the split named above attains the bound and has μ₋ + ‖J₋‖ = max(0, ‖J‖ − μ).
* `minimal_split`: the bound is the least element of the set of admissible μ₋.
-/

namespace NegativeMass.MinimalSplit

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- (μ₊, μ₋, J₊, J₋) splits (μ, J) into two parts with μ± ≥ ‖J±‖. -/
def IsSplit (μ : ℝ) (J : E) (μp μm : ℝ) (Jp Jm : E) : Prop :=
  μ = μp - μm ∧ J = Jp - Jm ∧ ‖Jp‖ ≤ μp ∧ ‖Jm‖ ≤ μm

omit [NormedSpace ℝ E] in
/-- Every split has μ₋ ≥ max(0, −μ, (‖J‖ − μ)/2). -/
theorem split_lower {μ μp μm : ℝ} {J Jp Jm : E} (h : IsSplit μ J μp μm Jp Jm) :
    max 0 (max (-μ) ((‖J‖ - μ) / 2)) ≤ μm := by
  obtain ⟨hμ, hJ, hp, hm⟩ := h
  have h1 : 0 ≤ ‖Jm‖ := norm_nonneg _
  have h2 : 0 ≤ ‖Jp‖ := norm_nonneg _
  have h3 : ‖J‖ ≤ ‖Jp‖ + ‖Jm‖ := by rw [hJ]; exact norm_sub_le Jp Jm
  exact max_le (by linarith) (max_le (by linarith) (by linarith))

/-- The split J₋ = −μ₋ Ĵ (‖J‖ ≥ −μ) or J₋ = −J (‖J‖ < −μ) attains the bound, with μ₋ + ‖J₋‖ = max(0, ‖J‖ − μ). -/
theorem split_attained (μ : ℝ) (J : E) :
    ∃ μp μm : ℝ, ∃ Jp Jm : E, IsSplit μ J μp μm Jp Jm ∧
      μm = max 0 (max (-μ) ((‖J‖ - μ) / 2)) ∧ μm + ‖Jm‖ = max 0 (‖J‖ - μ) := by
  have hn : 0 ≤ ‖J‖ := norm_nonneg J
  by_cases hc : -μ ≤ ‖J‖
  · set m := max 0 ((‖J‖ - μ) / 2) with hm
    have hm0 : 0 ≤ m := le_max_left _ _
    have hmJ : (‖J‖ - μ) / 2 ≤ m := le_max_right _ _
    -- ‖(m/‖J‖) • J‖ = m, and ‖J − (m/‖J‖) • J‖ = |‖J‖ − m|; at J = 0 both sides vanish since then m = 0
    have hnorms : ‖(m / ‖J‖) • J‖ = m ∧ ‖J - (m / ‖J‖) • J‖ = |‖J‖ - m| := by
      rcases eq_or_ne J 0 with hJ0 | hJ0
      · have hz : ‖J‖ = 0 := by simp [hJ0]
        have hm00 : m = 0 := by rw [hm, hz]; exact max_eq_left (by linarith)
        simp [hJ0, hm00]
      · have hpos : 0 < ‖J‖ := norm_pos_iff.mpr hJ0
        have hc1 : m / ‖J‖ * ‖J‖ = m := by field_simp
        constructor
        · rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg (div_nonneg hm0 hn), hc1]
        · have hs : J - (m / ‖J‖) • J = (1 - m / ‖J‖) • J := by rw [sub_smul, one_smul]
          rw [hs, norm_smul, Real.norm_eq_abs, ← abs_of_nonneg hn, ← abs_mul, abs_of_nonneg hn]
          congr 1
          rw [sub_mul, one_mul, hc1]
    obtain ⟨hJm, hJp⟩ := hnorms
    refine ⟨μ + m, m, J - (m / ‖J‖) • J, -((m / ‖J‖) • J), ⟨by ring, by abel, ?_, ?_⟩, ?_, ?_⟩
    · rw [hJp, abs_le]; constructor <;> linarith
    · rw [norm_neg, hJm]
    · rw [max_eq_right (show -μ ≤ (‖J‖ - μ) / 2 by linarith)]
    · rw [norm_neg, hJm]
      rcases le_total (‖J‖ - μ) 0 with h | h
      · have e1 : m = 0 := by rw [hm]; exact max_eq_left (by linarith)
        rw [e1, max_eq_left h]; ring
      · have e1 : m = (‖J‖ - μ) / 2 := by rw [hm]; exact max_eq_right (by linarith)
        rw [e1, max_eq_right h]; ring
  · replace hc : ‖J‖ < -μ := not_le.mp hc
    refine ⟨0, -μ, 0, -J, ⟨by ring, by simp, by simp, by rw [norm_neg]; linarith⟩, ?_, ?_⟩
    · rw [max_eq_left (show (‖J‖ - μ) / 2 ≤ -μ by linarith), max_eq_right (by linarith)]
    · rw [norm_neg, max_eq_right (by linarith)]; ring

/-- max(0, −μ, (‖J‖ − μ)/2) is the least μ₋ over all splits of (μ, J). -/
theorem minimal_split (μ : ℝ) (J : E) :
    IsLeast {μm : ℝ | ∃ μp : ℝ, ∃ Jp Jm : E, IsSplit μ J μp μm Jp Jm}
      (max 0 (max (-μ) ((‖J‖ - μ) / 2))) := by
  refine ⟨?_, fun x ⟨_, _, _, h⟩ => split_lower h⟩
  obtain ⟨μp, μm, Jp, Jm, h, he, -⟩ := split_attained μ J
  rw [← he]
  exact ⟨μp, Jp, Jm, h⟩

/-- The note's case: J ∈ ℝ³ with the Euclidean norm. -/
theorem minimal_split_R3 (μ : ℝ) (J : EuclideanSpace ℝ (Fin 3)) :
    IsLeast {μm : ℝ | ∃ μp : ℝ, ∃ Jp Jm : EuclideanSpace ℝ (Fin 3), IsSplit μ J μp μm Jp Jm}
      (max 0 (max (-μ) ((‖J‖ - μ) / 2))) :=
  minimal_split μ J

end NegativeMass.MinimalSplit
