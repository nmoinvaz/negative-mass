# Lean proofs for negative-mass

Lean 4 proofs of the project's lemmas, checked by the kernel. Toolchain and library are pinned as in choptuik-large-d's
`lean/`: Lean v4.34.1 (`lean-toolchain`) and Mathlib v4.34.1 (`lakefile.toml`, `lake-manifest.json`). The build sets
`warningAsError = true`, so a proof left as `sorry` fails it.

## Build

From this directory, with elan installed:

    lake exe cache get    # once, after a fresh clone: Mathlib's prebuilt files
    lake build            # exit 0: every theorem checked
    lake env lean scripts/Axioms.lean    # the axioms each theorem uses
    python3 scripts/mutate.py            # controls: each mutated statement must fail to build

## Contents

| file | theorem | status in the notes |
|---|---|---|
| `NegativeMass/TanhWall.lean` | #132's result (#134): for the tanh profile f with σ, R_b > 0, every r > 0 with (r f′)′ = 0 has r f′ < 0 and 1 − f + r f′ < 0 (`tanh_wall_no_typeI`); with `m_pos`, `sign_at_root` and `profile_cosh` | **certified** (Lean kernel, Mathlib v4.34.1), [notes/warp/he_tanh_bound.md](../notes/warp/he_tanh_bound.md) |

Each theorem depends only on `propext`, `Classical.choice` and `Quot.sound` (`scripts/Axioms.lean`). `scripts/mutate.py`
changes one statement at a time and expects the build of the copy to fail: four mutations make a statement false, and
`m-plus3` makes a true statement whose proof no longer fits.
