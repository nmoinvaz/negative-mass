"""
Controls for the Lean files: each copies a file, changes one statement, and expects `lake env lean` to reject the
copy. Each unmutated file must compile. Run from lean/: python3 scripts/mutate.py [name ...]
Exit 1 if a baseline fails or any mutated copy compiles.
"""
import os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TANH = os.path.join(HERE, "NegativeMass", "TanhWall.lean")
SPLIT = os.path.join(HERE, "NegativeMass", "MinimalSplit.lean")
MUTATIONS = [
    # (file, name, what it claims, old text, new text)
    (TANH, "m-plus3", "m_pos with + 3 sinh w: a true claim; fails only because the proof no longer fits (issue's control)",
     "theorem m_pos {w : ℝ} (hw : 0 < w) : 0 < w * (cosh w + 2) - 3 * sinh w := by",
     "theorem m_pos {w : ℝ} (hw : 0 < w) : 0 < w * (cosh w + 2) + 3 * sinh w := by"),
    (TANH, "m-minus4", "m_pos with - 4 sinh w: false (about -w near 0)",
     "theorem m_pos {w : ℝ} (hw : 0 < w) : 0 < w * (cosh w + 2) - 3 * sinh w := by",
     "theorem m_pos {w : ℝ} (hw : 0 < w) : 0 < w * (cosh w + 2) - 4 * sinh w := by"),
    (TANH, "profile-two", "profile_cosh with 2 + cosh 2p in the numerator: false",
     "(tanh (y + p) - tanh (y - p)) / (2 * tanh p) = (1 + cosh (2 * p)) / (cosh (2 * y) + cosh (2 * p)) := by",
     "(tanh (y + p) - tanh (y - p)) / (2 * tanh p) = (2 + cosh (2 * p)) / (cosh (2 * y) + cosh (2 * p)) := by"),
    (TANH, "root-positive", "sign_at_root claiming 0 < F + g: false",
     "∧ (c - 1) / (c + C) + -(1 + C) * w * s / (c + C) ^ 2 < 0 := by",
     "∧ 0 < (c - 1) / (c + C) + -(1 + C) * w * s / (c + C) ^ 2 := by"),
    (TANH, "wall-positive", "tanh_wall_no_typeI claiming 0 < 1 - f r + r f' r: false",
     "r * deriv f r < 0 ∧ 1 - f r + r * deriv f r < 0 := by",
     "r * deriv f r < 0 ∧ 0 < 1 - f r + r * deriv f r := by"),
    (SPLIT, "split-58", "split_attained with #58's first formula max(0, (|J| - mu)/2): false where -mu > |J|",
     "μm = max 0 (max (-μ) ((‖J‖ - μ) / 2)) ∧",
     "μm = max 0 ((‖J‖ - μ) / 2) ∧"),
    (SPLIT, "lower-nohalf", "split_lower without the /2: false (mu = 0, |J| = 1 has a split with mu- = 1/2)",
     "max 0 (max (-μ) ((‖J‖ - μ) / 2)) ≤ μm := by",
     "max 0 (max (-μ) (‖J‖ - μ)) ≤ μm := by"),
    (SPLIT, "sum-plus", "split_attained with mu- + |J-| = max(0, |J| + mu): false (mu = 1, J = 0)",
     "μm + ‖Jm‖ = max 0 (‖J‖ - μ) := by",
     "μm + ‖Jm‖ = max 0 (‖J‖ + μ) := by"),
    (SPLIT, "least-58", "minimal_split with #58's first formula as the least element: false where -mu > |J|",
     "    IsLeast {μm : ℝ | ∃ μp : ℝ, ∃ Jp Jm : E, IsSplit μ J μp μm Jp Jm}\n      (max 0 (max (-μ) ((‖J‖ - μ) / 2))) := by",
     "    IsLeast {μm : ℝ | ∃ μp : ℝ, ∃ Jp Jm : E, IsSplit μ J μp μm Jp Jm}\n      (max 0 ((‖J‖ - μ) / 2)) := by"),
]


def check(text):
    with tempfile.NamedTemporaryFile("w", suffix=".lean", delete=False, dir=tempfile.gettempdir()) as fh:
        fh.write(text)
        path = fh.name
    try:
        r = subprocess.run(["lake", "env", "lean", path], cwd=HERE, capture_output=True, text=True)
        first = next((ln for ln in (r.stdout + r.stderr).splitlines() if "error" in ln), "")
        return r.returncode, first.replace(path, "<copy>")
    finally:
        os.remove(path)


only = set(sys.argv[1:])
todo = [m for m in MUTATIONS if not only or m[1] in only]
assert not only or len(todo) == len(only), "unknown mutation name"
srcs = {f: open(f).read() for f in dict.fromkeys(m[0] for m in todo)}
bad = []
for f, src in srcs.items():
    code, msg = check(src)
    print(f"baseline {os.path.basename(f)}: exit {code} {'(compiles)' if code == 0 else 'FAILED: ' + msg}")
    if code != 0:
        bad.append("baseline " + os.path.basename(f))
for f, name, what, old, new in todo:
    src = srcs[f]
    assert src.count(old) == 1, name
    code, msg = check(src.replace(old, new))
    fired = code != 0
    print(f"{'FIRED ' if fired else 'SILENT'} {name}: {what}; exit {code}; {msg[:160]}")
    if not fired:
        bad.append(name)
print("RESULT", "FAIL: " + ", ".join(bad) if bad else "PASS")
sys.exit(1 if bad else 0)
