"""
Controls for NegativeMass/TanhWall.lean (#134): each copies the file, changes one statement, and expects
`lake env lean` to reject the copy. The unmutated copy must compile. Run from lean/: python3 scripts/mutate.py
Exit 1 if the baseline fails or any mutated copy compiles.
"""
import os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(HERE, "NegativeMass", "TanhWall.lean")
MUTATIONS = [
    # (name, what it claims, old text, new text)
    ("m-plus3", "m_pos with + 3 sinh w: a true claim; fails only because the proof no longer fits (issue's control)",
     "theorem m_pos {w : ℝ} (hw : 0 < w) : 0 < w * (cosh w + 2) - 3 * sinh w := by",
     "theorem m_pos {w : ℝ} (hw : 0 < w) : 0 < w * (cosh w + 2) + 3 * sinh w := by"),
    ("m-minus4", "m_pos with - 4 sinh w: false (about -w near 0)",
     "theorem m_pos {w : ℝ} (hw : 0 < w) : 0 < w * (cosh w + 2) - 3 * sinh w := by",
     "theorem m_pos {w : ℝ} (hw : 0 < w) : 0 < w * (cosh w + 2) - 4 * sinh w := by"),
    ("profile-two", "profile_cosh with 2 + cosh 2p in the numerator: false",
     "(tanh (y + p) - tanh (y - p)) / (2 * tanh p) = (1 + cosh (2 * p)) / (cosh (2 * y) + cosh (2 * p)) := by",
     "(tanh (y + p) - tanh (y - p)) / (2 * tanh p) = (2 + cosh (2 * p)) / (cosh (2 * y) + cosh (2 * p)) := by"),
    ("root-positive", "sign_at_root claiming 0 < F + g: false",
     "∧ (c - 1) / (c + C) + -(1 + C) * w * s / (c + C) ^ 2 < 0 := by",
     "∧ 0 < (c - 1) / (c + C) + -(1 + C) * w * s / (c + C) ^ 2 := by"),
    ("wall-positive", "tanh_wall_no_typeI claiming 0 < 1 - f r + r f' r: false",
     "r * deriv f r < 0 ∧ 1 - f r + r * deriv f r < 0 := by",
     "r * deriv f r < 0 ∧ 0 < 1 - f r + r * deriv f r := by"),
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


src = open(SRC).read()
bad = []
code, msg = check(src)
print(f"baseline: exit {code} {'(compiles)' if code == 0 else 'FAILED: ' + msg}")
if code != 0:
    bad.append("baseline")
for name, what, old, new in MUTATIONS:
    assert src.count(old) == 1, name
    code, msg = check(src.replace(old, new))
    fired = code != 0
    print(f"{'FIRED ' if fired else 'SILENT'} {name}: {what}; exit {code}; {msg[:160]}")
    if not fired:
        bad.append(name)
print("RESULT", "FAIL: " + ", ".join(bad) if bad else "PASS")
sys.exit(1 if bad else 0)
