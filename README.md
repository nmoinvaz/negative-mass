# negative-mass

Working notes and code on negative (active) mass in general relativity: Bondi's 1957 runaway, the
stability of negative-mass Schwarzschild, regular negative-mass sources, and what remains open.

- `PLAN.md` — phases, status, next steps
- `notes/literature.md` — papers read, with what each establishes
- `phase0/` — verification of the published perturbation results (sympy / mpmath)

## Setup
```
pip install sympy mpmath numpy scipy
python phase0/verify_zerilli.py   # several minutes (symbolic limits)
python phase0/mode_width.py
```
