# Setup

## Software

| what | needed for | version |
|---|---|---|
| Python | every script | 3.14 |
| numpy, scipy, sympy, mpmath, matplotlib, numba | the scripts, figures and evolutions | pinned in [requirements.txt](requirements.txt) |
| git, with an SSH key on GitHub | cloning, and pushing over SSH | any recent |
| gh CLI | agents posting to issues and the project | any recent |
| a TeX distribution with `latexmk` and `pdflatex` | building the papers | not pinned |

The figures need no TeX, since matplotlib's bundled Computer Modern is used.

## Environment

Clone the repository, and the raw-output repository beside it if you will push run output:

```
git clone git@github.com:nmoinvaz/negative-mass.git
git clone git@github.com:nmoinvaz/negative-mass-raw.git
cd negative-mass
```

Create a virtual environment in `.venv` (ignored by git) and install the pinned packages.

macOS and Linux:

```
python3.14 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
```

Windows:

```
py -3.14 -m venv .venv
.venv\Scripts\python -m pip install --upgrade pip
.venv\Scripts\python -m pip install -r requirements.txt
```

The commands below write `python` for the environment's interpreter. For long runs, set `OPENBLAS_NUM_THREADS`,
`OMP_NUM_THREADS` and `MKL_NUM_THREADS` to 1 and size the processes to the machine's budget in
[AGENTS.md](AGENTS.md#machines).

## Check

Each `*_output.txt` is the stdout of the script beside it. This compares one script's output with the stored file and
exits 1 on any difference, in about 15 s. Run it from `phase3`:

```
python -c "import subprocess,sys; out=subprocess.run([sys.executable,'circular_check.py'],capture_output=True,text=True).stdout; sys.exit(out != open('circular_check_output.txt').read())"
```

On macOS (arm64, Python 3.14.5) these packages reproduce the stored outputs of `phase3/circular_check.py`,
`collinear.py`, `point_particle.py`, `test_body.py` and `phase5/spherical/make_stars.py` exactly, as well as
`paper/figures/thresholds.csv`. win2 runs numpy, scipy, sympy and mpmath at these versions on Python 3.14.10.

## Before some scripts

- `phase5/spherical`: run `python make_stars.py` once (about 2 min). It writes the stars that `run.py` loads into
  `stars/`, which git ignores.
- `paper`: the build steps are in [paper/README.md](paper/README.md).
- `notes/wormhole_next/ghost_mu.py` runs only under choptuik-large-d's `.venv`, which has hpkernel and the
  `modecount_cx` build, from a read-only checkout beside this one ([notes/wormhole_next.md](notes/wormhole_next.md)).

Agents also need gh logged in with the project scope. The steps are in [AGENTS.md](AGENTS.md#tasks-the-github-project).
