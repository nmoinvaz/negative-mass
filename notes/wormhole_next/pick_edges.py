#!/usr/bin/env python3
"""Left edge aLL for a winding count that starts just right of the real-axis zero mu (notes/wormhole_next.md section 6).

    python3 pick_edges.py --hp=<choptuik-large-d>/hp STATE.npz MU [--lo=0.02] [--hi=0.08]

Picks aLL in (MU + lo, MU + hi) farthest from the real parts of the candidate poles of det M_psi (CLD's
hp/modecount_common.horizon_poles) with |Im| < 1.5 omega + 3 omega, as CLD's modecount_wind.far does, and prints it with
its distance to the nearest pole real part and to MU. Read-only on CLD.
"""
import sys

import numpy as np


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    sys.path.insert(0, opts["hp"])
    from modecount_common import setup, horizon_poles  # noqa: E402
    st = setup(args[0], threads=1, log=lambda m: None)
    mu = float(args[1])
    lo, hi = float(opts.get("lo", 0.02)), float(opts.get("hi", 0.08))
    y2 = 1.5 * st.omega
    pre = sorted(set(p.real for _, p in horizon_poles(st) if abs(p.imag) < y2 + 3 * st.omega))
    cands = np.linspace(mu + lo, mu + hi, 201)
    d = [min([abs(x - a) for a in pre] + [9.0]) for x in cands]
    i = int(np.argmax(d))
    print(f"{args[0]} D {st.D:.4f} omega {st.omega:.6f}: aLL {cands[i]:.4f} (pole real parts {d[i]:.4f} away, "
          f"{cands[i] - mu:.4f} right of mu {mu:+.6f}); nearest pole real parts {[round(a, 4) for a in pre if abs(a - cands[i]) < 0.3]}")


if __name__ == "__main__":
    main()
