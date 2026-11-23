#!/usr/bin/env python3
"""Post-hoc control of the review (issue 104), written AFTER the control run showed `accept_zero_enclosure` silent in common, roots and brackets.
Why it was silent: the mutation accepted an enclosure with an end exactly at 0 (v.a >= 0), but the planted enclosure of P(a) straddles 0
(the coefficients of P are not binary numbers), so the mutated test still refused it.  The pre-registered target list is not changed.
Here the intended bug is planted instead: a sign test that takes the midpoint's sign when the enclosure contains 0.
The planted check "endpoint exactly on a root" must then fail (verdict CERTIFIED instead of REFUSED).  Exit 0 means the control fired on target.
"""
import sys

import review_large_l_common as C


def sgn_midpoint(v):
    s = C.sgn(v)
    if s != 0:
        return s
    m = v.mid
    return 1 if m.a >= 0 else -1


def main():
    C.sgn_checked = sgn_midpoint  # certify_bracket looks the function up in the module namespace
    chk = C.Checker(failfast=False)
    print("review_large_l_posthoc: planted problem with a sign test that accepts enclosures containing 0 (midpoint sign)")
    ok = C.run_planted(chk)
    fired = [n for n in chk.fails]
    target = "planted: endpoint exactly on a root"
    print(f"planted checks that failed under the bug: {fired}")
    on_target = (not ok) and any(target in n for n in fired)
    print(f"{'OK' if on_target else 'NOT ON TARGET'}: the bug is {'caught' if on_target else 'missed'} by the planted problem (target: {target})")
    sys.exit(0 if on_target else 1)


if __name__ == "__main__":
    main()
