"""Tiny check harness: every check prints PASS/FAIL, and the script exits non-zero if any failed."""
import sys


class Checker:
    def __init__(self):
        self.fail = 0
        self.n = 0

    def __call__(self, name, ok, detail=''):
        self.n += 1
        print(('PASS  ' if ok else 'FAIL  ') + name + (('  [' + detail + ']') if detail else ''), flush=True)
        if not ok:
            self.fail += 1
        return ok

    def done(self):
        print(f"\n{self.n - self.fail}/{self.n} checks passed" + ('' if not self.fail else f", {self.fail} FAILED"))
        sys.exit(1 if self.fail else 0)
