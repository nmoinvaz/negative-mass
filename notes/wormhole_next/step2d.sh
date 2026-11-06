#!/bin/bash
# Step 2, queue T: D 4's count to the right of CLD's boxes (CLD counted the 11-node layout only for Re mu < 0.193),
# in the same layout and cache as box_D4, with a planted pair. Same conventions as step2.sh.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
CLD=/Users/nathan/Source/choptuik-large-d
PY=$CLD/.venv/bin/python
export PYTHONDONTWRITEBYTECODE=1
L=$HERE/logs2; C=$HERE/counts
done_() { [ -f "$L/$1.log" ] && tail -1 "$L/$1.log" | grep -q "exit 0"; }
job() { local tag=$1 dir=$2; shift 2; done_ "$tag" && return 0
        echo "== $tag $(date -u '+%Y-%m-%d %H:%M') UTC: $*" >> "$L/$tag.log"
        (cd "$dir" && nice -n 10 $PY -u "$@" >> "$L/$tag.log" 2>&1)
        echo "== $tag exit $? $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$L/$tag.log"; }
count() { local tag=$1 st=$2; shift 2; job "$tag" "$CLD/hp" modecount_test.py "$st" --threads=2 --data=$C --allow_unraw "$@"; }
S4=$CLD/hp/data/lowD_states/lowDn_D4.0000_M48.npz
count full_D4 $S4 --boxes=full --aLL=-0.3451 --tag=D4.0000_M48n
count full_D4_plant $S4 --boxes=full --aLL=-0.3451 --tag=D4.0000_M48n --plant=2.0+1.0j
