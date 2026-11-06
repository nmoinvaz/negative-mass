#!/bin/bash
# Step 2, queue S: queue Q's last four jobs (the other parity sector at D 20.04, 12, 10, 7), run in reverse order
# beside Q once P and R had finished; a job done here is skipped by Q (step2.sh's done_ check).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
CLD=/Users/nathan/Source/choptuik-large-d
RAW=/Users/nathan/Source/choptuik-large-d-raw
PY=$CLD/.venv/bin/python
export PYTHONDONTWRITEBYTECODE=1
L=$HERE/logs2; S=$CLD/hp/data/states
done_() { [ -f "$L/$1.log" ] && tail -1 "$L/$1.log" | grep -q "exit 0"; }
job() { local tag=$1 dir=$2; shift 2; done_ "$tag" && return 0; [ -f "$L/$tag.log" ] && return 0
        echo "== $tag $(date -u '+%Y-%m-%d %H:%M') UTC: $*" >> "$L/$tag.log"
        (cd "$dir" && nice -n 10 $PY -u "$@" >> "$L/$tag.log" 2>&1)
        echo "== $tag exit $? $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$L/$tag.log"; }
axis()  { local tag=$1 st=$2; shift 2; job "$tag" "$CLD/hp" modecount_test_axis.py "$st" --s=1 --threads=2 \
            --out=$HERE/axis_s1.jsonl --allow_unraw "$@"; }
axis ax_D20.04 $S/up2_D20.0437.npz --lo=-1.6 --hi=-0.02 --h=0.04
axis ax_D12 $S/up_D12.0000.npz --lo=-1.6 --hi=-0.02 --h=0.05
axis ax_D10 $S/up_D10.0000.npz --lo=-1.6 --hi=-0.02 --h=0.05
axis ax_D7 $RAW/T79-win/down_D7.0000.npz --lo=-1.6 --hi=-0.02 --h=0.05
