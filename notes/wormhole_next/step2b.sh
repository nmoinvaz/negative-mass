#!/bin/bash
# Step 2, queue R (after queue Q): mu(D) on the large-D pinned states D 68.5-78.1 (choptuik-large-d-raw T84-mac, T75,
# T85-mac), a second layout at D 74-76 (T85's x85), the other parity sector at D 52 and 60, and M 48 at D 30.07.
# Same conventions as step2.sh. usage: step2b.sh R PID
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
CLD=/Users/nathan/Source/choptuik-large-d
RAW=/Users/nathan/Source/choptuik-large-d-raw
PY=$CLD/.venv/bin/python
export PYTHONDONTWRITEBYTECODE=1
L=$HERE/logs2
mkdir -p "$L"
done_() { [ -f "$L/$1.log" ] && tail -1 "$L/$1.log" | grep -q "exit 0"; }
job() { local tag=$1 dir=$2; shift 2; done_ "$tag" && return 0
        echo "== $tag $(date -u '+%Y-%m-%d %H:%M') UTC: $*" >> "$L/$tag.log"
        (cd "$dir" && nice -n 10 $PY -u "$@" >> "$L/$tag.log" 2>&1)
        echo "== $tag exit $? $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$L/$tag.log"; }
mu()    { local tag=$1; shift; job "$tag" "$HERE" ghost_mu.py --hp=$CLD/hp --threads=2 --tag=$tag --out=$HERE/ghost_mu2.jsonl "$@"; }
axis()  { local tag=$1 st=$2; shift 2; job "$tag" "$CLD/hp" modecount_test_axis.py "$st" --s=1 --threads=2 \
            --out=$HERE/axis_s1.jsonl --allow_unraw "$@"; }
while kill -0 "$2" 2>/dev/null; do sleep 30; done
W="--lo=-1.3 --hi=-0.04 --h=0.06"
mu D68.50_m1 $RAW/T84-mac/pin7_D68.5000_full.npz $W
mu D70.50_m1 $RAW/T84-mac/pin7_D70.5000_full.npz $W
mu D72.00_m1 $RAW/T75/pin7_D72.0000.npz $W
mu D73.27_m1 $RAW/T75/pin7_D73.2657.npz $W
mu D74.00_m1 $RAW/T85-mac/pin7_x80_D74.0000.npz $W
mu D76.00_m1 $RAW/T85-mac/pin7_x80_D76.0000.npz $W
mu D78.10_m1 $RAW/T85-mac/pin7_x80_D78.1049.npz $W
mu D74.00_x85_m1 $RAW/T85-mac/pin7_x85_D74.0000.npz --lo=-1.05 --hi=-0.93 --h=0.06
mu D76.00_x85_m1 $RAW/T85-mac/pin7_x85_D76.0000.npz --lo=-1.05 --hi=-0.93 --h=0.06
mu D30.07_M48_m1 $CLD/hp/data/modecount_state_D30.0688_M48.npz --lo=-0.98 --hi=-0.92 --h=0.06
axis ax_D52.04 $CLD/hp/data/states/up9_D52.0414.npz --lo=-1.6 --hi=-0.02 --h=0.04
axis ax_D60.00 $CLD/hp/data/states/up13_D60.0020.npz --lo=-1.6 --hi=-0.02 --h=0.04
axis ax_D78.10 $RAW/T85-mac/pin7_x80_D78.1049.npz --lo=-1.6 --hi=-0.02 --h=0.04
