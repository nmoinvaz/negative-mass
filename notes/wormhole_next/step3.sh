#!/bin/bash
# Step 3 of notes/wormhole_next.md (issue #30, notes/wormhole_next/step3.md): winding counts of the ghost sector at
# D = 10, 12 and 20.04, and the other parity sector's zero at D = 10 and 12 refined on the line Im mu = omega, on
# choptuik-large-d's code and states, read-only. Two queues of 2 threads at nice 10; every job under a 3 h limit.
# Criteria, boxes, edges and planted pairs are in RAW/README.md, written before the first run. A job whose log ends in
# "exit 0" is skipped, so a restart repeats only unfinished jobs; modecount_test.py keeps every evaluation in
# counts/*_raw.jsonl and reuses it. A log is never overwritten (a rerun appends, or gets a number suffix).
# usage: step3.sh A|B|R2     (R2: rung 2 of the refinement, run only if rung 1 failed C4)
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
CLD=${CLD:-$(cd "$HERE/../../.." && pwd)/choptuik-large-d}
RAW=${RAW:-$(cd "$HERE/../../.." && pwd)/negative-mass-raw/I30-mac}
PY=$CLD/.venv/bin/python
SO=$CLD/hp/modecount_cx/build/modecount_cx.cpython-313-darwin.so
export PYTHONDONTWRITEBYTECODE=1
L=$RAW/logs; C=$RAW/counts
mkdir -p "$L" "$C"
done_() { [ -f "$L/$1.log" ] && tail -1 "$L/$1.log" | grep -q "exit 0"; }
# CLD's uncommitted changes to the scripts or the kernel stop a job (exit 99): ask the coordinator
dirty() { git -C "$CLD" status --short hp | grep -E "modecount_(test|test_axis|common|wind)\.py|modecount_cx|hpkernel"; }
envrec() { echo "-- env: CLD $(git -C "$CLD" log -1 --format='%h %s' | cut -c1-70); hp status [$(git -C "$CLD" status --short hp | tr '\n' ';')]; kernel sha1 $(shasum "$SO" | cut -c1-12)"; }
job() { local tag=$1 dir=$2; shift 2; done_ "$tag" && return 0
        echo "== $tag $(date -u '+%Y-%m-%d %H:%M') UTC: $*" >> "$L/$tag.log"; envrec >> "$L/$tag.log"
        if [ -n "$(dirty)" ]; then echo "== $tag exit 99 (uncommitted change in CLD hp/: $(dirty | tr '\n' ';')) $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$L/$tag.log"; return 99; fi
        (cd "$dir" && gtimeout 10800 nice -n 10 $PY -u "$@" >> "$L/$tag.log" 2>&1)
        local rc=$?
        echo "== $tag exit $rc$([ $rc = 124 ] && echo ' (stopped at the 3 h limit)') $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$L/$tag.log"; }
# E1: the real-axis zero re-found in this environment
mu()    { local tag=$1; shift; job "$tag" "$HERE" ghost_mu.py --hp=$CLD/hp --threads=2 --tag=$tag --out=$RAW/ghost_mu3.jsonl "$@"; }
# the line Im mu = omega
axis()  { local tag=$1 st=$2; shift 2; job "$tag" "$CLD/hp" modecount_test_axis.py "$st" --s=1 --threads=2 \
            --out=$RAW/axis_s1.jsonl --allow_unraw "$@"; }
count() { local tag=$1 st=$2; shift 2; job "$tag" "$CLD/hp" modecount_test.py "$st" --threads=2 --data=$C --allow_unraw "$@"; }
S=$CLD/hp/data/states
# aLL (pick_edges.py on the stored real-axis zero, 2026-11-26) and the planted pairs are fixed in RAW/README.md
case $1 in
  A)  mu D10_M48_m1_i30 $S/up_D10.0000.npz --lo=-0.76 --hi=-0.60 --h=0.04
      count cnt_D10 $S/up_D10.0000.npz --boxes=B0,U0,full --aLL=-0.5961 --tag=D10.0000_M48
      count cnt_D10_plant $S/up_D10.0000.npz --boxes=B0,U0,full --aLL=-0.5961 --tag=D10.0000_M48 --plant=-0.45+0.6j
      mu D20.04_M32_m1_i30 $S/up2_D20.0437.npz --lo=-0.96 --hi=-0.80 --h=0.04
      count cnt_D20.04 $S/up2_D20.0437.npz --boxes=B0,U0,full --aLL=-0.7943 --tag=D20.0437_M32
      count cnt_D20.04_plant $S/up2_D20.0437.npz --boxes=B0,U0,full --aLL=-0.7943 --tag=D20.0437_M32 --plant=-0.65+0.8j
      ;;
  B)  mu D12_M48_m1_i30 $S/up_D12.0000.npz --lo=-0.82 --hi=-0.66 --h=0.04
      count cnt_D12 $S/up_D12.0000.npz --boxes=B0,U0,full --aLL=-0.7145 --tag=D12.0000_M48
      count cnt_D12_plant $S/up_D12.0000.npz --boxes=B0,U0,full --aLL=-0.7145 --tag=D12.0000_M48 --plant=-0.55+0.7j
      # rung 0: the start of section 6.2's refinement (control; X1), rung 1: the record, then the controls of C5
      axis ax0_D10 $S/up_D10.0000.npz --lo=-1.05 --hi=-1.00 --h=0.05
      axis ax0_D12 $S/up_D12.0000.npz --lo=-1.05 --hi=-1.00 --h=0.05
      axis ax_D10 $S/up_D10.0000.npz --lo=-1.05 --hi=-1.00 --h=0.005
      axis ax_D12 $S/up_D12.0000.npz --lo=-1.05 --hi=-1.00 --h=0.005
      axis axc_D16.76 $S/up2_D16.7600.npz --lo=-1.04 --hi=-1.00 --h=0.005
      axis axc_D20.04 $S/up2_D20.0437.npz --lo=-1.04 --hi=-1.00 --h=0.005
      ;;
  R2) axis ax2_D10 $S/up_D10.0000.npz --lo=-1.05 --hi=-1.00 --h=0.002
      axis ax2_D12 $S/up_D12.0000.npz --lo=-1.05 --hi=-1.00 --h=0.002
      ;;
esac
