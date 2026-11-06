#!/bin/bash
# Step 2 of notes/wormhole_next.md (section 6): the ghost sector finished on choptuik-large-d's code and states,
# read-only. Two queues of 2 threads at nice 10. A job whose log ends in "exit 0" is skipped, so a restart repeats only
# unfinished jobs; modecount_test.py keeps every evaluation in counts/*_raw.jsonl and reuses it.
# usage: step2.sh P|Q
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
CLD=/Users/nathan/Source/choptuik-large-d
RAW=/Users/nathan/Source/choptuik-large-d-raw
PY=$CLD/.venv/bin/python
export PYTHONDONTWRITEBYTECODE=1
L=$HERE/logs2; C=$HERE/counts
mkdir -p "$L" "$C"
done_() { [ -f "$L/$1.log" ] && tail -1 "$L/$1.log" | grep -q "exit 0"; }
job() { local tag=$1 dir=$2; shift 2; done_ "$tag" && return 0
        echo "== $tag $(date -u '+%Y-%m-%d %H:%M') UTC: $*" >> "$L/$tag.log"
        (cd "$dir" && nice -n 10 $PY -u "$@" >> "$L/$tag.log" 2>&1)
        echo "== $tag exit $? $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$L/$tag.log"; }
mu()    { local tag=$1; shift; job "$tag" "$HERE" ghost_mu.py --hp=$CLD/hp --threads=2 --tag=$tag --out=$HERE/ghost_mu2.jsonl "$@"; }
axis()  { local tag=$1 st=$2; shift 2; job "$tag" "$CLD/hp" modecount_test_axis.py "$st" --s=1 --threads=2 \
            --out=$HERE/axis_s1.jsonl --allow_unraw "$@"; }
count() { local tag=$1 st=$2; shift 2; job "$tag" "$CLD/hp" modecount_test.py "$st" --threads=2 --data=$C --allow_unraw "$@"; }
# aLL for a count: just right of the real zero recorded in ghost_mu2.jsonl / ghost_mu.jsonl, away from candidate poles
edge()  { $PY "$HERE/pick_edges.py" --hp=$CLD/hp "$1" "$($PY -c "
import json,sys
z=[max(x['mu'] for x in r['zeros'] if x['ok']) for f in ('$HERE/ghost_mu.jsonl','$HERE/ghost_mu2.jsonl')
   for r in map(json.loads, open(f)) if r['tag']=='$2' and r['mode']==1 and r['plant'] is None]
print(z[-1])")" | sed 's/.*aLL \([-0-9.]*\) .*/\1/'; }
S=$CLD/hp/data/states
case $1 in
  P)  # winding counts at D 5, 6, 7 from just right of the real zero to aR, |Im| < 1.5 omega, each with a planted pair;
      # then D 4's box between its real zero and the origin
      count cnt_D5 $RAW/T79-win/down_D5.0000.npz --boxes=B0,U0,full --aLL=-0.4230 --tag=D5.0000_M48
      count cnt_D5_plant $RAW/T79-win/down_D5.0000.npz --boxes=B0,U0,full --aLL=-0.4230 --tag=D5.0000_M48 --plant=-0.30+0.6j
      count cnt_D6 $RAW/T79-win/down_D6.0000.npz --boxes=B0,U0,full --aLL=-0.4252 --tag=D6.0000_M48
      count cnt_D6_plant $RAW/T79-win/down_D6.0000.npz --boxes=B0,U0,full --aLL=-0.4252 --tag=D6.0000_M48 --plant=-0.35+0.6j
      count cnt_D7 $RAW/T79-win/down_D7.0000.npz --boxes=B0,U0,full --aLL=-0.4773 --tag=D7.0000_M48
      count cnt_D7_plant $RAW/T79-win/down_D7.0000.npz --boxes=B0,U0,full --aLL=-0.4773 --tag=D7.0000_M48 --plant=-0.40+0.6j
      count box_D4 $CLD/hp/data/lowD_states/lowDn_D4.0000_M48.npz --boxes=B0,U0 --aLL=-0.3451 --tag=D4.0000_M48n
      count box_D4_plant $CLD/hp/data/lowD_states/lowDn_D4.0000_M48.npz --boxes=B0,U0 --aLL=-0.3451 --tag=D4.0000_M48n --plant=-0.25+0.5j
      ;;
  Q)  # the other parity sector (Im mu = omega) at the pilot's D; boxes left of CLD's counts at D 8.25-45; mu(D) at
      # D 52-60; counts at D 52, 60, 67
      axis ax_D6 $RAW/T79-win/down_D6.0000.npz --lo=-1.6 --hi=-0.02 --h=0.05
      axis ax_D8.25 $CLD/hp/data/modecount_state_D8.2500_M32.npz --lo=-1.6 --hi=-0.02 --h=0.04
      axis ax_D16.76 $S/up2_D16.7600.npz --lo=-1.6 --hi=-0.02 --h=0.04
      axis ax_D30.07 $S/up6_D30.0688.npz --lo=-1.6 --hi=-0.02 --h=0.04
      axis ax_D45.05 $S/up7_D45.0507.npz --lo=-1.6 --hi=-0.02 --h=0.04
      axis ax_D67.01 $CLD/hp/data/large_d/up70_D67.0115.npz --lo=-1.6 --hi=-0.02 --h=0.04
      count box_D8.25 $CLD/hp/data/modecount_state_D8.2500_M32.npz --boxes=B0,U0 --aLL=-0.5325 --tag=D8.2500_M32
      count box_D8.25_plant $CLD/hp/data/modecount_state_D8.2500_M32.npz --boxes=B0,U0 --aLL=-0.5325 --tag=D8.2500_M32 --plant=-0.45+0.6j
      count box_D16.76 $S/up2_D16.7600.npz --boxes=B0,U0 --aLL=-0.7513 --tag=D16.7600_M32
      count box_D16.76_plant $S/up2_D16.7600.npz --boxes=B0,U0 --aLL=-0.7513 --tag=D16.7600_M32 --plant=-0.60+0.8j
      count box_D30.07 $S/up6_D30.0688.npz --boxes=B0,U0 --aLL=-0.9224 --tag=D30.0688_M32
      count box_D30.07_plant $S/up6_D30.0688.npz --boxes=B0,U0 --aLL=-0.9224 --tag=D30.0688_M32 --plant=-0.70+0.8j
      count box_D45.05 $S/up7_D45.0507.npz --boxes=B0,U0 --aLL=-0.8956 --tag=D45.0507_M32
      count box_D45.05_plant $S/up7_D45.0507.npz --boxes=B0,U0 --aLL=-0.8956 --tag=D45.0507_M32 --plant=-0.70+0.8j
      mu D52.04_M32_m1 $S/up9_D52.0414.npz --lo=-1.4 --hi=-0.04 --h=0.06
      mu D55.04_m1 $S/up10_D55.0351.npz --lo=-1.4 --hi=-0.04 --h=0.06
      mu D58.01_m1 $S/up12_D58.0060.npz --lo=-1.4 --hi=-0.04 --h=0.06
      mu D60.00_m1 $S/up13_D60.0020.npz --lo=-1.4 --hi=-0.04 --h=0.06
      count cnt_D67 $CLD/hp/data/large_d/up70_D67.0115.npz --boxes=B0,U0,full --aLL=-0.9096 --tag=D67.0115_M32
      count cnt_D67_plant $CLD/hp/data/large_d/up70_D67.0115.npz --boxes=B0,U0,full --aLL=-0.9096 --tag=D67.0115_M32 --plant=-0.70+0.8j
      a=$(edge $S/up9_D52.0414.npz D52.04_M32_m1)
      count cnt_D52 $S/up9_D52.0414.npz --boxes=B0,U0,full --aLL=$a --tag=D52.0414
      count cnt_D52_plant $S/up9_D52.0414.npz --boxes=B0,U0,full --aLL=$a --tag=D52.0414 --plant=-0.70+0.8j
      a=$(edge $S/up13_D60.0020.npz D60.00_m1)
      count cnt_D60 $S/up13_D60.0020.npz --boxes=B0,U0,full --aLL=$a --tag=D60.0020
      count cnt_D60_plant $S/up13_D60.0020.npz --boxes=B0,U0,full --aLL=$a --tag=D60.0020 --plant=-0.70+0.8j
      axis ax_D7 $RAW/T79-win/down_D7.0000.npz --lo=-1.6 --hi=-0.02 --h=0.05
      axis ax_D10 $S/up_D10.0000.npz --lo=-1.6 --hi=-0.02 --h=0.05
      axis ax_D12 $S/up_D12.0000.npz --lo=-1.6 --hi=-0.02 --h=0.05
      axis ax_D20.04 $S/up2_D20.0437.npz --lo=-1.6 --hi=-0.02 --h=0.04
      ;;
esac
