#!/bin/bash
# The pilot of wormhole_next.md section 3: ghost-sector decay rates on choptuik-large-d's stored states.
# Two queues of 2 threads each at nice 10 (python and states from the paths below; nothing is written there).
# usage: run_pilot.sh a|b
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
CLD=/Users/nathan/Source/choptuik-large-d
RAW=/Users/nathan/Source/choptuik-large-d-raw
PY=$CLD/.venv/bin/python
export PYTHONDONTWRITEBYTECODE=1
run() { local tag=$1; shift; echo "== $tag $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$HERE/logs/$tag.log"
        (cd "$HERE" && nice -n 10 $PY -u ghost_mu.py --hp=$CLD/hp --threads=2 --tag=$tag "$@" >> "$HERE/logs/$tag.log" 2>&1)
        echo "== $tag exit $? $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$HERE/logs/$tag.log"; }
mkdir -p "$HERE/logs"
case $1 in
  a)  # D 4: reproduce -0.3651202372 (choptuik-large-d) and GMG96's -0.362 +- 0.012 (the controls moved to queue c of
      # run_pilot2.sh, at D 5, where an evaluation costs a quarter of D 4's at M 64)
      run D4_M64_m1 $CLD/hp/data/lowD_states/lowDn_D4.0000_M64.npz --lo=-0.40 --hi=-0.30 --h=0.05
      ;;
  b)  # D 5 against MGM26's -0.448, two resolutions; D 6 (not in the notes)
      run D5_M48_m1 $RAW/T79-win/down_D5.0000.npz --lo=-1.0 --hi=-0.04 --h=0.04
      run D5_M64_m1 $RAW/T79-win/m64_D5.0000.npz --lo=-0.48 --hi=-0.40 --h=0.04
      run D6_M48_m1 $RAW/T79-win/down_D6.0000.npz --lo=-1.0 --hi=-0.04 --h=0.04
      ;;
esac
