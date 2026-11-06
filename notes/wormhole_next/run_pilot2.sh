#!/bin/bash
# Queues c and d of the pilot (see run_pilot.sh). usage: run_pilot2.sh c|d [pid to wait for]
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
if [ $# -ge 2 ]; then while kill -0 "$2" 2>/dev/null; do sleep 20; done; fi
case $1 in
  c)  # controls at D 5 (planted zero; the wrong operator TEST_MUT), then the larger D on the mode counts' M 32 states
      run D5_M48_plant $RAW/T79-win/down_D5.0000.npz --lo=-0.30 --hi=-0.06 --h=0.04 --plant=-0.17
      run D5_M48_m2 $RAW/T79-win/down_D5.0000.npz --lo=-1.0 --hi=-0.04 --h=0.04 --mode=2
      run D8.25_M32_m1 $CLD/hp/data/modecount_state_D8.2500_M32.npz --lo=-2.0 --hi=-0.04 --h=0.06
      run D16.76_M32_m1 $CLD/hp/data/states/up2_D16.7600.npz --lo=-2.0 --hi=-0.04 --h=0.06
      run D30.07_M32_m1 $CLD/hp/data/states/up6_D30.0688.npz --lo=-2.0 --hi=-0.04 --h=0.06
      run D45.05_M32_m1 $CLD/hp/data/states/up7_D45.0507.npz --lo=-2.0 --hi=-0.04 --h=0.06
      ;;
  d)  # D 6 at the second resolution, D 7, and D 8.25 at M 48
      run D6_M64_m1 $RAW/T79-win/m64_D6.0000.npz --lo=-0.70 --hi=-0.40 --h=0.06
      run D7_M48_m1 $RAW/T79-win/down_D7.0000.npz --lo=-1.2 --hi=-0.04 --h=0.06
      run D8.25_M48_m1 $CLD/hp/data/states/up_D8.2500.npz --lo=-1.2 --hi=-0.04 --h=0.06
      ;;
esac
