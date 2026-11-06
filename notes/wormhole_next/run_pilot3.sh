#!/bin/bash
# Queue e of the pilot: D 5 in a second node layout (the 11-node layout that removed D 4's artifacts), made by
# choptuik-large-d's relayout.py into the scratch folder given as $3, then its real-axis zero. usage: run_pilot3.sh e PID DIR
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
CLD=/Users/nathan/Source/choptuik-large-d
RAW=/Users/nathan/Source/choptuik-large-d-raw
PY=$CLD/.venv/bin/python
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$HERE/logs" "$3"
while kill -0 "$2" 2>/dev/null; do sleep 20; done
echo "== relayout $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$HERE/logs/D5_M48_n11_relayout.log"
(cd $CLD/hp && nice -n 10 $PY -u relayout.py $RAW/T79-win/down_D5.0000.npz "$3/n11_D5.0000_M48.npz" \
   --nodes 0.01,0.02,0.05,0.1,0.16,0.23,0.3,0.4,0.5,0.6,0.7 --kc 80 --ks 200 --threads 2 --allow_unraw) \
   >> "$HERE/logs/D5_M48_n11_relayout.log" 2>&1
echo "== relayout exit $? $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$HERE/logs/D5_M48_n11_relayout.log"
tag=D5_M48_n11_m1
echo "== $tag $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$HERE/logs/$tag.log"
(cd "$HERE" && nice -n 10 $PY -u ghost_mu.py --hp=$CLD/hp --threads=2 --tag=$tag "$3/n11_D5.0000_M48.npz" \
   --lo=-0.48 --hi=-0.40 --h=0.04 >> "$HERE/logs/$tag.log" 2>&1)
echo "== $tag exit $? $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$HERE/logs/$tag.log"
