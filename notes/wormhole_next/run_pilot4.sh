#!/bin/bash
# Queue f of the pilot: the other half-period parity sector (the line Im mu = omega) at D 5, with choptuik-large-d's
# hp/modecount_test_axis.py, output into this folder. usage: run_pilot4.sh f PID
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
CLD=/Users/nathan/Source/choptuik-large-d
RAW=/Users/nathan/Source/choptuik-large-d-raw
PY=$CLD/.venv/bin/python
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$HERE/logs"
while kill -0 "$2" 2>/dev/null; do sleep 20; done
tag=D5_M48_s1_axis
echo "== $tag $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$HERE/logs/$tag.log"
(cd $CLD/hp && nice -n 10 $PY -u modecount_test_axis.py $RAW/T79-win/down_D5.0000.npz --lo=-1.0 --hi=-0.02 --h=0.04 \
   --s=1 --threads=2 --out="$HERE/modecount_test_axis.jsonl" --allow_unraw) >> "$HERE/logs/$tag.log" 2>&1
echo "== $tag exit $? $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$HERE/logs/$tag.log"
