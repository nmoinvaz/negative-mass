#!/bin/bash
# Queue g of the pilot: more D for mu(D) (D 10, 12, 20.04, 60.00, 67.01). usage: run_pilot5.sh g PID
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
CLD=/Users/nathan/Source/choptuik-large-d
PY=$CLD/.venv/bin/python
export PYTHONDONTWRITEBYTECODE=1
run() { local tag=$1; shift; echo "== $tag $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$HERE/logs/$tag.log"
        (cd "$HERE" && nice -n 10 $PY -u ghost_mu.py --hp=$CLD/hp --threads=2 --tag=$tag "$@" >> "$HERE/logs/$tag.log" 2>&1)
        echo "== $tag exit $? $(date -u '+%Y-%m-%d %H:%M') UTC" >> "$HERE/logs/$tag.log"; }
mkdir -p "$HERE/logs"
while kill -0 "$2" 2>/dev/null; do sleep 20; done
run D20.04_M32_m1 $CLD/hp/data/states/up2_D20.0437.npz --lo=-1.4 --hi=-0.04 --h=0.06
run D67.01_M32_m1 $CLD/hp/data/large_d/up70_D67.0115.npz --lo=-1.4 --hi=-0.04 --h=0.06
run D10_M48_m1 $CLD/hp/data/states/up_D10.0000.npz --lo=-1.4 --hi=-0.04 --h=0.06
run D12_M48_m1 $CLD/hp/data/states/up_D12.0000.npz --lo=-1.4 --hi=-0.04 --h=0.06
