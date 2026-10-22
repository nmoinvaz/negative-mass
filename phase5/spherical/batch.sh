#!/bin/bash
# usage: batch.sh jobfile ; each line: NAME STAR h R T [k=v ...]; skips finished runs (final.npz present)
cd /home/claude/negative-mass/phase5/spherical
while read -r line; do
  [ -z "$line" ] && continue
  name=$(echo $line | cut -d' ' -f1)
  [ -f runs/$name.final.npz ] && continue
  python3 run.py $line >> runs/batch.log 2>&1
done < $1
