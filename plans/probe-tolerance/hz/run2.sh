#!/bin/sh
# Unseen starts: plain and rc2 optimized, then judged, then more crystal draws.
for n in $(seq 101 118); do echo "p none 2 $n"; echo "rc2 corr 2 $n"; done | xargs -P 11 -L1 sh -c 'ROB=$1 ROB_NM=$2 node genrob.mjs $3 $0 >> log/rob/$0-u.log 2>&1'
for n in $(seq 101 118); do echo "p $n"; echo "rc2 $n"; done | xargs -P 11 -L1 sh -c 'test -f des/$0-$1.json && node judge.mjs $0 $1 0.1 4 > log/judge2/$0-$1.log 2>&1'
for k in $(seq 7 18); do echo "2 $k"; echo "5 $k"; done | xargs -P 11 -L1 sh -c 'node xdraws.mjs $0 $1 >> log/xdraws.log 2>&1'
echo done > log/judge2/DONE
