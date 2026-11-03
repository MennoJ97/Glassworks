#!/bin/sh
# A move schedule by stage, with coefficient mixing and with the page's step, and lost designs kept for a check.
cat <<'L' | while read t n e; do echo "$t $n $e"; done | xargs -P 11 -L1 sh -c 't=$0; n=$1; shift; env TRACE=1 "$@" node genfix3.mjs $n $t >> log/mix/$t.log 2> log/mix/trace-$t-$n.txt'
keepr4 4 MIX=1 REACH=4 KEEP=1
keepr4 5 MIX=1 REACH=4 KEEP=1
keepr4 8 MIX=1 REACH=4 KEEP=1
keepr4 1 MIX=1 REACH=4 KEEP=1
keepbase 1 KEEP=1
L
cat <<'L' | while read t e; do for n in 1 2 3 4 5 6 7 8 9; do echo "$t $n $e"; done; done | xargs -P 11 -L1 sh -c 't=$0; n=$1; shift; env TRACE=1 "$@" node genfix3.mjs $n $t >> log/mix/$t.log 2> log/mix/trace-$t-$n.txt'
schedmix MIX=1 SCHEDULE=1,1,4,4,4,4,2,1
sched SCHEDULE=1,1,4,4,4,4,2,1
L
echo done > log/mix/DONE4
