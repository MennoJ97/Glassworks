#!/bin/sh
# Coefficient mixing with larger move limits, on starts 1 to 9.
cat <<'L' | while read t e; do for n in 1 2 3 4 5 6 7 8 9; do echo "$t $n $e"; done; done | xargs -P 11 -L1 sh -c 't=$0; n=$1; shift; env TRACE=1 "$@" node genfix3.mjs $n $t >> log/mix/$t.log 2> log/mix/trace-$t-$n.txt'
mixr2 MIX=1 REACH=2
mixr4 MIX=1 REACH=4
L
echo done > log/mix/DONE2
