#!/bin/sh
# Coefficient mixing and stage-change prediction, against the page's own step, on starts 1 to 9.
cat <<'L' | while read t e; do for n in 1 2 3 4 5 6 7 8 9; do echo "$t $n $e"; done; done | xargs -P 11 -L1 sh -c 't=$0; n=$1; shift; env TRACE=1 "$@" node genfix3.mjs $n $t >> log/mix/$t.log 2> log/mix/trace-$t-$n.txt'
base3 X=0
mix MIX=1
stage STAGE=1
mixstage MIX=1 STAGE=1
L
echo done > log/mix/DONE1
