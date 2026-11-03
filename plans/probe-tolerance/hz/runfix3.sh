#!/bin/sh
until [ -f log/fix/DONE2 ]; do sleep 5; done
cat <<'L' | while read t e; do for n in 1 2 3 4 5 6 7 8 9; do echo "$t $n $e"; done; done | xargs -P 11 -L1 sh -c 't=$0; n=$1; shift; env TRACE=1 "$@" node genfix.mjs $n $t >> log/fix/$t.log 2> log/fix/trace-$t-$n.txt'
m CONSTRAIN=1 PREDICT=1
mr5 CONSTRAIN=1 PREDICT=1 REACH=0.5
mr25 CONSTRAIN=1 PREDICT=1 REACH=0.25
sp SSP=1
spinf SSP=inf
spm SSP=1 CONSTRAIN=1 PREDICT=1
spinfm SSP=inf CONSTRAIN=1 PREDICT=1
L
echo done > log/fix/DONE3
