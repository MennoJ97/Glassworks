#!/bin/sh
# Second-order SSP, queued after the unseen check, on starts 1 to 9.
until [ -f log/fix/DONE4 ]; do sleep 5; done
cat <<'L' | while read t e; do for n in 1 2 3 4 5 6 7 8 9; do echo "$t $n $e"; done; done | xargs -P 11 -L1 sh -c 't=$0; n=$1; shift; env TRACE=1 "$@" node genfix.mjs $n $t >> log/fix/$t.log 2> log/fix/trace-$t-$n.txt'
sp2 SSP=1 SSP2=1
sp2inf SSP=inf SSP2=1
sp2m SSP=1 SSP2=1 CONSTRAIN=1 PREDICT=1
sp2infm SSP=inf SSP2=1 CONSTRAIN=1 PREDICT=1
sp2allm SSP=allinf SSP2=1 CONSTRAIN=1 PREDICT=1
L
echo done > log/fix/DONE5
