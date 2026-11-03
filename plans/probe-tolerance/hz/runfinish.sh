#!/bin/sh
# Whole-pixel projection in the last stages, with the page's step and with mixing at four times the move.
cat <<'L' | while read t e; do for n in 1 2 3 4 5 6 7 8 9; do echo "$t $n $e"; done; done | xargs -P 11 -L1 sh -c 't=$0; n=$1; shift; env TRACE=1 "$@" node genfix3.mjs $n $t >> log/mix/$t.log 2> log/mix/trace-$t-$n.txt'
fin2 FINISH=2
fin2mixr4 FINISH=2 MIX=1 REACH=4
L
echo done > log/mix/DONE5
