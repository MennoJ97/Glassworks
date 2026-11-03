#!/bin/sh
# A trust-region move limit, with coefficient mixing and with the page's step, on starts 1 to 9.
cat <<'L' | while read t e; do for n in 1 2 3 4 5 6 7 8 9; do echo "$t $n $e"; done; done | xargs -P 11 -L1 sh -c 't=$0; n=$1; shift; env TRACE=1 "$@" node genfix3.mjs $n $t >> log/mix/$t.log 2> log/mix/trace-$t-$n.txt'
trustmix MIX=1 TRUST=0.75,0.25,1.5,0.5,0.25,4
trust TRUST=0.75,0.25,1.5,0.5,0.25,4
L
echo done > log/mix/DONE3
