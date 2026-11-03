#!/bin/sh
# Starts 1 to 9 again, keeping each run's last optimized variables, for the anisotropic-averaging probe.
cat <<'L' | while read t e; do for n in 1 2 3 4 5 6 7 8 9; do echo "$t $n $e"; done; done | xargs -P 11 -L1 sh -c 't=$0; n=$1; shift; env "$@" node genfix3.mjs $n $t >> log/mix/$t.log 2>> log/mix/err-$t-$n.txt'
kbase KEEP=1
kmixr4 KEEP=1 MIX=1 REACH=4
L
echo done > log/mix/DONE7
