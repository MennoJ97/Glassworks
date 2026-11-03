#!/bin/sh
# The page's step with smaller pixels over the same area: six starts at 128, three at 256, one at 384.
cat <<'L' | xargs -P 10 -L1 sh -c 't=$0; n=$1; shift; env TRACE=1 "$@" node genfix3.mjs $n $t >> log/mix/$t.log 2> log/mix/trace-$t-$n.txt'
px384 1 PIXELS=384
px256 1 PIXELS=256
px256 2 PIXELS=256
px256 3 PIXELS=256
px128 1 PIXELS=128
px128 2 PIXELS=128
px128 3 PIXELS=128
px128 4 PIXELS=128
px128 5 PIXELS=128
px128 6 PIXELS=128
L
echo done > log/mix/DONE8
