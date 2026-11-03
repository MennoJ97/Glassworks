#!/bin/sh
# More starts with smaller pixels and the dipole's spacing fixed, beside the 384 run.
cat <<'L' | xargs -P 5 -L1 sh -c 't=$0; n=$1; shift; env TRACE=1 "$@" node --max-old-space-size=8000 genfix3a.mjs $n $t >> log/mix/$t.log 2> log/mix/trace-$t-$n.txt'
d256 4 PIXELS=256 DSTEP=1
d256 5 PIXELS=256 DSTEP=1
d128 7 PIXELS=128 DSTEP=1
d128 8 PIXELS=128 DSTEP=1
d128 9 PIXELS=128 DSTEP=1
L
