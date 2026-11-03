#!/bin/sh
# Smaller pixels again, with the dipole's two points a fixed λ/48 apart, as on the page's mesh.
cat <<'L' | xargs -P 9 -L1 sh -c 't=$0; n=$1; shift; env TRACE=1 "$@" node --max-old-space-size=12000 genfix3a.mjs $n $t >> log/mix/$t.log 2> log/mix/trace-$t-$n.txt'
d256 1 PIXELS=256 DSTEP=1
d256 2 PIXELS=256 DSTEP=1
d256 3 PIXELS=256 DSTEP=1
d128 1 PIXELS=128 DSTEP=1
d128 2 PIXELS=128 DSTEP=1
d128 3 PIXELS=128 DSTEP=1
d128 4 PIXELS=128 DSTEP=1
d128 5 PIXELS=128 DSTEP=1
d128 6 PIXELS=128 DSTEP=1
L
echo done > log/mix/DONE9
