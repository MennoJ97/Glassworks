#!/bin/sh
# Finer meshes on the top designs, then m and spm on unseen starts, judged by drawn disorder.
cp des2/m-*.json des2/spm-*.json des/
for id in m-8 m-2 m-3 spm-3 spm-4 m-1; do echo $id; done | xargs -P 6 -I{} sh -c 'node check.mjs {} finer > log/fix/chk-{}.json 2>&1'
for t in "m CONSTRAIN=1 PREDICT=1" "spm SSP=1 CONSTRAIN=1 PREDICT=1"; do echo "$t"; done | while read t e; do for n in $(seq 101 118); do echo "$t $n $e"; done; done > /dev/null
cat <<'L' | while read t e; do for n in $(seq 101 118); do echo "$t $n $e"; done; done | xargs -P 11 -L1 sh -c 't=$0; n=$1; shift; env "$@" node genfix.mjs $n $t >> log/fix/$t-u.log 2>&1'
m CONSTRAIN=1 PREDICT=1
spm SSP=1 CONSTRAIN=1 PREDICT=1
L
cp des2/m-1??.json des2/spm-1??.json des/
for n in $(seq 101 118); do echo "m $n"; echo "spm $n"; done | xargs -P 11 -L1 sh -c 'test -f des/$0-$1.json && node judge.mjs $0 $1 0.1 4 > log/judge2/$0-$1.log 2>&1'
for n in $(seq 101 118); do echo "m-$n"; echo "spm-$n"; done | xargs -P 11 -I{} sh -c 'node check.mjs {} > log/chk-{}.s 2>&1'
echo done > log/fix/DONE4
