#!/bin/sh
# Second-order SSP with the held step, on unseen starts, judged by drawn disorder.
for n in $(seq 101 118); do echo "$n"; done | xargs -P 11 -I{} sh -c 'env SSP=1 SSP2=1 CONSTRAIN=1 PREDICT=1 node genfix.mjs {} sp2m >> log/fix/sp2m-u.log 2>&1'
cp des2/sp2m-1??.json des/
for n in $(seq 101 118); do echo "$n"; done | xargs -P 11 -I{} sh -c 'test -f des/sp2m-{}.json && node judge.mjs sp2m {} 0.1 4 > log/judge2/sp2m-{}.log 2>&1'
cp des2/sp2m-[1-9].json des/
for n in $(seq 101 118); do echo "sp2m-$n"; done | xargs -P 11 -I{} sh -c 'node check.mjs {} > log/chk-{}.s 2>&1'
for id in sp2m-1 sp2m-6 sp2m-7; do echo $id; done | xargs -P 3 -I{} sh -c 'node check.mjs {} finer > log/fix/chk-{}.json 2>&1'
echo done > log/fix/DONE6
