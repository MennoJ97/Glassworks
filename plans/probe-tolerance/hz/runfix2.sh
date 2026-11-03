#!/bin/sh
until [ -f log/fix/DONE ]; do sleep 5; done
for v in "o 1 '' 1" "p '' 1 1" "op 1 1 1" "opr5 1 1 0.5"; do :; done
printf '%s\n' "o 1 0 1" "p 0 1 1" "op 1 1 1" "opr5 1 1 0.5" | while read t o pr re; do for n in 1 2 3 4 5 6 7 8 9; do echo "$t $o $pr $re $n"; done; done | xargs -P 11 -L1 sh -c 'env TRACE=1 $( [ "$1" = 1 ] && echo ORTH=1 ) $( [ "$2" = 1 ] && echo PREDICT=1 ) REACH=$3 node genfix.mjs $4 $0 >> log/fix/$0.log 2> log/fix/trace-$0-$4.txt'
echo done > log/fix/DONE2
