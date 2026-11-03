#!/bin/sh
# Variants on starts 1 to 9: tag settle reach seed.
for v in "base 0 1" "s1 1 1" "s05 0.5 1" "r5 0 0.5" "r25 0 0.25" "s1r5 1 0.5"; do for n in 1 2 3 4 5 6 7 8 9; do echo "$v $n"; done; done | xargs -P 11 -L1 sh -c 'TRACE=1 SETTLE=$1 REACH=$2 node genfix.mjs $3 $0 >> log/fix/$0.log 2> log/fix/trace-$0-$3.txt'
echo done > log/fix/DONE
