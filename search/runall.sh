#!/bin/bash
# complete tri-Cayley sweep over Z_q: all sorted size triples a<=b<=c with a+b+c<=q+1
q=$1
for a in $(seq 1 $q); do for b in $(seq $a $q); do for c in $(seq $b $q); do
  if [ $((a+b+c)) -le $((q+1)) ]; then python3 fsearch.py $q $a $b $c; fi
done; done; done
