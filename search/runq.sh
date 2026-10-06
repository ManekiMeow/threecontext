#!/bin/bash
for q in $@; do
for sz in "3 3 3" "3 3 4" "2 4 4" "3 4 4" "3 3 5" "2 3 6" "2 4 5" "4 4 4" "3 4 5" "3 3 6" "2 5 5"; do
  timeout 5400 python3 fsearch.py $q $sz
done; done
