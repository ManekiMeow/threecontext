#!/bin/bash
q=$1
for sz in "4 4 5" "3 5 5" "4 5 5" "3 4 6" "4 4 6" "5 5 5" "3 5 6" "2 5 6" "2 6 6"; do
  timeout 5400 python3 fsearch.py $q $sz
done
