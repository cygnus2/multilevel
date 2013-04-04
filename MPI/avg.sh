#!/bin/bash
for i in 40 41 42
do
  awk '{print $1}' fort.$i > test
#  awk '{print $2}' fort.$i > test
  ./avg test 1 1 -1 
  rm test
done
exit 0
