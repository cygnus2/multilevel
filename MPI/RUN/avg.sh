#!/bin/bash
echo "Pol Loop:"
for i in 40 41 42
do
  awk '{print $1}' fort.$i > test
#  awk '{print $2}' fort.$i > test
  ./avg test 1 1 -1 
  rm test
done

echo "G(1):"
for i in 51 61 71
do
  awk '{print $1}' fort.$i > test
#  awk '{print $2}' fort.$i > test
  ./avg test 1 1 -1
  rm test
done

echo "G(2):"
for i in 52 62 72
do
  awk '{print $1}' fort.$i > test
#  awk '{print $2}' fort.$i > test
  ./avg test 1 1 -1
  rm test
done


exit 0
