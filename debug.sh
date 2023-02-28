#!/bin/bash

set -e

cmake .
make

DATA=data/rongeur.fa
OUT=data/rongeur
# DATA=data/debug/test.fa
# OUT=data/debug/test
K=63


# echo ">Test" > data/test.fa
# head $DATA -n 200 | tail -n 100 >> data/test.fa
# head $DATA -n 200 > data/test.fa
DATA=data/test.fa


echo "kmc -k$K -ci0 -okff -fm $DATA ${OUT}_kmc /tmp"
kmc -k$K -fm -ci0 -okff $DATA ${OUT}_kmc /tmp > ${OUT}_kmc.stdout
./counter -k 63 -m 21 -b 9 -f $DATA -o ${OUT}_brisk.kff --mode 2 > ${OUT}_brisk.stdout

grep "No. of unique counted k-mers" ${OUT}_kmc.stdout
grep "nb kmers" ${OUT}_brisk.stdout
