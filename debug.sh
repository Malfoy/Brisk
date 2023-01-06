#!/bin/bash

set -e

cmake .
make

echo ">test" > data/test.fa
tail data/rongeur.fa -n 4000000 >> data/test.fa
./counter -k 63 -m 21 -b 9 -f data/test.fa -o data/test.kff --mode 2
# ./counter -k 63 -m 21 -b 9 -f data/rongeur.fa -o data/rongeur.fa --mode 2
# ./counter -k 63 -m 21 -b 9 -f data/bacterie.fa -o data/bacterie.fa --mode 2
# ./counter -k 21 -m 7 -b 5 -f data/bacterie.fa -o data/bacterie.fa --mode 2
# valgrind ./counter -k 63 -m 21 -b 9 -f data/debug/test.fa -o data/debug/test_count.kff --mode 2
# kfftools validate -v -i data/debug/test_count.kff