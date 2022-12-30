#!/bin/bash

set -e

cmake .
make

valgrind ./counter -k 63 -m 11 -b 7 -f data/debug/test.fa -o data/debug/test_count.kff --mode 2
# kfftools validate -v -i data/debug/test_count.kff