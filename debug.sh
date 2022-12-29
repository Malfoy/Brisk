#!/bin/bash

set -e

cmake .
make

valgrind ./counter -k 63 -m 5 -b 4 -f data/debug/test.fa -o data/debug/test_count.kff
# kfftools validate -v -i data/debug/test_count.kff