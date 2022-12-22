#!/bin/bash

set -e

cmake .
make

valgrind ./counter -k 17 -m 11 -f data/debug/test.fa -o data/debug/test_count.kff
kfftools validate -v -i data/debug/test_count.kff