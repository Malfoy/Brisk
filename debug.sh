#!/bin/bash

cmake . && make
valgrind ./counter -k 17 -f data/debug/test.fa -o data/debug/test_count.kff