#!/bin/bash

SCRPT=/home/michael_science/git/partitionfinder/PartitionFinder.py
IND=/home/michael_science/Desktop/TESTFELD/Datasets_partitioning/03_testing
FILESTEM=test3

# WRAP CMD IN VALGRIND/MASSIF AND EXECUTE CMD
valgrind --tool=massif --trace-children=yes --pages-as-heap=yes --time-unit=ms --massif-out-file=${FILESTEM}_massif.out python2 $SCRPT -p 1 $IND 1>${FILESTEM}.log 2>&1

# SUMMARIZE VALGRIND/MASSIF OUTPUT
ms_print --threshold=50.0 ${FILESTEM}_massif.out > ${FILESTEM}_massif.summary
# For interpretation of graph, see: http://valgrind.org/docs/manual/ms-manual.html
