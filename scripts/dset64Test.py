#!/usr/bin/python3

helpMessage = """
This runs a unit test for dset64.hpp.

Invoke with 5 arguments:
- The number of items (vertices).
- The number of union operations (edges).
- The number of threads.
- The number of union operations per batch.
- The random seed.
"""


import Nanopore2
import sys

if not len(sys.argv)==6:
    print(helpMessage)
    exit(1);
    
n = int(sys.argv[1])
m = int(sys.argv[2])
threadCount = int(sys.argv[3])
batchSize = int(sys.argv[4])
seed = int(sys.argv[5])


Nanopore2.dset64Test(
    n = n,
    m = m,
    threadCount = threadCount,
    batchSize = batchSize,
    seed = seed)

