#!/usr/bin/python3

import Nanopore2
import sys

helpMessage="""
This uses the MinHash method to find overlapping reads.

Invoke with one argument, log2MinHashBucketCount.

To avoid too many collisions, this should be greater
than the base 2 log of the number of reads plus 3.
"""

if not len(sys.argv)==2:
    print(helpMessage)
    exit(1)
    
log2MinHashBucketCount = int(sys.argv[1])    

a = Nanopore2.Assembler()
a.accessKmers()
a.accessMarkers()
a.findOverlaps(
    m=5, 
    minHashIterationCount=100, 
    log2MinHashBucketCount=log2MinHashBucketCount,
    maxBucketSize = 30,
    minFrequency = 1)

