#!/usr/bin/python3

import Nanopore2
import sys

helpMessage = """
This computes a local marker graph for a set of oriented reads
specified as (ReadId, Strand) pairs (Python tuples)..

The (ReadId, Strand) pairs are read one per line from file 
OrientedReads.txt.

Invoke with one argument: the minimum coverage for a leaf to be kept.
"""

if not len(sys.argv) == 2:
    print(helpMessage)
    exit(1)
minCoverage = int(sys.argv[1])

# Read the oriented reads.
orientedReads = []
for line in open('OrientedReads.txt', 'r'):
    tokens = line.split(' ')
    if not len(tokens) == 2:
        print('OrientedReads.txt must have two space separated tokens per line.')
        exit(1)
    readId = int(tokens[0])
    strand = int(tokens[1])
    if not (strand==0 or strand==1):
        print("Strand must be 0 or 1.")
        exit(2)
    orientedReads.append((readId, strand))    


a = Nanopore2.Assembler()
a.accessKmers()
a.accessReadsReadOnly()
a.accessMarkers()
a.createLocalMarkerGraph(
    readIdsAndStrands = orientedReads, 
    alignAllPairs = True,
    alignmentMaxSkip = 30,
    alignmentMaxVertexCountPerKmer = 100,
    minAlignmentLength = 40,
    minCoverage = minCoverage,
    minConsensus = 3
    )


        

