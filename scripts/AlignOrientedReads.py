#!/usr/bin/python3

import Nanopore2
import sys

helpMessage = """
This computes a marker alignment of two oriented reads.

Invoke with four arguments: readId0, strand0, readId1, strand1.
"""

if not len(sys.argv) == 5:
    print(helpMessage)
    exit(1)

readId0 = int(sys.argv[1]);
strand0 = int(sys.argv[2]);
readId1 = int(sys.argv[3]);
strand1 = int(sys.argv[4]);


a = Nanopore2.Assembler()
a.accessKmers()
a.accessReadsReadOnly()
a.accessMarkers()

# For convenience, write the markers sorted by position.
a.writeMarkers(readId=readId0, strand=strand0, 
    fileName = 'Markers-ByPosition-0.csv')
a.writeMarkers(readId=readId1, strand=strand1, 
    fileName = 'Markers-ByPosition-1.csv')

a.alignOrientedReads(
    readId0=readId0, strand0=strand0,
    readId1=readId1, strand1=strand1,
    maxSkip=30)

