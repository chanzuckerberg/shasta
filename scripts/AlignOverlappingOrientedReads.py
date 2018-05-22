#!/usr/bin/python3

import Nanopore2
import sys

helpMessage = """
This computes alignments of an oriented read with all reads
for which we have an Overlap.

Invoke with two arguments: readId, strand.
"""

if not len(sys.argv) == 3:
    print(helpMessage)
    exit(1)

readId = int(sys.argv[1]);
strand = int(sys.argv[2]);


a = Nanopore2.Assembler()
a.accessKmers()
a.accessReadsReadOnly()
a.accessReadNamesReadOnly()
a.accessMarkers()
a.accessOverlaps()


a.alignOverlappingOrientedReads(
    readId=readId, strand=strand,
    maxSkip=30,
    minAlignedMarkerCount = 40,
    maxTrim = 200)

