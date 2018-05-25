#!/usr/bin/python3

import Nanopore2
import sys

helpMessage = """
Invoke with three arguments, the id and strand of the read to start from,
and the distance to which to local read graph should extend.
"""

if not len(sys.argv) == 4:
    print(helpMessage)
    exit(1)

readId = int(sys.argv[1])
strand = int(sys.argv[2])
distance = int(sys.argv[3])
if not (strand==0 or strand==1):
    print('Invalid strand')
    print(helpMessage)
    exit(1)   

a = Nanopore2.Assembler()

a.accessKmers()
a.accessReadsReadOnly()
a.accessReadNamesReadOnly()
a.accessMarkers()
a.accessOverlaps()
a.accessAlignmentInfos()

a.createLocalReadGraph(
    readId = readId,
    strand = strand,
    minFrequency = 3,
    minAlignedMarkerCount = 100,
    maxTrim = 200,
    distance = distance
    )

