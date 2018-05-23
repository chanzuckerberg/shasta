#!/usr/bin/python3

import Nanopore2
import sys

   

a = Nanopore2.Assembler()
a.accessKmers()
a.accessReadsReadOnly()
a.accessMarkers()
a.accessOverlaps()
a.accessAlignmentInfos()
a.computeOverlapGraphComponents(
    minFrequency = 1,
    minComponentSize = 100,
    minAlignedMarkerCount = 40,
    maxTrim = 200
    )

