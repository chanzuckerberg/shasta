#!/usr/bin/python3

import Nanopore2
import sys

a = Nanopore2.Assembler()
a.accessReadsReadOnly()
a.accessKmers()
a.accessMarkers()
a.accessOverlaps()
a.computeAllAlignments(maxSkip=30, maxVertexCountPerKmer=100)

