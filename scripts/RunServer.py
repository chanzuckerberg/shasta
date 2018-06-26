#!/usr/bin/python3

import Nanopore2


# Initialize the assembler and access what we need.
a = Nanopore2.Assembler()
a.accessReadsReadOnly()
a.accessReadNamesReadOnly()
a.accessKmers()
a.accessMarkers()
a.accessOverlaps()
a.accessAlignmentInfos()
a.accessGlobalMarkerGraph()
a.explore()


