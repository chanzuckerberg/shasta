#!/usr/bin/python3

import shasta

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices(readWriteAccess = True)
a.cleanupDuplicateMarkers(duplicateCoverageRatioThreshold = 0.5)
