#!/usr/bin/python3

import shasta

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices(readWriteAccess = True)
a.accessMarkerGraphReverseComplementVertex(readWriteAccess = True)
a.cleanupDuplicateMarkers(
    duplicateMarkersPattern1Threshold = 0.5,
    pattern1CreateNewVertices = False)
