#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessAlignmentData()
a.accessCompressedAlignments()
a.accessMarkerGraphVertices(readWriteAccess = True)
a.accessMarkerGraphReverseComplementVertex(readWriteAccess = True)
a.cleanupDuplicateMarkers(
    minCoverage = int(config['MarkerGraph']['minCoverage']),
    minCoveragePerStrand = int(config['MarkerGraph']['minCoveragePerStrand']),
    duplicateMarkersPattern1Threshold = float(config['MarkerGraph']['duplicateMarkersPattern1Threshold']),
    pattern1CreateNewVertices = False,
    pattern2CreateNewVertices = False)
