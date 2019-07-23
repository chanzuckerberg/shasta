#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()


# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(accessEdgesReadWrite=True)
a.accessMarkerGraphReverseComplementEdge()
a.transitiveReduction(
    lowCoverageThreshold = int(config['MarkerGraph']['lowCoverageThreshold']),
    highCoverageThreshold = int(config['MarkerGraph']['highCoverageThreshold']),
    maxDistance = int(config['MarkerGraph']['maxDistance']),
    edgeMarkerSkipThreshold = int(config['MarkerGraph']['edgeMarkerSkipThreshold'])
    )


