#!/usr/bin/python3

import shasta
import GetConfig

config = GetConfig.getConfig()

# To get meaningful results from this, use
# the following options when running the assembly,
# to make sure all vertices are generated:
# --MarkerGraph.allowDuplicateMarkers 
# --MarkerGraph.minCoverage 1
# --MarkerGraph.minCoverage 1000000000


a = shasta.Assembler()
a.accessKmers()
a.accessMarkers()
a.accessMarkerGraphVertices()

a.vertexCoverageStatisticsByKmerId()


