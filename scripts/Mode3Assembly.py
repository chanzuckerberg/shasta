#!/usr/bin/python3

"""

This run the final portion of Mode 3 assembly.
It assumes that the marker graph has already been created.

"""

import ast
import shasta
import GetConfig

config = GetConfig.getConfig()

shasta.openPerformanceLog('Mode3Assembly.log')

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphEdges()
a.accessMarkerGraphReverseComplementEdge()
a.accessMarkerGraphConsensus()

a.mode3Assembly()




