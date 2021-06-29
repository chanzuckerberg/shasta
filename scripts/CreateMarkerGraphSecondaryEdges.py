#!/usr/bin/python3

import shasta
import GetConfig

config = GetConfig.getConfig()

a = shasta.Assembler()

a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphEdges(True, True)

# ********************** Expose this as an expotion when code stabilizes
a.createMarkerGraphSecondaryEdges(secondaryEdgeMaxSkip = 100)


