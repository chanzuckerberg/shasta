#!/usr/bin/python3

import shasta
import GetConfig

config = GetConfig.getConfig()

a = shasta.Assembler()

a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphEdges(True, True)

a.createMarkerGraphSecondaryEdges(
    secondaryEdgeMaxSkip = int(config['MarkerGraph']['secondaryEdges.maxSkip']))


