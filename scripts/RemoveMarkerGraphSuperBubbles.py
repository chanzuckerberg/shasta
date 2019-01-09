#!/usr/bin/python3

import shasta
import GetConfig

config = GetConfig.getConfig()

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(accessEdgesReadWrite = True)
a.removeMarkerGraphSuperBubbles(
    maxLength = int(config['MarkerGraph']['superBubbleLengthThreshold']))


