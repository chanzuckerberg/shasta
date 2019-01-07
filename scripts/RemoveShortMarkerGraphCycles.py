#!/usr/bin/python3

import shasta
import GetConfig

config = GetConfig.getConfig()

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(accessEdgesReadWrite = True)
a.removeShortMarkerGraphCycles(
    maxLength = int(config['MarkerGraph']['shortCycleLengthThreshold']))


