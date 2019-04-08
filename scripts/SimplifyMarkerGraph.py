#!/usr/bin/python3

import shasta
import GetConfig

config = GetConfig.getConfig()

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(accessEdgesReadWrite = True)
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphReverseComplementEdge()
a.simplifyMarkerGraph(
    maxLength = [int(s) for s in config['MarkerGraph']['simplifyMaxLength'].split(',')],
    debug = True)


