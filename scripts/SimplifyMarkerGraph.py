#!/usr/bin/python3

import shasta
import GetConfig

config = GetConfig.getConfig()

simplifyString = config['MarkerGraph']['simplifyMaxLength']
if simplifyString:
	simplifyList = [int(s) for s in simplifyString.split(',')]
else:
    simplifyList = []

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(accessEdgesReadWrite = True)
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphReverseComplementEdge()
a.simplifyMarkerGraph(
    maxLength = simplifyList,
    debug = True)


