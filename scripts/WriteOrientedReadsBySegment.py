#!/usr/bin/python3

import shasta

a = shasta.Assembler()

a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphEdgeLists()

a.gatherOrientedReadsByAssemblyGraphEdge()
a.writeOrientedReadsByAssemblyGraphEdge()




