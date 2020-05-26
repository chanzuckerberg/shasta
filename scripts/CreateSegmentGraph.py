#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.createSegmentGraph()

while True:
	chainId = int(input("Enter chain id for GFA coloring:\n"))
	a.colorGfaBySegmentGraphChain(chainId)
	





