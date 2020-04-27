#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphEdgeLists()
a.createCompressedAssemblyGraph()

while True:
	s = input("Enter GFA segment id or oriented read id for coloring: ")
	a.colorCompressedAssemblyGraph(s)
	



