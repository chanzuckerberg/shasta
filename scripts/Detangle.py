#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphEdgeLists()
a.detangle()

