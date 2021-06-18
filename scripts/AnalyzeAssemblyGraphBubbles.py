#!/usr/bin/python3

import shasta

import GetConfig
config = GetConfig.getConfig()

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphEdges()
a.accessMarkerGraphVertices()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdgeLists()
a.accessAssemblyGraphEdges()
a.analyzeAssemblyGraphBubbles()


