#!/usr/bin/python3

import shasta

import GetConfig
config = GetConfig.getConfig()

a = shasta.Assembler()
a.accessMarkers()
a.accessAlignmentDataReadWrite()
a.accessMarkerGraphEdges()
a.accessMarkerGraphVertices()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdgeLists()
a.accessAssemblyGraphEdges()

a.analyzeAssemblyGraphBubbles(debug=True)

a.createReadGraphMode1(
    maxAlignmentCount = int(config['ReadGraph']['maxAlignmentCount']))


