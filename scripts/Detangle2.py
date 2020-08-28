#!/usr/bin/python3

import shasta

import GetConfig
config = GetConfig.getConfig()

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphEdgeLists()
a.detangle2(
	diagonalReadCountMin = int(config['Assembly']['detangle.diagonalReadCountMin']),
	offDiagonalReadCountMax = int(config['Assembly']['detangle.offDiagonalReadCountMax']),
	offDiagonalRatio = float(config['Assembly']['detangle.offDiagonalRatio']))

