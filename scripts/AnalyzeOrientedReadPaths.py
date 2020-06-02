#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

a = shasta.Assembler()
a.accessMarkers()
a.accessReadGraph()
a.accessReadFlags()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphEdgeLists()
a.analyzeOrientedReadPaths(
	readGraphCreationMethod = int(config['ReadGraph']['creationMethod']))

	





