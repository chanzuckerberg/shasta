#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

a = shasta.Assembler()
a.accessKmers()
a.accessReadsReadOnly()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphEdgeLists()
a.setupConsensusCaller(config['Assembly']['consensusCaller'])
a.assemble()



