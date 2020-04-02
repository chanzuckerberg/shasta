#!/usr/bin/python3

import shasta
import GetConfig
import ast

# Read the config file.
config = GetConfig.getConfig()

# Create the Assembler.
a = shasta.Assembler()

# Set up the consensus caller.
a.setupConsensusCaller(config['Assembly']['consensusCaller'])

a.accessKmers()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphEdgeLists()
a.accessMarkerGraphConsensus()
a.assemble()



