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

# Figure out if we should use marginPhase, and if so set it up.
useMarginPhase = ast.literal_eval(config['Assembly']['useMarginPhase'])
if useMarginPhase:
    a.setupMarginPhase()

a.accessKmers()
a.accessReadsReadOnly()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphEdgeLists()
a.accessMarkerGraphVertexRepeatCounts()
a.accessMarkerGraphEdgeConsensus()
a.assemble()



