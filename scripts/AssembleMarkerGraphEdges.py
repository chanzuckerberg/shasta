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

# Access what we need.
a.accessKmers()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(accessEdgesReadWrite = True)
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()

# Do it.
a.assembleMarkerGraphEdges(
    markerGraphEdgeLengthThresholdForConsensus =
    int(config['Assembly']['markerGraphEdgeLengthThresholdForConsensus']),
    useMarginPhase = useMarginPhase,
    storeCoverageData = ast.literal_eval(config['Assembly']['storeCoverageData']))



