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

# Access what we need.
a.accessKmers()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(accessEdgesReadWrite = True)
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()

# Do it.
assembleAllEdges = (int(config['Assembly']['mode']) == 2)
a.assembleMarkerGraphEdges(
    markerGraphEdgeLengthThresholdForConsensus =
    int(config['Assembly']['markerGraphEdgeLengthThresholdForConsensus']),
    storeCoverageData = ast.literal_eval(config['Assembly']['storeCoverageData']),
    assembleAllEdges = assembleAllEdges)



