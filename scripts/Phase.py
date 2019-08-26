#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Create the Assembler.
a = shasta.Assembler()

# Access what we need.
a.accessMarkers();
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphEdgeLists()

# Do the work.
a.createPhasingGraph(
    phasingSimilarityThreshold = float(config['Phasing']['phasingSimilarityThreshold']),
    maxNeighborCount = int(config['Phasing']['maxNeighborCount']))



