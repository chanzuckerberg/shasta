#!/usr/bin/python3

import shasta
import GetConfig



# Read the config file.
config = GetConfig.getConfig()

# Create the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessDirectedReadGraphReadOnly()

# Create the conflict read graph.
# Use hardwired parameter values for now.
a.createConflictReadGraph(
    maxOffsetSigma = 50,
    maxTrim = 100,
    maxSkip = 100)


