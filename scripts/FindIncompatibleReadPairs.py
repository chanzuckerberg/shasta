#!/usr/bin/python3

import shasta
import GetConfig

readId = 0  # For now

# Read the config file.
config = GetConfig.getConfig()

# Create the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()

# Find incompatible read pairs
incompatiblePairs = a.findIncompatibleReadPairs(
    readId = readId,
    onlyConsiderLowerReadIds = False,
    skipReadGraphEdges = True)
    
print('Found the following incompatible read pairs:')
for p in incompatiblePairs:
    print(p.readIds[0], p.readIds[1], p.isSameStrand)

