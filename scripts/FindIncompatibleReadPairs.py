#!/usr/bin/python3

import argparse
import shasta
import GetConfig

# Get the argument.
parser = argparse.ArgumentParser(description=
    'Find reads that have bad induced alignments with a given read.')    
parser.add_argument('--readId', type=int, required=True)
arguments = parser.parse_args()
readId = arguments.readId


# Read the config file.
config = GetConfig.getConfig()

# Create the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessDirectedReadGraphReadOnly()

# Find incompatible read pairs
incompatiblePairs = a.findIncompatibleReadPairs(
    readId = readId,
    onlyConsiderLowerReadIds = False,
    skipReadGraphEdges = True)
    
print('Found the following incompatible read pairs:')
for p in incompatiblePairs:
    print(p.readIds[0], p.readIds[1], p.isSameStrand)

