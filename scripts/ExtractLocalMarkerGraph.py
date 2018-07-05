#!/usr/bin/python3

import shasta
import GetConfig
import sys

helpMessage = """
This extracts a local marker graph from the local marker graph.

Invoke with 4 arguments:
- The read id, strand, and ordinal of the marker that defines the start vertex.
- The distance to move away from the start vertex (number of edges on the 
  global marker graph).
"""


if not len(sys.argv) == 5:
    print(helpMessage)
    exit(1)
readId = int(sys.argv[1])
strand = int(sys.argv[2])
ordinal = int(sys.argv[3])
distance = int(sys.argv[4])


# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessReadsReadOnly()
a.accessKmers()
a.accessMarkers()
a.accessOverlaps()
a.accessGlobalMarkerGraph()

# Do it.
a.extractLocalMarkerGraph(
    readId = readId,
    strand = strand,
    ordinal = ordinal,
    distance = distance,
    minCoverage = int(config['MarkerGraph']['minCoverage']),
    minConsensus = int(config['MarkerGraph']['minConsensus']))
    
