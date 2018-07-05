#!/usr/bin/python3

import shasta
import GetConfig
import sys

helpMessage = """
This computes alignments of an oriented read with all reads
for which we have an Overlap.

Invoke with two arguments: readId, strand.
"""

# Get the arguments.
if not len(sys.argv) == 3:
    print(helpMessage)
    exit(1)
readId = int(sys.argv[1]);
strand = int(sys.argv[2]);

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessKmers()
a.accessReadsReadOnly()
a.accessReadNamesReadOnly()
a.accessMarkers()
a.accessOverlaps()

# Compute the alignments.
a.alignOverlappingOrientedReads(
    readId=readId, strand=strand,
    maxSkip = int(config['Align']['maxSkip']),
    maxVertexCountPerKmer = int(config['Align']['maxVertexCountPerKmer']),
    minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount']),
    maxTrim = int(config['Align']['maxTrim'])
    )

