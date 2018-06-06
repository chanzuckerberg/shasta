#!/usr/bin/python3

import Nanopore2
import Nanopore2GetConfig
import sys

# Read the config file.
config = Nanopore2GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = Nanopore2.Assembler()
a.accessReadsReadOnly()
a.accessKmers()
a.accessMarkers()
a.accessOverlaps()

# Do the computation.
a.computeAllAlignments(
    minFrequency = int(config['MinHash']['minFrequency']),
    maxVertexCountPerKmer = int(config['Align']['maxVertexCountPerKmer']),
    maxSkip = int(config['Align']['maxSkip']),
    minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount']),
    maxTrim = int(config['Align']['minAlignedMarkerCount']))



