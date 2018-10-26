#!/usr/bin/python3

import shasta
import GetConfig
import sys

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessReadsReadOnly()
a.accessKmers()
a.accessMarkers()
a.accessAlignmentData()
a.accessReadGraph()

# Do the computation.
a.createMarkerGraphVertices(
    maxVertexCountPerKmer = int(config['Align']['maxVertexCountPerKmer']),
    maxSkip = int(config['Align']['maxSkip']),
    minCoverage = int(config['MarkerGraph']['minCoverage']))

