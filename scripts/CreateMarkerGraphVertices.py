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
a.accessChimericReadsFlags()

# Do the computation.
a.createMarkerGraphVertices(
    maxMarkerFrequency = int(config['Align']['maxMarkerFrequency']),
    maxSkip = int(config['Align']['maxSkip']),
    minCoverage = int(config['MarkerGraph']['minCoverage']),
    maxCoverage = int(config['MarkerGraph']['maxCoverage']))

