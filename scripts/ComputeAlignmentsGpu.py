#!/usr/bin/python3

import shasta
import GetConfig
import sys

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessKmers()
a.accessMarkers()
a.accessAlignmentCandidates()

# Do the computation.
a.computeAlignmentsGpu(
    maxMarkerFrequency = int(config['Align']['maxMarkerFrequency']),
    maxSkip = int(config['Align']['maxSkip']),
    minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount']),
    maxTrim = int(config['Align']['maxTrim']))

