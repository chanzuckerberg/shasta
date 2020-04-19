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
a.accessAlignmentData()
a.accessReadGraph()
a.accessReadFlags()

# Do the computation.
a.createMarkerGraphVertices(
    alignMethod = int(config['Align']['alignMethod']),
    maxMarkerFrequency = int(config['Align']['maxMarkerFrequency']),
    maxSkip = int(config['Align']['maxSkip']),
    matchScore = int(config['Align']['matchScore']),,
    mismatchScore = int(config['Align']['mismatchScore']),,
    gapScore = int(config['Align']['gapScore']),
    downsamplingFactor = float(config['Align']['downsamplingFactor']),
    bandExtend = int(config['Align']['bandExtend']),
    readGraphCreationMethod = int(config['ReadGraph']['creationMethod']),
    minCoverage = int(config['MarkerGraph']['minCoverage']),
    maxCoverage = int(config['MarkerGraph']['maxCoverage']))

