#!/usr/bin/python3

import shasta
import GetConfig
import sys

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessReadFlags()
a.accessKmers()
a.accessMarkers()
a.accessAlignmentData()
a.accessReadGraph()

# Create vertices of the marker graph.
a.createMarkerGraphVertices(
    alignMethod = int(config['Align']['alignMethod']),
    maxMarkerFrequency = int(config['Align']['maxMarkerFrequency']),
    maxSkip = int(config['Align']['maxSkip']),
    maxDrift = int(config['Align']['maxDrift']),
    matchScore = int(config['Align']['matchScore']),
    mismatchScore = int(config['Align']['mismatchScore']),
    gapScore = int(config['Align']['gapScore']),
    downsamplingFactor = float(config['Align']['downsamplingFactor']),
    bandExtend = int(config['Align']['bandExtend']),
    readGraphCreationMethod = int(config['ReadGraph']['creationMethod']),
    minCoverage = int(config['MarkerGraph']['minCoverage']),
    maxCoverage = int(config['MarkerGraph']['maxCoverage']),
    minCoveragePerStrand = int(config['MarkerGraph']['minCoveragePerStrand']))

# Create edges of the marker graph.
a.createMarkerGraphEdges()

# Approximate transitive reduction.
a.transitiveReduction(
    lowCoverageThreshold = int(config['MarkerGraph']['lowCoverageThreshold']),
    highCoverageThreshold = int(config['MarkerGraph']['highCoverageThreshold']),
    maxDistance = int(config['MarkerGraph']['maxDistance']),
)

# Prune the strong subgraph of the marker graph.
a.pruneMarkerGraphStrongSubgraph(
    iterationCount = int(config['MarkerGraph']['pruneIterationCount']))

