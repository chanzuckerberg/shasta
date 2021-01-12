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

# Create vertices of the marker graph.
a.createMarkerGraphVertices(
    minCoverage = int(config['MarkerGraph']['minCoverage']),
    maxCoverage = int(config['MarkerGraph']['maxCoverage']),
    minCoveragePerStrand = int(config['MarkerGraph']['minCoveragePerStrand']),
    peakFinderMinAreaFraction = float(config['MarkerGraph']['peakFinder.minAreaFraction']),
    peakFinderAreaStartIndex = int(config['MarkerGraph']['peakFinder.areaStartIndex']))
a.findMarkerGraphReverseComplementVertices()


# Create edges of the marker graph.
a.createMarkerGraphEdges()
a.findMarkerGraphReverseComplementEdges()

# Approximate transitive reduction.
a.transitiveReduction(
    lowCoverageThreshold = int(config['MarkerGraph']['lowCoverageThreshold']),
    highCoverageThreshold = int(config['MarkerGraph']['highCoverageThreshold']),
    maxDistance = int(config['MarkerGraph']['maxDistance']),
    edgeMarkerSkipThreshold = int(config['MarkerGraph']['edgeMarkerSkipThreshold']),
)

# Prune the strong subgraph of the marker graph.
a.pruneMarkerGraphStrongSubgraph(
    iterationCount = int(config['MarkerGraph']['pruneIterationCount']))

