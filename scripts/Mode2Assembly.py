#!/usr/bin/python3

"""
Mode 2 assembly experiments.
This assumes that everything up to alignment computation was done.
"""

import shasta
import GetConfig

import ast
import os

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessKmers()
a.accessMarkers()
a.accessAlignmentDataReadWrite()
a.accessCompressedAlignments()
a.setupConsensusCaller(config['Assembly']['consensusCaller'])

# Create the read graph.
readGraphCreationMethod = int(config['ReadGraph']['creationMethod'])
if readGraphCreationMethod == 0:
    a.createReadGraph(
        maxAlignmentCount = int(config['ReadGraph']['maxAlignmentCount']),
        maxTrim = int(config['Align']['maxTrim']))

elif readGraphCreationMethod == 2:
    a.createReadGraph2(
        maxAlignmentCount = int(config['ReadGraph']['maxAlignmentCount']),
        markerCountPercentile = int(config['ReadGraph']['markerCountPercentile']),
        alignedFractionPercentile = int(config['ReadGraph']['alignedFractionPercentile']),
        maxSkipPercentile = int(config['ReadGraph']['maxSkipPercentile']),
        maxDriftPercentile = int(config['ReadGraph']['maxDriftPercentile']),
        maxTrimPercentile = int(config['ReadGraph']['maxTrimPercentile']))
else:
    raise ValueError('Invalid value for --ReadGraph.creationMethod.')
    
a.flagCrossStrandReadGraphEdges(
    maxDistance = int(config['ReadGraph']['crossStrandMaxDistance']))
a.flagChimericReads(
    maxChimericReadDistance = int(config['ReadGraph']['maxChimericReadDistance']))
    


# Flag inconsistent alignments, if requested.
flagInconsistentAlignments = ast.literal_eval(config['ReadGraph']['flagInconsistentAlignments'])
if flagInconsistentAlignments:
    a.flagInconsistentAlignments(
        triangleErrorThreshold = config['ReadGraph']['flagInconsistentAlignmentsTriangleErrorThreshold'],
        leastSquareErrorThreshold = config['ReadGraph']['flagInconsistentAlignmentsLeastSquareErrorThreshold'],
        leastSquareMaxDistance = config['ReadGraph']['flagInconsistentAlignmentsLeastSquareMaxDistance'])

# Compute connected components of the read graph.
# These are currently not used.
a.computeReadGraphConnectedComponents(
    int(config['ReadGraph']['minComponentSize']))
    
# Create the marker graph using strict edge creation.    
a.createMarkerGraphVertices(
    minCoverage = int(config['MarkerGraph']['minCoverage']),
    maxCoverage = int(config['MarkerGraph']['maxCoverage']),
    minCoveragePerStrand = int(config['MarkerGraph']['minCoveragePerStrand']),
    allowDuplicateMarkers = ast.literal_eval(config['MarkerGraph']['allowDuplicateMarkers']),
    peakFinderMinAreaFraction = float(config['MarkerGraph']['peakFinder.minAreaFraction']),
    peakFinderAreaStartIndex = int(config['MarkerGraph']['peakFinder.areaStartIndex']))    
a.findMarkerGraphReverseComplementVertices()
a.createMarkerGraphEdgesStrict(
    minEdgeCoverage = int(config['MarkerGraph']['minEdgeCoverage']),
    minEdgeCoveragePerStrand = int(config['MarkerGraph']['minEdgeCoveragePerStrand']))
a.findMarkerGraphReverseComplementEdges()
a.computeMarkerGraphCoverageHistogram()

# Add secondary edges.
a.createMarkerGraphSecondaryEdges(
    secondaryEdgeMaxSkip = 1000000)

# Missing: 
# - Create AssemblyGraph2 (including, later, bubble detection and phasing).
# - Assemble marker graph vertices, edges.
# - Write gfa, fasta.




