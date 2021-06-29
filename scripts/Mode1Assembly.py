#!/usr/bin/python3

"""
This runs Mode 1 assembly from initial creation of the read graph to the end.
It is useful for debugging.
It assumes that everything up to alignment computation was done.
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
    
# Create the marker graph.    
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

# Create the assembly graph.
a.createAssemblyGraphEdges()
a.createAssemblyGraphVertices()

# Analyze bubbles in the assembly graph.
a.analyzeAssemblyGraphBubbles()

# Create the final read graph.
a.createReadGraphMode1(
    maxAlignmentCount = int(config['ReadGraph']['maxAlignmentCount']))
    
# Remove the initial marker graph and assembly graph.
a.removeMarkerGraph()
a.removeAssemblyGraph()

# ******* EXPOSE THESE WHEN CODE STABILIZES (ALSO in main.cpp)
minVertexCoverageFinal = 4
minVertexCoveragePerStrandFinal = 1
minEdgeCoverageFinal = 4
minEdgeCoveragePerStrandFinal = 1
secondaryEdgeMaxSkip = 100

# Create the final marker graph.
a.createMarkerGraphVertices(
    minCoverage = minVertexCoverageFinal,
    maxCoverage = int(config['MarkerGraph']['maxCoverage']),
    minCoveragePerStrand = minVertexCoveragePerStrandFinal,
    allowDuplicateMarkers = ast.literal_eval(config['MarkerGraph']['allowDuplicateMarkers']),
    peakFinderMinAreaFraction = float(config['MarkerGraph']['peakFinder.minAreaFraction']),
    peakFinderAreaStartIndex = int(config['MarkerGraph']['peakFinder.areaStartIndex']))    
a.findMarkerGraphReverseComplementVertices()
a.createMarkerGraphEdgesStrict(
    minEdgeCoverage = minEdgeCoverageFinal,
    minEdgeCoveragePerStrand = minEdgeCoveragePerStrandFinal)
a.findMarkerGraphReverseComplementEdges()
a.computeMarkerGraphCoverageHistogram()

# Add secondary edges.
a.createMarkerGraphSecondaryEdges(
    secondaryEdgeMaxSkip = secondaryEdgeMaxSkip)

# Bubble/superbubble removal.
simplifyString = config['MarkerGraph']['simplifyMaxLength']
if simplifyString:
	simplifyList = [int(s) for s in simplifyString.split(',')]
else:
    simplifyList = []
a.simplifyMarkerGraph(
    maxLength = simplifyList,
    debug = True)

# Create the final assembly graph.
a.createAssemblyGraphEdges()
a.createAssemblyGraphVertices()

# Assemble.
a.assembleMarkerGraphVertices()
a.assembleMarkerGraphEdges(
    markerGraphEdgeLengthThresholdForConsensus =
    int(config['Assembly']['markerGraphEdgeLengthThresholdForConsensus']),
    storeCoverageData = ast.literal_eval(config['Assembly']['storeCoverageData']))
a.assemble()
a.computeAssemblyStatistics()

# Wtite out the assembly.
a.writeGfa1('Assembly.gfa')
a.writeGfa1BothStrands('Assembly-BothStrands.gfa')
a.writeFasta('Assembly.fasta')




