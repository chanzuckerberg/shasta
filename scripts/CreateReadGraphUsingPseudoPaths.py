#!/usr/bin/python3

import os
import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessAlignmentData()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphEdgeLists()

# Create the read graph.
a.createReadGraphUsingPseudoPaths(
     matchScore = int(config['Assembly']['iterative.pseudoPathAlignMatchScore']),
     mismatchScore = int(config['Assembly']['iterative.pseudoPathAlignMismatchScore']),
     gapScore = int(config['Assembly']['iterative.pseudoPathAlignGapScore']),
     mismatchSquareFactor = int(config['Assembly']['iterative.mismatchSquareFactor']),
     minScore = int(config['Assembly']['iterative.minScore']),
     maxAlignmentCount = int(config['Assembly']['iterative.maxAlignmentCount']))

