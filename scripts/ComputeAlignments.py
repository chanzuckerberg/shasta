#!/usr/bin/python3

import shasta
import GetConfig
import sys
import ast

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessKmers()
a.accessMarkers()
a.accessAlignmentCandidates()

# Do the computation.
a.computeAlignments(
    alignmentMethod = int(config['Align']['alignMethodForReadGraph']
    maxMarkerFrequency = int(config['Align']['maxMarkerFrequency']),
    maxSkip = int(config['Align']['maxSkip']),
    minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount']),
    minAlignedFraction = int(config['Align']['minAlignedFraction']),
    maxTrim = int(config['Align']['maxTrim']),
    matchScore = int(config['Align']['matchScore']),
    mismatchScore = int(config['Align']['mismatchScore']),
    gapScore = int(config['Align']['gapScore']),
    suppressContainments = ast.literal_eval(config['Align']['suppressContainments'])
    )
    
    
    

