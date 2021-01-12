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

alignOptions = shasta.AlignOptions()
alignOptions.alignMethod = int(config['Align']['alignMethod'])
alignOptions.maxSkip = int(config['Align']['maxSkip'])
alignOptions.maxDrift = int(config['Align']['maxDrift'])
alignOptions.maxTrim = int(config['Align']['maxTrim'])
alignOptions.maxMarkerFrequency = int(config['Align']['maxMarkerFrequency'])
alignOptions.minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount'])
alignOptions.minAlignedFraction = float(config['Align']['minAlignedFraction'])
alignOptions.matchScore = int(config['Align']['matchScore'])
alignOptions.mismatchScore = int(config['Align']['mismatchScore'])
alignOptions.gapScore = int(config['Align']['gapScore'])
alignOptions.downsamplingFactor = float(config['Align']['downsamplingFactor'])
alignOptions.bandExtend = int(config['Align']['bandExtend'])
alignOptions.maxBand = int(config['Align']['maxBand'])
alignOptions.suppressContainments = ast.literal_eval(config['Align']['suppressContainments'])
alignOptions.sameChannelReadAlignmentSuppressDeltaThreshold = \
    int(config['Align']['sameChannelReadAlignment.suppressDeltaThreshold'])
alignOptions.align4DeltaX = int(config['Align']['align4.deltaX'])
alignOptions.align4DeltaY = int(config['Align']['align4.deltaY'])
alignOptions.align4MinEntryCountPerCell = int(config['Align']['align4.minEntryCountPerCell'])
alignOptions.align4MaxDistanceFromBoundary = int(config['Align']['align4.maxDistanceFromBoundary'])

# Do the computation.
a.computeAlignments(alignOptions, 0)

    
    

