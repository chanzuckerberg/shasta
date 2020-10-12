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
a.accessCompressedAlignments()

# Do the computation.
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
    maxBand = int(config['Align']['maxBand']),
    minCoverage = int(config['MarkerGraph']['minCoverage']),
    maxCoverage = int(config['MarkerGraph']['maxCoverage']),
    minCoveragePerStrand = int(config['MarkerGraph']['minCoveragePerStrand']),
    peakFinderMinAreaFraction = float(config['MarkerGraph']['peakFinder.minAreaFraction']),
    peakFinderAreaStartIndex = int(config['MarkerGraph']['peakFinder.areaStartIndex']))
    

