#!/usr/bin/python3

import shasta
import GetConfig
import sys

helpMessage = """
This computes a marker alignment of two oriented reads 
using alignment method 4.

Invoke with four arguments: readId0, strand0, readId1, strand1.
"""

# Get the arguments.
import argparse
parser = argparse.ArgumentParser(
    description='This computes a marker alignment of two oriented reads using alignment method 4.')
parser.add_argument('readId0', type=int)
parser.add_argument('strand0', type=int, choices=range(2))
parser.add_argument('readId1', type=int)
parser.add_argument('strand1', type=int, choices=range(2))
arguments = parser.parse_args()

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()

# Compute the alignment.
a.alignOrientedReads4(
    readId0 = arguments.readId0, strand0 = arguments.strand0,
    readId1 = arguments.readId1, strand1 = arguments.strand1,
    deltaX = int(config['Align']['align4.deltaX']),
    deltaY = int(config['Align']['align4.deltaY']),
    minEntryCountPerCell = int(config['Align']['align4.minEntryCountPerCell']),
    maxDistanceFromBoundary = int(config['Align']['align4.maxDistanceFromBoundary']),
    minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount']),
    minAlignedFraction = float(config['Align']['minAlignedFraction']),
    maxSkip = int(config['Align']['maxSkip']),
    maxDrift = int(config['Align']['maxDrift']),
    maxTrim = int(config['Align']['maxTrim']),
    matchScore = int(config['Align']['matchScore']),
    mismatchScore = int(config['Align']['mismatchScore']),
    gapScore = int(config['Align']['gapScore']))
    
    
    
    

