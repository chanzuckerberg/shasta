#!/usr/bin/python3

import shasta
import GetConfig
import sys


# Check that there are no arguments.
if not len(sys.argv)==1:
    print(helpMessage)
    exit(1)
    
# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessKmers()
a.accessMarkers()

# Do the computation.
a.flagPalindromicReads(
    maxSkip = int(config['Reads']['palindromicReads.maxSkip']),
    maxMarkerFrequency = int(config['Reads']['palindromicReads.maxMarkerFrequency']),
    alignedFractionThreshold = float(config['Reads']['palindromicReads.alignedFractionThreshold']),
    nearDiagonalFractionThreshold = float(config['Reads']['palindromicReads.nearDiagonalFractionThreshold']),
    deltaThreshold = int(config['Reads']['palindromicReads.deltaThreshold']))

