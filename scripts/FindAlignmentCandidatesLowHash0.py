#!/usr/bin/python3

import shasta
import GetConfig
import sys

helpMessage="""
This uses the LowHash method to find alignment candidates.

Invoke without arguments.
"""

# Check that there are no arguments.
if not len(sys.argv)==1:
    print(helpMessage)
    exit(1)
    
# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessKmers()
a.accessReadFlags()
a.accessMarkers()

# Do the computation.
a.findAlignmentCandidatesLowHash0(
    m = int(config['MinHash']['m']), 
    hashFraction = float(config['MinHash']['hashFraction']),
    minHashIterationCount = int(config['MinHash']['minHashIterationCount']), 
    alignmentCandidatesPerRead = float(config['MinHash']['alignmentCandidatesPerRead']), 
    minBucketSize = int(config['MinHash']['minBucketSize']),
    maxBucketSize = int(config['MinHash']['maxBucketSize']),
    minFrequency = int(config['MinHash']['minFrequency']))

