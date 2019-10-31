#!/usr/bin/python3

import shasta
import GetConfig
import sys

helpMessage = """
This computes a marker alignment of two oriented reads.

Invoke with four arguments: readId0, strand0, readId1, strand1.
"""

# Get the arguments.
if not len(sys.argv) == 5:
    print(helpMessage)
    exit(1)
readId0 = int(sys.argv[1]);
strand0 = int(sys.argv[2]);
readId1 = int(sys.argv[3]);
strand1 = int(sys.argv[4]);

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessKmers()
a.accessMarkers()

# For convenience, write the markers sorted by position.
a.writeMarkers(readId=readId0, strand=strand0, 
    fileName = 'Markers-ByPosition-0.csv')
a.writeMarkers(readId=readId1, strand=strand1, 
    fileName = 'Markers-ByPosition-1.csv')

# Compute the alignment.
a.alignOrientedReads1(
    readId0 = readId0, strand0 = strand0,
    readId1 = readId1, strand1 = strand1)

