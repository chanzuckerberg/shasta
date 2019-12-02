#!/usr/bin/python3

import argparse
import GetConfig
import os
import shasta

# Get the arguments.
parser = argparse.ArgumentParser(description=
    'Analyze a vertex of the directed read graph.')    
parser.add_argument('--readId', type=int, required=True, metavar='readId')
parser.add_argument('--strand', type=int, choices=range(2), required=True, metavar='strand')
arguments = parser.parse_args()
readId = arguments.readId
strand = arguments.strand

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessDirectedReadGraphReadOnly()
a.accessMarkers()

# Analyze the vertex corresponding to this readId and strand.
a.analyzeDirectedReadGraphVertex(
    readId = readId,
    strand = strand)


