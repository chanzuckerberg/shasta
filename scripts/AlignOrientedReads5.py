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
a.alignOrientedReads5(
    readId0 = arguments.readId0, strand0 = arguments.strand0,
    readId1 = arguments.readId1, strand1 = arguments.strand1,
    deltaX = 200,
    deltaY = 10,
    matchScore = int(config['Align']['matchScore']),
    mismatchScore = int(config['Align']['mismatchScore']),
    gapScore = int(config['Align']['gapScore']))
    
    
    
    

