#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()

# Generate the k-mers and write them out.
a.selectKmersBasedOnFrequency(
    k = int(config['Kmers']['k']), 
    markerDensity = float(config['Kmers']['probability']),
    enrichmentThreshold = float(config['Kmers']['enrichmentThreshold']))

