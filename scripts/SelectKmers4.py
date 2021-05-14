#!/usr/bin/python3

import shasta

# Read the config file.
import GetConfig
config = GetConfig.getConfig()

# Create the assembler.
a = shasta.Assembler()

# select k-mers.
a.selectKmers4(
    k = int(config['Kmers']['k']), 
    markerDensity = float(config['Kmers']['probability']),
    distanceThreshold = int(config['Kmers']['distanceThreshold']))

