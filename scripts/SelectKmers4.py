#!/usr/bin/python3

import shasta

# Read the config file.
import GetConfig
config = GetConfig.getConfig()

a = shasta.Assembler()

# For now, use a hardwired distance threshold.
a.selectKmers4(
    k = int(config['Kmers']['k']), 
    markerDensity = float(config['Kmers']['probability']),
    distanceThreshold = 1000)

