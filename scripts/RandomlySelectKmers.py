#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()

# Generate the k-mers and write them out.
a.randomlySelectKmers(
    k = int(config['Kmers']['k']), 
    probability = float(config['Kmers']['probability']))
a.writeKmers()
