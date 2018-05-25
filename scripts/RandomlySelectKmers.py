#!/usr/bin/python3

import Nanopore2
import Nanopore2GetConfig

# Read the config file.
config = Nanopore2GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = Nanopore2.Assembler()

# Generate the k-mers and write them out.
a.randomlySelectKmers(
    k = int(config['Kmers']['k']), 
    probability = float(config['Kmers']['probability']))
a.writeKmers()
