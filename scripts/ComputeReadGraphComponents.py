#!/usr/bin/python3

import Nanopore2
import Nanopore2GetConfig
import sys
 
# Read the config file.
config = Nanopore2GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = Nanopore2.Assembler()
a.accessKmers()
a.accessReadsReadOnly()
a.accessMarkers()
a.accessOverlaps()
a.accessAlignmentInfos()

# Do the computation.
a.computeReadGraphComponents(
    minFrequency = int(config['MinHash']['minFrequency']),
    minComponentSize = int(config['ReadGraph']['minComponentSize']),
    minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount']),
    maxTrim = int(config['Align']['maxTrim'])
    )

