#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Create the Assembler.
a = shasta.Assembler()

# Set up the consensus caller.
a.setupConsensusCaller(config['Assembly']['consensusCaller'])

# Access what we need.
a.accessKmers()
a.accessReadsReadOnly()
a.accessMarkers()
a.accessMarkerGraphVertices()

# Do it.
a.assembleMarkerGraphVertices()



