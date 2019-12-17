#!/usr/bin/python3

import shasta
import GetConfig



# Read the config file.
config = GetConfig.getConfig()

# Create the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessDirectedReadGraphReadOnly()

# Find incompatible read pairs
a.findAllIncompatibleReadPairs()
    


