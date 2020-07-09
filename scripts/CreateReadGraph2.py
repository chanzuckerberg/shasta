#!/usr/bin/python3

import os
import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessAlignmentData()
a.accessCompressedAlignments()

# Create the global read graph using creation method 2.
a.createReadGraph2()


