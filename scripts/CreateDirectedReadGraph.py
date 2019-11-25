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

# Create the global read graph.
a.createDirectedReadGraph(maxTrim = int(config['Align']['maxTrim']))


