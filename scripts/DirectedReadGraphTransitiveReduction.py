#!/usr/bin/python3

import os
import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessDirectedReadGraphReadWrite()

# Do the transitive reduction.
# Use hardwired parameters for now.
# These should be added as additional command line options
# in the ReadGraph section when this code stabilizes.
a.directedReadGraphTransitiveReduction(
    offsetTolerance0 = 100,
    offsetTolerance1 = 0.1)


