#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()


# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkerGraphConnectivity(accessEdgesReadWrite=True)
a.flagMarkerGraphEdges(
    minCoverage = int(config['MarkerGraph']['minCoverage']),
    maxPathLength = int(config['MarkerGraph']['maxPathLength']),
    )


