#!/usr/bin/python3

import os
import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessAlignmentDataReadWrite()

# Create the global read graph.
a.createReadGraph2(
    maxAlignmentCount = int(config['ReadGraph']['maxAlignmentCount']),
    markerCountPercentile = float(config['ReadGraph']['markerCountPercentile']),
    alignedFractionPercentile = float(config['ReadGraph']['alignedFractionPercentile']),
    maxSkipPercentile = float(config['ReadGraph']['maxSkipPercentile']),
    maxDriftPercentile = float(config['ReadGraph']['maxDriftPercentile']),
    maxTrimPercentile = float(config['ReadGraph']['maxTrimPercentile']))


