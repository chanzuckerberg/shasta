#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

a = shasta.Assembler()
a.accessAlignmentData()
a.accessReadGraph()
a.removeReadGraphBridges(
	maxDistance = int(config['Assembly']['iterative.bridgeRemovalMaxDistance']))

