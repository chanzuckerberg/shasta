#!/usr/bin/python3

import shasta
import GetConfig

# Read the config file.
config = GetConfig.getConfig()


# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessAlignmentData()
a.accessReadGraph()
a.flagChimericReads(
    maxChimericReadDistance = int(config['ReadGraph']['maxChimericReadDistance']))


