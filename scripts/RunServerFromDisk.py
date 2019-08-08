#!/usr/bin/python3

import os
import shasta
import GetConfig

# Find the path to the docs directory.
thisScriptPath = os.path.realpath(__file__)
thisScriptDirectory = os.path.dirname(thisScriptPath)
thisScriptParentDirectory = os.path.dirname(thisScriptDirectory)
docsDirectory = thisScriptParentDirectory + '/docs'

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler.
a = shasta.Assembler(
    largeDataFileNamePrefix='DataOnDisk/')
a.accessAllSoft()
a.setupConsensusCaller(config['Assembly']['consensusCaller'])

a.setDocsDirectory(docsDirectory)
a.explore(port=17100, localOnly=False, sameUserOnly=False)


