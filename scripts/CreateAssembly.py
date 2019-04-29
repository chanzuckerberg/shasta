#!/usr/bin/python3

import shasta
import GetConfig
import os
import sys

# Check that we have what we need.
if not os.path.lexists('Data'):
    raise Exception('Missing: Data. Use SetupRunDirectory.py to set up the run directory.')
if not os.path.lexists('shasta.conf'):
    raise Exception('Missing: configuration file shasta.conf. Sample available in shasta-install/conf.')


# Read the config file.
config = GetConfig.getConfig()

# Create the assembler.
a = shasta.Assembler(createNew=True)


