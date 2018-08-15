#!/usr/bin/python3

import shasta
import GetConfig
import os
import sys

# Check that we have what we need.
if not os.path.lexists('data'):
    raise Exception('Missing: data. Use SetupRunDirectory.py to set up the run directory.')
if not os.path.lexists('Data'):
    raise Exception('Missing: Data. Use SetupRunDirectory.py to set up the run directory.')
if not os.path.lexists('threadLogs'):
    raise Exception('Missing: threadLogs. Use SetupRunDirectory.py to set up the run directory.')
if not os.path.lexists('shasta.conf'):
    raise Exception('Missing: configuration file shasta.conf. Sample available in shasta-install/conf.')


# Read the config file.
config = GetConfig.getConfig()

# Create the assembler.
useRunLengthReadsString = config['Reads']['useRunLengthReads']
if useRunLengthReadsString == 'True':
    useRunLengthReads = True
elif useRunLengthReadsString == 'False':
    useRunLengthReads = False
else:
    raise RuntimeError("Configuration parameter useRunLengthReads in section Reads must be True or False.")
a = shasta.Assembler(useRunLengthReads = useRunLengthReads)


