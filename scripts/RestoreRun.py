#!/usr/bin/python3

import glob
import os
import shutil
import sys
import shasta

# Get from the arguments the list of input fasta files and check that they all exist.
helpMessage = "This script restores a run to memory from directory DataOnDisk."
if not len(sys.argv)==1:
    print(helpMessage)
    exit(1)

if not os.path.lexists('Data'):
    raise Exception('Missing: Data. Use SetupRunDirectory.py to set up the run directory.')

if not os.path.lexists('DataOnDisk'):
    raise Exception('Missing: DataOnDisk.')
    
if glob.glob('Data/*'):
    raise Exception('The Data directory is not empty. Nothing was restored.')

# Copy the Data directory.
# We cannot use regular copy commands because
# this is on the huge page filesystem.
inputFileNames = glob.glob('DataOnDisk/*')
for inputFileName in inputFileNames:
    lastSlashPosition = inputFileName.rfind('/')
    outputFileName = 'Data/' + inputFileName[lastSlashPosition+1:]
    shasta.mappedCopy(inputFileName, outputFileName)

