#!/usr/bin/python3

import glob
import os
import shutil
import sys
import shasta

# Get from the arguments the list of input fasta files and check that they all exist.
helpMessage = "This script restores a run from directories dataOnDisk and DataOnDisk."
if not len(sys.argv)==1:
    print(helpMessage)
    exit(1)

if not os.path.lexists('data'):
    raise Exception('Missing: data. Use SetupRunDirectory.py to set up the run directory.')
if not os.path.lexists('Data'):
    raise Exception('Missing: Data. Use SetupRunDirectory.py to set up the run directory.')

if not os.path.lexists('dataOnDisk'):
    raise Exception('Missing: dataOnDisk.')
if not os.path.lexists('DataOnDisk'):
    raise Exception('Missing: DataOnDisk.')
    
if glob.glob('data/*'):
    raise Exception('The data directory is not empty. Nothing was restored.')
if glob.glob('Data/*'):
    raise Exception('The Data directory is not empty. Nothing was restored.')

# Restore the data directory.    
inputFileNames = glob.glob('dataOnDisk/*')
for inputFileName in inputFileNames:
    lastSlashPosition = inputFileName.rfind('/')
    outputFileName = 'data/' + inputFileName[lastSlashPosition+1:]
    shasta.mappedCopy(inputFileName, outputFileName)

# Copy the Data directory.
# We cannot use regular copy commands because
# this is on the huge page filesystem.
inputFileNames = glob.glob('DataOnDisk/*')
for inputFileName in inputFileNames:
    lastSlashPosition = inputFileName.rfind('/')
    outputFileName = 'Data/' + inputFileName[lastSlashPosition+1:]
    shasta.mappedCopy(inputFileName, outputFileName)

