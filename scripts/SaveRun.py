#!/usr/bin/python3

import os
import shutil
import sys


# Get from the arguments the list of input fasta files and check that they all exist.
helpMessage = "This script copies a run to directories dataOnDisk and DataOnDisk."
if not len(sys.argv)==1:
    print(helpMessage)
    exit(1)

if not os.path.lexists('data'):
    raise Exception('data does not exist.')
if not os.path.lexists('Data'):
    raise Exception('Data does not exist.')

if os.path.lexists('dataOnDisk'):
    raise Exception('dataOnDisk already exists. Remove before running this script.')
if os.path.lexists('DataOnDisk'):
    raise Exception('DataOnDisk already exists. Remove before running this script.')
    
shutil.copytree('data', 'dataOnDisk')
shutil.copytree('Data', 'DataOnDisk')

