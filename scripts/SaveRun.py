#!/usr/bin/python3

import os
import shutil
import sys


# Get from the arguments the list of input fasta files and check that they all exist.
helpMessage = "This script copies a run to directory DataOnDisk."
if not len(sys.argv)==1:
    print(helpMessage)
    exit(1)

if not os.path.lexists('Data'):
    raise Exception('Data does not exist.')

if os.path.lexists('DataOnDisk'):
    raise Exception('DataOnDisk already exists. Remove before running this script.')
    
shutil.copytree('Data', 'DataOnDisk')

