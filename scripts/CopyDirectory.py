#!/usr/bin/python3

import glob
import shasta
import sys

helpMessage = """
This can be used to copy all files in a directory
to the huge page filesystem.
The regular cp command does not work (but it works to copy
the other way around, from the huge page filesystem).

All files in the input directory must have size equal to a multiple of the page
size of the huge page filesystem. This will always be the case
if the files were oroginally copied from a huge page filesystem
with the same page size.

Invoke with two arguments:
- The path for the input directory.
- The path for the output directory.
Both directories must exist. 
The second one will be a symbolic link to the huge page filesystem.

"""

if not len(sys.argv) == 3:
    print(helpMessage)
    exit(1)
    
inputName = sys.argv[1]
outputName = sys.argv[2]

inputFileNames = glob.glob(inputName + '/*')

for inputFileName in inputFileNames:
    lastSlashPosition = inputFileName.rfind('/')
    outputFileName = outputName + '/' + inputFileName[lastSlashPosition+1:]
    shasta.mappedCopy(inputFileName, outputFileName)


