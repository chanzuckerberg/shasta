#!/usr/bin/python3

import Nanopore2
import sys

helpMessage = """
This can be used to copy a file to the huge page filesystem.
The regular cp command does not work (but it works to copy
the other way around, from the huge page filesystem).

The input file must have size equal to a multiple of the page
size of the huge page filesystem. This will always be the case
if the file was oroginally copied from a huge page filesystem
with the same page size.

Invoke with two arguments:
- The path for the input file.
- The path for the output file.

"""

if not len(sys.argv) == 3:
    print(helpMessage)
    exit(1)
    
inputName = sys.argv[1]
outputName = sys.argv[2]
Nanopore2.mappedCopy(inputName, outputName)

