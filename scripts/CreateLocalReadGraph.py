#!/usr/bin/python3

import shasta
import GetConfig
import sys

helpMessage = """
Invoke with three arguments, the id and strand of the read to start from,
and the distance to which to local read graph should extend.
"""

# Get the arguments.
if not len(sys.argv) == 4:
    print(helpMessage)
    exit(1)
readId = int(sys.argv[1])
strand = int(sys.argv[2])
distance = int(sys.argv[3])
if not (strand==0 or strand==1):
    print('Invalid strand')
    print(helpMessage)
    exit(1)   

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessKmers()
a.accessReadsReadOnly()
a.accessReadNamesReadOnly()
a.accessMarkers()
a.accessOverlaps()
a.accessAlignmentData()

# Do the computation.
a.createLocalReadGraph(
    readId = readId,
    strand = strand,
    minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount']),
    maxTrim = int(config['Align']['maxTrim']),
    distance = distance
    )

