#!/usr/bin/python3

import Nanopore2
import Nanopore2GetConfig
import sys

helpMessage = """
This computes a local marker graph in one of two ways:

- Case A: If called without arguments, reads a set of oriented reads
  and creates a local marker graph for them, 
  computing alignments for all possible pairs. 
  The (ReadId, Strand) pairs are read one per line from file 
  OrientedReads.txt.
- Case B: If called with three arguments, the three arguments are
  a read id, strand, and distance. It first creates a local read
  graph starting from the given oriented read and extending out
  up to the given distance (just like CreateLocalReadGraph.py). 
  Then, it creates the local marker graph
  corresponding to this local read graph.  

"""

# Read the config file.
config = Nanopore2GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = Nanopore2.Assembler()
a.accessKmers()
a.accessReadsReadOnly()
a.accessReadNamesReadOnly()
a.accessMarkers()
a.accessOverlaps()
a.accessAlignmentInfos()



if len(sys.argv) == 1:

    # Case A.

    # Read the oriented reads that will be used to create
    # the local marker graph.
    orientedReads = []
    for line in open('OrientedReads.txt', 'r'):
        tokens = line.split(' ')
        if not len(tokens) == 2:
            print('OrientedReads.txt must have two space separated tokens per line.')
            exit(1)
        readId = int(tokens[0])
        strand = int(tokens[1])
        if not (strand==0 or strand==1):
            print("Strand must be 0 or 1.")
            exit(2)
        orientedReads.append((readId, strand))    


    # Do the computation.
    a.createLocalMarkerGraph(
        readIdsAndStrands = orientedReads, 
        alignAllPairs = True,
        alignmentMaxSkip = int(config['Align']['maxSkip']),
        alignmentMaxVertexCountPerKmer = int(config['Align']['maxVertexCountPerKmer']),
        minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount']),
        minCoverage = int(config['MarkerGraph']['minCoverage']),
        minConsensus = int(config['MarkerGraph']['minConsensus'])
        )



elif len(sys.argv) == 4:

    # Case B.

    readId = int(sys.argv[1])
    strand = int(sys.argv[2])
    distance = int(sys.argv[3])
    if not (strand==0 or strand==1):
        print('Invalid strand')
        print(helpMessage)
        exit(1)   

    a.createLocalMarkerGraph(
        readId = readId,
        strand = strand,
        minFrequency = int(config['MinHash']['minFrequency']),
        minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount']),
        maxTrim = int(config['Align']['maxTrim']),
        distance = distance,
        alignmentMaxSkip = int(config['Align']['maxSkip']),
        alignmentMaxVertexCountPerKmer = int(config['Align']['maxVertexCountPerKmer']),
        minCoverage = int(config['MarkerGraph']['minCoverage']),
        minConsensus = int(config['MarkerGraph']['minConsensus']))


    
else:
    
    # If getting here, the number of arguments is not right
    # for either case A or case B.
    print(helpMessage)
    exit(1)

        

